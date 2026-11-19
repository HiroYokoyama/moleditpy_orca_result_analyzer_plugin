"""Thermochemistry dialog: enthalpy/entropy/free-energy table with CSV export."""

from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QPushButton,
    QApplication,
    QCheckBox,
    QComboBox,
    QLabel,
    QAbstractItemView,
    QFileDialog,
    QMessageBox,
)
import os
import csv
import logging
from .settings import load_section, save_section
from .utils import notify

HEADERS = ["Property", "Value", "Unit"]
_SECTION = "section"

#: Energy units offered for display: hartree multiplier and decimals shown.
#: Factors are CODATA 2018.
ENERGY_UNITS = {
    "Eh": (1.0, 8),
    "kcal/mol": (627.5094740631, 4),
    "kJ/mol": (2625.4996394799, 4),
    "eV": (27.211386245988, 6),
}
DEFAULT_UNIT = "Eh"


def _energy(name, value, unit=DEFAULT_UNIT):
    """A row for an energy given in hartree, shown in ``unit``; '-' when absent."""
    if value is None:
        return (name, "-", unit)
    factor, decimals = ENERGY_UNITS[unit]
    return (name, f"{value * factor:.{decimals}f}", unit)


def selection_rows(table):
    """The selected cells as rows of text, ordered by row then column."""
    rows = {}
    for it in table.selectedItems():
        rows.setdefault(it.row(), {})[it.column()] = it.text()
    return [[cols[c] for c in sorted(cols)] for _, cols in sorted(rows.items())]


class _ThermalTable(QTableWidget):
    """Read-only table whose Ctrl+C copies every selected cell as TSV."""

    def __init__(self, on_copy, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._on_copy = on_copy

    def keyPressEvent(self, event):
        """Copy the selection on Ctrl+C, otherwise defer to Qt."""
        if event.matches(QKeySequence.StandardKey.Copy):
            self._on_copy()
            return
        super().keyPressEvent(event)


class ThermalTableDialog(QDialog):
    """Dialog listing thermochemistry values (enthalpy, entropy, free energy)."""

    def __init__(self, parent, data):
        super().__init__(parent)
        self.setWindowTitle("Thermochemistry")
        self.resize(480, 480)
        self.data = data
        self.settings_file = os.path.join(os.path.dirname(__file__), "settings.json")

        layout = QVBoxLayout(self)

        top_row = QHBoxLayout()
        self.chk_details = QCheckBox("Show Detailed Values")
        self.chk_details.setChecked(False)
        self.chk_details.stateChanged.connect(lambda: self.update_table())
        top_row.addWidget(self.chk_details)
        top_row.addStretch()
        top_row.addWidget(QLabel("Energy unit:"))
        self.cmb_unit = QComboBox()
        self.cmb_unit.addItems(list(ENERGY_UNITS))
        self.cmb_unit.setCurrentText(DEFAULT_UNIT)
        self.cmb_unit.setToolTip("Unit for every energy row; temperature stays in K.")
        self.cmb_unit.currentTextChanged.connect(lambda: self.update_table())
        top_row.addWidget(self.cmb_unit)
        layout.addLayout(top_row)

        self.table = _ThermalTable(self.copy_selection)
        self.table.setColumnCount(len(HEADERS))
        self.table.setHorizontalHeaderLabels(HEADERS)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectItems)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.table.setAlternatingRowColors(True)
        self.table.verticalHeader().setVisible(False)
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        header.setHighlightSections(False)
        self.table.setToolTip(
            "Select cells (Shift/Ctrl+click or drag) and press Ctrl+C to copy them."
        )

        act_copy_sel = QAction("Copy Selection", self.table)
        act_copy_sel.triggered.connect(lambda: self.copy_selection())
        act_copy_all = QAction("Copy All", self.table)
        act_copy_all.triggered.connect(lambda: self.copy_all())
        self.table.addAction(act_copy_sel)
        self.table.addAction(act_copy_all)
        self.table.setContextMenuPolicy(Qt.ContextMenuPolicy.ActionsContextMenu)
        layout.addWidget(self.table)

        btn_row = QHBoxLayout()
        btn_copy = QPushButton("Copy to Clipboard")
        btn_copy.setToolTip(
            "Copies the selected cells, or the whole table when at most one "
            "cell is selected."
        )
        btn_copy.clicked.connect(lambda: self.copy_table())
        btn_row.addWidget(btn_copy)

        btn_csv = QPushButton("Export to CSV")
        btn_csv.clicked.connect(lambda: self.export_csv())
        btn_row.addWidget(btn_csv)

        btn_row.addStretch()
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)  # accept() skips closeEvent (Qt >= 6.3)
        btn_row.addWidget(btn_close)
        layout.addLayout(btn_row)

        self.update_table()
        self.load_settings()

    def reject(self):
        """Route Esc through close() so closeEvent cleanup always runs."""
        self.close()

    def closeEvent(self, event):
        """Save settings before the dialog closes."""
        self.save_settings()
        # accept() not super().closeEvent(): QDialog.closeEvent calls reject(),
        # which is routed back through close() and would recurse.
        event.accept()

    def load_settings(self):
        """Restore the saved 'show detailed values' and energy-unit preferences."""
        if os.path.exists(self.settings_file):
            settings = load_section(self.settings_file, "thermal_settings")
            try:
                if "show_details" in settings:
                    self.chk_details.setChecked(bool(settings["show_details"]))
                if settings.get("energy_unit") in ENERGY_UNITS:
                    self.cmb_unit.setCurrentText(settings["energy_unit"])

            except (
                RuntimeError,
                AttributeError,
                KeyError,
                IndexError,
                ValueError,
            ) as e:
                logging.warning("Error loading thermal settings: %s", e)

    def save_settings(self):
        """Persist the 'show detailed values' and energy-unit preferences."""
        thermal_settings = {
            "show_details": self.chk_details.isChecked(),
            "energy_unit": self.energy_unit(),
        }
        save_section(self.settings_file, "thermal_settings", thermal_settings)

    def energy_unit(self):
        """The energy unit picked in the combo box, falling back to hartree."""
        unit = self.cmb_unit.currentText()
        return unit if unit in ENERGY_UNITS else DEFAULT_UNIT

    def build_rows(self):
        """The table rows as (property, value, unit) or (title, _SECTION) tuples."""
        data = self.data
        unit = self.energy_unit()
        rows = []
        e_el = data.get("electronic_energy")
        if e_el is not None:
            rows.append(_energy("Electronic Energy (SP)", e_el, unit))
        temp = data.get("temperature")
        if temp is not None:
            rows.append(("Temperature", f"{temp:.2f}", "K"))

        rows.extend(
            _energy(name, data.get(key), unit)
            for name, key in (
                ("Zero Point Energy", "zpe"),
                ("Total Thermal Energy (U)", "thermal_energy"),
                ("Total Enthalpy (H)", "enthalpy"),
                ("Enthalpy Correction (H - E_el)", "enthalpy_corr"),
                ("Entropy Term (T*S)", "entropy"),
                ("Gibbs Free Energy (G)", "gibbs"),
                ("Gibbs Correction (G - E_el)", "gibbs_corr"),
            )
        )

        imag_count = data.get("imaginary_freq_count")
        if imag_count is not None:
            rows.append(("Imaginary Frequencies", str(imag_count), ""))

        if self.chk_details.isChecked():
            rows.append(("Detailed Corrections", _SECTION))
            rows.extend(
                _energy(name, data.get(key), unit)
                for name, key in (
                    ("Thermal Vib Correction", "corr_vib"),
                    ("Thermal Rot Correction", "corr_rot"),
                    ("Thermal Trans Correction", "corr_trans"),
                    ("Total Thermal Correction", "corr_thermal_total"),
                    ("Non-thermal (ZPE) Correction", "corr_zpe"),
                    ("Total Correction (ZPE + Thermal)", "corr_total"),
                )
            )
            rows.append(("Entropy Breakdown", _SECTION))
            rows.extend(
                _energy(name, data.get(key), unit)
                for name, key in (
                    ("Electronic Entropy (T*S_el)", "s_el"),
                    ("Vibrational Entropy (T*S_vib)", "s_vib"),
                    ("Rotational Entropy (T*S_rot)", "s_rot"),
                    ("Translational Entropy (T*S_trans)", "s_trans"),
                )
            )
        return rows

    def update_table(self):
        """Rebuild the table, including the detail rows when enabled."""
        rows = self.build_rows()
        self.table.clearSpans()
        self.table.setRowCount(len(rows))
        right = Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        center = Qt.AlignmentFlag.AlignCenter
        for r, row in enumerate(rows):
            if row[1] == _SECTION:
                item = QTableWidgetItem(row[0])
                font = item.font()
                font.setBold(True)
                item.setFont(font)
                self.table.setItem(r, 0, item)
                self.table.setItem(r, 1, QTableWidgetItem(""))
                self.table.setItem(r, 2, QTableWidgetItem(""))
                self.table.setSpan(r, 0, 1, len(HEADERS))
                continue
            name, value, unit = row
            self.table.setItem(r, 0, QTableWidgetItem(name))
            val_item = QTableWidgetItem(value)
            val_item.setTextAlignment(right)
            self.table.setItem(r, 1, val_item)
            unit_item = QTableWidgetItem(unit)
            unit_item.setTextAlignment(center)
            self.table.setItem(r, 2, unit_item)

    def table_rows(self):
        """Every row of the table as [property, value, unit] text."""
        out = []
        for r in range(self.table.rowCount()):
            cells = [self.table.item(r, c) for c in range(len(HEADERS))]
            out.append([c.text() if c is not None else "" for c in cells])
        return out

    @staticmethod
    def _to_clipboard(rows):
        QApplication.clipboard().setText("".join("\t".join(r) + "\n" for r in rows))

    def copy_all(self):
        """Put the whole table, with its header, on the clipboard as TSV."""
        self._to_clipboard([HEADERS] + self.table_rows())

    def copy_selection(self):
        """Put the selected cells on the clipboard as TSV; nothing if none."""
        rows = selection_rows(self.table)
        if rows:
            self._to_clipboard(rows)

    def copy_table(self):
        """Copy the selected cells, or the whole table for a single-cell selection."""
        if len(self.table.selectedItems()) > 1:
            self.copy_selection()
        else:
            self.copy_all()

    def export_csv(self):
        """Prompt for a path and write the property/value/unit table to CSV."""
        path, _ = QFileDialog.getSaveFileName(self, "Save CSV", "", "CSV Files (*.csv)")
        if path:
            try:
                with open(path, "w", newline="", encoding="utf-8") as f:
                    writer = csv.writer(f)
                    writer.writerow(HEADERS)
                    writer.writerows(self.table_rows())
                notify(self, f"Data exported to {path}", 5000)
            except (OSError, IndexError, ValueError) as e:
                logging.warning("Thermochemistry: exporting to %s failed: %s", path, e)
                QMessageBox.critical(
                    self, "Export Error", f"Failed to export CSV:\n{e}"
                )
