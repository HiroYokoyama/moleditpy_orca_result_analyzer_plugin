"""Properties panel: scalar / global quantities from an ORCA result.

Collects values that do not warrant a dedicated analysis panel (energy,
charge/multiplicity, dispersion correction, <S**2> spin contamination, ...).
"""

from PyQt6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QPushButton,
)
from PyQt6.QtCore import Qt


_TABLE_STYLE = """
    QTableWidget { gridline-color: #e6e6e6; background: #ffffff; }
    QTableWidget::item { padding: 7px 16px; }
    QTableWidget::item:selected { background: #cfe5ff; color: #000; }
    QHeaderView::section {
        background-color: #f3f3f3;
        padding: 7px 16px;
        border: none;
        border-bottom: 1px solid #cccccc;
        font-weight: bold;
    }
"""


class PropertiesDialog(QDialog):
    def __init__(self, parent, data):
        super().__init__(parent)
        self.setWindowTitle("Properties")
        self.resize(540, 460)

        rows = []
        energy = data.get("scf_energy")
        if energy is not None:
            rows.append(("Final Single Point Energy (Eh)", f"{energy:.8f}"))
        rows.append(("Charge", str(data.get("charge", 0))))
        rows.append(("Multiplicity", str(data.get("mult", 1))))
        if data.get("version"):
            rows.append(("ORCA Version", str(data["version"])))
        rows.append(("Converged", "Yes" if data.get("converged") else "No"))

        disp = data.get("dispersion")
        if disp is not None:
            rows.append(("Dispersion correction (Eh)", f"{disp:.8f}"))

        s2 = data.get("spin_s2")
        if s2 and s2.get("actual") is not None:
            rows.append(("⟨S²⟩", f"{s2['actual']:.4f}"))
            if s2.get("ideal") is not None:
                rows.append(("Ideal S(S+1)", f"{s2['ideal']:.4f}"))
            if s2.get("contamination") is not None:
                rows.append(("Spin contamination", f"{s2['contamination']:+.4f}"))

        layout = QVBoxLayout(self)

        table = QTableWidget(len(rows), 2)
        table.setHorizontalHeaderLabels(["Property", "Value"])
        table.verticalHeader().setVisible(False)
        table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        table.setAlternatingRowColors(True)
        table.setShowGrid(False)
        table.setStyleSheet(_TABLE_STYLE)
        for r, (name, value) in enumerate(rows):
            name_item = QTableWidgetItem(name)
            name_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            table.setItem(r, 0, name_item)
            val_item = QTableWidgetItem(value)
            val_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            table.setItem(r, 1, val_item)
        table.verticalHeader().setDefaultSectionSize(32)
        header = table.horizontalHeader()
        header.setHighlightSections(False)
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        layout.addWidget(table)

        btn_row = QHBoxLayout()
        btn_row.addStretch()
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        btn_row.addWidget(btn_close)
        layout.addLayout(btn_row)
