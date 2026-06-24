"""Bond Analysis panel: Mayer bond orders and NBO results.

Tabs are shown only for the data present in the file:
  - Mayer bond orders
  - NBO orbitals (NATURAL BOND ORBITALS summary)
  - NBO second-order perturbation (donor -> acceptor E(2))

Selecting a row highlights the relevant bond / atoms in the host 3D editor:
a Mayer bond row highlights that bond; an NBO orbital row highlights the
atom(s) it spans.
"""

import logging

from PyQt6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QLabel,
    QPushButton,
)
from PyQt6.QtCore import Qt


_TABLE_STYLE = """
    QTableWidget { gridline-color: #e6e6e6; background: #ffffff; }
    QTableWidget::item { padding: 4px 8px; }
    QTableWidget::item:selected { background: #cfe5ff; color: #000; }
    QHeaderView::section {
        background-color: #f3f3f3;
        padding: 5px 8px;
        border: none;
        border-bottom: 1px solid #cccccc;
        font-weight: bold;
    }
"""


def _make_table(headers, rows):
    table = QTableWidget(len(rows), len(headers))
    table.setHorizontalHeaderLabels(headers)
    table.verticalHeader().setVisible(False)
    table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
    table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
    table.setSelectionMode(QTableWidget.SelectionMode.SingleSelection)
    table.setAlternatingRowColors(True)
    table.setShowGrid(False)
    table.setStyleSheet(_TABLE_STYLE)
    for r, row in enumerate(rows):
        for c, value in enumerate(row):
            item = QTableWidgetItem(str(value))
            if c > 0:
                item.setTextAlignment(
                    Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
                )
            table.setItem(r, c, item)
    table.verticalHeader().setDefaultSectionSize(26)
    table.resizeColumnsToContents()
    header = table.horizontalHeader()
    header.setHighlightSections(False)
    header.setStretchLastSection(True)
    header.setMinimumSectionSize(60)
    return table


class BondAnalysisDialog(QDialog):
    def __init__(self, parent, data):
        super().__init__(parent)
        self.parent_dlg = parent  # OrcaResultAnalyzerDialog (has .mw, .parser)
        self.setWindowTitle("Bond Analysis")
        self.resize(640, 500)

        self._actors = []
        self._mbo = data.get("mayer_bond_orders", [])
        self._nbo = data.get("nbo_orbitals", [])
        pert = data.get("nbo_perturbation", [])

        layout = QVBoxLayout(self)
        tabs = QTabWidget()

        if self._mbo:
            rows = [
                [
                    f"{b['atom_idx1']} {b['atom_sym1']}",
                    f"{b['atom_idx2']} {b['atom_sym2']}",
                    f"{b['order']:.4f}",
                ]
                for b in self._mbo
            ]
            t = _make_table(["Atom 1", "Atom 2", "Bond Order"], rows)
            t.itemSelectionChanged.connect(lambda tbl=t: self._on_mayer_selected(tbl))
            tabs.addTab(t, f"Mayer Bond Orders ({len(self._mbo)})")

        if self._nbo:
            rows = [
                [
                    str(o["index"]),
                    o["type"],
                    o["atoms"],
                    f"{o['occupancy']:.5f}",
                    f"{o['energy']:.5f}",
                ]
                for o in self._nbo
            ]
            t = _make_table(["#", "Type", "Atoms", "Occupancy", "Energy (Eh)"], rows)
            t.itemSelectionChanged.connect(lambda tbl=t: self._on_nbo_selected(tbl))
            tabs.addTab(t, f"NBO Orbitals ({len(self._nbo)})")

        if pert:
            rows = [
                [
                    p["donor"],
                    p["acceptor"],
                    f"{p['e2_kcal']:.2f}",
                    f"{p['e_diff']:.2f}",
                    f"{p['fock']:.3f}",
                ]
                for p in pert
            ]
            tabs.addTab(
                _make_table(
                    ["Donor NBO", "Acceptor NBO", "E(2) kcal/mol", "ΔE a.u.", "F a.u."],
                    rows,
                ),
                f"NBO E(2) ({len(pert)})",
            )

        if tabs.count() == 0:
            layout.addWidget(QLabel("No bond-analysis data found."))
        else:
            hint = QLabel("Select a row to highlight it in the 3D view.")
            hint.setStyleSheet("color:#777; font-size:9pt;")
            layout.addWidget(hint)
            layout.addWidget(tabs)

        btn_row = QHBoxLayout()
        btn_row.addStretch()
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        btn_row.addWidget(btn_close)
        layout.addLayout(btn_row)

    # ----- 3D highlighting -------------------------------------------------

    def _coords(self):
        parser = getattr(self.parent_dlg, "parser", None)
        if parser is None:
            return []
        return parser.data.get("coords", [])

    def _plotter(self):
        mw = getattr(self.parent_dlg, "mw", None)
        if mw is None or not hasattr(mw, "plotter"):
            return None
        return mw.plotter

    def _clear_highlight(self):
        plotter = self._plotter()
        if plotter is None:
            self._actors = []
            return
        for actor in self._actors:
            try:
                plotter.remove_actor(actor)
            except Exception as _e:
                logging.warning("silenced: %s", _e)
        self._actors = []
        try:
            plotter.render()
        except Exception as _e:
            logging.warning("silenced: %s", _e)

    def _highlight_atoms(self, indices):
        self._clear_highlight()
        plotter = self._plotter()
        coords = self._coords()
        if plotter is None or not coords:
            return
        try:
            import pyvista as pv

            for idx in indices:
                if 0 <= idx < len(coords):
                    sphere = pv.Sphere(radius=0.45, center=coords[idx])
                    self._actors.append(
                        plotter.add_mesh(sphere, color="yellow", opacity=0.45)
                    )
            plotter.render()
        except Exception as _e:
            logging.warning("silenced: %s", _e)

    def _highlight_bond(self, i, j):
        self._clear_highlight()
        plotter = self._plotter()
        coords = self._coords()
        if plotter is None or not coords:
            return
        if not (0 <= i < len(coords) and 0 <= j < len(coords)):
            return
        try:
            import pyvista as pv

            p1, p2 = coords[i], coords[j]
            tube = pv.Line(p1, p2).tube(radius=0.12)
            self._actors.append(plotter.add_mesh(tube, color="yellow", opacity=0.7))
            for p in (p1, p2):
                self._actors.append(
                    plotter.add_mesh(
                        pv.Sphere(radius=0.3, center=p), color="orange", opacity=0.6
                    )
                )
            plotter.render()
        except Exception as _e:
            logging.warning("silenced: %s", _e)

    def _on_mayer_selected(self, table):
        row = table.currentRow()
        if 0 <= row < len(self._mbo):
            b = self._mbo[row]
            self._highlight_bond(b["atom_idx1"], b["atom_idx2"])

    def _on_nbo_selected(self, table):
        row = table.currentRow()
        if 0 <= row < len(self._nbo):
            self._highlight_atoms(self._nbo[row].get("atom_indices", []))

    def closeEvent(self, event):
        self._clear_highlight()
        super().closeEvent(event)
