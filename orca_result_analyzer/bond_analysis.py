"""Bond Analysis panel: Mayer bond orders and NBO results.

Tabs are shown only for the data present in the file:
  - Mayer bond orders
  - NBO orbitals (NATURAL BOND ORBITALS summary)
  - NBO second-order perturbation (donor -> acceptor E(2))
"""

from PyQt6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QLabel,
    QPushButton,
)
from PyQt6.QtCore import Qt


def _make_table(headers, rows):
    table = QTableWidget(len(rows), len(headers))
    table.setHorizontalHeaderLabels(headers)
    table.verticalHeader().setVisible(False)
    table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
    table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
    for r, row in enumerate(rows):
        for c, value in enumerate(row):
            item = QTableWidgetItem(str(value))
            if c > 0:
                item.setTextAlignment(
                    Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
                )
            table.setItem(r, c, item)
    table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
    for c in range(1, len(headers)):
        table.horizontalHeader().setSectionResizeMode(
            c, QHeaderView.ResizeMode.ResizeToContents
        )
    return table


class BondAnalysisDialog(QDialog):
    def __init__(self, parent, data):
        super().__init__(parent)
        self.setWindowTitle("Bond Analysis")
        self.resize(620, 480)

        layout = QVBoxLayout(self)
        tabs = QTabWidget()

        mbo = data.get("mayer_bond_orders", [])
        if mbo:
            rows = [
                [
                    f"{b['atom_idx1']} {b['atom_sym1']}",
                    f"{b['atom_idx2']} {b['atom_sym2']}",
                    f"{b['order']:.4f}",
                ]
                for b in mbo
            ]
            tabs.addTab(
                _make_table(["Atom 1", "Atom 2", "Bond Order"], rows),
                f"Mayer Bond Orders ({len(mbo)})",
            )

        nbo = data.get("nbo_orbitals", [])
        if nbo:
            rows = [
                [
                    str(o["index"]),
                    o["type"],
                    o["atoms"],
                    f"{o['occupancy']:.5f}",
                    f"{o['energy']:.5f}",
                ]
                for o in nbo
            ]
            tabs.addTab(
                _make_table(["#", "Type", "Atoms", "Occupancy", "Energy (Eh)"], rows),
                f"NBO Orbitals ({len(nbo)})",
            )

        pert = data.get("nbo_perturbation", [])
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
            layout.addWidget(tabs)

        btn_row = QHBoxLayout()
        btn_row.addStretch()
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        btn_row.addWidget(btn_close)
        layout.addLayout(btn_row)
