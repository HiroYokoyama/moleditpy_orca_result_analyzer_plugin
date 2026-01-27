
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QTableWidget, QTableWidgetItem, 
                             QHeaderView, QPushButton, QApplication)

class ThermalTableDialog(QDialog):
    def __init__(self, parent, data):
        super().__init__(parent)
        self.setWindowTitle("Thermochemistry")
        self.resize(400, 300)
        
        layout = QVBoxLayout(self)
        
        self.table = QTableWidget()
        self.table.setColumnCount(2)
        self.table.setHorizontalHeaderLabels(["Property", "Value (Eh)"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        
        items = [
            ("Zero Point Energy", data.get("zpe")),
            ("Total Thermal Energy", data.get("thermal_energy")),
            ("Enthalpy (H)", data.get("enthalpy")),
            ("Entropy (S)", data.get("entropy")),
            ("Gibbs Free Energy (G)", data.get("gibbs"))
        ]
        
        self.table.setRowCount(len(items))
        
        for i, (name, val) in enumerate(items):
            self.table.setItem(i, 0, QTableWidgetItem(name))
            if val is not None:
                self.table.setItem(i, 1, QTableWidgetItem(f"{val:.6f}"))
            else:
                self.table.setItem(i, 1, QTableWidgetItem("-"))
                
        layout.addWidget(self.table)
        
        btn_copy = QPushButton("Copy to Clipboard")
        btn_copy.clicked.connect(self.copy_table)
        layout.addWidget(btn_copy)
        
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        layout.addWidget(btn_close)
        
    def copy_table(self):
        text = ""
        for r in range(self.table.rowCount()):
            p = self.table.item(r, 0).text()
            v = self.table.item(r, 1).text()
            text += f"{p}\t{v}\n"
        QApplication.clipboard().setText(text)
