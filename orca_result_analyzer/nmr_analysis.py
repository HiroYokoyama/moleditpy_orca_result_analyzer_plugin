
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, QComboBox, 
                             QDoubleSpinBox, QTableWidget, QTableWidgetItem, QHeaderView, 
                             QPushButton, QApplication)

class NMRDialog(QDialog):

    def __init__(self, parent, data):
        super().__init__(parent)
        self.setWindowTitle("NMR Chemical Shielding")
        self.resize(500, 600)
        self.data = data
        self.displayed_data = list(data) # copy
        
        layout = QVBoxLayout(self)
        
        # Controls
        ctrl = QHBoxLayout()
        ctrl.addWidget(QLabel("Reference Info:"))
        
        # Preset references
        self.combo_ref = QComboBox()
        self.combo_ref.addItems(["Custom", "TMS (C) ~182.4", "TMS (H) ~31.8", "Benzene (C) ~57.6", "Water (H) ~30.8"])
        self.combo_ref.currentIndexChanged.connect(self.on_preset_change)
        ctrl.addWidget(self.combo_ref)
        
        ctrl.addWidget(QLabel("Ref Shielding (ppm):"))
        self.spin_ref = QDoubleSpinBox()
        self.spin_ref.setRange(-1000, 20000)
        self.spin_ref.setValue(0.0)
        self.spin_ref.valueChanged.connect(self.recalc)
        ctrl.addWidget(self.spin_ref)
        
        layout.addLayout(ctrl)
        
        # Filter by Element
        flt_layout = QHBoxLayout()
        flt_layout.addWidget(QLabel("Filter Element:"))
        self.combo_elem = QComboBox()
        elements = sorted(list(set([d['atom_sym'] for d in self.data])))
        self.combo_elem.addItem("All")
        self.combo_elem.addItems(elements)
        self.combo_elem.currentTextChanged.connect(self.apply_filter)
        flt_layout.addWidget(self.combo_elem)
        flt_layout.addStretch()
        layout.addLayout(flt_layout)
        
        # Table
        self.table = QTableWidget()
        self.table.setColumnCount(4)
        self.table.setHorizontalHeaderLabels(["Idx", "Element", "Shielding (ppm)", "Shift (Ref-S) (ppm)"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.table)
        
        btn_copy = QPushButton("Copy Table")
        btn_copy.clicked.connect(self.copy_table)
        layout.addWidget(btn_copy)
        
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        layout.addWidget(btn_close)
        
        self.apply_filter() # Populate
        
    def on_preset_change(self, idx):
        # Rough values
        if idx == 1: self.spin_ref.setValue(182.4) # TMS C
        elif idx == 2: self.spin_ref.setValue(31.8) # TMS H
        elif idx == 3: self.spin_ref.setValue(57.6) # Benzene C (Gas)
        elif idx == 4: self.spin_ref.setValue(30.6) # Water H ~30-31? Not standard ref usually but ok.
        
    def apply_filter(self):
        elem = self.combo_elem.currentText()
        if elem == "All":
            self.displayed_data = self.data
        else:
            self.displayed_data = [d for d in self.data if d['atom_sym'] == elem]
        self.recalc()
        
    def recalc(self):
        ref = self.spin_ref.value()
        self.table.setRowCount(len(self.displayed_data))
        
        for r, item in enumerate(self.displayed_data):
            # Index
            self.table.setItem(r, 0, QTableWidgetItem(str(item.get("atom_idx", ""))))
            self.table.setItem(r, 1, QTableWidgetItem(item.get("atom_sym", "")))
            
            val = item.get("shielding", 0.0)
            self.table.setItem(r, 2, QTableWidgetItem(f"{val:.2f}"))
            
            shift = ref - val
            self.table.setItem(r, 3, QTableWidgetItem(f"{shift:.2f}"))
            
    def copy_table(self):
        text = "Idx\tElem\tShielding\tShift\n"
        for r in range(self.table.rowCount()):
             cols = []
             for c in range(self.table.columnCount()):
                 it = self.table.item(r, c)
                 cols.append(it.text() if it else "")
             text += "\t".join(cols) + "\n"
        QApplication.clipboard().setText(text)
