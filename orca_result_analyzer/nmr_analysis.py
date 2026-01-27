
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, QComboBox, 
                             QDoubleSpinBox, QTableWidget, QTableWidgetItem, QHeaderView, 
                             QPushButton, QApplication, QGroupBox)

class NMRDialog(QDialog):

    def __init__(self, parent, data):
        super().__init__(parent)
        self.setWindowTitle("NMR Chemical Shielding")
        self.resize(550, 650)
        self.data = data
        self.displayed_data = list(data) # copy
        
        main_layout = QVBoxLayout(self)
        
        # 1. Reference Shielding Section
        ref_group = QGroupBox("Chemical Shift Reference")
        ref_layout = QVBoxLayout(ref_group)
        
        ref_row = QHBoxLayout()
        ref_row.addWidget(QLabel("Preset:"))
        self.combo_ref = QComboBox()
        self.combo_ref.addItems(["Custom", "TMS (C) ~182.4", "TMS (H) ~31.8", "Benzene (C) ~57.6", "Water (H) ~30.8"])
        self.combo_ref.currentIndexChanged.connect(self.on_preset_change)
        ref_row.addWidget(self.combo_ref)
        
        ref_row.addWidget(QLabel(" | Ref Shielding (ppm):"))
        self.spin_ref = QDoubleSpinBox()
        self.spin_ref.setRange(-10000, 20000)
        self.spin_ref.setValue(0.0)
        self.spin_ref.setDecimals(2)
        self.spin_ref.valueChanged.connect(self.recalc)
        ref_row.addWidget(self.spin_ref)
        ref_layout.addLayout(ref_row)
        main_layout.addWidget(ref_group)
        
        # 2. Filtering Section
        filter_group = QGroupBox("Data Filter")
        filter_layout = QHBoxLayout(filter_group)
        filter_layout.addWidget(QLabel("Show Element:"))
        self.combo_elem = QComboBox()
        elements = sorted(list(set([d['atom_sym'] for d in self.data])))
        self.combo_elem.addItem("All")
        self.combo_elem.addItems(elements)
        self.combo_elem.currentTextChanged.connect(self.apply_filter)
        filter_layout.addWidget(self.combo_elem)
        filter_layout.addStretch()
        main_layout.addWidget(filter_group)
        
        # 3. Table Section
        table_group = QGroupBox("Computed Shifts")
        table_layout = QVBoxLayout(table_group)
        
        self.table = QTableWidget()
        self.table.setColumnCount(4)
        self.table.setHorizontalHeaderLabels(["Idx", "Element", "Shielding (ppm)", "Shift (ppm)"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.verticalHeader().setVisible(False)
        table_layout.addWidget(self.table)
        
        btn_copy = QPushButton("Copy to Clipboard")
        btn_copy.clicked.connect(self.copy_table)
        table_layout.addWidget(btn_copy)
        main_layout.addWidget(table_group)
        
        # 4. Actions
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        main_layout.addWidget(btn_close)
        
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
