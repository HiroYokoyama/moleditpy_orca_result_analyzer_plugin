
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QTableWidget, QTableWidgetItem, 
                             QHeaderView, QPushButton, QApplication, QCheckBox, QAbstractItemView)

class ThermalTableDialog(QDialog):
    def __init__(self, parent, data):
        super().__init__(parent)
        self.setWindowTitle("Thermochemistry")
        self.resize(550, 460)
        self.data = data
        
        layout = QVBoxLayout(self)
        
        # Checkbox for details
        self.chk_details = QCheckBox("Show Detailed Values")
        self.chk_details.setChecked(False) # Reduced by default
        self.chk_details.stateChanged.connect(self.update_table)
        layout.addWidget(self.chk_details)
        
        self.table = QTableWidget()
        self.table.setColumnCount(2)
        self.table.setHorizontalHeaderLabels(["Property", "Value"])
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.table)
        
        btn_copy = QPushButton("Copy to Clipboard")
        btn_copy.clicked.connect(self.copy_table)
        layout.addWidget(btn_copy)
        
        btn_csv = QPushButton("Export to CSV")
        btn_csv.clicked.connect(self.export_csv)
        layout.addWidget(btn_csv)
        
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        layout.addWidget(btn_close)
        
        self.update_table()
        
    def update_table(self):
        show_details = self.chk_details.isChecked()
        data = self.data
        
        # Extract values
        temp = data.get("temperature")
        e_el = data.get("electronic_energy")
        zpe = data.get("zpe")
        u = data.get("thermal_energy") 
        h = data.get("enthalpy")
        h_corr = data.get("enthalpy_corr")
        s = data.get("entropy") 
        g = data.get("gibbs")
        g_corr = data.get("gibbs_corr")
        
        # Detailed
        corr_vib = data.get("corr_vib")
        corr_rot = data.get("corr_rot")
        corr_trans = data.get("corr_trans")
        corr_thermal_total = data.get("corr_thermal_total")
        corr_zpe = data.get("corr_zpe")
        corr_total = data.get("corr_total")
        
        s_el = data.get("s_el")
        s_vib = data.get("s_vib")
        s_rot = data.get("s_rot")
        s_trans = data.get("s_trans")
        
        items = []
        if e_el is not None:
             items.append(("Electronic Energy (SP)", f"{e_el:.8f} Eh"))
        
        if temp is not None:
            items.append(("Temperature (K)", f"{temp:.2f} K"))
             
        # Basic Items
        items.extend([
            ("Zero Point Energy", zpe),
            ("Total Thermal Energy (U)", u),
            ("Total Enthalpy (H)", h),
            ("Enthalpy Correction (H - E_el)", h_corr),
            ("Entropy Term (T*S)", s),
            ("Gibbs Free Energy (G)", g),
            ("Gibbs Correction (G - E_el)", g_corr)
        ])
        
        # Imaginary Frequencies
        imag_count = data.get("imaginary_freq_count")
        if imag_count is not None:
             items.append(("Imaginary Frequencies", str(imag_count)))
        
        if show_details:
            # Add Separator and Detailed Items
            items.extend([
                ("", ""), # Separator
                ("--- Detailed Corrections ---", ""),
                ("Thermal Vib Correction", corr_vib),
                ("Thermal Rot Correction", corr_rot),
                ("Thermal Trans Correction", corr_trans),
                ("Total Thermal Correction", corr_thermal_total),
                ("Non-thermal (ZPE) Correction", corr_zpe),
                ("Total Correction (ZPE + Thermal)", corr_total),
                ("", ""),
                ("--- Entropy Breakdown ---", ""),
                ("Electronic Entropy (T*S_el)", s_el),
                ("Vibrational Entropy (T*S_vib)", s_vib),
                ("Rotational Entropy (T*S_rot)", s_rot),
                ("Translational Entropy (T*S_trans)", s_trans)
            ])
            
        self.table.setRowCount(len(items))
        
        for i, (name, val) in enumerate(items):
            self.table.setItem(i, 0, QTableWidgetItem(name))
            if isinstance(val, str):
                 self.table.setItem(i, 1, QTableWidgetItem(val))
            elif val is not None:
                self.table.setItem(i, 1, QTableWidgetItem(f"{val:.8f} Eh"))
            else:
                self.table.setItem(i, 1, QTableWidgetItem("-"))
        
    def copy_table(self):
        text = ""
        for r in range(self.table.rowCount()):
            p = self.table.item(r, 0).text()
            v = self.table.item(r, 1).text()
            text += f"{p}\t{v}\n"
        QApplication.clipboard().setText(text)

    def export_csv(self):
        from PyQt6.QtWidgets import QFileDialog
        import csv
        path, _ = QFileDialog.getSaveFileName(self, "Save CSV", "", "CSV Files (*.csv)")
        if path:
            try:
                with open(path, 'w', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow(["Property", "Value"])
                    for r in range(self.table.rowCount()):
                        p = self.table.item(r, 0).text()
                        v = self.table.item(r, 1).text()
                        writer.writerow([p, v])
            except Exception as e:
                pass # Simple error handling
