from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, 
                             QRadioButton, QDoubleSpinBox, QCheckBox, QPushButton, 
                             QFileDialog, QMessageBox, QGroupBox)

try:
    from .spectrum_widget import SpectrumWidget
except ImportError:
    try:
        from spectrum_widget import SpectrumWidget
    except:
        SpectrumWidget = None

class TDDFTDialog(QDialog):
    def __init__(self, parent, excitations):
        super().__init__(parent)
        self.setWindowTitle("TDDFT Spectrum")
        self.resize(900, 700)
        self.excitations = excitations
        
        main_layout = QVBoxLayout(self)
        
        # 0. Spectrum Display (Existing Widget)
        if SpectrumWidget:
            # Data format for widget: List of dicts. 
            # Our excitations have 'energy_nm' and 'osc'.
            self.spectrum = SpectrumWidget(self.excitations, x_key='energy_nm', y_key='osc', x_unit='Wavelength (nm)', sigma=5.0)
            self.spectrum.show_legend = False
            main_layout.addWidget(self.spectrum)
        else:
            main_layout.addWidget(QLabel("SpectrumWidget not available."))
            self.spectrum = None
            
        # 1. Spectrum Settings Section
        settings_group = QGroupBox("Spectrum Settings")
        settings_layout = QHBoxLayout(settings_group)
        
        settings_layout.addWidget(QLabel("Spectrum Type:"))
        self.radio_abs = QRadioButton("Absorption")
        self.radio_abs.setChecked(True)
        self.radio_abs.toggled.connect(self.switch_spectrum_type)
        settings_layout.addWidget(self.radio_abs)
        
        self.radio_cd = QRadioButton("CD")
        self.radio_cd.toggled.connect(self.switch_spectrum_type)
        settings_layout.addWidget(self.radio_cd)
        
        settings_layout.addWidget(QLabel(" | Sigma (nm):"))
        self.spin_sigma = QDoubleSpinBox()
        self.spin_sigma.setRange(1.0, 100.0)
        self.spin_sigma.setValue(5.0)
        if self.spectrum:
            self.spin_sigma.valueChanged.connect(self.spectrum.set_sigma)
        settings_layout.addWidget(self.spin_sigma)
        
        self.chk_sticks = QCheckBox("Show Sticks")
        self.chk_sticks.setChecked(True)
        if self.spectrum:
            self.chk_sticks.stateChanged.connect(self.spectrum.set_sticks)
        settings_layout.addWidget(self.chk_sticks)
        
        settings_layout.addStretch()
        main_layout.addWidget(settings_group)
        
        # 2. Axis Controls Section
        axis_group = QGroupBox("Axis Scale")
        axis_layout = QVBoxLayout(axis_group)
        
        # X-Range
        x_row = QHBoxLayout()
        self.chk_auto_x = QCheckBox("Auto Wavelength Range")
        self.chk_auto_x.setChecked(True)
        self.chk_auto_x.stateChanged.connect(self.toggle_auto_x)
        x_row.addWidget(self.chk_auto_x)
        
        x_row.addWidget(QLabel("Range (nm):"))
        self.spin_x_min = QDoubleSpinBox()
        self.spin_x_min.setRange(0, 50000)
        self.spin_x_min.setValue(200)
        self.spin_x_min.setEnabled(False)
        self.spin_x_min.setFixedWidth(80)
        self.spin_x_min.valueChanged.connect(self.update_x_range)
        x_row.addWidget(self.spin_x_min)
        
        x_row.addWidget(QLabel("-"))
        self.spin_x_max = QDoubleSpinBox()
        self.spin_x_max.setRange(0, 50000)
        self.spin_x_max.setValue(800)
        self.spin_x_max.setEnabled(False)
        self.spin_x_max.setFixedWidth(80)
        self.spin_x_max.valueChanged.connect(self.update_x_range)
        x_row.addWidget(self.spin_x_max)
        x_row.addStretch()
        axis_layout.addLayout(x_row)
        
        # Y-Range
        y_row = QHBoxLayout()
        self.chk_auto_y = QCheckBox("Auto Intensity Range  ")
        self.chk_auto_y.setChecked(True)
        self.chk_auto_y.stateChanged.connect(self.toggle_auto_y)
        y_row.addWidget(self.chk_auto_y)
        
        y_row.addWidget(QLabel("Range:"))
        self.spin_y_min = QDoubleSpinBox()
        self.spin_y_min.setRange(-10000, 10000)
        self.spin_y_min.setValue(0)
        self.spin_y_min.setEnabled(False)
        self.spin_y_min.setFixedWidth(80)
        self.spin_y_min.valueChanged.connect(self.update_range)
        y_row.addWidget(self.spin_y_min)
        
        y_row.addWidget(QLabel("-"))
        self.spin_y_max = QDoubleSpinBox()
        self.spin_y_max.setRange(-10000, 10000)
        self.spin_y_max.setValue(1.0)
        self.spin_y_max.setEnabled(False)
        self.spin_y_max.setFixedWidth(80)
        self.spin_y_max.valueChanged.connect(self.update_range)
        y_row.addWidget(self.spin_y_max)
        y_row.addStretch()
        axis_layout.addLayout(y_row)
        
        main_layout.addWidget(axis_group)
        
        # 3. Actions Section
        action_layout = QHBoxLayout()
        
        self.btn_png = QPushButton("Export Image (PNG)")
        self.btn_png.clicked.connect(self.save_png)
        action_layout.addWidget(self.btn_png)
        
        self.btn_csv = QPushButton("Export Data (CSV)")
        self.btn_csv.clicked.connect(self.save_csv)
        action_layout.addWidget(self.btn_csv)
        
        self.btn_sticks = QPushButton("Export Sticks (CSV)")
        self.btn_sticks.clicked.connect(self.save_sticks)
        action_layout.addWidget(self.btn_sticks)
        
        action_layout.addStretch()
        
        self.btn_close = QPushButton("Close")
        self.btn_close.setFixedWidth(100)
        self.btn_close.clicked.connect(self.accept)
        action_layout.addWidget(self.btn_close)
        
        main_layout.addLayout(action_layout)

    def toggle_auto_y(self):
        is_auto = self.chk_auto_y.isChecked()
        self.spin_y_min.setEnabled(not is_auto)
        self.spin_y_max.setEnabled(not is_auto)
        if is_auto and self.spectrum:
            self.spectrum.set_auto_range()
        else:
            self.update_range()
            
    def update_range(self):
        if self.chk_auto_y.isChecked() or not self.spectrum: return
        ymin = self.spin_y_min.value()
        ymax = self.spin_y_max.value()
        self.spectrum.set_y_range(ymin, ymax)
        
    def toggle_auto_x(self):
        is_auto = self.chk_auto_x.isChecked()
        self.spin_x_min.setEnabled(not is_auto)
        self.spin_x_max.setEnabled(not is_auto)
        if is_auto and self.spectrum:
            self.spectrum.set_auto_x_range()
        else:
            self.update_x_range()
            
    def update_x_range(self):
        if self.chk_auto_x.isChecked() or not self.spectrum: return
        xmin = self.spin_x_min.value()
        xmax = self.spin_x_max.value()
        self.spectrum.set_x_range(xmin, xmax)
        
    def save_png(self):
        if not self.spectrum: return
        path, _ = QFileDialog.getSaveFileName(self, "Save Graph", "", "Images (*.png)")
        if path:
            self.spectrum.save_png(path)
            
    def save_csv(self):
        if not self.spectrum: return
        path, _ = QFileDialog.getSaveFileName(self, "Save Data", "", "CSV Files (*.csv)")
        if path:
            success = self.spectrum.save_csv(path)
            if success:
                QMessageBox.information(self, "Saved", f"Data saved to:\n{path}")
            else:
                QMessageBox.warning(self, "Error", "Failed to save CSV.")
                
    def save_sticks(self):
        if not self.spectrum: return
        path, _ = QFileDialog.getSaveFileName(self, "Export Sticks", "", "CSV Files (*.csv)")
        if path:
            success = self.spectrum.save_sticks_csv(path)
            if success:
                QMessageBox.information(self, "Exported", f"Stick data saved to:\n{path}")
            else:
                QMessageBox.warning(self, "Error", "Failed to export stick data.")
                
    def switch_spectrum_type(self):
        if not self.spectrum: return
        is_cd = self.radio_cd.isChecked()
        if is_cd:
            # Switch to CD (Rotatory Strength)
            self.spectrum.y_key = 'rotatory_strength'
            self.spectrum.y_unit = 'CD Intensity (10⁻⁴⁰ esu² cm²)'
        else:
            # Switch to Absorption (Oscillator Strength)
            self.spectrum.y_key = 'osc'
            self.spectrum.y_unit = 'Oscillator Strength'
        self.spectrum.update()
