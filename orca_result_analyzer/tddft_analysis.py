
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, 
                             QRadioButton, QDoubleSpinBox, QCheckBox, QPushButton, 
                             QFileDialog, QMessageBox)

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
        self.resize(900, 650)
        self.excitations = excitations
        
        layout = QVBoxLayout(self)
        
        # Spectrum Widget
        if SpectrumWidget:
            # Data format for widget: List of dicts. 
            # Our excitations have 'energy_nm' and 'osc'.
            self.spectrum = SpectrumWidget(self.excitations, x_key='energy_nm', y_key='osc', x_unit='Wavelength (nm)', sigma=5.0)
            self.spectrum.show_legend = False
            layout.addWidget(self.spectrum)
        else:
            layout.addWidget(QLabel("SpectrumWidget not available."))
            self.spectrum = None
        
        # Controls organized into rows
        # Row 1: Spectrum Type and Basic Controls
        ctrl_row1 = QHBoxLayout()
        
        # Spectrum Type
        ctrl_row1.addWidget(QLabel("Type:"))
        self.radio_abs = QRadioButton("Absorption")
        self.radio_abs.setChecked(True)
        self.radio_abs.toggled.connect(self.switch_spectrum_type)
        ctrl_row1.addWidget(self.radio_abs)
        
        self.radio_cd = QRadioButton("CD")
        self.radio_cd.toggled.connect(self.switch_spectrum_type)
        ctrl_row1.addWidget(self.radio_cd)
        
        ctrl_row1.addWidget(QLabel(" | Sigma (nm):"))
        self.spin_sigma = QDoubleSpinBox()
        self.spin_sigma.setRange(1.0, 100.0)
        self.spin_sigma.setValue(5.0)
        if self.spectrum:
            self.spin_sigma.valueChanged.connect(self.spectrum.set_sigma)
        ctrl_row1.addWidget(self.spin_sigma)
        
        self.chk_sticks = QCheckBox("Sticks")
        self.chk_sticks.setChecked(True)
        if self.spectrum:
            self.chk_sticks.stateChanged.connect(self.spectrum.set_sticks)
        ctrl_row1.addWidget(self.chk_sticks)
        
        ctrl_row1.addStretch()
        layout.addLayout(ctrl_row1)
        
        # Row 2: Axis Range Controls
        ctrl_row2 = QHBoxLayout()
        
        # X-Range
        self.chk_auto_x = QCheckBox("Auto X")
        self.chk_auto_x.setChecked(True)
        self.chk_auto_x.stateChanged.connect(self.toggle_auto_x)
        ctrl_row2.addWidget(self.chk_auto_x)
        
        ctrl_row2.addWidget(QLabel("X:"))
        self.spin_x_min = QDoubleSpinBox()
        self.spin_x_min.setRange(0, 10000)
        self.spin_x_min.setValue(200)
        self.spin_x_min.setSuffix(" nm")
        self.spin_x_min.setEnabled(False)
        self.spin_x_min.setMaximumWidth(100)
        self.spin_x_min.valueChanged.connect(self.update_x_range)
        ctrl_row2.addWidget(self.spin_x_min)
        
        ctrl_row2.addWidget(QLabel("-"))
        self.spin_x_max = QDoubleSpinBox()
        self.spin_x_max.setRange(0, 10000)
        self.spin_x_max.setValue(800)
        self.spin_x_max.setSuffix(" nm")
        self.spin_x_max.setEnabled(False)
        self.spin_x_max.setMaximumWidth(100)
        self.spin_x_max.valueChanged.connect(self.update_x_range)
        ctrl_row2.addWidget(self.spin_x_max)
        
        ctrl_row2.addWidget(QLabel(" | "))
        
        # Y-Range
        self.chk_auto_y = QCheckBox("Auto Y")
        self.chk_auto_y.setChecked(True)
        self.chk_auto_y.stateChanged.connect(self.toggle_auto_y)
        ctrl_row2.addWidget(self.chk_auto_y)
        
        ctrl_row2.addWidget(QLabel("Y:"))
        self.spin_y_min = QDoubleSpinBox()
        self.spin_y_min.setRange(-1000, 10000)
        self.spin_y_min.setValue(0)
        self.spin_y_min.setEnabled(False)
        self.spin_y_min.setMaximumWidth(100)
        self.spin_y_min.valueChanged.connect(self.update_range)
        ctrl_row2.addWidget(self.spin_y_min)
        
        ctrl_row2.addWidget(QLabel("-"))
        self.spin_y_max = QDoubleSpinBox()
        self.spin_y_max.setRange(-1000, 10000)
        self.spin_y_max.setValue(1.0)
        self.spin_y_max.setEnabled(False)
        self.spin_y_max.setMaximumWidth(100)
        self.spin_y_max.valueChanged.connect(self.update_range)
        ctrl_row2.addWidget(self.spin_y_max)
        
        ctrl_row2.addStretch()
        layout.addLayout(ctrl_row2)
        
        # Row 3: Export Controls
        ctrl_row3 = QHBoxLayout()
        
        btn_png = QPushButton("Save PNG")
        btn_png.clicked.connect(self.save_png)
        ctrl_row3.addWidget(btn_png)
        
        btn_csv = QPushButton("Save CSV")
        btn_csv.clicked.connect(self.save_csv)
        ctrl_row3.addWidget(btn_csv)
        
        ctrl_row3.addStretch()
        
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        ctrl_row3.addWidget(btn_close)
        
        layout.addLayout(ctrl_row3)

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
