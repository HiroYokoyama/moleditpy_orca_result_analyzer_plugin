from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, 
                             QRadioButton, QDoubleSpinBox, QCheckBox, QPushButton, 
                             QFileDialog, QMessageBox, QGroupBox)
import os
import json

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
        self.settings_file = os.path.join(os.path.dirname(__file__), "settings.json")
        
        main_layout = QVBoxLayout(self)
        
        # 0. Spectrum Display (Existing Widget)
        if SpectrumWidget:
            # Data format for widget: List of dicts. 
            # Our excitations have 'energy_nm' and 'osc'.
            self.spectrum = SpectrumWidget(self.excitations, x_key='energy_nm', y_key='osc', x_unit='Wavelength (nm)', sigma=5.0)
            self.spectrum.show_legend = False
            main_layout.addWidget(self.spectrum)
            self.spectrum.clicked.connect(self.on_spectrum_click)
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
        self.load_settings()
        
    def closeEvent(self, event):
        self.save_settings()
        super().closeEvent(event)

    def load_settings(self):
        if os.path.exists(self.settings_file):
            try:
                with open(self.settings_file, 'r') as f:
                    all_settings = json.load(f)
                
                settings = all_settings.get("tddft_settings", {})
                
                if "sigma" in settings:
                    self.spin_sigma.setValue(float(settings["sigma"]))
                
                if "show_sticks" in settings:
                    self.chk_sticks.setChecked(bool(settings["show_sticks"]))
                    
                if "type" in settings:
                    if settings["type"] == "cd":
                        self.radio_cd.setChecked(True)
                    else:
                        self.radio_abs.setChecked(True)
                        
            except Exception as e:
                print(f"Error loading TDDFT settings: {e}")

    def save_settings(self):
        all_settings = {}
        if os.path.exists(self.settings_file):
            try:
                with open(self.settings_file, 'r') as f:
                    all_settings = json.load(f)
            except: pass
            
        tddft_settings = {
            "sigma": self.spin_sigma.value(),
            "show_sticks": self.chk_sticks.isChecked(),
            "type": "cd" if self.radio_cd.isChecked() else "abs"
        }
        
        all_settings["tddft_settings"] = tddft_settings
        
        try:
            with open(self.settings_file, 'w') as f:
                json.dump(all_settings, f, indent=2)
        except Exception as e:
            print(f"Error saving TDDFT settings: {e}")

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
                # print(f"Data saved to {path}")
                # QMessageBox.information(self, "Saved", f"Data saved to:\n{path}")
                if hasattr(self.parent(), 'mw') and self.parent().mw:
                    self.parent().mw.statusBar().showMessage(f"Data saved to {path}", 5000)
            else:
                QMessageBox.warning(self, "Error", "Failed to save CSV.")
                
    def save_sticks(self):
        if not self.spectrum: return
        path, _ = QFileDialog.getSaveFileName(self, "Export Sticks", "", "CSV Files (*.csv)")
        if path:
            success = self.spectrum.save_sticks_csv(path)
            if success:
                # print(f"Stick data saved to {path}")
                # QMessageBox.information(self, "Exported", f"Stick data saved to:\n{path}")
                if hasattr(self.parent(), 'mw') and self.parent().mw:
                    self.parent().mw.statusBar().showMessage(f"Stick data saved to {path}", 5000)
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

    def on_spectrum_click(self, item):
        """Handle click on spectrum peak"""
        state = item.get('state', '?')
        energy_ev = item.get('energy_ev', 0.0)
        energy_nm = item.get('energy_nm', 0.0)
        osc = item.get('osc', 0.0)
        
        msg = f"Transition to Excited State {state}\n"
        msg += f"Energy: {energy_ev:.4f} eV ({energy_nm:.2f} nm)\n"
        msg += f"Oscillator Strength: {osc:.6f}\n"
        
        if 'rotatory_strength' in item:
             msg += f"Rotatory Strength: {item['rotatory_strength']:.6f}\n"



        # Use instantiated QMessageBox for selectable text
        from PyQt6.QtCore import Qt
        msg_box = QMessageBox(self)
        msg_box.setWindowTitle(f"Transition Details - State {state}")
        
        # Add scroll area if message is too long? QMessageBox doesn't support scrolling easily.
        # But for 20 lines it's fine.
        # Let's increase the limit or make it nicer.
        
        if 'transitions' in item:
            trans_data = item['transitions']
            if trans_data and len(trans_data) > 0:
                msg += "\nOrbital Contributions:\n"
                # Show ALL transitions
                msg += "\n".join(trans_data)
            else:
                msg += f"\n(No detailed orbital contributions found - data type: {type(trans_data)}, len: {len(trans_data) if trans_data else 0})"
        else:
            msg += "\n(No 'transitions' key in item)"

        msg_box.setText(msg)
        msg_box.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        msg_box.exec()
