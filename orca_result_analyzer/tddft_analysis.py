from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, 
                             QRadioButton, QDoubleSpinBox, QCheckBox, QPushButton, 
                             QFileDialog, QMessageBox, QGroupBox, QComboBox)
import os
import json
import csv

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
        self.resize(600, 700)
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
        settings_vbox = QVBoxLayout(settings_group)
        
        # Row 1: Type & Gauge
        row1_layout = QHBoxLayout()
        row1_layout.addWidget(QLabel("Spectrum Type:"))
        self.radio_abs = QRadioButton("Absorption")
        self.radio_abs.setChecked(True)
        self.radio_abs.toggled.connect(self.switch_spectrum_type)
        row1_layout.addWidget(self.radio_abs)
        
        self.radio_cd = QRadioButton("CD")
        self.radio_cd.toggled.connect(self.switch_spectrum_type)
        row1_layout.addWidget(self.radio_cd)
        
        row1_layout.addWidget(QLabel(" | Gauge:"))
        self.combo_gauge = QComboBox()
        self.combo_gauge.addItems(["Length (Electric)", "Velocity"])
        self.combo_gauge.currentIndexChanged.connect(self.switch_spectrum_type)
        self.combo_gauge.setEnabled(False)
        self.combo_gauge.setFixedWidth(130)
        row1_layout.addWidget(self.combo_gauge)
        row1_layout.addStretch()
        settings_vbox.addLayout(row1_layout)
        
        # Row 2: Sigma & Sticks
        row2_layout = QHBoxLayout()
        row2_layout.addWidget(QLabel("Gaussian Sigma:"))
        self.spin_sigma = QDoubleSpinBox()
        self.spin_sigma.setRange(0.01, 5000.0)
        self.spin_sigma.setValue(5.0)
        self.spin_sigma.valueChanged.connect(self.update_spectrum_sigma)
        row2_layout.addWidget(self.spin_sigma)
        
        self.combo_sigma_unit = QComboBox()
        self.combo_sigma_unit.addItems(["nm", "cm⁻¹"])
        self.combo_sigma_unit.currentIndexChanged.connect(self.update_spectrum_sigma)
        row2_layout.addWidget(self.combo_sigma_unit)
        
        row2_layout.addSpacing(20)
        self.chk_sticks = QCheckBox("Show Transitions (Sticks)")
        self.chk_sticks.setChecked(True)
        if self.spectrum:
            self.chk_sticks.stateChanged.connect(self.spectrum.set_sticks)
        row2_layout.addWidget(self.chk_sticks)
        
        row2_layout.addStretch()
        settings_vbox.addLayout(row2_layout)
        
        # Row 3: Physical Broadening
        row3_layout = QHBoxLayout()
        self.chk_physical = QCheckBox("Physical (Area-preserving) Broadening")
        self.chk_physical.setChecked(False)
        self.chk_physical.toggled.connect(self.switch_spectrum_type)
        row3_layout.addWidget(self.chk_physical)
        row3_layout.addStretch()
        settings_vbox.addLayout(row3_layout)
        
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
        self.spin_x_min.setValue(100)
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

        self.btn_report = QPushButton("Export Full Data (.txt)")
        self.btn_report.clicked.connect(self.save_orca_report)
        action_layout.addWidget(self.btn_report)

        action_layout.addStretch()
        
        self.btn_close = QPushButton("Close")
        self.btn_close.setFixedWidth(100)
        self.btn_close.clicked.connect(self.accept)
        action_layout.addWidget(self.btn_close)
        
        main_layout.addLayout(action_layout)
        self.load_settings()
        self.switch_spectrum_type()
        
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
                
                if "sigma_unit_idx" in settings:
                    self.combo_sigma_unit.setCurrentIndex(int(settings["sigma_unit_idx"]))

                if "show_sticks" in settings:
                    self.chk_sticks.setChecked(bool(settings["show_sticks"]))
                
                if "physical" in settings:
                    self.chk_physical.setChecked(bool(settings["physical"]))
                        
            except Exception as e:
                print(f"Error loading TDDFT settings: {e}")

    def update_spectrum_sigma(self):
        if not self.spectrum: return
        val = self.spin_sigma.value()
        unit_idx = self.combo_sigma_unit.currentIndex() # 0: nm, 1: cm-1
        
        self.spectrum.sigma = val
        self.spectrum.broaden_in_energy = (unit_idx == 1)
        self.spectrum.update()

    def save_settings(self):
        all_settings = {}
        if os.path.exists(self.settings_file):
            try:
                with open(self.settings_file, 'r') as f:
                    all_settings = json.load(f)
            except: pass
            
        tddft_settings = {
            "sigma": self.spin_sigma.value(),
            "sigma_unit_idx": self.combo_sigma_unit.currentIndex(),
            "show_sticks": self.chk_sticks.isChecked(),
            "physical": self.chk_physical.isChecked()
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
                
    def save_orca_report(self):
        if not self.excitations:
            QMessageBox.warning(self, "No Data", "No excitation data to export.")
            return

        path, _ = QFileDialog.getSaveFileName(self, "Save Report", "", "Text Files (*.txt)")
        if not path:
            return

        import datetime

        try:
            with open(path, 'w', encoding='utf-8') as f:
                # --- Header ---
                f.write("=" * 80 + "\n")
                f.write("             ORCA TD-DFT / TDA EXCITATION SPECTRUM ANALYSIS\n")
                f.write("=" * 80 + "\n")
                f.write(f"Generated on: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write("-" * 80 + "\n")
                f.write("UNITS & DEFINITIONS:\n")
                f.write("  Energy      : eV, nm, cm^-1\n")
                f.write("  f (Osc)     : Oscillator Strength (Dimensionless)\n")
                f.write("  R (Rot)     : Rotatory Strength   (10^-40 esu^2 cm^2)\n")
                f.write("  Gauge       : L = Length (Dipole), V = Velocity\n")
                f.write("-" * 80 + "\n\n")

                for item in self.excitations:
                    state = item.get('state', '?')
                    
                    # Skip Ground State if labeled '0'
                    if str(state) == '0':
                        continue

                    # --- Prepare Energy Values ---
                    energy_ev = item.get('energy_ev', 0.0)
                    energy_nm = item.get('energy_nm', 0.0)
                    energy_cm = item.get('energy_cm', 0.0)
                    
                    # Auto-fill missing energy units
                    if energy_cm == 0.0 and energy_nm > 0:
                        energy_cm = 10000000.0 / energy_nm
                    if energy_ev == 0.0 and energy_nm > 0:
                        energy_ev = 1239.84193 / energy_nm

                    # --- Prepare Strength Values (Length vs Velocity) ---
                    # Oscillator Strength
                    # 変更: デフォルトを 0.0 ではなく None にして区別する
                    f_len = item.get('osc_len', item.get('osc'))
                    f_vel = item.get('osc_vel') 
                    
                    # Rotatory Strength
                    r_len = item.get('rot_len', item.get('rotatory_strength'))
                    r_vel = item.get('rot_vel')

                    # Helper for formatting (None -> "N/A")
                    def fmt_val(v):
                        return f"{v:>9.6f}" if v is not None else "      N/A"

                    # --- Writing the Block ---
                    f.write(f"STATE {state:>3}  --------------------------------------------------------\n")
                    f.write(f"  Energy : {energy_ev:>8.4f} eV  |  {energy_nm:>8.2f} nm  |  {energy_cm:>8.1f} cm^-1\n")
                    
                    f.write(f"  f (Osc): L={fmt_val(f_len)}  |  V={fmt_val(f_vel)}\n")

                    # Rotatory Strength (CD計算データがある場合のみ行自体を表示)
                    if r_len is not None or r_vel is not None:
                        f.write(f"  R (Rot): L={fmt_val(r_len)}  |  V={fmt_val(r_vel)}\n")
                    
                    f.write("\n  Major Transitions:\n")
                    
                    # Transition Details
                    transitions = item.get('transitions', [])
                    if isinstance(transitions, list):
                        if not transitions:
                            f.write("    (No transition details available)\n")
                        for trans in transitions:
                            # Clean up formatting if needed
                            f.write(f"    {trans}\n")
                    else:
                        f.write(f"    {transitions}\n")

                    f.write("\n") # Empty line between states

            if hasattr(self.parent(), 'mw') and self.parent().mw:
                self.parent().mw.statusBar().showMessage(f"Report saved to {path}", 5000)
            else:
                QMessageBox.information(self, "Exported", f"Report saved to:\n{path}")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save report:\n{e}")
    
    def switch_spectrum_type(self):
        """
        Switch between Absorption (Oscillator Strength) and CD (Rotatory Strength).
        Robustly handles missing specific gauge keys by falling back to generic keys.
        Apply physical conversion factors if requested.
        """
        if not self.spectrum: return
        
        is_cd = self.radio_cd.isChecked()
        is_physical = self.chk_physical.isChecked()
        
        self.combo_gauge.setEnabled(True)
        gauge_mode = self.combo_gauge.currentIndex() # 0: Length, 1: Velocity

        # Sample data to check keys
        sample_data = self.excitations[0] if self.excitations and len(self.excitations) > 0 else {}

        # 1. Determine stick data key and units
        if is_cd:
            target_key = 'rot_len' if gauge_mode == 0 else 'rot_vel'
            if target_key not in sample_data and 'rotatory_strength' in sample_data:
                target_key = 'rotatory_strength'
            
            y_unit_main = 'Molar Circular Dichroism (Δε) [L mol⁻¹ cm⁻¹]' if is_physical else 'CD Intensity (arb. units)'
            y_unit_sticks = 'Rotatory Strength (R) [10⁻⁴⁰ esu² cm²]'
        else:
            target_key = 'osc_len' if gauge_mode == 0 else 'osc_vel'
            if target_key not in sample_data and 'osc' in sample_data:
                target_key = 'osc'
            
            y_unit_main = 'Molar Absorbance (ε) [L mol⁻¹ cm⁻¹]' if is_physical else 'Intensity (arb. units)'
            y_unit_sticks = 'Oscillator Strength (f)'

        # 2. Prepare processed data list for the widget
        processed_data = []
        
        # Physical Constants
        # f -> epsilon: integral(eps) d_nu = 2.315e8 * f
        EPS_FACTOR = 2.315e8
        # R -> delta_eps: Delta_epsilon = (1/22.96) * nu * R * G(nu)
        CD_FACTOR = 1.0 / 22.96

        for item in self.excitations:
            new_item = item.copy()
            val = item.get(target_key, 0.0)
            
            if is_physical:
                if is_cd:
                    # R (10^-40 cgs) -> Delta_epsilon
                    # Need wavenumber (cm-1)
                    wn = item.get('energy_cm', 0.0)
                    if wn == 0 and item.get('energy_nm', 0) > 0:
                        wn = 1e7 / item.get('energy_nm')
                    
                    # Target "Area" is val * wn * CD_FACTOR
                    new_item['processed_y'] = val * wn * CD_FACTOR
                else:
                    # f -> epsilon
                    # Target "Area" is f * EPS_FACTOR
                    new_item['processed_y'] = val * EPS_FACTOR
            else:
                new_item['processed_y'] = val
            
            processed_data.append(new_item)

        # 3. Update Widget
        self.spectrum.data_list = processed_data
        self.spectrum.y_key = 'processed_y'
        self.spectrum.y_unit = y_unit_main
        self.spectrum.y_unit_sticks = y_unit_sticks
        self.spectrum.normalization_mode = 'area' if is_physical else 'height'
        
        # Ensure sigma settings are sync'd
        self.update_spectrum_sigma()
        
        self.spectrum.update()

    def on_spectrum_click(self, item):
        """Handle click on spectrum peak"""
        if item is None:
            return
        state = item.get('state', '?')
        energy_ev = item.get('energy_ev', 0.0)
        energy_nm = item.get('energy_nm', 0.0)
        
        is_cd = self.radio_cd.isChecked()
        gauge_mode = self.combo_gauge.currentIndex() # 0: Length, 1: Velocity
        
        msg = f"Transition to Excited State {state}\n"
        msg += f"Energy: {energy_ev:.4f} eV ({energy_nm:.2f} nm)\n"
        msg += "-" * 30 + "\n"

        # Show Oscillator Strengths
        osc_len = item.get('osc_len', item.get('osc', 0.0))
        osc_vel = item.get('osc_vel', 0.0)
        
        # Show Rotatory Strengths
        rot_len = item.get('rot_len', item.get('rotatory_strength', 0.0))
        rot_vel = item.get('rot_vel', 0.0)

        if not is_cd:
            msg += f"Oscillator Strength (f):\n"
            msg += f"  - Length   : {osc_len:.6f}" + (" *" if gauge_mode == 0 else "") + "\n"
            msg += f"  - Velocity : {osc_vel:.6f}" + (" *" if gauge_mode == 1 else "") + "\n"
        else:
            msg += f"Rotatory Strength (R) [10⁻⁴⁰ esu² cm²]:\n"
            msg += f"  - Length   : {rot_len:.6f}" + (" *" if gauge_mode == 0 else "") + "\n"
            msg += f"  - Velocity : {rot_vel:.6f}" + (" *" if gauge_mode == 1 else "") + "\n"
        
        msg += "-" * 30 + "\n"



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
