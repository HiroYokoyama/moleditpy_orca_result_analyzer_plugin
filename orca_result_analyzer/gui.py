
import os
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, 
                             QWidget, QGridLayout, QMessageBox, QMenuBar, QMenu, QFileDialog)
from PyQt6.QtGui import QAction

try:
    from rdkit import Chem
    from rdkit.Geometry import Point3D
    from rdkit.Chem import rdDetermineBonds
except ImportError:
    Chem = None
    Point3D = None
    rdDetermineBonds = None

# Imported Modules for Analysis
from .mo_analysis import MODialog
from .freq_analysis import FrequencyDialog
from .traj_analysis import TrajectoryResultDialog
from .force_analysis import ForceViewerDialog
from .charge_analysis import ChargeDialog
from .dipole_analysis import DipoleDialog
from .nmr_analysis import NMRDialog
from .tddft_analysis import TDDFTDialog
from .thermal_analysis import ThermalTableDialog

from . import PLUGIN_VERSION

class OrcaResultAnalyzerDialog(QDialog):
    def __init__(self, parent, parser, file_path, context=None):
        super().__init__(parent)
        self.mw = parent
        self.parser = parser
        self.file_path = file_path
        self.context = context
        
        self.setWindowTitle(f"ORCA Result Analyzer (v{PLUGIN_VERSION})")
        self.resize(450, 600)
        
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout(self)
        
        # Menu Bar (added as widget since QDialog doesn't have native menu bar)
        menu_bar = QMenuBar(self)
        menu_bar.setStyleSheet("QMenuBar { background-color: #f0f0f0; padding: 5px; font-size: 11pt; }")
        layout.addWidget(menu_bar)
        
        # File Menu
        file_menu = menu_bar.addMenu("&File")
        
        open_action = QAction("&Open File...", self)
        open_action.setShortcut("Ctrl+O")
        open_action.triggered.connect(self.open_file)
        file_menu.addAction(open_action)
        
        reload_action = QAction("&Reload File", self)
        reload_action.setShortcut("Ctrl+R")
        reload_action.triggered.connect(self.reload_file)
        file_menu.addAction(reload_action)
        
        file_menu.addSeparator()
        
        close_action = QAction("&Close", self)
        close_action.setShortcut("Ctrl+W")
        close_action.triggered.connect(self.close)
        file_menu.addAction(close_action)
        
        # Current File Display with Open Button
        file_frame = QWidget()
        file_frame.setStyleSheet("""
            QWidget {
                background-color: #e8f4f8;
                border: 2px solid #0066cc;
                border-radius: 5px;
                padding: 10px;
            }
        """)
        file_frame_layout = QHBoxLayout(file_frame)
        file_frame_layout.setContentsMargins(10, 10, 10, 10)
        
        # File info section
        file_info_layout = QVBoxLayout()
        lbl_current = QLabel("<b>Current File:</b>")
        lbl_current.setStyleSheet("font-size: 10pt; background: transparent; border: none; padding: 0;")
        file_info_layout.addWidget(lbl_current)
        
        self.lbl_file_path = QLabel(self.file_path)
        self.lbl_file_path.setStyleSheet("color: #0066cc; font-size: 9pt; background: transparent; border: none; padding: 0;")
        self.lbl_file_path.setWordWrap(True)
        file_info_layout.addWidget(self.lbl_file_path)
        
        file_frame_layout.addLayout(file_info_layout, 1)
        
        # Buttons Layout (Open + Reload)
        btns_top_layout = QVBoxLayout()
        
        # Large Open File Button
        btn_open_large = QPushButton("Open File...")
        btn_open_large.setStyleSheet("""
            QPushButton {
                background-color: #0066cc;
                color: white;
                font-size: 10pt;
                font-weight: bold;
                padding: 8px 15px;
                border-radius: 5px;
                border: none;
            }
            QPushButton:hover {
                background-color: #0052a3;
            }
            QPushButton:pressed {
                background-color: #003d7a;
            }
        """)
        btn_open_large.clicked.connect(self.open_file)
        btns_top_layout.addWidget(btn_open_large)
        
        # Reload Button
        btn_reload = QPushButton("Reload")
        btn_reload.setStyleSheet("""
            QPushButton {
                background-color: #28a745;
                color: white;
                font-size: 10pt;
                font-weight: bold;
                padding: 8px 15px;
                border-radius: 5px;
                border: none;
            }
            QPushButton:hover {
                background-color: #218838;
            }
            QPushButton:pressed {
                background-color: #1e7e34;
            }
        """)
        btn_reload.clicked.connect(self.reload_file)
        btns_top_layout.addWidget(btn_reload)
        
        file_frame_layout.addLayout(btns_top_layout)
        
        layout.addWidget(file_frame)
        
        # Grid for buttons
        grid = QGridLayout()
        grid.setSpacing(12)
        layout.addLayout(grid)
        
        # Button style for all analysis buttons
        button_style = """
            QPushButton {
                background-color: #ffffff;
                border: 2px solid #cccccc;
                border-radius: 8px;
                padding: 12px 15px;
                font-size: 10pt;
                font-weight: bold;
                text-align: left;
                min-height: 40px;
            }
            QPushButton:hover {
                background-color: #f0f8ff;
                border: 2px solid #0066cc;
            }
            QPushButton:pressed {
                background-color: #e0e8ff;
            }
            QPushButton:disabled {
                background-color: #f0f0f0;
                border: 2px solid #dddddd;
                color: #999999;
            }
        """
        
        # MO Analysis
        self.btn_mo = QPushButton("MO Coefficients & Cube Gen")
        self.btn_mo.setStyleSheet(button_style)
        self.btn_mo.clicked.connect(self.show_mo_analyzer)
        grid.addWidget(self.btn_mo, 0, 0)
        
        # Frequencies
        self.btn_freq = QPushButton("Frequency Analysis")
        self.btn_freq.setStyleSheet(button_style)
        self.btn_freq.clicked.connect(self.show_freq)
        grid.addWidget(self.btn_freq, 0, 1)
        
        # Scan / Optimization
        self.btn_traj = QPushButton("Scan / Optimization Results")
        self.btn_traj.setStyleSheet(button_style)
        self.btn_traj.clicked.connect(self.show_trajectory)
        grid.addWidget(self.btn_traj, 1, 0)
        
        # Forces
        self.btn_forces = QPushButton("Forces (Gradients)")
        self.btn_forces.setStyleSheet(button_style)
        self.btn_forces.clicked.connect(self.show_forces)
        grid.addWidget(self.btn_forces, 1, 1)
        
        # Thermochemistry
        self.btn_therm = QPushButton("Thermochemistry")
        self.btn_therm.setStyleSheet(button_style)
        self.btn_therm.clicked.connect(self.show_thermal)
        grid.addWidget(self.btn_therm, 2, 0)

        # TDDFT
        self.btn_tddft = QPushButton("TDDFT Spectra")
        self.btn_tddft.setStyleSheet(button_style)
        self.btn_tddft.clicked.connect(self.show_tddft)
        grid.addWidget(self.btn_tddft, 2, 1)

        # Dipole
        self.btn_dipole = QPushButton("Dipole Moment")
        self.btn_dipole.setStyleSheet(button_style)
        self.btn_dipole.clicked.connect(self.show_dipole)
        grid.addWidget(self.btn_dipole, 3, 0)
        
        # Charges
        self.btn_charge = QPushButton("Atomic Charges")
        self.btn_charge.setStyleSheet(button_style)
        self.btn_charge.clicked.connect(self.show_charges)
        grid.addWidget(self.btn_charge, 3, 1)
        
        # NMR
        self.btn_nmr = QPushButton("NMR Shielding")
        self.btn_nmr.setStyleSheet(button_style)
        self.btn_nmr.clicked.connect(self.show_nmr)
        grid.addWidget(self.btn_nmr, 4, 0)
        
        layout.addStretch()
        
        # Close button only (Open is now in menu)
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()
        
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        btn_layout.addWidget(btn_close)
        
        layout.addLayout(btn_layout)
        
        self.update_button_states()
        
    def open_file(self):
        # Get last directory from current file
        start_dir = os.path.dirname(self.file_path) if self.file_path else ""
        
        path, _ = QFileDialog.getOpenFileName(self, "Open ORCA Output", start_dir, "ORCA Output (*.out *.log)")
        if not path: return
        
        self.load_file(path)
    
    def load_file(self, path):
        """Load a file and update UI"""
        try:
            with open(path, 'r', encoding='utf-8', errors='replace') as f:
                content = f.read()
                
            import importlib
            from . import parser as parser_mod
            importlib.reload(parser_mod)
            from .parser import OrcaParser
            
            new_parser = OrcaParser()
            new_parser.load_from_memory(content, path)
            
            self.parser = new_parser
            self.file_path = path
            self.setWindowTitle(f"ORCA Analyzer - {os.path.basename(path)}")
            self.lbl_file_path.setText(path)
            
            # Auto-load 3D structure
            self.load_structure_3d()
            self.update_button_states()
            
            QMessageBox.information(self, "Loaded", f"Successfully loaded:\n{os.path.basename(path)}")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load file:\n{e}")

    def reload_file(self):
        if self.file_path and os.path.exists(self.file_path):
             self.load_file(self.file_path)
        else:
             QMessageBox.warning(self, "Error", "No file currently loaded or file not found.")

    def update_button_states(self):
        data = self.parser.data
        
        # Enable MO button if MO coefficients OR orbital energies exist
        # Strict check: must be non-empty dict or list
        mo_coeffs = data.get("mo_coeffs", {})
        orb_energies = data.get("orbital_energies", [])
        
        has_mo = bool(mo_coeffs) or (bool(orb_energies) and len(orb_energies) > 0)
        self.btn_mo.setEnabled(has_mo)
        tooltip = ""
        if not has_mo:
            tooltip = "No MO data found"
        elif orb_energies and not mo_coeffs:
            tooltip = "Orbital energies available (no coefficients for visualization)"
        self.btn_mo.setToolTip(tooltip)
        
        freqs = data.get("frequencies", [])
        has_freq = bool(freqs and len(freqs) > 0)
        self.btn_freq.setEnabled(has_freq)
        self.btn_freq.setToolTip("" if has_freq else "No frequency data found")
        
        has_traj = bool(data.get("scan_steps") or data.get("opt_steps")) 
        has_scan = bool(data.get("scan_steps"))
        self.btn_traj.setEnabled(has_scan)
        self.btn_traj.setToolTip("" if has_scan else "No trajectory/scan steps found")
       
        # Enable Forces button if gradients OR frequencies exist (can load Hessian)
        grads = data.get("gradients", [])
        has_gradients = bool(grads and len(grads) > 0)
        has_forces = has_gradients or has_freq
        self.btn_forces.setEnabled(has_forces)
        tooltip = ""
        if not has_forces:
             tooltip = "No gradients or frequency data found"
        elif has_freq and not has_gradients:
             tooltip = "Load Hessian file for force constants"
        elif has_gradients:
             tooltip = "View Forces (Gradients)"
        self.btn_forces.setToolTip(tooltip)
        
        has_thermal = bool(data.get("thermal") or (data.get("frequencies") and "thermo" in str(data))) 
        self.btn_therm.setEnabled(bool(data.get("thermal")))
        
        has_tddft = bool(data.get("tddft"))
        self.btn_tddft.setEnabled(has_tddft)
        
        has_dipole = bool(data.get("dipoles"))
        self.btn_dipole.setEnabled(has_dipole)
        
        has_charges = bool(data.get("charges"))
        self.btn_charge.setEnabled(has_charges)
        
        has_nmr = bool(data.get("nmr_shielding"))
        self.btn_nmr.setEnabled(has_nmr)
        
    def load_structure_3d(self):
        # Helper to ensure 3D structure is loaded
        atoms = self.parser.data.get("atoms", [])
        coords = self.parser.data.get("coords", [])
        
        if not atoms or not coords: return
        
        if not Chem: return
        
        try:
            mol = Chem.RWMol()
            conf = Chem.Conformer(len(atoms))
            
            for i, sym in enumerate(atoms):
                idx = mol.AddAtom(Chem.Atom(sym))
                x, y, z = coords[i]
                conf.SetAtomPosition(idx, Point3D(x, y, z))
                
            mol.AddConformer(conf)
            final_mol = mol.GetMol()
            
            try:
                rdDetermineBonds.DetermineBonds(final_mol)
            except: pass
            
            # Set as current molecule for export functionality
            if hasattr(self.mw, 'current_mol'):
                self.mw.current_mol = final_mol
            
            if hasattr(self.mw, 'draw_molecule_3d'):
                self.mw.draw_molecule_3d(final_mol)
                
            # Reset 3D camera to fit molecule
            if hasattr(self.mw, 'plotter') and self.mw.plotter:
                try:
                    self.mw.plotter.reset_camera()
                except: pass
                
            # Minimize 2D editor to focus on 3D view
            if hasattr(self.mw, 'splitter'):
                try:
                    # splitter has 2D editor at index 0, 3D view at index 1
                    # Completely hide 2D, show only 3D
                    total = self.mw.splitter.width()
                    self.mw.splitter.setSizes([0, total])
                except: pass
        except Exception as e:
            print(f"Error loading 3D: {e}")

    def show_mo_analyzer(self):
        mo_coeffs = self.parser.data.get("mo_coeffs")
        orb_energies = self.parser.data.get("orbital_energies")
        
        data_to_show = mo_coeffs if mo_coeffs else orb_energies
        
        if not data_to_show:
            QMessageBox.warning(self, "No Data", "No Molecular Orbital coefficients or energies found.")
            return
            
        if hasattr(self, 'mo_dlg') and self.mo_dlg is not None:
            self.mo_dlg.close()
            
        self.mo_dlg = MODialog(self, data_to_show)
        self.mo_dlg.show()
        
    def show_freq(self):
        freqs = self.parser.data.get("frequencies", [])
        if not freqs:
            QMessageBox.warning(self, "No Data", "No frequency data found.")
            return
        atoms = self.parser.data.get("atoms", [])
        coords = self.parser.data.get("coords", [])
        
        # Ensure only one instance is open
        if hasattr(self, 'freq_dlg') and self.freq_dlg is not None:
            self.freq_dlg.close()
            
        self.freq_dlg = FrequencyDialog(self.mw, freqs, atoms, coords)
        self.freq_dlg.show()
        
    def show_trajectory(self):
        data = self.parser.data.get("scan_steps", [])
        if not data:
             QMessageBox.warning(self, "No Info", "No trajectory steps (Scan/Optimization) found.")
             return
        if hasattr(self, 'traj_dlg') and self.traj_dlg is not None:
             self.traj_dlg.close()
        charge = self.parser.data.get("charge", 0)
        self.traj_dlg = TrajectoryResultDialog(self.mw, data, charge=charge, title="Trajectory / Scan Analysis")
        self.traj_dlg.show()
        
    def show_forces(self):
        grads = self.parser.data.get("gradients", [])
        has_freq = bool(self.parser.data.get("frequencies"))
        
        # If no gradients but has frequency data, allow opening for Hessian loading
        if not grads and not has_freq:
             QMessageBox.warning(self, "No Info", "No cartesian gradients or frequency data found.")
             return
        
        if hasattr(self, 'force_dlg') and self.force_dlg is not None:
             self.force_dlg.close()
        # Pass parser to enable Hessian file loading
        self.force_dlg = ForceViewerDialog(self, grads, parser=self.parser)
        self.force_dlg.show()
        
    def show_thermal(self):
        data = self.parser.data.get("thermal", {})
        if not data:
            QMessageBox.warning(self, "No Info", "No thermochemistry section found.")
            return
        dlg = ThermalTableDialog(self, data)
        dlg.exec()

    def show_tddft(self):
        excitations = self.parser.data.get("tddft", [])
        if not excitations:
            QMessageBox.warning(self, "No Analysis", "No TDDFT/TDA excitation energies found.")
            return
        dlg = TDDFTDialog(self, excitations)
        dlg.exec()

    def show_dipole(self):
        d = self.parser.data.get("dipoles")
        if not d:
            QMessageBox.warning(self, "No Info", "No dipole moment found.")
            return
        if hasattr(self, 'dipole_dlg') and self.dipole_dlg is not None:
             self.dipole_dlg.close()
        self.dipole_dlg = DipoleDialog(self, d)
        self.dipole_dlg.show()

    def show_charges(self):
        charges = self.parser.data.get("charges", {})
        if not charges:
            QMessageBox.warning(self, "No Info", "No atomic charges found.")
            return
        if hasattr(self, 'charge_dlg') and self.charge_dlg is not None:
             self.charge_dlg.close()
        self.charge_dlg = ChargeDialog(self, charges)
        self.charge_dlg.show()

    def show_nmr(self):
        data = self.parser.data.get("nmr_shielding", [])
        if not data:
            QMessageBox.warning(self, "No Info", "No NMR chemical shielding data found.")
            return
        dlg = NMRDialog(self, data)
        dlg.exec()
