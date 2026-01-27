import os
import matplotlib
matplotlib.use('QtAgg')
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, 
                             QWidget, QGridLayout, QMessageBox, QScrollArea, QHeaderView, QDoubleSpinBox, 
                             QCheckBox, QComboBox, QAbstractItemView, QFileDialog, QFormLayout, QDialogButtonBox, QSpinBox,
                             QTableWidget, QTableWidgetItem, QSlider, QRadioButton, QMenuBar, QMenu)
from PyQt6.QtGui import QColor, QPalette, QPainter, QPen, QAction
from PyQt6.QtCore import Qt, QTimer, QPointF, QSettings
try:
    from rdkit import Chem
    from rdkit.Geometry import Point3D
except ImportError:
    Chem = None
    Point3D = None

try:
    from PIL import Image
    HAS_PIL = True
except ImportError:
    HAS_PIL = False

class OrcaResultAnalyzerDialog(QDialog):
    def __init__(self, parent, parser, file_path, context=None):
        super().__init__(parent)
        self.mw = parent
        self.parser = parser
        self.file_path = file_path
        self.context = context
        
        self.setWindowTitle(f"ORCA Analyzer - {self.parser.filename}")
        self.resize(800, 600)
        
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
        
        # Large Open File Button
        btn_open_large = QPushButton("Open File...")
        btn_open_large.setStyleSheet("""
            QPushButton {
                background-color: #0066cc;
                color: white;
                font-size: 11pt;
                font-weight: bold;
                padding: 12px 20px;
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
        file_frame_layout.addWidget(btn_open_large)
        
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
        """
        
        # MO Analysis
        btn_mo = QPushButton("MO Coefficients & Cube Gen")
        btn_mo.setStyleSheet(button_style)
        btn_mo.clicked.connect(self.show_mo_analyzer)
        grid.addWidget(btn_mo, 0, 0)
        
        # Frequencies
        btn_freq = QPushButton("Frequency Analysis")
        btn_freq.setStyleSheet(button_style)
        btn_freq.clicked.connect(self.show_freq)
        grid.addWidget(btn_freq, 0, 1)
        
        # Scan / Optimization
        btn_traj = QPushButton("Scan / Optimization Results")
        btn_traj.setStyleSheet(button_style)
        btn_traj.clicked.connect(self.show_trajectory)
        grid.addWidget(btn_traj, 1, 0)
        
        # Forces
        btn_forces = QPushButton("Forces (Gradients)")
        btn_forces.setStyleSheet(button_style)
        btn_forces.clicked.connect(self.show_forces)
        grid.addWidget(btn_forces, 1, 1)
        
        # Thermochemistry
        btn_therm = QPushButton("Thermochemistry")
        btn_therm.setStyleSheet(button_style)
        btn_therm.clicked.connect(self.show_thermal)
        grid.addWidget(btn_therm, 2, 0)

        # TDDFT
        btn_tddft = QPushButton("TDDFT Spectra")
        btn_tddft.setStyleSheet(button_style)
        btn_tddft.clicked.connect(self.show_tddft)
        grid.addWidget(btn_tddft, 2, 1)

        # Dipole
        btn_dipole = QPushButton("Dipole Moment")
        btn_dipole.setStyleSheet(button_style)
        btn_dipole.clicked.connect(self.show_dipole)
        grid.addWidget(btn_dipole, 3, 0)
        
        # Charges
        btn_charge = QPushButton("Atomic Charges")
        btn_charge.setStyleSheet(button_style)
        btn_charge.clicked.connect(self.show_charges)
        grid.addWidget(btn_charge, 3, 1)
        
        # NMR
        btn_nmr = QPushButton("NMR Shielding")
        btn_nmr.setStyleSheet(button_style)
        btn_nmr.clicked.connect(self.show_nmr)
        grid.addWidget(btn_nmr, 4, 0)
        
        layout.addStretch()
        
        # Close button only (Open is now in menu)
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()
        
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        btn_layout.addWidget(btn_close)
        
        layout.addLayout(btn_layout)
        
    def open_file(self):
        from PyQt6.QtWidgets import QFileDialog
        
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
            
            QMessageBox.information(self, "Loaded", f"Successfully loaded:\n{os.path.basename(path)}")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load file:\n{e}")
    
    
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
                from rdkit.Chem import rdDetermineBonds
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
        if not self.parser.data.get("mo_coeffs"):
            QMessageBox.warning(self, "No Data", "No Molecular Orbital coefficients found.")
            return
        dlg = MODialog(self, self.parser.data["mo_coeffs"])
        dlg.exec()
        
    def show_freq(self):
        freqs = self.parser.data.get("frequencies", [])
        if not freqs:
            QMessageBox.warning(self, "No Data", "No frequency data found.")
            return
        atoms = self.parser.data.get("atoms", [])
        coords = self.parser.data.get("coords", [])
        
        from .freq_analysis import FrequencyDialog
        dlg = FrequencyDialog(self.mw, freqs, atoms, coords)
        dlg.exec()
        
    def show_trajectory(self):
        data = self.parser.data.get("scan_steps", [])
        if not data:
             QMessageBox.warning(self, "No Info", "No trajectory steps (Scan/Optimization) found.")
             return
        if hasattr(self, 'traj_dlg') and self.traj_dlg is not None:
             self.traj_dlg.close()
        self.traj_dlg = TrajectoryResultDialog(self.mw, data, title="Trajectory / Scan Analysis")
        self.traj_dlg.show()
        
    def show_forces(self):
        grads = self.parser.data.get("gradients", [])
        if not grads:
             QMessageBox.warning(self, "No Info", "No cartesian gradients found.")
             return
        if hasattr(self, 'force_dlg') and self.force_dlg is not None:
             self.force_dlg.close()
        self.force_dlg = ForceViewerDialog(self, grads)
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
        d = self.parser.data.get("dipole")
        if not d:
            QMessageBox.warning(self, "No Info", "No dipole moment found.")
            return
        QMessageBox.information(self, "Dipole Moment", 
                                f"Total Dipole Moment:\n"
                                f"X: {d['vector'][0]:.4f}\n"
                                f"Y: {d['vector'][1]:.4f}\n"
                                f"Z: {d['vector'][2]:.4f}\n\n"
                                f"Magnitude: {d['magnitude']:.4f} Debye")

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

class ChargeDialog(QDialog):
    def __init__(self, parent, all_charges):
        super().__init__(parent)
        self.parent_dlg = parent # OrcaResultAnalyzerDialog
        self.setWindowTitle("Atomic Charges")
        self.resize(500, 600)
        self.all_charges = all_charges 
        self.current_type = next(iter(self.all_charges)) if self.all_charges else ""
        
        # Color Schemes - default schemes
        self.schemes = {
            "Red(-) - White - Blue(+)": ["red", "white", "blue"],
            "Blue(-) - White - Red(+)": ["blue", "white", "red"],
            "Red(-) - Blue(+)": ["red", "blue"],
            "Green(-) - White - Purple(+)": ["green", "white", "purple"]
        }
        
        # Load custom schemes from settings.json
        settings_file = os.path.join(os.path.dirname(__file__), "settings.json")
        self.current_scheme = "Red(-) - White - Blue(+)"
        
        if os.path.exists(settings_file):
            try:
                import json
                with open(settings_file, 'r') as f:
                    settings_data = json.load(f)
                
                # Load custom schemes
                if "custom_color_schemes" in settings_data:
                    for scheme_data in settings_data["custom_color_schemes"]:
                        name = scheme_data.get("name", "")
                        colors = scheme_data.get("colors", [])
                        if name and colors:
                            self.schemes[f"Custom: {name}"] = colors
                
                # Load last used scheme
                if "last_charge_scheme" in settings_data:
                    self.current_scheme = settings_data["last_charge_scheme"]
            except Exception as e:
                print(f"Error loading settings: {e}")
        
        layout = QVBoxLayout(self)
        
        # Type Selector
        head_layout = QHBoxLayout()
        head_layout.addWidget(QLabel("Charge Type:"))
        self.combo_type = QComboBox()
        self.combo_type.addItems(list(self.all_charges.keys()))
        self.combo_type.currentTextChanged.connect(self.on_type_change)
        head_layout.addWidget(self.combo_type)
        layout.addLayout(head_layout)
        
        # Table
        from PyQt6.QtWidgets import QTableWidget, QTableWidgetItem
        self.table = QTableWidget()
        self.table.setColumnCount(3)
        self.table.setHorizontalHeaderLabels(["Idx", "Atom", "Charge"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.table)
        
        # Coloring Controls
        grp_color = QWidget()
        grp_layout = QVBoxLayout(grp_color)
        
        # Scheme Selector
        scheme_layout = QHBoxLayout()
        scheme_layout.addWidget(QLabel("Color Scheme:"))
        self.combo_scheme = QComboBox()
        self.combo_scheme.addItems(list(self.schemes.keys()))
        self.combo_scheme.setCurrentText(self.current_scheme)
        self.combo_scheme.currentTextChanged.connect(self.on_scheme_change)
        scheme_layout.addWidget(self.combo_scheme)
        
        # Custom scheme button
        btn_custom = QPushButton("+ Custom")
        btn_custom.clicked.connect(self.edit_custom_scheme)
        scheme_layout.addWidget(btn_custom)
        
        grp_layout.addLayout(scheme_layout)
        
        # Gradient Bar
        self.grad_bar = GradientBar(self, self.schemes[self.current_scheme])
        grp_layout.addWidget(self.grad_bar)
        
        # Labels for gradient
        lbl_layout = QHBoxLayout()
        lbl_layout.addWidget(QLabel("Negative"))
        lbl_mid = QLabel("Neutral")
        lbl_mid.setAlignment(Qt.AlignmentFlag.AlignCenter)
        lbl_layout.addWidget(lbl_mid)
        lbl_max = QLabel("Positive")
        lbl_max.setAlignment(Qt.AlignmentFlag.AlignRight)
        lbl_layout.addWidget(lbl_max)
        grp_layout.addLayout(lbl_layout)
        
        # Show labels checkbox
        self.chk_show_labels = QCheckBox("Show charge labels in 3D")
        self.chk_show_labels.setChecked(False)
        self.chk_show_labels.stateChanged.connect(self.toggle_labels)
        grp_layout.addWidget(self.chk_show_labels)
        
        btn_colorize = QPushButton("Colorize Atoms in 3D View")
        btn_colorize.setStyleSheet("font-weight: bold; background-color: #2196F3; color: white;")
        btn_colorize.clicked.connect(self.apply_colors)
        grp_layout.addWidget(btn_colorize)
        
        # Reset color button
        btn_reset = QPushButton("Reset Colors")
        btn_reset.setStyleSheet("background-color: #f44336; color: white;")
        btn_reset.clicked.connect(self.reset_colors)
        grp_layout.addWidget(btn_reset)
        
        layout.addWidget(grp_color)
        
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        layout.addWidget(btn_close)
        
        self.update_table()
        
    def on_type_change(self, text):
        self.current_type = text
        self.update_table()
    
    def on_scheme_change(self, text):
        self.current_scheme = text
        colors = self.schemes.get(text, ["red", "white", "blue"])
        self.grad_bar.set_colors(colors)
        
        # Save last used scheme to settings.json
        self.save_settings()
    
    def edit_custom_scheme(self):
        """Create a custom color scheme"""
        from PyQt6.QtWidgets import QInputDialog, QColorDialog
        
        # Get scheme name
        name, ok = QInputDialog.getText(self, "Custom Scheme", "Enter scheme name:")
        if not ok or not name:
            return
        
        # Get number of colors
        num_colors, ok = QInputDialog.getInt(self, "Custom Scheme", "Number of colors (2-5):", 3, 2, 5)
        if not ok:
            return
        
        # Collect colors
        colors = []
        for i in range(num_colors):
            label = ["Negative", "Mid-Negative", "Neutral", "Mid-Positive", "Positive"][i] if num_colors <= 5 else f"Color {i+1}"
            color = QColorDialog.getColor(QColor("white"), self, f"Select {label} color")
            if not color.isValid():
                return
            colors.append(color.name())
        
        # Add to schemes
        scheme_name = f"Custom: {name}"
        self.schemes[scheme_name] = colors
        
        # Update combo box
        self.combo_scheme.addItem(scheme_name)
        self.combo_scheme.setCurrentText(scheme_name)
        
        # Save to settings
        self.save_custom_schemes()
    
    def save_custom_schemes(self):
        """Save custom schemes to settings.json"""
        self.save_settings()
    
    def save_settings(self):
        """Save all settings to settings.json"""
        import json
        settings_file = os.path.join(os.path.dirname(__file__), "settings.json")
        
        # Load existing settings or create new
        settings_data = {}
        if os.path.exists(settings_file):
            try:
                with open(settings_file, 'r') as f:
                    settings_data = json.load(f)
            except:
                pass
        
        # Get custom schemes only
        custom_schemes = []
        for name, colors in self.schemes.items():
            if name.startswith("Custom: "):
                custom_schemes.append({
                    "name": name.replace("Custom: ", ""),
                    "colors": colors
                })
        
        # Update settings
        settings_data["custom_color_schemes"] = custom_schemes
        settings_data["last_charge_scheme"] = self.current_scheme
        
        # Save to file
        try:
            with open(settings_file, 'w') as f:
                json.dump(settings_data, f, indent=2)
        except Exception as e:
            print(f"Error saving settings: {e}")
    
    def toggle_labels(self):
        """Toggle charge value labels in 3D view"""
        show = self.chk_show_labels.isChecked()
        
        # Remove old labels if they exist
        if hasattr(self, '_charge_labels'):
            for actor in self._charge_labels:
                try:
                    self.parent_dlg.mw.plotter.remove_actor(actor)
                except: pass
            self._charge_labels = []
        
        if not show:
            if hasattr(self.parent_dlg.mw, 'plotter'):
                self.parent_dlg.mw.plotter.render()
            return
        
        # Add labels
        data = self.all_charges.get(self.current_type, [])
        coords = self.parent_dlg.parser.data.get("coords", [])
        
        if not coords or len(coords) != len(data):
            QMessageBox.warning(self, "Error", "No coordinates available for labels.")
            self.chk_show_labels.setChecked(False)
            return
        
        self._charge_labels = []
        try:
            for item in data:
                idx = item['atom_idx']
                pos = coords[idx]
                label_pos = [pos[0], pos[1], pos[2] + 0.3]
                
                actor = self.parent_dlg.mw.plotter.add_point_labels(
                    [label_pos],
                    [f"{item['charge']:.2f}"],
                    font_size=10,
                    text_color='yellow',
                    point_size=0,
                    always_visible=True
                )
                self._charge_labels.append(actor)
            
            self.parent_dlg.mw.plotter.render()
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Could not add labels: {e}")
            self.chk_show_labels.setChecked(False)
    
    def reset_colors(self):
        """Reset atom colors to default CPK colors"""
        context = getattr(self.parent_dlg, "context", None)
        if not context:
            QMessageBox.warning(self, "Error", "Plugin Context not found.")
            return
        
        controller = context.get_3d_controller()
        if not controller:
            QMessageBox.warning(self, "Error", "3D Controller not available.")
            return
        
        try:
            # Reset all atom colors to default CPK
            data = self.all_charges.get(self.current_type, [])
            for item in data:
                controller.reset_atom_color(item['atom_idx'])
            
            # Remove scalar bar if exists
            if hasattr(self, '_charge_scalar_bar'):
                try:
                    self.parent_dlg.mw.plotter.remove_actor(self._charge_scalar_bar)
                    delattr(self, '_charge_scalar_bar')
                except: pass
            
            # Remove labels if exist
            if hasattr(self, '_charge_labels'):
                for actor in self._charge_labels:
                    try:
                        self.parent_dlg.mw.plotter.remove_actor(actor)
                    except: pass
                self._charge_labels = []
                self.chk_show_labels.setChecked(False)
            
            # Redraw molecule with default CPK colors
            if hasattr(self.parent_dlg.mw, 'current_mol') and self.parent_dlg.mw.current_mol:
                if hasattr(self.parent_dlg.mw, 'draw_molecule_3d'):
                    self.parent_dlg.mw.draw_molecule_3d(self.parent_dlg.mw.current_mol)
            
            # Render
            if hasattr(self.parent_dlg.mw, 'plotter'):
                self.parent_dlg.mw.plotter.render()
            
            QMessageBox.information(self, "Done", "Colors reset to CPK default.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to reset colors:\n{e}")
        
    def update_table(self):
        data = self.all_charges.get(self.current_type, [])
        self.table.setRowCount(len(data))
        for r, item in enumerate(data):
            self.table.setItem(r, 0, QTableWidgetItem(str(item['atom_idx'])))
            self.table.setItem(r, 1, QTableWidgetItem(item['atom_sym']))
            self.table.setItem(r, 2, QTableWidgetItem(f"{item['charge']:.4f}"))
            
    def apply_colors(self):
        data = self.all_charges.get(self.current_type, [])
        if not data: return
        
        charges = [d['charge'] for d in data]
        if not charges: return
        
        max_c = max(abs(min(charges)), abs(max(charges)))
        if max_c == 0: max_c = 1.0
        
        context = getattr(self.parent_dlg, "context", None)
        if not context:
            QMessageBox.warning(self, "Error", "Plugin Context not found. Cannot access 3D controller.")
            return

        controller = context.get_3d_controller()
        if not controller:
            QMessageBox.warning(self, "Error", "3D Controller not available.")
            return
        
        # Get current color scheme
        colors = self.schemes.get(self.current_scheme, ["red", "white", "blue"])
        
        # Helper to interpolate color between multiple colors
        def get_color(q):
            # Normalize charge to 0.0 (negative) to 1.0 (positive)
            norm = (q + max_c) / (2 * max_c)
            norm = max(0.0, min(1.0, norm))
            
            n = len(colors)
            if n < 2: 
                return QColor(colors[0]).name()
            
            # Find which segment we're in
            seg_len = 1.0 / (n - 1)
            idx = int(norm / seg_len)
            if idx >= n - 1: 
                return QColor(colors[-1]).name()
            
            # Interpolate within segment
            local_t = (norm - (idx * seg_len)) / seg_len
            
            c1 = QColor(colors[idx])
            c2 = QColor(colors[idx + 1])
            
            r = int(c1.red() + (c2.red() - c1.red()) * local_t)
            g = int(c1.green() + (c2.green() - c1.green()) * local_t)
            b = int(c1.blue() + (c2.blue() - c1.blue()) * local_t)
            
            return f"#{r:02x}{g:02x}{b:02x}"

        try:
            for item in data:
                idx = item['atom_idx']
                q = item['charge']
                color = get_color(q)
                controller.set_atom_color(idx, color)
            
            # Redraw molecule once after all colors are set
            if hasattr(self.parent_dlg.mw, 'current_mol') and self.parent_dlg.mw.current_mol:
                if hasattr(self.parent_dlg.mw, 'draw_molecule_3d'):
                    self.parent_dlg.mw.draw_molecule_3d(self.parent_dlg.mw.current_mol)
                
            # Add charge labels if checkbox is enabled
            if self.chk_show_labels.isChecked():
                coords = self.parent_dlg.parser.data.get("coords", [])
                if coords and len(coords) == len(data):
                    # Remove old labels if exist
                    if hasattr(self, '_charge_labels'):
                        for actor in self._charge_labels:
                            try:
                                self.parent_dlg.mw.plotter.remove_actor(actor)
                            except: pass
                    
                    self._charge_labels = []
                    for i, item in enumerate(data):
                        pos = coords[item['atom_idx']]
                        label_pos = [pos[0], pos[1], pos[2] + 0.3]
                        actor = self.parent_dlg.mw.plotter.add_point_labels(
                            [label_pos],
                            [f"{item['charge']:.2f}"],
                            font_size=10,
                            text_color='yellow',
                            point_size=0,
                            always_visible=True
                        )
                        self._charge_labels.append(actor)
            
            # Add scalar bar legend showing charge scale with color gradient
            try:
                import pyvista as pv
                import numpy as np
                
                # Create a mesh with charge values to display scalar bar
                vmin = min(charges)
                vmax = max(charges)
                
                # Convert scheme colors to colormap
                import matplotlib.colors as mcolors
                cmap_colors = [mcolors.to_rgb(QColor(c).name()) for c in colors]
                from matplotlib.colors import LinearSegmentedColormap
                cmap = LinearSegmentedColormap.from_list("charge_cmap", cmap_colors, N=256)
                
                # Remove old scalar bar if exists
                if hasattr(self, '_charge_scalar_bar'):
                    try:
                        self.parent_dlg.mw.plotter.remove_actor(self._charge_scalar_bar)
                    except: pass
                
                # Create dummy mesh for scalar bar
                dummy = pv.Box()
                dummy.point_data['charges'] = np.linspace(vmin, vmax, dummy.n_points)
                
                # Add mesh with scalar bar (invisible mesh, visible bar)
                self._charge_scalar_bar = self.parent_dlg.mw.plotter.add_mesh(
                    dummy,
                    scalars='charges',
                    cmap=cmap,
                    clim=[vmin, vmax],
                    opacity=0.0,  # Invisible mesh
                    show_scalar_bar=True,
                    scalar_bar_args={
                        'title': f'{self.current_type}',
                        'title_font_size': 14,
                        'label_font_size': 12,
                        'n_labels': 5,
                        'vertical': True,
                        'height': 0.3,
                        'width': 0.08,
                        'position_x': 0.88,
                        'position_y': 0.35,
                        'color': 'white'
                    }
                )
            except Exception as e:
                print(f"Error adding scalar bar: {e}")
            
            # Trigger update
            if hasattr(self.parent_dlg.mw, 'plotter'):
                self.parent_dlg.mw.plotter.render()
                
            QMessageBox.information(self, "Done", f"Applied '{self.current_scheme}' coloring to 3D view.")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to color atoms:\n{e}")

        

class MODialog(QDialog):
    def __init__(self, parent, mos):
        super().__init__(parent)
        self.setWindowTitle("MO Energy Levels")
        self.resize(600, 500)
        self.mos = mos
        
        layout = QVBoxLayout(self)
        
        # Simple List for now, maybe diagram later
        from PyQt6.QtWidgets import QTreeWidget, QTreeWidgetItem
        
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["No.", "Occ", "Energy (eV)", "Energy (Eh)"])
        self.tree.setColumnWidth(0, 60)
        self.tree.setColumnWidth(1, 60)
        self.tree.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection) # Allow multiple
        
        for mo in self.mos:
            item = QTreeWidgetItem([str(mo['id']), f"{mo['occ']:.2f}", f"{mo['energy_ev']:.3f}", f"{mo['energy_eh']:.5f}"])
            
            # Color code by occupancy
            if mo['occ'] > 1.0:
                item.setForeground(0, QColor("blue")) # Occupied
                item.setForeground(1, QColor("blue"))
            elif mo['occ'] > 0.0:
                item.setForeground(0, QColor("green")) # SOMO
            else:
                item.setForeground(0, QColor("gray")) # Virtual
                
            self.tree.addTopLevelItem(item)
            
        layout.addWidget(self.tree)
        
        # Buttons
        btn_layout = QHBoxLayout()
        
        self.btn_coeffs = QPushButton("View Composition")
        self.btn_coeffs.clicked.connect(self.view_coeffs)
        btn_layout.addWidget(self.btn_coeffs)
        
        # Native 3D generator
        self.btn_gen_3d = QPushButton("Generate Cube(s)")
        self.btn_gen_3d.clicked.connect(self.generate_cubes)
        btn_layout.addWidget(self.btn_gen_3d)

        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        btn_layout.addWidget(btn_close)
        
        layout.addLayout(btn_layout)

    def generate_cubes(self):
        # 1. Get Selection
        items = self.tree.selectedItems()
        if not items:
            QMessageBox.warning(self, "Selection", "Please select at least one MO.")
            return
            
        # 2. Check Data Availability
        parser = self.parent().parser
        shells = parser.data.get("basis_set_shells", [])
        if not shells:
            QMessageBox.warning(self, "Error", "No Basis Set definitions found.\nMake sure output contains 'BASIS SET IN INPUT FORMAT'.\n(Add ! PrintBasis to input)")
            return
            
        coeffs_map = parser.data.get("mo_coeffs", {})
        if not coeffs_map:
             QMessageBox.warning(self, "Error", "No MO coefficients parsed.")
             return
             
        # 3. Directory Selection
        # If multiple, ask for directory. If single, ask for file?
        # Let's standardize: Ask for directory to save MO_X.cube, MO_Y.cube etc.
        out_dir = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if not out_dir: return
        
        # 4. Prepare Batch
        from .mo_engine import BasisSetEngine, CalcWorker
        import os
        import numpy as np
        
        # Init Engine once
        try:
            engine = BasisSetEngine(shells)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to initialize Basis Engine: {e}")
            return
            
        tasks = []
        for item in items:
            try:
                mo_id = int(item.text(0))
            except: continue
            
            # Get coefficients from parser
            mo_data = coeffs_map.get(mo_id)
            if not mo_data: continue
            
            # Extract values list (assuming parser kept order)
            raw_coeffs = [x['coeff'] for x in mo_data['coeffs']]
            
            # Build dense vector (pad/truncate if needed)
            dense_vec = np.zeros(engine.n_basis)
            n = min(len(raw_coeffs), engine.n_basis)
            dense_vec[:n] = raw_coeffs[:n]
            
            fname = f"orbitals_mo_{mo_id}.cube"
            fpath = os.path.join(out_dir, fname)
            
            tasks.append({
                'mo_id': mo_id,
                'coeffs': dense_vec,
                'path': fpath
            })
            
        if not tasks: return

        # 5. Run Batch
        self.run_cube_batch(engine, tasks, out_dir)

    def run_cube_batch(self, engine, tasks, out_dir):
        # Interactive progress dialog
        from PyQt6.QtWidgets import QProgressDialog
        
        self.cube_pd = QProgressDialog("Generating Cubes...", "Cancel", 0, len(tasks)*100, self)
        self.cube_pd.setWindowModality(Qt.WindowModality.WindowModal)
        self.cube_pd.setMinimumDuration(0)
        
        self.cube_tasks = tasks
        self.cube_engine = engine
        self.cube_current_idx = 0
        self.cube_out_dir = out_dir
        
        self.start_next_cube_task()
        
    def start_next_cube_task(self):
        if self.cube_current_idx >= len(self.cube_tasks):
            self.cube_pd.setValue(len(self.cube_tasks)*100)
            self.cube_pd.close()
            
            # Ask to open Folder or load one?
            reply = QMessageBox.question(self, "Done", 
                f"Generated {len(self.cube_tasks)} cube files in:\n{self.cube_out_dir}\n\nOpen output directory?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
                
            if reply == QMessageBox.StandardButton.Yes and os.path.exists(self.cube_out_dir):
                import webbrowser
                webbrowser.open(self.cube_out_dir)
            return
            
        if self.cube_pd.wasCanceled():
            return
            
        task = self.cube_tasks[self.cube_current_idx]
        self.cube_pd.setLabelText(f"Generating MO {task['mo_id']}...")
        
        from .mo_engine import CalcWorker
        
        # Params
        n_points = 40 # Configurable?
        margin = 3.0
        
        parser = self.parent().parser
        
        worker = CalcWorker(
            self.cube_engine, 
            task['mo_id'], 
            n_points, 
            margin, 
            parser.data["atoms"], 
            parser.data["coords"], 
            task['coeffs'], 
            task['path']
        )
        
        worker.progress_sig.connect(self.on_cube_progress)
        worker.finished_sig.connect(self.on_cube_finished)
        
        self.current_cube_worker = worker
        worker.start()
        
    def on_cube_progress(self, pct):
        if self.cube_pd.wasCanceled():
            if self.current_cube_worker:
                self.current_cube_worker._is_cancelled = True
            return
            
        base = self.cube_current_idx * 100
        self.cube_pd.setValue(base + pct)
        
    def on_cube_finished(self, success, msg):
        if not success:
            QMessageBox.warning(self, "Error", f"Failed to generate cube: {msg}")
            
        self.cube_current_idx += 1
        self.start_next_cube_task()


    def view_coeffs(self):
        item = self.tree.currentItem()
        if not item: return
        
        # Column 0 is index
        try:
            mo_idx = int(item.text(0))
        except: return
        
        # Get coeffs
        mo_data = self.parent().parser.data.get("mo_coeffs", {}).get(mo_idx)
        if not mo_data or not mo_data.get('coeffs'):
            QMessageBox.information(self, "No Info", f"No coefficients found for MO {mo_idx}. \n(Check if 'MOLECULAR ORBITALS' block exists in output)")
            return
            
        dlg = MOCoeffDialog(self, mo_idx, mo_data)
        dlg.exec()

class MOCoeffDialog(QDialog):
    def __init__(self, parent, mo_idx, mo_data):
        super().__init__(parent)
        energy_str = f"E={mo_data.get('energy', 0.0):.4f} Eh"
        occ_str = f"Occ={mo_data.get('occ', 0.0):.2f}"
        spin_str = mo_data.get('spin', '').capitalize()
        
        self.setWindowTitle(f"MO {mo_idx} ({spin_str}) - {energy_str}, {occ_str}")
        self.resize(500, 600)
        
        coeffs = mo_data.get('coeffs', [])
        
        layout = QVBoxLayout(self)
        
        # Sort by abs(coeff) descending
        coeffs_sorted = sorted(coeffs, key=lambda x: abs(x['coeff']), reverse=True)
        
        from PyQt6.QtWidgets import QTableWidget, QTableWidgetItem
        self.table = QTableWidget()
        self.table.setColumnCount(4)
        self.table.setHorizontalHeaderLabels(["Atom", "Sym", "Orb", "Coeff"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.setRowCount(len(coeffs_sorted))
        
        for r, c in enumerate(coeffs_sorted):
            self.table.setItem(r, 0, QTableWidgetItem(str(c['atom_idx'])))
            self.table.setItem(r, 1, QTableWidgetItem(c['sym']))
            self.table.setItem(r, 2, QTableWidgetItem(c['orb']))
            coeff_val = c['coeff']
            item_c = QTableWidgetItem(f"{coeff_val:.5f}")
            # Colorize for visibility
            if coeff_val > 0:
                item_c.setForeground(QColor("blue"))
            else:
                item_c.setForeground(QColor("red"))
            self.table.setItem(r, 3, item_c)
            
        layout.addWidget(self.table)
        
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        layout.addWidget(btn_close)

    def launch_3d_viewer(self):
        # file_path comes from parent dialog
        # Access parent's file_path
        # self.parent() is the Dialog, but might need casting or just store logic
        
        # In __init__, parent is passed.
        # OrcaResultAnalyzerDialog has self.file_path
        
        if not hasattr(self.parent(), "file_path"):
            return
            
        base_path = os.path.splitext(self.parent().file_path)[0]
        fchk_path = base_path + ".fchk"
        
        if os.path.exists(fchk_path):
            # Ask confirmation if we should open it? 
            # Opening might close current if logic dictates, but usually plugins overlap or replace.
            # Best is to just request load.
            reply = QMessageBox.question(self, "Open FCHK", 
                                         f"Found {os.path.basename(fchk_path)}.\nOpen it for 3D orbital visualization?",
                                         QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
            
            if reply == QMessageBox.StandardButton.Yes:
                self.parent().mw.load_file(fchk_path)
                self.close() # Close MO list to avoid clutter
        else:
            QMessageBox.information(self, "Not Found", 
                                    f"Could not find corresponding FCHK file:\n{fchk_path}\n\n"
                                    "To generate it, run: orca_2mkl {base} -molden\n"
                                    "(Note: Gaussian MO Analyzer generates cubes from FCHK)")

    def show_thermal(self):
        data = self.parser.data.get("thermal", {})
        if not data:
            QMessageBox.warning(self, "No Info", "No thermochemistry section found.")
            return
            
        dlg = ThermalTableDialog(self, data)
        dlg.exec()

class ThermalTableDialog(QDialog):
    def __init__(self, parent, data):
        super().__init__(parent)
        self.setWindowTitle("Thermochemistry")
        self.resize(400, 300)
        
        layout = QVBoxLayout(self)
        
        from PyQt6.QtWidgets import QTableWidget, QTableWidgetItem
        
        self.table = QTableWidget()
        self.table.setColumnCount(2)
        self.table.setHorizontalHeaderLabels(["Property", "Value (Eh)"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        
        items = [
            ("Zero Point Energy", data.get("zpe")),
            ("Total Thermal Energy", data.get("thermal_energy")),
            ("Enthalpy (H)", data.get("enthalpy")),
            ("Entropy (S)", data.get("entropy")), # Note: S is usually T*S or S? In ORCA T*S is implicit in G=H-TS
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
        from PyQt6.QtWidgets import QApplication
        text = ""
        for r in range(self.table.rowCount()):
            p = self.table.item(r, 0).text()
            v = self.table.item(r, 1).text()
            text += f"{p}\t{v}\n"
        QApplication.clipboard().setText(text)

    def show_tddft(self):
        excitations = self.parser.data.get("tddft", [])
        if not excitations:
            QMessageBox.warning(self, "No Analysis", "No TDDFT/TDA excitation energies found.")
            return
            
        dlg = TDDFTDialog(self, excitations)
        dlg.exec()

class TDDFTDialog(QDialog):
    def __init__(self, parent, excitations):
        super().__init__(parent)
        self.setWindowTitle("TDDFT Spectrum")
        self.resize(900, 650)
        self.excitations = excitations
        
        layout = QVBoxLayout(self)
        
        # Spectrum Widget
        from .spectrum_widget import SpectrumWidget
        # Data format for widget: List of dicts. 
        # Our excitations have 'energy_nm' and 'osc'.
        self.spectrum = SpectrumWidget(self.excitations, x_key='energy_nm', y_key='osc', x_unit='Wavelength (nm)', sigma=20.0)
        layout.addWidget(self.spectrum)
        
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
        self.spin_sigma.setValue(20.0)
        self.spin_sigma.valueChanged.connect(self.spectrum.set_sigma)
        ctrl_row1.addWidget(self.spin_sigma)
        
        self.chk_sticks = QCheckBox("Sticks")
        self.chk_sticks.setChecked(True)
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
        if is_auto:
            self.spectrum.set_auto_range()
        else:
            self.update_range()
            
    def update_range(self):
        if self.chk_auto_y.isChecked(): return
        ymin = self.spin_y_min.value()
        ymax = self.spin_y_max.value()
        self.spectrum.set_y_range(ymin, ymax)
        
    def toggle_auto_x(self):
        is_auto = self.chk_auto_x.isChecked()
        self.spin_x_min.setEnabled(not is_auto)
        self.spin_x_max.setEnabled(not is_auto)
        if is_auto:
            self.spectrum.set_auto_x_range()
        else:
            self.update_x_range()
            
    def update_x_range(self):
        if self.chk_auto_x.isChecked(): return
        xmin = self.spin_x_min.value()
        xmax = self.spin_x_max.value()
        self.spectrum.set_x_range(xmin, xmax)
        
    def save_png(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save Graph", "", "Images (*.png)")
        if path:
            self.spectrum.save_png(path)
            
    def save_csv(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save Data", "", "CSV Files (*.csv)")
        if path:
            success = self.spectrum.save_csv(path)
            if success:
                QMessageBox.information(self, "Saved", f"Data saved to:\n{path}")
            else:
                QMessageBox.warning(self, "Error", "Failed to save CSV.")
                
    def switch_spectrum_type(self):
        is_cd = self.radio_cd.isChecked()
        if is_cd:
            # Switch to CD (Rotatory Strength)
            self.spectrum.y_key = 'rotatory_strength'
            self.spectrum.y_unit = 'R (10⁻⁴⁰ esu² cm²)'
        else:
            # Switch to Absorption (Oscillator Strength)
            self.spectrum.y_key = 'osc'
            self.spectrum.y_unit = 'Oscillator Strength'
        self.spectrum.update()

    def show_dipole(self):
        d = self.parser.data.get("dipole")
        if not d:
            QMessageBox.warning(self, "No Info", "No dipole moment found.")
            return
        
        # Create dipole dialog
        dlg = QDialog(self)
        dlg.setWindowTitle("Dipole Moment")
        dlg.resize(400, 200)
        
        layout = QVBoxLayout(dlg)
        
        # Dipole information
        info_text = f"""<b>Dipole Moment</b><br><br>
X: {d['vector'][0]:.4f} Debye<br>
Y: {d['vector'][1]:.4f} Debye<br>
Z: {d['vector'][2]:.4f} Debye<br>
<br>
<b>Magnitude: {d['magnitude']:.4f} Debye</b>"""
        
        lbl_info = QLabel(info_text)
        layout.addWidget(lbl_info)
        
        # Toggle 3D arrow
        chk_show_arrow = QCheckBox("Show dipole vector in 3D")
        chk_show_arrow.setChecked(False)
        
        dipole_arrow_actor = [None]  # Container for actor reference
        
        def toggle_arrow(state):
            # Remove old arrow
            if dipole_arrow_actor[0]:
                try:
                    self.mw.plotter.remove_actor(dipole_arrow_actor[0])
                    self.mw.plotter.render()
                except: pass
                dipole_arrow_actor[0] = None
            
            # Add arrow if checked
            if state:
                try:
                    import numpy as np
                    
                    # Get center of mass
                    coords = self.parser.data.get("coords", [])
                    if coords:
                        center = np.mean(coords, axis=0)
                    else:
                        center = np.array([0.0, 0.0, 0.0])
                    
                    # Dipole vector (normalized and scaled to 2 Angstrom)
                    dipole_vec = np.array(d['vector'])
                    magnitude = np.linalg.norm(dipole_vec)
                    if magnitude > 1e-6:
                        dipole_vec = (dipole_vec / magnitude) * 2.0
                    
                    # Add arrow
                    dipole_arrow_actor[0] = self.mw.plotter.add_arrows(
                        center.reshape(1, 3),
                        dipole_vec.reshape(1, 3),
                        mag=1.0,
                        color='cyan',
                        name='dipole_vector'
                    )
                    self.mw.plotter.render()
                except Exception as e:
                    print(f"Error adding dipole arrow: {e}")
        
        chk_show_arrow.stateChanged.connect(toggle_arrow)
        layout.addWidget(chk_show_arrow)
        
        # Close button
        btn_close = QPushButton("Close")
        
        def on_close():
            # Clean up arrow
            if dipole_arrow_actor[0]:
                try:
                    self.mw.plotter.remove_actor(dipole_arrow_actor[0])
                    self.mw.plotter.render()
                except: pass
            dlg.accept()
        
        btn_close.clicked.connect(on_close)
        layout.addWidget(btn_close)
        
        dlg.exec()

    def show_charges(self):
        charges = self.parser.data.get("charges", {})
        if not charges:
            QMessageBox.warning(self, "No Info", "No atomic charges found (Mulliken/Loewdin).")
            return
            
        if hasattr(self, 'charge_dlg') and self.charge_dlg is not None:
             self.charge_dlg.close()
             
        self.charge_dlg = ChargeDialog(self, charges)
        self.charge_dlg.show()
        # Note: ChargeDialog relies on parent (self) for context access.
        # Since self is now modeless (global ref), this is safe from GC.

    def show_nmr(self):
        data = self.parser.data.get("nmr_shielding", [])
        if not data:
            QMessageBox.warning(self, "No Info", "No NMR chemical shielding data found.")
        dlg = NMRDialog(self, data)
        dlg.exec() # NMR is just table, modal is OK

    def show_trajectory(self):
        data = self.parser.data.get("scan_steps", [])
        if not data:
             QMessageBox.warning(self, "No Info", "No trajectory steps (Scan/Optimization) found.")
             return
             
        if hasattr(self, 'traj_dlg') and self.traj_dlg is not None:
             self.traj_dlg.close()
             
        # self.mw is the main window instance passed in init
        self.traj_dlg = TrajectoryResultDialog(self.mw, data, title="Trajectory / Scan Analysis")
        self.traj_dlg.show()

    def show_forces(self):
        grads = self.parser.data.get("gradients", [])
        if not grads:
             QMessageBox.warning(self, "No Info", "No cartesian gradients found.")
             return
             
        dlg = ForceViewerDialog(self, grads)
        dlg.show() # Modeless relative to parent? Or modal? 
        # ForceViewer uses context from parent.
        # It handles closeEvent to clean up.
        # exec() blocks, so user can't rotate 3D view easily if not careful.
        # But for interactors, modeless is better.
        # However, dialogs usually default to modal or blocking unless Show() is used.
        # If I use show(), I need to keep reference.
        # Let's make it modal for safety (simpler) OR keep ref.
        # Reference plugin 'vector_viewer' used modeless.
        # But here, as a sub-sub-window, let's try modal first to ensure cleanup, 
        # OR just .exec(). 
        # Actually, user wants to see 3D view while tuning scale.
        # So .show() is better. 
        # But local variable 'dlg' will die. 
        # I need to store self.force_dlg
        
        if hasattr(self, 'force_dlg') and self.force_dlg is not None:
            self.force_dlg.close()
            
        self.force_dlg = ForceViewerDialog(self, grads)
        self.force_dlg.show()

class GradientBar(QWidget):
    def __init__(self, parent=None, colors=["red", "white", "blue"]):
        super().__init__(parent)
        self.colors = colors
        self.setFixedHeight(30)
        
    def set_colors(self, colors):
        self.colors = colors
        self.update()
        
    def paintEvent(self, event):
        painter = QPainter(self)
        grad =  self.get_gradient()
        painter.fillRect(self.rect(), grad)
        
        # Draw border
        painter.setPen(Qt.GlobalColor.black)
        painter.drawRect(0, 0, self.width()-1, self.height()-1)
        
    def get_gradient(self):
        from PyQt6.QtGui import QLinearGradient
        grad = QLinearGradient(0, 0, self.width(), 0)
        
        n = len(self.colors)
        if n < 2: return grad
        
        for i, c in enumerate(self.colors):
            pos = i / (n - 1)
            grad.setColorAt(pos, QColor(c))
        return grad



class ForceViewerDialog(QDialog):
    def __init__(self, parent_dlg, gradients):
        super().__init__(parent_dlg)
        self.setWindowTitle("Gradient (Force) Viewer")
        self.resize(300, 400)
        self.parent_dlg = parent_dlg
        self.gradients = gradients # List of {atom_idx, atom_sym, vector}
        self.actors = []
        
        layout = QVBoxLayout(self)
        
        layout.addWidget(QLabel(f"Visualizing gradients for {len(gradients)} atoms."))
        
        # Scale
        scale_layout = QHBoxLayout()
        scale_layout.addWidget(QLabel("Scale:"))
        self.spin_scale = QDoubleSpinBox()
        self.spin_scale.setRange(0.1, 100.0)
        self.spin_scale.setValue(5.0) # Forces are small usually
        self.spin_scale.setSingleStep(0.5)
        self.spin_scale.valueChanged.connect(self.update_vectors)
        scale_layout.addWidget(self.spin_scale)
        layout.addLayout(scale_layout)
        
        # Color
        self.btn_color = QPushButton("Color: Red")
        self.btn_color.setStyleSheet("background-color: red; color: white;")
        # No color picker for now, just red default
        layout.addWidget(self.btn_color) 
        
        # Close
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        layout.addWidget(btn_close)
        
        # Init
        self.update_vectors()
        
    def update_vectors(self):
        # Clear old
        self.clear_vectors()
        
        if not self.parent_dlg.context: return
        try:
            # We need direct plotter access efficiently.
            # Using context.get_main_window().plotter is best.
            mw = self.parent_dlg.context.get_main_window()
            if not hasattr(mw, 'plotter'): return
            
            import pyvista as pv
            import numpy as np
            
            # Helper to get atom pos
            mol = mw.current_mol
            if not mol: return
            conf = mol.GetConformer()
            
            scale = self.spin_scale.value()
            
            for item in self.gradients:
                idx = item['atom_idx']
                if idx < mol.GetNumAtoms():
                    pos = conf.GetAtomPosition(idx)
                    start = np.array([pos.x, pos.y, pos.z])
                    vec = np.array(item['vector'])
                    
                    mag = np.linalg.norm(vec)
                    if mag < 1e-6: continue
                    
                    # Gradient is derivative of Energy. Force is -Gradient.
                    # Usually we visualize Force (-Gradient) to show where atoms want to go.
                    # Let's assume user wants Force.
                    force = -vec * scale # Reverse direction
                    
                    # Centered on atom? Yes.
                    
                    actor = mw.plotter.add_arrow(start, force, mag=np.linalg.norm(force), color='red', resolution=10)
                    self.actors.append(actor)
            
            mw.plotter.render()
            
        except Exception as e:
            print(f"Error drawing vectors: {e}")
            
    def clear_vectors(self):
        if not self.parent_dlg.context: return
        mw = self.parent_dlg.context.get_main_window()
        if not hasattr(mw, 'plotter'): return
        
        for actor in self.actors:
            mw.plotter.remove_actor(actor)
        self.actors = []
        mw.plotter.render()
        
    def closeEvent(self, event):
        self.clear_vectors()
        super().closeEvent(event)


class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MplCanvas, self).__init__(self.fig)

class TrajectoryResultDialog(QDialog):

    def __init__(self, gl_widget, steps, title="Trajectory Analysis"):
        super().__init__()
        self.setWindowTitle(title)
        self.resize(900, 700)
        self.gl_widget = gl_widget
        self.steps = steps # List of {step, energy, atoms, coords}
        
        # Extract energies
        self.energies = [s['energy'] for s in self.steps]
        # Calculate relative energies (kcal/mol)
        min_e = min(self.energies) if self.energies else 0
        self.rel_energies = [(e - min_e)*627.509 for e in self.energies] 
        
        self.init_ui()
        self.plot_data()
        
    def init_ui(self):
        layout = QVBoxLayout(self)
        
        # 1. Matplotlib Canvas
        self.canvas = MplCanvas(self, width=5, height=4, dpi=100)
        layout.addWidget(self.canvas)
        
        # Connect events
        self.canvas.mpl_connect('pick_event', self.on_pick)
        self.canvas.mpl_connect('motion_notify_event', self.on_hover)
        
        # Tooltip annotation
        self.annot = self.canvas.axes.annotate("", xy=(0,0), xytext=(20,20),
                            textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="w", alpha=0.9),
                            arrowprops=dict(arrowstyle="->"))
        self.annot.set_visible(False)
        
        # 2. Controls
        ctrl_layout = QHBoxLayout()
        
        # Slider
        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.setRange(0, len(self.steps) - 1)
        self.slider.valueChanged.connect(self.on_step_changed)
        ctrl_layout.addWidget(QLabel("Step:"))
        ctrl_layout.addWidget(self.slider)
        
        self.lbl_info = QLabel("Step 0")
        ctrl_layout.addWidget(self.lbl_info)
        
        layout.addLayout(ctrl_layout)
        
        # 3. Buttons (Playback & Export)
        btn_layout = QHBoxLayout()
        
        self.btn_play = QPushButton("Play")
        self.btn_play.clicked.connect(self.toggle_play)
        btn_layout.addWidget(self.btn_play)
        
        btn_layout.addStretch()
        
        self.btn_save_img = QPushButton("Save Graph")
        self.btn_save_img.clicked.connect(self.save_graph)
        btn_layout.addWidget(self.btn_save_img)
        
        self.btn_save_csv = QPushButton("Save CSV")
        self.btn_save_csv.clicked.connect(self.save_csv)
        btn_layout.addWidget(self.btn_save_csv)
        
        self.btn_save_gif = QPushButton("Save GIF")
        self.btn_save_gif.clicked.connect(self.save_gif)
        self.btn_save_gif.setEnabled(HAS_PIL)
        btn_layout.addWidget(self.btn_save_gif)
        
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        btn_layout.addWidget(btn_close)
        
        layout.addLayout(btn_layout)
        
        # Timer
        self.timer = QTimer()
        self.timer.timeout.connect(self.next_frame)
        self.is_playing = False
        
    def plot_data(self):
        self.canvas.axes.clear()
        
        x = list(range(len(self.rel_energies)))
        y = self.rel_energies
        
        # Draw Line and Scatter
        self.canvas.axes.plot(x, y, 'b-', label='Energy', picker=5)
        self.scatter = self.canvas.axes.scatter(x, y, c='red', s=40, picker=5, zorder=5)
        
        self.canvas.axes.set_xlabel("Step")
        self.canvas.axes.set_ylabel("Relative Energy (kcal/mol)")
        self.canvas.axes.set_title("Energy Profile")
        self.canvas.axes.grid(True)
        
        self.highlight_point(0)
        self.canvas.draw()
        
    def highlight_point(self, idx):
        # Remove old markers
        if hasattr(self, '_highlight_marker'):
            try: self._highlight_marker.remove()
            except: pass
        if hasattr(self, '_highlight_line'):
            try: self._highlight_line.remove()
            except: pass
            
        x = idx
        y = self.rel_energies[idx]
        
        # Red circle
        self._highlight_marker, = self.canvas.axes.plot(x, y, 'ro', markersize=12, markeredgecolor='black', markeredgewidth=2, zorder=10)
        # Vertical Line
        self._highlight_line = self.canvas.axes.axvline(x=x, color='gray', linestyle='--', alpha=0.6)
        
        self.canvas.draw()
        
    def on_step_changed(self, idx):
        self.highlight_point(idx)
        step = self.steps[idx]
        val = self.rel_energies[idx]
        abs_e = self.energies[idx]
        self.lbl_info.setText(f"Step {idx+1}/{len(self.steps)} | E={abs_e:.6f} Eh ({val:.2f} kcal/mol)")
        
        self.update_structure(step['atoms'], step['coords'])
        
    def update_structure(self, atoms, coords):
        # RDKit build
        mol = Chem.RWMol()
        conf = Chem.Conformer()
        pt = Chem.GetPeriodicTable()
        try:
            for i, sym in enumerate(atoms):
               if ":" in sym: sym = sym.split(":")[0]
               an = pt.GetAtomicNumber(sym)
               mol.AddAtom(Chem.Atom(an))
               conf.SetAtomPosition(i, Point3D(coords[i][0], coords[i][1], coords[i][2]))
        except: return
        
        mol.AddConformer(conf)
        final_mol = mol.GetMol()
        
        if hasattr(self.gl_widget, 'draw_molecule_3d'):
             self.gl_widget.draw_molecule_3d(final_mol)
             
    def on_pick(self, event):
        if event.artist and hasattr(event, 'ind'):
            idx = event.ind[0] # Index of point
            self.slider.setValue(idx)
            
    def on_hover(self, event):
        vis = self.annot.get_visible()
        if event.inaxes == self.canvas.axes:
            cont, ind = self.scatter.contains(event)
            if cont:
                idx = ind['ind'][0]
                pos = self.scatter.get_offsets()[idx]
                self.annot.xy = pos
                val = self.rel_energies[idx]
                self.annot.set_text(f"Step {idx+1}\nRel E={val:.2f} kcal/mol")
                self.annot.set_visible(True)
                self.canvas.draw_idle()
                return
        
        if vis:
            self.annot.set_visible(False)
            self.canvas.draw_idle()
            
    def toggle_play(self):
        if self.is_playing:
            self.timer.stop()
            self.btn_play.setText("Play")
            self.is_playing = False
        else:
            self.timer.start(500)
            self.btn_play.setText("Pause")
            self.is_playing = True
            
    def next_frame(self):
        idx = self.slider.value() + 1
        if idx >= len(self.steps):
            idx = 0
        self.slider.setValue(idx)
        
    def save_graph(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save Graph", "", "Images (*.png *.jpg *.svg)")
        if path:
            self.canvas.fig.savefig(path, dpi=300)
            QMessageBox.information(self, "Saved", f"Graph saved to:\n{path}")

    def save_csv(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save CSV", "", "CSV Files (*.csv)")
        if path:
            import csv
            try:
                with open(path, 'w', newline='', encoding='utf-8') as f:
                    writer = csv.writer(f)
                    writer.writerow(["Step", "Energy_Hartree", "Rel_Energy_kcal_mol"])
                    for i, step in enumerate(self.steps):
                         writer.writerow([i+1, step['energy'], self.rel_energies[i]])
                QMessageBox.information(self, "Saved", f"Data saved to:\n{path}")
            except Exception as e:
                QMessageBox.critical(self, "Error", str(e))

    def save_gif(self):
        if not HAS_PIL:
            QMessageBox.warning(self, "Error", "PIL (Pillow) not installed.")
            return

        # Settings Dialog
        dialog = QDialog(self)
        dialog.setWindowTitle("GIF Settings")
        form = QFormLayout(dialog)
        
        spin_fps = QSpinBox()
        spin_fps.setRange(1, 60)
        spin_fps.setValue(10)
        form.addRow("FPS:", spin_fps)
        
        chk_trans = QCheckBox()
        chk_trans.setChecked(True)
        form.addRow("Transparent:", chk_trans)
        
        btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        btns.accepted.connect(dialog.accept)
        btns.rejected.connect(dialog.reject)
        form.addRow(btns)
        
        if dialog.exec() != QDialog.DialogCode.Accepted: return
        
        fps = spin_fps.value()
        transparent = chk_trans.isChecked()
        
        path, _ = QFileDialog.getSaveFileName(self, "Save GIF", "", "GIF Files (*.gif)")
        if not path: return
        if not path.lower().endswith('.gif'): path += '.gif'
        
        # Stop playback if running
        was_playing = self.is_playing
        if self.is_playing: self.toggle_play()
        
        # Capture logic
        import time
        images = []
        original_idx = self.slider.value()
        
        # Access Plotter from gl_widget? 
        # Ideally gl_widget IS the main window or has controller.
        # The constructor passed 'gl_widget' which in show_scan is 'self.mw'.
        mw = self.gl_widget
        
        if not hasattr(mw, 'plotter'):
            QMessageBox.warning(self, "Error", "Cannot access 3D plotter for capture.")
            return
            
        try:
            self.setCursor(Qt.CursorShape.WaitCursor)
            for i in range(len(self.steps)):
                self.slider.setValue(i)
                # Force update
                QApplication.processEvents()
                mw.plotter.render()
                
                img_array = mw.plotter.screenshot(transparent_background=transparent, return_img=True)
                if img_array is not None:
                     img = Image.fromarray(img_array)
                     if transparent:
                         img = img.convert("RGBA")
                     else:
                         img = img.convert("RGB").quantize(colors=256)
                     images.append(img)
            
            # Save GIF
            if images:
                duration = int(1000 / fps)
                # Basic save
                images[0].save(path, save_all=True, append_images=images[1:], duration=duration, loop=0, disposal=2)
                QMessageBox.information(self, "Success", f"GIF saved to:\n{path}")
                
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save GIF:\n{e}")
        finally:
            self.setCursor(Qt.CursorShape.ArrowCursor)
            self.slider.setValue(original_idx)
            if was_playing: self.toggle_play()

    def closeEvent(self, event):
        self.timer.stop()
        super().closeEvent(event)


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
        from PyQt6.QtWidgets import QTableWidget, QTableWidgetItem
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
        from PyQt6.QtWidgets import QApplication
        text = "Idx\tElem\tShielding\tShift\n"
        for r in range(self.table.rowCount()):
             cols = []
             for c in range(self.table.columnCount()):
                 it = self.table.item(r, c)
                 cols.append(it.text() if it else "")
             text += "\t".join(cols) + "\n"
        QApplication.clipboard().setText(text)

    def show_vib(self):
        freqs = self.parser.data.get("frequencies", [])
        if not freqs:
            QMessageBox.warning(self, "No Data", "No vibrational frequencies found.")
            return
            
        from .freq_analysis import FrequencyDialog
        dlg = FrequencyDialog(self.mw, freqs, self.parser.data["atoms"], self.parser.data["coords"])
        dlg.exec()

    def show_scan(self): QMessageBox.information(self, "Info", "Scan not implemented yet.")
