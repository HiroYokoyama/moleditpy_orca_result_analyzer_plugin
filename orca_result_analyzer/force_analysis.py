
import numpy as np
import pyvista as pv
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
                             QDoubleSpinBox, QGroupBox, QSpinBox, QColorDialog,
                             QFileDialog, QMessageBox, QTableWidget, QTableWidgetItem,
                             QTabWidget, QWidget, QFormLayout, QCheckBox)
from PyQt6.QtGui import QColor
import os

class ForceViewerDialog(QDialog):
    def __init__(self, parent_dlg, gradients, parser=None):
        super().__init__(parent_dlg)
        self.setWindowTitle("Geometry Stability Analysis")
        self.resize(600, 700)
        self.parent_dlg = parent_dlg
        self.gradients = gradients # From .out file
        self.parser = parser
        self.actors = []
        self.force_color = "red"
        self.force_res = 20
        self.hessian_data = None
        self.force_res = 20
        self.hessian_data = None
        
        main_layout = QVBoxLayout(self)
        
        # 1. Stability Summary Dashboard
        self.summary_grp = QGroupBox("Stability Summary")
        self.summary_layout = QFormLayout(self.summary_grp)
        self.lbl_status = QLabel("Load Hessian to analyze...")
        self.lbl_max_force = QLabel("-")
        self.lbl_rms_force = QLabel("-")
        self.lbl_imag_freq = QLabel("-")
        
        self.summary_layout.addRow("Stationary Point:", self.lbl_status)
        self.summary_layout.addRow("Imaginary Freqs:", self.lbl_imag_freq)
        self.summary_layout.addRow("Max Force:", self.lbl_max_force)
        self.summary_layout.addRow("RMS Force:", self.lbl_rms_force)
        
        main_layout.addWidget(self.summary_grp)
        
        # Tabs for Gradients and Hessian
        self.tabs = QTabWidget()
        main_layout.addWidget(self.tabs)
        
        # Tab 1: Hessian/Force Constants (Prioritized for Stability Analysis)
        hess_widget = QWidget()
        hess_layout = QVBoxLayout(hess_widget)
        self._setup_hessian_tab(hess_layout)
        self.tabs.addTab(hess_widget, "Frequencies & Constants")
        
        # Tab 2: Gradient Visualization
        grad_widget = QWidget()
        grad_layout = QVBoxLayout(grad_widget)
        self._setup_gradient_tab(grad_layout)
        self.tabs.addTab(grad_widget, "3D Forces")
             
        # Close button
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        main_layout.addWidget(btn_close)
        
        # Init
        self.update_vectors()
        
        # Auto-load Hessian if exists
        if self.parser and hasattr(self.parser, "filename") and self.parser.filename:
            base, _ = os.path.splitext(self.parser.filename)
            hess_path = base + ".hess"
            if os.path.exists(hess_path):
                 #print(f"Auto-loading Hessian: {hess_path}")
                 self.load_hessian_from_path(hess_path)
    
    def _setup_gradient_tab(self, layout):
        """Setup the gradient visualization tab"""
        # 3D Appearance Section
        view_group = QGroupBox("3D Visualization Controls")
        view_layout = QVBoxLayout(view_group)
        
        # Visibility Toggle
        self.chk_vectors = QCheckBox("Show Force Vectors")
        self.chk_vectors.setChecked(True)
        self.chk_vectors.stateChanged.connect(self.update_vectors)
        view_layout.addWidget(self.chk_vectors)
        
        # Scale row
        scale_row = QHBoxLayout()
        scale_row.addWidget(QLabel("Scale:"))
        self.spin_scale = QDoubleSpinBox()
        self.spin_scale.setRange(0.01, 1000.0)
        self.spin_scale.setValue(1.0) 
        self.spin_scale.setSingleStep(0.1)
        self.spin_scale.valueChanged.connect(self.update_vectors)
        scale_row.addWidget(self.spin_scale)
        
        scale_row.addWidget(QLabel(" Res:"))
        self.spin_res = QSpinBox()
        self.spin_res.setRange(3, 100)
        self.spin_res.setValue(20)
        self.spin_res.valueChanged.connect(self.on_res_changed)
        scale_row.addWidget(self.spin_res)
        view_layout.addLayout(scale_row)
        
        # Color row
        color_row = QHBoxLayout()
        color_row.addWidget(QLabel("Color:"))
        self.btn_color = QPushButton()
        self.btn_color.setFixedWidth(60)
        self.btn_color.setStyleSheet(f"background-color: {self.force_color}; border: 1px solid gray; height: 20px;")
        self.btn_color.clicked.connect(self.pick_color)
        color_row.addWidget(self.btn_color)
        color_row.addStretch()
        view_layout.addLayout(color_row)
        
        layout.addWidget(view_group)
        
        layout.addStretch()
        
        # Legend / Hints
        self.hint_label = QLabel("Displaying force vectors (Negative Gradient).\nArrows point where atoms want to move.")
        self.hint_label.setWordWrap(True)
        layout.addWidget(self.hint_label)
        layout.addStretch()
    
    def _setup_hessian_tab(self, layout):
        """Setup the Hessian/Force Constants tab"""
        # Load button
        btn_load = QPushButton("Load Hessian File (.hess)")
        btn_load.clicked.connect(self.load_hessian_file)
        layout.addWidget(btn_load)
        
        # Frequencies Table
        self.freq_label = QLabel("Frequencies:")
        layout.addWidget(self.freq_label)
        
        self.freq_table = QTableWidget()
        self.freq_table.setColumnCount(2)
        self.freq_table.setHorizontalHeaderLabels(["Mode", "Frequency (cm⁻¹)"])
        self.freq_table.setAlternatingRowColors(True)
        self.freq_table.setFixedHeight(200)
        header = self.freq_table.header() if hasattr(self.freq_table, "header") else self.freq_table.horizontalHeader()
        header.setStretchLastSection(True)
        self.freq_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.freq_table.setSelectionMode(QTableWidget.SelectionMode.SingleSelection)
        self.freq_table.itemSelectionChanged.connect(self.update_vectors)
        layout.addWidget(self.freq_table)
        
        # Force Constants Table
        layout.addWidget(QLabel("Force Constants (Diagonal):"))
        self.hessian_table = QTableWidget()
        self.hessian_table.setColumnCount(3)
        self.hessian_table.setHorizontalHeaderLabels(["Atom Pair", "Force Costant", "Unit"])
        h2 = self.hessian_table.header() if hasattr(self.hessian_table, "header") else self.hessian_table.horizontalHeader()
        h2.setStretchLastSection(True)
        layout.addWidget(self.hessian_table)
    
    def load_hessian_file(self):
        """Open dialog to load a Hessian file"""
        # Get directory from parent
        start_dir = ""
        if hasattr(self.parent_dlg, 'file_path'):
            start_dir = os.path.dirname(self.parent_dlg.file_path)
        
        filepath, _ = QFileDialog.getOpenFileName(
            self, "Open Hessian File", start_dir, "Hessian Files (*.hess);;All Files (*.*)"
        )
        
        if filepath:
            self.load_hessian_from_path(filepath)

    def load_hessian_from_path(self, filepath):
        """Load and parse a Hessian file from path"""
        if not self.parser:
            # Should have parser, but if not, can't fallback properly
            pass
            
        # Parse the file
        self.hessian_data = self.parser.parse_hessian_file(filepath)
        
        # Validate and Fallback
        if not self.hessian_data:
            QMessageBox.warning(self, "Parse Error", f"Failed to parse Hessian file:\n{filepath}")
            return
            
        # Fallback for Frequencies (common in Opt-only jobs)
        if not self.hessian_data.get("frequencies") and self.parser:
             p_freqs = self.parser.data.get("frequencies", [])
             if p_freqs:
                 # parser freqs are dicts: {'freq': val, ...}
                 # hessian freqs are floats: [val, ...]
                 try:
                     # Check if it's list of dicts or floats (parser might change)
                     if isinstance(p_freqs[0], dict):
                         self.hessian_data["frequencies"] = [d['freq'] for d in p_freqs]
                     else:
                         self.hessian_data["frequencies"] = p_freqs
                 except: pass

        # Analyze Stability
        self.analyze_stability()
        
        # Display Data
        self.display_force_constants()
        self.update_vectors() # Update vectors using Hessian gradients if available
        
        # Only show success message if called interactively? 
        # For auto-load maybe silent? But helpful to know.
        # Let's show it but maybe smaller or relying on the UI update.
        # Actually, if auto-load, user sees the table filled.
        if self.isVisible():
            QMessageBox.information(self, "Success", f"Loaded Hessian file successfully!\n{os.path.basename(filepath)}")
    
    def analyze_stability(self):
        if not self.hessian_data: return
        
        # 1. Stationary Point Classification
        freqs = self.hessian_data.get("frequencies", [])
        imag_count = 0
        for f in freqs:
            if f < 0: imag_count += 1
            elif f == 0: pass # Rotation/Translation usually near 0, but exactly 0 is rare in output
            
        self.lbl_imag_freq.setText(str(imag_count))
        
        if imag_count == 0:
            self.lbl_status.setText("Local Minimum (Stable)")
            self.lbl_status.setStyleSheet("color: green; font-weight: bold;")
        elif imag_count == 1:
            self.lbl_status.setText("Transition State (First-Order Saddle Point)")
            self.lbl_status.setStyleSheet("color: orange; font-weight: bold;")
        else:
            self.lbl_status.setText(f"Higher-Order Saddle Point ({imag_count} imag)")
            self.lbl_status.setStyleSheet("color: red; font-weight: bold;")
            
        # 2. Force Metrics
        # Check if we have gradient vector in hessian data
        grad_vec = self.hessian_data.get("gradient_vec", [])
        if grad_vec:
            # list of 3N floats
            forces = np.array(grad_vec)
            # Reshape to (N, 3)
            if len(forces) % 3 == 0:
                forces_3d = forces.reshape(-1, 3)
                # Max Force (max component? or max norm?)
                # ORCA usually reports Max Force as max(|grad component|)
                max_force = np.max(np.abs(forces))
                rms_force = np.sqrt(np.mean(forces**2))
                
                self.lbl_max_force.setText(f"{max_force:.6f} Eh/Bohr")
                self.lbl_rms_force.setText(f"{rms_force:.6f} Eh/Bohr")
            else:
                self.lbl_max_force.setText("Invalid Dim")
        else:
            self.lbl_max_force.setText("No Gradient Data")
            self.lbl_rms_force.setText("-")

    def display_force_constants(self):
        """Display force constants from Hessian matrix"""
        if not self.hessian_data: return
        
        matrix = self.hessian_data.get("matrix", [])
        atoms = self.hessian_data.get("atoms", [])
        
        # Fallback: Use parser atoms if Hessian has none
        if not atoms and self.parser and "atoms" in self.parser.data:
            atoms = self.parser.data["atoms"]
            
        # Display frequencies
        freqs = self.hessian_data.get("frequencies", [])
        self.freq_table.setRowCount(0)
        self.freq_table.setRowCount(len(freqs))
        if freqs:
            self.freq_label.setText(f"Frequencies ({len(freqs)} modes):")
            for i, freq in enumerate(freqs):
                self.freq_table.setItem(i, 0, QTableWidgetItem(f"Mode {i+1}")) # 1-based mode
                
                item = QTableWidgetItem(f"{freq:.2f}")
                if freq < 0:
                     item.setForeground(QColor("red"))
                self.freq_table.setItem(i, 1, item)
        else:
             self.freq_label.setText("Frequencies: None found")
        
        if not matrix or not atoms: return
        
        # Clear table
        self.hessian_table.setRowCount(0)
        
        # Display diagonal elements (force constants for each coordinate)
        n_atoms = len(atoms)
        row = 0
        
        for i in range(min(len(matrix), n_atoms * 3)):
            atom_idx = i // 3
            coord_idx = i % 3
            coord_name = ['X', 'Y', 'Z'][coord_idx]
            
            if atom_idx < len(atoms) and i < len(matrix) and i < len(matrix[i]):
                atom_label = f"{atoms[atom_idx]}{atom_idx + 1}" # 1-based atom index
                force_const = matrix[i][i] if i < len(matrix[i]) else 0.0
                
                self.hessian_table.insertRow(row)
                self.hessian_table.setItem(row, 0, QTableWidgetItem(f"{atom_label} - {coord_name}"))
                self.hessian_table.setItem(row, 1, QTableWidgetItem(f"{force_const:.6f}"))
                self.hessian_table.setItem(row, 2, QTableWidgetItem("Eh/Bohr²"))
                row += 1
        
        self.hessian_table.resizeColumnsToContents()
        
    def update_vectors(self):
        """Update the force vectors in the visualizer"""
        if not hasattr(self, 'spin_scale'): return
            
        try:
            mw = None
            if hasattr(self.parent_dlg, 'context') and self.parent_dlg.context:
                 mw = self.parent_dlg.context.get_main_window()
            elif hasattr(self.parent_dlg, 'mw'):
                 mw = self.parent_dlg.mw
            
            if not mw or not hasattr(mw, 'plotter'): return
            
            # Clear old
            self.clear_vectors()
            
            # Check visibility
            if hasattr(self, 'chk_vectors') and not self.chk_vectors.isChecked():
                return
            
            # Get settings
            scale = self.spin_scale.value()
            
            target_data = [] # List of (start_pos, vector)
            is_mode = False
            
            # 1. Check if a Normal Mode is selected in the table
            selected_rows = self.freq_table.selectedItems()
            if selected_rows and self.hessian_data and self.hessian_data.get("normal_modes"):
                mode_idx = selected_rows[0].row()
                if mode_idx < len(self.hessian_data["normal_modes"]):
                    is_mode = True
                    mode_vecs = self.hessian_data["normal_modes"][mode_idx]
                    
                    # Determine Atoms/Coords
                    coords = []
                    if self.parser and "coords" in self.parser.data:
                         coords = self.parser.data["coords"]
                    
                    if coords:
                        for i, vec in enumerate(mode_vecs):
                            if i < len(coords):
                                target_data.append((coords[i], vec))
                    
                    self.hint_label.setText(f"Displaying Normal Mode {mode_idx + 1} ({self.hessian_data['frequencies'][mode_idx]:.2f} cm⁻¹)")
            
            # 2. Fallback to Forces if no mode selected
            if not is_mode:
                # Decide Source: Gradients (Output) vs Hessian
                use_hessian = False
                grad_vecs = []
                
                if self.hessian_data and self.hessian_data.get("gradient_vec"):
                     use_hessian = True
                     raw_vecs = self.hessian_data.get("gradient_vec")
                     # Chunk into (x,y,z)
                     for i in range(0, len(raw_vecs), 3):
                         if i+2 < len(raw_vecs):
                             grad_vecs.append([raw_vecs[i], raw_vecs[i+1], raw_vecs[i+2]])
                
                # Determine Atoms
                coords = []
                if self.parser and "coords" in self.parser.data:
                     coords = self.parser.data["coords"]
                
                if coords:
                    if use_hessian:
                        # Map using index
                        for i, vec in enumerate(grad_vecs):
                            if i < len(coords):
                                target_data.append((coords[i], vec))
                    else:
                        # Use current gradients
                        if self.gradients:
                            for item in self.gradients:
                                atom_idx = item.get('atom_idx')
                                if atom_idx is not None and atom_idx < len(coords):
                                    vec = item.get('grad', item.get('vector'))
                                    if vec:
                                        target_data.append((coords[atom_idx], vec))
                
                self.hint_label.setText("Displaying force vectors (Negative Gradient).\nArrows point where atoms want to move.")

            # Draw
            for idx, (start, raw_vec) in enumerate(target_data):
                vec = np.array(raw_vec)
                length = np.linalg.norm(vec)
                if length < 1e-12: continue
                
                # For modes, we just scale the displacement. For forces, it's -grad
                if is_mode:
                    disp = vec * scale
                else:
                    disp = -vec * scale # Force = -Gradient
                
                disp_mag = np.linalg.norm(disp)
                if disp_mag < 1e-12: continue
                
                disp_dir = disp / disp_mag
                
                # Use add_mesh with pv.Arrow for robustness
                arrow = pv.Arrow(start=start, direction=disp_dir, scale=disp_mag,
                                 shaft_resolution=self.force_res, tip_resolution=self.force_res)
                actor = mw.plotter.add_mesh(arrow, color=self.force_color, name=f'force_{idx}')
                self.actors.append(actor)
            
            mw.plotter.render()
            
        except Exception as e:
            print(f"Error drawing vectors: {e}")

    def pick_color(self):
        color = QColorDialog.getColor(QColor(self.force_color), self, "Select Force Vector Color")
        if color.isValid():
            self.force_color = color.name()
            self.btn_color.setStyleSheet(f"background-color: {self.force_color}; border: 1px solid gray; height: 20px;")
            self.update_vectors()

    def on_res_changed(self, val):
        self.force_res = val
        self.update_vectors()
            
    def clear_vectors(self):
        mw = None
        if hasattr(self.parent_dlg, 'context') and self.parent_dlg.context:
             mw = self.parent_dlg.context.get_main_window()
        elif hasattr(self.parent_dlg, 'mw'):
             mw = self.parent_dlg.mw
             
        if not mw or not hasattr(mw, 'plotter'): return
        
        for actor in self.actors:
            try:
                mw.plotter.remove_actor(actor)
            except: pass
        self.actors = []
        mw.plotter.render()
        
    def closeEvent(self, event):
        self.clear_vectors()
        super().closeEvent(event)
