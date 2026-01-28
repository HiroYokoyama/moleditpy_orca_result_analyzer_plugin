
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
                 # print(f"Auto-loading Hessian: {hess_path}")
                 self.load_hessian_from_path(hess_path)
        
        # If still no Hessian data, try to populate from parser (Output file)
        if not self.hessian_data and self.parser:
             # Create dummy hessian_data from parser info
             # Gradients
             grads = self.parser.data.get("gradients", [])
             freqs = self.parser.data.get("frequencies", [])
             
             if grads or freqs:
                 self.hessian_data = {}
                 
                 # Gradients
                 if grads:
                     try:
                         flat_grads = []
                         sorted_grads = sorted(grads, key=lambda x: x.get('atom_idx', 0))
                         for g in sorted_grads:
                             flat_grads.extend(g['vector'])
                         self.hessian_data["gradient_vec"] = flat_grads
                     except: pass
                     
                 # Frequencies
                 if freqs:
                     # Check format (list of dicts vs list of floats)
                     try:
                         if isinstance(freqs[0], dict):
                             self.hessian_data["frequencies"] = [d['freq'] for d in freqs]
                             # Extract normal modes if available
                             # parser stores vectors in freq dicts
                             if 'vector' in freqs[0]:
                                 # hessian_data["normal_modes"] expects list of lists of values [x,y,z, x,y,z...]?
                                 # OR list of lists of tuples?
                                 # ForceViewerDialog.update_vectors expects:
                                 # mode_vecs = self.hessian_data["normal_modes"][mode_idx]
                                 # then iterates: for i, vec in enumerate(mode_vecs)
                                 # So it expects a LIST OF VECTORS (tuples or lists [x,y,z])
                                 
                                 # Parser stores: "vector": [(x,y,z), (x,y,z)...]
                                 # So we can just map it directly.
                                 self.hessian_data["normal_modes"] = [d['vector'] for d in freqs]
                         else:
                             self.hessian_data["frequencies"] = freqs
                     except: pass
                 
                 # Analyze and Update UI
                 self.analyze_stability()
                 if self.hessian_data.get("frequencies"):
                     self.populate_frequency_table()
                 self.update_vectors()

    def closeEvent(self, event):
        """Clean up vectors when dialog is closed"""
        if self.parent_dlg and hasattr(self.parent_dlg, 'mw') and hasattr(self.parent_dlg.mw, 'plotter'):
            plotter = self.parent_dlg.mw.plotter
            if plotter:
                for actor in self.actors:
                    try:
                        plotter.remove_actor(actor)
                    except: pass
                self.actors = []
                try:
                    plotter.render()
                except: pass
        super().closeEvent(event)
    
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
        # Hessian Load Buttons
        btn_load = QPushButton("Load Hessian File (.hess)")
        btn_load.clicked.connect(self.load_hessian_file)
        layout.addWidget(btn_load)
        
        # Frequencies Section
        freq_layout = QHBoxLayout()
        self.freq_label = QLabel("Frequencies:")
        freq_layout.addWidget(self.freq_label)
        
        freq_layout.addStretch()
        
        # Reset/Visualize Forces Button (Relocated here)
        btn_reset = QPushButton("Reset View (Show Forces)")
        btn_reset.setToolTip("Clear normal mode selection and show static force vectors")
        btn_reset.clicked.connect(self.visualize_hessian_forces)
        freq_layout.addWidget(btn_reset)
        
        layout.addLayout(freq_layout)
        

        
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
        layout.addWidget(QLabel("Force Constants (Diagonal of Hessian):"))
        self.hessian_table = QTableWidget()
        self.hessian_table.setColumnCount(3)
        self.hessian_table.setHorizontalHeaderLabels(["Atom Pair", "Force Constant (k)", "Unit"])
        h2 = self.hessian_table.header() if hasattr(self.hessian_table, "header") else self.hessian_table.horizontalHeader()
        h2.setStretchLastSection(True)
        layout.addWidget(self.hessian_table)
        
        # Add hint label to Hessian tab too (so user can see status updates)
        layout.addStretch()
        self.hessian_hint_label = QLabel("")
        self.hessian_hint_label.setWordWrap(True)
        self.hessian_hint_label.setStyleSheet("color: #0066cc; font-style: italic;")
        layout.addWidget(self.hessian_hint_label)
    
        if self.hessian_data:
             self.update_stability_summary()
             self.populate_frequency_table()
             self.populate_hessian_table()

    def visualize_hessian_forces(self):
         """Visualize forces (negated gradients)"""
         try:
             # Show hint immediately to confirm click
             if self.hessian_data and self.hessian_data.get("gradient_vec"):
                 count = len(self.hessian_data.get("gradient_vec"))
                 msg = f"Displaying Forces from Hessian File ({count} values)"
                 self.hint_label.setText(msg)
                 self.hessian_hint_label.setText(msg)
             elif self.gradients:
                 count = len(self.gradients)
                 msg = f"Displaying Forces from Output File ({count} atoms)"
                 self.hint_label.setText(msg)
                 self.hessian_hint_label.setText(msg)
             else:
                 msg = "No gradient data found (Structure likely optimized).\nSelect a frequency to see normal modes."
                 self.hint_label.setText(msg)
                 self.hessian_hint_label.setText(msg)
                 QMessageBox.information(self, "No Gradient Data", "No gradient forces found.\nThis is normal for an optimized structure (Forces ≈ 0).\n\nPlease select a Frequency in the table to visualize vibrational modes.")
                 return

             # Clear any normal mode selection to force 'Force' mode in update_vectors
             self.freq_table.blockSignals(True) # Prevent double update
             self.freq_table.clearSelection()
             self.freq_table.blockSignals(False)
             
             # Force update
             self.update_vectors()
             
         except Exception as e:
             import traceback
             QMessageBox.critical(self, "Error", f"An error occurred in visualize_hessian_forces:\n{str(e)}\n\n{traceback.format_exc()}")
 
         # Let's reuse update_vectors by adding a mode? 
         # The existing update_vectors checks chk_vectors and tables.
         
         # BETTER: Just call render directly here since it's a specific action
         scale = self.spin_scale.value()
         
         # Clear existing
         for actor in self.actors:
             try:
                 self.parent_dlg.plotter.remove_actor(actor)
             except: pass
         self.actors.clear()
         
         # Draw
         # Convert to pyvista vectors
         # We need start positions (current geometry)
         # If existing geometry is loaded in main window, use that.
         if not hasattr(self.parent_dlg, "atom_coords") or not self.parent_dlg.atom_coords:
             return
             
         coords = self.parent_dlg.atom_coords
         
         target_data = [] # (start, vec)
         
         if using_fallback:
             # Use self.gradients
             for item in self.gradients:
                 atom_idx = item.get('atom_idx')
                 if atom_idx is not None and atom_idx < len(coords):
                     vec = item.get('grad', item.get('vector'))
                     if vec:
                         target_data.append((coords[atom_idx], vec))
         else:
             # Use Hessian gradients
             # Reshape flat list to list of (x,y,z)
             forces = gradients
             n_atoms = len(coords)
             for i in range(n_atoms):
                 idx = i * 3
                 if idx + 2 < len(forces):
                     gx, gy, gz = forces[idx], forces[idx+1], forces[idx+2]
                     target_data.append((coords[i], [gx, gy, gz]))

         # Create vectors
         for i, (start, raw_vec) in enumerate(target_data):
             fx, fy, fz = raw_vec
             # Force is NEGATIVE gradient
             vector = np.array([-fx, -fy, -fz])
             mag = np.linalg.norm(vector)
             
             if mag < 1e-6: continue
             
             # Scale
             scaled_vec = vector * scale
             
             arrow = pv.Arrow(start=start, direction=scaled_vec, scale=np.linalg.norm(scaled_vec), 
                              shaft_resolution=self.force_res, tip_resolution=self.force_res,
                              tip_length=0.25, tip_radius=0.1, shaft_radius=0.05)
             
             actor = self.parent_dlg.plotter.add_mesh(arrow, color=self.force_color)
             self.actors.append(actor)
             
         self.parent_dlg.plotter.render()
         msg = f"Visualizing forces.\nScale: {scale}"
         if using_fallback:
             msg += "\n(Source: Output File Gradients)"
         else:
             msg += "\n(Source: Hessian $gradient block)"
         QMessageBox.information(self, "Visualization", msg)

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

        # Fallback for Gradients (Crucial for Opt outputs)
        if not self.hessian_data.get("gradient_vec") and self.gradients:
             # self.gradients is list of dicts: [{'atom_idx': 0, 'vector': [x,y,z]}, ...]
             # hessian gradient_vec is flat list [x0,y0,z0, x1,y1,z1...]
             try:
                 flat_grads = []
                 # Sort by atom_idx just in case
                 sorted_grads = sorted(self.gradients, key=lambda x: x.get('atom_idx', 0))
                 for g in sorted_grads:
                     flat_grads.extend(g['vector'])
                 self.hessian_data["gradient_vec"] = flat_grads
             except Exception as e:
                 print(f"Error converting gradients fallback: {e}")

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
        """Display force constants and frequencies"""
        self.populate_frequency_table()
        self.populate_hessian_table()

    def populate_frequency_table(self):
        """Populate the frequency table"""
        if not self.hessian_data: return
        
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

    def populate_hessian_table(self):
        """Populate the hessian force constants table"""
        if not self.hessian_data: return

        matrix = self.hessian_data.get("matrix", [])
        atoms = self.hessian_data.get("atoms", [])
        
        # Fallback: Use parser atoms if Hessian has none
        if not atoms and self.parser and "atoms" in self.parser.data:
            atoms = self.parser.data["atoms"]
            
        if not matrix or not atoms: return
        
        # Clear table
        self.hessian_table.setRowCount(0)
        
        # Display diagonal elements (force constants for each coordinate)
        n_atoms = len(atoms)
        row = 0
        
        # Limit rows to avoid freezing if massive
        max_rows = 300 # arbitrary limit for display
        
        for i in range(min(len(matrix), n_atoms * 3)):
            if row >= max_rows: break
            
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
