import os
import numpy as np
import pyvista as pv
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
                             QDoubleSpinBox, QGroupBox, QSlider, QTableWidget, 
                             QTableWidgetItem, QCheckBox, QAbstractItemView)
from PyQt6.QtCore import Qt, QTimer

class ForceViewerDialog(QDialog):
    def __init__(self, parent_dlg, gradients, parser=None):
        super().__init__(parent_dlg)
        self.setWindowTitle("Force Analysis")
        self.resize(700, 500)
        self.parent_dlg = parent_dlg
        self.gradients = gradients  # Current gradients from .out file
        self.parser = parser
        self.actors = []
        self.force_color = "red"
        
        # Get trajectory steps if available
        self.traj_steps = []
        self.current_step_idx = -1  # -1 means current/final structure
        if self.parser and "scan_steps" in self.parser.data:
            self.traj_steps = self.parser.data["scan_steps"]
        
        main_layout = QVBoxLayout(self)
        
        # Visualization Controls (Restoring deleted code)
        view_group = QGroupBox("Visualization Settings")
        view_layout = QVBoxLayout(view_group)
        
        # Row 1: Visualize Button and Scale
        row1 = QHBoxLayout()
        
        self.btn_visualize = QPushButton("Visualize Forces")
        self.btn_visualize.setCheckable(True) # Use checkable for toggle state
        self.btn_visualize.clicked.connect(self.toggle_visualization)
        row1.addWidget(self.btn_visualize)
        
        row1.addSpacing(20)
        row1.addWidget(QLabel("Vector Scale:"))
        self.spin_scale = QDoubleSpinBox()
        self.spin_scale.setRange(0.01, 100000.0) # Increased range for small gradients
        self.spin_scale.setValue(1.0)
        self.spin_scale.setSingleStep(0.1)
        self.spin_scale.valueChanged.connect(self.update_vectors)
        row1.addWidget(self.spin_scale)
        
        # Auto Scale Button
        btn_autoscale = QPushButton("Auto Scale")
        btn_autoscale.clicked.connect(self.auto_scale)
        row1.addWidget(btn_autoscale)
        
        self.chk_auto_scale = QCheckBox("Auto Apply")
        self.chk_auto_scale.setToolTip("Automatically auto-scale vectors when step changes")
        row1.addWidget(self.chk_auto_scale)
        
        row1.addStretch()
        
        view_layout.addLayout(row1)
        main_layout.addWidget(view_group)
        
        # Force and Gradient Table
        table_label = QLabel("Force and Gradient Data (Eh/Bohr):")
        main_layout.addWidget(table_label)
        
        self.force_table = QTableWidget()
        self.force_table.setColumnCount(8)
        self.force_table.setHorizontalHeaderLabels([
            "Atom", "Grad X", "Grad Y", "Grad Z", "Force X", "Force Y", "Force Z", "Force Mag"
        ])
        self.force_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        header = self.force_table.horizontalHeader()
        header.setStretchLastSection(True)
        main_layout.addWidget(self.force_table)

        # Trajectory Navigation (if available) - Moved AFTER table creation
        # (becuase controls trigger update which updates table)
        if self.traj_steps:
            self._setup_trajectory_controls(main_layout)
        
        # Buttons row
        buttons_layout = QHBoxLayout()
        
        btn_reload = QPushButton("Reload from File")
        btn_reload.setToolTip("Reload data from the ORCA output file")
        btn_reload.clicked.connect(self.reload_data)
        buttons_layout.addWidget(btn_reload)
        
        buttons_layout.addStretch()
        
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        buttons_layout.addWidget(btn_close)
        
        main_layout.addLayout(buttons_layout)
        
        # Initialize
        self.populate_force_table()
        # self.auto_scale() # Manual only per user request

    def toggle_visualization(self):
        """Toggle force vector visualization"""
        if self.btn_visualize.isChecked():
            self.btn_visualize.setText("Clear Forces")
            self.update_vectors()
        else:
            self.btn_visualize.setText("Visualize Forces")
            self.clear_vectors()

    def auto_scale(self):
        """Calculate and set an appropriate scale factor"""
        if not self.gradients: 
            return
            
        try:
            # Find max magnitude
            max_mag = 0.0
            for g in self.gradients:
                vec = g.get('grad', g.get('vector'))
                if vec:
                    mag = np.linalg.norm(vec)
                    if mag > max_mag: max_mag = mag
            
            if max_mag > 1e-12:
                # Target max vector length ~ 1.5 - 2.0 units (Angstrom)
                target = 2.0
                new_scale = target / max_mag
                
                # Round nicely
                if new_scale > 100:
                    new_scale = round(new_scale, -int(np.log10(new_scale))+1)
                else:
                    new_scale = round(new_scale, 2)
                    
                self.spin_scale.blockSignals(True)
                self.spin_scale.setValue(new_scale)
                self.spin_scale.blockSignals(False)
                
                # If currently visualizing, update
                if self.btn_visualize.isChecked():
                    self.update_vectors()
        except: pass
    
    def _setup_trajectory_controls(self, layout):
        """Setup trajectory navigation controls"""
        traj_group = QGroupBox(f"Optimization Trajectory ({len(self.traj_steps)} steps)")
        traj_layout = QVBoxLayout(traj_group)
        
        # Slider
        slider_layout = QHBoxLayout()
        slider_layout.addWidget(QLabel("Step:"))
        
        self.traj_slider = QSlider(Qt.Orientation.Horizontal)
        self.traj_slider.setRange(0, len(self.traj_steps)) # 0 to N-1 (traj), N = final

        self.traj_label = QLabel("Current (Final)")
        
        # Info labels (define early)
        self.traj_info = QLabel("Showing final optimized structure")
        self.traj_conv = QLabel("")
        self.traj_conv.setStyleSheet("color: #444; font-size: 9pt;")
        self.traj_conv.setWordWrap(True)

        # Navigation buttons
        btn_prev = QPushButton("<")
        btn_prev.setFixedWidth(30)
        btn_prev.clicked.connect(lambda: self.traj_slider.setValue(max(0, self.traj_slider.value() - 1)))
        
        btn_next = QPushButton(">")
        btn_next.setFixedWidth(30)
        btn_next.clicked.connect(lambda: self.traj_slider.setValue(min(self.traj_slider.maximum(), self.traj_slider.value() + 1)))
        
        slider_layout.addWidget(btn_prev)
        slider_layout.addWidget(self.traj_slider)
        slider_layout.addWidget(btn_next)
        slider_layout.addWidget(self.traj_label)
        
        # Jump to Final button
        btn_final = QPushButton("Jump to Final")
        btn_final.setToolTip("Jump to final structure with force data")
        btn_final.clicked.connect(lambda: self.traj_slider.setValue(len(self.traj_steps)))
        slider_layout.addWidget(btn_final)
        
        traj_layout.addLayout(slider_layout)
        
        traj_layout.addWidget(self.traj_info)
        traj_layout.addWidget(self.traj_conv)
        
        # Connect signals AFTER UI construction
        self.traj_slider.valueChanged.connect(self.on_trajectory_change)
        
        # Set initial value (triggers signal)
        self.traj_slider.setValue(len(self.traj_steps)) 
        
        # Initial trigger if value didn't change (e.g. if we started at 0 and desired is 0)
        if self.traj_slider.value() == 0 and len(self.traj_steps) == 0:
             self.on_trajectory_change(0)
        
        # Insert at the top (index 0) so it appears above other controls
        # even though it is initialized last.
        layout.insertWidget(0, traj_group)
    
    def on_trajectory_change(self, val):
        """Handle trajectory step change"""
        self.current_step_idx = val
        num_steps = len(self.traj_steps)
        
        if val == num_steps:
            # Final/current structure
            self.traj_label.setText("Current (Final)")
            self.traj_info.setText("Showing final optimized structure")
            
            # Show convergence from the LAST trajectory step if it represents the same state
            # (which it usually does in a completed optimization)
            last_conv = {}
            if self.traj_steps:
                last_conv = self.traj_steps[-1].get('convergence', {})
            
            if last_conv:
                self._update_conv_label(last_conv)
            else:
                self.traj_conv.setText("")
            
            # Use current gradients
            self.gradients = self.parser.data.get("gradients", []) if self.parser else []
            
            # Auto-Scale if requested
            if hasattr(self, 'chk_auto_scale') and self.chk_auto_scale.isChecked():
                self.auto_scale()
                
            # Populate table
            self.populate_force_table()

            # Update structure first
            if self.parser and "atoms" in self.parser.data and "coords" in self.parser.data:
                self.update_structure(self.parser.data["atoms"], self.parser.data["coords"])
            
            # Then update vectors with a slight delay to ensure structure is drawn
            QTimer.singleShot(100, self.update_vectors)
        else:
            # Historical step
            step = self.traj_steps[val]
            energy = step.get('energy', 0.0)
            self.traj_label.setText(f"Step {val + 1}/{num_steps}")
            self.traj_info.setText(f"Energy: {energy:.8f} Eh")
            
            # Update gradients for this step
            self.gradients = step.get('gradients', [])
            
            # Auto-Scale if requested
            if hasattr(self, 'chk_auto_scale') and self.chk_auto_scale.isChecked():
                self.auto_scale()
                
            # Populate table
            self.populate_force_table()
            
            # Update structure first
            atoms = step.get('atoms', [])
            coords = step.get('coords', [])
            if atoms and coords:
                self.update_structure(atoms, coords)

            # Then update vectors with delay
            QTimer.singleShot(100, self.update_vectors)

    def _update_conv_label(self, conv):
        """Helper to update the convergence info label with rich text and normalization"""
        if not conv:
            self.traj_conv.setText("No convergence data available")
            return

        # Map specific keys to display names (keys are now lowercased from parser)
        display_keys = {
            "rms gradient": "RMS Grad",
            "max gradient": "MAX Grad",
            "rms step": "RMS Step",
            "max step": "MAX Step",
            "energy change": "Energy Change"
        }
        
        items = []
        for k, v in conv.items():
            if isinstance(v, dict):
                dn = display_keys.get(k.lower(), k.title())
                status = v.get('converged', '??').upper()
                # Status determines color of the whole line
                color = "#28a745" if status == "YES" else "#dc3545"
                items.append(f"<span style='color: {color};'><b>{dn}:</b> {v.get('value')} ({status})</span>")
            else:
                items.append(f"<b>{k}:</b> {v:.6f}")
        
        self.traj_conv.setText("<br>".join(items))
    
    def reload_data(self):
        """Reload data from the output file"""
        if not self.parser or not self.parser.filename:
            QMessageBox.warning(self, "Error", "No file associated with this parser.")
            return
            
        if not os.path.exists(self.parser.filename):
            QMessageBox.warning(self, "Error", f"File not found: {self.parser.filename}")
            return
            
        try:
            with open(self.parser.filename, 'r', encoding='utf-8', errors='replace') as f:
                content = f.read()
            
            # Re-parse (using the already initialized parser)
            self.parser.load_from_memory(content, self.parser.filename)
            
            # Update local state
            self.gradients = self.parser.data.get("gradients", [])
            self.traj_steps = self.parser.data.get("scan_steps", [])
            
            # Update slider if it exists
            if hasattr(self, 'traj_slider'):
                self.traj_slider.setRange(-1, len(self.traj_steps) - 1)
                # Ensure current_step_idx is still valid
                if self.current_step_idx >= len(self.traj_steps):
                    self.current_step_idx = -1
                
                # Update label/text
                self.on_trajectory_change(self.current_step_idx)
            
            # Update force table and vectors
            self.populate_force_table()
            # self.auto_scale() # Use button only
            if self.btn_visualize.isChecked():
                self.update_vectors()
            
            # print(f"Force Viewer: Reloaded from {os.path.basename(self.parser.filename)}")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to reload data: {e}")
            
    def update_structure(self, atoms, coords):
        """Update the 3D structure visualization using RDKit"""
        try:
            from rdkit import Chem
            from rdkit.Geometry import Point3D
            from rdkit.Chem import rdDetermineBonds
        except ImportError:
            return
        
        mol = Chem.RWMol()
        conf = Chem.Conformer()
        pt = Chem.GetPeriodicTable()
        
        try:
            for i, sym in enumerate(atoms):
                if ":" in sym:
                    sym = sym.split(":")[0]
                an = pt.GetAtomicNumber(sym)
                mol.AddAtom(Chem.Atom(an))
                conf.SetAtomPosition(i, Point3D(coords[i][0], coords[i][1], coords[i][2]))
        except:
            return
        
        mol.AddConformer(conf)
        
        # Determine bonds
        if rdDetermineBonds:
            try:
                rdDetermineBonds.DetermineConnectivity(mol)
                charge = self.parser.data.get('charge', 0) if self.parser else 0
                rdDetermineBonds.DetermineBondOrders(mol, charge=charge)
            except:
                pass
        
        final_mol = mol.GetMol()
        
        # Update in main window
        mw = None
        if hasattr(self.parent_dlg, 'context') and self.parent_dlg.context:
            mw = self.parent_dlg.context.get_main_window()
        elif hasattr(self.parent_dlg, 'mw'):
            mw = self.parent_dlg.mw
        
        if mw and hasattr(mw, 'draw_molecule_3d'):
            mw.draw_molecule_3d(final_mol)
    
    def populate_force_table(self):
        """Populate the force and gradient table from current gradient data"""
        self.force_table.setRowCount(0)
        
        # Get atoms from parser
        atoms = []
        if self.parser and "atoms" in self.parser.data:
            atoms = self.parser.data["atoms"]
        
        if not atoms or not self.gradients:
            # Show "No Data" message if gradients are missing but we are supposed to show something
            self.force_table.setRowCount(1)
            msg_item = QTableWidgetItem("No gradient data available for this step in output file")
            msg_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            self.force_table.setItem(0, 0, msg_item)
            self.force_table.setSpan(0, 0, 1, 8)
            return
        
        # Populate with gradient data (force = -gradient)
        for grad_item in self.gradients:
            atom_idx = grad_item.get('atom_idx')
            if atom_idx is None or atom_idx >= len(atoms):
                continue
            
            vec = grad_item.get('grad', grad_item.get('vector'))
            if not vec or len(vec) != 3:
                continue
            
            # Gradient components
            gx, gy, gz = vec[0], vec[1], vec[2]
            
            # Force is negative gradient
            fx, fy, fz = -gx, -gy, -gz
            magnitude = np.sqrt(fx**2 + fy**2 + fz**2)
            
            atom_label = f"{atoms[atom_idx]}{atom_idx + 1}"  # 1-based
            
            row = self.force_table.rowCount()
            self.force_table.insertRow(row)
            self.force_table.setItem(row, 0, QTableWidgetItem(atom_label))
            self.force_table.setItem(row, 1, QTableWidgetItem(f"{gx:.8f}"))
            self.force_table.setItem(row, 2, QTableWidgetItem(f"{gy:.8f}"))
            self.force_table.setItem(row, 3, QTableWidgetItem(f"{gz:.8f}"))
            self.force_table.setItem(row, 4, QTableWidgetItem(f"{fx:.8f}"))
            self.force_table.setItem(row, 5, QTableWidgetItem(f"{fy:.8f}"))
            self.force_table.setItem(row, 6, QTableWidgetItem(f"{fz:.8f}"))
            self.force_table.setItem(row, 7, QTableWidgetItem(f"{magnitude:.8f}"))
        
        self.force_table.resizeColumnsToContents()
    
    def update_vectors(self):
        """Update the force vectors in the 3D visualizer"""
        try:
            mw = None
            if hasattr(self.parent_dlg, 'context') and self.parent_dlg.context:
                mw = self.parent_dlg.context.get_main_window()
            elif hasattr(self.parent_dlg, 'mw'):
                mw = self.parent_dlg.mw
            
            if not mw or not hasattr(mw, 'plotter'):
                return
            
            # Clear old vectors
            self.clear_vectors()
            
            # Only show vectors if checked (visualizing)
            if not self.btn_visualize.isChecked():
                return
            
            # Get coordinates for the current view
            current_coords = []
            if self.current_step_idx == len(self.traj_steps):
                # Final structure
                if self.parser and "coords" in self.parser.data:
                    current_coords = self.parser.data["coords"]
            else:
                # Trajectory step
                current_coords = self.traj_steps[self.current_step_idx].get("coords", [])
            
            # Get settings
            scale = self.spin_scale.value()
            
            if not current_coords or not self.gradients:
                return
            
            # Draw force vectors
            for grad_item in self.gradients:
                atom_idx = grad_item.get('atom_idx')
                if atom_idx is None or atom_idx >= len(current_coords):
                    continue
                
                vec = grad_item.get('grad', grad_item.get('vector'))
                if not vec or len(vec) != 3:
                    continue
                
                # Force = -Gradient
                force = np.array([-vec[0], -vec[1], -vec[2]])
                magnitude = np.linalg.norm(force)
                
                if magnitude < 1e-12:
                    continue
                
                # Scale the vector
                scaled_force = force * scale
                scaled_mag = np.linalg.norm(scaled_force)
                
                if scaled_mag < 1e-12:
                    continue
                
                direction = scaled_force / scaled_mag
                
                # Create arrow
                arrow = pv.Arrow(
                    start=current_coords[atom_idx],
                    direction=direction,
                    scale=scaled_mag,
                    shaft_resolution=20,
                    tip_resolution=20
                )
                
                actor = mw.plotter.add_mesh(arrow, color=self.force_color, name=f'force_{atom_idx}')
                self.actors.append(actor)
            
            mw.plotter.render()
            
        except Exception as e:
            print(f"Error drawing force vectors: {e}")
    
    def clear_vectors(self):
        """Clear all force vector actors"""
        mw = None
        if hasattr(self.parent_dlg, 'context') and self.parent_dlg.context:
            mw = self.parent_dlg.context.get_main_window()
        elif hasattr(self.parent_dlg, 'mw'):
            mw = self.parent_dlg.mw
        
        if not mw or not hasattr(mw, 'plotter'):
            return
        
        for actor in self.actors:
            try:
                mw.plotter.remove_actor(actor)
            except:
                pass
        
        self.actors = []
        mw.plotter.render()
    
    def closeEvent(self, event):
        """Clean up when dialog closes"""
        self.clear_vectors()
        super().closeEvent(event)
