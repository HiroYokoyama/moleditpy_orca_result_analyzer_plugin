
import numpy as np
import pyvista as pv
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
                             QDoubleSpinBox, QGroupBox, QSpinBox, QColorDialog,
                             QFileDialog, QMessageBox, QTableWidget, QTableWidgetItem,
                             QTabWidget, QWidget)
from PyQt6.QtGui import QColor
import os

class ForceViewerDialog(QDialog):
    def __init__(self, parent_dlg, gradients, parser=None):
        super().__init__(parent_dlg)
        self.setWindowTitle("Force Viewer")
        self.resize(500, 600)
        self.parent_dlg = parent_dlg
        self.gradients = gradients # List of {atom_idx, atom_sym, vector}
        self.parser = parser
        self.actors = []
        self.force_color = "red"
        self.force_res = 20
        self.hessian_data = None
        
        main_layout = QVBoxLayout(self)
        
        # Tabs for Gradients and Hessian
        self.tabs = QTabWidget()
        main_layout.addWidget(self.tabs)
        
        # Tab 1: Gradient Visualization
        # Tab 1: Gradient Visualization
        grad_widget = QWidget()
        grad_layout = QVBoxLayout(grad_widget)
        
        has_gradients = bool(self.gradients and len(self.gradients) > 0)
        
        if has_gradients:
            self._setup_gradient_tab(grad_layout)
            self.tabs.addTab(grad_widget, "Gradients")
        else:
            grad_layout.addWidget(QLabel("No gradient (force) data found in output file.\nThis is normal for successful frequency jobs."))
            grad_layout.addStretch()
            self.tabs.addTab(grad_widget, "Gradients (Empty)")
        
        # Tab 2: Hessian/Force Constants
        hess_widget = QWidget()
        hess_layout = QVBoxLayout(hess_widget)
        self._setup_hessian_tab(hess_layout)
        self.tabs.addTab(hess_widget, "Force Constants")
        
        # Auto-switch to Hessian tab if no gradients
        if not has_gradients:
             self.tabs.setCurrentIndex(1)
             
        # Close button
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        main_layout.addWidget(btn_close)
        
        # Init
        self.update_vectors()
    
    def _setup_gradient_tab(self, layout):
        """Setup the gradient visualization tab"""
        # 1. Info Section
        info_group = QGroupBox("Gradient Info")
        info_layout = QVBoxLayout(info_group)
        info_layout.addWidget(QLabel(f"Visualizing forces for <b>{len(self.gradients)}</b> atoms."))
        info_layout.addWidget(QLabel("<i>Note: Force = -Gradient</i>"))
        layout.addWidget(info_group)
        
        # 2. 3D Appearance Section
        view_group = QGroupBox("3D Visualization")
        view_layout = QVBoxLayout(view_group)
        
        # Scale row
        scale_row = QHBoxLayout()
        scale_row.addWidget(QLabel("Scale:"))
        self.spin_scale = QDoubleSpinBox()
        self.spin_scale.setRange(0.1, 100.0)
        self.spin_scale.setValue(5.0) 
        self.spin_scale.setSingleStep(0.5)
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
    
    def _setup_hessian_tab(self, layout):
        """Setup the Hessian/Force Constants tab"""
        info_lbl = QLabel("Load a Hessian file (.hess) to view force constants and frequencies")
        layout.addWidget(info_lbl)
        
        # Load button
        btn_load = QPushButton("Load Hessian File")
        btn_load.clicked.connect(self.load_hessian_file)
        layout.addWidget(btn_load)
        
        # Frequencies Table
        self.freq_label = QLabel("Frequencies:")
        layout.addWidget(self.freq_label)
        
        self.freq_table = QTableWidget()
        self.freq_table.setColumnCount(2)
        self.freq_table.setHorizontalHeaderLabels(["Mode", "Frequency (cm⁻¹)"])
        self.freq_table.setAlternatingRowColors(True)
        # Limit height
        self.freq_table.setFixedHeight(150)
        layout.addWidget(self.freq_table)
        
        # Force Constants Table
        layout.addWidget(QLabel("Force Constants (Diagonal):"))
        self.hessian_table = QTableWidget()
        self.hessian_table.setColumnCount(3)
        self.hessian_table.setHorizontalHeaderLabels(["Atom Pair", "Force Constant", "Unit"])
        layout.addWidget(self.hessian_table)
        
        layout.addStretch()
    
    def load_hessian_file(self):
        """Load and parse a Hessian file"""
        if not self.parser:
            QMessageBox.warning(self, "No Parser", "Parser not available")
            return
        
        # Get directory from parent
        start_dir = ""
        if hasattr(self.parent_dlg, 'file_path'):
            start_dir = os.path.dirname(self.parent_dlg.file_path)
        
        filepath, _ = QFileDialog.getOpenFileName(
            self, "Open Hessian File", start_dir, "Hessian Files (*.hess);;All Files (*.*)"
        )
        
        if not filepath:
            return
        
        # Parse the file
        self.hessian_data = self.parser.parse_hessian_file(filepath)
        
        if not self.hessian_data or not self.hessian_data.get("matrix"):
            QMessageBox.warning(self, "Parse Error", "Failed to parse Hessian file")
            return
        
        # Display force constants
        self.display_force_constants()
        
        QMessageBox.information(self, "Success", f"Loaded Hessian file successfully!\nFrequencies: {len(self.hessian_data.get('frequencies', []))}")
    
    def display_force_constants(self):
        """Display force constants from Hessian matrix"""
        if not self.hessian_data:
            return
        
        matrix = self.hessian_data.get("matrix", [])
        atoms = self.hessian_data.get("atoms", [])
        
        if not matrix:
            return
        
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
                atom_label = f"{atoms[atom_idx]}{atom_idx}"
                force_const = matrix[i][i] if i < len(matrix[i]) else 0.0
                
                self.hessian_table.insertRow(row)
                self.hessian_table.setItem(row, 0, QTableWidgetItem(f"{atom_label}-{coord_name}"))
                self.hessian_table.setItem(row, 1, QTableWidgetItem(f"{force_const:.6f}"))
                self.hessian_table.setItem(row, 2, QTableWidgetItem("Eh/Bohr²"))
                row += 1
        
        self.hessian_table.resizeColumnsToContents()
        
        # Display frequencies if available
        freqs = self.hessian_data.get("frequencies", [])
        self.freq_table.setRowCount(len(freqs))
        if freqs:
            self.freq_label.setText(f"Frequencies ({len(freqs)} modes):")
            for i, freq in enumerate(freqs):
                self.freq_table.setItem(i, 0, QTableWidgetItem(f"Mode {i}"))
                self.freq_table.setItem(i, 1, QTableWidgetItem(f"{freq:.2f}"))
        else:
             self.freq_label.setText("Frequencies: None found")
        
    def update_vectors(self):
        """Update the force vectors in the visualizer"""
        # Guard clause: if widgets not initialized (e.g. no gradients tab), return
        if not hasattr(self, 'spin_scale'):
            return
            
        try:
            mw = None
            if hasattr(self.parent_dlg, 'context') and self.parent_dlg.context:
                 mw = self.parent_dlg.context.get_main_window()
            elif hasattr(self.parent_dlg, 'mw'):
                 mw = self.parent_dlg.mw
            
            if not mw or not hasattr(mw, 'plotter'): return
            
            # Clear old
            self.clear_vectors()
            
            # Get settings
            scale = self.spin_scale.value()
            
            # Draw
            for idx, item in enumerate(self.gradients):
                # item: {atom_idx, atom_sym, vector: [dx, dy, dz]}
                # We need atom positions. 
                # Parser has data["coords"]. 
                # Let's hope self.parser exists or gradients include positions?
                # Gradients dict in parser usually: {'atom_idx': i, 'sym': s, 'grad': [x,y,z]}
                # We need coords separately.
                atom_idx = item.get('atom_idx', idx)
                
                # Fetch coord from parser
                if self.parser and "coords" in self.parser.data:
                    coords = self.parser.data["coords"]
                    if atom_idx < len(coords):
                        start = coords[atom_idx]
                    else: continue
                else: continue
                
                grad = item.get('grad', item.get('vector'))
                if grad:
                    vec = np.array(grad)
                    length = np.linalg.norm(vec)
                    if length < 1e-6: continue
                    
                    # Gradient is derivative of Energy. Force is -Gradient.
                    # Usually we visualize Force (-Gradient) to show where atoms want to go.
                    # Let's assume user wants Force.
                    force = -vec * scale # Reverse direction
                    
                    # Centered on atom? Yes.
                    # add_arrow(start, direction, mag=length, ...)
                    # direction needs to be normalized? pv docs say direction vector.
                    # If we use `mag`, it scales direction.
                    # Better to use normalized direction and set mag explicitly.
                    force_mag = np.linalg.norm(force)
                    if force_mag < 1e-9: continue
                    
                    force_dir = force / force_mag
                    
                    # Use add_mesh with pv.Arrow for robustness
                    arrow = pv.Arrow(start=start, direction=force_dir, scale=force_mag,
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
