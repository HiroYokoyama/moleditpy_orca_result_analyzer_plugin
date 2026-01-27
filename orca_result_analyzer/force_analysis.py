
import numpy as np
import pyvista as pv
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, 
                             QPushButton, QDoubleSpinBox, QGroupBox, QSpinBox, QColorDialog)
from PyQt6.QtGui import QColor

class ForceViewerDialog(QDialog):
    def __init__(self, parent_dlg, gradients):
        super().__init__(parent_dlg)
        self.setWindowTitle("Force Viewer")
        self.resize(350, 450)
        self.parent_dlg = parent_dlg
        self.gradients = gradients # List of {atom_idx, atom_sym, vector}
        self.actors = []
        self.force_color = "red"
        self.force_res = 20
        
        main_layout = QVBoxLayout(self)
        
        # 1. Info Section
        info_group = QGroupBox("Gradient Info")
        info_layout = QVBoxLayout(info_group)
        info_layout.addWidget(QLabel(f"Visualizing forces for <b>{len(gradients)}</b> atoms."))
        info_layout.addWidget(QLabel("<i>Note: Force = -Gradient</i>"))
        main_layout.addWidget(info_group)
        
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
        
        main_layout.addWidget(view_group)
        main_layout.addStretch()
        
        # 3. Actions
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        main_layout.addWidget(btn_close)
        
        # Init
        self.update_vectors()
        
    def update_vectors(self):
        # Clear old
        self.clear_vectors()
        
        # Access main window via context if available, or direct parent if not?
        # OrcaResultAnalyzerDialog uses context to get main_window().
        
        mw = None
        if hasattr(self.parent_dlg, 'context') and self.parent_dlg.context:
             mw = self.parent_dlg.context.get_main_window()
        elif hasattr(self.parent_dlg, 'mw'):
             mw = self.parent_dlg.mw
             
        if not mw or not hasattr(mw, 'plotter'): return
        
        try:
            # Helper to get atom pos
            mol = getattr(mw, 'current_mol', None)
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
