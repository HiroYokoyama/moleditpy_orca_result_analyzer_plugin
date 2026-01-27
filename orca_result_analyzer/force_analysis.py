
import numpy as np
import pyvista as pv
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, 
                             QPushButton, QDoubleSpinBox)

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
        
        # Color - Currently fixed to red
        self.btn_color = QPushButton("Color: Red")
        self.btn_color.setStyleSheet("background-color: red; color: white;")
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
                    
                    actor = mw.plotter.add_arrow(start, force_dir, mag=force_mag, color='red', resolution=10)
                    self.actors.append(actor)
            
            mw.plotter.render()
            
        except Exception as e:
            print(f"Error drawing vectors: {e}")
            
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
