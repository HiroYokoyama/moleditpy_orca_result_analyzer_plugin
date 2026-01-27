
import numpy as np
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, 
                             QPushButton, QCheckBox, QDoubleSpinBox)

class DipoleDialog(QDialog):
    def __init__(self, parent_dlg, dipole_data):
        super().__init__(parent_dlg)
        self.setWindowTitle("Dipole Moment")
        self.resize(300, 200)
        self.parent_dlg = parent_dlg
        self.dipole_data = dipole_data
        self.arrow_actor = None
        
        layout = QVBoxLayout(self)
        
        # Info
        vec = dipole_data.get('vector', [0.0, 0.0, 0.0])
        mag = dipole_data.get('magnitude', 0.0)
        
        info_text = f"<b>Dipole Moment Vector</b><br>" \
                    f"X: {vec[0]:.4f} Debye<br>" \
                    f"Y: {vec[1]:.4f} Debye<br>" \
                    f"Z: {vec[2]:.4f} Debye<br><br>" \
                    f"<b>Magnitude: {mag:.4f} Debye</b>"
        layout.addWidget(QLabel(info_text))
        
        # Controls
        ctrl_layout = QHBoxLayout()
        
        self.chk_show = QCheckBox("Show in 3D")
        self.chk_show.setChecked(False)
        self.chk_show.stateChanged.connect(self.update_view)
        ctrl_layout.addWidget(self.chk_show)
        
        ctrl_layout.addStretch()
        
        layout.addLayout(ctrl_layout)
        
        scale_layout = QHBoxLayout()
        scale_layout.addWidget(QLabel("Scale:"))
        self.spin_scale = QDoubleSpinBox()
        self.spin_scale.setRange(0.1, 10.0)
        self.spin_scale.setSingleStep(0.1)
        self.spin_scale.setValue(1.0)
        self.spin_scale.valueChanged.connect(self.update_view)
        scale_layout.addWidget(self.spin_scale)
        scale_layout.addStretch()
        layout.addLayout(scale_layout)
        
        # Close
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        layout.addWidget(btn_close)
        
    def update_view(self):
        # Clear old
        if self.arrow_actor:
            try:
                self.parent_dlg.mw.plotter.remove_actor(self.arrow_actor)
            except: pass
            self.arrow_actor = None
            
        if not self.chk_show.isChecked():
            self.parent_dlg.mw.plotter.render()
            return
            
        try:
            mw = self.parent_dlg.mw
            # Center of mass calculation
            coords = self.parent_dlg.parser.data.get("coords", [])
            if coords:
                center = np.mean(coords, axis=0)
            else:
                center = np.array([0.0, 0.0, 0.0])
                
            vec = np.array(self.dipole_data.get('vector', [0.0, 0.0, 0.0]))
            mag = np.linalg.norm(vec)
            
            if mag < 1e-6: return
            
            # Normalize and scale
            scale = self.spin_scale.value()
            # Default length 2.0 angstrom for visibility? Or proportional?
            # User usually wants to see direction.
            # Let's use scale * 2.0 as length
            
            direction = vec / mag
            length = 2.0 * scale
            
            # Add Arrow: add_arrow(start, heading, mag=..., color=...)
            self.arrow_actor = mw.plotter.add_arrow(center, direction, mag=length, color='cyan', name='dipole_vector')
            mw.plotter.render()
            
        except Exception as e:
            print(f"Error drawing dipole: {e}")

    def closeEvent(self, event):
        if self.arrow_actor:
             try:
                 self.parent_dlg.mw.plotter.remove_actor(self.arrow_actor)
                 self.parent_dlg.mw.plotter.render()
             except: pass
        # Clean up reference in parent
        if hasattr(self.parent_dlg, 'dipole_dlg'):
             self.parent_dlg.dipole_dlg = None
        event.accept()
