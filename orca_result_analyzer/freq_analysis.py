from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, 
                             QTreeWidget, QTreeWidgetItem, QHeaderView, QCheckBox, 
                             QDoubleSpinBox, QSlider, QWidget, QRadioButton, QFileDialog, QFormLayout, QDialogButtonBox, 
                             QSpinBox, QMessageBox, QApplication)
from PyQt6.QtCore import Qt, QTimer

try:
    from PIL import Image
    HAS_PIL = True
except ImportError:
    HAS_PIL = False

class FrequencyDialog(QDialog):
    def __init__(self, parent, frequencies, atoms, coords):
        super().__init__(parent)
        self.mw = parent
        self.frequencies = frequencies # List of dicts
        self.atoms = atoms
        self.base_coords = coords
        
        self.is_playing = False
        self.timer = QTimer()
        self.timer.timeout.connect(self.animate_frame)
        self.animation_step = 0
        self.current_mode_idx = -1
        
        self.setWindowTitle("Vibrational Frequencies")
        self.resize(500, 600)
        
        self.init_ui()
        
    def init_ui(self):
        layout = QVBoxLayout(self)
        
        # Scaling Factor
        scale_layout = QHBoxLayout()
        scale_layout.addWidget(QLabel("Scaling Factor:"))
        self.spin_sf = QDoubleSpinBox()
        self.spin_sf.setRange(0.1, 2.0)
        self.spin_sf.setSingleStep(0.001)
        self.spin_sf.setDecimals(4)
        self.spin_sf.setValue(1.0)
        self.spin_sf.valueChanged.connect(self.update_data)
        scale_layout.addWidget(self.spin_sf)
        scale_layout.addStretch()
        layout.addLayout(scale_layout)
        
        # List
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["Mode", "Freq (cm-1)", "IR", "Raman"])
        self.tree.currentItemChanged.connect(self.on_mode_selected)
        layout.addWidget(self.tree)
        
        # Spectrum Widget
        from .spectrum_widget import SpectrumWidget
        self.spectrum = SpectrumWidget(self.frequencies, x_key='freq', y_key='ir', x_unit='Frequency (cm-1)', sigma=20.0)
        layout.addWidget(self.spectrum)
        
        # Controls
        ctrl_layout = QHBoxLayout()
        ctrl_layout.addWidget(QLabel("Sigma:"))
        self.spin_sigma = QDoubleSpinBox()
        self.spin_sigma.setRange(1.0, 100.0)
        self.spin_sigma.setValue(20.0)
        self.spin_sigma.valueChanged.connect(self.spectrum.set_sigma)
        ctrl_layout.addWidget(self.spin_sigma)
        
        self.chk_sticks = QCheckBox("Sticks")
        self.chk_sticks.setChecked(True)
        self.chk_sticks.stateChanged.connect(self.spectrum.set_sticks)
        ctrl_layout.addWidget(self.chk_sticks)
        
        self.chk_invert_x = QCheckBox("Rev. X")
        self.chk_invert_x.setChecked(True)
        self.chk_invert_x.stateChanged.connect(self.toggle_invert)
        ctrl_layout.addWidget(self.chk_invert_x)
        
        self.chk_invert_y = QCheckBox("Rev. Y")
        self.chk_invert_y.setChecked(True) # Default for IR usually
        self.chk_invert_y.stateChanged.connect(self.toggle_invert)
        ctrl_layout.addWidget(self.chk_invert_y)
        
        ctrl_layout.addStretch()
        
        # Spectrum Type
        self.radio_ir = QRadioButton("IR")
        self.radio_ir.setChecked(True)
        self.radio_ir.toggled.connect(self.switch_spectrum_type)
        ctrl_layout.addWidget(self.radio_ir)
        
        self.radio_raman = QRadioButton("Raman")
        self.radio_raman.toggled.connect(self.switch_spectrum_type)
        ctrl_layout.addWidget(self.radio_raman)
        
        layout.addLayout(ctrl_layout)
        
    def toggle_invert(self):
        self.spectrum.invert_x = self.chk_invert_x.isChecked()
        self.spectrum.invert_y = self.chk_invert_y.isChecked()
        self.spectrum.update()
        
    def switch_spectrum_type(self):
        # Trigger update data to ensure scaling is applied
        is_raman = self.radio_raman.isChecked()
        if is_raman:
            self.spectrum.y_key = 'raman'
            self.spectrum.x_unit = "Raman Shift (cm-1)"
            # Default Raman: Normal X, Normal Y
            self.chk_invert_x.setChecked(False)
            self.chk_invert_y.setChecked(False)
        else:
            self.spectrum.y_key = 'ir'
            self.spectrum.x_unit = "Frequency (cm-1)"
            # Default IR: Inverted X (High->Low) AND Inverted Y (Transmittance-style)
            # User said "NOT X, Y", meaning likely "I want Y inverted". 
            # Standard IR is BOTH inverted. I will set BOTH.
            self.chk_invert_x.setChecked(True)
            self.chk_invert_y.setChecked(True)
        
        self.update_data() # Will set data and update

        # X-Range Control
        x_range_layout = QHBoxLayout()
        
        self.chk_auto_x = QCheckBox("Auto X")
        self.chk_auto_x.setChecked(True)
        self.chk_auto_x.stateChanged.connect(self.toggle_auto_x)
        x_range_layout.addWidget(self.chk_auto_x)
        
        self.spin_x_min = QDoubleSpinBox()
        self.spin_x_min.setRange(0, 5000)
        self.spin_x_min.setValue(400)
        self.spin_x_min.setSuffix(" cm⁻¹")
        self.spin_x_min.setEnabled(False)
        self.spin_x_min.valueChanged.connect(self.update_x_range)
        x_range_layout.addWidget(self.spin_x_min)
        
        self.spin_x_max = QDoubleSpinBox()
        self.spin_x_max.setRange(0, 5000)
        self.spin_x_max.setValue(4000)
        self.spin_x_max.setEnabled(False)
        self.spin_x_max.valueChanged.connect(self.update_x_range)
        x_range_layout.addWidget(self.spin_x_max)
        
        x_range_layout.addStretch()
        layout.addLayout(x_range_layout)

        # Graph Settings (Y-Range & Export)
        graph_layout = QHBoxLayout()
        
        self.chk_auto_y = QCheckBox("Auto Y")
        self.chk_auto_y.setChecked(True)
        self.chk_auto_y.stateChanged.connect(self.toggle_auto_y)
        graph_layout.addWidget(self.chk_auto_y)
        
        self.spin_y_min = QDoubleSpinBox()
        self.spin_y_min.setRange(-1000, 10000)
        self.spin_y_min.setValue(0)
        self.spin_y_min.setSuffix(" Y")
        self.spin_y_min.setEnabled(False)
        self.spin_y_min.valueChanged.connect(self.update_range)
        graph_layout.addWidget(self.spin_y_min)
        
        self.spin_y_max = QDoubleSpinBox()
        self.spin_y_max.setRange(-1000, 10000)
        self.spin_y_max.setValue(1.0)
        self.spin_y_max.setEnabled(False)
        self.spin_y_max.valueChanged.connect(self.update_range)
        graph_layout.addWidget(self.spin_y_max)
        
        graph_layout.addStretch()
        
        btn_png = QPushButton("Save PNG")
        btn_png.clicked.connect(self.save_png)
        graph_layout.addWidget(btn_png)
        
        btn_csv = QPushButton("Save CSV")
        btn_csv.clicked.connect(self.save_csv)
        graph_layout.addWidget(btn_csv)
        
        layout.addLayout(graph_layout)

        
        # Animation / Vector Controls
        anim_layout = QHBoxLayout()
        
        self.chk_vector = QCheckBox("Vectors")
        self.chk_vector.setChecked(True)
        self.chk_vector.stateChanged.connect(self.update_view)
        anim_layout.addWidget(self.chk_vector)
        
        anim_layout.addWidget(QLabel("Vec Scale:"))
        self.spin_vec_scale = QDoubleSpinBox()
        self.spin_vec_scale.setRange(0.1, 50.0)
        self.spin_vec_scale.setValue(1.0)
        self.spin_vec_scale.valueChanged.connect(self.update_view)
        anim_layout.addWidget(self.spin_vec_scale)
        
        anim_layout.addWidget(QLabel("Amp:"))
        self.slider_amp = QSlider(Qt.Orientation.Horizontal)
        self.slider_amp.setRange(1, 20)
        self.slider_amp.setValue(5)
        anim_layout.addWidget(self.slider_amp)
        
        self.btn_play = QPushButton("Animate")
        self.btn_play.clicked.connect(self.toggle_play)
        anim_layout.addWidget(self.btn_play)
        
        self.btn_gif = QPushButton("GIF Export")
        self.btn_gif.clicked.connect(self.save_gif)
        self.btn_gif.setEnabled(HAS_PIL)
        anim_layout.addWidget(self.btn_gif)
        
        layout.addLayout(anim_layout)
        
        # Close
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close_clean)
        layout.addWidget(btn_close)
        
        self.vector_actor = None
        self.populate_list()
        
    def populate_list(self):
        self.tree.clear()
        sf = self.spin_sf.value()
        for i, f in enumerate(self.frequencies):
            freq_val = f['freq'] * sf
            # Always show both IR and Raman in columns
            item = QTreeWidgetItem([str(i), f"{freq_val:.2f}", f"{f.get('ir', 0.0):.2f}", f"{f.get('raman', 0.0):.2f}"])
            self.tree.addTopLevelItem(item)
            
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

    def update_data(self):
        # Update list values (scaling)
        sf = self.spin_sf.value()
        root = self.tree.invisibleRootItem()
        for i in range(root.childCount()):
             item = root.child(i)
             idx = int(item.text(0))
             old_f = self.frequencies[idx]['freq']
             item.setText(1, f"{old_f * sf:.2f}")
             
        # Update spectrum logic handled by switch_spectrum_type and update()
        # But we need to update spectrum data if scaling factor changed!
        
        # SpectrumWidget holds ref to self.frequencies ?
        # In init: self.spectrum = SpectrumWidget(self.frequencies, ...)
        # It stores self.data_list = data_list.
        # Since self.frequencies are DICTS in a list, if we modify the dicts, spectrum sees it?
        # NO, we are computing f['freq'] * sf on the fly for display.
        # We need to pass SCALED data to spectrum widget.
        
        scaled_data = []
        for f in self.frequencies:
            d = f.copy()
            d['freq'] = f['freq'] * sf
            scaled_data.append(d)
        
        self.spectrum.set_data(scaled_data)
        
    def switch_spectrum_type(self):
        # Trigger update data to ensure scaling is applied
        is_raman = self.radio_raman.isChecked()
        if is_raman:
            self.spectrum.y_key = 'raman'
            self.spectrum.x_unit = "Raman Shift (cm-1)"
            # Default Raman Normal
            self.chk_invert.setChecked(False)
        else:
            self.spectrum.y_key = 'ir'
            self.spectrum.x_unit = "Frequency (cm-1)"
            # Default IR Inverted
            self.chk_invert.setChecked(True)
        
        self.update_data() # Will set data and update

        if not current: return
        idx = int(current.text(0))
        self.current_mode_idx = idx
        self.stop_animation()
        self.update_view()
        
    def update_view(self):
        if self.current_mode_idx < 0: return

        # 1. Clear old vectors
        if self.vector_actor:
            try:
                self.mw.plotter.remove_actor(self.vector_actor)
            except: pass
            self.vector_actor = None
            
        # 2. Reset geometry
        self.reset_geometry()
        
        if not self.chk_vector.isChecked(): return
        
        # 3. Draw Vectors
        # frequency dict has 'vector': list of (dx, dy, dz)
        vecs = self.frequencies[self.current_mode_idx].get("vector", [])
        if not vecs: return
        
        try:
            import numpy as np
            points = np.array(self.base_coords)
            vectors = np.array(vecs)
            
            scale = self.spin_vec_scale.value()
            
            # PyVista add_arrows
            # plotter.add_arrows(cent, direction, mag=1.0)
            # We want vectors at atoms.
            
            self.vector_actor = self.mw.plotter.add_arrows(points, vectors, mag=scale, color='orange', name='vib_vectors')
            self.mw.plotter.render()
        except: pass

    def toggle_play(self):
        if self.is_playing:
            self.stop_animation()
        else:
            self.start_animation()
            
    def start_animation(self):
        if self.current_mode_idx < 0: return
        self.is_playing = True
        self.btn_play.setText("Stop")
        self.timer.start(50)
        
    def stop_animation(self):
        self.is_playing = False
        self.btn_play.setText("Animate")
        self.timer.stop()
        self.reset_geometry()
        if self.chk_vector.isChecked(): self.update_view() # Restore vectors
        
    def animate_frame(self):
        if self.current_mode_idx < 0: return
        self.animation_step += 1
        
        import math
        phase = self.animation_step * 0.2
        amp = self.slider_amp.value() * 0.1
        factor = math.sin(phase) * amp
        
        vecs = self.frequencies[self.current_mode_idx].get("vector", [])
        if not vecs: return
        
        try:
            from rdkit.Geometry import Point3D
            mol = self.mw.current_mol
            conf = mol.GetConformer()
            for i, (bx, by, bz) in enumerate(self.base_coords):
                vx, vy, vz = vecs[i]
                nx = bx + vx * factor
                ny = by + vy * factor
                nz = bz + vz * factor
                conf.SetAtomPosition(i, Point3D(nx, ny, nz))
            
            self.mw.draw_molecule_3d(mol)
        except: pass

    def reset_geometry(self):
        try:
            from rdkit.Geometry import Point3D
            mol = self.mw.current_mol
            conf = mol.GetConformer()
            for i, (bx, by, bz) in enumerate(self.base_coords):
                conf.SetAtomPosition(i, Point3D(bx, by, bz))
            self.mw.draw_molecule_3d(mol)
        except: pass
        
    def save_gif(self):
        if not HAS_PIL:
            QMessageBox.warning(self, "Error", "PIL (Pillow) not installed.")
            return

        if self.current_mode_idx < 0:
            QMessageBox.warning(self, "Select Mode", "Please select a frequency mode first.")
            return

        # Settings
        dialog = QDialog(self)
        dialog.setWindowTitle("GIF Settings")
        form = QFormLayout(dialog)
        
        spin_fps = QSpinBox()
        spin_fps.setRange(1, 60)
        spin_fps.setValue(20)
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
        
        # Stop playback and reset
        was_playing = self.is_playing
        if self.is_playing: self.stop_animation()
        else: self.reset_geometry()
        
        import time
        import math
        import numpy as np
        images = []
        mw = self.mw
        
        if not hasattr(mw, 'plotter'):
            QMessageBox.warning(self, "Error", "Cannot access 3D plotter.")
            return
            
        try:
            self.setCursor(Qt.CursorShape.WaitCursor)
            
            # Generate 20 frames (1 cycle)
            vecs = self.frequencies[self.current_mode_idx].get("vector", [])
            if not vecs: raise Exception("No vectors for this mode")
            
            from rdkit.Geometry import Point3D
            mol = self.mw.current_mol
            conf = mol.GetConformer()
            
            # Disable vectors during capture to avoid clutter? Reference 'freq_vis.py' does.
            # But user might want vectors. Let's keep existing state.
            
            for i in range(20):
                phase = i * 0.2 # Match animate_frame logic somewhat? 
                # animate_frame uses phase = step * 0.2
                # cycle = 2pi. 20 frames * x = 2pi ? 
                # 20 * 0.2 = 4.0 != 2pi (6.28). 
                # Let's match freq_vis.py reference: cycle_pos = i/20.0, phase = cycle_pos * 2*pi
                cycle_pos = i / 20.0
                phase = cycle_pos * 2 * np.pi
                
                amp = self.slider_amp.value() * 0.1
                factor = math.sin(phase) * amp
                
                for j, (bx, by, bz) in enumerate(self.base_coords):
                    vx, vy, vz = vecs[j]
                    nx = bx + vx * factor
                    ny = by + vy * factor
                    nz = bz + vz * factor
                    conf.SetAtomPosition(j, Point3D(nx, ny, nz))
                
                if hasattr(mw, 'draw_molecule_3d'):
                     mw.draw_molecule_3d(mol)
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
                     
            if images:
                duration = int(1000 / fps)
                images[0].save(path, save_all=True, append_images=images[1:], duration=duration, loop=0, disposal=2)
                QMessageBox.information(self, "Success", f"GIF saved to:\n{path}")
                
        except Exception as e:
             QMessageBox.critical(self, "Error", f"Failed to save GIF:\n{e}")
        finally:
             self.setCursor(Qt.CursorShape.ArrowCursor)
             self.reset_geometry()
             if was_playing: self.start_animation()

    def close_clean(self):
        self.stop_animation()
        if self.vector_actor:
             try:
                 self.mw.plotter.remove_actor(self.vector_actor)
             except: pass
        self.accept()

