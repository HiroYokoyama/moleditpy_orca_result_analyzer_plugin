from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, 
                             QTreeWidget, QTreeWidgetItem, QHeaderView, QCheckBox, 
                             QDoubleSpinBox, QSlider, QWidget, QRadioButton, QFileDialog, QFormLayout, QDialogButtonBox, 
                             QSpinBox, QMessageBox, QApplication, QColorDialog, QGroupBox)
from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QColor
import time
import math
import numpy as np
import pyvista as pv

try:
    from PIL import Image
    HAS_PIL = True
except ImportError:
    HAS_PIL = False

class FreqSpectrumWindow(QWidget):
    """
    Separate window for displaying the spectrum.
    """
    def __init__(self, parent_dialog, frequencies):
        super().__init__()
        # We don't verify parent strictly to allow independent testing if needed, 
        # but logically it belongs to FrequencyDialog
        self.freq_dialog = parent_dialog 
        self.frequencies = frequencies
        self.scaling_factor = 1.0
        
        self.setWindowTitle("IR/Raman Spectrum")
        self.resize(600, 500)
        
        self.init_ui()
        
    def init_ui(self):
        layout = QVBoxLayout(self)
        
        # Spectrum Widget
        from .spectrum_widget import SpectrumWidget
        self.spectrum = SpectrumWidget(self.frequencies, x_key='freq', y_key='ir', x_unit='Frequency (cm-1)', y_unit='Intensity (a.u.)', sigma=20.0)
        self.spectrum.show_legend = False
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

        # X-Range Control
        x_range_layout = QHBoxLayout()
        
        self.chk_auto_x = QCheckBox("Auto X")
        self.chk_auto_x.setChecked(False)
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
        
        self.toggle_auto_x() # Initialize state (enable spinners, set range)
        
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
        
        # Initial Update
        self.switch_spectrum_type()
        
    def set_scaling_factor(self, sf):
        self.scaling_factor = sf
        self.update_data()
        
    def update_data(self):
        # Re-calc scaled data
        scaled_data = []
        for f in self.frequencies:
            d = f.copy()
            if 'freq' in d:
                d['freq'] = d['freq'] * self.scaling_factor
            scaled_data.append(d)
        self.spectrum.set_data(scaled_data)
        
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
            self.spectrum.y_unit = "Raman Intensity (a.u.)"
            # Default Raman: Normal X, Normal Y
            self.chk_invert_x.setChecked(False)
            self.chk_invert_y.setChecked(False)
        else:
            self.spectrum.y_key = 'ir'
            self.spectrum.x_unit = "Frequency (cm-1)"
            self.spectrum.y_unit = "Intensity (a.u.)"
            # Default IR: Inverted X (High->Low) AND Inverted Y (Transmittance-style)
            self.chk_invert_x.setChecked(True)
            self.chk_invert_y.setChecked(True)
        
        self.toggle_invert() # Force update of spectrum widget properties
        self.update_data() # Will set data and update
            
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
    
    def closeEvent(self, event):
        # Notify parent that we are closing (optional, effectively done by parent checking .isVisible())
        # Or we can nullify reference in parent?
        if self.freq_dialog and hasattr(self.freq_dialog, 'spectrum_win'):
            # Don't nullify, just let it stay? 
            # Better to let parent know we can reopen.
            # But parent just creates new one or shows existing?
            # We'll just hide basically.
            pass
        super().closeEvent(event)


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
        self.vector_color = "orange"
        self.vector_res = 20
        
        self.spectrum_win = None
        
        self.setWindowTitle("Vibrational Frequencies")
        self.resize(450, 650) 
        
        self.init_ui()
        
    def init_ui(self):
        main_layout = QVBoxLayout(self)
        
        # 1. Frequency List Section
        list_group = QGroupBox("Frequency Modes")
        list_layout = QVBoxLayout(list_group)
        
        # Scaling Factor
        sf_layout = QHBoxLayout()
        sf_layout.addWidget(QLabel("Frequency Scaling:"))
        self.spin_sf = QDoubleSpinBox()
        self.spin_sf.setRange(0.1, 2.0)
        self.spin_sf.setSingleStep(0.001)
        self.spin_sf.setDecimals(4)
        self.spin_sf.setValue(1.0)
        self.spin_sf.valueChanged.connect(self.update_data)
        sf_layout.addWidget(self.spin_sf)
        sf_layout.addStretch()
        list_layout.addLayout(sf_layout)
        
        # List
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["Mode", "Freq (cm⁻¹)", "IR", "Raman"])
        self.tree.currentItemChanged.connect(self.on_mode_selected)
        list_layout.addWidget(self.tree)
        
        # Spectrum Button
        btn_spectrum = QPushButton("Open IR/Raman Spectrum...")
        btn_spectrum.clicked.connect(self.open_spectrum)
        list_layout.addWidget(btn_spectrum)
        
        main_layout.addWidget(list_group)
        
        # 2. 3D Appearance (Vectors)
        vec_group = QGroupBox("3D Vector Appearance")
        vec_layout = QVBoxLayout(vec_group)
        
        vec_top_row = QHBoxLayout()
        self.chk_vector = QCheckBox("Show Vectors")
        self.chk_vector.setChecked(True)
        self.chk_vector.stateChanged.connect(self.update_view)
        vec_top_row.addWidget(self.chk_vector)
        
        vec_top_row.addWidget(QLabel("Scale:"))
        self.spin_vec_scale = QDoubleSpinBox()
        self.spin_vec_scale.setRange(0.1, 50.0)
        self.spin_vec_scale.setValue(2.0)
        self.spin_vec_scale.valueChanged.connect(self.update_view)
        vec_top_row.addWidget(self.spin_vec_scale)
        vec_layout.addLayout(vec_top_row)
        
        vec_bot_row = QHBoxLayout()
        vec_bot_row.addWidget(QLabel("Resolution:"))
        self.spin_vec_res = QSpinBox()
        self.spin_vec_res.setRange(3, 100)
        self.spin_vec_res.setValue(20)
        self.spin_vec_res.valueChanged.connect(self.on_res_changed)
        vec_bot_row.addWidget(self.spin_vec_res)
        
        vec_bot_row.addWidget(QLabel(" Color:"))
        self.btn_vec_color = QPushButton()
        self.btn_vec_color.setFixedWidth(60)
        self.btn_vec_color.setStyleSheet(f"background-color: {self.vector_color}; border: 1px solid gray; height: 20px;")
        self.btn_vec_color.clicked.connect(self.pick_color)
        vec_bot_row.addWidget(self.btn_vec_color)
        vec_bot_row.addStretch()
        vec_layout.addLayout(vec_bot_row)
        
        main_layout.addWidget(vec_group)
        
        # 3. Animation Section
        anim_group = QGroupBox("Animation Parameters")
        anim_layout = QVBoxLayout(anim_group)
        
        anim_row1 = QHBoxLayout()
        anim_row1.addWidget(QLabel("Displacement Scale:"))
        self.spin_amp = QDoubleSpinBox()
        self.spin_amp.setRange(0.1, 10.0)
        self.spin_amp.setSingleStep(0.1)
        self.spin_amp.setValue(1.0)
        anim_row1.addWidget(self.spin_amp)
        
        anim_row1.addWidget(QLabel(" | FPS:"))
        self.spin_fps = QSpinBox()
        self.spin_fps.setRange(1, 120)
        self.spin_fps.setValue(30)
        self.spin_fps.valueChanged.connect(self.update_fps)
        anim_row1.addWidget(self.spin_fps)
        anim_layout.addLayout(anim_row1)
        
        # Playback row
        action_row = QHBoxLayout()
        self.btn_play = QPushButton("Play")
        self.btn_play.clicked.connect(self.start_animation)
        action_row.addWidget(self.btn_play)
        
        self.btn_pause = QPushButton("Pause")
        self.btn_pause.clicked.connect(self.pause_animation)
        self.btn_pause.setEnabled(False)
        action_row.addWidget(self.btn_pause)
        
        self.btn_stop = QPushButton("Stop")
        self.btn_stop.clicked.connect(self.stop_animation)
        self.btn_stop.setEnabled(False)
        action_row.addWidget(self.btn_stop)
        
        action_row.addStretch()
        
        self.btn_gif = QPushButton("GIF")
        self.btn_gif.clicked.connect(self.save_gif)
        self.btn_gif.setEnabled(HAS_PIL)
        action_row.addWidget(self.btn_gif)
        
        anim_layout.addLayout(action_row)
        main_layout.addWidget(anim_group)
        
        # Close
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close_clean)
        main_layout.addWidget(btn_close)
        
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

    def update_data(self):
        # Update list values (scaling)
        sf = self.spin_sf.value()
        root = self.tree.invisibleRootItem()
        for i in range(root.childCount()):
             item = root.child(i)
             idx = int(item.text(0))
             old_f = self.frequencies[idx]['freq']
             item.setText(1, f"{old_f * sf:.2f}")
             
        # Update spectrum window if it's open
        if self.spectrum_win is not None:
             self.spectrum_win.set_scaling_factor(sf)

    def open_spectrum(self):
        if self.spectrum_win is None:
            self.spectrum_win = FreqSpectrumWindow(self, self.frequencies)
            # Apply current scaling factor
            self.spectrum_win.set_scaling_factor(self.spin_sf.value())
            self.spectrum_win.show()
        else:
            self.spectrum_win.show()
            self.spectrum_win.activateWindow()
            self.spectrum_win.raise_()

    def on_mode_selected(self, current, previous):
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
            
            # Use add_mesh with glyphs for robustness across plotter implementations
            poly = pv.PolyData(points)
            poly.point_data['vectors'] = vectors
            geom = pv.Arrow(shaft_resolution=self.vector_res, tip_resolution=self.vector_res)
            arrows = poly.glyph(orient=True, scale=True, factor=scale, geom=geom)
            self.vector_actor = self.mw.plotter.add_mesh(arrows, color=self.vector_color, name='vib_vectors')
            self.mw.plotter.render()
        except Exception as e:
            print(f"Error in FrequencyDialog.update_view: {e}")

    def start_animation(self):
        if self.current_mode_idx < 0: return
        if self.is_playing: return
        
        self.is_playing = True
        self.btn_play.setEnabled(False)
        self.btn_pause.setEnabled(True)
        self.btn_stop.setEnabled(True)
        
        fps = self.spin_fps.value()
        self.timer.start(int(1000/fps))
        
    def pause_animation(self):
        self.is_playing = False
        self.timer.stop()
        self.btn_play.setEnabled(True)
        self.btn_pause.setEnabled(False)
        self.btn_stop.setEnabled(True) # Can still stop to reset
        
    def stop_animation(self):
        self.is_playing = False
        self.timer.stop()
        self.animation_step = 0 # Reset phase
        self.reset_geometry()
        if self.chk_vector.isChecked(): self.update_view() # Restore vectors
        
        self.btn_play.setEnabled(True)
        self.btn_pause.setEnabled(False)
        self.btn_stop.setEnabled(False)
        
    def update_fps(self):
        if self.is_playing:
            self.timer.setInterval(int(1000/self.spin_fps.value()))

    def animate_frame(self):
        if self.current_mode_idx < 0: return
        
        self.animation_step += 1
        phase = self.animation_step * 0.2
        
        amp = self.spin_amp.value()
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
            
            # Update vectors if enabled
            if self.chk_vector.isChecked():
                self.update_vectors_at_displaced_position()
                
        except Exception as e:
            print(f"Error in animate_frame: {e}")
            import traceback
            traceback.print_exc()
    
    def update_vectors_at_displaced_position(self):
        """Redraw vectors at current displaced atomic positions"""
        if self.vector_actor:
            try:
                self.mw.plotter.remove_actor(self.vector_actor)
            except: pass
            self.vector_actor = None
        
        vecs = self.frequencies[self.current_mode_idx].get("vector", [])
        if not vecs or not hasattr(self.mw, 'current_mol'):
            return
            
        try:
            mol = self.mw.current_mol
            conf = mol.GetConformer()
            
            # Get current (displaced) coordinates
            coords = []
            for i in range(conf.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                coords.append([pos.x, pos.y, pos.z])
            
            coords_array = np.array(coords)
            vecs_array = np.array(vecs)
            scale = self.spin_vec_scale.value()
            
            # Use add_arrows if available, otherwise fall back to glyph method
            if hasattr(self.mw.plotter, 'add_arrows'):
                self.vector_actor = self.mw.plotter.add_arrows(
                    coords_array, vecs_array, mag=scale, color=self.vector_color, show_scalar_bar=False
                )
            else:
                # Fallback to glyph method
                poly = pv.PolyData(coords_array)
                poly.point_data['vectors'] = vecs_array
                geom = pv.Arrow(shaft_resolution=self.vector_res, tip_resolution=self.vector_res)
                arrows = poly.glyph(orient=True, scale=True, factor=scale, geom=geom)
                self.vector_actor = self.mw.plotter.add_mesh(arrows, color=self.vector_color, name='vib_vectors')
        except Exception as e:
            print(f"Error updating vectors: {e}")

    def reset_geometry(self):
        try:
            from rdkit.Geometry import Point3D
            mol = self.mw.current_mol
            conf = mol.GetConformer()
            for i, (bx, by, bz) in enumerate(self.base_coords):
                conf.SetAtomPosition(i, Point3D(bx, by, bz))
            self.mw.draw_molecule_3d(mol)
        except Exception as e:
            print(f"Error in reset_geometry: {e}")
            import traceback
            traceback.print_exc()
        
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
        
        chk_hq = QCheckBox()
        chk_hq.setChecked(True)
        form.addRow("High Quality (Adaptive):", chk_hq)
        
        btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        btns.accepted.connect(dialog.accept)
        btns.rejected.connect(dialog.reject)
        form.addRow(btns)
        
        if dialog.exec() != QDialog.DialogCode.Accepted: return
        
        fps = spin_fps.value()
        transparent = chk_trans.isChecked()
        use_hq = chk_hq.isChecked()
        
        path, _ = QFileDialog.getSaveFileName(self, "Save GIF", "", "GIF Files (*.gif)")
        if not path: return
        if not path.lower().endswith('.gif'): path += '.gif'
        
        # Stop playback and reset
        was_playing = self.is_playing
        if self.is_playing: self.stop_animation()
        else: self.reset_geometry()
        
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
            
            for i in range(20):
                cycle_pos = i / 20.0
                phase = cycle_pos * 2 * np.pi
                
                amp = self.spin_amp.value()
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
                     images.append(img)
                     
            if images:
                duration = int(1000 / fps)
                processed_images = []
                for img in images:
                    if use_hq:
                        if transparent:
                            # Alpha preservation with adaptive palette
                            alpha = img.split()[3]
                            img_rgb = img.convert("RGB")
                            # Quantize to 255 colors to leave room for transparency
                            img_p = img_rgb.convert('P', palette=Image.Palette.ADAPTIVE, colors=255)
                            # Set transparency
                            mask = Image.eval(alpha, lambda a: 255 if a <= 128 else 0)
                            img_p.paste(255, mask)
                            img_p.info['transparency'] = 255
                            processed_images.append(img_p)
                        else:
                            processed_images.append(img.convert("P", palette=Image.Palette.ADAPTIVE, colors=256))
                    else:
                        if transparent:
                            processed_images.append(img.convert("RGBA"))
                        else:
                            processed_images.append(img.convert("RGB"))
                
                processed_images[0].save(path, save_all=True, append_images=processed_images[1:], duration=duration, loop=0, disposal=2)
                QMessageBox.information(self, "Success", f"GIF saved to:\n{path}")
                
        except Exception as e:
             QMessageBox.critical(self, "Error", f"Failed to save GIF:\n{e}")
        finally:
             self.setCursor(Qt.CursorShape.ArrowCursor)
             self.reset_geometry()
             if was_playing: self.start_animation()

    def pick_color(self):
        color = QColorDialog.getColor(QColor(self.vector_color), self, "Select Vector Color")
        if color.isValid():
            self.vector_color = color.name()
            self.btn_vec_color.setStyleSheet(f"background-color: {self.vector_color};")
            self.update_view()

    def on_res_changed(self, val):
        self.vector_res = val
        self.update_view()

    def _dummy_import_pv(self):
        import pyvista as pv # Ensure pv is available in local scope if needed, though usually global
              
    def closeEvent(self, event):
        self.stop_animation()
        if self.vector_actor:
             try:
                 self.mw.plotter.remove_actor(self.vector_actor)
             except: pass
        if self.spectrum_win:
             self.spectrum_win.close()
        event.accept()

    def close_clean(self):
        self.close()
