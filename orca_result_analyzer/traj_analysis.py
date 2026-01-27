
import csv
import numpy as np
import matplotlib
matplotlib.use('QtAgg')
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure

from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, QSlider, 
                             QRadioButton, QComboBox, QPushButton, QSpinBox, 
                             QFormLayout, QDialogButtonBox, QCheckBox, QFileDialog, 
                             QMessageBox, QApplication)
from PyQt6.QtCore import Qt, QTimer

try:
    from rdkit import Chem
    from rdkit.Geometry import Point3D
    from rdkit.Chem import rdDetermineBonds
except ImportError:
    Chem = None
    Point3D = None
    rdDetermineBonds = None

try:
    from PIL import Image
    HAS_PIL = True
except ImportError:
    HAS_PIL = False

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MplCanvas, self).__init__(self.fig)

class TrajectoryResultDialog(QDialog):

    def __init__(self, gl_widget, steps, charge=0, title="Trajectory Analysis"):
        super().__init__()
        self.setWindowTitle(title)
        self.resize(800, 600)
        self.gl_widget = gl_widget
        self.steps = steps # List of {step, energy, atoms, coords}
        self.charge = charge
        self.show_relative = False
        
        # Extract energies
        self.energies = [s['energy'] for s in self.steps]
        # Absolute and Relative energies in current unit
        self.min_e = min(self.energies) if self.energies else 0
        self.current_unit = "kJ/mol" 
        self.display_energies = []
        self.update_display_values()        
        layout = QVBoxLayout(self)
        # 1. Matplotlib Canvas
        self.canvas = MplCanvas(self, width=5, height=4, dpi=100)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        
        # Connect events
        self.canvas.mpl_connect('pick_event', self.on_pick)
        self.canvas.mpl_connect('motion_notify_event', self.on_hover)
        self.canvas.mpl_connect('scroll_event', self.on_scroll)
        
        # Tooltip annotation
        self.annot = self.canvas.axes.annotate("", xy=(0,0), xytext=(20,20),
                            textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="w", alpha=0.9),
                            arrowprops=dict(arrowstyle="->"))
        self.annot.set_visible(False)
        
        self.init_ui()
        self.plot_data()
        
    def init_ui(self):
        
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
        
        self.layout().addLayout(ctrl_layout)
        
        # Toggle Layout
        toggle_layout = QHBoxLayout()
        self.radio_abs = QRadioButton("Absolute")
        self.radio_abs.setChecked(True)
        self.radio_rel = QRadioButton("Relative")
        self.radio_abs.toggled.connect(self.on_toggle_mode)
        
        self.combo_unit = QComboBox()
        self.combo_unit.addItems(["kJ/mol", "kcal/mol", "eV", "Eh"])
        self.combo_unit.currentTextChanged.connect(self.on_unit_changed)
        
        toggle_layout.addWidget(QLabel("Display Mode:"))
        toggle_layout.addWidget(self.radio_abs)
        toggle_layout.addWidget(self.radio_rel)
        toggle_layout.addWidget(QLabel("Unit:"))
        toggle_layout.addWidget(self.combo_unit)
        toggle_layout.addStretch()
        self.layout().addLayout(toggle_layout)
        
        # 3. Buttons (Playback & Export)
        btn_layout = QHBoxLayout()
        
        self.btn_play = QPushButton("Play")
        self.btn_play.clicked.connect(self.toggle_play)
        btn_layout.addWidget(self.btn_play)
        
        btn_layout.addWidget(QLabel("FPS:"))
        self.spin_fps = QSpinBox()
        self.spin_fps.setRange(1, 60)
        self.spin_fps.setValue(5)
        self.spin_fps.valueChanged.connect(self.on_fps_changed)
        btn_layout.addWidget(self.spin_fps)
        
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
        
        btn_clear_sel = QPushButton("Clear Selection")
        # User might mean clear the highlight in graph.
        btn_clear_sel.clicked.connect(self.clear_selection)
        btn_layout.addWidget(btn_clear_sel)
        
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        btn_layout.addWidget(btn_close)
        
        self.layout().addLayout(btn_layout)
        
        # Timer
        self.timer = QTimer()
        self.timer.timeout.connect(self.next_frame)
        self.is_playing = False
        
    def on_toggle_mode(self):
        self.show_relative = self.radio_rel.isChecked()
        self.update_display_values()
        self.plot_data()
        self.on_step_changed(self.slider.value())
        
    def on_unit_changed(self, unit):
        self.current_unit = unit
        self.update_display_values()
        self.plot_data()
        self.on_step_changed(self.slider.value())

    def update_display_values(self):
        if self.current_unit == "kcal/mol":
            factor = 627.509
        elif self.current_unit == "kJ/mol":
            factor = 2625.50
        elif self.current_unit == "eV":
            factor = 27.2114
        else: # Eh
            factor = 1.0
            
        if self.show_relative:
            self.display_energies = [(e - self.min_e) * factor for e in self.energies]
        else:
            self.display_energies = [e * factor for e in self.energies]
        
    def plot_data(self):
        self.canvas.axes.clear()
        
        x = list(range(len(self.energies)))
        y = self.display_energies
        
        if self.show_relative:
            ylabel = f"Relative Energy ({self.current_unit})"
        else:
            ylabel = f"Absolute Energy ({self.current_unit})"
        
        # Draw Line and Scatter
        self.canvas.axes.plot(x, y, 'b-', label='Energy', picker=5)
        self.scatter = self.canvas.axes.scatter(x, y, c='red', s=40, picker=5, zorder=5)
        
        self.canvas.axes.set_xlabel("Step")
        self.canvas.axes.set_ylabel(ylabel)
        self.canvas.axes.set_title("Energy Profile")
        self.canvas.axes.grid(True)
        
        # Re-add tooltip annotation to cleared axis
        self.annot = self.canvas.axes.annotate("", xy=(0,0), xytext=(20,20),
                            textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="w", alpha=0.9),
                            arrowprops=dict(arrowstyle="->"))
        self.annot.set_visible(False)
        
        self.highlight_point(self.slider.value())
        self.canvas.draw()
        
    def highlight_point(self, idx):
        # Remove old markers
        if hasattr(self, '_highlight_marker'):
            try: self._highlight_marker.remove()
            except: pass
        if hasattr(self, '_highlight_line'):
            try: self._highlight_line.remove()
            except: pass
            
        if self.show_relative:
            y = self.display_energies[idx]
        else:
            y = self.display_energies[idx]
            
        x_val = idx # Use idx for x-coordinate
        # Red circle
        self._highlight_marker, = self.canvas.axes.plot(x_val, y, 'ro', markersize=12, markeredgecolor='black', markeredgewidth=2, zorder=10)
        # Vertical Line
        self._highlight_line = self.canvas.axes.axvline(x=x_val, color='gray', linestyle='--', alpha=0.6)
        
        self.canvas.draw()
        
    def on_step_changed(self, idx):
        self.highlight_point(idx)
        step = self.steps[idx]
        val = self.display_energies[idx]
        self.lbl_info.setText(f"Step {idx+1}/{len(self.steps)} | {val:.4f} {self.current_unit}")
        
        self.update_structure(step['atoms'], step['coords'])
        
    def update_structure(self, atoms, coords):
        # RDKit build
        if not Chem: return
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
        
        # Determine Bonds and Bond Orders for each step
        if rdDetermineBonds:
            try:
                rdDetermineBonds.DetermineConnectivity(mol)
                rdDetermineBonds.DetermineBondOrders(mol, charge=self.charge)
            except Exception as e:
                print(f"Bond determination failed at step: {e}")
                
        final_mol = mol.GetMol()
        
        if hasattr(self.gl_widget, 'draw_molecule_3d'):
             self.gl_widget.draw_molecule_3d(final_mol)
             
    def on_scroll(self, event):
        if event.button == 'up':
            self.slider.setValue(min(self.slider.value() + 1, len(self.steps) - 1))
        elif event.button == 'down':
            self.slider.setValue(max(self.slider.value() - 1, 0))

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
                val = self.display_energies[idx]
                text = f"Step {idx}\n{val:.4f} {self.current_unit}"
                self.annot.set_text(text)
                self.annot.set_visible(True)
                self.canvas.draw_idle()
            else:
                if vis:
                    self.annot.set_visible(False)
                    self.canvas.draw_idle()
            
    def toggle_play(self):
        if self.is_playing:
            self.timer.stop()
            self.btn_play.setText("Play")
            self.is_playing = False
        else:
            fps = self.spin_fps.value()
            self.timer.start(int(1000/fps))
            self.btn_play.setText("Pause")
            self.is_playing = True
            
    def on_fps_changed(self):
        if self.is_playing:
            fps = self.spin_fps.value()
            self.timer.start(int(1000/fps))
            
    def next_frame(self):
        idx = self.slider.value() + 1
        if idx >= len(self.steps):
            idx = 0
        self.slider.setValue(idx)
        
    def save_graph(self):
        # Hide annotation before saving
        was_visible = self.annot.get_visible()
        self.annot.set_visible(False)
        self.canvas.draw()
        
        path, _ = QFileDialog.getSaveFileName(self, "Save Graph", "", "Images (*.png *.jpg *.svg)")
        if path:
            self.canvas.fig.savefig(path, dpi=300)
            QMessageBox.information(self, "Saved", f"Graph saved to:\n{path}")
            
        # Restore annotation visibility
        self.annot.set_visible(was_visible)
        self.canvas.draw()

    def clear_selection(self):
        # Remove highlight markers and line
        if hasattr(self, '_highlight_marker'):
            try: self._highlight_marker.remove()
            except: pass
            del self._highlight_marker
        if hasattr(self, '_highlight_line'):
            try: self._highlight_line.remove()
            except: pass
            del self._highlight_line
        
        self.lbl_info.setText("Selection Cleared")
        self.canvas.draw()

    def save_csv(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save CSV", "", "CSV Files (*.csv)")
        if path:
            try:
                with open(path, 'w', newline='', encoding='utf-8') as f:
                    writer = csv.writer(f)
                    writer.writerow(["Step", f"Energy_{self.current_unit.replace('/', '_')}", "Mode"])
                    mode = "Relative" if self.show_relative else "Absolute"
                    for i, step in enumerate(self.steps):
                         writer.writerow([i+1, self.display_energies[i], mode])
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
