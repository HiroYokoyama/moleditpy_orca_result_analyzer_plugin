
import os
import json
import numpy as np
import pyvista as pv
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, 
                             QComboBox, QTableWidget, QTableWidgetItem, QHeaderView, 
                             QWidget, QCheckBox, QInputDialog, QColorDialog, QMessageBox)
from PyQt6.QtGui import QColor, QPainter, QLinearGradient
from PyQt6.QtCore import Qt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap

# GradientBar Widget
class GradientBar(QWidget):
    def __init__(self, parent=None, colors=["red", "white", "blue"]):
        super().__init__(parent)
        self.colors = colors
        self.setFixedHeight(30)
        
    def set_colors(self, colors):
        self.colors = colors
        self.update()
        
    def paintEvent(self, event):
        painter = QPainter(self)
        grad =  self.get_gradient()
        painter.fillRect(self.rect(), grad)
        
        # Draw border
        painter.setPen(Qt.GlobalColor.black)
        painter.drawRect(0, 0, self.width()-1, self.height()-1)
        
    def get_gradient(self):
        grad = QLinearGradient(0, 0, self.width(), 0)
        
        n = len(self.colors)
        if n < 2: return grad
        
        for i, c in enumerate(self.colors):
            pos = i / (n - 1)
            grad.setColorAt(pos, QColor(c))
        return grad

class ChargeDialog(QDialog):
    def __init__(self, parent, all_charges):
        super().__init__(parent)
        self.parent_dlg = parent # OrcaResultAnalyzerDialog
        self.setWindowTitle("Atomic Charges")
        self.resize(500, 600)
        self.all_charges = all_charges 
        self.current_type = next(iter(self.all_charges)) if self.all_charges else ""
        
        # Color Schemes - default schemes
        self.schemes = {
            "Red(-) - White - Blue(+)": ["red", "white", "blue"],
            "Blue(-) - White - Red(+)": ["blue", "white", "red"],
            "Red(-) - Blue(+)": ["red", "blue"],
            "Green(-) - White - Purple(+)": ["green", "white", "purple"]
        }
        
        # Load custom schemes from settings.json
        settings_file = os.path.join(os.path.dirname(__file__), "settings.json")
        self.current_scheme = "Red(-) - White - Blue(+)"
        
        if os.path.exists(settings_file):
            try:
                with open(settings_file, 'r') as f:
                    settings_data = json.load(f)
                
                # Load custom schemes
                if "custom_color_schemes" in settings_data:
                    for scheme_data in settings_data["custom_color_schemes"]:
                        name = scheme_data.get("name", "")
                        colors = scheme_data.get("colors", [])
                        if name and colors:
                            self.schemes[f"Custom: {name}"] = colors
                
                # Load last used scheme
                if "last_charge_scheme" in settings_data:
                    self.current_scheme = settings_data["last_charge_scheme"]
            except Exception as e:
                print(f"Error loading settings: {e}")
        
        layout = QVBoxLayout(self)
        
        # Type Selector
        head_layout = QHBoxLayout()
        head_layout.addWidget(QLabel("Charge Type:"))
        self.combo_type = QComboBox()
        self.combo_type.addItems(list(self.all_charges.keys()))
        self.combo_type.currentTextChanged.connect(self.on_type_change)
        head_layout.addWidget(self.combo_type)
        layout.addLayout(head_layout)
        
        # Table
        self.table = QTableWidget()
        self.table.setColumnCount(3)
        self.table.setHorizontalHeaderLabels(["Idx", "Atom", "Charge"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.table)
        
        # Coloring Controls
        grp_color = QWidget()
        grp_layout = QVBoxLayout(grp_color)
        
        # Scheme Selector
        scheme_layout = QHBoxLayout()
        scheme_layout.addWidget(QLabel("Color Scheme:"))
        self.combo_scheme = QComboBox()
        self.combo_scheme.addItems(list(self.schemes.keys()))
        self.combo_scheme.setCurrentText(self.current_scheme)
        self.combo_scheme.currentTextChanged.connect(self.on_scheme_change)
        scheme_layout.addWidget(self.combo_scheme)
        
        # Custom scheme button
        btn_custom = QPushButton("+ Custom")
        btn_custom.clicked.connect(self.edit_custom_scheme)
        scheme_layout.addWidget(btn_custom)
        
        grp_layout.addLayout(scheme_layout)
        
        # Gradient Bar
        self.grad_bar = GradientBar(self, self.schemes[self.current_scheme])
        grp_layout.addWidget(self.grad_bar)
        
        # Labels for gradient
        lbl_layout = QHBoxLayout()
        lbl_layout.addWidget(QLabel("Negative"))
        lbl_mid = QLabel("Neutral")
        lbl_mid.setAlignment(Qt.AlignmentFlag.AlignCenter)
        lbl_layout.addWidget(lbl_mid)
        lbl_max = QLabel("Positive")
        lbl_max.setAlignment(Qt.AlignmentFlag.AlignRight)
        lbl_layout.addWidget(lbl_max)
        grp_layout.addLayout(lbl_layout)
        
        # Show labels checkbox
        self.chk_show_labels = QCheckBox("Show charge labels in 3D")
        self.chk_show_labels.setChecked(False)
        self.chk_show_labels.stateChanged.connect(self.toggle_labels)
        grp_layout.addWidget(self.chk_show_labels)
        
        btn_colorize = QPushButton("Colorize Atoms in 3D View")
        btn_colorize.setStyleSheet("font-weight: bold; background-color: #2196F3; color: white;")
        btn_colorize.clicked.connect(self.apply_colors)
        grp_layout.addWidget(btn_colorize)
        
        # Reset color button
        btn_reset = QPushButton("Reset Colors")
        btn_reset.setStyleSheet("background-color: #f44336; color: white;")
        btn_reset.clicked.connect(self.reset_colors)
        grp_layout.addWidget(btn_reset)
        
        layout.addWidget(grp_color)
        
        # Bottom Buttons
        bottom_layout = QHBoxLayout()
        
        btn_clear = QPushButton("Clear Selection")
        btn_clear.clicked.connect(self.table.clearSelection)
        bottom_layout.addWidget(btn_clear)
        
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        bottom_layout.addWidget(btn_close)
        
        layout.addLayout(bottom_layout)
        
        self.update_table()
        
    def on_type_change(self, text):
        self.current_type = text
        self.update_table()
    
    def on_scheme_change(self, text):
        self.current_scheme = text
        colors = self.schemes.get(text, ["red", "white", "blue"])
        self.grad_bar.set_colors(colors)
        
        # Save last used scheme to settings.json
        self.save_settings()
    
    def edit_custom_scheme(self):
        """Create a custom color scheme"""
        
        # Get scheme name
        name, ok = QInputDialog.getText(self, "Custom Scheme", "Enter scheme name:")
        if not ok or not name:
            return
        
        # Get number of colors
        num_colors, ok = QInputDialog.getInt(self, "Custom Scheme", "Number of colors (2-5):", 3, 2, 5)
        if not ok:
            return
        
        # Collect colors
        colors = []
        for i in range(num_colors):
            label = ["Negative", "Mid-Negative", "Neutral", "Mid-Positive", "Positive"][i] if num_colors <= 5 else f"Color {i+1}"
            color = QColorDialog.getColor(QColor("white"), self, f"Select {label} color")
            if not color.isValid():
                return
            colors.append(color.name())
        
        # Add to schemes
        scheme_name = f"Custom: {name}"
        self.schemes[scheme_name] = colors
        
        # Update combo box
        self.combo_scheme.addItem(scheme_name)
        self.combo_scheme.setCurrentText(scheme_name)
        
        # Save to settings
        self.save_custom_schemes()
    
    def save_custom_schemes(self):
        """Save custom schemes to settings.json"""
        self.save_settings()
    
    def save_settings(self):
        """Save all settings to settings.json"""
        settings_file = os.path.join(os.path.dirname(__file__), "settings.json")
        
        # Load existing settings or create new
        settings_data = {}
        if os.path.exists(settings_file):
            try:
                with open(settings_file, 'r') as f:
                    settings_data = json.load(f)
            except:
                pass
        
        # Get custom schemes only
        custom_schemes = []
        for name, colors in self.schemes.items():
            if name.startswith("Custom: "):
                custom_schemes.append({
                    "name": name.replace("Custom: ", ""),
                    "colors": colors
                })
        
        # Update settings
        settings_data["custom_color_schemes"] = custom_schemes
        settings_data["last_charge_scheme"] = self.current_scheme
        
        # Save to file
        try:
            with open(settings_file, 'w') as f:
                json.dump(settings_data, f, indent=2)
        except Exception as e:
            print(f"Error saving settings: {e}")
    
    def toggle_labels(self):
        """Toggle charge value labels in 3D view"""
        show = self.chk_show_labels.isChecked()
        
        # Remove old labels if they exist
        if hasattr(self, '_charge_labels'):
            for actor in self._charge_labels:
                try:
                    self.parent_dlg.mw.plotter.remove_actor(actor)
                except: pass
            self._charge_labels = []
        
        if not show:
            if hasattr(self.parent_dlg.mw, 'plotter'):
                self.parent_dlg.mw.plotter.render()
            return
        
        # Add labels
        data = self.all_charges.get(self.current_type, [])
        coords = self.parent_dlg.parser.data.get("coords", [])
        
        if not coords or len(coords) != len(data):
            QMessageBox.warning(self, "Error", "No coordinates available for labels.")
            self.chk_show_labels.setChecked(False)
            return
        
        self._charge_labels = []
        try:
            for item in data:
                idx = item['atom_idx']
                pos = coords[idx]
                label_pos = [pos[0], pos[1], pos[2] + 0.3]
                
                actor = self.parent_dlg.mw.plotter.add_point_labels(
                    [label_pos],
                    [f"{item['charge']:.2f}"],
                    font_size=10,
                    text_color='yellow',
                    point_size=0,
                    always_visible=True
                )
                self._charge_labels.append(actor)
            
            self.parent_dlg.mw.plotter.render()
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Could not add labels: {e}")
            self.chk_show_labels.setChecked(False)
    
    def reset_colors(self):
        """Reset atom colors to default CPK colors"""
        try:
            mw = self.parent_dlg.mw
            if hasattr(mw, 'main_window_view_3d'):
                data = self.all_charges.get(self.current_type, [])
                for item in data:
                    mw.main_window_view_3d.update_atom_color_override(item['atom_idx'], None)
            
            # Remove scalar bar if exists
            if hasattr(self, '_charge_scalar_bar'):
                try:
                    self.parent_dlg.mw.plotter.remove_actor(self._charge_scalar_bar)
                    delattr(self, '_charge_scalar_bar')
                except: pass
            
            # Remove labels if exist
            if hasattr(self, '_charge_labels'):
                for actor in self._charge_labels:
                    try:
                        self.parent_dlg.mw.plotter.remove_actor(actor)
                    except: pass
                self._charge_labels = []
                self.chk_show_labels.setChecked(False)
            
            # Redraw molecule with default CPK colors
            if hasattr(self.parent_dlg.mw, 'current_mol') and self.parent_dlg.mw.current_mol:
                if hasattr(self.parent_dlg.mw, 'draw_molecule_3d'):
                    self.parent_dlg.mw.draw_molecule_3d(self.parent_dlg.mw.current_mol)
            
            # Render
            if hasattr(self.parent_dlg.mw, 'plotter'):
                self.parent_dlg.mw.plotter.render()
            
            QMessageBox.information(self, "Done", "Colors reset to CPK default.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to reset colors:\n{e}")
        
    def update_table(self):
        data = self.all_charges.get(self.current_type, [])
        self.table.setRowCount(len(data))
        for r, item in enumerate(data):
            self.table.setItem(r, 0, QTableWidgetItem(str(item['atom_idx'])))
            self.table.setItem(r, 1, QTableWidgetItem(item['atom_sym']))
            self.table.setItem(r, 2, QTableWidgetItem(f"{item['charge']:.4f}"))
            
    def apply_colors(self):
        data = self.all_charges.get(self.current_type, [])
        if not data: return
        
        charges = [d['charge'] for d in data]
        if not charges: return
        
        max_c = max(abs(min(charges)), abs(max(charges)))
        if max_c == 0: max_c = 1.0
        
        # Get current color scheme
        colors = self.schemes.get(self.current_scheme, ["red", "white", "blue"])
        
        # Helper to interpolate color between multiple colors
        def get_color(q):
            # Normalize charge to 0.0 (negative) to 1.0 (positive)
            norm = (q + max_c) / (2 * max_c)
            norm = max(0.0, min(1.0, norm))
            
            n = len(colors)
            if n < 2: 
                return QColor(colors[0]).name()
            
            # Find which segment we're in
            seg_len = 1.0 / (n - 1)
            idx = int(norm / seg_len)
            if idx >= n - 1: 
                return QColor(colors[-1]).name()
            
            # Interpolate within segment
            local_t = (norm - (idx * seg_len)) / seg_len
            
            c1 = QColor(colors[idx])
            c2 = QColor(colors[idx + 1])
            
            r = int(c1.red() + (c2.red() - c1.red()) * local_t)
            g = int(c1.green() + (c2.green() - c1.green()) * local_t)
            b = int(c1.blue() + (c2.blue() - c1.blue()) * local_t)
            
            return f"#{r:02x}{g:02x}{b:02x}"

        try:
            mw = self.parent_dlg.mw
            if hasattr(mw, 'main_window_view_3d'):
                for item in data:
                    idx = item['atom_idx']
                    q = item['charge']
                    color = get_color(q)
                    mw.main_window_view_3d.update_atom_color_override(idx, color)
            
            # Redraw molecule once after all colors are set
            if hasattr(self.parent_dlg.mw, 'current_mol') and self.parent_dlg.mw.current_mol:
                if hasattr(self.parent_dlg.mw, 'draw_molecule_3d'):
                    self.parent_dlg.mw.draw_molecule_3d(self.parent_dlg.mw.current_mol)
                
            # Add charge labels if checkbox is enabled
            if self.chk_show_labels.isChecked():
                coords = self.parent_dlg.parser.data.get("coords", [])
                if coords and len(coords) == len(data):
                    # Remove old labels if exist
                    if hasattr(self, '_charge_labels'):
                        for actor in self._charge_labels:
                            try:
                                self.parent_dlg.mw.plotter.remove_actor(actor)
                            except: pass
                    
                    self._charge_labels = []
                    for i, item in enumerate(data):
                        pos = coords[item['atom_idx']]
                        label_pos = [pos[0], pos[1], pos[2] + 0.3]
                        actor = self.parent_dlg.mw.plotter.add_point_labels(
                            [label_pos],
                            [f"{item['charge']:.2f}"],
                            font_size=10,
                            text_color='yellow',
                            point_size=0,
                            always_visible=True
                        )
                        self._charge_labels.append(actor)
            
            # Add scalar bar legend showing charge scale with color gradient
            try:
                # Create a mesh with charge values to display scalar bar
                vmin = min(charges)
                vmax = max(charges)
                
                # Convert scheme colors to colormap
                cmap_colors = [mcolors.to_rgb(QColor(c).name()) for c in colors]
                cmap = LinearSegmentedColormap.from_list("charge_cmap", cmap_colors, N=256)
                
                # Remove old scalar bar if exists
                if hasattr(self, '_charge_scalar_bar'):
                    try:
                        self.parent_dlg.mw.plotter.remove_actor(self._charge_scalar_bar)
                    except: pass
                
                # Create dummy mesh for scalar bar
                dummy = pv.Box()
                dummy.point_data['charges'] = np.linspace(vmin, vmax, dummy.n_points)
                
                # Add mesh with scalar bar (invisible mesh, visible bar)
                self._charge_scalar_bar = self.parent_dlg.mw.plotter.add_mesh(
                    dummy,
                    scalars='charges',
                    cmap=cmap,
                    clim=[vmin, vmax],
                    opacity=0.0,  # Invisible mesh
                    show_scalar_bar=True,
                    scalar_bar_args={
                        'title': f'{self.current_type}',
                        'title_font_size': 14,
                        'label_font_size': 12,
                        'n_labels': 5,
                        'vertical': True,
                        'height': 0.3,
                        'width': 0.08,
                        'position_x': 0.88,
                        'position_y': 0.35,
                        'color': 'white'
                    }
                )
            except Exception as e:
                print(f"Error adding scalar bar: {e}")
            
            # Trigger update
            if hasattr(self.parent_dlg.mw, 'plotter'):
                self.parent_dlg.mw.plotter.render()
                
            QMessageBox.information(self, "Done", f"Applied '{self.current_scheme}' coloring to 3D view.")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to color atoms:\n{e}")
