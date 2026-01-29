
import os
import tempfile
import numpy as np
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, 
                             QTreeWidget, QTreeWidgetItem, QAbstractItemView, QMessageBox, 
                             QFileDialog, QProgressDialog, QTableWidget, QTableWidgetItem, 
                             QHeaderView, QGroupBox, QSpinBox, QDoubleSpinBox, QSplitter, QWidget,
                             QFormLayout, QTreeWidgetItemIterator, QApplication, QColorDialog, QInputDialog, QComboBox)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor, QBrush
import json

try:
    from .mo_engine import BasisSetEngine, CalcWorker
    from .vis import CubeVisualizer
except ImportError:
    try:
        from mo_engine import BasisSetEngine, CalcWorker
        from vis import CubeVisualizer
    except:
        BasisSetEngine = None
        CalcWorker = None
        CubeVisualizer = None

class MODialog(QDialog):
    def __init__(self, parent, mos):
        super().__init__(parent)
        self.setWindowTitle("MO Analysis & Visualization")
        self.resize(550, 750) 
        self.mos = mos
        self.parent_dlg = parent
        self.last_cube_path = None
        self.setup_ui()

    def get_cube_path(self, display_id):
        if not hasattr(self.parent_dlg, 'parser') or not self.parent_dlg.parser:
            return None
            
        parser = self.parent_dlg.parser
        if hasattr(parser, 'filename') and parser.filename:
            fpath = parser.filename
            base_dir = os.path.dirname(fpath)
            filename_base = os.path.splitext(os.path.basename(fpath))[0]
            out_dir = os.path.join(base_dir, f"{filename_base}_cubes")
            # Ensure dir exists? Logic doesn't seem to create it here, but usually caller handles or it's fine.
            # Standardizing path return.
            return os.path.join(out_dir, f"{filename_base}_MO_{display_id}.cube")
        return None

    def setup_ui(self):
        # Use simpler Vertical Layout to fill the window
        layout = QVBoxLayout(self)
        
        # 1. Visualization Settings (Grouping "Grid" vs "Calc")
        vis_grp = QGroupBox("Visualization Controls")
        vis_layout = QVBoxLayout(vis_grp)
        
        # --- Presets ---
        h_preset = QHBoxLayout()
        h_preset.addWidget(QLabel("Preset:"))
        self.combo_presets = QComboBox()
        self.combo_presets.addItem("Default")
        h_preset.addWidget(self.combo_presets)
        
        btn_save_preset = QPushButton("Save")
        btn_save_preset.setFixedWidth(50)
        btn_save_preset.clicked.connect(self.save_preset)
        h_preset.addWidget(btn_save_preset)
        
        btn_del_preset = QPushButton("Del")
        btn_del_preset.setFixedWidth(40)
        btn_del_preset.clicked.connect(self.delete_preset)
        h_preset.addWidget(btn_del_preset)
        vis_layout.addLayout(h_preset)
        
        # --- Style & Colors ---
        h_style = QHBoxLayout()
        h_style.addWidget(QLabel("Style:"))
        self.combo_style = QComboBox()
        self.combo_style.addItems(["Surface", "Wireframe", "Points"])
        h_style.addWidget(self.combo_style)
        vis_layout.addLayout(h_style)
        
        h_colors = QHBoxLayout()
        self.btn_color_p = QPushButton("Pos (+)")
        self.btn_color_p.setStyleSheet("background-color: red; color: white; font-weight: bold;")
        self.btn_color_p.clicked.connect(lambda: self.pick_color('p'))
        h_colors.addWidget(self.btn_color_p)
        
        self.btn_color_n = QPushButton("Neg (-)")
        self.btn_color_n.setStyleSheet("background-color: blue; color: white; font-weight: bold;")
        self.btn_color_n.clicked.connect(lambda: self.pick_color('n'))
        h_colors.addWidget(self.btn_color_n)
        vis_layout.addLayout(h_colors)
        
        # Isovalue / Opacity
        h_iso = QHBoxLayout()
        h_iso.addWidget(QLabel("Isovalue:"))
        self.spin_iso = QDoubleSpinBox()
        self.spin_iso.setRange(0.001, 1.0)
        self.spin_iso.setSingleStep(0.005)
        self.spin_iso.setDecimals(3)
        self.spin_iso.setValue(0.02)
        h_iso.addWidget(self.spin_iso)
        
        h_iso.addWidget(QLabel("Opacity:"))
        self.spin_opacity = QDoubleSpinBox()
        self.spin_opacity.setRange(0.0, 1.0)
        self.spin_opacity.setSingleStep(0.1)
        self.spin_opacity.setValue(0.5)
        h_iso.addWidget(self.spin_opacity)
        vis_layout.addLayout(h_iso)
        
        # Connect changes to update view immediately if exists
        self.combo_style.currentTextChanged.connect(self.update_vis_only)
        self.spin_iso.valueChanged.connect(self.update_vis_only)
        self.spin_opacity.valueChanged.connect(self.update_vis_only)
        
        # Load Settings must be called after UI init
        self.load_settings()
        self.combo_presets.currentTextChanged.connect(self.apply_preset)
        
        # Calculation Settings
        calc_grp = QGroupBox("Calculation Parameters")
        calc_layout = QFormLayout(calc_grp)
        
        self.spin_pts = QSpinBox()
        self.spin_pts.setRange(10, 200)
        self.spin_pts.setValue(40)
        self.spin_pts.setSuffix(" pts")
        calc_layout.addRow("Grid Resolution (x,y,z):", self.spin_pts)
        
        self.spin_margin = QDoubleSpinBox()
        self.spin_margin.setRange(1.0, 15.0)
        self.spin_margin.setValue(4.0)
        self.spin_margin.setSuffix(" Bohr") # Correct unit
        calc_layout.addRow("Calc Boundary Margin:", self.spin_margin)
        
        vis_layout.addWidget(calc_grp)
        
        
        # Warning Layout (Label + Copy Button)
        warn_layout = QHBoxLayout()
        
        self.lbl_warning = QLabel("")
        self.lbl_warning.setStyleSheet("color: #e65100; font-weight: bold; font-size: 9pt;")
        self.lbl_warning.setWordWrap(True)
        self.lbl_warning.setVisible(False)
        warn_layout.addWidget(self.lbl_warning)
        
        self.btn_copy_input = QPushButton("Copy Input")
        self.btn_copy_input.setFixedWidth(80)
        self.btn_copy_input.setVisible(False)
        self.btn_copy_input.clicked.connect(self.copy_orca_input)
        warn_layout.addWidget(self.btn_copy_input)
        
        vis_layout.addLayout(warn_layout)
        
        layout.addWidget(vis_grp)
        
        # 2. MO List
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["MO", "Occ", "E (eV)", "E (Eh)"])
        # Columns: Let them fill
        self.tree.setColumnWidth(0, 150)
        self.tree.setColumnWidth(1, 80)
        self.tree.setColumnWidth(2, 80)
        # Last column stretches
        header = self.tree.header()
        if header:
            header.setSectionResizeMode(3, QHeaderView.ResizeMode.Stretch)
            
        self.tree.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.tree.itemDoubleClicked.connect(self.on_double_click)
        layout.addWidget(self.tree)
        
        # 3. Action Buttons
        btn_layout = QHBoxLayout()
        self.btn_vis = QPushButton("Visualize Selected")
        self.btn_vis.setStyleSheet("font-weight: bold; background-color: #d0f0c0;")
        self.btn_vis.clicked.connect(self.visualize_current_mo)
        self.btn_vis.setEnabled(False) # Default disabled until selection
        btn_layout.addWidget(self.btn_vis)
        
        # Connect selection
        self.tree.itemSelectionChanged.connect(self.on_selection_changed)
        
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        btn_layout.addWidget(btn_close)
        
        # Add Stretch to push buttons to right or keep centered? 
        # Standard dialogs often have buttons on right.
        # But centering is fine.
        
        # IMPORTANT: Add the button layout to the main layout!
        layout.addLayout(btn_layout)
        # If we want interaction with 3D view, we need Non-Modal (show()).
        
        # For now, we keep it Modal, but 'Visualize' updates the background window?
        # If modal, user can't rotate 3D view easily without closing dialog.
        # Let's add a "Apply/Update" button or just let it update.
        
        # To make it better:
        # We'll just show info on right side.
        # info_panel = QWidget()
        # info_layout = QVBoxLayout(info_panel)
        # info_layout.addWidget(QLabel("<b>Selected MO Info</b>"))
        # self.lbl_info = QLabel("Select an MO to view details.")
        # self.lbl_info.setWordWrap(True)
        # info_layout.addWidget(self.lbl_info)
        # info_layout.addStretch()
        # main_layout.addWidget(info_panel)
        
        # Populate
        self.normalize_and_populate()

    def normalize_and_populate(self):
        self.tree.clear() # Verify clear first
        self.mo_list = []
        if isinstance(self.mos, dict):
            for k in sorted(self.mos.keys()):
                val = self.mos[k]
                if 'id' not in val: val['id'] = k
                self.mo_list.append(val)
        elif isinstance(self.mos, list):
            self.mo_list = self.mos
            
        # Estimate HOMO/LUMO
        # Count electrons? Or check occ drop
        n_occ = 0
        for mo in self.mo_list:
            occ = mo.get('occ', mo.get('occupation', 0.0))
            if occ > 0.1: n_occ += 1
            
        homo_idx = n_occ - 1
        
        # Display order: Reversed (High Energy Top, Low Energy Bottom) - "Orbital Diagram" style
        # Iterate backwards
        count = len(self.mo_list)
        for i in range(count - 1, -1, -1):
            mo = self.mo_list[i]
            
            # Internal mo_id is 0-based
            mo_id = mo.get('id', mo.get('index', i))
            occ = mo.get('occ', mo.get('occupation', 0.0))
            e_eh = mo.get('energy_eh', mo.get('energy'))
            e_ev = mo.get('energy_ev')
            
            if e_eh is None and e_ev is not None: e_eh = e_ev / 27.2114
            elif e_ev is None and e_eh is not None: e_ev = e_eh * 27.2114
            if e_eh is None: e_eh = 0.0
            if e_ev is None: e_ev = 0.0
            
            # Display MO ID as 1-based (User Request)
            display_id = mo_id + 1
            label_id = str(display_id)
            if i == homo_idx: label_id += " (HOMO)"
            elif i == homo_idx + 1: label_id += " (LUMO)"
            
            item = QTreeWidgetItem([label_id, f"{occ:.2f}", f"{e_ev:.3f}", f"{e_eh:.5f}"])
            
            # Text Color based on Occ
            if occ > 1.9:
                item.setForeground(0, QColor("blue")) # HOMO/Occupied
            elif occ > 0.1:
                item.setForeground(0, QColor("green"))
            else:
                item.setForeground(0, QColor("gray")) # Virtual
            
            # Background Color if Generated (User Request)
            path = self.get_cube_path(display_id)
            if path and os.path.exists(path):
                # Light Green Background (Lighter shade requested)
                bg = QBrush(QColor(240, 255, 240))
                for c in range(4):
                    item.setBackground(c, bg)

            self.tree.addTopLevelItem(item)
            
        # Scroll to HOMO (Default View)
        # Since list is reversed, index of HOMO in tree is: (N - 1) - homo_idx
        # homo_idx is 0-based index in original sorted list.
        # Tree index: 0 is highest energy (last item in original list).
        if homo_idx >= 0:
            tree_homo_idx = (count - 1) - homo_idx
            if tree_homo_idx >= 0 and tree_homo_idx < self.tree.topLevelItemCount():
                item = self.tree.topLevelItem(tree_homo_idx)
                self.tree.scrollToItem(item, QAbstractItemView.ScrollHint.PositionAtCenter)
                self.tree.setCurrentItem(item)

    def on_double_click(self, item, col):
        self.visualize_current_mo()

    def on_selection_changed(self):
        items = self.tree.selectedItems()
        has_coeffs = False
        if items:
            try:
                mo_id_txt = items[0].text(0).split()[0]
                display_id = int(mo_id_txt)
                internal_id = display_id - 1
                
                # Check coeffs in parser data
                # Safety check for parser and data
                if hasattr(self.parent_dlg, 'parser') and self.parent_dlg.parser:
                    if hasattr(self.parent_dlg.parser, 'data') and self.parent_dlg.parser.data:
                        if "mo_coeffs" in self.parent_dlg.parser.data:
                            if internal_id in self.parent_dlg.parser.data["mo_coeffs"]:
                                has_coeffs = True
            except: pass
            
        self.btn_vis.setEnabled(has_coeffs)
        if items and not has_coeffs:
            self.btn_vis.setToolTip(f"No coefficients available for MO {items[0].text(0)}")
            self.lbl_warning.setText(
                "<b>Coefficients Missing.</b> Required Input:"
                "<pre style='margin-top:0; margin-bottom:0;'>%output<br>  Print[P_Basis] 2<br>  Print[P_Mos] 1<br>end</pre>"
            )
            self.lbl_warning.setVisible(True)
            self.btn_copy_input.setVisible(True)
        else:
            self.btn_vis.setToolTip("")
            self.lbl_warning.setVisible(False)
            self.btn_copy_input.setVisible(False)

    def copy_orca_input(self):
        text = "%output\n  Print[P_Basis] 2\n  Print[P_Mos] 1\nend"
        QApplication.clipboard().setText(text)
        QMessageBox.information(self, "Copied", "ORCA Input block copied to clipboard.")

    def get_engine(self):
        if not BasisSetEngine:
            QMessageBox.critical(self, "Error", "BasisSetEngine not available")
            return None
        try:
            # Safety checks
            if not hasattr(self.parent_dlg, 'parser') or not self.parent_dlg.parser:
                return None
            if not hasattr(self.parent_dlg.parser, 'data') or not self.parent_dlg.parser.data:
                return None
                
            shells = self.parent_dlg.parser.data.get("basis_set_shells", [])
            
            # Convert lists to numpy arrays for the engine AND map keys
            # IMPORTANT: Engine expects Bohr for centers (matching internal grid)
            # Parser provides Angstroms.
            BOHR_TO_ANG = 0.529177249
            
            clean_shells = []
            for s in shells:
                center_ang = np.array(s.get('origin', s.get('center', [0,0,0])))
                # Convert Angstrom -> Bohr
                center_bohr = center_ang / BOHR_TO_ANG
                
                d = {
                    'type': s.get('l', s.get('type', 0)),
                    'center': center_bohr, 
                    'exps': np.array(s['exps']),
                    'coeffs': np.array(s['coeffs'])
                }
                clean_shells.append(d)
                
            return BasisSetEngine(clean_shells)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Engine Init Failed: {e}")
            return None

    def visualize_current_mo(self):
        items = self.tree.selectedItems()
        if not items: return
        try:
            mo_id_txt = items[0].text(0).split()[0]
            display_id = int(mo_id_txt)
            internal_mo_id = display_id - 1
        except: return
        
        # Check coefficients with INTERNAL ID
        # Safety checks
        if not hasattr(self.parent_dlg, 'parser') or not self.parent_dlg.parser:
            QMessageBox.warning(self, "Error", "Parser not available")
            return
        if not hasattr(self.parent_dlg.parser, 'data') or not self.parent_dlg.parser.data:
            QMessageBox.warning(self, "Error", "No parser data available")
            return

        coeffs_map = self.parent_dlg.parser.data.get("mo_coeffs", {})
        mo_data = coeffs_map.get(internal_mo_id)
        if not mo_data:
            QMessageBox.warning(self, "Error", f"No coefficients for MO {display_id}")
            return
            
        engine = self.get_engine()
        if not engine: return
        
        # Prepare Data
        # Flatten coeffs
        raw_coeffs = [x['coeff'] for x in mo_data['coeffs']]
        dense_vec = np.zeros(engine.n_basis)
        n = min(len(raw_coeffs), engine.n_basis)
        dense_vec[:n] = raw_coeffs[:n]
        
        # Path Logic via Helper
        out_path = self.get_cube_path(display_id)
        if not out_path:
             # Fallback if no filename
             out_path = os.path.join(tempfile.gettempdir(), f"orca_mo_{display_id}.cube")
        
        self.last_cube_path = out_path
        
        # Check if exists
        if os.path.exists(out_path):
             # Reuse!
             #print(f"Loading cached cube: {out_path}")
             self.show_cube(out_path)
             return

        # Worker
        dialog = QProgressDialog("Generating Cube...", "Cancel", 0, 100, self)
        dialog.setWindowModality(Qt.WindowModality.WindowModal)
        
        self.worker = CalcWorker(
            engine, internal_mo_id, self.spin_pts.value(), self.spin_margin.value(),
            self.parent_dlg.parser.data["atoms"], self.parent_dlg.parser.data["coords"], dense_vec, out_path
        )
        
        self.worker.progress_sig.connect(dialog.setValue)
        
        def on_finished(success, res):
            dialog.close()
            if success:
                self.show_cube(res)
                # Mark as Generated (Green Background)
                # We need to find the item in tree. 
                # Since tree is reversed, iterating to find it is easiest or calculating index.
                # display_id is unique.
                # Just iterate all items?
                # Just iterate all items?
                it = QTreeWidgetItemIterator(self.tree)
                bg = QBrush(QColor(240, 255, 240)) # Light Green (Lighter)
                while it.value():
                    item = it.value()
                    txt = item.text(0).split()[0]
                    if txt == str(display_id):
                        for c in range(4):
                            item.setBackground(c, bg)
                        break
                    it += 1
            else:
                QMessageBox.critical(self, "Error", f"Calculation failed: {res}")
                
        self.worker.finished_sig.connect(on_finished)
        self.worker.start()
        dialog.exec()

    def pick_color(self, which):
        current_col = QColor("red") if which == 'p' else QColor("blue")
        # Try to parse from button style
        try:
            style = self.btn_color_p.styleSheet() if which == 'p' else self.btn_color_n.styleSheet()
            if "background-color:" in style:
                c_str = style.split("background-color:")[1].split(";")[0].strip()
                current_col = QColor(c_str)
        except: pass
        
        col = QColorDialog.getColor(current_col, self, "Select Color")
        if col.isValid():
            hex_c = col.name()
            # Determine contrasting text color
            # Simple brightness check
            brightness = (col.red() * 299 + col.green() * 587 + col.blue() * 114) / 1000
            text_c = "black" if brightness > 128 else "white"
            
            style_sheet = f"background-color: {hex_c}; color: {text_c}; font-weight: bold;"
            if which == 'p':
                self.btn_color_p.setStyleSheet(style_sheet)
            else:
                self.btn_color_n.setStyleSheet(style_sheet)
            self.update_vis_only()

    def load_settings(self):
        self.settings_file = os.path.join(os.path.dirname(__file__), "settings.json")
        self.presets = {"Default": {
            "iso": 0.02, "opacity": 0.5, "style": "Surface", "color_p": "#ff0000", "color_n": "#0000ff"
        }}
        
        if os.path.exists(self.settings_file):
            try:
                with open(self.settings_file, 'r') as f:
                    all_settings = json.load(f)
                    mo_settings = all_settings.get("mo_settings", {})
                    
                    # Load presets
                    saved_presets = mo_settings.get("presets", {})
                    for name, data in saved_presets.items():
                        self.presets[name] = data
                    
                    # Last used
                    last_preset = mo_settings.get("last_preset", "Default")
                    
                    # Populate combo
                    self.combo_presets.blockSignals(True)
                    self.combo_presets.clear()
                    self.combo_presets.addItems(list(self.presets.keys()))
                    
                    if last_preset in self.presets:
                        self.combo_presets.setCurrentText(last_preset)
                        self.apply_preset(last_preset)
                    else:
                        self.combo_presets.setCurrentText("Default")
                        self.apply_preset("Default")
                        
                    self.combo_presets.blockSignals(False)
            except Exception as e:
                print(f"Error loading settings: {e}")

    def save_settings(self):
        # Save current presets and selection
        all_settings = {}
        if os.path.exists(self.settings_file):
            try:
                with open(self.settings_file, 'r') as f:
                    all_settings = json.load(f)
            except: pass
            
        mo_settings = {
            "presets": {k:v for k,v in self.presets.items() if k != "Default"},
            "last_preset": self.combo_presets.currentText()
        }
        all_settings["mo_settings"] = mo_settings
        
        try:
            with open(self.settings_file, 'w') as f:
                json.dump(all_settings, f, indent=2)
        except Exception as e:
            print(f"Error saving settings: {e}")

    def save_preset(self):
        name, ok = QInputDialog.getText(self, "Save Preset", "Preset Name:")
        if not ok or not name: return
        
        # Get current state
        data = {
            "iso": self.spin_iso.value(),
            "opacity": self.spin_opacity.value(),
            "style": self.combo_style.currentText(),
            "color_p": self.get_color_hex('p'),
            "color_n": self.get_color_hex('n')
        }
        
        self.presets[name] = data
        
        # Update combo
        curr = self.combo_presets.currentText()
        self.combo_presets.blockSignals(True)
        self.combo_presets.clear()
        self.combo_presets.addItems(list(self.presets.keys()))
        self.combo_presets.setCurrentText(name) # switch to new
        self.combo_presets.blockSignals(False)
        
        self.save_settings()

    def delete_preset(self):
        curr = self.combo_presets.currentText()
        if curr == "Default":
            QMessageBox.warning(self, "Error", "Cannot delete Default preset.")
            return
            
        if curr in self.presets:
            del self.presets[curr]
            
        self.combo_presets.blockSignals(True)
        self.combo_presets.clear()
        self.combo_presets.addItems(list(self.presets.keys()))
        self.combo_presets.setCurrentText("Default")
        self.apply_preset("Default")
        self.combo_presets.blockSignals(False)
        
        self.save_settings()

    def apply_preset(self, name):
        if name not in self.presets: return
        data = self.presets[name]
        
        # Block signals to avoid partial updates? No, updates are fine.
        self.spin_iso.setValue(data.get("iso", 0.02))
        self.spin_opacity.setValue(data.get("opacity", 0.5))
        
        style = data.get("style", "Surface")
        idx = self.combo_style.findText(style)
        if idx >= 0: self.combo_style.setCurrentIndex(idx)
        
        cp = data.get("color_p", "#ff0000")
        cn = data.get("color_n", "#0000ff")
        self.set_btn_color(self.btn_color_p, cp)
        self.set_btn_color(self.btn_color_n, cn)
        
        self.update_vis_only()
        self.save_settings() # Save last used

    def get_color_hex(self, which):
        # Extract from stylesheet
        btn = self.btn_color_p if which == 'p' else self.btn_color_n
        style = btn.styleSheet()
        if "background-color:" in style:
            return style.split("background-color:")[1].split(";")[0].strip()
        return "#ff0000" if which == 'p' else "#0000ff"

    def set_btn_color(self, btn, hex_c):
        col = QColor(hex_c)
        brightness = (col.red() * 299 + col.green() * 587 + col.blue() * 114) / 1000
        text_c = "black" if brightness > 128 else "white"
        btn.setStyleSheet(f"background-color: {hex_c}; color: {text_c}; font-weight: bold;")

    def update_vis_only(self):
        if self.last_cube_path and os.path.exists(self.last_cube_path):
            self.show_cube(self.last_cube_path)

    def show_cube(self, path):
        if not CubeVisualizer: return
        
        # Access Main Window
        mw = None
        if hasattr(self.parent_dlg, 'mw'): mw = self.parent_dlg.mw
        elif hasattr(self.parent_dlg, 'context'): mw = self.parent_dlg.context.get_main_window()
        
        if not mw: return
        
        cp = self.get_color_hex('p')
        cn = self.get_color_hex('n')
        style = self.combo_style.currentText().lower() # surface, wireframe, points
        
        vis = CubeVisualizer(mw)
        if vis.load_file(path):
            vis.show_iso(self.spin_iso.value(), opacity=self.spin_opacity.value(), 
                         color_p=cp, color_n=cn, style=style)
            mw.plotter.render()
