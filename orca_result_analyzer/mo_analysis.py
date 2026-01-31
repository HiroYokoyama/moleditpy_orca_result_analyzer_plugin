
import os
import tempfile
import numpy as np
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, 
                             QTreeWidget, QTreeWidgetItem, QAbstractItemView, QMessageBox, 
                             QFileDialog, QProgressDialog, QTableWidget, QTableWidgetItem, 
                             QHeaderView, QGroupBox, QSpinBox, QDoubleSpinBox, QSplitter, QWidget,
                             QFormLayout, QTreeWidgetItemIterator, QApplication, QColorDialog, QInputDialog, QComboBox, QCheckBox)
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
            # Sanitize display ID (which might contain "1 (a)")
            safe_id = str(display_id).replace(" ", "_").replace("(", "").replace(")", "").replace(":", "")
            return os.path.join(out_dir, f"{filename_base}_MO_{safe_id}.cube")
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
        
        self.check_smooth = QCheckBox("Smooth Shading")
        self.check_smooth.setChecked(True)
        h_style.addWidget(self.check_smooth)
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
        self.check_smooth.toggled.connect(self.update_vis_only)
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
        self.tree.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
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
        
        btn_csv = QPushButton("Export CSV")
        btn_csv.clicked.connect(self.export_csv)
        btn_layout.addWidget(btn_csv)
                
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
        
        # Sort keys logic:
        # Keys are likely strings "0_alpha", "1_alpha" etc. or integers (if old/other format)
        # We want to sort by Spin then Index, or just by Energy if available?
        # Standard: Sort by Spin (Alpha first), then Index.
        
        keys = []
        if isinstance(self.mos, dict):
            raw_keys = list(self.mos.keys())
            # Check key type
            if raw_keys and isinstance(raw_keys[0], str) and "_" in raw_keys[0]:
                # Sort by spin, then index
                def sort_key(k):
                    if "_" in str(k):
                        idx_s, spin_s = str(k).split("_")
                        # alpha < beta < restricted ???
                        # "alpha" < "beta" (alphabetical ok).
                        return (spin_s, int(idx_s))
                    return (str(k), 0)
                keys = sorted(raw_keys, key=sort_key)
            else:
                # Fallback numeric sort
                keys = sorted(raw_keys, key=lambda x: int(x) if str(x).isdigit() else str(x))
                
            for k in keys:
                val = self.mos[k]
                if 'id' not in val: val['id'] = k
                # Store the Access Key used in self.mos dictionary
                val['_access_key'] = k 
                self.mo_list.append(val)
        elif isinstance(self.mos, list):
            self.mo_list = self.mos
            # Add index as access key
            for i, mo in enumerate(self.mo_list):
                 mo['_access_key'] = i
            
        # Separate HOMO/LUMO logic for each spin
        spin_mos = {'restricted': [], 'alpha': [], 'beta': []}
        for mo in self.mo_list:
            s = mo.get('spin', 'restricted')
            if s not in spin_mos: s = 'restricted' # fallback using restricted logic for unknown strings if any
            spin_mos[s].append(mo)
            
        spin_homo_idx = {}
        for s in spin_mos:
             n_occ = 0
             for mo in spin_mos[s]:
                 occ = mo.get('occ', mo.get('occupation', 0.0))
                 if occ > 0.1: n_occ += 1
             spin_homo_idx[s] = n_occ - 1

        # Display order: Reversed (High Energy Top)
        # We need to sort global list carefully if mixing spins?
        # Or just show them in the order we sorted above (Alpha then Beta)
        # Our self.mo_list is already sorted by Spin, Index.
        # But we iterate backwards! so Beta (high idx -> low idx) -> Alpha (high idx -> low idx)
        
        count = len(self.mo_list)
        for i in range(count - 1, -1, -1):
            mo = self.mo_list[i]
            
            # Internal unique key
            access_key = mo.get('_access_key')
            
            # Display ID
            mo_idx_val = mo.get('id', mo.get('index', i))
            spin = mo.get('spin', 'restricted')
            
            # Determine if HOMO/LUMO for this specific spin channel
            # We need the index OF THIS MO within its spin group
            # Since self.mo_list is sorted, we can find it?
            # Or just use the ID? Usually ID corresponds to index.
            # Let's use the ID if available, assuming 0-based index.
            
            local_idx = -1
            try:
                local_idx = int(mo_idx_val)
            except: pass
            
            is_homo = False
            is_lumo = False
            if spin in spin_homo_idx:
                h = spin_homo_idx[spin]
                if local_idx == h: is_homo = True
                elif local_idx == h + 1: is_lumo = True
            
            display_id = str(mo_idx_val)
            try:
                 display_id = str(int(mo_idx_val) + 1)
            except: pass
            
            label_id = display_id
            if spin == 'alpha':
                label_id += " (a)"
            elif spin == 'beta':
                label_id += " (b)"
                
            if is_homo: label_id += " (HOMO)"
            elif is_lumo: label_id += " (LUMO)"
            
            occ = mo.get('occ', mo.get('occupation', 0.0))
            e_eh = mo.get('energy_eh', mo.get('energy'))
            e_ev = mo.get('energy_ev')
            
            if e_eh is None and e_ev is not None: e_eh = e_ev / 27.2114
            elif e_ev is None and e_eh is not None: e_ev = e_eh * 27.2114
            if e_eh is None: e_eh = 0.0
            if e_ev is None: e_ev = 0.0
            
            item = QTreeWidgetItem([label_id, f"{occ:.2f}", f"{e_ev:.3f}", f"{e_eh:.5f}"])
            
            # Store the unique lookup key
            item.setData(0, Qt.ItemDataRole.UserRole, access_key)
            
            # Text Color based on Occ
            if occ > 1.9:
                item.setForeground(0, QColor("blue")) 
            elif occ > 0.1:
                item.setForeground(0, QColor("green"))
            else:
                item.setForeground(0, QColor("gray"))
            
            # Background Color if Generated
            # Need to match path logic with key? 
            # get_cube_path usually uses display_id (int).
            # If keys are strings, cache might break if we just use INT.
            # We should probably incorporate spin into cache filename.
            
            # Quick fix: Pass full label_id hash or just access_key to path
            # But get_cube_path expects display_id.
            # Let's keep existing check but maybe use access_key for path generation in future.
            
            
            self.tree.addTopLevelItem(item)
            
        # Scroll to HOMO of first available spin branch
        # Since we just fill them, scrolling to center might be ambiguous.
        # Let's scroll to the first 'Occupied' item found from top of tree (High energy)
        # That would be the HOMO of the last spin block added (e.g. Beta HOMO if present, or Alpha HOMO).
        # Actually tree has High Energy at top. So First Occupied Item is HOMO.
        
        iterator = QTreeWidgetItemIterator(self.tree)
        while iterator.value():
            item = iterator.value()
            if "HOMO" in item.text(0):
                 self.tree.scrollToItem(item, QAbstractItemView.ScrollHint.PositionAtCenter)
                 self.tree.setCurrentItem(item)
                 break
            iterator += 1

    def on_double_click(self, item, col):
        self.visualize_current_mo()

    def on_selection_changed(self):
        items = self.tree.selectedItems()
        has_coeffs = False
        if items:
            try:
                # Use UserRole data for robust lookup
                key = items[0].data(0, Qt.ItemDataRole.UserRole)
                
                # Check coeffs in parser data
                if hasattr(self.parent_dlg, 'parser') and self.parent_dlg.parser:
                    if hasattr(self.parent_dlg.parser, 'data') and self.parent_dlg.parser.data:
                        if "mo_coeffs" in self.parent_dlg.parser.data:
                            # Use key directly
                            if key in self.parent_dlg.parser.data["mo_coeffs"]:
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
        if self.mw and hasattr(self.mw, 'statusBar'):
            self.mw.statusBar().showMessage("ORCA Input block copied to clipboard.", 5000)
        else:
            print("ORCA Input block copied to clipboard.")

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
        try:
            items = self.tree.selectedItems()
            if not items: return
            
            key = items[0].data(0, Qt.ItemDataRole.UserRole)
            display_id = items[0].text(0) # Just for caching/display
            
        except: return
        
        # Check coefficients with KEY
        # Safety checks
        if not hasattr(self.parent_dlg, 'parser') or not self.parent_dlg.parser:
            QMessageBox.warning(self, "Error", "Parser not available")
            return
        if not hasattr(self.parent_dlg.parser, 'data') or not self.parent_dlg.parser.data:
            QMessageBox.warning(self, "Error", "No parser data available")
            return
            
        coeffs_map = self.parent_dlg.parser.data.get("mo_coeffs", {})
        mo_data = coeffs_map.get(key)
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
            engine, key, self.spin_pts.value(), self.spin_margin.value(),
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
            "iso": 0.02, "opacity": 0.5, "style": "Surface", "color_p": "#ff0000", "color_n": "#0000ff", "smooth_shading": True
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
            "last_preset": self.combo_presets.currentText(),
            "smooth_shading": self.check_smooth.isChecked()
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
            "color_n": self.get_color_hex('n'),
            "smooth_shading": self.check_smooth.isChecked()
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
        
        self.check_smooth.setChecked(data.get("smooth_shading", True))
        
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
                         color_p=cp, color_n=cn, style=style, smooth_shading=self.check_smooth.isChecked())
            mw.plotter.render()

    def closeEvent(self, event):
        """Clean up 3D actors when closing"""
        if hasattr(self.parent_dlg, 'mw'):
             plotter = self.parent_dlg.mw.plotter
             plotter.remove_actor("mo_iso_p")
             plotter.remove_actor("mo_iso_n")
             plotter.render()
        elif hasattr(self.parent_dlg, 'context'):
             plotter = self.parent_dlg.context.get_main_window().plotter
             plotter.remove_actor("mo_iso_p")
             plotter.remove_actor("mo_iso_n")
             plotter.render()
             
        super().closeEvent(event)

    def export_csv(self):
        filename, _ = QFileDialog.getSaveFileName(self, "Export MO Data", "", "CSV Files (*.csv)")
        if not filename: return
        
        try:
            import csv
            with open(filename, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(["ID", "Occupation", "Energy (eV)", "Energy (Eh)"])
                
                # Iterate tree items to respect current sort/display?
                # Using the data list is safer and complete.
                # However, the tree has formatted text. Let's use the list self.mo_list
                # But self.mo_list is raw. Let's replicate the display logic or use tree items.
                
                # Using tree items ensures we export what is seen (including calculated Energy if missing)
                it = QTreeWidgetItemIterator(self.tree)
                while it.value():
                    item = it.value()
                    # ID, Occ, eV, Eh
                    row = [item.text(0), item.text(1), item.text(2), item.text(3)]
                    writer.writerow(row)
                    it += 1
            # print(f"Data exported to {filename}")
            # QMessageBox.information(self, "Success", f"Data exported to {filename}")
            if hasattr(self.parent_dlg, 'mw') and self.parent_dlg.mw:
                self.parent_dlg.mw.statusBar().showMessage(f"Data exported to {filename}", 5000)
            elif hasattr(self.parent_dlg, 'context'):
                 self.parent_dlg.context.get_main_window().statusBar().showMessage(f"Data exported to {filename}", 5000)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to export CSV: {e}")
