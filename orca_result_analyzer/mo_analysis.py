
import os
import tempfile
import numpy as np
import webbrowser
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, 
                             QTreeWidget, QTreeWidgetItem, QAbstractItemView, QMessageBox, 
                             QFileDialog, QProgressDialog, QTableWidget, QTableWidgetItem, 
                             QHeaderView, QGroupBox, QSpinBox, QDoubleSpinBox, QSplitter, QWidget,
                             QFormLayout, QTreeWidgetItemIterator)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor, QBrush

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
        parser = self.parent_dlg.parser
        if hasattr(parser, 'filename') and parser.filename:
            fpath = parser.filename
            base_dir = os.path.dirname(fpath)
            filename_base = os.path.splitext(os.path.basename(fpath))[0]
            out_dir = os.path.join(base_dir, f"{filename_base}_cubes")
            return os.path.join(out_dir, f"{filename_base}_MO_{display_id}.cube")
        return None

    def setup_ui(self):
        # Use simpler Vertical Layout to fill the window
        layout = QVBoxLayout(self)
        
        # 1. Visualization Settings (Grouping "Grid" vs "Calc")
        vis_grp = QGroupBox("Visualization Controls")
        vis_layout = QVBoxLayout(vis_grp)
        
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
                if self.parent_dlg.parser and "mo_coeffs" in self.parent_dlg.parser.data:
                    if internal_id in self.parent_dlg.parser.data["mo_coeffs"]:
                        has_coeffs = True
            except: pass
            
        self.btn_vis.setEnabled(has_coeffs)
        if items and not has_coeffs:
            self.btn_vis.setToolTip(f"No coefficients available for MO {items[0].text(0)}")
        else:
            self.btn_vis.setToolTip("")

    def get_engine(self):
        if not BasisSetEngine:
            QMessageBox.critical(self, "Error", "BasisSetEngine not available")
            return None
        try:
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

    def show_cube(self, path):
        if not CubeVisualizer: return
        
        # Access Main Window
        mw = None
        if hasattr(self.parent_dlg, 'mw'): mw = self.parent_dlg.mw
        elif hasattr(self.parent_dlg, 'context'): mw = self.parent_dlg.context.get_main_window()
        
        if not mw: return
        
        vis = CubeVisualizer(mw)
        if vis.load_file(path):
            vis.show_iso(self.spin_iso.value(), opacity=self.spin_opacity.value())
            mw.plotter.render()
            # self.lbl_info.setText(f"Visualizing MO {path} \nIso={self.spin_iso.value()}")
            # Bring main window to focus? No, dialog is modal.
            # Ideally dialog should be modeless.
        
    def update_vis_only(self):
        if self.last_cube_path and os.path.exists(self.last_cube_path):
            self.show_cube(self.last_cube_path)

    def generate_cubes_batch(self):
        # (Reuse previous batch logic but simplified or just call it)
        # For brevity, implementing a simpler version or reusing logic if methods existed.
        QMessageBox.information(self, "TODO", "Batch export not fully refactored yet. Use individual visualize for now.")

    def view_coeffs(self):
        # Disabled
        pass
        
    def closeEvent(self, event):
        self.cleanup_vis()
        super().closeEvent(event)
        
    def cleanup_vis(self):
        if not CubeVisualizer: return
        
        # Access Main Window
        mw = None
        if hasattr(self.parent_dlg, 'mw'): mw = self.parent_dlg.mw
        elif hasattr(self.parent_dlg, 'context'): mw = self.parent_dlg.context.get_main_window()
        
        if mw:
            vis = CubeVisualizer(mw)
            if hasattr(vis, 'clear'):
                vis.clear()
            else:
                # Fallback if logic mismatch
                if hasattr(mw.plotter, 'remove_actor'):
                    mw.plotter.remove_actor("mo_iso_p")
                    mw.plotter.remove_actor("mo_iso_n")
                    mw.plotter.render()
