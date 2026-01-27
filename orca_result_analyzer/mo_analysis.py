
import os
import numpy as np
import webbrowser
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, 
                             QTreeWidget, QTreeWidgetItem, QAbstractItemView, QMessageBox, 
                             QFileDialog, QProgressDialog, QTableWidget, QTableWidgetItem, 
                             QHeaderView)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor

try:
    from .mo_engine import BasisSetEngine, CalcWorker
except ImportError:
    # If running as script or in flat struct
    try:
        from mo_engine import BasisSetEngine, CalcWorker
    except:
        BasisSetEngine = None
        CalcWorker = None

class MODialog(QDialog):
    def __init__(self, parent, mos):
        super().__init__(parent)
        self.setWindowTitle("MO Energy Levels")
        self.resize(600, 500)
        self.mos = mos
        self.parent_dlg = parent
        
        layout = QVBoxLayout(self)
        
        # Simple List for now, maybe diagram later
        
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["No.", "Occ", "Energy (eV)", "Energy (Eh)"])
        self.tree.setColumnWidth(0, 60)
        self.tree.setColumnWidth(1, 60)
        self.tree.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection) # Allow multiple
        
        for mo in self.mos:
            item = QTreeWidgetItem([str(mo['id']), f"{mo['occ']:.2f}", f"{mo['energy_ev']:.3f}", f"{mo['energy_eh']:.5f}"])
            
            # Color code by occupancy
            if mo['occ'] > 1.0:
                item.setForeground(0, QColor("blue")) # Occupied
                item.setForeground(1, QColor("blue"))
            elif mo['occ'] > 0.0:
                item.setForeground(0, QColor("green")) # SOMO
            else:
                item.setForeground(0, QColor("gray")) # Virtual
                
            self.tree.addTopLevelItem(item)
            
        layout.addWidget(self.tree)
        
        # Buttons
        btn_layout = QHBoxLayout()
        
        self.btn_coeffs = QPushButton("View Composition")
        self.btn_coeffs.clicked.connect(self.view_coeffs)
        btn_layout.addWidget(self.btn_coeffs)
        
        # Native 3D generator
        self.btn_gen_3d = QPushButton("Generate Cube(s)")
        self.btn_gen_3d.clicked.connect(self.generate_cubes)
        btn_layout.addWidget(self.btn_gen_3d)

        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        btn_layout.addWidget(btn_close)
        
        layout.addLayout(btn_layout)

    def generate_cubes(self):
        # 1. Get Selection
        items = self.tree.selectedItems()
        if not items:
            QMessageBox.warning(self, "Selection", "Please select at least one MO.")
            return
            
        # 2. Check Data Availability
        # parent is OrcaResultAnalyzerDialog
        parser = self.parent_dlg.parser
        shells = parser.data.get("basis_set_shells", [])
        if not shells:
            QMessageBox.warning(self, "Error", "No Basis Set definitions found.\nMake sure output contains 'BASIS SET IN INPUT FORMAT'.\n(Add ! PrintBasis to input)")
            return
            
        coeffs_map = parser.data.get("mo_coeffs", {})
        if not coeffs_map:
             QMessageBox.warning(self, "Error", "No MO coefficients parsed.")
             return
             
        # 3. Directory Selection
        # If multiple, ask for directory. If single, ask for file?
        # Let's standardize: Ask for directory to save MO_X.cube, MO_Y.cube etc.
        out_dir = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if not out_dir: return
        
        # 4. Prepare Batch
        if BasisSetEngine is None:
             QMessageBox.critical(self, "Error", "BasisSetEngine not available. Check mo_engine.py.")
             return
        
        # Init Engine once
        try:
            engine = BasisSetEngine(shells)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to initialize Basis Engine: {e}")
            return
            
        tasks = []
        for item in items:
            try:
                mo_id = int(item.text(0))
            except: continue
            
            # Get coefficients from parser
            mo_data = coeffs_map.get(mo_id)
            if not mo_data: continue
            
            # Extract values list (assuming parser kept order)
            raw_coeffs = [x['coeff'] for x in mo_data['coeffs']]
            
            # Build dense vector (pad/truncate if needed)
            dense_vec = np.zeros(engine.n_basis)
            n = min(len(raw_coeffs), engine.n_basis)
            dense_vec[:n] = raw_coeffs[:n]
            
            fname = f"orbitals_mo_{mo_id}.cube"
            fpath = os.path.join(out_dir, fname)
            
            tasks.append({
                'mo_id': mo_id,
                'coeffs': dense_vec,
                'path': fpath
            })
            
        if not tasks: return

        # 5. Run Batch
        self.run_cube_batch(engine, tasks, out_dir)

    def run_cube_batch(self, engine, tasks, out_dir):
        # Interactive progress dialog
        
        self.cube_pd = QProgressDialog("Generating Cubes...", "Cancel", 0, len(tasks)*100, self)
        self.cube_pd.setWindowModality(Qt.WindowModality.WindowModal)
        self.cube_pd.setMinimumDuration(0)
        
        self.cube_tasks = tasks
        self.cube_engine = engine
        self.cube_current_idx = 0
        self.cube_out_dir = out_dir
        
        self.start_next_cube_task()
        
    def start_next_cube_task(self):
        if self.cube_current_idx >= len(self.cube_tasks):
            self.cube_pd.setValue(len(self.cube_tasks)*100)
            self.cube_pd.close()
            
            # Ask to open Folder or load one?
            reply = QMessageBox.question(self, "Done", 
                f"Generated {len(self.cube_tasks)} cube files in:\n{self.cube_out_dir}\n\nOpen output directory?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
                
            if reply == QMessageBox.StandardButton.Yes and os.path.exists(self.cube_out_dir):
                webbrowser.open(self.cube_out_dir)
            return
            
        if self.cube_pd.wasCanceled():
            return
            
        task = self.cube_tasks[self.cube_current_idx]
        self.cube_pd.setLabelText(f"Generating MO {task['mo_id']}...")
        
        # Params
        n_points = 40 # Configurable?
        margin = 3.0
        
        parser = self.parent_dlg.parser
        
        worker = CalcWorker(
            self.cube_engine, 
            task['mo_id'], 
            n_points, 
            margin, 
            parser.data["atoms"], 
            parser.data["coords"], 
            task['coeffs'], 
            task['path']
        )
        
        worker.progress_sig.connect(self.on_cube_progress)
        worker.finished_sig.connect(self.on_cube_finished)
        
        self.current_cube_worker = worker
        worker.start()
        
    def on_cube_progress(self, pct):
        if self.cube_pd.wasCanceled():
            if self.current_cube_worker:
                self.current_cube_worker._is_cancelled = True
            return
            
        base = self.cube_current_idx * 100
        self.cube_pd.setValue(base + pct)
        
    def on_cube_finished(self, success, msg):
        if not success:
            QMessageBox.warning(self, "Error", f"Failed to generate cube: {msg}")
            
        self.cube_current_idx += 1
        self.start_next_cube_task()


    def view_coeffs(self):
        item = self.tree.currentItem()
        if not item: return
        
        # Column 0 is index
        try:
            mo_idx = int(item.text(0))
        except: return
        
        # Get coeffs
        mo_data = self.parent_dlg.parser.data.get("mo_coeffs", {}).get(mo_idx)
        if not mo_data or not mo_data.get('coeffs'):
            QMessageBox.information(self, "No Info", f"No coefficients found for MO {mo_idx}. \n(Check if 'MOLECULAR ORBITALS' block exists in output)")
            return
            
        dlg = MOCoeffDialog(self, mo_idx, mo_data)
        dlg.exec()

class MOCoeffDialog(QDialog):
    def __init__(self, parent, mo_idx, mo_data):
        super().__init__(parent)
        self.parent_mo_dlg = parent # MODialog
        
        energy_str = f"E={mo_data.get('energy', 0.0):.4f} Eh"
        occ_str = f"Occ={mo_data.get('occ', 0.0):.2f}"
        spin_str = mo_data.get('spin', '').capitalize()
        
        self.setWindowTitle(f"MO {mo_idx} ({spin_str}) - {energy_str}, {occ_str}")
        self.resize(500, 600)
        
        coeffs = mo_data.get('coeffs', [])
        
        layout = QVBoxLayout(self)
        
        # Sort by abs(coeff) descending
        coeffs_sorted = sorted(coeffs, key=lambda x: abs(x['coeff']), reverse=True)
        
        self.table = QTableWidget()
        self.table.setColumnCount(4)
        self.table.setHorizontalHeaderLabels(["Atom", "Sym", "Orb", "Coeff"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.setRowCount(len(coeffs_sorted))
        
        for r, c in enumerate(coeffs_sorted):
            self.table.setItem(r, 0, QTableWidgetItem(str(c['atom_idx'])))
            self.table.setItem(r, 1, QTableWidgetItem(c['sym']))
            self.table.setItem(r, 2, QTableWidgetItem(c['orb']))
            coeff_val = c['coeff']
            item_c = QTableWidgetItem(f"{coeff_val:.5f}")
            # Colorize for visibility
            if coeff_val > 0:
                item_c.setForeground(QColor("blue"))
            else:
                item_c.setForeground(QColor("red"))
            self.table.setItem(r, 3, item_c)
            
        layout.addWidget(self.table)
        
        # Button layout
        btn_layout = QHBoxLayout()
        
        # Link to 3D Viewer if possible
        # Requires access to file path.
        # MODialog -> OrcaResultAnalyzerDialog -> file_path
        
        btn_view_3d = QPushButton("Open FCHK in 3D Viewer")
        btn_view_3d.clicked.connect(self.launch_3d_viewer)
        btn_layout.addWidget(btn_view_3d)

        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        btn_layout.addWidget(btn_close)
        
        layout.addLayout(btn_layout)

    def launch_3d_viewer(self):
        # Access parent's file_path via chain
        # self.parent_mo_dlg is MODialog
        # self.parent_mo_dlg.parent_dlg is OrcaResultAnalyzerDialog
        
        if not hasattr(self.parent_mo_dlg, "parent_dlg"): return
        main_dlg = self.parent_mo_dlg.parent_dlg
        
        if not hasattr(main_dlg, "file_path"):
            return
            
        base_path = os.path.splitext(main_dlg.file_path)[0]
        fchk_path = base_path + ".fchk"
        
        if os.path.exists(fchk_path):
            reply = QMessageBox.question(self, "Open FCHK", 
                                         f"Found {os.path.basename(fchk_path)}.\nOpen it for 3D orbital visualization?",
                                         QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
            
            if reply == QMessageBox.StandardButton.Yes:
                main_dlg.mw.load_file(fchk_path)
                self.close() # Close MO list to avoid clutter
        else:
            QMessageBox.information(self, "Not Found", 
                                    f"Could not find corresponding FCHK file:\n{fchk_path}\n\n"
                                    "To generate it, run: orca_2mkl {base} -molden\n"
                                    "(Note: Gaussian MO Analyzer generates cubes from FCHK)")
