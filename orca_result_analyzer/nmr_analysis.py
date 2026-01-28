import os
import json
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, QComboBox, 
                             QDoubleSpinBox, QTableWidget, QTableWidgetItem, QHeaderView, 
                             QPushButton, QApplication, QGroupBox, QMessageBox,
                              QFileDialog, QCheckBox, QButtonGroup)
from PyQt6.QtCore import Qt, QTimer
import pyvista as pv
import numpy as np

# Import RDKit for VDW radii calculation
try:
    from rdkit import Chem
    _pt = Chem.GetPeriodicTable()
    # Base VDW radii (scaled by 0.3 as in moledit core)
    VDW_RADII = {_pt.GetElementSymbol(i): _pt.GetRvdw(i) * 0.3 for i in range(1, 119)}
except ImportError:
    VDW_RADII = {'H': 1.2 * 0.3, 'C': 1.7 * 0.3, 'N': 1.55 * 0.3, 'O': 1.52 * 0.3}

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT
from matplotlib.figure import Figure
from .nmr_custom_ref_dialog import CustomReferenceDialog


class NMRDialog(QDialog):
    """Enhanced NMR Chemical Shielding Dialog with Spectrum"""
    
    def __init__(self, parent, data, file_path=None):
        super().__init__(parent)
        self.setWindowTitle("NMR Chemical Shielding & Spectrum")
        self.resize(600, 800)  # More compact width
        
        # Make dialog modeless (non-blocking)
        self.setWindowModality(Qt.WindowModality.NonModal)
        
        # Store parent dialog for 3D viewer access
        self.parent_dlg = parent
        
        self.data = data
        self.displayed_data = list(data)
        
        # Track atom labels in 3D viewer
        self._atom_labels = []
        
        # Track selected peaks for highlighting
        self.selected_peak_indices = set()
        self.highlight_artists = []
        self.show_all_mode = False  # Track if showing all labels without highlights
        
        # Reference standards database (delta = reference position, sigma = isotropic shielding)
        # Chemical shift formula: δ_sample = δ_ref + (σ_ref - σ_sample)
        self.reference_standards = {
            "1H": {
        "No Reference": {"delta_ref": 0.0, "sigma_ref": 0.0},
        "TMS": {"delta_ref": 0.0, "sigma_ref": 31.8},
        "CDCl3": {"delta_ref": 7.26, "sigma_ref": 24.5},
        "DMSO-d6": {"delta_ref": 2.50, "sigma_ref": 29.3},
        "Custom": {"delta_ref": 0.0, "sigma_ref": 0.0}
    },
    "13C": {
        "No Reference": {"delta_ref": 0.0, "sigma_ref": 0.0},
        "TMS": {"delta_ref": 0.0, "sigma_ref": 182.4},
        "CDCl3": {"delta_ref": 77.16, "sigma_ref": 105.2},
        "DMSO-d6": {"delta_ref": 39.52, "sigma_ref": 142.9},
        "Custom": {"delta_ref": 0.0, "sigma_ref": 0.0}
    },
    "15N": {
        "No Reference": {"delta_ref": 0.0, "sigma_ref": 0.0},
        "CH3NO2": {"delta_ref": 0.0, "sigma_ref": -135.8},
        "NH3": {"delta_ref": -381.9, "sigma_ref": 244.4},
        "Custom": {"delta_ref": 0.0, "sigma_ref": 0.0}
    },
    "31P": {
        "No Reference": {"delta_ref": 0.0, "sigma_ref": 0.0},
        "H3PO4 (85%)": {"delta_ref": 0.0, "sigma_ref": 328.4},
        "Custom": {"delta_ref": 0.0, "sigma_ref": 0.0}
    },
    "19F": {
        "No Reference": {"delta_ref": 0.0, "sigma_ref": 0.0},
        "CFCl3": {"delta_ref": 0.0, "sigma_ref": 188.5},
        "Custom": {"delta_ref": 0.0, "sigma_ref": 0.0}
    }        }
        
        # Current reference values
        self.delta_ref = 0.0
        self.sigma_ref = 0.0
        
        # Spectrum settings
        self.linewidth = 1.0  # ppm for spectrum
        self.peak_intensity = 1.0
        
        # Settings file
        self.file_path = file_path
        if file_path:
            base_name = os.path.splitext(os.path.basename(file_path))[0]
            file_dir = os.path.dirname(file_path)
            self.merged_peaks_file = os.path.join(file_dir, f"{base_name}-nmr_peak_info.json")
        else:
            # Fallback if no file path provided
            self.merged_peaks_file = os.path.join(os.path.dirname(__file__), "nmr_merged_peaks.json")
        
        self.settings_file = os.path.join(os.path.dirname(__file__), "settings.json")
        
        # Track manually merged peaks: [{"indices": [0, 1, 2], "avg_delta": 7.5, ...}]
        self.merged_peaks = []
        self.load_merged_peaks()
        
        self.load_settings()
        
        # Initialize current nucleus before UI setup
        self.current_nucleus = "All"
        
        self.setup_ui()
        
        # Timer for polling main window selection (sync 3D -> NMR)
        self.sel_timer = QTimer(self)
        self.sel_timer.timeout.connect(self._check_external_selection)
        self.sel_timer.start(200) # Check every 200ms
        
        # Custom 3D highlight actors and names
        self._nmr_sphere_actors = []
        self._nmr_label_names = [] # Explicitly track label names for removal

    def _check_external_selection(self):
        """Poll main window for 3D selection changes"""
        if not hasattr(self.parent_dlg, 'mw'):
            return
            
        mw = self.parent_dlg.mw
        
        # 0. Check if Selection Mode is active (Optimization: Don't hijack selection if user is doing something else)
        # Assuming mw has a 'current_mode' or we check if 'SelectionTool' is active if that structure exists.
        # Fallback: If mw has 'selection_enabled' flag.
        # 0. Check if Selection Mode is active
        # We removed the mw.scene.mode check because it refers to the 2D editor mode.
        # 3D selection should be allowed unless we are in a specific conflicting mode like Measurement.

        if getattr(mw, 'measurement_mode', False):
            # Measurement mode uses a different selection list
            pass

        indices = set()
        
        # Check standard 3D selection
        if hasattr(mw, 'selected_atoms_3d') and mw.selected_atoms_3d:
            indices.update(mw.selected_atoms_3d)
            
        # Check measurement selection
        if hasattr(mw, 'selected_atoms_for_measurement') and mw.selected_atoms_for_measurement:
            for item in mw.selected_atoms_for_measurement:
                if isinstance(item, int):
                    indices.add(item)
                    
        # 1. State Tracking for Stability
        # Compare current 3D selection with what we LAST knew about/set.
        # If they are identical, NO CHANGE has happened, so we do nothing.
        # This prevents the "Echo" loop where we clear the selection (setting it to empty)
        # and then this poller sees "Empty" vs "My Internal Selection" and clears the graph.
        
        current_mw_selection_frozen = frozenset(indices)
        
        # Initialize tracker if missing
        if not hasattr(self, '_last_synced_mw_selection'):
            self._last_synced_mw_selection = frozenset()
            
        # If the 3D selection HAS NOT CHANGED from what we last saw/set, STOP.
        if current_mw_selection_frozen == self._last_synced_mw_selection:
            return

        # Update our tracker to the new state
        self._last_synced_mw_selection = current_mw_selection_frozen

        # Calculate what the NMR selection SHOULD be based on the current 3D selection
        # Note: We use the "Any Member" rule for selecting, but by not expanding the 3D set,
        # unselection of the specific clicked atom correctly clears the indices.
        new_peak_selection = self._calculate_peak_selection_from_atoms(indices)
        
        # Only update if the peak selection itself has changed
        if new_peak_selection != self.selected_peak_indices:
            # If 3D selection is empty, force clear
            if not indices:
                self.clear_peak_selection()
            else:
                self.selected_peak_indices = new_peak_selection
                self.highlight_selected_peaks()
                # Update visual labels and spheres (Yellow), 
                # but tell it NOT to sync back to the main window's selection set (Green).
                self.update_selected_labels(is_external_sync=True)

    def _calculate_peak_selection_from_atoms(self, target_atoms):
        """Helper to determine which peaks should be selected based on atom set"""
        if not hasattr(self, 'peaks_metadata') or not self.peaks_metadata:
            return set()
            
        new_selection = set()
        target_atoms = {int(i) for i in target_atoms}
        
        for peak_idx, metadata in enumerate(self.peaks_metadata):
            # metadata: (shift, intensity, is_merged, atom_indices)
            _, _, _, peak_atoms = metadata
            peak_atoms_set = {int(i) for i in peak_atoms}
            
            # Selection Rule: A peak is selected if ANY of its atoms are in the 3D selection set.
            # This is stable and prevents the "flicker" caused by switching between ANY and SUBSET.
            if not peak_atoms_set.isdisjoint(target_atoms):
                new_selection.add(peak_idx)
        return new_selection

    def select_peaks_by_atom_indices(self, atom_indices):
        """Deprecated/Legacy: Now uses _check_external_selection logic directly"""
        # Kept for potential internal calls, but redirected to robust logic
        new_peaks = self._calculate_peak_selection_from_atoms(atom_indices)
        if new_peaks != self.selected_peak_indices:
            self.selected_peak_indices = new_peaks
            self.highlight_selected_peaks()
            self.update_selected_labels()
    
    def get_nucleus_key(self, atom_sym):
        """Map atom symbol to nucleus key for reference standards"""
        mapping = {
            "H": "1H", "D": "2H", "T": "3H",
            "Li": "7Li", "Be": "9Be", "B": "11B", "C": "13C", "N": "15N", "O": "17O", "F": "19F",
            "Na": "23Na", "Mg": "25Mg", "Al": "27Al", "Si": "29Si", "P": "31P", "S": "33S",
            "Cl": "35Cl", "K": "39K", "Ca": "43Ca", "Sc": "45Sc", "Ti": "47Ti", "V": "51V",
            "Cr": "53Cr", "Mn": "55Mn", "Fe": "57Fe", "Co": "59Co", "Ni": "61Ni", "Cu": "63Cu",
            "Zn": "67Zn", "Ga": "71Ga", "Ge": "73Ge", "As": "75As", "Se": "77Se", "Br": "81Br",
            "Kr": "83Kr", "Rb": "87Rb", "Sr": "87Sr", "Y": "89Y", "Zr": "91Zr", "Nb": "93Nb",
            "Mo": "95Mo", "Tc": "99Tc", "Ru": "99Ru", "Rh": "103Rh", "Pd": "105Pd", "Ag": "109Ag",
            "Cd": "113Cd", "In": "115In", "Sn": "119Sn", "Sb": "121Sb", "Te": "125Te", "I": "127I",
            "Xe": "129Xe", "Cs": "133Cs", "Ba": "137Ba", "La": "139La", "W": "183W", "Os": "187Os",
            "Pt": "195Pt", "Au": "197Au", "Hg": "199Hg", "Tl": "205Tl", "Pb": "207Pb"
        }
        return mapping.get(atom_sym, atom_sym)
        
    def load_settings(self):
        """Load NMR settings from JSON"""
        if os.path.exists(self.settings_file):
            try:
                with open(self.settings_file, 'r') as f:
                    settings = json.load(f)
                
                nmr_settings = settings.get("nmr_settings", {})
                self.linewidth = nmr_settings.get("spectrum_linewidth", 1.0)
                self.peak_intensity = nmr_settings.get("peak_intensity", 1.0)
                
                # Load custom references
                custom_refs = nmr_settings.get("custom_references", {})
                for nucleus, refs in custom_refs.items():
                    if nucleus not in self.reference_standards:
                        self.reference_standards[nucleus] = {}
                    for ref_name, ref_val in refs.items():
                        self.reference_standards[nucleus][ref_name] = ref_val
            except Exception as e:
                print(f"Error loading NMR settings: {e}")
    
    def merge_selected_peaks(self):
        """Merge selected peaks into a single entry"""
        if len(self.selected_peak_indices) < 2:
            QMessageBox.warning(self, "Invalid Selection", "Please select at least 2 peaks to merge.")
            return
        
        # Get selected atom indices
        selected_indices = []
        
        # Use peaks_metadata which maps plot indices to atoms/groups
        if hasattr(self, 'peaks_metadata') and self.peaks_metadata:
            for peak_idx in self.selected_peak_indices:
                if peak_idx < len(self.peaks_metadata):
                    # Metadata format: (shift, intensity, is_merged, atom_indices)
                    _, _, _, atom_indices = self.peaks_metadata[peak_idx]
                    selected_indices.extend(atom_indices)
        else:
            # Fallback (though this shouldn't happen if selection exists and spectrum plotted)
            for peak_idx in sorted(self.selected_peak_indices):
                if peak_idx < len(self.displayed_data):
                    atom_idx = self.displayed_data[peak_idx].get('atom_idx', peak_idx)
                    selected_indices.append(atom_idx)
        
        # Ensure unique indices
        selected_indices = sorted(list(set(selected_indices)))
        
        # Calculate averages (for display/confirmation only - not saved)
        total_sigma = 0.0
        total_delta = 0.0
        for atom_idx in selected_indices:
            # Find item in original data
            item = next((d for d in self.data if d.get('atom_idx') == atom_idx), None)
            if item:
                sigma = item.get('shielding', 0.0)
                delta = self.delta_ref + (self.sigma_ref - sigma)
                total_sigma += sigma
                total_delta += delta
        
        count = len(selected_indices)
        avg_sigma = total_sigma / count
        avg_delta = total_delta / count
        
        # Create merge group - ONLY save indices, not averages
        merge_group = {
            'indices': selected_indices
        }
        
        # Check if any indices are already merged
        for existing_group in self.merged_peaks:
            if any(idx in existing_group['indices'] for idx in selected_indices):
                reply = QMessageBox.question(
                    self, "Already Merged",
                    "Some selected peaks are already in a merged group. Replace existing merge?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
                )
                if reply == QMessageBox.StandardButton.Yes:
                    # Remove old group
                    self.merged_peaks.remove(existing_group)
                else:
                    return
        
        # Add new merge group
        self.merged_peaks.append(merge_group)
        
        # Save to JSON
        self.save_merged_peaks()
        
        # Clear selection and refresh
        self.clear_peak_selection()
        
        # Also clear 3D selection in main window
        if hasattr(self.parent_dlg, 'mw'):
            mw = self.parent_dlg.mw
            if hasattr(mw, 'selected_atoms_3d'):
                mw.selected_atoms_3d.clear()
            if hasattr(mw, 'update_3d_selection_display'):
                mw.update_3d_selection_display()
                
        self.recalc()
    
    def unmerge_selected_peaks(self):
        """Separate previously merged peaks back into individuals"""
        if not self.selected_peak_indices:
            return
            
        groups_to_remove = []
        for peak_idx in self.selected_peak_indices:
            if peak_idx < len(self.peaks_metadata):
                _, _, is_merged, atom_indices = self.peaks_metadata[peak_idx]
                if is_merged:
                    # Find which group in self.merged_peaks contains these indices
                    # Since groups are unique per atom, we can just match any index
                    for group in self.merged_peaks:
                        if all(idx in group['indices'] for idx in atom_indices) and len(group['indices']) == len(atom_indices):
                            if group not in groups_to_remove:
                                groups_to_remove.append(group)
                            break
        
        if not groups_to_remove:
            return
            
        for group in groups_to_remove:
            self.merged_peaks.remove(group)
            
        self.save_merged_peaks()
        self.clear_peak_selection()
        self.recalc()
    
    def save_merged_peaks(self):
        """Save merged peaks to JSON file"""
        try:
            with open(self.merged_peaks_file, 'w') as f:
                json.dump(self.merged_peaks, f, indent=2)
        except Exception as e:
            print(f"Error saving merged peaks: {e}")
    
    def load_merged_peaks(self):
        """Load merged peaks from JSON file"""
        if os.path.exists(self.merged_peaks_file):
            try:
                with open(self.merged_peaks_file, 'r') as f:
                    self.merged_peaks = json.load(f)
            except Exception as e:
                print(f"Error loading merged peaks: {e}")
                self.merged_peaks = []
        else:
            self.merged_peaks = []
    
    def save_settings(self):
        """Save NMR settings to JSON"""
        settings = {}
        if os.path.exists(self.settings_file):
            try:
                with open(self.settings_file, 'r') as f:
                    settings = json.load(f)
            except:
                pass
        
        # Extract custom references only (non-default)
        # Extract custom references only (non-default)
        default_standards = {
            "1H": {
                "No Reference": {"delta_ref": 0.0, "sigma_ref": 0.0},
                "TMS": {"delta_ref": 0.0, "sigma_ref": 31.8},
                "CDCl3": {"delta_ref": 7.26, "sigma_ref": 24.5},
                "DMSO-d6": {"delta_ref": 2.50, "sigma_ref": 29.3}
            },
            "13C": {
                "No Reference": {"delta_ref": 0.0, "sigma_ref": 0.0},
                "TMS": {"delta_ref": 0.0, "sigma_ref": 182.4},
                "CDCl3": {"delta_ref": 77.16, "sigma_ref": 105.2},
                "DMSO-d6": {"delta_ref": 39.52, "sigma_ref": 142.9}
            },
            "15N": {
                "No Reference": {"delta_ref": 0.0, "sigma_ref": 0.0},
                "CH3NO2": {"delta_ref": 0.0, "sigma_ref": -135.8},
                "NH3": {"delta_ref": -381.9, "sigma_ref": 244.4}
            },
            "31P": {
                "No Reference": {"delta_ref": 0.0, "sigma_ref": 0.0},
                "H3PO4 (85%)": {"delta_ref": 0.0, "sigma_ref": 328.4}
            },
            "19F": {
                "No Reference": {"delta_ref": 0.0, "sigma_ref": 0.0},
                "CFCl3": {"delta_ref": 0.0, "sigma_ref": 188.5}
            }
        }
        
        custom_refs = {}
        for nucleus, refs in self.reference_standards.items():
            for ref_name, ref_val in refs.items():
                if ref_name == "Custom":
                    continue  # Skip custom placeholder
                default_dict = default_standards.get(nucleus, {})
                if ref_name not in default_dict:
                    # Completely custom reference
                    if nucleus not in custom_refs:
                        custom_refs[nucleus] = {}
                    custom_refs[nucleus][ref_name] = ref_val
                elif default_dict.get(ref_name) != ref_val:
                    # Modified default reference
                    if nucleus not in custom_refs:
                        custom_refs[nucleus] = {}
                    custom_refs[nucleus][ref_name] = ref_val
        
        settings["nmr_settings"] = {
            "custom_references": custom_refs
        }
        
        try:
            with open(self.settings_file, 'w') as f:
                json.dump(settings, f, indent=2)
        except Exception as e:
            print(f"Error saving NMR settings: {e}")
    
    def setup_ui(self):
        main_layout = QVBoxLayout(self)
        
        # 1. Reference & Element Filter Row
        top_row = QHBoxLayout()
        
        # Nucleus filter with toggle buttons
        nucleus_box = QGroupBox("Nucleus Filter")
        nucleus_layout = QHBoxLayout(nucleus_box)
        
        # Get available nuclei
        nuclei = ["All"] + sorted(list(set([d['atom_sym'] for d in self.data])))
        
        # Create button group for exclusive selection
        self.nucleus_button_group = QButtonGroup()
        self.nucleus_buttons = {}
        
        for nucleus in nuclei:
            btn = QPushButton(nucleus)
            btn.setCheckable(True)
            btn.setStyleSheet("""
                QPushButton {
                    background-color: #f8f8f8;
                    border: 1px solid #ccc;
                    border-radius: 4px;
                    padding: 6px 15px;
                    font-size: 10pt;
                    font-weight: 500;
                }
                QPushButton:checked {
                    background-color: #0066cc;
                    color: white;
                    border-color: #004d99;
                    font-weight: bold;
                }
                QPushButton:hover {
                    border-color: #0066cc;
                    background-color: #eef6ff;
                }
                QPushButton:checked:hover {
                    background-color: #0059b3;
                }
            """)
            btn.toggled.connect(lambda checked, n=nucleus: self.on_nucleus_changed(n) if checked else None)
            self.nucleus_button_group.addButton(btn)
            self.nucleus_buttons[nucleus] = btn
            nucleus_layout.addWidget(btn)
        
        
        nucleus_layout.addStretch()
        top_row.addWidget(nucleus_box)
        
        # Reference selection and values (merged into single group)
        ref_box = QGroupBox("Reference Standard")
        ref_layout = QVBoxLayout(ref_box)
        
        # Reference selection row
        ref_sel_row = QHBoxLayout()
        ref_sel_row.addWidget(QLabel("Standard:"))
        self.combo_ref = QComboBox()
        self.combo_ref.setMinimumWidth(250)
        from PyQt6.QtWidgets import QSizePolicy
        self.combo_ref.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.combo_ref.currentIndexChanged.connect(self.on_ref_change)
        ref_sel_row.addWidget(self.combo_ref)
        
        btn_add_ref = QPushButton("+ Custom")
        btn_add_ref.setFixedWidth(80)
        btn_add_ref.clicked.connect(self.add_custom_reference)
        ref_sel_row.addWidget(btn_add_ref)
        
        btn_del_ref = QPushButton("Delete")
        btn_del_ref.setFixedWidth(60)
        btn_del_ref.setToolTip("Delete selected custom reference")
        btn_del_ref.clicked.connect(self.delete_custom_reference)
        ref_sel_row.addWidget(btn_del_ref)
        
        # Duplicate combo_ref removed
        
        ref_layout.addLayout(ref_sel_row)
        
        # Delta ref (reference peak position)
        delta_row = QHBoxLayout()
        delta_row.addWidget(QLabel("δ_ref (ppm):"))
        self.spin_delta_ref = QDoubleSpinBox()
        self.spin_delta_ref.setRange(-500, 500)
        self.spin_delta_ref.setValue(self.delta_ref)
        self.spin_delta_ref.setDecimals(2)
        self.spin_delta_ref.setToolTip("Reference peak position on chemical shift scale")
        self.spin_delta_ref.valueChanged.connect(self.on_ref_value_change)
        delta_row.addWidget(self.spin_delta_ref)
        ref_layout.addLayout(delta_row)
        
        # Sigma ref (isotropic shielding)
        sigma_row = QHBoxLayout()
        sigma_row.addWidget(QLabel("σ_ref (ppm):"))
        self.spin_sigma_ref = QDoubleSpinBox()
        self.spin_sigma_ref.setRange(-10000, 20000)
        self.spin_sigma_ref.setValue(self.sigma_ref)
        self.spin_sigma_ref.setDecimals(2)
        self.spin_sigma_ref.setToolTip("Isotropic shielding value of reference compound")
        self.spin_sigma_ref.valueChanged.connect(self.on_ref_value_change)
        sigma_row.addWidget(self.spin_sigma_ref)
        ref_layout.addLayout(sigma_row)
        
        top_row.addWidget(ref_box)
        
        main_layout.addLayout(top_row)
        
        # 2. Spectrum Plot
        spec_group = QGroupBox("NMR Stick Spectrum (δ = δ_ref + σ_ref - σ)")
        spec_layout = QVBoxLayout(spec_group)
        
        # Spectrum settings row
        spec_settings = QHBoxLayout()
        
        # Add checkbox to show all labels
        self.chk_show_all_labels = QCheckBox("Show all atom labels")
        self.chk_show_all_labels.stateChanged.connect(self.toggle_all_labels)
        spec_settings.addWidget(self.chk_show_all_labels)
        
        # Add button to clear selection
        btn_clear_selection = QPushButton("Clear Selection")
        btn_clear_selection.setFixedWidth(120)
        btn_clear_selection.setToolTip("Clear all selected peaks and labels")
        btn_clear_selection.clicked.connect(self.clear_peak_selection)
        spec_settings.addWidget(btn_clear_selection)
        
        # Add merge selected button
        btn_merge = QPushButton("Merge Selected")
        btn_merge.setFixedWidth(120)
        btn_merge.setToolTip("Merge selected peaks into one entry")
        btn_merge.clicked.connect(self.merge_selected_peaks)
        spec_settings.addWidget(btn_merge)
        
        btn_unmerge = QPushButton("Unmerge")
        btn_unmerge.setFixedWidth(120)
        btn_unmerge.setToolTip("Separate merged peaks back to individuals")
        btn_unmerge.clicked.connect(self.unmerge_selected_peaks)
        spec_settings.addWidget(btn_unmerge)
        
        spec_settings.addStretch()
        
        btn_export = QPushButton("Export Spectrum")
        btn_export.clicked.connect(self.export_spectrum)
        spec_settings.addWidget(btn_export)
        
        spec_layout.addLayout(spec_settings)
        
        # Matplotlib canvas - adjusted for narrower dialog
        self.figure = Figure(figsize=(5.5, 4))  # Narrower to fit 600px width
        self.canvas = FigureCanvas(self.figure)
        
        # Add navigation toolbar for zoom/pan
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        spec_layout.addWidget(self.toolbar)
        spec_layout.addWidget(self.canvas)
        
        main_layout.addWidget(spec_group)
        
        # 3. Data Table
        table_group = QGroupBox("Chemical Shift Data")
        table_layout = QVBoxLayout(table_group)
        
        self.table = QTableWidget()
        self.table.setColumnCount(4)
        self.table.setHorizontalHeaderLabels(["Idx", "Nucleus", "σ (ppm)", "δ (ppm)"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.verticalHeader().setVisible(False)
        self.table.setMinimumHeight(100)  # Half size for compact layout
        table_layout.addWidget(self.table)
        
        table_btn_row = QHBoxLayout()
        btn_copy = QPushButton("Copy to Clipboard")
        btn_copy.clicked.connect(self.copy_table)
        table_btn_row.addWidget(btn_copy)
        table_btn_row.addStretch()
        table_layout.addLayout(table_btn_row)
        
        main_layout.addWidget(table_group)
        
        # 4. Bottom Buttons
        bottom_row = QHBoxLayout()
        bottom_row.addStretch()
        btn_close = QPushButton("Close")
        btn_close.setFixedWidth(100)
        btn_close.clicked.connect(self.close)
        bottom_row.addWidget(btn_close)
        main_layout.addLayout(bottom_row)
        
        # Finally select "All" by default to trigger population
        if "All" in self.nucleus_buttons:
            self.nucleus_buttons["All"].setChecked(True)
    
    def update_reference_combo(self):
        """Update reference combo box based on current nucleus"""
        if not hasattr(self, 'combo_ref') or self.combo_ref is None:
            return
            
        current_nucleus = self.current_nucleus
        if current_nucleus == "All":
            # Special case for "All": Only allow "No Reference"
            # We want to visualize raw values (no shifting) across different nuclei
            self.combo_ref.blockSignals(True)
            self.combo_ref.clear()
            self.combo_ref.addItems(["No Reference"])
            self.combo_ref.setCurrentText("No Reference")
            self.combo_ref.blockSignals(False)
            
            # Set values to 0,0
            self.delta_ref = 0.0
            self.sigma_ref = 0.0
            
            # Update spinboxes
            self.spin_delta_ref.blockSignals(True)
            self.spin_sigma_ref.blockSignals(True)
            self.spin_delta_ref.setValue(0.0)
            self.spin_sigma_ref.setValue(0.0)
            self.spin_delta_ref.blockSignals(False)
            self.spin_sigma_ref.blockSignals(False)
            
            # Disable inputs
            self.spin_delta_ref.setEnabled(False)
            self.spin_sigma_ref.setEnabled(False)
            return

        # Map atom symbol to nucleus key (e.g., "H" -> "1H")
        current_nucleus = self.get_nucleus_key(current_nucleus)
        
        # Save current selection to preserve it if possible
        current_ref = self.combo_ref.currentText() if self.combo_ref.count() > 0 else None
        
        # Block signals to prevent triggering on_ref_change during population
        self.combo_ref.blockSignals(True)
        
        # Clear and repopulate
        self.combo_ref.clear()
        refs = self.reference_standards.get(current_nucleus, {})
        
        # Create list of items
        items = list(refs.keys())
        
        # Always ensure "Custom" is in the list
        if "Custom" not in items:
            items.append("Custom")
            # Ensure it exists in logic too (transiently if not in storage)
            if "Custom" not in refs:
                 # We don't save this to self.reference_standards to avoid polluting it with empty Customs,
                 # but we need to handle its selection.
                 pass

        self.combo_ref.addItems(items)
        
        # Try to restore previous selection, otherwise set default to TMS or first item
        if current_ref and current_ref in refs:
            self.combo_ref.setCurrentText(current_ref)
            ref_data = refs[current_ref]
        elif current_ref == "Custom":
            # Transient custom selected
            self.combo_ref.setCurrentText("Custom")
            ref_data = {"delta_ref": 0.0, "sigma_ref": 0.0}
        elif "TMS" in refs:
            self.combo_ref.setCurrentText("TMS")
            ref_data = refs["TMS"]
        elif refs:
            first_ref = list(refs.keys())[0]
            self.combo_ref.setCurrentText(first_ref)
            ref_data = refs[first_ref]
        else:
            # Should be Custom default
            self.combo_ref.setCurrentText("Custom")
            ref_data = {"delta_ref": 0.0, "sigma_ref": 0.0}
        
        # Update internal values
        self.delta_ref = ref_data.get("delta_ref", 0.0)
        self.sigma_ref = ref_data.get("sigma_ref", 0.0)
        
        # Update spinboxes (these updates won't trigger recalc since we blocked combo signals)
        self.spin_delta_ref.blockSignals(True)
        self.spin_sigma_ref.blockSignals(True)
        self.spin_delta_ref.setValue(self.delta_ref)
        self.spin_sigma_ref.setValue(self.sigma_ref)
        self.spin_delta_ref.blockSignals(False)
        self.spin_sigma_ref.blockSignals(False)

        # Only allow editing if "Custom" is selected
        is_custom = (self.combo_ref.currentText() == "Custom")
        self.spin_delta_ref.setEnabled(is_custom)
        self.spin_sigma_ref.setEnabled(is_custom)
        
        # Re-enable combo signals
        self.combo_ref.blockSignals(False)
    
    def on_ref_change(self):
        """Handle reference standard change"""
        current_nucleus = self.current_nucleus
        if current_nucleus == "All":
            # Force 0,0 for All view
            self.delta_ref = 0.0
            self.sigma_ref = 0.0
            self.recalc()
            return
        
        # Map atom symbol to nucleus key
        nucleus_key = self.get_nucleus_key(current_nucleus)
        
        ref_name = self.combo_ref.currentText()
        refs = self.reference_standards.get(nucleus_key, {})
        ref_data = refs.get(ref_name, {"delta_ref": 0.0, "sigma_ref": 0.0})
        
        self.delta_ref = ref_data["delta_ref"]
        self.sigma_ref = ref_data["sigma_ref"]
        
        # Block signals only during setValue to prevent triggering on_ref_value_change
        self.spin_delta_ref.blockSignals(True)
        self.spin_sigma_ref.blockSignals(True)
        self.spin_delta_ref.setValue(self.delta_ref)
        self.spin_sigma_ref.setValue(self.sigma_ref)
        self.spin_delta_ref.blockSignals(False)
        self.spin_sigma_ref.blockSignals(False)
        
        # Only allow editing if "Custom" is selected
        is_custom = (ref_name == "Custom")
        self.spin_delta_ref.setEnabled(is_custom)
        self.spin_sigma_ref.setEnabled(is_custom)
        
        self.recalc()
    
    def on_ref_value_change(self):
        """Handle manual reference value changes"""
        self.delta_ref = self.spin_delta_ref.value()
        self.sigma_ref = self.spin_sigma_ref.value()
        
        # Update the reference standard dict
        current_nucleus = self.current_nucleus
        ref_name = self.combo_ref.currentText()
        if current_nucleus != "All":
            # Map atom symbol to nucleus key
            nucleus_key = self.get_nucleus_key(current_nucleus)
            if nucleus_key not in self.reference_standards:
                self.reference_standards[nucleus_key] = {}
            self.reference_standards[nucleus_key][ref_name] = {
                "delta_ref": self.delta_ref,
                "sigma_ref": self.sigma_ref
            }
        self.recalc()
    
    def add_custom_reference(self):
        """Add a custom reference standard using custom dialog"""
        # Get available nuclei from data
        available_nuclei = sorted(list(set([d['atom_sym'] for d in self.data])))
        if not available_nuclei:
            available_nuclei = ["1H", "13C", "15N", "31P", "19F"]
        
        # Show custom dialog
        dialog = CustomReferenceDialog(self, available_nuclei)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        
        ref_name, nucleus_data = dialog.get_reference_data()
        
        # Add to standards for each nucleus
        for raw_nucleus, values in nucleus_data.items():
            # Ensure we use the standardized key (e.g. "H" -> "1H")
            nucleus = self.get_nucleus_key(raw_nucleus)
            
            if nucleus not in self.reference_standards:
                self.reference_standards[nucleus] = {}
            self.reference_standards[nucleus][ref_name] = values
        
        self.save_settings()
        
        # If we are in "All" mode, switch to the specific nucleus of the added reference
        # This prevents the reference from being hidden by the "All" restriction ("No Reference" only)
        if self.current_nucleus == "All" and nucleus_data:
             # Pick the first nucleus added (e.g. "1H" or "H")
             first_raw = list(nucleus_data.keys())[0]
             # Get the key used for buttons (e.g. "H" or "1H") - keys in nucleus_buttons usually match atom symbols
             # We should try to find the matching button.
             target_btn_key = None
             
             # Try standardized key first
             std_key = self.get_nucleus_key(first_raw) 
             if std_key in self.nucleus_buttons:
                 target_btn_key = std_key
             elif first_raw in self.nucleus_buttons:
                 target_btn_key = first_raw
             
             if target_btn_key:
                 self.nucleus_buttons[target_btn_key].setChecked(True)
                 # Manually trigger the mode change handler since setChecked might not if via code (depends on signal)
                 # But usually the group handles it. Let's explicitly call the update to be safe.
                 self.current_nucleus = target_btn_key
        
        # Force UI update
        self.update_reference_combo()
        
        # Select the new reference
        # Now that we switched inputs (if we were in All), this should work
        self.combo_ref.setCurrentText(ref_name)
        
        QMessageBox.information(self, "Success", 
                              f"Added reference '{ref_name}' for {len(nucleus_data)} nucleus/nuclei.")

    def delete_custom_reference(self):
        """Delete currently selected custom reference"""
        if not hasattr(self, 'current_nucleus') or not self.current_nucleus:
            return

        current_nucleus = self.get_nucleus_key(self.current_nucleus) if self.current_nucleus != "All" else "1H"
        current_ref = self.combo_ref.currentText()
        
        if not current_ref:
            return

        # Check if it is a built-in standard
        default_standards = {
            "1H": ["TMS", "CDCl3", "DMSO-d6"],
            "13C": ["TMS", "CDCl3", "DMSO-d6"],
            "15N": ["CH3NO2", "NH3"],
            "31P": ["H3PO4 (85%)"],
            "19F": ["CFCl3"]
        }
        
        # Also check hardcoded defaults in save_settings to be safe
        is_default = False
        if current_nucleus in default_standards and current_ref in default_standards[current_nucleus]:
            is_default = True
        
        # Prevent deletion of "Custom" placeholder AND "No Reference"
        if current_ref in ["Custom", "No Reference"]:
            is_default = True
            
        if is_default:
            QMessageBox.warning(self, "Cannot Delete", 
                              f"'{current_ref}' is a built-in standard and cannot be deleted.")
            return

        # Confirm deletion
        reply = QMessageBox.question(self, "Confirm Deletion",
                                   f"Are you sure you want to delete custom reference '{current_ref}'?",
                                   QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
        
        if reply == QMessageBox.StandardButton.Yes:
            # Remove from storage
            if current_nucleus in self.reference_standards:
                if current_ref in self.reference_standards[current_nucleus]:
                    del self.reference_standards[current_nucleus][current_ref]
                    self.save_settings()
                    self.update_reference_combo()
                    QMessageBox.information(self, "Deleted", f"Reference '{current_ref}' removed.")
                else:
                    # Should not happen if UI is consistent
                    QMessageBox.warning(self, "Error", "Reference not found in storage.")
    
    def on_spectrum_settings_change(self):
        """Handle spectrum visualization parameter changes"""
        self.linewidth = self.spin_linewidth.value()
        self.peak_intensity = self.spin_intensity.value()
        self.save_settings()
        self.plot_spectrum()
    
    def on_nucleus_changed(self, nucleus):
        """Handle nucleus button toggle"""
        self.current_nucleus = nucleus
        self.apply_filter()
    
    def apply_filter(self):
        """Filter data by nucleus and update UI"""
        if not hasattr(self, 'table') or self.table is None:
            return
            
        nucleus = self.current_nucleus
        if nucleus == "All":
            self.displayed_data = self.data
        else:
            self.displayed_data = [d for d in self.data if d['atom_sym'] == nucleus]
        
        # Update reference combo for this nucleus
        self.update_reference_combo()
        
        # Force recalculation with new reference values
        self.recalc()
    
    def recalc(self):
        """Recalculate chemical shifts and update table + spectrum"""
        # Get merged peak groups
        merged_indices = set()
        for group in self.merged_peaks:
            merged_indices.update(group['indices'])
        
        # Track which rows to display
        rows_to_display = []
        
        # Add merged peaks first
        for group in self.merged_peaks:
            # Calculate averages dynamically based on current reference
            total_sigma = 0.0
            total_delta = 0.0
            count = 0
            
            for atom_idx in group['indices']:
                item = next((d for d in self.data if d.get('atom_idx') == atom_idx), None)
                if item:
                    sigma = item.get('shielding', 0.0)
                    delta = self.delta_ref + (self.sigma_ref - sigma)
                    total_sigma += sigma
                    total_delta += delta
                    count += 1
            
            if count > 0:
                avg_sigma = total_sigma / count
                avg_delta = total_delta / count
                
                # Check if any atoms in this group are in displayed_data
                group_items = [item for item in self.displayed_data if item.get('atom_idx', -1) in group['indices']]
                if group_items:
                    rows_to_display.append(('merged', {
                        'indices': group['indices'],
                        'avg_sigma': avg_sigma,
                        'avg_delta': avg_delta
                    }, group_items))
        
        # Add individual non-merged items
        for item in self.displayed_data:
            if item.get('atom_idx', -1) not in merged_indices:
                rows_to_display.append(('individual', item, None))
        
        # Update table
        self.table.setRowCount(len(rows_to_display))
        
        for r, row_data in enumerate(rows_to_display):
            row_type, data, group_items = row_data
            
            if row_type == 'merged':
                # Display merged group
                indices_str = ", ".join([str(idx) for idx in data['indices']])
                self.table.setItem(r, 0, QTableWidgetItem(f"[{indices_str}]"))
                self.table.setItem(r, 1, QTableWidgetItem(f"{len(data['indices'])}{group_items[0].get('atom_sym', '')}"))
                self.table.setItem(r, 2, QTableWidgetItem(f"{data['avg_sigma']:.2f}"))
                self.table.setItem(r, 3, QTableWidgetItem(f"{data['avg_delta']:.2f}"))
            else:
                # Display individual item
                self.table.setItem(r, 0, QTableWidgetItem(str(data.get("atom_idx", ""))))
                self.table.setItem(r, 1, QTableWidgetItem(data.get("atom_sym", "")))
                
                sigma_sample = data.get("shielding", 0.0)
                self.table.setItem(r, 2, QTableWidgetItem(f"{sigma_sample:.2f}"))
                
                # Chemical shift: δ = δ_ref + (σ_ref - σ_sample)
                delta_sample = self.delta_ref + (self.sigma_ref - sigma_sample)
                self.table.setItem(r, 3, QTableWidgetItem(f"{delta_sample:.2f}"))
        
        self.plot_spectrum()
    
    def plot_spectrum(self):
        """Plot NMR stick spectrum using stem plot for better visibility"""
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        
        if not self.displayed_data:
            ax.text(0.5, 0.5, "No data to display",
                   ha='center', va='center', transform=ax.transAxes, fontsize=14)
            self.canvas.draw()
            return
        
        # Build list of peaks to plot (merged + individual)
        peaks_to_plot = []  # List of (chemical_shift, intensity, is_merged, indices)
        
        # Collect merged indices
        merged_indices_set = set()
        for group in self.merged_peaks:
            merged_indices_set.update(group['indices'])
        
        # Add merged peaks first
        for group in self.merged_peaks:
            # Calculate averages dynamically based on current reference
            total_sigma = 0.0
            total_delta = 0.0
            count = 0
            
            for atom_idx in group['indices']:
                item = next((d for d in self.data if d.get('atom_idx') == atom_idx), None)
                if item:
                    sigma = item.get('shielding', 0.0)
                    delta = self.delta_ref + (self.sigma_ref - sigma)
                    total_sigma += sigma
                    total_delta += delta
                    count += 1
            
            if count > 0:
                avg_delta = total_delta / count
                
                # Check if any atoms in this group are in displayed_data
                group_items = [item for item in self.displayed_data if item.get('atom_idx', -1) in group['indices']]
                if group_items:
                    intensity = len(group['indices'])  # Integration
                    peaks_to_plot.append((avg_delta, intensity, True, group['indices']))
        
        # Add individual non-merged peaks
        for item in self.displayed_data:
            atom_idx = item.get('atom_idx', -1)
            if atom_idx not in merged_indices_set:
                sigma_sample = item.get("shielding", 0.0)
                delta_sample = self.delta_ref + (self.sigma_ref - sigma_sample)
                peaks_to_plot.append((delta_sample, 1.0, False, [atom_idx]))
        
        if not peaks_to_plot:
            ax.text(0.5, 0.5, "No data", ha='center', va='center', transform=ax.transAxes)
            self.canvas.draw()
            return
        
        # Extract shifts and intensities
        shifts = [p[0] for p in peaks_to_plot]
        intensities = [p[1] for p in peaks_to_plot]
        
        # Store peak metadata for click handling (shift, intensity, is_merged, atom_indices)
        self.peaks_metadata = peaks_to_plot
        
        # Store shifts for click detection
        self.current_shifts = shifts
        
        # Generate x-axis range
        min_shift = min(shifts)
        max_shift = max(shifts)
        padding = max(0.5, (max_shift - min_shift) * 0.15) if max_shift > min_shift else 5.0
        
        # Adjust y-axis limit based on max intensity
        max_intensity = max(intensities) if intensities else 1.0
        y_limit = max(1.5, max_intensity * 1.3)  # Add 30% headroom
        
        # Plot stick spectrum using stem plot
        markerline, stemlines, baseline = ax.stem(shifts, intensities, 
                                                   linefmt='b-', markerfmt='None', 
                                                   basefmt='k-')
        # Style the stem plot (no markers, just stems)
        stemlines.set_linewidth(2.5)
        stemlines.set_alpha(0.8)
        baseline.set_linewidth(1)
        baseline.set_alpha(0.3)
        
        # Set axis limits - larger y-axis for traditional NMR appearance
        ax.set_xlim(max_shift + padding, min_shift - padding)  # Inverted for NMR convention
        ax.set_ylim(0, y_limit)  # Dynamic limit based on integration
        
        # Formatting
        ax.set_xlabel('Chemical Shift δ (ppm)', fontsize=11, fontweight='bold')
        ax.set_ylabel('Intensity (normalized)', fontsize=11, fontweight='bold')
        ax.set_yticks([0, 0.5, 1.0])
        ax.grid(True, alpha=0.25, linestyle='--', axis='x', linewidth=0.8)
        ax.grid(True, alpha=0.15, linestyle=':', axis='y', linewidth=0.5)
        
        # Title with nucleus and reference info
        current_nucleus = self.current_nucleus
        ref_name = self.combo_ref.currentText() if hasattr(self, 'combo_ref') else "Custom"
        ax.set_title(f'{current_nucleus} NMR Stick Spectrum\n'
                    f'Ref: {ref_name} (δ_ref = {self.delta_ref:.2f} ppm, σ_ref = {self.sigma_ref:.1f} ppm)', 
                    fontsize=11, fontweight='bold')
        
        # Add click event for peak selection
        self.canvas.mpl_connect('button_press_event', self.on_peak_click)
        
        self.figure.tight_layout()
        self.canvas.draw()
    
    def highlight_selected_peaks(self):
        """Add red highlights and labels to selected peaks"""
        if not hasattr(self, 'current_shifts') or not self.selected_peak_indices:
            # Clear highlights if no selection
            for artist in self.highlight_artists:
                try:
                    artist.remove()
                except:
                    pass
            self.highlight_artists = []
            self.canvas.draw()
            return
        
        ax = self.figure.axes[0] if self.figure.axes else None
        if not ax:
            return
        
        # Clear old highlights and labels
        for artist in self.highlight_artists:
            try:
                artist.remove()
            except:
                pass
        self.highlight_artists = []
        
        # Add red highlights and text labels for selected peaks
        for idx in self.selected_peak_indices:
            if idx < len(self.current_shifts):
                shift = self.current_shifts[idx]
                
                # Draw red line over selected peak (only if not in show_all_mode)
                if not self.show_all_mode:
                    line = ax.axvline(shift, ymin=0, ymax=1, color='red', linewidth=3.5, alpha=0.7, zorder=10)
                    self.highlight_artists.append(line)
                
                # Add text label above the peak
                if idx < len(self.peaks_metadata):
                    # Get peak metadata
                    _, _, is_merged, atom_indices = self.peaks_metadata[idx]
                    
                    # Build label text from all atoms in this peak
                    label_parts = []
                    for atom_idx in atom_indices:
                        atom_item = next((d for d in self.data if d.get('atom_idx') == atom_idx), None)
                        if atom_item:
                            atom_sym = atom_item.get('atom_sym', '?')
                            label_parts.append(f"{atom_sym}{atom_idx}")
                    
                    label_text = ",".join(label_parts) if label_parts else "?"
                    
                    # Position label above the peak
                    label_color = 'black' if self.show_all_mode else 'red'
                    text = ax.text(shift, 1.1, label_text, 
                                  ha='center', va='bottom',
                                  fontsize=10, fontweight='bold',
                                  color=label_color, zorder=11)
                    self.highlight_artists.append(text)
        
        self.canvas.draw()
    
    def export_spectrum(self):
        """Export spectrum to file"""
        current_nucleus = self.current_nucleus
        default_name = f"nmr_spectrum_{current_nucleus}.png"
        
        filename, _ = QFileDialog.getSaveFileName(
            self, "Export Spectrum", default_name,
            "PNG Image (*.png);;PDF (*.pdf);;SVG (*.svg)"
        )
        
        if filename:
            try:
                self.figure.savefig(filename, dpi=300, bbox_inches='tight')
                QMessageBox.information(self, "Success", f"Spectrum exported to:\n{os.path.basename(filename)}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Export failed:\n{e}")
    
    def copy_table(self):
        """Copy table data to clipboard"""
        text = "Idx\tNucleus\tShielding\tShift\n"
        for r in range(self.table.rowCount()):
            cols = []
            for c in range(self.table.columnCount()):
                it = self.table.item(r, c)
                cols.append(it.text() if it else "")
            text += "\t".join(cols) + "\n"
        QApplication.clipboard().setText(text)
        QMessageBox.information(self, "Copied", "Table data copied to clipboard!")
    
    def on_peak_click(self, event):
        """Handle click on spectrum peak"""
        if event.inaxes is None or not hasattr(self, 'current_shifts'):
            return
        
        # Disable selection when show all mode is active
        if self.show_all_mode:
            return

        # Remove blocking check for measurement mode
        # User confirmed graph should be clickable and syncing should work.
        
        # Get click position
        click_x = event.xdata
        if click_x is None:
            return
        
        # Find nearest peak within tolerance
        tolerance = 0.5  # ppm
        distances = [abs(shift - click_x) for shift in self.current_shifts]
        min_distance = min(distances)
        
        if min_distance > tolerance:
            return  # Click too far from any peak
        
        # Find the clicked peak index
        peak_idx = distances.index(min_distance)
        
        # Check for Shift or Ctrl key (using Qt modifiers)
        modifiers = QApplication.keyboardModifiers()
        is_multi = bool(modifiers & (Qt.KeyboardModifier.ShiftModifier | Qt.KeyboardModifier.ControlModifier))
        
        if peak_idx < len(self.current_shifts):
            if not is_multi:
                # Normal click
                if len(self.selected_peak_indices) == 1 and peak_idx in self.selected_peak_indices:
                    # User clicked the ONLY selected peak -> Toggle OFF (Deselect)
                    self.selected_peak_indices = set()
                else:
                    # New peak, or switching from multi-selection -> Select ONLY this peak
                    self.selected_peak_indices = {peak_idx}
            else:
                # Shift or Ctrl + Click: Toggle selection
                if peak_idx in self.selected_peak_indices:
                    self.selected_peak_indices.remove(peak_idx)
                else:
                    self.selected_peak_indices.add(peak_idx)
            
            # Force a fresh draw on next poll if desired
            self._last_highlight_atoms = set()
            
            # Update highlights
            self.highlight_selected_peaks()
            
            # Update 3D labels for all selected peaks
            self.update_selected_labels()
    
    def update_selected_labels(self, is_external_sync=False):
        """Update 3D labels and spheres for all selected peaks"""
        # 1. Clear existing labels
        self.clear_atom_labels()
        
        # 2. Add labels for each selected peak
        all_peak_indices = set()
        
        for peak_idx in sorted(self.selected_peak_indices):
            if peak_idx < len(self.peaks_metadata):
                # Get metadata for this peak (shift, intensity, is_merged, atom_indices)
                _, _, is_merged, atom_indices = self.peaks_metadata[peak_idx]
                
                # Add label for each atom in this peak (handles both merged and individual)
                for atom_idx in atom_indices:
                    all_peak_indices.add(atom_idx)
                    # Find atom symbol from data
                    atom_item = next((d for d in self.data if d.get('atom_idx') == atom_idx), None)
                    if atom_item:
                        atom_sym = atom_item.get('atom_sym', '?')
                        self.add_atom_label(atom_idx, atom_sym)
        
        # 3. Synchronize with Main Window
        if hasattr(self.parent_dlg, 'mw'):
            mw = self.parent_dlg.mw
            
            # If this is a sync FROM 3D (user clicked in viewer), we should NOT clear the 3D selection!
            # We only clear it if the user clicked in the Graph (internal sync), to replace Green with Yellow.
            if not is_external_sync:
                # Ensure we don't have double spheres (Green + Yellow).
                # We clear the global selection so ONLY our internal Yellow spheres are visible.
                if hasattr(mw, 'selected_atoms_3d'):
                    mw.selected_atoms_3d.clear()
                
                # CRITICAL: Update our sync tracker so the polling loop knows WE did this
                # and doesn't interpret the empty set as a user unselection.
                self._last_synced_mw_selection = frozenset()

                # Sync to MW if we are the originator (internal sync)
                if hasattr(mw, 'update_3d_selection_display'):
                    mw.update_3d_selection_display()
                elif hasattr(mw, 'update_selection_visuals'):
                    mw.update_selection_visuals()

            # Draw yellow highlights for NMR selection
            self.draw_custom_nmr_highlights_3d(all_peak_indices)

            # Debug print to confirming highlighting path is taken
            # print(f"Highlighting {len(self.selected_peak_indices)} peaks with {len(atom_coords)} atoms")
            
            # Render once after all labels added
            if hasattr(self.parent_dlg.mw, 'plotter'):
                self.parent_dlg.mw.plotter.render()
    
    def add_atom_label(self, atom_idx, atom_sym):
        """Add a single atom label to 3D viewer"""
        # Check if parent has plotter
        if not hasattr(self.parent_dlg, 'mw') or not hasattr(self.parent_dlg.mw, 'plotter'):
            return
        
        # Get coordinates
        coords = self.parent_dlg.parser.data.get('coords', [])
        if not coords or atom_idx >= len(coords):
            return
        
        try:
            pos = coords[atom_idx]
            label_pos = [pos[0], pos[1], pos[2] + 0.4]  # Offset above atom
            
            label_name = f"nmr_label_{atom_idx}"
            actor = self.parent_dlg.mw.plotter.add_point_labels(
                [label_pos],
                [f"{atom_sym}{atom_idx}"],
                font_size=12,
                text_color='cyan',
                point_size=0,
                always_visible=True,
                bold=True,
                name=label_name
            )
            self._atom_labels.append(actor)
            self._nmr_label_names.append(label_name)
        except Exception as e:
            print(f"Error adding atom label: {e}")
    
    def highlight_atom_in_3d(self, atom_idx, atom_sym):
        """Highlight selected atom with a label in 3D viewer (legacy - now uses update_selected_labels)"""
        # This is now handled by update_selected_labels
        pass
    
    def clear_peak_selection(self):
        """Clear all selected peaks and their labels"""
        # Clear selected peaks
        self.selected_peak_indices.clear()
        
        # Clear highlights
        for artist in self.highlight_artists:
            try:
                artist.remove()
            except:
                pass
        self.highlight_artists = []
        
        # Clear 3D labels
        self.clear_atom_labels()
        
        # [Commented out to avoid doubled spheres]
        if hasattr(self.parent_dlg, 'mw'):
            mw = self.parent_dlg.mw
            # We clear mw.selected_atoms_3d to avoid greenish spheres from "polluting" the view
            # The NMR dialog handles its own yellow highlights (internal logic)
            if hasattr(mw, 'selected_atoms_3d'):
                mw.selected_atoms_3d.clear()
            
            # Re-enable MW visual update so its internal highlights are also cleared
            if hasattr(mw, 'update_3d_selection_display'):
                mw.update_3d_selection_display()
            elif hasattr(mw, 'update_selection_visuals'):
                mw.update_selection_visuals()
        
        # Redraw spectrum
        if hasattr(self, 'canvas'):
            self.canvas.draw()
    
    def toggle_all_labels(self):
        """Toggle showing all atom labels on the spectrum graph"""
        show_all = self.chk_show_all_labels.isChecked()
        
        if show_all:
            # Enable show all mode (labels without red highlights)
            self.show_all_mode = True
            self.selected_peak_indices.clear()
            
            # Select all peaks to show labels
            # Use peaks_metadata length if available, otherwise displayed_data as fallback
            count = len(self.peaks_metadata) if hasattr(self, 'peaks_metadata') else len(self.displayed_data)
            for i in range(count):
                self.selected_peak_indices.add(i)
            
            # Update graph with labels only (no red highlights)
            self.highlight_selected_peaks()
            
            # For "Show All", we also want to update the 3D view to reflect "All"
            # or at least clear the specific "red" selection we had.
            # If we want to show labels for ALL atoms in 3D:
            self.update_selected_labels()
        else:
            # Disable show all mode
            self.show_all_mode = False
            # Clear all selections (graph and 3D)
            self.clear_peak_selection()
    
    def clear_atom_labels(self):
        """Remove all atom labels and custom selection spheres from 3D viewer"""
        if not hasattr(self.parent_dlg, 'mw') or not hasattr(self.parent_dlg.mw, 'plotter'):
            return
            
        plotter = self.parent_dlg.mw.plotter
        
        # 1. Clear custom NMR selection spheres by name (most reliable in PyVista)
        try:
            plotter.remove_actor('nmr_selection_highlights')
        except:
            pass
            
        # 2. Clear labels by tracked name
        if hasattr(self, '_nmr_label_names'):
            for name in self._nmr_label_names:
                try:
                    plotter.remove_actor(name)
                except:
                    pass
            self._nmr_label_names = []
            
        # 3. Fallback: Clear labels by list reference
        for actor in self._atom_labels:
            try:
                plotter.remove_actor(actor)
            except:
                pass
        self._atom_labels = []
        
        # 4. Clean up private spheres actor list
        if hasattr(self, '_nmr_sphere_actors'):
            for actor in self._nmr_sphere_actors:
                try:
                    plotter.remove_actor(actor)
                except:
                    pass
            self._nmr_sphere_actors = []
        
        try:
            plotter.render()
        except:
            pass
            
    def draw_custom_nmr_highlights_3d(self, atom_indices):
        """Draw yellow highlight spheres for selected atoms in 3D viewer"""
        if not hasattr(self.parent_dlg, 'mw') or not hasattr(self.parent_dlg.mw, 'plotter'):
            return
            
        plotter = self.parent_dlg.mw.plotter
        mw = self.parent_dlg.mw
        
        # ALWAYS clear existing custom highlights first to prevent stacking/phantom spheres
        try:
            plotter.remove_actor('nmr_selection_highlights')
        except:
            pass
            
        # Clear tracker list to prevent stale references
        self._nmr_sphere_actors = []
        
        # If no indices provided, just render the cleared state and return
        if not atom_indices or not hasattr(mw, 'atom_positions_3d'):
             try:
                 plotter.render()
             except:
                 pass
             return

        indices = list(atom_indices)
        valid_indices = [i for i in indices if i < len(mw.atom_positions_3d)]
        
        if not valid_indices:
            return
            
        try:
            # Get positions
            selected_positions = mw.atom_positions_3d[valid_indices]
            
            # Highlight sphere size: 40% (1.4x) relative to VDW radii per user request
            radii = []
            for i in valid_indices:
                # Find atom symbol from parser data
                atom_item = next((d for d in self.data if i == d.get('atom_idx')), None)
                sym = atom_item.get('atom_sym', 'C') if atom_item else 'C'
                # Use 1.4x scaling factor (40% larger)
                r = VDW_RADII.get(sym, 0.4) * 1.4
                radii.append(r)

            # Create glyphs for highlights
            highlight_source = pv.PolyData(selected_positions)
            highlight_source['radii'] = np.array(radii)
            
            highlight_glyphs = highlight_source.glyph(
                scale='radii',
                geom=pv.Sphere(radius=1.0, theta_resolution=16, phi_resolution=16),
                orient=False
            )
            
            # Add to plotter and track actor
            actor = plotter.add_mesh(
                highlight_glyphs,
                color='yellow',
                opacity=0.3,
                name='nmr_selection_highlights' 
            )
            self._nmr_sphere_actors.append(actor)
            plotter.render()
            
        except Exception as e:
            print(f"Error drawing custom NMR highlights: {e}")
    
    def closeEvent(self, event):
        """Clean up labels when dialog closes"""
        self.clear_atom_labels()
        super().closeEvent(event)
