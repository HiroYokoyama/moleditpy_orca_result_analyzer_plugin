import os
import json
import math
import numpy as np
import pyvista as pv
from PyQt6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QDoubleSpinBox,
    QGroupBox,
    QSlider,
    QTableWidget,
    QTableWidgetItem,
    QCheckBox,
    QAbstractItemView,
    QMessageBox,
)
from PyQt6.QtCore import Qt, QTimer
import logging

try:
    from matplotlib.backends.backend_qtagg import (
        FigureCanvasQTAgg,
        NavigationToolbar2QT,
    )
    from matplotlib.figure import Figure

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


class ConvergenceGraphDialog(QDialog):
    def __init__(self, parent, traj_steps, current_idx=None):
        super().__init__(parent)
        self.setWindowTitle("Convergence Thresholds")
        self.resize(800, 1000)

        layout = QVBoxLayout(self)

        if not HAS_MATPLOTLIB:
            layout.addWidget(QLabel("Matplotlib is required to view the graph."))
            return

        self.figure = Figure()
        self.canvas = FigureCanvasQTAgg(self.figure)
        layout.addWidget(self.canvas)

        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        layout.addWidget(self.toolbar)

        self.plot_data(traj_steps, current_idx)

    def plot_data(self, traj_steps, current_idx):
        display_keys = {
            "rms gradient": "RMS Grad",
            "max gradient": "MAX Grad",
            "rms step": "RMS Step",
            "max step": "MAX Step",
            "energy change": "Energy Change",
        }

        # Gather data
        data = {k: [] for k in display_keys.keys()}
        targets = {k: None for k in display_keys.keys()}
        steps = []

        def safe_float(v):
            try:
                return float(v)
            except (ValueError, TypeError):
                return np.nan

        for i, step in enumerate(traj_steps):
            conv = step.get("convergence")
            if not conv:
                continue
            steps.append(i + 1)
            for k in display_keys.keys():
                val_obj = conv.get(k)
                if val_obj and isinstance(val_obj, dict):
                    data[k].append(safe_float(val_obj.get("value", np.nan)))
                    if targets[k] is None and "target" in val_obj:
                        targets[k] = safe_float(val_obj["target"])
                elif val_obj is not None:
                    data[k].append(safe_float(val_obj))
                else:
                    data[k].append(np.nan)

        keys_with_data = [k for k in display_keys.keys() if any(not np.isnan(v) for v in data[k])]

        if not steps or not keys_with_data:
            ax = self.figure.add_subplot(111)
            ax.text(0.5, 0.5, "No convergence data available", ha="center", va="center")
            return

        ax1 = self.figure.add_subplot(111)
        axes = [ax1]
        for _ in range(1, len(keys_with_data)):
            axes.append(ax1.twinx())

        # Adjust spines to prevent overlap when there are more than 2 axes
        if len(axes) > 2:
            axes[2].spines["right"].set_position(("axes", 1.15))
        if len(axes) > 3:
            axes[3].spines["left"].set_position(("axes", -0.15))
            axes[3].yaxis.set_label_position("left")
            axes[3].yaxis.set_ticks_position("left")
        if len(axes) > 4:
            axes[4].spines["right"].set_position(("axes", 1.30))

        # Add padding to figure to fit the extra axes
        self.figure.subplots_adjust(left=0.2, right=0.75, bottom=0.1, top=0.9)

        colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]
        lines = []
        labels = []

        for idx, (ax, k) in enumerate(zip(axes, keys_with_data)):
            name = display_keys[k]
            color = colors[idx % len(colors)]
            
            (line,) = ax.plot(steps, data[k], marker="o", color=color, label=name)
            lines.append(line)
            labels.append(name)
            
            # Plot the threshold target as a dashed line
            if targets[k] is not None:
                ax.axhline(
                    y=targets[k],
                    color=color,
                    linestyle="--",
                    linewidth=2,
                    alpha=0.7,
                )

            ax.set_ylabel(name, color=color)
            ax.tick_params(axis="y", colors=color)
            ax.set_yscale("symlog", linthresh=1e-6)

        ax1.set_xlabel("Optimization Step")
        ax1.grid(True, which="both", ls="-", alpha=0.2)
        
        # Create a single legend for all lines
        ax1.legend(lines, labels, loc="upper right")


class ForceViewerDialog(QDialog):
    def __init__(self, parent_dlg, gradients, parser=None):
        super().__init__(parent_dlg)
        self.setWindowTitle("Force Analysis")
        self.resize(700, 750)
        self.parent_dlg = parent_dlg
        self.gradients = gradients  # Current gradients from .out file
        self.parser = parser
        self.actors = []
        self.force_color = "red"
        self.settings_file = os.path.join(os.path.dirname(__file__), "settings.json")

        # Get trajectory steps if available
        self.traj_steps = []
        self.current_step_idx = -1  # -1 means current/final structure
        if self.parser and "scan_steps" in self.parser.data:
            self.traj_steps = self.parser.data["scan_steps"]

        main_layout = QVBoxLayout(self)

        # Visualization Controls (Restoring deleted code)
        view_group = QGroupBox("Visualization Settings")
        view_layout = QVBoxLayout(view_group)

        # Row 1: Visualize Button and Scale
        row1 = QHBoxLayout()

        self.btn_visualize = QPushButton("Visualize Forces")
        self.btn_visualize.setCheckable(True)  # Use checkable for toggle state
        self.btn_visualize.clicked.connect(self.toggle_visualization)
        row1.addWidget(self.btn_visualize)

        row1.addSpacing(20)
        row1.addWidget(QLabel("Vector Scale:"))
        self.spin_scale = QDoubleSpinBox()
        self.spin_scale.setRange(0.01, 100000.0)  # Increased range for small gradients
        self.spin_scale.setValue(1.0)
        self.spin_scale.setSingleStep(0.1)
        self.spin_scale.valueChanged.connect(self.update_vectors)
        row1.addWidget(self.spin_scale)

        # Auto Scale Button
        btn_autoscale = QPushButton("Auto Scale")
        btn_autoscale.clicked.connect(self.auto_scale)
        row1.addWidget(btn_autoscale)

        self.chk_auto_scale = QCheckBox("Auto Apply")
        self.chk_auto_scale.setToolTip(
            "Automatically auto-scale vectors when step changes"
        )
        row1.addWidget(self.chk_auto_scale)

        row1.addStretch()

        view_layout.addLayout(row1)
        main_layout.addWidget(view_group)

        # Force and Gradient Table
        table_label = QLabel("Force and Gradient Data (Eh/Bohr):")
        main_layout.addWidget(table_label)

        self.force_table = QTableWidget()
        self.force_table.setColumnCount(8)
        self.force_table.setHorizontalHeaderLabels(
            [
                "Atom",
                "Grad X",
                "Grad Y",
                "Grad Z",
                "Force X",
                "Force Y",
                "Force Z",
                "Force Mag",
            ]
        )
        self.force_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        header = self.force_table.horizontalHeader()
        header.setStretchLastSection(True)
        main_layout.addWidget(self.force_table)

        # Trajectory Navigation (if available) - Moved AFTER table creation
        # (becuase controls trigger update which updates table)
        if self.traj_steps:
            self._setup_trajectory_controls(main_layout)

        # Buttons row
        buttons_layout = QHBoxLayout()

        btn_reload = QPushButton("Reload from File")
        btn_reload.setToolTip("Reload data from the ORCA output file")
        btn_reload.clicked.connect(self.reload_data)
        buttons_layout.addWidget(btn_reload)

        buttons_layout.addStretch()

        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        buttons_layout.addWidget(btn_close)

        main_layout.addLayout(buttons_layout)

        # Initialize
        self.populate_force_table()

        self.load_settings()

        # Vector visualization is OFF by default as per user request
        self.btn_visualize.setChecked(False)

    def toggle_visualization(self):
        """Toggle force vector visualization"""
        if self.btn_visualize.isChecked():
            self.btn_visualize.setText("Clear Forces")
            self.update_vectors()
        else:
            self.btn_visualize.setText("Visualize Forces")
            self.clear_vectors()

    def auto_scale(self):
        """Calculate and set an appropriate scale factor"""
        if not self.gradients:
            return

        try:
            # Find max magnitude
            max_mag = 0.0
            for g in self.gradients:
                vec = g.get("grad", g.get("vector", None))
                if vec:
                    mag = np.linalg.norm(vec)
                    if mag > max_mag:
                        max_mag = mag

            if max_mag > 1e-12:
                # Target max vector length ~ 1.5 - 2.0 units (Angstrom)
                target = 2.0
                new_scale = target / max_mag

                # Round nicely
                if new_scale > 100:
                    new_scale = round(new_scale, -int(np.log10(new_scale)) + 1)
                else:
                    new_scale = round(new_scale, 2)

                self.spin_scale.blockSignals(True)
                self.spin_scale.setValue(new_scale)
                self.spin_scale.blockSignals(False)

                # If currently visualizing, update
                if self.btn_visualize.isChecked():
                    self.update_vectors()
        except Exception as _e:
            logging.warning("silenced: %s", _e)

    def _setup_trajectory_controls(self, layout):
        """Setup trajectory navigation controls"""
        traj_group = QGroupBox(
            f"Optimization Trajectory ({len(self.traj_steps)} steps)"
        )
        traj_layout = QVBoxLayout(traj_group)

        # Slider
        slider_layout = QHBoxLayout()
        slider_layout.addWidget(QLabel("Step:"))

        self.traj_slider = QSlider(Qt.Orientation.Horizontal)
        self.traj_slider.setRange(0, len(self.traj_steps))  # 0 to N-1 (traj), N = final

        self.traj_label = QLabel("Current (Final)")

        # Info labels (define early)
        self.traj_info = QLabel("Showing final optimized structure")
        self.traj_conv = QLabel("")
        self.traj_conv.setStyleSheet("color: #444; font-size: 9pt;")
        self.traj_conv.setWordWrap(True)

        # Navigation buttons
        btn_prev = QPushButton("<")
        btn_prev.setFixedWidth(30)
        btn_prev.clicked.connect(
            lambda: self.traj_slider.setValue(max(0, self.traj_slider.value() - 1))
        )

        btn_next = QPushButton(">")
        btn_next.setFixedWidth(30)
        btn_next.clicked.connect(
            lambda: self.traj_slider.setValue(
                min(self.traj_slider.maximum(), self.traj_slider.value() + 1)
            )
        )

        slider_layout.addWidget(btn_prev)
        slider_layout.addWidget(self.traj_slider)
        slider_layout.addWidget(btn_next)
        slider_layout.addWidget(self.traj_label)

        # Jump to Final button
        btn_final = QPushButton("Jump to Final")
        btn_final.setToolTip("Jump to final structure with force data")
        btn_final.clicked.connect(
            lambda: self.traj_slider.setValue(len(self.traj_steps))
        )
        slider_layout.addWidget(btn_final)

        traj_layout.addLayout(slider_layout)

        traj_layout.addWidget(self.traj_info)
        traj_layout.addWidget(self.traj_conv)

        # Show Convergence Graph Button
        btn_graph = QPushButton("Show Threshold Graph")
        btn_graph.setToolTip(
            "Show convergence thresholds plotted across all optimization steps"
        )
        btn_graph.clicked.connect(self.show_convergence_graph)
        traj_layout.addWidget(btn_graph)

        # Connect signals AFTER UI construction
        self.traj_slider.valueChanged.connect(self.on_trajectory_change)

        # Block signals while setting initial value to avoid firing before data is ready
        self.traj_slider.blockSignals(True)

        initial_val = len(self.traj_steps)
        if (
            self.parser
            and not self.parser.data.get("converged", False)
            and self.traj_steps
        ):
            # If running, search backwards for the last step that has convergence data AND coords
            found = False
            for idx in range(len(self.traj_steps) - 1, -1, -1):
                step = self.traj_steps[idx]
                if step.get("convergence") and step.get("coords") and step.get("atoms"):
                    initial_val = idx
                    found = True
                    break

            # If none found (very early in calculation), just default to the latest step
            if not found:
                initial_val = len(self.traj_steps)

        self.traj_slider.setValue(initial_val)
        self.traj_slider.blockSignals(False)

        # Manually trigger initial update
        self.on_trajectory_change(self.traj_slider.value())

        # Insert at the top (index 0) so it appears above other controls
        # even though it is initialized last.
        layout.insertWidget(0, traj_group)

    def show_convergence_graph(self):
        if getattr(self, "traj_steps", None) is None or not self.traj_steps:
            QMessageBox.warning(
                self, "No Data", "No trajectory convergence data available."
            )
            return
        
        # pass current_step_idx so it can draw a vertical line for the current frame
        current_idx = getattr(self, "current_step_idx", None)
        dlg = ConvergenceGraphDialog(self, self.traj_steps, current_idx)
        dlg.exec()

    def on_trajectory_change(self, val):
        """Handle trajectory step change"""
        self.current_step_idx = val
        num_steps = len(self.traj_steps)

        if val == num_steps:
            # Final/current structure
            self.traj_label.setText("Current (Final)")
            self.traj_info.setText("Showing final optimized structure")

            # Show convergence from the LAST trajectory step if it represents the same state
            # (which it usually does in a completed optimization)
            last_conv = {}
            if self.traj_steps:
                last_conv = self.traj_steps[-1].get("convergence", {})

            if last_conv:
                self._update_conv_label(last_conv)
            else:
                self.traj_conv.setText("")

            # Use current gradients
            self.gradients = (
                self.parser.data.get("gradients", []) if self.parser else []
            )

            # Auto-Scale if requested
            if (
                getattr(self, "chk_auto_scale", None) is not None
                and self.chk_auto_scale.isChecked()
            ):
                self.auto_scale()

            # Populate table
            self.populate_force_table()

            # Update structure first
            if (
                self.parser
                and "atoms" in self.parser.data
                and "coords" in self.parser.data
            ):
                self.update_structure(
                    self.parser.data["atoms"], self.parser.data["coords"]
                )

            # Then update vectors with a slight delay to ensure structure is drawn
            QTimer.singleShot(100, self.update_vectors)
        else:
            # Historical step
            step = self.traj_steps[val]
            energy = step.get("energy", 0.0)
            self.traj_label.setText(f"Step {val + 1}/{num_steps}")
            self.traj_info.setText(f"Energy: {energy:.8f} Eh")

            # Update convergence info
            self._update_conv_label(step.get("convergence", {}))

            # Update gradients for this step
            self.gradients = step.get("gradients", [])

            # Auto-Scale if requested
            if (
                getattr(self, "chk_auto_scale", None) is not None
                and self.chk_auto_scale.isChecked()
            ):
                self.auto_scale()

            # Populate table
            self.populate_force_table()

            # Update structure first
            atoms = step.get("atoms", [])
            coords = step.get("coords", [])
            if atoms and coords:
                self.update_structure(atoms, coords)

            # Then update vectors with delay
            QTimer.singleShot(100, self.update_vectors)

    def _update_conv_label(self, conv):
        """Helper to update the convergence info label with rich text and normalization"""
        if not conv:
            self.traj_conv.setText("No convergence data available")
            return

        # Map specific keys to display names (keys are now lowercased from parser)
        display_keys = {
            "rms gradient": "RMS Grad",
            "max gradient": "MAX Grad",
            "rms step": "RMS Step",
            "max step": "MAX Step",
            "energy change": "Energy Change",
        }

        items = []
        for k, v in conv.items():
            if isinstance(v, dict):
                dn = display_keys.get(k.lower(), k.title())
                status = v.get("converged", "??").upper()
                # Status determines color of the whole line
                if status == "YES":
                    color = "#28a745"  # Green
                    status_text = f" ({status})"
                elif status == "NO":
                    color = "#dc3545"  # Red
                    status_text = f" ({status})"
                elif status == "INFO":
                    color = "#000000"  # Black
                    status_text = ""
                else:
                    color = "#000000"
                    status_text = f" ({status})"

                items.append(
                    f"<span style='color: {color};'><b>{dn}:</b> {v.get('value')}{status_text}</span>"
                )
            else:
                items.append(f"<b>{k}:</b> {v:.6f}")

        num_items = len(items)
        if num_items == 0:
            self.traj_conv.setText("")
            return

        # Split into 2 columns
        half = math.ceil(num_items / 2)
        col1_items = items[:half]
        col2_items = items[half:]

        html = "<table width='100%'><tr>"
        html += "<td valign='top'>" + "<br>".join(col1_items) + "</td>"
        html += "<td valign='top'>" + "<br>".join(col2_items) + "</td>"
        html += "</tr></table>"

        self.traj_conv.setText(html)

    def reload_data(self):
        """Reload data from the output file"""
        if getattr(self, "_reloading", False):
            return
        self._reloading = True
        if not self.parser or not self.parser.filename:
            QMessageBox.warning(self, "Error", "No file associated with this parser.")
            self._reloading = False
            return

        if not os.path.exists(self.parser.filename):
            QMessageBox.warning(
                self, "Error", f"File not found: {self.parser.filename}"
            )
            self._reloading = False
            return

        try:
            content = ""
            encodings = ["utf-8", "utf-16", "latin-1", "cp1252"]
            found = False
            for enc in encodings:
                try:
                    with open(self.parser.filename, "r", encoding=enc) as f:
                        content = f.read()
                    found = True
                    break
                except UnicodeError:
                    continue

            if not found:
                with open(
                    self.parser.filename, "r", encoding="utf-8", errors="replace"
                ) as f:
                    content = f.read()

            # Re-parse (using the already initialized parser)
            self.parser.load_from_memory(content, self.parser.filename)

            # Update local state
            self.gradients = self.parser.data.get("gradients", [])
            self.traj_steps = self.parser.data.get("scan_steps", [])

            # Update slider if it exists
            if getattr(self, "traj_slider", None) is not None:
                self.traj_slider.setRange(-1, len(self.traj_steps) - 1)
                # Ensure current_step_idx is still valid
                if self.current_step_idx >= len(self.traj_steps):
                    self.current_step_idx = -1

                # Update label/text
                self.on_trajectory_change(self.current_step_idx)

            # Update force table and vectors
            self.populate_force_table()
            # self.auto_scale() # Use button only
            if self.btn_visualize.isChecked():
                self.update_vectors()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to reload data: {e}")
        finally:
            self._reloading = False

    def update_structure(self, atoms, coords):
        """Update the 3D structure visualization using RDKit"""
        try:
            from rdkit import Chem
            from rdkit.Geometry import Point3D
            from rdkit.Chem import rdDetermineBonds
            from .utils import normalize_atom_symbol, determine_bonds_without_dummies
        except ImportError:
            return

        if not hasattr(self, "_pt"):
            self._pt = Chem.GetPeriodicTable()

        mol = Chem.RWMol()
        conf = Chem.Conformer()
        pt = self._pt

        try:
            for i, sym in enumerate(atoms):
                sym = normalize_atom_symbol(sym)
                an = pt.GetAtomicNumber(sym)
                mol.AddAtom(Chem.Atom(an))
                conf.SetAtomPosition(
                    i, Point3D(coords[i][0], coords[i][1], coords[i][2])
                )
        except Exception:
            return

        mol.AddConformer(conf)

        # Determine bonds and bond orders on every load (no animation in this view).
        if rdDetermineBonds:
            try:
                charge = self.parser.data.get("charge", 0) if self.parser else 0
                determine_bonds_without_dummies(mol, charge=charge, bond_orders=True)
            except Exception as _e:
                logging.warning("silenced: %s", _e)

        final_mol = mol.GetMol()

        # Update in main window
        if self.parent_dlg.context:
            self.parent_dlg.context.draw_molecule_3d(final_mol)

    def populate_force_table(self):
        """Populate the force and gradient table from current gradient data"""
        self.force_table.setRowCount(0)

        # Get atoms from parser
        atoms = []
        if self.parser and "atoms" in self.parser.data:
            atoms = self.parser.data["atoms"]

        if not atoms or not self.gradients:
            # Show "No Data" message if gradients are missing but we are supposed to show something
            self.force_table.setRowCount(1)
            msg_item = QTableWidgetItem(
                "No gradient data available for this step in output file"
            )
            msg_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            self.force_table.setItem(0, 0, msg_item)
            self.force_table.setSpan(0, 0, 1, 8)
            return

        # Populate with gradient data (force = -gradient)
        for grad_item in self.gradients:
            atom_idx = grad_item.get("atom_idx", None)
            if atom_idx is None or atom_idx >= len(atoms):
                continue

            vec = grad_item.get("grad", grad_item.get("vector", None))
            if not vec or len(vec) != 3:
                continue

            # Gradient components
            gx, gy, gz = vec[0], vec[1], vec[2]

            # Force is negative gradient
            fx, fy, fz = -gx, -gy, -gz
            magnitude = np.sqrt(fx**2 + fy**2 + fz**2)

            atom_label = f"{atoms[atom_idx]}{atom_idx + 1}"  # 1-based

            row = self.force_table.rowCount()
            self.force_table.insertRow(row)
            self.force_table.setItem(row, 0, QTableWidgetItem(atom_label))
            self.force_table.setItem(row, 1, QTableWidgetItem(f"{gx:.8f}"))
            self.force_table.setItem(row, 2, QTableWidgetItem(f"{gy:.8f}"))
            self.force_table.setItem(row, 3, QTableWidgetItem(f"{gz:.8f}"))
            self.force_table.setItem(row, 4, QTableWidgetItem(f"{fx:.8f}"))
            self.force_table.setItem(row, 5, QTableWidgetItem(f"{fy:.8f}"))
            self.force_table.setItem(row, 6, QTableWidgetItem(f"{fz:.8f}"))
            self.force_table.setItem(row, 7, QTableWidgetItem(f"{magnitude:.8f}"))

        self.force_table.resizeColumnsToContents()

    def update_vectors(self):
        """Update the force vectors in the 3D visualizer"""
        try:
            mw = None
            if hasattr(self.parent_dlg, "context") and self.parent_dlg.context:
                mw = self.parent_dlg.context.get_main_window()
            elif hasattr(self.parent_dlg, "mw"):
                mw = self.parent_dlg.mw

            if not mw or not hasattr(mw, "plotter"):
                return

            # Clear old vectors
            self.clear_vectors()

            # Only show vectors if checked (visualizing)
            if not self.btn_visualize.isChecked():
                return

            # Get coordinates for the current view
            current_coords = []
            if self.current_step_idx < 0 or self.current_step_idx == len(
                self.traj_steps
            ):
                # Final / current structure (sentinel -1 or slider at max)
                if self.parser and "coords" in self.parser.data:
                    current_coords = self.parser.data["coords"]
            else:
                # Trajectory step
                current_coords = self.traj_steps[self.current_step_idx].get(
                    "coords", []
                )

            # Get settings
            scale = self.spin_scale.value()

            if not current_coords or not self.gradients:
                return

            # Draw force vectors
            for grad_item in self.gradients:
                atom_idx = grad_item.get("atom_idx", None)
                if atom_idx is None or atom_idx >= len(current_coords):
                    continue

                vec = grad_item.get("grad", grad_item.get("vector", None))
                if not vec or len(vec) != 3:
                    continue

                # Force = -Gradient
                force = np.array([-vec[0], -vec[1], -vec[2]])

                # Reverse if requested
                # if hasattr(self, 'chk_reverse') and self.chk_reverse.isChecked():
                #     force = -force
                magnitude = np.linalg.norm(force)

                if magnitude < 1e-12:
                    continue

                # Scale the vector
                scaled_force = force * scale
                scaled_mag = np.linalg.norm(scaled_force)

                if scaled_mag < 1e-12:
                    continue

                direction = scaled_force / scaled_mag

                # Create arrow
                arrow = pv.Arrow(
                    start=current_coords[atom_idx],
                    direction=direction,
                    scale=scaled_mag,
                    shaft_resolution=20,
                    tip_resolution=20,
                )

                actor = mw.plotter.add_mesh(
                    arrow, color=self.force_color, name=f"force_{atom_idx}"
                )
                self.actors.append(actor)

            mw.plotter.render()

        except Exception as e:
            logging.warning("Error drawing force vectors: %s", e)

    def clear_vectors(self):
        """Clear all force vector actors"""
        mw = None
        if hasattr(self.parent_dlg, "context") and self.parent_dlg.context:
            mw = self.parent_dlg.context.get_main_window()
        elif hasattr(self.parent_dlg, "mw"):
            mw = self.parent_dlg.mw

        if not mw or not hasattr(mw, "plotter"):
            return

        for actor in self.actors:
            try:
                mw.plotter.remove_actor(actor)
            except Exception as _e:
                logging.warning("silenced: %s", _e)

        self.actors = []
        mw.plotter.render()

    def closeEvent(self, event):
        """Clean up when dialog closes"""
        self.clear_vectors()
        self.save_settings()
        super().closeEvent(event)

    def load_settings(self):
        if os.path.exists(self.settings_file):
            try:
                with open(self.settings_file, "r", encoding="utf-8") as f:
                    all_settings = json.load(f)

                settings = all_settings.get("force_settings", {})

                # Note: 'scale' is intentionally NOT loaded to allow auto-scaling based on specific molecule data

                if "force_color" in settings:
                    self.force_color = settings["force_color"]

            except Exception as e:
                logging.warning("Error loading force settings: %s", e)

    def save_settings(self):
        all_settings = {}
        if os.path.exists(self.settings_file):
            try:
                with open(self.settings_file, "r", encoding="utf-8") as f:
                    all_settings = json.load(f)
            except Exception as _e:
                logging.warning("silenced: %s", _e)

        force_settings = {
            # "scale": self.spin_scale.value(), # Do not save scale
            # "reverse_vector": self.chk_reverse.isChecked() if hasattr(self, 'chk_reverse') else True,
            "force_color": self.force_color
        }

        all_settings["force_settings"] = force_settings

        try:
            with open(self.settings_file, "w", encoding="utf-8") as f:
                json.dump(all_settings, f, indent=2)
        except Exception as e:
            logging.warning("Error saving force settings: %s", e)
