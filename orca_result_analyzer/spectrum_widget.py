from PyQt6.QtWidgets import QWidget, QVBoxLayout
from PyQt6.QtCore import pyqtSignal
import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        # We manually call tight_layout() in plot_spectrum() to avoid clashing with constrained_layout.
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super().__init__(fig)

class SpectrumWidget(QWidget):
    clicked = pyqtSignal(object)
    range_changed = pyqtSignal(float, float, float, float, bool) # xmin, xmax, ymin, ymax, is_manual

    def __init__(self, data_list, x_key='energy', y_key='intensity', x_unit='nm', y_unit='Intensity', sigma=20.0, invert_x=False, invert_y=False):
        super().__init__()
        self.data_list = data_list
        self.x_key = x_key
        self.y_key = y_key
        self.x_unit = x_unit
        self.y_unit = y_unit
        self.y_unit_sticks = "Intensity" # Default for dual axis stick labels
        self.sigma = sigma
        self.invert_x = invert_x
        self.invert_y = invert_y
        self.normalization_mode = 'height' # 'height' or 'area'
        self.broaden_in_energy = False    # If True, broaden in cm-1
        
        self.x_range = None
        self.y_range = None
        
        if 'freq' in x_unit.lower() or 'cm' in x_unit.lower():
            self.invert_x = True
            self.x_range = (400, 4000)
        
        self.show_sticks = True
        self.show_gaussian = True
        self.show_legend = True
        
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        
        self.canvas = MplCanvas(self, width=7, height=3, dpi=100)
        layout.addWidget(self.canvas)
        
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        layout.addWidget(self.toolbar)
        
        self.selected_item = None
        self._initial_plot_done = False
        self._is_plotting = False
        self.plot_spectrum()
        
        self.canvas.mpl_connect('button_press_event', self.on_click)
        # Connect to axis changes (zoom/pan)
        self.canvas.axes.callbacks.connect('xlim_changed', self._on_axes_changed)
        self.canvas.axes.callbacks.connect('ylim_changed', self._on_axes_changed)

    def _on_axes_changed(self, ax):
        if self._is_plotting or not self._initial_plot_done:
            return
        
        xlim = self.canvas.axes.get_xlim()
        ylim = self.canvas.axes.get_ylim()
        
        xmin, xmax = min(xlim), max(xlim)
        ymin, ymax = min(ylim), max(ylim)
        
        # Update internal range so redraws (from selection etc) preserve this zoom
        # We only update if it's not Auto? Actually better to update always.
        self.x_range = (xmin, xmax)
        self.y_range = (ymin, ymax)
        
        self.range_changed.emit(xmin, xmax, ymin, ymax, True)

    def set_selected_item(self, item):
        self.selected_item = item
        self.plot_spectrum()

    def set_data(self, data_list):
        self.data_list = data_list
        self.plot_spectrum()
        
    def set_x_range(self, x_min, x_max):
        self.x_range = (x_min, x_max)
        self.plot_spectrum()
        
    def set_auto_x_range(self):
        self.x_range = None
        self.plot_spectrum()
        
    def set_y_range(self, y_min, y_max):
        self.y_range = (y_min, y_max)
        self.plot_spectrum()
        
    def set_auto_range(self):
        self.y_range = None
        self.plot_spectrum()

    def set_sigma(self, val):
        self.sigma = val
        self.plot_spectrum()
        
    def set_sticks(self, state):
        from PyQt6.QtCore import Qt
        self.show_sticks = (state == Qt.CheckState.Checked.value or state == True)
        self.plot_spectrum()

    def set_gaussian(self, state):
        from PyQt6.QtCore import Qt
        self.show_gaussian = (state == Qt.CheckState.Checked.value or state == True)
        self.plot_spectrum()
        
    def save_png(self, path):
        self.canvas.figure.savefig(path, dpi=300, bbox_inches='tight')
        
    def save_csv(self, path):
        import csv
        
        # Filter valid data
        points = []
        for item in self.data_list:
            x = item.get(self.x_key, 0.0)
            y = item.get(self.y_key, 0.0)
            if abs(y) > 1e-12: 
                points.append((x, y))
            
        if not points: return False
        
        xs = [p[0] for p in points]
        min_x = min(xs) - self.sigma * 3
        max_x = max(xs) + self.sigma * 3
        if max_x - min_x < 1.0:
            min_x -= 10
            max_x += 10
            
        # 1000 points resolution for CSV
        display_x = np.linspace(min_x, max_x, 1000)
        curve_y = np.zeros_like(display_x)
        
        for x0, y0 in points:
            term = np.exp(-0.5 * ((display_x - x0) / self.sigma)**2)
            curve_y += y0 * term
            
        try:
            with open(path, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerow(["X_Value", "Broadened_Intensity"])
                for i in range(len(display_x)):
                    writer.writerow([display_x[i], curve_y[i]])
            return True
        except Exception as e:
            return False

    def save_sticks_csv(self, path):
        import csv
        
        # Filter valid data
        points = []
        for item in self.data_list:
            x = item.get(self.x_key, 0.0)
            y = item.get(self.y_key, 0.0)
            if abs(y) > 1e-12: 
                points.append((x, y))
            
        if not points: return False
        
        # Sort by X
        points.sort(key=lambda p: p[0])
            
        try:
            with open(path, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                # Header
                x_head = self.x_unit if self.x_unit else "X"
                y_head = self.y_unit if self.y_unit else "Y"
                writer.writerow([x_head, y_head])
                for x, y in points:
                    writer.writerow([x, y])
            return True
        except Exception as e:
            return False

    def set_scaling(self, factor):
        self.scaling_factor = factor
        self.plot_spectrum()
        
    def set_dual_axis(self, enable):
        self.use_dual_axis = enable
        # Clear twin axis if disabling
        if not enable and hasattr(self, 'ax2'):
            try: 
                self.ax2.remove()
                del self.ax2
            except: pass
        self.plot_spectrum()

    def plot_spectrum(self):
        self._is_plotting = True
        try:
            self.canvas.axes.clear()
            if hasattr(self, 'ax2'):
                 try: self.ax2.clear()
                 except: pass
            
            if not hasattr(self, 'scaling_factor'): self.scaling_factor = 1.0
            if not hasattr(self, 'use_dual_axis'): self.use_dual_axis = False
            
            # Filter valid data
            all_points = []
            for item in self.data_list:
                x = item.get(self.x_key, 0.0)
                y = item.get(self.y_key, 0.0)
                if abs(y) > 1e-12: 
                     all_points.append((x, y))
                     
            if not all_points:
                self.canvas.axes.text(0.5, 0.5, "No Data", 
                                    ha='center', va='center', 
                                    transform=self.canvas.axes.transAxes)
                self.canvas.draw()
                return
            
            # Adaptive filtering
            points = all_points
                
            xs = [p[0] for p in points]
            ys = [p[1] for p in points]
            
            # Determine X Range (for plotting calculation)
            if self.x_range:
                new_xmin, new_xmax = min(self.x_range), max(self.x_range)
            else:
                new_xmin = min(xs) - self.sigma * 3
                new_xmax = max(xs) + self.sigma * 3
                if new_xmax - new_xmin < 1.0:
                    new_xmin -= 10
                    new_xmax += 10

            # Determine Plotting Range for Gaussian (Full data range + margins)
            # This ensures the curve remains visible when zooming out or clicking HOME
            data_min_x = min(xs)
            data_max_x = max(xs)
            plot_min_x = data_min_x - max(self.sigma * 10, 500)
            plot_max_x = data_max_x + max(self.sigma * 10, 500)
            
            # Use wavenumber grid if broadening in energy
            is_energy = 'cm' in self.x_unit.lower() or 'freq' in self.x_unit.lower()
            
            # Decide interpolation grid
            if self.broaden_in_energy and not is_energy:
                # Broaden in cm-1, but display is in nm
                # Convert data points to cm-1
                all_pts_cm = []
                for x0, y0 in all_points:
                    if x0 > 1e-6:
                        all_pts_cm.append((1e7 / x0, y0))
                
                if not all_pts_cm: return
                
                # Create wavenumber grid
                cms = [p[0] for p in all_pts_cm]
                cm_min_x = min(cms) - max(self.sigma * 10, 500)
                cm_max_x = max(cms) + max(self.sigma * 10, 500)
                grid_cm = np.linspace(cm_min_x, cm_max_x, 5000)
                curve_y_cm = np.zeros_like(grid_cm)
                
                # Broaden in cm-1
                norm_factor = 1.0
                if self.normalization_mode == 'area':
                    norm_factor = 1.0 / (np.sqrt(np.pi) * self.sigma)
                
                for x0, y0 in all_pts_cm:
                    term = np.exp(-((grid_cm - x0) / self.sigma)**2)
                    curve_y_cm += y0 * term * norm_factor
                
                # Interpolate back to nm display grid
                display_x = np.linspace(plot_min_x, plot_max_x, 3000)
                # Filter grid_cm to avoid divide by zero if 0 in grid_cm? 
                # (handled by min_x range usually)
                valid_grid = grid_cm[grid_cm > 1.0]
                valid_curve = curve_y_cm[grid_cm > 1.0]
                # x_nm = 1e7 / valid_grid
                # np.interp expects sorted x. x_nm will be sorted descending if valid_grid is ascending.
                curve_y = np.interp(display_x, 1e7 / valid_grid[::-1], valid_curve[::-1])
            else:
                # Standard broadening in current axis space
                display_x = np.linspace(plot_min_x, plot_max_x, 3000)
                curve_y = np.zeros_like(display_x)
                
                norm_factor = 1.0
                if self.normalization_mode == 'area':
                    norm_factor = 1.0 / (np.sqrt(np.pi) * self.sigma)
                
                if self.show_gaussian:
                    for x0, y0 in all_points:
                        term = np.exp(-((display_x - x0) / self.sigma)**2)
                        curve_y += y0 * term * norm_factor
            
            # APPLY SCALING
            curve_y *= self.scaling_factor
            
            # Plot Gaussian on Primary Axis
            if self.show_gaussian:
                self.canvas.axes.plot(display_x, curve_y, 'r-', linewidth=2, label='Spectrum')
                
            # Determine Axes for Sticks
            ax_sticks = self.canvas.axes
            if self.use_dual_axis:
                 if not hasattr(self, 'ax2'):
                     self.ax2 = self.canvas.axes.twinx()
                 ax_sticks = self.ax2
                 ax_sticks.set_ylabel(self.y_unit_sticks)
            elif hasattr(self, 'ax2'):
                 self.ax2.axis('off')
            
            # Plot Sticks
            if self.show_sticks:
                all_xs = [p[0] for p in all_points]
                all_ys = [p[1] for p in all_points]
                ax_sticks.vlines(all_xs, 0, all_ys, colors='black', alpha=0.6, linewidth=1.5, label='Transitions')
                
            # Highlight Selected Item
            if self.selected_item:
                sel_x = self.selected_item.get(self.x_key)
                sel_y = self.selected_item.get(self.y_key)
                if sel_x is not None and sel_y is not None:
                    ax_sticks.vlines([sel_x], 0, [sel_y], colors='blue', linewidth=2.5, zorder=5)
                    self.canvas.axes.axvline(sel_x, color='blue', linestyle='--', alpha=0.5, zorder=4)
            
            # Determine Y scale for Primary Axis (Curve)
            if self.y_range:
                 new_ymin, new_ymax = min(self.y_range), max(self.y_range)
            else:
                 g_max = np.max(curve_y) if self.show_gaussian and len(curve_y) > 0 else 1.0
                 g_min = np.min(curve_y) if self.show_gaussian and len(curve_y) > 0 else 0.0
                 if not self.use_dual_axis and self.show_sticks:
                      scale_max_sticks = max(ys) if ys else 1.0
                      g_max = max(g_max, scale_max_sticks)
                 if g_min >= 0:
                     new_ymin, new_ymax = 0, g_max * 1.1
                 else:
                     margin = (g_max - g_min) * 0.1
                     new_ymin, new_ymax = g_min - margin, g_max + margin

            # === FINAL AXIS AND CALLBACK CONFIGURATION ===
            # We do this at the very end to override any internal auto-scaling from plotting calls
            
            if self.invert_x:
                self.canvas.axes.set_xlim(new_xmax, new_xmin)
            else:
                self.canvas.axes.set_xlim(new_xmin, new_xmax)
                
            if self.invert_y:
                self.canvas.axes.set_ylim(new_ymax, new_ymin)
            else:
                self.canvas.axes.set_ylim(new_ymin, new_ymax)

            # Re-connect callbacks (ax.clear() removes them)
            self.canvas.axes.callbacks.connect('xlim_changed', self._on_axes_changed)
            self.canvas.axes.callbacks.connect('ylim_changed', self._on_axes_changed)

            # Set axis properties
            self.canvas.axes.set_xlabel(self.x_unit)
            self.canvas.axes.set_ylabel(self.y_unit)
            self.canvas.axes.grid(True, alpha=0.3)
            
            self._initial_plot_done = True
            
            # Use tight_layout to prevent clipping of labels/units
            # This is applied just before drawing
            try:
                self.canvas.figure.tight_layout()
            except: pass
                
            self.canvas.draw()
            
            # Emit range changed to sync UI
            xlim = self.canvas.axes.get_xlim()
            ylim = self.canvas.axes.get_ylim()
            self.range_changed.emit(min(xlim), max(xlim), min(ylim), max(ylim), False)
        finally:
            self._is_plotting = False
        
    def update(self):
        """Override to trigger plot update"""
        self.plot_spectrum()

    def on_click(self, event):
        if event.inaxes != self.canvas.axes: return
        
        # Double-click to reset zoom and selection
        if event.dblclick:
            self.selected_item = None
            self.clicked.emit(None)
            if 'freq' in self.x_unit.lower() or 'cm' in self.x_unit.lower():
                self.set_x_range(400, 4000)
            else:
                self.set_auto_x_range()
            self.set_auto_range()
            return
            
        if not self.data_list: return
        
        click_x = event.xdata
        if click_x is None: return

        # Relative tolerance: 1% of current view range
        xlim = self.canvas.axes.get_xlim()
        x_range = abs(xlim[1] - xlim[0])
        tolerance = x_range * 0.01

        # Find nearest point
        best_item = None
        min_dist = float('inf')

        for item in self.data_list:
            x = item.get(self.x_key, 0.0)
            y = item.get(self.y_key, 0.0)
            if abs(y) < 1e-12: continue

            dist = abs(x - click_x)
            if dist < min_dist:
                min_dist = dist
                best_item = item
        
        if best_item and min_dist <= tolerance:
            self.clicked.emit(best_item)
        else:
            # Clicked on empty space: clear selection
            self.selected_item = None
            self.clicked.emit(None)
            self.plot_spectrum()
