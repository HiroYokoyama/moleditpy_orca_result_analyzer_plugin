from PyQt6.QtWidgets import QWidget, QVBoxLayout
from PyQt6.QtCore import pyqtSignal
import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super().__init__(fig)

class SpectrumWidget(QWidget):
    clicked = pyqtSignal(object)

    def __init__(self, data_list, x_key='energy', y_key='intensity', x_unit='nm', y_unit='Intensity', sigma=20.0, invert_x=False, invert_y=False):
        """
        data_list: List of dicts containing the data.
        x_key: Key for X-axis value (e.g., 'energy_nm', 'freq').
        y_key: Key for Y-axis value (e.g., 'osc', 'ir', 'raman', 'rotatory_strength').
        x_unit: Unit string for X-axis label.
        y_unit: Unit string for Y-axis label.
        sigma: Initial Gaussian broadening width.
        invert_x: If True, plot X axis from High to Low (e.g. for IR).
        invert_y: If True, plot Y axis Downwards (e.g. for Transmittance-style).
        """
        super().__init__()
        self.data_list = data_list
        self.x_key = x_key
        self.y_key = y_key
        self.x_unit = x_unit
        self.y_unit = y_unit
        self.sigma = sigma
        self.invert_x = invert_x
        self.invert_y = invert_y
        
        self.x_range = None # (min, max) or None for auto
        self.y_range = None # (min, max) or None for auto
        
        # Auto-detect IR convention
        if 'freq' in x_unit.lower() or 'cm' in x_unit.lower():
            self.invert_x = True # Standard IR X
        
        self.show_sticks = True
        self.show_gaussian = True
        self.show_legend = True
        
        # Setup matplotlib canvas
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        
        self.canvas = MplCanvas(self, width=7, height=4, dpi=100)
        layout.addWidget(self.canvas)
        
        # Enable interactive toolbar features
        from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        layout.addWidget(self.toolbar)
        
        self.plot_spectrum()
        
        self.canvas.mpl_connect('button_press_event', self.on_click)

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
        self.canvas.axes.clear()
        if hasattr(self, 'ax2'): # Cleanup previous twinx to avoid overlap/duplication
             try: self.ax2.clear()
             except: pass
             # We might need to completely remove and recreate if axes property is persistent?
             # Matplotlib clear() usually clears content but twinx persists.
             # We'll re-checkout ax2 inside logic.
        
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
        
        # Adaptive filtering: If one axis is manual, filter points for the other axis
        points = all_points
        
        # If Y is manual, filter X to only visible Y range
        if self.y_range and not self.x_range:
            y_min_filter, y_max_filter = self.y_range
            # Note: filtering based on UNSCALED Y for sticking? Or Scaled?
            # Usually users filter based on what they see. 
            # If scaling is active, y_range is likely scaled.
            # Let's apply scaling filter logic carefully.
            # Usually filtering is for data reduction.
            # Let's just use all_points for robust plotting unless huge.
            points = all_points # Simplified
                
        # If X is manual, filter Y to only visible X range  
        elif self.x_range and not self.y_range:
            x_min_filter, x_max_filter = self.x_range
            points = [(x, y) for x, y in all_points if x_min_filter <= x <= x_max_filter]
            if not points:  # No data in X range, use all
                points = all_points
            
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        
        # Determine X Range
        if self.x_range:
            min_x, max_x = self.x_range
        else:
            # Auto range based on filtered data
            min_x = min(xs) - self.sigma * 3
            max_x = max(xs) + self.sigma * 3
            
            if max_x - min_x < 1.0:
                min_x -= 10
                max_x += 10
            
        # Calculate Gaussian broadening (use all_points, not filtered)
        display_x = np.linspace(min_x, max_x, 1000)
        curve_y = np.zeros_like(display_x)
        
        if self.show_gaussian:
            for x0, y0 in all_points:  # Use all data for broadening
                term = np.exp(-0.5 * ((display_x - x0) / self.sigma)**2)
                curve_y += y0 * term
        
        # APPLY SCALING
        curve_y *= self.scaling_factor
        
        # Plot Gaussian on Primary Axis
        color_curve = 'r-'
        if self.show_gaussian:
            self.canvas.axes.plot(display_x, curve_y, color_curve, linewidth=2, label='Spectrum')
            
        # Determine Axes for Sticks
        ax_sticks = self.canvas.axes
        if self.use_dual_axis:
             if not hasattr(self, 'ax2'):
                 self.ax2 = self.canvas.axes.twinx()
             ax_sticks = self.ax2
             ax_sticks.set_ylabel("Oscillator Strength") # Fixed Label for now
        elif hasattr(self, 'ax2'):
             # If switching off dual axis, hide it
             self.ax2.axis('off')
        
        # Plot Sticks
        if self.show_sticks:
            all_xs = [p[0] for p in all_points]
            all_ys = [p[1] for p in all_points]
            # Sticks are NOT scaled by scaling_factor usually (they represent discrete f)
            ax_sticks.vlines(all_xs, 0, all_ys, colors='black', alpha=0.6, linewidth=1.5, label='Transitions')
        
        # Determine Y scale for Primary Axis (Curve)
        g_max = np.max(curve_y) if self.show_gaussian and len(curve_y) > 0 else 1.0
        g_min = np.min(curve_y) if self.show_gaussian and len(curve_y) > 0 else 0.0
        
        # If dual axis, Primary only cares about Curve
        if not self.use_dual_axis and self.show_sticks:
             # Mix of Curve and Sticks
             # Check scaled sticks? No sticks are unscaled f. 
             # If scaling factor is huge (Epsilon), mixture is broken. 
             # Hence Dual Axis is required for Epsilon.
             scale_max_sticks = max(ys) if ys else 1.0
             g_max = max(g_max, scale_max_sticks)
        
        if self.y_range:
             self.canvas.axes.set_ylim(self.y_range)
        else:
             # Auto scale Primary
             if g_min >= 0:
                 self.canvas.axes.set_ylim(0, g_max * 1.1)
             else:
                 margin = (g_max - g_min) * 0.1
                 self.canvas.axes.set_ylim(g_min - margin, g_max + margin)

        # Set axis properties
        self.canvas.axes.set_xlabel(self.x_unit)
        self.canvas.axes.set_ylabel(self.y_unit)
        
        if self.invert_x:
            self.canvas.axes.set_xlim(max_x, min_x)
        else:
            self.canvas.axes.set_xlim(min_x, max_x)
            
        if self.invert_y:
            self.canvas.axes.invert_yaxis()
        
        self.canvas.axes.grid(True, alpha=0.3)
        self.canvas.draw()
        
    def update(self):
        """Override to trigger plot update"""
        self.plot_spectrum()

    def on_click(self, event):
        if event.inaxes != self.canvas.axes: return
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
