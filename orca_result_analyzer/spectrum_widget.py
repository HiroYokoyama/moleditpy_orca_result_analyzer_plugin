from PyQt6.QtWidgets import QWidget, QVBoxLayout
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super().__init__(fig)

class SpectrumWidget(QWidget):
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
        from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        layout.addWidget(self.toolbar)
        
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

    def plot_spectrum(self):
        self.canvas.axes.clear()
        
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
            points = [(x, y) for x, y in all_points if y_min_filter <= y <= y_max_filter]
            if not points:  # No data in Y range, use all
                points = all_points
                
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
        
        # Determine Y scale
        if self.y_range:
             min_y_scale, max_y_scale = self.y_range
        else:
             # Auto range based on filtered data
             max_peak_y = np.max(curve_y) if self.show_gaussian and len(curve_y) > 0 else max(ys)
             min_peak_y = np.min(curve_y) if self.show_gaussian and len(curve_y) > 0 else min(ys)
             
             max_stick_y = max(ys) if ys else 1.0
             min_stick_y = min(ys) if ys else 0.0
             
             g_max = max_peak_y if self.show_gaussian else max_stick_y
             g_min = min_peak_y if self.show_gaussian else min_stick_y
             
             if self.show_gaussian and self.show_sticks:
                 g_max = max(g_max, max_stick_y)
                 g_min = min(g_min, min_stick_y)
             
             # If purely positive (Absorption), min usually 0
             if g_min >= 0:
                 min_y_scale = 0.0
                 max_y_scale = g_max * 1.1
                 if max_y_scale == 0: max_y_scale = 1.0
             else:
                 # CD: bounds with margin
                 margin = (g_max - g_min) * 0.1
                 if margin == 0: margin = 1.0
                 min_y_scale = g_min - margin
                 max_y_scale = g_max + margin
        
        # Plot Gaussian curve (use all_points for complete curve)
        if self.show_gaussian:
            self.canvas.axes.plot(display_x, curve_y, 'r-', linewidth=2, label='Broadened')
            
        # Plot Sticks (use all_points, not filtered)
        if self.show_sticks:
            all_xs = [p[0] for p in all_points]
            all_ys = [p[1] for p in all_points]
            self.canvas.axes.vlines(all_xs, 0, all_ys, colors='black', alpha=0.6, linewidth=1.5, label='Transitions')
        
        # Set axis properties
        self.canvas.axes.set_xlabel(self.x_unit)
        self.canvas.axes.set_ylabel(self.y_unit)
        self.canvas.axes.set_ylim(min_y_scale, max_y_scale)
        
        if self.invert_x:
            self.canvas.axes.set_xlim(max_x, min_x)
        else:
            self.canvas.axes.set_xlim(min_x, max_x)
            
        if self.invert_y:
            self.canvas.axes.invert_yaxis()
        
        self.canvas.axes.grid(True, alpha=0.3)
        
        # Add legend if both shown and legend enabled
        if self.show_gaussian and self.show_sticks and self.show_legend:
            self.canvas.axes.legend(loc='best')
        
        self.canvas.draw()
        
    def update(self):
        """Override to trigger plot update"""
        self.plot_spectrum()
