"""
Tests for force_analysis.py including ConvergenceGraphDialog.
"""

import sys
import unittest
from unittest.mock import MagicMock

# Stub Qt
if "PyQt6" not in sys.modules:
    pyqt6 = MagicMock()
    qtw = MagicMock()

    class _QDialog:
        def __init__(self, *a, **k):
            pass

        def setWindowTitle(self, *a):
            pass

        def resize(self, *a):
            pass

        def close(self):
            pass

        def setLayout(self, *a):
            pass

    qtw.QDialog = _QDialog

    sys.modules["PyQt6"] = pyqt6
    sys.modules["PyQt6.QtWidgets"] = qtw
    sys.modules["PyQt6.QtCore"] = MagicMock()
    sys.modules["PyQt6.QtGui"] = MagicMock()

    # Stub matplotlib and pyvista
    sys.modules["matplotlib.backends.backend_qtagg"] = MagicMock()
    sys.modules["matplotlib.figure"] = MagicMock()
    sys.modules["pyvista"] = MagicMock()

# Import the module to test
from orca_result_analyzer.force_analysis import ConvergenceGraphDialog


class TestForceAnalysis(unittest.TestCase):
    def test_convergence_graph_dialog_plot_data(self):
        """Test the data parsing logic inside ConvergenceGraphDialog plot_data"""
        # Create a mock figure to capture ax.plot
        mock_figure = MagicMock()
        mock_ax = MagicMock()
        mock_figure.add_subplot.return_value = mock_ax

        mock_line = MagicMock()
        mock_ax.plot.return_value = [mock_line]

        # Override the Figure class in the module temporarily
        import orca_result_analyzer.force_analysis as fa

        original_figure = getattr(fa, "Figure", None)
        fa.Figure = MagicMock(return_value=mock_figure)

        try:
            # Traj steps containing convergence data
            traj_steps = [
                {
                    "convergence": {
                        "rms gradient": {
                            "value": 0.001,
                            "target": 0.0001,
                            "converged": "NO",
                        },
                        "energy change": {
                            "value": -0.5,
                            "target": 0.00001,
                            "converged": "NO",
                        },
                    }
                },
                {
                    "convergence": {
                        "rms gradient": {
                            "value": 0.00005,
                            "target": 0.0001,
                            "converged": "YES",
                        },
                        "energy change": {
                            "value": -0.000001,
                            "target": 0.00001,
                            "converged": "YES",
                        },
                    }
                },
            ]

            dlg = ConvergenceGraphDialog(None, traj_steps)

            # Verify ax.plot and ax.axhline were called
            self.assertTrue(mock_ax.plot.called)
            self.assertTrue(mock_ax.axhline.called)

            # The plot should be called for RMS Grad and Energy Change
            plot_calls = mock_ax.plot.call_args_list
            self.assertEqual(len(plot_calls), 2)

            # Check steps values (x-axis)
            x_vals = plot_calls[0][0][0]
            self.assertEqual(list(x_vals), [1, 2])

        finally:
            if original_figure:
                fa.Figure = original_figure


if __name__ == "__main__":
    unittest.main()
