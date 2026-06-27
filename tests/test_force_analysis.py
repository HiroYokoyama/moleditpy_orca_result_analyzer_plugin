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
        mock_ax1 = MagicMock()
        mock_ax2 = MagicMock()
        mock_figure.add_subplot.return_value = mock_ax1
        mock_ax1.twinx.return_value = mock_ax2

        mock_line1 = MagicMock()
        mock_ax1.plot.return_value = [mock_line1]
        mock_line2 = MagicMock()
        mock_ax2.plot.return_value = [mock_line2]

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
                            "tolerance": 0.0001,
                            "converged": "NO",
                        },
                        "energy change": {
                            "value": -0.5,
                            "tolerance": 0.00001,
                            "converged": "NO",
                        },
                    }
                },
                {
                    "convergence": {
                        "rms gradient": {
                            "value": 0.00005,
                            "tolerance": 0.0001,
                            "converged": "YES",
                        },
                        "energy change": {
                            "value": -0.000001,
                            "tolerance": 0.00001,
                            "converged": "YES",
                        },
                    }
                },
            ]

            dlg = ConvergenceGraphDialog(None, traj_steps)

            # Verify ax.plot and ax.axhline were called
            self.assertTrue(mock_ax1.plot.called)
            self.assertTrue(mock_ax1.axhline.called)
            self.assertTrue(mock_ax2.plot.called)
            self.assertTrue(mock_ax2.axhline.called)

            # The plot should be called for RMS Grad and Energy Change
            plot_calls1 = mock_ax1.plot.call_args_list
            plot_calls2 = mock_ax2.plot.call_args_list
            self.assertEqual(len(plot_calls1), 1)
            self.assertEqual(len(plot_calls2), 1)

            # Check steps values (x-axis)
            x_vals = plot_calls1[0][0][0]
            self.assertEqual(list(x_vals), [1, 2])

        finally:
            if original_figure:
                fa.Figure = original_figure

    def _make_traj_steps(self):
        """Helper: two steps with tolerance-keyed convergence data."""
        return [
            {
                "convergence": {
                    "rms gradient": {
                        "value": 0.001,
                        "tolerance": 0.0001,
                        "converged": "NO",
                    },
                    "max gradient": {
                        "value": 0.002,
                        "tolerance": 0.0003,
                        "converged": "NO",
                    },
                    "energy change": {
                        "value": -0.5,
                        "tolerance": 0.00001,
                        "converged": "NO",
                    },
                }
            },
            {
                "convergence": {
                    "rms gradient": {
                        "value": 0.00005,
                        "tolerance": 0.0001,
                        "converged": "YES",
                    },
                    "max gradient": {
                        "value": 0.0001,
                        "tolerance": 0.0003,
                        "converged": "YES",
                    },
                    "energy change": {
                        "value": -0.000001,
                        "tolerance": 0.00001,
                        "converged": "YES",
                    },
                }
            },
        ]

    def test_tolerance_key_draws_axhline(self):
        """Threshold line must appear when parser stores 'tolerance' (not 'target')."""
        import orca_result_analyzer.force_analysis as fa

        mock_figure = MagicMock()
        mock_ax = MagicMock()
        mock_figure.add_subplot.return_value = mock_ax
        mock_ax.twinx.return_value = MagicMock()
        mock_ax.plot.return_value = [MagicMock()]
        mock_ax.twinx.return_value.plot.return_value = [MagicMock()]

        original_figure = getattr(fa, "Figure", None)
        fa.Figure = MagicMock(return_value=mock_figure)
        try:
            ConvergenceGraphDialog(None, self._make_traj_steps())
            # axhline must have been called (threshold line visible)
            self.assertTrue(
                mock_ax.axhline.called, "axhline not called – threshold line missing"
            )
            # Verify the y value matches the tolerance
            y_val = mock_ax.axhline.call_args_list[0][1]["y"]
            self.assertAlmostEqual(y_val, 0.0001)
        finally:
            if original_figure:
                fa.Figure = original_figure

    def test_target_key_compat_draws_axhline(self):
        """'target' key (legacy) must also produce a threshold line."""
        import orca_result_analyzer.force_analysis as fa

        mock_figure = MagicMock()
        mock_ax = MagicMock()
        mock_figure.add_subplot.return_value = mock_ax
        mock_ax.twinx.return_value = MagicMock()
        mock_ax.plot.return_value = [MagicMock()]
        mock_ax.twinx.return_value.plot.return_value = [MagicMock()]

        original_figure = getattr(fa, "Figure", None)
        fa.Figure = MagicMock(return_value=mock_figure)
        try:
            legacy_steps = [
                {
                    "convergence": {
                        "rms gradient": {
                            "value": 0.001,
                            "target": 0.0001,
                            "converged": "NO",
                        }
                    }
                }
            ]
            ConvergenceGraphDialog(None, legacy_steps)
            self.assertTrue(
                mock_ax.axhline.called, "axhline not called for 'target' key"
            )
        finally:
            if original_figure:
                fa.Figure = original_figure

    def test_empty_tolerance_no_axhline(self):
        """Empty tolerance string (Max(Bonds) info rows) must not draw an axhline."""
        import orca_result_analyzer.force_analysis as fa

        mock_figure = MagicMock()
        mock_ax = MagicMock()
        mock_figure.add_subplot.return_value = mock_ax
        mock_ax.twinx.return_value = MagicMock()
        mock_ax.plot.return_value = [MagicMock()]

        original_figure = getattr(fa, "Figure", None)
        fa.Figure = MagicMock(return_value=mock_figure)
        try:
            no_tol_steps = [
                {
                    "convergence": {
                        "rms gradient": {
                            "value": 0.001,
                            "tolerance": "",
                            "converged": "INFO",
                        }
                    }
                }
            ]
            ConvergenceGraphDialog(None, no_tol_steps)
            self.assertFalse(
                mock_ax.axhline.called,
                "axhline should not be called for empty tolerance",
            )
        finally:
            if original_figure:
                fa.Figure = original_figure

    def test_metric_filter_single(self):
        """Selecting a single metric must call plot only once on ax1."""
        import orca_result_analyzer.force_analysis as fa

        mock_figure = MagicMock()
        mock_ax = MagicMock()
        mock_figure.add_subplot.return_value = mock_ax
        mock_ax.plot.return_value = [MagicMock()]  # must return a list with one line

        original_figure = getattr(fa, "Figure", None)
        fa.Figure = MagicMock(return_value=mock_figure)
        try:
            dlg = ConvergenceGraphDialog(None, self._make_traj_steps())
            # Reset call counts, then call plot_data with a single metric
            mock_ax.reset_mock()
            dlg.plot_data(self._make_traj_steps(), None, selection="RMS Grad")
            self.assertEqual(mock_ax.plot.call_count, 1)
            # twinx should NOT be called for a single metric
            mock_ax.twinx.assert_not_called()
        finally:
            if original_figure:
                fa.Figure = original_figure

    def test_no_convergence_data_no_crash(self):
        """Steps with no convergence key must not crash."""
        import orca_result_analyzer.force_analysis as fa

        mock_figure = MagicMock()
        mock_ax = MagicMock()
        mock_figure.add_subplot.return_value = mock_ax

        original_figure = getattr(fa, "Figure", None)
        fa.Figure = MagicMock(return_value=mock_figure)
        try:
            empty_steps = [{"atoms": ["C"], "coords": [[0, 0, 0]]}]
            ConvergenceGraphDialog(None, empty_steps)  # must not raise
            # Should display the "No convergence data" text
            self.assertTrue(mock_ax.text.called)
        finally:
            if original_figure:
                fa.Figure = original_figure


if __name__ == "__main__":
    unittest.main()
