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

    def test_convergence_graph_dialog_spines_and_colors(self):
        """Test that top spine is visible and twin axes are colored in ConvergenceGraphDialog"""
        import orca_result_analyzer.force_analysis as fa

        mock_figure = MagicMock()
        mock_ax1 = MagicMock()
        mock_ax2 = MagicMock()
        mock_figure.add_subplot.return_value = mock_ax1
        mock_ax1.twinx.return_value = mock_ax2

        # Create mocks for spines
        spines1 = {k: MagicMock() for k in ["left", "bottom", "right", "top"]}
        spines2 = {k: MagicMock() for k in ["left", "bottom", "right", "top"]}
        mock_ax1.spines = spines1
        mock_ax2.spines = spines2

        mock_ax1.plot.return_value = [MagicMock()]
        mock_ax2.plot.return_value = [MagicMock()]

        original_figure = getattr(fa, "Figure", None)
        fa.Figure = MagicMock(return_value=mock_figure)

        try:
            # Traj steps containing convergence data to trigger multi axis (is_multi = True)
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
                }
            ]

            ConvergenceGraphDialog(None, traj_steps)

            # Check that axes[0] (ax1) spines left, bottom, and top are set to visible
            spines1["left"].set_visible.assert_called_with(True)
            spines1["bottom"].set_visible.assert_called_with(True)
            spines1["top"].set_visible.assert_called_with(True)

            # Check that twin axis spine is colored to match the metric color
            # Since index of energy change is > 0, it should color mock_ax2's right spine
            spines2["right"].set_color.assert_called()
            # The first axis (ax1) spine should NOT be colored (keep black axis)
            spines1["left"].set_color.assert_not_called()

        finally:
            if original_figure:
                fa.Figure = original_figure

    def test_force_viewer_dialog_get_last_force_containing_step_idx(self):
        """Test get_last_force_containing_step_idx method of ForceViewerDialog"""
        from orca_result_analyzer.force_analysis import ForceViewerDialog

        # Mock parser and parent
        mock_parent = MagicMock()
        mock_parent.context = MagicMock()

        mock_parser = MagicMock()
        mock_parser.data = {
            "gradients": [{"atom_idx": 0, "vector": [0.1, 0.2, 0.3]}],
            "scan_steps": [
                {
                    "step": 1,
                    "gradients": [{"atom_idx": 0, "vector": [0.01, 0.02, 0.03]}],
                    "atoms": ["C"],
                    "coords": [[0.0, 0.0, 0.0]],
                },
                {
                    "step": 2,
                    "gradients": [],
                    "atoms": ["C"],
                    "coords": [[0.1, 0.1, 0.1]],
                },
            ],
            "termination_status": "Running",
        }

        # Gradients are present in final/current structure
        dlg = ForceViewerDialog(
            mock_parent, mock_parser.data["gradients"], parser=mock_parser
        )

        # Test index when final gradients are non-empty
        self.assertEqual(dlg.get_last_force_containing_step_idx(), 2)

        # Test index when final gradients are empty but first step has gradients
        mock_parser.data["gradients"] = []
        dlg.gradients = []
        self.assertEqual(dlg.get_last_force_containing_step_idx(), 0)

        # Test index when no gradients are found anywhere
        mock_parser.data["scan_steps"][0]["gradients"] = []
        self.assertEqual(dlg.get_last_force_containing_step_idx(), 2)

    def test_force_viewer_dialog_reload_data_running(self):
        """Test reload_data in ForceViewerDialog when job is running"""
        from orca_result_analyzer.force_analysis import ForceViewerDialog

        mock_parent = MagicMock()
        mock_parent.context = MagicMock()

        mock_parser = MagicMock()
        mock_parser.filename = "dummy.out"
        mock_parser.data = {
            "gradients": [],
            "scan_steps": [
                {
                    "step": 1,
                    "gradients": [{"atom_idx": 0, "vector": [0.01, 0.02, 0.03]}],
                    "atoms": ["C"],
                    "coords": [[0.0, 0.0, 0.0]],
                },
                {
                    "step": 2,
                    "gradients": [],
                    "atoms": ["C"],
                    "coords": [[0.1, 0.1, 0.1]],
                },
            ],
            "termination_status": "Running",
        }

        # Mock file reading inside reload_data
        import builtins

        original_open = builtins.open

        # Helper mock for open
        mock_open = MagicMock()
        mock_open.__enter__.return_value.read.return_value = "dummy content"

        # Stub os.path.exists
        import os

        original_exists = os.path.exists
        os.path.exists = MagicMock(return_value=True)

        try:
            builtins.open = MagicMock(return_value=mock_open)

            dlg = ForceViewerDialog(mock_parent, [], parser=mock_parser)

            # Mock UI slider
            mock_slider = MagicMock()
            dlg.traj_slider = mock_slider

            # Run reload_data
            dlg.reload_data()

            # Verify slider range is set from 0 to len(traj_steps) -> 0 to 2
            mock_slider.setRange.assert_called_with(0, 2)

            # Since the job is running and only step 0 has gradients, self.current_step_idx should be 0
            self.assertEqual(dlg.current_step_idx, 0)
            mock_slider.setValue.assert_called_with(0)

        finally:
            builtins.open = original_open
            os.path.exists = original_exists


if __name__ == "__main__":
    unittest.main()
