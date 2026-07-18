"""
tests/test_scf_analysis.py
Unit tests for SCFTraceDialog's pure logic: CSV export (all-blocks
concatenation vs per-block) and the update_plot guard branches.

PyQt6 and the matplotlib Qt backend are stubbed so the module imports
headlessly; the dialog's Qt constructor is bypassed (object.__new__) and
lightweight fakes stand in for the combo box, axes, canvas and parent.
"""

import os
import sys
import csv
import types
import tempfile
import importlib.util
import unittest
from unittest.mock import MagicMock


# ---------------------------------------------------------------------------
# Stub PyQt6 + the matplotlib Qt backend before importing scf_analysis
# ---------------------------------------------------------------------------


def _install_stubs():
    for name in ["PyQt6", "PyQt6.QtWidgets", "PyQt6.QtCore", "PyQt6.QtGui"]:
        sys.modules.setdefault(name, types.ModuleType(name))

    class _Base:
        def __init__(self, *a, **k):
            pass

        def __getattr__(self, n):
            return MagicMock()

    qtw = sys.modules["PyQt6.QtWidgets"]
    for cls in [
        "QDialog",
        "QVBoxLayout",
        "QHBoxLayout",
        "QPushButton",
        "QLabel",
        "QComboBox",
        "QGroupBox",
        "QFileDialog",
    ]:
        setattr(qtw, cls, type(cls, (_Base,), {}))

    # matplotlib's qtagg backend needs a real Qt binding; stub it out.
    sys.modules["matplotlib.backends.backend_qtagg"] = MagicMock()


_install_stubs()

_PKG_DIR = os.path.normpath(
    os.path.join(os.path.dirname(__file__), "..", "orca_result_analyzer")
)

# Load under a throwaway package so `from .utils import get_default_export_path`
# resolves to the real utils without executing the Qt-heavy package __init__.
_PKG = "_scf_test_pkg"
if _PKG not in sys.modules:
    _pkg = types.ModuleType(_PKG)
    _pkg.__path__ = [_PKG_DIR]
    sys.modules[_PKG] = _pkg
if f"{_PKG}.utils" not in sys.modules:
    _uspec = importlib.util.spec_from_file_location(
        f"{_PKG}.utils", os.path.join(_PKG_DIR, "utils.py")
    )
    _umod = importlib.util.module_from_spec(_uspec)
    _uspec.loader.exec_module(_umod)
    sys.modules[f"{_PKG}.utils"] = _umod

_spec = importlib.util.spec_from_file_location(
    f"{_PKG}.scf_analysis", os.path.join(_PKG_DIR, "scf_analysis.py")
)
S = importlib.util.module_from_spec(_spec)
S.__package__ = _PKG
_spec.loader.exec_module(S)


# ---------------------------------------------------------------------------
# Fakes
# ---------------------------------------------------------------------------


class _FakeCombo:
    def __init__(self, data):
        self._data = data

    def currentData(self):
        return self._data


def _traces():
    return [
        {
            "step": "Geom 1",
            "iterations": [
                {"iter": 1, "energy": -76.0},
                {"iter": 2, "energy": -76.1},
            ],
        },
        {
            "step": "Geom 2",
            "iterations": [
                {"iter": 1, "energy": -76.2},
            ],
        },
    ]


def _bare_dialog(traces, current):
    inst = S.SCFTraceDialog.__new__(S.SCFTraceDialog)
    inst.scf_traces = traces
    inst.combo_steps = _FakeCombo(current)
    inst.ax = MagicMock()
    inst.figure = MagicMock()
    inst.canvas = MagicMock()
    return inst


# ---------------------------------------------------------------------------
# export_csv
# ---------------------------------------------------------------------------


class TestExportCsv(unittest.TestCase):
    def _run_export(self, dlg, out):
        parent = types.SimpleNamespace(file_path="job.out", context=None)
        dlg.parent = MagicMock(return_value=parent)
        saved = S.QFileDialog
        S.QFileDialog = MagicMock()
        S.QFileDialog.getSaveFileName.return_value = (out, "CSV Files (*.csv)")
        try:
            dlg.export_csv()
        finally:
            S.QFileDialog = saved

    def test_export_all_blocks_concatenates_with_cumulative_index(self):
        dlg = _bare_dialog(_traces(), current=-1)
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "all.csv")
            self._run_export(dlg, out)
            with open(out, newline="", encoding="utf-8") as f:
                rows = list(csv.reader(f))
        self.assertEqual(
            rows[0], ["Cumulative Iteration", "Block", "Internal Iter", "Energy (Eh)"]
        )
        # 3 data rows total (2 + 1), cumulative index 1..3
        self.assertEqual([r[0] for r in rows[1:]], ["1", "2", "3"])
        self.assertEqual(rows[1][1], "Geom 1")
        self.assertEqual(rows[3][1], "Geom 2")
        self.assertEqual(rows[3][3], "-76.2")

    def test_export_single_block(self):
        dlg = _bare_dialog(_traces(), current=0)
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "one.csv")
            self._run_export(dlg, out)
            with open(out, newline="", encoding="utf-8") as f:
                rows = list(csv.reader(f))
        self.assertEqual(rows[0], ["Iteration", "Energy (Eh)"])
        self.assertEqual([r[0] for r in rows[1:]], ["1", "2"])
        self.assertEqual(rows[1][1], "-76.0")

    def test_export_none_selection_is_noop(self):
        dlg = _bare_dialog(_traces(), current=None)
        # getSaveFileName must never be reached; if it were, this would still
        # be safe, but assert the early return by checking no file is written.
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "never.csv")
            self._run_export(dlg, out)
            self.assertFalse(os.path.exists(out))

    def test_export_cancelled_writes_nothing(self):
        dlg = _bare_dialog(_traces(), current=0)
        dlg.parent = MagicMock(
            return_value=types.SimpleNamespace(file_path="job.out", context=None)
        )
        saved = S.QFileDialog
        S.QFileDialog = MagicMock()
        S.QFileDialog.getSaveFileName.return_value = ("", "")
        try:
            dlg.export_csv()  # must not raise
        finally:
            S.QFileDialog = saved


# ---------------------------------------------------------------------------
# update_plot guards
# ---------------------------------------------------------------------------


class TestUpdatePlot(unittest.TestCase):
    def test_none_selection_returns_before_clearing(self):
        dlg = _bare_dialog(_traces(), current=None)
        dlg.update_plot()
        dlg.ax.clear.assert_not_called()

    def test_all_blocks_plots_concatenated_series(self):
        dlg = _bare_dialog(_traces(), current=-1)
        dlg.update_plot()
        dlg.ax.plot.assert_called()
        xs, ys = dlg.ax.plot.call_args[0][0], dlg.ax.plot.call_args[0][1]
        self.assertEqual(xs, [1, 2, 3])
        self.assertEqual(ys, [-76.0, -76.1, -76.2])

    def test_single_block_plots_that_block(self):
        dlg = _bare_dialog(_traces(), current=1)
        dlg.update_plot()
        xs, ys = dlg.ax.plot.call_args[0][0], dlg.ax.plot.call_args[0][1]
        self.assertEqual(xs, [1])
        self.assertEqual(ys, [-76.2])

    def test_out_of_range_index_returns(self):
        dlg = _bare_dialog(_traces(), current=99)
        dlg.update_plot()
        dlg.ax.plot.assert_not_called()

    def test_empty_block_shows_message(self):
        dlg = _bare_dialog([{"step": "Empty", "iterations": []}], current=0)
        dlg.update_plot()
        dlg.ax.plot.assert_not_called()
        dlg.ax.text.assert_called()


if __name__ == "__main__":
    unittest.main()
