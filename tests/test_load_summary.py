"""
tests/test_load_summary.py

Tests for build_load_summary() (status-bar summary shown after a file is
loaded) and its wiring into OrcaResultAnalyzerDialog.load_file().
"""

import sys
import unittest
from unittest.mock import MagicMock, patch

# Stub Qt (mirrors tests/test_gui_reset_on_load.py's bootstrap so this file is
# self-sufficient regardless of test collection order).
if "PyQt6.QtWidgets" in sys.modules:
    qtw = sys.modules["PyQt6.QtWidgets"]
else:
    qtw = MagicMock()
    sys.modules["PyQt6.QtWidgets"] = qtw

if not hasattr(qtw, "QLabel") or "Mock" in type(getattr(qtw, "QLabel")).__name__:

    class _QLabel:
        def __init__(self, *a, **k):
            pass

        def setText(self, *a):
            pass

        def fontMetrics(self):
            fm = MagicMock()
            fm.elidedText.return_value = ""
            return fm

        def width(self):
            return 100

        def resizeEvent(self, *a):
            pass

        def __getattr__(self, name):
            if name.startswith("_"):
                raise AttributeError(name)
            return MagicMock()

    qtw.QLabel = _QLabel

if (
    not hasattr(qtw, "QPushButton")
    or "Mock" in type(getattr(qtw, "QPushButton")).__name__
):

    class _QPushButton:
        def __init__(self, *a, **k):
            self.clicked = MagicMock()

        def setIcon(self, *a):
            pass

        def setIconSize(self, *a):
            pass

        def setStyleSheet(self, *a):
            pass

        def setEnabled(self, *a):
            pass

        def setToolTip(self, *a):
            pass

        def installEventFilter(self, *a):
            pass

        def __getattr__(self, name):
            if name.startswith("_"):
                raise AttributeError(name)
            return MagicMock()

    qtw.QPushButton = _QPushButton

if not hasattr(qtw, "QDialog") or "Mock" in type(getattr(qtw, "QDialog")).__name__:

    class _QDialog:
        class DialogCode:
            Accepted = 1
            Rejected = 0

        def __init__(self, *a, **k):
            self.menu_bar_mock = MagicMock()

        def menuBar(self):
            return self.menu_bar_mock

        def setWindowTitle(self, *a):
            pass

        def resize(self, *a):
            pass

        def setAcceptDrops(self, *a):
            pass

        def close(self):
            pass

        def __getattr__(self, name):
            if name.startswith("_"):
                raise AttributeError(name)
            return MagicMock()

    qtw.QDialog = _QDialog

if "PyQt6" not in sys.modules:
    sys.modules["PyQt6"] = MagicMock()
if "PyQt6.QtCore" not in sys.modules:
    sys.modules["PyQt6.QtCore"] = MagicMock()
qt_core = sys.modules["PyQt6.QtCore"]
if (
    not hasattr(qt_core, "QObject")
    or "Mock" in type(getattr(qt_core, "QObject")).__name__
):

    class _QObject:
        def __init__(self, *a, **k):
            pass

        def __getattr__(self, name):
            if name.startswith("_"):
                raise AttributeError(name)
            return MagicMock()

    qt_core.QObject = _QObject
if "PyQt6.QtGui" not in sys.modules:
    sys.modules["PyQt6.QtGui"] = MagicMock()

if "matplotlib.backends.backend_qtagg" not in sys.modules:
    sys.modules["matplotlib.backends.backend_qtagg"] = MagicMock()
if "pyvista" not in sys.modules:
    sys.modules["pyvista"] = MagicMock()

from orca_result_analyzer.gui import (  # noqa: E402
    OrcaResultAnalyzerDialog,
    build_load_summary,
)
from orca_result_analyzer.parser import OrcaParser  # noqa: E402


class TestBuildLoadSummary(unittest.TestCase):
    def test_no_termination_status(self):
        msg, count = build_load_summary({}, "job.out")
        self.assertEqual(msg, "Loaded: job.out")
        self.assertEqual(count, 0)

    def test_missing_frequencies_key(self):
        data = {"termination_status": "Terminated normally"}
        msg, count = build_load_summary(data, "sp.out")
        self.assertEqual(msg, "Loaded: sp.out | Terminated normally")
        self.assertEqual(count, 0)

    def test_empty_frequencies_list(self):
        data = {"termination_status": "Terminated normally", "frequencies": []}
        msg, count = build_load_summary(data, "sp.out")
        self.assertEqual(msg, "Loaded: sp.out | Terminated normally")
        self.assertEqual(count, 0)

    def test_all_positive_freqs_no_imaginary(self):
        data = {
            "termination_status": "Terminated normally",
            "frequencies": [{"freq": 100.0}, {"freq": 200.0}],
        }
        msg, count = build_load_summary(data, "job.out")
        self.assertEqual(
            msg, "Loaded: job.out | Terminated normally | no imaginary modes"
        )
        self.assertEqual(count, 0)

    def test_zero_valued_modes_not_counted(self):
        # ORCA freq lists often carry leading 0.00 rot/trans entries.
        data = {
            "termination_status": "Terminated normally",
            "frequencies": [
                {"freq": 0.0},
                {"freq": 0.0},
                {"freq": 0.0},
                {"freq": 0.0},
                {"freq": 0.0},
                {"freq": 0.0},
                {"freq": 50.0},
            ],
        }
        msg, count = build_load_summary(data, "job.out")
        self.assertIn("no imaginary modes", msg)
        self.assertEqual(count, 0)

    def test_single_imaginary_mode(self):
        data = {
            "termination_status": "Terminated normally",
            "frequencies": [{"freq": 0.0}, {"freq": -123.45}, {"freq": 300.0}],
        }
        msg, count = build_load_summary(data, "ts_search.out")
        self.assertEqual(
            msg,
            "Loaded: ts_search.out | Terminated normally | "
            "1 imaginary mode (-123.45 cm-1)",
        )
        self.assertEqual(count, 1)

    def test_multiple_imaginary_modes_lowest_reported(self):
        data = {
            "termination_status": "Error termination",
            "frequencies": [
                {"freq": -450.12},
                {"freq": -20.0},
                {"freq": 300.0},
            ],
        }
        msg, count = build_load_summary(data, "broken.out")
        self.assertEqual(
            msg,
            "Loaded: broken.out | Error termination | "
            "2 imaginary modes (lowest: -450.12 cm-1)",
        )
        self.assertEqual(count, 2)

    def test_no_frequency_part_when_missing(self):
        data = {"termination_status": "Terminated normally"}
        msg, _ = build_load_summary(data, "sp.out")
        self.assertNotIn("imaginary", msg)
        self.assertNotIn("cm-1", msg)

    def test_freq_key_missing_defaults_to_zero(self):
        data = {
            "termination_status": "Terminated normally",
            "frequencies": [{}, {"freq": 100.0}],
        }
        msg, count = build_load_summary(data, "job.out")
        self.assertIn("no imaginary modes", msg)
        self.assertEqual(count, 0)


class TestLoadFileStatusMessage(unittest.TestCase):
    """Drives the real OrcaResultAnalyzerDialog.load_file() success path
    against a real temp file, mocking only the parser instantiation and the
    3D/UI update methods -- following the mocking style used in
    tests/test_gui_reset_on_load.py."""

    def _make_dialog(self, ctx):
        parser = OrcaParser()
        dlg = OrcaResultAnalyzerDialog(None, parser, "", ctx)
        dlg.close_all_sub_dialogs = MagicMock()
        dlg.update_file_info_labels = MagicMock()
        dlg.load_structure_3d = MagicMock()
        dlg.update_button_states = MagicMock()
        return dlg

    def _load_with_data(self, data):
        import os
        import tempfile

        ctx = MagicMock()
        dlg = self._make_dialog(ctx)

        fake_parser_instance = MagicMock()
        fake_parser_instance.data = data
        fake_parser_instance.load_from_memory = MagicMock()

        fd, path = tempfile.mkstemp(suffix=".out")
        os.close(fd)
        try:
            with patch(
                "orca_result_analyzer.gui.clear_atom_color_overrides"
            ), patch(
                "orca_result_analyzer.gui.OrcaParser",
                return_value=fake_parser_instance,
            ):
                dlg.load_file(path)
        finally:
            os.remove(path)
        return dlg, ctx, os.path.basename(path)

    def test_status_message_includes_imaginary_and_long_duration(self):
        data = {
            "termination_status": "Terminated normally",
            "frequencies": [{"freq": 0.0}, {"freq": -123.45}],
        }
        dlg, ctx, basename = self._load_with_data(data)
        ctx.show_status_message.assert_called_once()
        args, kwargs = ctx.show_status_message.call_args
        self.assertIn("imaginary", args[0])
        self.assertIn(basename, args[0])
        self.assertEqual(args[1], 10000)

    def test_status_message_no_imaginary_uses_default_duration(self):
        data = {
            "termination_status": "Terminated normally",
            "frequencies": [{"freq": 100.0}],
        }
        dlg, ctx, basename = self._load_with_data(data)
        ctx.show_status_message.assert_called_once()
        args, kwargs = ctx.show_status_message.call_args
        self.assertIn("no imaginary modes", args[0])
        self.assertEqual(args[1], 5000)


if __name__ == "__main__":
    unittest.main()
