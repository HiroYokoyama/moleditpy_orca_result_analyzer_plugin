"""
tests/test_nmr_mo_fixes.py
Unit tests for the bug-fixes applied in v2.6.0:

  NMR:
    - sel_timer is stopped in closeEvent so it cannot fire on a dead widget.
    - reset_selection() exists and stops/restarts the timer while clearing state.

  MO:
    - apply_preset() blocks signals on all vis widgets so only ONE show_cube()
      call is made per preset application (no intermediate phantom renders).

PyQt6, pyvista, rdkit, numpy, matplotlib, PIL and all related sub-modules are
fully stubbed so tests run headlessly.
"""

import os
import sys
import types
import importlib.util
import unittest
from unittest.mock import MagicMock

_SRC_DIR = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))


# ---------------------------------------------------------------------------
# Shared stub helpers
# ---------------------------------------------------------------------------


def _make_qt_base():
    """Return a base class stub for QDialog/QWidget with noop Qt methods."""

    class _Base:
        class DialogCode:
            Accepted = 1
            Rejected = 0

        def __init__(self, *a, **kw):
            pass

        def closeEvent(self, event):
            pass

        def __getattr__(self, name):
            return MagicMock()

    return _Base


class _QColor:
    def __init__(self, name=None):
        self._name = name

    def name(self):
        return self._name or ""

    def isValid(self):
        return self._name is not None

    def red(self):
        return 0

    def green(self):
        return 0

    def blue(self):
        return 0


class _SpinBox:
    def __init__(self, val=0.0):
        self._val = val
        self._blocked = False
        self._signals_blocked = False

    def value(self):
        return self._val

    def setValue(self, v):
        self._val = v

    def blockSignals(self, v):
        self._signals_blocked = v

    def setRange(self, *a):
        pass

    def setSingleStep(self, *a):
        pass

    def setDecimals(self, *a):
        pass

    def setEnabled(self, *a):
        pass

    def setFixedWidth(self, *a):
        pass

    def setSuffix(self, *a):
        pass

    def setChecked(self, v):
        self._val = bool(v)

    def isChecked(self):
        return bool(self._val)

    def stateChanged(self):
        return MagicMock()

    def valueChanged(self):
        return MagicMock()

    def toggled(self):
        return MagicMock()

    def currentTextChanged(self):
        return MagicMock()


class _QComboBox:
    def __init__(self):
        self._items = []
        self._idx = 0
        self._blocked = False

    def addItem(self, t):
        self._items.append(t)

    def addItems(self, ts):
        self._items.extend(ts)

    def findText(self, t):
        try:
            return self._items.index(t)
        except ValueError:
            return -1

    def setCurrentIndex(self, i):
        self._idx = i

    def currentText(self):
        return self._items[self._idx] if self._items else ""

    def currentIndex(self):
        return self._idx

    def blockSignals(self, v):
        self._blocked = v

    def currentTextChanged(self):
        return MagicMock()

    def currentIndexChanged(self):
        return MagicMock()

    def clear(self):
        self._items = []
        self._idx = 0


class _RecordingTimer:
    """QTimer stub that records stop()/start() calls."""

    def __init__(self, *a, **kw):
        self.calls = []
        self._active = False

    def start(self, ms=200):
        self.calls.append(("start", ms))
        self._active = True

    def stop(self):
        self.calls.append(("stop",))
        self._active = False

    def isActive(self):
        return self._active

    def timeout(self):
        return MagicMock()


# ---------------------------------------------------------------------------
# NMR stubs + loader
# ---------------------------------------------------------------------------


def _install_nmr_stubs():
    """Install stubs for headless NMR loading.

    When the full suite runs, real PyQt6 may already be in sys.modules from
    earlier tests.  Replacing it wholesale breaks the test runner.  Instead:
      - If QtWidgets is already present, patch only the names nmr_analysis
        needs that are missing (MagicMock stand-ins).
      - If PyQt6 is absent, install a complete stub module.
    numpy / matplotlib / pyvista stubs are always installed under their real
    names so nmr_analysis can import them without a display.
    """
    marker = "nmr_fix_stubs_installed"
    if marker in sys.modules:
        return

    _Base = _make_qt_base()

    # All names nmr_analysis.py imports from QtWidgets
    _WIDGET_NAMES = [
        "QDialog",
        "QVBoxLayout",
        "QHBoxLayout",
        "QLabel",
        "QComboBox",
        "QDoubleSpinBox",
        "QTableWidget",
        "QTableWidgetItem",
        "QHeaderView",
        "QPushButton",
        "QApplication",
        "QGroupBox",
        "QMessageBox",
        "QFileDialog",
        "QCheckBox",
        "QButtonGroup",
        "QAbstractItemView",
        "QTreeWidget",
        "QTreeWidgetItem",
        "QFormLayout",
        "QDialogButtonBox",
        "QSplitter",
        "QFrame",
        "QScrollArea",
        "QSizePolicy",
        "QWidget",
        "QSpinBox",
    ]

    existing_widgets = sys.modules.get("PyQt6.QtWidgets")
    existing_core = sys.modules.get("PyQt6.QtCore")

    if existing_widgets is not None:
        # Real PyQt6 already loaded — only patch missing names
        for name in _WIDGET_NAMES:
            if not hasattr(existing_widgets, name):
                setattr(existing_widgets, name, MagicMock)
        if existing_core is not None and not hasattr(existing_core, "QTimer"):
            existing_core.QTimer = _RecordingTimer
    else:
        # No PyQt6 yet — build and install full stub
        _pyqt6 = types.ModuleType("PyQt6")
        _core = types.ModuleType("PyQt6.QtCore")
        _core.Qt = MagicMock()
        _core.QTimer = _RecordingTimer
        _widgets = types.ModuleType("PyQt6.QtWidgets")
        for name in _WIDGET_NAMES:
            setattr(_widgets, name, MagicMock)
        _widgets.QDialog = _Base
        _pyqt6.QtWidgets = _widgets
        _pyqt6.QtCore = _core
        _gui = types.ModuleType("PyQt6.QtGui")
        _gui.QColor = _QColor
        _gui.QBrush = MagicMock
        _pyqt6.QtGui = _gui
        sys.modules.update(
            {
                "PyQt6": _pyqt6,
                "PyQt6.QtWidgets": _widgets,
                "PyQt6.QtCore": _core,
                "PyQt6.QtGui": _gui,
            }
        )

    # Always install these stubs (they don't conflict with a live Qt session)
    _np = types.ModuleType("numpy")
    _np.array = list
    _np.zeros_like = lambda a: [0] * len(a)

    _mpl = types.ModuleType("matplotlib")
    _mpl.use = MagicMock()
    _mpl_fig = types.ModuleType("matplotlib.figure")
    _mpl_fig.Figure = MagicMock
    _mpl_back = types.ModuleType("matplotlib.backends")
    _mpl_back_qt = types.ModuleType("matplotlib.backends.backend_qtagg")
    _mpl_back_qt.FigureCanvasQTAgg = _Base
    _mpl_back_qt.NavigationToolbar2QT = MagicMock
    _mpl_tick = types.ModuleType("matplotlib.ticker")
    _mpl_tick.MaxNLocator = MagicMock

    _pv = types.ModuleType("pyvista")
    _pv.PolyData = MagicMock
    _pv.Sphere = MagicMock

    _utils_mod = types.ModuleType("orca_result_analyzer.utils")
    _utils_mod.get_default_export_path = MagicMock(return_value="")

    _custom_ref = types.ModuleType("orca_result_analyzer.nmr_custom_ref_dialog")
    _custom_ref.CustomReferenceDialog = MagicMock

    _ver_mod = types.ModuleType("orca_result_analyzer")
    _ver_mod.PLUGIN_VERSION = "2.6.0"

    sys.modules.update(
        {
            "numpy": _np,
            "matplotlib": _mpl,
            "matplotlib.figure": _mpl_fig,
            "matplotlib.backends": _mpl_back,
            "matplotlib.backends.backend_qtagg": _mpl_back_qt,
            "matplotlib.ticker": _mpl_tick,
            "pyvista": _pv,
            "orca_result_analyzer": _ver_mod,
            "orca_result_analyzer.utils": _utils_mod,
            "orca_result_analyzer.nmr_custom_ref_dialog": _custom_ref,
            marker: types.ModuleType(marker),
        }
    )


def _load_nmr():
    _install_nmr_stubs()
    path = os.path.join(_SRC_DIR, "orca_result_analyzer", "nmr_analysis.py")
    spec = importlib.util.spec_from_file_location("nmr_analysis_fix_mod", path)
    mod = importlib.util.module_from_spec(spec)
    mod.__package__ = "orca_result_analyzer"
    sys.modules["nmr_analysis_fix_mod"] = mod
    spec.loader.exec_module(mod)
    return mod


_nmr_mod = _load_nmr()
NMRDialog = _nmr_mod.NMRDialog


# ---------------------------------------------------------------------------
# NMR test helpers
# ---------------------------------------------------------------------------


def _make_nmr_dialog():
    """Build a minimal NMRDialog instance without running Qt __init__."""
    parent_dlg = MagicMock()
    parent_dlg.mw = MagicMock()

    dlg = NMRDialog.__new__(NMRDialog)
    dlg.parent_dlg = parent_dlg
    dlg.data = [{"atom_idx": 0, "atom_sym": "H", "shielding": 30.0}]
    dlg.couplings = []
    dlg.displayed_data = list(dlg.data)
    dlg.selected_peak_indices = set()
    dlg.highlight_artists = []
    dlg.show_all_mode = False
    dlg._atom_labels = []
    dlg._nmr_sphere_actors = []
    dlg._nmr_label_names = []
    dlg._last_synced_mw_selection = frozenset()
    dlg.peaks_metadata = None

    # Real recording timer
    timer = _RecordingTimer()
    dlg.sel_timer = timer

    return dlg


# ---------------------------------------------------------------------------
# NMR Tests
# ---------------------------------------------------------------------------


class TestNMRCloseEventStopsTimer(unittest.TestCase):
    """closeEvent() must stop sel_timer before it fires on a dead dialog."""

    def test_close_event_stops_timer(self):
        dlg = _make_nmr_dialog()
        dlg.sel_timer.start(200)  # simulate running timer
        self.assertTrue(dlg.sel_timer.isActive())

        # Patch save_settings and clear_atom_labels (they touch Qt)
        dlg.save_settings = MagicMock()
        dlg.clear_atom_labels = MagicMock()

        event = MagicMock()
        # Call the real closeEvent
        NMRDialog.closeEvent(dlg, event)

        # Timer must be stopped
        self.assertFalse(dlg.sel_timer.isActive())
        stop_calls = [c for c in dlg.sel_timer.calls if c[0] == "stop"]
        self.assertGreater(len(stop_calls), 0, "sel_timer.stop() was never called")

    def test_close_event_stop_before_save(self):
        """stop() must precede save_settings() to avoid mid-close timer fires."""
        dlg = _make_nmr_dialog()
        call_order = []
        dlg.sel_timer.stop = lambda: call_order.append("stop")
        dlg.save_settings = lambda: call_order.append("save")
        dlg.clear_atom_labels = MagicMock()
        NMRDialog.closeEvent(dlg, MagicMock())
        self.assertLess(call_order.index("stop"), call_order.index("save"))


class TestNMRResetSelection(unittest.TestCase):
    """reset_selection() stops the timer, clears state, then restarts."""

    def test_reset_selection_exists(self):
        self.assertTrue(hasattr(NMRDialog, "reset_selection"))
        self.assertTrue(callable(NMRDialog.reset_selection))

    def test_reset_selection_stops_and_restarts_timer(self):
        dlg = _make_nmr_dialog()
        dlg.sel_timer.start(200)
        dlg.clear_peak_selection = MagicMock()

        NMRDialog.reset_selection(dlg)

        ops = [c[0] for c in dlg.sel_timer.calls]
        self.assertIn("stop", ops)
        self.assertIn("start", ops)
        # stop must come before start
        last_start = len(ops) - 1 - list(reversed(ops)).index("start")
        self.assertLess(ops.index("stop"), last_start)

    def test_reset_selection_clears_sync_tracker(self):
        dlg = _make_nmr_dialog()
        dlg._last_synced_mw_selection = frozenset({1, 2, 3})
        dlg.clear_peak_selection = MagicMock()
        NMRDialog.reset_selection(dlg)
        self.assertEqual(dlg._last_synced_mw_selection, frozenset())

    def test_reset_selection_calls_clear_peak_selection(self):
        dlg = _make_nmr_dialog()
        dlg.clear_peak_selection = MagicMock()
        NMRDialog.reset_selection(dlg)
        dlg.clear_peak_selection.assert_called_once()


# ---------------------------------------------------------------------------
# MO stubs + loader
# ---------------------------------------------------------------------------


def _install_mo_stubs():
    marker = "mo_fix_stubs_installed"
    if marker in sys.modules:
        return

    _Base = _make_qt_base()

    _MO_WIDGET_NAMES = [
        "QVBoxLayout",
        "QHBoxLayout",
        "QPushButton",
        "QLabel",
        "QTreeWidget",
        "QTreeWidgetItem",
        "QAbstractItemView",
        "QTreeWidgetItemIterator",
        "QFileDialog",
        "QMessageBox",
        "QGroupBox",
        "QHeaderView",
        "QProgressDialog",
        "QFormLayout",
        "QInputDialog",
        "QApplication",
        "QColorDialog",
        "QWidget",
        "QDialog",
        "QDoubleSpinBox",
        "QSpinBox",
        "QCheckBox",
        "QComboBox",
    ]

    existing_widgets = sys.modules.get("PyQt6.QtWidgets")

    if existing_widgets is not None:
        # Patch whatever is already there
        for name in _MO_WIDGET_NAMES:
            if not hasattr(existing_widgets, name):
                if name in ["QWidget", "QDialog"]:
                    setattr(existing_widgets, name, _Base)
                else:
                    setattr(existing_widgets, name, MagicMock)
    else:
        _pyqt6_mo = types.ModuleType("PyQt6")
        _core_mo = types.ModuleType("PyQt6.QtCore")
        _core_mo.Qt = MagicMock()
        _widgets_mo = types.ModuleType("PyQt6.QtWidgets")
        for name in _MO_WIDGET_NAMES:
            if name in ["QWidget", "QDialog"]:
                setattr(_widgets_mo, name, _Base)
            else:
                setattr(_widgets_mo, name, MagicMock)
        _pyqt6_mo.QtWidgets = _widgets_mo
        _pyqt6_mo.QtCore = _core_mo
        _gui_mo = types.ModuleType("PyQt6.QtGui")
        _gui_mo.QColor = _QColor
        _gui_mo.QBrush = MagicMock
        _pyqt6_mo.QtGui = _gui_mo
        sys.modules.update(
            {
                "PyQt6": _pyqt6_mo,
                "PyQt6.QtWidgets": _widgets_mo,
                "PyQt6.QtCore": _core_mo,
                "PyQt6.QtGui": _gui_mo,
            }
        )

    _utils_mo = types.ModuleType("orca_result_analyzer.utils")
    _utils_mo.get_default_export_path = MagicMock(return_value="")
    sys.modules.setdefault("orca_result_analyzer.utils", _utils_mo)

    sys.modules[marker] = types.ModuleType(marker)


def _load_mo():
    _install_mo_stubs()
    # Stub heavy deps before loading mo_analysis
    for mod_name in [
        "orca_result_analyzer.mo_engine",
        "orca_result_analyzer.vis",
        "orca_result_analyzer.energy_diag",
    ]:
        sys.modules.setdefault(mod_name, types.ModuleType(mod_name))
    path = os.path.join(_SRC_DIR, "orca_result_analyzer", "mo_analysis.py")
    spec = importlib.util.spec_from_file_location("mo_analysis_fix_mod", path)
    mod = importlib.util.module_from_spec(spec)
    mod.__package__ = "orca_result_analyzer"
    sys.modules["mo_analysis_fix_mod"] = mod
    spec.loader.exec_module(mod)
    return mod


_mo_mod = _load_mo()
MODialog = _mo_mod.MODialog


# ---------------------------------------------------------------------------
# MO test helpers
# ---------------------------------------------------------------------------


def _make_mo_dialog():
    dlg = MODialog.__new__(MODialog)
    dlg.mw = MagicMock()
    dlg.parent_dlg = MagicMock()
    dlg.mo_data = {}
    dlg.mos = {}
    dlg.last_cube_path = None
    dlg.generation_queue = []
    dlg.energy_dlg = None
    dlg.settings_file = os.path.join(_SRC_DIR, "orca_result_analyzer", "settings.json")

    # Vis control widgets (recording spins/combos/checkboxes)
    dlg.spin_iso = _SpinBox(0.02)
    dlg.spin_opacity = _SpinBox(0.5)
    dlg.combo_style = _QComboBox()
    dlg.combo_style.addItems(["Surface", "Wireframe", "Points"])
    dlg.check_smooth = _SpinBox(1.0)  # isChecked() → True

    dlg.btn_color_p = MagicMock()
    dlg.btn_color_n = MagicMock()
    dlg.btn_color_p.styleSheet.return_value = "background-color: #ff0000;"
    dlg.btn_color_n.styleSheet.return_value = "background-color: #0000ff;"

    dlg.combo_presets = _QComboBox()
    dlg.combo_presets.addItem("Default")

    dlg.presets = {
        "Default": {
            "iso": 0.02,
            "opacity": 0.5,
            "style": "Surface",
            "color_p": "#ff0000",
            "color_n": "#0000ff",
            "smooth_shading": True,
        },
        "MyPreset": {
            "iso": 0.05,
            "opacity": 0.7,
            "style": "Wireframe",
            "color_p": "#00ff00",
            "color_n": "#ff00ff",
            "smooth_shading": False,
        },
    }
    dlg.combo_presets.addItem("MyPreset")

    return dlg


# ---------------------------------------------------------------------------
# MO Tests
# ---------------------------------------------------------------------------


class TestMOApplyPresetSignalBlock(unittest.TestCase):
    """apply_preset() must block all widget signals so only one
    update_vis_only() is fired at the end, not one per widget change."""

    def test_update_vis_only_called_exactly_once(self):
        dlg = _make_mo_dialog()
        dlg.update_vis_only = MagicMock()
        dlg.save_settings = MagicMock()
        dlg.set_btn_color = MagicMock()

        MODialog.apply_preset(dlg, "MyPreset")

        dlg.update_vis_only.assert_called_once()

    def test_spin_iso_blocked_during_apply(self):
        """spin_iso.blockSignals(True) must be called before setValue."""
        dlg = _make_mo_dialog()
        dlg.update_vis_only = MagicMock()
        dlg.save_settings = MagicMock()
        dlg.set_btn_color = MagicMock()

        block_calls = []

        def recording_bs(v):
            block_calls.append(v)

        dlg.spin_iso.blockSignals = recording_bs

        MODialog.apply_preset(dlg, "MyPreset")

        self.assertIn(True, block_calls, "blockSignals(True) never called on spin_iso")
        self.assertIn(
            False, block_calls, "blockSignals(False) never called on spin_iso"
        )

    def test_values_applied_from_preset(self):
        dlg = _make_mo_dialog()
        dlg.update_vis_only = MagicMock()
        dlg.save_settings = MagicMock()
        dlg.set_btn_color = MagicMock()

        MODialog.apply_preset(dlg, "MyPreset")

        self.assertAlmostEqual(dlg.spin_iso.value(), 0.05)
        self.assertAlmostEqual(dlg.spin_opacity.value(), 0.7)

    def test_signals_unblocked_after_apply(self):
        """Even if an exception occurs during apply, signals must be unblocked."""
        dlg = _make_mo_dialog()
        dlg.update_vis_only = MagicMock()
        dlg.save_settings = MagicMock()
        dlg.set_btn_color = MagicMock()

        # Force an error mid-apply by making combo_style.findText raise
        dlg.combo_style.findText = MagicMock(side_effect=RuntimeError("oops"))

        try:
            MODialog.apply_preset(dlg, "MyPreset")
        except RuntimeError:
            pass

        # After the try/finally, signals must be unblocked
        self.assertFalse(dlg.spin_iso._signals_blocked)
        self.assertFalse(dlg.spin_opacity._signals_blocked)
        self.assertFalse(dlg.check_smooth._signals_blocked)

    def test_unknown_preset_is_noop(self):
        dlg = _make_mo_dialog()
        dlg.update_vis_only = MagicMock()
        MODialog.apply_preset(dlg, "DoesNotExist")
        dlg.update_vis_only.assert_not_called()

    def test_show_cube_updates_last_cube_path(self):
        dlg = _make_mo_dialog()
        dlg.mw = MagicMock()
        dlg.get_color_hex = MagicMock(return_value="#ffffff")
        dlg.combo_style = MagicMock()
        dlg.combo_style.currentText = MagicMock(return_value="Surface")
        dlg.spin_iso = MagicMock()
        dlg.spin_iso.value = MagicMock(return_value=0.02)
        dlg.spin_opacity = MagicMock()
        dlg.spin_opacity.value = MagicMock(return_value=0.5)
        dlg.check_smooth = MagicMock()
        dlg.check_smooth.isChecked = MagicMock(return_value=True)

        mock_vis = MagicMock()
        mock_vis.load_file = MagicMock(return_value=True)
        _mo_mod.CubeVisualizer = MagicMock(return_value=mock_vis)

        dlg.last_cube_path = "old_path.cube"
        MODialog.show_cube(dlg, "new_path.cube")

        self.assertEqual(dlg.last_cube_path, "new_path.cube")


if __name__ == "__main__":
    unittest.main()
