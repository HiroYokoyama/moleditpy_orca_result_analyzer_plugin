"""
tests/test_freq_analysis.py
Unit tests for FrequencyDialog frequency list population.

Covers:
  - populate_list(): imaginary modes (raw freq < 0) appear red; real modes use
    default colour; the new Unscaled column always shows the raw value.
  - update_data(): re-applies (or removes) red colour based on the unscaled
    value whenever scaling coefficients change; only the Scaled column (col 1)
    is updated; the Unscaled column (col 2) stays fixed.

PyQt6, pyvista, rdkit, numpy, and PIL are stubbed so the tests run headlessly.
"""

import os
import sys
import types
import importlib.util
import unittest
from unittest.mock import MagicMock

_SRC_DIR = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))


# ---------------------------------------------------------------------------
# Stubs
# ---------------------------------------------------------------------------


def _install_stubs():
    if "freq_analysis_test_stubs_installed" in sys.modules:
        return

    class _QBase:
        def __init__(self, *a, **kw):
            pass

        def closeEvent(self, event):
            pass

        def __getattr__(self, name):
            return MagicMock()

    class _QDialog(_QBase):
        pass

    class _QWidget(_QBase):
        pass

    class _QTimer:
        def __init__(self, *a, **kw):
            self._interval = 0

        def timeout(self):
            pass

        def start(self, ms=0):
            self._interval = ms

        def stop(self):
            pass

        def isActive(self):
            return False

    class _QColor:
        """Minimal QColor stub that tracks name/invalid state."""

        def __init__(self, name=None):
            self._name = name
            self._valid = name is not None

        def name(self):
            return self._name or ""

        def isValid(self):
            return self._valid

        def red(self):
            return 0

        def green(self):
            return 0

        def blue(self):
            return 0

    # QTreeWidgetItem: stores column texts and foreground colours
    class _QTreeWidgetItem:
        def __init__(self, texts=None):
            self._texts = list(texts) if texts else []
            self._fgs = {}

        def text(self, col):
            return self._texts[col] if col < len(self._texts) else ""

        def setText(self, col, val):
            while len(self._texts) <= col:
                self._texts.append("")
            self._texts[col] = val

        def setForeground(self, col, color):
            self._fgs[col] = color

        def foreground(self, col):
            return self._fgs.get(col, _QColor())

    class _QTreeWidget:
        def __init__(self, *a, **kw):
            self._items = []
            self._header_labels = []
            self._current_item = None

        def setHeaderLabels(self, labels):
            self._header_labels = list(labels)

        def setEditTriggers(self, *a):
            pass

        def currentItemChanged(self):
            return MagicMock()

        def addTopLevelItem(self, item):
            self._items.append(item)

        def clear(self):
            self._items = []

        def invisibleRootItem(self):
            class _Root:
                def __init__(self, items):
                    self._items = items

                def childCount(self):
                    return len(self._items)

                def child(self, i):
                    return self._items[i]

            return _Root(self._items)

        def setCurrentItem(self, item):
            self._current_item = item

        def scrollToItem(self, *a):
            pass

        def clearSelection(self):
            pass

        def selectedItems(self):
            return []

    class _QDoubleSpinBox:
        def __init__(self, *a, **kw):
            self._val = 0.0
            self._enabled = True

        def value(self):
            return self._val

        def setValue(self, v):
            self._val = v

        def setRange(self, *a):
            pass

        def setSingleStep(self, *a):
            pass

        def setDecimals(self, *a):
            pass

        def setEnabled(self, v):
            self._enabled = v

        def setFixedWidth(self, *a):
            pass

        def setSuffix(self, *a):
            pass

        def valueChanged(self):
            return MagicMock()

        def blockSignals(self, *a):
            pass

    class _QSpinBox(_QDoubleSpinBox):
        def __init__(self, *a, **kw):
            super().__init__()
            self._val = 0

    class _QCheckBox:
        def __init__(self, *a, **kw):
            self._checked = False

        def isChecked(self):
            return self._checked

        def setChecked(self, v):
            self._checked = v

        def stateChanged(self):
            return MagicMock()

    existing_widgets = sys.modules.get("PyQt6.QtWidgets")
    existing_core = sys.modules.get("PyQt6.QtCore")
    existing_gui = sys.modules.get("PyQt6.QtGui")
    existing_pyqt6 = sys.modules.get("PyQt6")

    _WIDGET_MOCKS = [
        "QVBoxLayout",
        "QHBoxLayout",
        "QPushButton",
        "QLabel",
        "QFileDialog",
        "QMessageBox",
        "QFormLayout",
        "QDialogButtonBox",
        "QGroupBox",
        "QAbstractItemView",
        "QTreeWidgetItemIterator",
        "QComboBox",
        "QApplication",
        "QColorDialog",
        "QRadioButton",
    ]
    _WIDGET_CLASSES = {
        "QDialog": _QDialog,
        "QWidget": _QWidget,
        "QTreeWidget": _QTreeWidget,
        "QTreeWidgetItem": _QTreeWidgetItem,
        "QDoubleSpinBox": _QDoubleSpinBox,
        "QSpinBox": _QSpinBox,
        "QCheckBox": _QCheckBox,
    }

    if existing_widgets is not None:
        # Cooperative patching
        for name in _WIDGET_MOCKS:
            if not hasattr(existing_widgets, name):
                setattr(existing_widgets, name, MagicMock)
        for name, cls in _WIDGET_CLASSES.items():
            if not hasattr(existing_widgets, name):
                setattr(existing_widgets, name, cls)

        if existing_core is not None:
            if not hasattr(existing_core, "QTimer"):
                existing_core.QTimer = _QTimer
            if not hasattr(existing_core, "Qt"):
                existing_core.Qt = MagicMock()

        if existing_gui is not None:
            if not hasattr(existing_gui, "QColor"):
                existing_gui.QColor = _QColor
            if not hasattr(existing_gui, "QBrush"):
                existing_gui.QBrush = MagicMock
    else:
        # Build new stubs
        _pyqt6 = existing_pyqt6 or types.ModuleType("PyQt6")
        _core = existing_core or types.ModuleType("PyQt6.QtCore")
        _core.Qt = MagicMock()
        _core.QTimer = _QTimer
        _widgets = types.ModuleType("PyQt6.QtWidgets")
        for name in _WIDGET_MOCKS:
            setattr(_widgets, name, MagicMock)
        for name, cls in _WIDGET_CLASSES.items():
            setattr(_widgets, name, cls)
        _pyqt6.QtWidgets = _widgets
        _pyqt6.QtCore = _core
        _gui = existing_gui or types.ModuleType("PyQt6.QtGui")
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

    # pyvista stub
    _pv = types.ModuleType("pyvista")
    _pv.PolyData = MagicMock
    _pv.Arrow = MagicMock

    # spectrum_widget stub
    _sw = types.ModuleType("orca_result_analyzer.spectrum_widget")
    _sw.SpectrumWidget = MagicMock

    # utils stub
    _ut = types.ModuleType("orca_result_analyzer.utils")
    _ut.get_default_export_path = MagicMock(return_value="")

    # PIL stub
    _pil = types.ModuleType("PIL")
    _pil.Image = MagicMock

    # numpy stub (keep real numpy if available, else stub)
    try:
        import numpy  # noqa: F401
    except ImportError:
        _np = types.ModuleType("numpy")
        _np.array = list
        sys.modules["numpy"] = _np

    sys.modules.update(
        {
            "pyvista": _pv,
            "orca_result_analyzer.spectrum_widget": _sw,
            "orca_result_analyzer.utils": _ut,
            "PIL": _pil,
            "freq_analysis_test_stubs_installed": types.ModuleType("marker"),
        }
    )


def _load_freq_analysis():
    _install_stubs()
    path = os.path.join(_SRC_DIR, "orca_result_analyzer", "freq_analysis.py")
    spec = importlib.util.spec_from_file_location("freq_analysis_mod", path)
    mod = importlib.util.module_from_spec(spec)
    mod.__package__ = "orca_result_analyzer"
    sys.modules["freq_analysis_mod"] = mod
    spec.loader.exec_module(mod)
    return mod


_fa_mod = _load_freq_analysis()
FrequencyDialog = _fa_mod.FrequencyDialog


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_REAL_COLOR = "#cc0000"  # the red used for imaginary modes


def _make_dialog(frequencies):
    """Instantiate FrequencyDialog with minimal mocked parent."""
    mw = MagicMock()
    mw.plotter = MagicMock()
    dlg = FrequencyDialog.__new__(FrequencyDialog)
    # Minimal attribute init (bypass Qt __init__)
    dlg.mw = mw
    dlg.frequencies = frequencies
    dlg.atoms = []
    dlg.base_coords = []
    dlg.is_playing = False
    dlg.vector_color = "orange"
    dlg.vector_res = 20
    dlg.scaling_a = 1.0
    dlg.scaling_b = 0.0
    dlg.vector_actor = None
    dlg.spectrum_win = None
    dlg.default_presets = {"Unscaled": {"a": 1.0, "b": 0.0}}
    dlg.custom_presets = {}

    # Spin boxes
    dlg.spin_sf_a = sys.modules["PyQt6.QtWidgets"].QDoubleSpinBox()
    dlg.spin_sf_a._val = 1.0
    dlg.spin_sf_b = sys.modules["PyQt6.QtWidgets"].QDoubleSpinBox()
    dlg.spin_sf_b._val = 0.0

    # Tree widget (real stub class)
    dlg.tree = sys.modules["PyQt6.QtWidgets"].QTreeWidget()

    return dlg


# ---------------------------------------------------------------------------
# Tests: populate_list
# ---------------------------------------------------------------------------


class TestPopulateListColouring(unittest.TestCase):
    """populate_list() colours imaginary modes red using the unscaled value."""

    def _get_items(self, dlg):
        root = dlg.tree.invisibleRootItem()
        return [root.child(i) for i in range(root.childCount())]

    def test_positive_freq_not_red(self):
        freqs = [{"freq": 3000.0, "ir": 10.0, "raman": 1.0}]
        dlg = _make_dialog(freqs)
        dlg.populate_list()
        items = self._get_items(dlg)
        self.assertEqual(len(items), 1)
        fg = items[0].foreground(0)
        self.assertNotEqual(fg.name(), _REAL_COLOR)

    def test_negative_freq_is_red(self):
        freqs = [{"freq": -50.0, "ir": 0.0, "raman": 0.0}]
        dlg = _make_dialog(freqs)
        dlg.populate_list()
        items = self._get_items(dlg)
        self.assertEqual(len(items), 1)
        for col in range(5):
            self.assertEqual(
                items[0].foreground(col).name(),
                _REAL_COLOR,
                f"col {col} not red for imaginary mode",
            )

    def test_zero_freq_not_red(self):
        """Zero is not negative → not imaginary."""
        freqs = [{"freq": 0.0, "ir": 0.0, "raman": 0.0}]
        dlg = _make_dialog(freqs)
        dlg.populate_list()
        items = self._get_items(dlg)
        fg = items[0].foreground(0)
        self.assertNotEqual(fg.name(), _REAL_COLOR)

    def test_mixed_modes_colouring(self):
        freqs = [
            {"freq": -120.0, "ir": 0.0, "raman": 0.0},
            {"freq": 500.0, "ir": 5.0, "raman": 0.5},
            {"freq": -30.0, "ir": 0.0, "raman": 0.0},
            {"freq": 1200.0, "ir": 20.0, "raman": 2.0},
        ]
        dlg = _make_dialog(freqs)
        dlg.populate_list()
        items = self._get_items(dlg)
        self.assertEqual(len(items), 4)
        # imaginary
        self.assertEqual(items[0].foreground(0).name(), _REAL_COLOR)
        self.assertEqual(items[2].foreground(0).name(), _REAL_COLOR)
        # real
        self.assertNotEqual(items[1].foreground(0).name(), _REAL_COLOR)
        self.assertNotEqual(items[3].foreground(0).name(), _REAL_COLOR)

    def test_all_five_columns_coloured(self):
        freqs = [{"freq": -75.0, "ir": 0.0, "raman": 0.0}]
        dlg = _make_dialog(freqs)
        dlg.populate_list()
        items = self._get_items(dlg)
        for col in range(5):
            self.assertEqual(
                items[0].foreground(col).name(),
                _REAL_COLOR,
                f"column {col} should be red",
            )


class TestPopulateListUnscaledColumn(unittest.TestCase):
    """Column index 2 always shows the raw (unscaled) frequency value."""

    def _item0(self, dlg):
        return dlg.tree.invisibleRootItem().child(0)

    def test_unscaled_col_equals_raw_freq(self):
        freqs = [{"freq": 3500.5, "ir": 1.0, "raman": 0.0}]
        dlg = _make_dialog(freqs)
        dlg.populate_list()
        item = self._item0(dlg)
        self.assertAlmostEqual(float(item.text(2)), 3500.5, places=1)

    def test_unscaled_col_negative_freq(self):
        freqs = [{"freq": -100.25, "ir": 0.0, "raman": 0.0}]
        dlg = _make_dialog(freqs)
        dlg.populate_list()
        item = self._item0(dlg)
        self.assertAlmostEqual(float(item.text(2)), -100.25, places=1)

    def test_scaled_col_differs_from_unscaled(self):
        """When a!=1, col 1 (scaled) should differ from col 2 (unscaled)."""
        freqs = [{"freq": 1000.0, "ir": 5.0, "raman": 0.0}]
        dlg = _make_dialog(freqs)
        dlg.spin_sf_a._val = 0.96
        dlg.spin_sf_b._val = 10.0
        dlg.populate_list()
        item = self._item0(dlg)
        scaled = float(item.text(1))
        unscaled = float(item.text(2))
        self.assertAlmostEqual(unscaled, 1000.0, places=1)
        self.assertAlmostEqual(scaled, 1000.0 * 0.96 + 10.0, places=1)
        self.assertNotAlmostEqual(scaled, unscaled, places=1)

    def test_column_order_mode_scaled_unscaled_ir_raman(self):
        freqs = [{"freq": 2000.0, "ir": 7.5, "raman": 3.3}]
        dlg = _make_dialog(freqs)
        dlg.spin_sf_a._val = 0.90
        dlg.spin_sf_b._val = 0.0
        dlg.populate_list()
        item = self._item0(dlg)
        # col 0: mode index
        self.assertEqual(item.text(0), "0")
        # col 1: scaled
        self.assertAlmostEqual(float(item.text(1)), 2000.0 * 0.90, places=1)
        # col 2: unscaled
        self.assertAlmostEqual(float(item.text(2)), 2000.0, places=1)
        # col 3: IR
        self.assertAlmostEqual(float(item.text(3)), 7.5, places=1)
        # col 4: Raman
        self.assertAlmostEqual(float(item.text(4)), 3.3, places=1)


# ---------------------------------------------------------------------------
# Tests: update_data
# ---------------------------------------------------------------------------


class TestUpdateDataColouring(unittest.TestCase):
    """update_data() must re-colour rows using the unscaled value."""

    def _populate_then_update(self, freqs, new_a=1.0, new_b=0.0):
        dlg = _make_dialog(freqs)
        dlg.populate_list()
        dlg.spin_sf_a._val = new_a
        dlg.spin_sf_b._val = new_b
        dlg.update_data()
        root = dlg.tree.invisibleRootItem()
        return [root.child(i) for i in range(root.childCount())]

    def test_imaginary_stays_red_after_scaling(self):
        """Even if scaled value becomes positive, raw<0 → stays red."""
        freqs = [{"freq": -50.0, "ir": 0.0, "raman": 0.0}]
        # With a=-2, b=0 the scaled value would be +100, but unscaled=-50 → red
        items = self._populate_then_update(freqs, new_a=-2.0, new_b=0.0)
        for col in range(5):
            self.assertEqual(items[0].foreground(col).name(), _REAL_COLOR)

    def test_real_mode_not_red_after_scaling(self):
        freqs = [{"freq": 1000.0, "ir": 5.0, "raman": 0.0}]
        items = self._populate_then_update(freqs, new_a=0.96, new_b=10.0)
        for col in range(5):
            self.assertNotEqual(items[0].foreground(col).name(), _REAL_COLOR)

    def test_scaled_column_updates(self):
        freqs = [{"freq": 1000.0, "ir": 0.0, "raman": 0.0}]
        items = self._populate_then_update(freqs, new_a=0.95, new_b=5.0)
        self.assertAlmostEqual(float(items[0].text(1)), 1000.0 * 0.95 + 5.0, places=1)

    def test_unscaled_column_unchanged(self):
        """Col 2 must not be altered by update_data()."""
        freqs = [{"freq": 3000.0, "ir": 0.0, "raman": 0.0}]
        dlg = _make_dialog(freqs)
        dlg.populate_list()
        original_unscaled = dlg.tree.invisibleRootItem().child(0).text(2)
        # Apply a very different scaling
        dlg.spin_sf_a._val = 0.5
        dlg.spin_sf_b._val = -100.0
        dlg.update_data()
        updated_unscaled = dlg.tree.invisibleRootItem().child(0).text(2)
        self.assertEqual(original_unscaled, updated_unscaled)


if __name__ == "__main__":
    unittest.main()
