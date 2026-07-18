"""
tests/test_dipole_analysis.py
Unit tests for DipoleDialog's pure logic: the update_view arrow builder
(center-of-mass, reverse toggle, zero-magnitude guard, actor clearing) and the
settings load/save round-trip.

PyQt6 is stubbed so the module imports headlessly; numpy/pyvista are real.
The dialog's Qt constructor is bypassed (object.__new__) with lightweight
fakes for the checkboxes/spin boxes and a recording plotter. pyvista.Arrow is
patched to capture the direction/scale it is handed.
"""

import os
import sys
import json
import types
import tempfile
import importlib.util
import unittest
from unittest.mock import MagicMock

import numpy as np


# ---------------------------------------------------------------------------
# Stub PyQt6 before importing dipole_analysis (numpy/pyvista are real)
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
        "QLabel",
        "QPushButton",
        "QCheckBox",
        "QDoubleSpinBox",
        "QColorDialog",
        "QSpinBox",
        "QGroupBox",
    ]:
        setattr(qtw, cls, type(cls, (_Base,), {}))
    sys.modules["PyQt6.QtGui"].QColor = MagicMock()


_install_stubs()

_PKG_DIR = os.path.normpath(
    os.path.join(os.path.dirname(__file__), "..", "orca_result_analyzer")
)

_PKG = "_dipole_test_pkg"
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
    f"{_PKG}.dipole_analysis", os.path.join(_PKG_DIR, "dipole_analysis.py")
)
D = importlib.util.module_from_spec(_spec)
D.__package__ = _PKG
_spec.loader.exec_module(D)


# ---------------------------------------------------------------------------
# Fakes
# ---------------------------------------------------------------------------


class _FakeCheck:
    def __init__(self, checked):
        self._checked = checked

    def isChecked(self):
        return self._checked

    def setChecked(self, v):
        self._checked = bool(v)


class _FakeSpin:
    def __init__(self, value):
        self._value = value

    def value(self):
        return self._value

    def setValue(self, v):
        self._value = v


class _RecordingPlotter:
    def __init__(self):
        self.added = []
        self.removed = []
        self.rendered = 0

    def add_mesh(self, mesh, **k):
        actor = types.SimpleNamespace(mesh=mesh, kwargs=k)
        self.added.append(actor)
        return actor

    def remove_actor(self, actor):
        self.removed.append(actor)

    def render(self):
        self.rendered += 1


def _make_host(coords):
    plotter = _RecordingPlotter()
    return types.SimpleNamespace(
        mw=types.SimpleNamespace(plotter=plotter),
        parser=types.SimpleNamespace(data={"coords": coords}),
    )


def _bare_dialog(dipole_data, coords, show=True, reverse=False, scale=2.0, res=20):
    inst = D.DipoleDialog.__new__(D.DipoleDialog)
    inst.dipole_data = dipole_data
    inst.parent_dlg = _make_host(coords)
    inst.arrow_actor = None
    inst.arrow_color = "cyan"
    inst.arrow_res = res
    inst.chk_show = _FakeCheck(show)
    inst.chk_reverse = _FakeCheck(reverse)
    inst.spin_scale = _FakeSpin(scale)
    inst.spin_res = _FakeSpin(res)
    return inst


# ---------------------------------------------------------------------------
# update_view
# ---------------------------------------------------------------------------


class TestUpdateView(unittest.TestCase):
    def setUp(self):
        self._saved_arrow = D.pv.Arrow
        self.arrow_calls = []

        def _fake_arrow(*a, **k):
            self.arrow_calls.append(k)
            return object()

        D.pv.Arrow = _fake_arrow

    def tearDown(self):
        D.pv.Arrow = self._saved_arrow

    def test_arrow_placed_at_center_of_mass(self):
        dlg = _bare_dialog(
            {"vector": [1.0, 0.0, 0.0]}, coords=[[0.0, 0, 0], [2.0, 0, 0]]
        )
        dlg.update_view()
        self.assertEqual(len(self.arrow_calls), 1)
        np.testing.assert_allclose(self.arrow_calls[0]["start"], [1.0, 0.0, 0.0])

    def test_reverse_flips_direction(self):
        fwd = _bare_dialog({"vector": [0.0, 0.0, 3.0]}, coords=[[0, 0, 0]], reverse=False)
        fwd.update_view()
        rev = _bare_dialog({"vector": [0.0, 0.0, 3.0]}, coords=[[0, 0, 0]], reverse=True)
        rev.update_view()
        np.testing.assert_allclose(self.arrow_calls[0]["direction"], [0, 0, 1])
        np.testing.assert_allclose(self.arrow_calls[1]["direction"], [0, 0, -1])

    def test_scale_sets_arrow_length(self):
        dlg = _bare_dialog({"vector": [0.0, 0.0, 4.0]}, coords=[[0, 0, 0]], scale=2.0)
        dlg.update_view()
        # length = magnitude(4) * scale(2) = 8
        self.assertAlmostEqual(self.arrow_calls[0]["scale"], 8.0)

    def test_zero_magnitude_draws_no_arrow(self):
        dlg = _bare_dialog({"vector": [0.0, 0.0, 0.0]}, coords=[[0, 0, 0]])
        dlg.update_view()
        self.assertEqual(len(self.arrow_calls), 0)
        self.assertIsNone(dlg.arrow_actor)

    def test_no_coords_uses_origin(self):
        dlg = _bare_dialog({"vector": [1.0, 0.0, 0.0]}, coords=[])
        dlg.update_view()
        np.testing.assert_allclose(self.arrow_calls[0]["start"], [0.0, 0.0, 0.0])

    def test_show_adds_actor_to_plotter(self):
        dlg = _bare_dialog({"vector": [1.0, 0.0, 0.0]}, coords=[[0, 0, 0]], show=True)
        dlg.update_view()
        self.assertEqual(len(dlg.parent_dlg.mw.plotter.added), 1)
        self.assertIsNotNone(dlg.arrow_actor)

    def test_unchecked_clears_existing_actor_and_returns(self):
        dlg = _bare_dialog({"vector": [1.0, 0.0, 0.0]}, coords=[[0, 0, 0]], show=True)
        dlg.update_view()  # draws
        prev = dlg.arrow_actor
        dlg.chk_show.setChecked(False)
        dlg.update_view()  # should remove and not re-add
        self.assertIn(prev, dlg.parent_dlg.mw.plotter.removed)
        self.assertIsNone(dlg.arrow_actor)

    def test_redraw_removes_previous_actor_first(self):
        dlg = _bare_dialog({"vector": [1.0, 0.0, 0.0]}, coords=[[0, 0, 0]], show=True)
        dlg.update_view()
        first = dlg.arrow_actor
        dlg.update_view()
        self.assertIn(first, dlg.parent_dlg.mw.plotter.removed)


# ---------------------------------------------------------------------------
# settings load / save round-trip
# ---------------------------------------------------------------------------


class TestSettings(unittest.TestCase):
    def test_save_then_load_round_trips(self):
        with tempfile.TemporaryDirectory() as d:
            settings = os.path.join(d, "settings.json")

            dlg = _bare_dialog(
                {"vector": [1, 0, 0]}, coords=[[0, 0, 0]],
                show=True, reverse=False, res=42,
            )
            dlg.arrow_color = "#ff0000"
            dlg.settings_file = settings
            dlg.save_settings()

            with open(settings, encoding="utf-8") as f:
                on_disk = json.load(f)["dipole_settings"]
            self.assertEqual(on_disk["res"], 42)
            self.assertEqual(on_disk["color"], "#ff0000")
            self.assertTrue(on_disk["show"])
            self.assertFalse(on_disk["reverse"])

            other = _bare_dialog(
                {"vector": [1, 0, 0]}, coords=[[0, 0, 0]],
                show=False, reverse=True, res=20,
            )
            other.btn_color = MagicMock()
            other.settings_file = settings
            other.load_settings()
            self.assertEqual(other.spin_res.value(), 42)
            self.assertEqual(other.arrow_color, "#ff0000")
            self.assertTrue(other.chk_show.isChecked())
            self.assertFalse(other.chk_reverse.isChecked())

    def test_save_preserves_unrelated_sections(self):
        with tempfile.TemporaryDirectory() as d:
            settings = os.path.join(d, "settings.json")
            with open(settings, "w", encoding="utf-8") as f:
                json.dump({"thermal_settings": {"show_details": True}}, f)

            dlg = _bare_dialog({"vector": [1, 0, 0]}, coords=[[0, 0, 0]])
            dlg.settings_file = settings
            dlg.save_settings()

            with open(settings, encoding="utf-8") as f:
                on_disk = json.load(f)
        self.assertTrue(on_disk["thermal_settings"]["show_details"])
        self.assertIn("dipole_settings", on_disk)

    def test_load_missing_file_is_noop(self):
        dlg = _bare_dialog({"vector": [1, 0, 0]}, coords=[[0, 0, 0]])
        dlg.btn_color = MagicMock()
        dlg.settings_file = os.path.join(
            tempfile.gettempdir(), "dipole_does_not_exist_xyz.json"
        )
        dlg.load_settings()  # must not raise


if __name__ == "__main__":
    unittest.main()
