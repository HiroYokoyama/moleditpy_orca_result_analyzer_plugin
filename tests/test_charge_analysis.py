"""
tests/test_charge_analysis.py
Instantiation + method coverage for charge_analysis (GradientBar, ChargeDialog).

Drives the whole dialog under the headless Qt harness: construction (which runs
update_table), scheme/type switching, apply/reset coloring, 3D labels, CSV
export and close cleanup, against a fake OrcaResultAnalyzerDialog host.
"""

import os
import sys
import types
import tempfile
import unittest
from unittest.mock import MagicMock

sys.path.insert(0, os.path.dirname(__file__))
import gui_harness  # noqa: E402

_C = gui_harness.load_isolated("charge_analysis")
GradientBar = _C.GradientBar
ChargeDialog = _C.ChargeDialog


def _make_host(coords, n_atoms):
    v3d = types.SimpleNamespace(_plugin_color_overrides={})
    plotter = MagicMock()
    plotter.add_point_labels.return_value = object()
    plotter.add_mesh.return_value = object()
    mw = types.SimpleNamespace(
        plotter=plotter,
        view_3d_manager=v3d,
        current_mol=object(),
    )
    context = MagicMock()
    return types.SimpleNamespace(
        mw=mw,
        context=context,
        file_path=os.path.join(tempfile.gettempdir(), "job.out"),
        parser=types.SimpleNamespace(data={"coords": coords}),
    )


def _charges():
    return {
        "Mulliken": [
            {"atom_idx": 0, "atom_sym": "O", "charge": -0.42, "spin": 0.0},
            {"atom_idx": 1, "atom_sym": "H", "charge": 0.21, "spin": 0.0},
            {"atom_idx": 2, "atom_sym": "H", "charge": 0.21, "spin": 0.0},
        ],
        "Mayer": [
            {
                "atom_idx": 0,
                "atom_sym": "O",
                "charge": -0.4,
                "valency": 2.0,
                "bonded_valency": 1.9,
                "free_valency": 0.1,
            },
            {
                "atom_idx": 1,
                "atom_sym": "H",
                "charge": 0.2,
                "valency": 1.0,
                "bonded_valency": 1.0,
                "free_valency": 0.0,
            },
            {
                "atom_idx": 2,
                "atom_sym": "H",
                "charge": 0.2,
                "valency": 1.0,
                "bonded_valency": 1.0,
                "free_valency": 0.0,
            },
        ],
        "NBO": [
            {"atom_idx": 0, "atom_sym": "O", "charge": -0.5, "core": 2.0,
             "valence": 6.4, "rydberg": 0.1, "total": 8.5},
            {"atom_idx": 1, "atom_sym": "H", "charge": 0.25, "core": 0.0,
             "valence": 0.74, "rydberg": 0.01, "total": 0.75},
            {"atom_idx": 2, "atom_sym": "H", "charge": 0.25, "core": 0.0,
             "valence": 0.74, "rydberg": 0.01, "total": 0.75},
        ],
    }


COORDS = [[0.0, 0.0, 0.0], [0.76, 0.59, 0.0], [-0.76, 0.59, 0.0]]


class TestGradientBar(unittest.TestCase):
    def test_construct_and_set_colors(self):
        bar = GradientBar(None, ["red", "white", "blue"])
        bar.set_colors(["blue", "red"])
        # paintEvent / get_gradient touch MagicMock painters; just ensure no raise
        bar.get_gradient()
        bar.paintEvent(MagicMock())


class TestChargeDialog(unittest.TestCase):
    def _dialog(self, charges=None, coords=None, n=3):
        host = _make_host(coords if coords is not None else COORDS, n)
        dlg = ChargeDialog(host, charges if charges is not None else _charges())
        return dlg, host

    def test_construct_runs_update_table(self):
        dlg, _ = self._dialog()
        self.assertEqual(dlg.current_type, "Mulliken")

    def test_empty_charges(self):
        dlg, _ = self._dialog(charges={})
        self.assertEqual(dlg.current_type, "")

    def test_type_switch_rebuilds_table(self):
        dlg, _ = self._dialog()
        dlg.on_type_change("Mayer")
        self.assertEqual(dlg.current_type, "Mayer")
        dlg.on_type_change("NBO")
        self.assertEqual(dlg.current_type, "NBO")

    def test_scheme_change_updates_gradient_and_saves(self):
        dlg, _ = self._dialog()
        with tempfile.TemporaryDirectory() as d:
            dlg.settings_file = os.path.join(d, "settings.json")  # unused by save
            dlg.on_scheme_change("Blue(-) - White - Red(+)")
        self.assertEqual(dlg.current_scheme, "Blue(-) - White - Red(+)")

    def test_apply_colors_writes_overrides(self):
        dlg, host = self._dialog()
        dlg.apply_colors()
        overrides = host.mw.view_3d_manager._plugin_color_overrides
        self.assertEqual(set(overrides.keys()), {0, 1, 2})
        # each override is a hex color string
        self.assertTrue(all(v.startswith("#") for v in overrides.values()))
        host.context.draw_molecule_3d.assert_called()

    def test_apply_colors_two_color_scheme(self):
        dlg, host = self._dialog()
        dlg.on_scheme_change("Red(-) - Blue(+)")
        dlg.apply_colors()
        self.assertEqual(
            len(host.mw.view_3d_manager._plugin_color_overrides), 3
        )

    def test_apply_colors_all_zero_charges(self):
        charges = {"Z": [{"atom_idx": 0, "atom_sym": "C", "charge": 0.0}]}
        dlg, host = self._dialog(charges=charges, coords=[[0, 0, 0]], n=1)
        dlg.apply_colors()  # max_c==0 branch
        self.assertIn(0, host.mw.view_3d_manager._plugin_color_overrides)

    def test_reset_colors_clears_overrides(self):
        dlg, host = self._dialog()
        dlg.apply_colors()
        dlg.reset_colors()
        self.assertEqual(host.mw.view_3d_manager._plugin_color_overrides, {})
        host.context.show_status_message.assert_called()

    def test_toggle_labels_on_then_off(self):
        dlg, host = self._dialog()
        dlg.chk_show_labels = _Check(True)
        dlg.toggle_labels()
        self.assertTrue(host.mw.plotter.add_point_labels.called)
        dlg.chk_show_labels = _Check(False)
        dlg.toggle_labels()

    def test_toggle_labels_coord_mismatch_warns(self):
        dlg, host = self._dialog(coords=[[0, 0, 0]])  # 1 coord vs 3 atoms
        dlg.chk_show_labels = _Check(True)
        dlg.toggle_labels()  # hits the "no coordinates" warning branch

    def test_apply_colors_with_labels_enabled(self):
        dlg, host = self._dialog()
        dlg.chk_show_labels = _Check(True)
        dlg.apply_colors()
        host.mw.plotter.add_point_labels.assert_called()

    def test_export_csv_writes_file(self):
        dlg, host = self._dialog()
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "charges.csv")
            saved = _patch_savedialog(out)
            try:
                dlg.export_csv()
            finally:
                saved()
            self.assertTrue(os.path.exists(out))
            with open(out, encoding="utf-8") as f:
                head = f.readline()
        self.assertIn("Charge", head)

    def test_export_csv_cancelled(self):
        dlg, _ = self._dialog()
        saved = _patch_savedialog("")
        try:
            dlg.export_csv()
        finally:
            saved()

    def test_export_mayer_headers(self):
        dlg, _ = self._dialog()
        dlg.on_type_change("Mayer")
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "mayer.csv")
            saved = _patch_savedialog(out)
            try:
                dlg.export_csv()
            finally:
                saved()
            with open(out, encoding="utf-8") as f:
                head = f.readline()
        self.assertIn("Valency", head)

    def test_close_event_saves_and_cleans(self):
        dlg, host = self._dialog()
        dlg.apply_colors()
        dlg.closeEvent(MagicMock())
        # scalar bar + render invoked during cleanup
        host.mw.plotter.render.assert_called()

    def test_reject_calls_close(self):
        dlg, _ = self._dialog()
        dlg.close = MagicMock()
        dlg.reject()
        dlg.close.assert_called_once()


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


class _Check:
    def __init__(self, v):
        self._v = v

    def isChecked(self):
        return self._v

    def setChecked(self, v):
        self._v = bool(v)


def _patch_savedialog(return_path):
    saved = _C.QFileDialog
    _C.QFileDialog = MagicMock()
    _C.QFileDialog.getSaveFileName.return_value = (return_path, "CSV Files (*.csv)")

    def restore():
        _C.QFileDialog = saved

    return restore


if __name__ == "__main__":
    unittest.main()
