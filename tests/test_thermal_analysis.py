"""
tests/test_thermal_analysis.py
Unit tests for ThermalTableDialog's pure logic: the update_table row builder,
copy_table, export_csv, and the settings load/save round-trip.

PyQt6 is stubbed so the module imports headlessly; the dialog's Qt constructor
is bypassed (object.__new__) and lightweight fakes stand in for the checkbox
and table so the row-generation logic can be exercised directly.
"""

import os
import sys
import csv
import json
import types
import tempfile
import importlib.util
import unittest
from unittest.mock import MagicMock


# ---------------------------------------------------------------------------
# Stub PyQt6 before importing thermal_analysis
# ---------------------------------------------------------------------------


class _FakeItem:
    """A QTableWidgetItem stand-in that actually stores its text."""

    def __init__(self, text=""):
        self._text = text

    def text(self):
        return self._text


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
        "QTableWidget",
        "QHeaderView",
        "QPushButton",
        "QCheckBox",
        "QAbstractItemView",
        "QFileDialog",
    ]:
        setattr(qtw, cls, type(cls, (_Base,), {}))
    qtw.QTableWidgetItem = _FakeItem
    qtw.QApplication = MagicMock()


_install_stubs()

_PKG_DIR = os.path.normpath(
    os.path.join(os.path.dirname(__file__), "..", "orca_result_analyzer")
)


# Load thermal_analysis under a throwaway package so the lazy
# `from .utils import save_json_atomic` in save_settings resolves to the real
# utils module WITHOUT executing orca_result_analyzer/__init__ (which imports
# the full Qt GUI, unavailable under these stubs) and without shadowing the
# real package for other test files in the suite.
_PKG = "_thermal_test_pkg"


def _preload():
    if _PKG not in sys.modules:
        pkg = types.ModuleType(_PKG)
        pkg.__path__ = [_PKG_DIR]
        sys.modules[_PKG] = pkg
    if f"{_PKG}.utils" not in sys.modules:
        uspec = importlib.util.spec_from_file_location(
            f"{_PKG}.utils", os.path.join(_PKG_DIR, "utils.py")
        )
        umod = importlib.util.module_from_spec(uspec)
        uspec.loader.exec_module(umod)
        sys.modules[f"{_PKG}.utils"] = umod


_preload()

_SRC = os.path.join(_PKG_DIR, "thermal_analysis.py")
_spec = importlib.util.spec_from_file_location(f"{_PKG}.thermal_analysis", _SRC)
T = importlib.util.module_from_spec(_spec)
T.__package__ = _PKG
_spec.loader.exec_module(T)
# Ensure the module's QTableWidgetItem is our capturing fake (in case another
# test module replaced the shared PyQt6 stub after import).
T.QTableWidgetItem = _FakeItem


# ---------------------------------------------------------------------------
# Fakes
# ---------------------------------------------------------------------------


class _FakeCheck:
    def __init__(self, checked=False):
        self._checked = checked

    def isChecked(self):
        return self._checked

    def setChecked(self, v):
        self._checked = bool(v)


class _FakeTable:
    def __init__(self):
        self._rows = 0
        self._cells = {}

    def setRowCount(self, n):
        self._rows = n

    def rowCount(self):
        return self._rows

    def setItem(self, r, c, item):
        self._cells[(r, c)] = item

    def item(self, r, c):
        return self._cells.get((r, c))

    def rows_as_pairs(self):
        return [
            (self._cells[(r, 0)].text(), self._cells[(r, 1)].text())
            for r in range(self._rows)
        ]


def _bare_dialog(data, show_details=False):
    inst = T.ThermalTableDialog.__new__(T.ThermalTableDialog)
    inst.data = data
    inst.chk_details = _FakeCheck(show_details)
    inst.table = _FakeTable()
    return inst


# ---------------------------------------------------------------------------
# update_table row generation
# ---------------------------------------------------------------------------


class TestUpdateTable(unittest.TestCase):
    def test_energy_and_temperature_formatting(self):
        dlg = _bare_dialog({"electronic_energy": -76.12345678, "temperature": 298.15})
        dlg.update_table()
        pairs = dict(dlg.table.rows_as_pairs())
        self.assertEqual(pairs["Electronic Energy (SP)"], "-76.12345678 Eh")
        self.assertEqual(pairs["Temperature (K)"], "298.15 K")

    def test_missing_values_render_as_dash(self):
        dlg = _bare_dialog({})  # nothing provided
        dlg.update_table()
        pairs = dict(dlg.table.rows_as_pairs())
        # Basic energy rows are always present; absent -> "-"
        self.assertEqual(pairs["Zero Point Energy"], "-")
        self.assertEqual(pairs["Gibbs Free Energy (G)"], "-")

    def test_float_value_gets_eh_suffix(self):
        dlg = _bare_dialog({"gibbs": -76.5})
        dlg.update_table()
        pairs = dict(dlg.table.rows_as_pairs())
        self.assertEqual(pairs["Gibbs Free Energy (G)"], "-76.50000000 Eh")

    def test_imaginary_freq_count_included_when_present(self):
        dlg = _bare_dialog({"imaginary_freq_count": 2})
        dlg.update_table()
        pairs = dict(dlg.table.rows_as_pairs())
        self.assertEqual(pairs["Imaginary Frequencies"], "2")

    def test_imaginary_freq_absent_when_none(self):
        dlg = _bare_dialog({})
        dlg.update_table()
        labels = [p for p, _ in dlg.table.rows_as_pairs()]
        self.assertNotIn("Imaginary Frequencies", labels)

    def test_details_off_hides_breakdown(self):
        dlg = _bare_dialog({"corr_vib": -0.001}, show_details=False)
        dlg.update_table()
        labels = [p for p, _ in dlg.table.rows_as_pairs()]
        self.assertNotIn("--- Detailed Corrections ---", labels)
        self.assertNotIn("Thermal Vib Correction", labels)

    def test_details_on_shows_breakdown(self):
        dlg = _bare_dialog({"corr_vib": -0.001, "s_vib": 0.02}, show_details=True)
        dlg.update_table()
        labels = [p for p, _ in dlg.table.rows_as_pairs()]
        self.assertIn("--- Detailed Corrections ---", labels)
        self.assertIn("Thermal Vib Correction", labels)
        self.assertIn("--- Entropy Breakdown ---", labels)

    def test_row_count_matches_items(self):
        dlg = _bare_dialog({"electronic_energy": -1.0})
        dlg.update_table()
        self.assertEqual(dlg.table.rowCount(), len(dlg.table.rows_as_pairs()))


# ---------------------------------------------------------------------------
# copy_table
# ---------------------------------------------------------------------------


class TestCopyTable(unittest.TestCase):
    def test_copy_builds_tab_separated_text(self):
        dlg = _bare_dialog({"gibbs": -76.5, "electronic_energy": -76.0})
        dlg.update_table()
        captured = {}
        clip = MagicMock()
        clip.setText.side_effect = lambda t: captured.setdefault("text", t)
        saved_app = T.QApplication
        T.QApplication = MagicMock()
        T.QApplication.clipboard.return_value = clip
        try:
            dlg.copy_table()
        finally:
            T.QApplication = saved_app
        text = captured["text"]
        self.assertIn("Gibbs Free Energy (G)\t-76.50000000 Eh\n", text)
        # every logical row becomes one line
        self.assertEqual(text.count("\n"), dlg.table.rowCount())


# ---------------------------------------------------------------------------
# export_csv
# ---------------------------------------------------------------------------


class TestExportCsv(unittest.TestCase):
    def test_export_writes_header_and_rows(self):
        dlg = _bare_dialog({"gibbs": -76.5})
        dlg.update_table()
        dlg.parent = MagicMock(return_value=None)  # skip status-message branch

        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "thermo.csv")
            saved = T.QFileDialog
            T.QFileDialog = MagicMock()
            T.QFileDialog.getSaveFileName.return_value = (out, "CSV Files (*.csv)")
            try:
                dlg.export_csv()
            finally:
                T.QFileDialog = saved

            with open(out, newline="", encoding="utf-8") as f:
                rows = list(csv.reader(f))
        self.assertEqual(rows[0], ["Property", "Value"])
        flat = {r[0]: r[1] for r in rows[1:] if len(r) == 2}
        self.assertEqual(flat["Gibbs Free Energy (G)"], "-76.50000000 Eh")

    def test_export_cancelled_writes_nothing(self):
        dlg = _bare_dialog({"gibbs": -76.5})
        dlg.update_table()
        saved = T.QFileDialog
        T.QFileDialog = MagicMock()
        T.QFileDialog.getSaveFileName.return_value = ("", "")  # user cancelled
        try:
            dlg.export_csv()  # must not raise
        finally:
            T.QFileDialog = saved


# ---------------------------------------------------------------------------
# settings load / save round-trip
# ---------------------------------------------------------------------------


class TestSettings(unittest.TestCase):
    def test_save_then_load_round_trips_show_details(self):
        with tempfile.TemporaryDirectory() as d:
            settings = os.path.join(d, "settings.json")

            dlg = _bare_dialog({}, show_details=True)
            dlg.settings_file = settings
            dlg.save_settings()

            with open(settings, encoding="utf-8") as f:
                on_disk = json.load(f)
            self.assertTrue(on_disk["thermal_settings"]["show_details"])

            other = _bare_dialog({}, show_details=False)
            other.settings_file = settings
            other.load_settings()
            self.assertTrue(other.chk_details.isChecked())

    def test_save_preserves_unrelated_sections(self):
        with tempfile.TemporaryDirectory() as d:
            settings = os.path.join(d, "settings.json")
            with open(settings, "w", encoding="utf-8") as f:
                json.dump({"dipole_settings": {"color": "red"}}, f)

            dlg = _bare_dialog({}, show_details=True)
            dlg.settings_file = settings
            dlg.save_settings()

            with open(settings, encoding="utf-8") as f:
                on_disk = json.load(f)
        self.assertEqual(on_disk["dipole_settings"]["color"], "red")
        self.assertIn("thermal_settings", on_disk)

    def test_load_missing_file_is_noop(self):
        dlg = _bare_dialog({}, show_details=False)
        dlg.settings_file = os.path.join(tempfile.gettempdir(), "does_not_exist_xyz.json")
        dlg.load_settings()  # must not raise
        self.assertFalse(dlg.chk_details.isChecked())

    def test_load_corrupt_file_is_swallowed(self):
        with tempfile.TemporaryDirectory() as d:
            settings = os.path.join(d, "settings.json")
            with open(settings, "w", encoding="utf-8") as f:
                f.write("{not valid json")
            dlg = _bare_dialog({}, show_details=False)
            dlg.settings_file = settings
            dlg.load_settings()  # must not raise
        self.assertFalse(dlg.chk_details.isChecked())


if __name__ == "__main__":
    unittest.main()
