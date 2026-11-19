"""
tests/test_thermal_analysis.py
Unit tests for ThermalTableDialog's pure logic: the update_table row builder,
copy_table, export_csv, and the settings load/save round-trip.

PySide6 is stubbed so the module imports headlessly; the dialog's Qt constructor
is bypassed (object.__new__) and lightweight fakes stand in for the checkbox
and table so the row-generation logic can be exercised directly.
"""

import os
import sys
import csv
import json
import tempfile
import unittest
from unittest.mock import MagicMock

sys.path.insert(0, os.path.dirname(__file__))
import gui_harness  # noqa: E402


class _FakeItem:
    """A QTableWidgetItem stand-in that actually stores its text."""

    def __init__(self, text=""):
        self._text = text
        self._row = self._col = 0

    def text(self):
        return self._text

    def row(self):
        return self._row

    def column(self):
        return self._col

    def font(self):
        return MagicMock()

    def setFont(self, _font):
        pass

    def setTextAlignment(self, _align):
        pass


# Load thermal_analysis in isolation (Qt stubs swapped in then restored, leaving
# the shared stubs other test modules rely on untouched).
T = gui_harness.load_isolated("thermal_analysis")
# The update_table / copy_table logic asserts on cell text, so give the loaded
# module a QTableWidgetItem that actually stores it.
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
        self.selected = []

    def clearSpans(self):
        pass

    def setSpan(self, *_args):
        pass

    def selectedItems(self):
        return [self._cells[rc] for rc in self.selected]

    def setRowCount(self, n):
        self._rows = n

    def rowCount(self):
        return self._rows

    def setItem(self, r, c, item):
        item._row, item._col = r, c
        self._cells[(r, c)] = item

    def item(self, r, c):
        return self._cells.get((r, c))

    def units(self):
        return {
            self._cells[(r, 0)].text(): self._cells[(r, 2)].text()
            for r in range(self._rows)
        }

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
        self.assertEqual(pairs["Electronic Energy (SP)"], "-76.12345678")
        self.assertEqual(pairs["Temperature"], "298.15")
        units = dlg.table.units()
        self.assertEqual(units["Electronic Energy (SP)"], "Eh")
        self.assertEqual(units["Temperature"], "K")

    def test_missing_values_render_as_dash(self):
        dlg = _bare_dialog({})  # nothing provided
        dlg.update_table()
        pairs = dict(dlg.table.rows_as_pairs())
        # Basic energy rows are always present; absent -> "-"
        self.assertEqual(pairs["Zero Point Energy"], "-")
        self.assertEqual(pairs["Gibbs Free Energy (G)"], "-")

    def test_value_and_unit_are_separate_cells(self):
        dlg = _bare_dialog({"gibbs": -76.5})
        dlg.update_table()
        pairs = dict(dlg.table.rows_as_pairs())
        self.assertEqual(pairs["Gibbs Free Energy (G)"], "-76.50000000")
        self.assertEqual(dlg.table.units()["Gibbs Free Energy (G)"], "Eh")

    def test_imaginary_freq_count_included_when_present(self):
        dlg = _bare_dialog({"imaginary_freq_count": 2})
        dlg.update_table()
        pairs = dict(dlg.table.rows_as_pairs())
        self.assertEqual(pairs["Imaginary Frequencies"], "2")
        self.assertEqual(dlg.table.units()["Imaginary Frequencies"], "")

    def test_imaginary_freq_absent_when_none(self):
        dlg = _bare_dialog({})
        dlg.update_table()
        labels = [p for p, _ in dlg.table.rows_as_pairs()]
        self.assertNotIn("Imaginary Frequencies", labels)

    def test_details_off_hides_breakdown(self):
        dlg = _bare_dialog({"corr_vib": -0.001}, show_details=False)
        dlg.update_table()
        labels = [p for p, _ in dlg.table.rows_as_pairs()]
        self.assertNotIn("Detailed Corrections", labels)
        self.assertNotIn("Thermal Vib Correction", labels)

    def test_details_on_shows_breakdown(self):
        dlg = _bare_dialog({"corr_vib": -0.001, "s_vib": 0.02}, show_details=True)
        dlg.update_table()
        labels = [p for p, _ in dlg.table.rows_as_pairs()]
        self.assertIn("Detailed Corrections", labels)
        self.assertIn("Thermal Vib Correction", labels)
        self.assertIn("Entropy Breakdown", labels)

    def test_row_count_matches_items(self):
        dlg = _bare_dialog({"electronic_energy": -1.0})
        dlg.update_table()
        self.assertEqual(dlg.table.rowCount(), len(dlg.table.rows_as_pairs()))


# ---------------------------------------------------------------------------
# copy_table
# ---------------------------------------------------------------------------


class TestCopyTable(unittest.TestCase):
    def _copy(self, dlg, method="copy_table"):
        captured = {}
        clip = MagicMock()
        clip.setText.side_effect = lambda t: captured.setdefault("text", t)
        saved_app = T.QApplication
        T.QApplication = MagicMock()
        T.QApplication.clipboard.return_value = clip
        try:
            getattr(dlg, method)()
        finally:
            T.QApplication = saved_app
        return captured.get("text")

    def test_copy_without_selection_copies_the_whole_table(self):
        dlg = _bare_dialog({"gibbs": -76.5, "electronic_energy": -76.0})
        dlg.update_table()
        text = self._copy(dlg)
        self.assertTrue(text.startswith("Property\tValue\tUnit\n"))
        self.assertIn("Gibbs Free Energy (G)\t-76.50000000\tEh\n", text)
        # header plus every logical row, one line each
        self.assertEqual(text.count("\n"), dlg.table.rowCount() + 1)

    def test_single_selected_cell_still_copies_the_whole_table(self):
        dlg = _bare_dialog({"gibbs": -76.5})
        dlg.update_table()
        dlg.table.selected = [(0, 1)]
        self.assertTrue(self._copy(dlg).startswith("Property\tValue\tUnit\n"))

    def test_multi_cell_selection_copies_only_those_cells(self):
        dlg = _bare_dialog({"electronic_energy": -76.0, "temperature": 298.15})
        dlg.update_table()
        # Values of the first two rows, selected out of order.
        dlg.table.selected = [(1, 1), (0, 1)]
        self.assertEqual(self._copy(dlg), "-76.00000000\n298.15\n")

    def test_selection_keeps_rows_and_columns_as_a_grid(self):
        dlg = _bare_dialog({"electronic_energy": -76.0, "temperature": 298.15})
        dlg.update_table()
        dlg.table.selected = [(0, 1), (0, 2), (1, 1), (1, 2)]
        self.assertEqual(
            self._copy(dlg, "copy_selection"), "-76.00000000\tEh\n298.15\tK\n"
        )

    def test_copy_selection_with_nothing_selected_leaves_the_clipboard(self):
        dlg = _bare_dialog({"gibbs": -76.5})
        dlg.update_table()
        self.assertIsNone(self._copy(dlg, "copy_selection"))


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
        self.assertEqual(rows[0], ["Property", "Value", "Unit"])
        flat = {r[0]: r[1:] for r in rows[1:]}
        self.assertEqual(flat["Gibbs Free Energy (G)"], ["-76.50000000", "Eh"])

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
        dlg.settings_file = os.path.join(
            tempfile.gettempdir(), "does_not_exist_xyz.json"
        )
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

    def test_save_preserves_several_other_dialogs_sections(self):
        with tempfile.TemporaryDirectory() as d:
            settings = os.path.join(d, "settings.json")
            others = {
                "mo_settings": {"iso": 0.02},
                "nmr_settings": {"ref": "TMS"},
                "dipole_settings": {"res": 20},
            }
            with open(settings, "w", encoding="utf-8") as f:
                json.dump(others, f)

            dlg = _bare_dialog({}, show_details=True)
            dlg.settings_file = settings
            dlg.save_settings()

            with open(settings, encoding="utf-8") as f:
                on_disk = json.load(f)
        for key, val in others.items():
            self.assertEqual(on_disk[key], val)
        self.assertTrue(on_disk["thermal_settings"]["show_details"])

    def test_round_trips_through_the_shared_helpers(self):
        S = gui_harness.load_isolated("settings")
        with tempfile.TemporaryDirectory() as d:
            settings = os.path.join(d, "settings.json")
            dlg = _bare_dialog({}, show_details=True)
            dlg.settings_file = settings
            dlg.save_settings()

            section = S.load_section(settings, "thermal_settings")
            self.assertTrue(section["show_details"])


if __name__ == "__main__":
    unittest.main()


class TestExportFailure(unittest.TestCase):
    def test_a_failed_export_is_reported_to_the_user(self):
        dlg = _bare_dialog({"gibbs": -76.5})
        dlg.update_table()
        with tempfile.TemporaryDirectory() as d:
            blocker = os.path.join(d, "a_file")
            with open(blocker, "w", encoding="utf-8") as fh:
                fh.write("x")
            saved_fd, saved_mb = T.QFileDialog, T.QMessageBox
            T.QFileDialog, T.QMessageBox = MagicMock(), MagicMock()
            T.QFileDialog.getSaveFileName.return_value = (
                os.path.join(blocker, "thermo.csv"),
                "",
            )
            try:
                dlg.export_csv()
                T.QMessageBox.critical.assert_called_once()
            finally:
                T.QFileDialog, T.QMessageBox = saved_fd, saved_mb
