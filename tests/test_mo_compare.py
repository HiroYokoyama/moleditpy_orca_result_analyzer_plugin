"""
tests/test_mo_compare.py
Coverage for the multi-orbital comparison window (mo_compare.MOCompareDialog)
and the MODialog plumbing it depends on: the silent cube-generation batch, the
launcher/teardown handshake, and the context-menu label.

The two modules are loaded through separate isolated harness loads, so their
``Qt`` stubs are distinct objects and a role constant taken from one does not
match the other. Tests that cross the boundary rebind ``mo_compare.Qt`` to the
MO dialog's copy, which is what a real single Qt binding gives you.
"""

import os
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch

sys.path.insert(0, os.path.dirname(__file__))
import gui_harness  # noqa: E402

M = gui_harness.load_isolated("mo_analysis")
CMP = gui_harness.load_isolated("mo_compare")


def _restricted():
    return {
        str(i): {"id": i, "energy": e, "occ": occ, "spin": "restricted"}
        for i, (e, occ) in enumerate([(-1.0, 2.0), (-0.5, 2.0), (0.2, 0.0), (0.6, 0.0)])
    }


class _FakeVisualizer:
    """Records what the compare dialog asked the 3D view to draw."""

    calls = []
    loadable = True

    def __init__(self, mw):
        self.mw = mw

    def load_file(self, path):
        self.path = path
        return self.loadable

    def show_iso(self, isovalue, **kw):
        _FakeVisualizer.calls.append({"iso": isovalue, "path": self.path, **kw})


class _CompareCase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.addCleanup(self._tmp.cleanup)
        self.out = os.path.join(self.tmp, "job.out")

        # Both modules recompute settings.json from their own directory, so
        # without this the suite writes into the package source tree.
        saved = M.__file__
        M.__file__ = os.path.join(self.tmp, "mo_analysis.py")
        self.addCleanup(lambda: setattr(M, "__file__", saved))
        saved_cmp = CMP.__file__
        CMP.__file__ = os.path.join(self.tmp, "mo_compare.py")
        self.addCleanup(lambda: setattr(CMP, "__file__", saved_cmp))

        self.host = MagicMock()
        self.host.context = MagicMock()
        self.host.parser.filename = self.out
        self.host.file_path = self.out

        self.dlg = M.MODialog(self.host, _restricted(), result_dir=self.tmp)
        self.dlg.parent_dlg = self.host
        self.dlg.settings_file = os.path.join(self.tmp, "settings.json")

        # One Qt enum, as a real binding would have.
        p = patch.object(CMP, "Qt", M.Qt)
        p.start()
        self.addCleanup(p.stop)

        _FakeVisualizer.calls = []
        _FakeVisualizer.loadable = True
        p2 = patch.object(CMP, "CubeVisualizer", _FakeVisualizer)
        p2.start()
        self.addCleanup(p2.stop)

    def _make_compare(self):
        return CMP.MOCompareDialog(self.dlg)

    def _write_cube(self, display_id):
        path = self.dlg.get_cube_path(display_id)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("cube\n")
        return path

    def _row(self, display_id):
        tree = self.dlg.tree
        for i in range(tree.topLevelItemCount()):
            if tree.topLevelItem(i).text(0) == display_id:
                return tree.topLevelItem(i)
        raise AssertionError(f"no row {display_id}")


# ---------------------------------------------------------------------------
# Context menu wording
# ---------------------------------------------------------------------------


class TestContextMenuLabel(_CompareCase):
    def _labels(self):
        captured = []

        class _Menu:
            def __init__(self, *a, **k):
                pass

            def addAction(self, text):
                captured.append(text)
                act = MagicMock()
                act.setEnabled = lambda v: None
                return act

            def exec(self, *a, **k):
                return None

        self._row("1").setSelected(True)
        with patch.object(M, "QMenu", _Menu):
            self.dlg.show_tree_context_menu(MagicMock())
        return captured

    def test_regenerate_entry_no_longer_says_overwrite(self):
        self.assertNotIn("overwrite", " ".join(self._labels()).lower())

    def test_regenerate_entry_still_names_the_action_and_count(self):
        self.assertEqual(self._labels()[1], "Regenerate Cube (1)")

    def test_the_overwrite_meaning_survives_in_the_tooltip(self):
        """The wording moved, the behaviour did not."""
        self.assertIn(
            "replacing it",
            open(
                os.path.join(
                    os.path.dirname(__file__),
                    "..",
                    "orca_result_analyzer",
                    "mo_analysis.py",
                ),
                encoding="utf-8",
            ).read(),
        )


# ---------------------------------------------------------------------------
# Silent batch generation on the MO dialog
# ---------------------------------------------------------------------------


class TestSilentGeneration(_CompareCase):
    def setUp(self):
        super().setUp()
        self.host.parser.data = {
            "mo_coeffs": {"1": {"coeffs": [{"coeff": 0.1}, {"coeff": 0.2}]}},
            "atoms": ["H", "H"],
            "coords": [(0.0, 0.0, 0.0), (0.0, 0.0, 0.74)],
        }
        engine = MagicMock()
        engine.n_basis = 2
        p = patch.object(self.dlg, "get_engine", return_value=engine)
        p.start()
        self.addCleanup(p.stop)

        self.shown = []
        p2 = patch.object(self.dlg, "show_cube", self.shown.append)
        p2.start()
        self.addCleanup(p2.stop)

    def test_a_silent_batch_writes_the_cube_without_redrawing(self):
        self._write_cube("1")
        self.dlg.generate_cubes(["1"])
        self.assertEqual(self.shown, [])

    def test_a_normal_batch_still_redraws(self):
        path = self._write_cube("1")
        self.dlg.generation_queue = ["1"]
        self.dlg.process_generation_queue()
        self.assertEqual(self.shown, [path])

    def test_the_callback_runs_once_the_batch_drains(self):
        self._write_cube("1")
        seen = []
        self.dlg.generate_cubes(["1"], on_done=lambda: seen.append(True))
        self.assertEqual(seen, [True])

    def test_an_empty_key_list_still_calls_back(self):
        seen = []
        self.dlg.generate_cubes([], on_done=lambda: seen.append(True))
        self.assertEqual(seen, [True])

    def test_none_keys_are_dropped(self):
        seen = []
        self.dlg.generate_cubes([None, None], on_done=lambda: seen.append(True))
        self.assertEqual(seen, [True])
        self.assertEqual(self.dlg.generation_queue, [])

    def test_silence_is_lifted_when_the_batch_ends(self):
        """A stuck flag would make the next plain Visualize draw nothing."""
        self._write_cube("1")
        self.dlg.generate_cubes(["1"])
        self.assertFalse(self.dlg.generation_silent)

    def test_the_callback_does_not_fire_twice(self):
        self._write_cube("1")
        seen = []
        self.dlg.generate_cubes(["1"], on_done=lambda: seen.append(True))
        self.dlg.process_generation_queue()
        self.assertEqual(seen, [True])

    def test_a_silent_batch_does_not_hijack_the_main_view(self):
        """last_cube_path drives update_vis_only; a silent batch stealing it
        made the next colour tweak redraw somebody else's orbital."""
        self.dlg.last_cube_path = "/somewhere/else.cube"
        self._write_cube("1")
        self.dlg.generate_cubes(["1"])
        self.assertEqual(self.dlg.last_cube_path, "/somewhere/else.cube")

    def test_a_normal_batch_does_take_over_the_main_view(self):
        path = self._write_cube("1")
        self.dlg.generation_queue = ["1"]
        self.dlg.process_generation_queue()
        self.assertEqual(self.dlg.last_cube_path, path)

    def test_forcing_is_honoured_by_a_silent_batch(self):
        self._write_cube("1")
        with patch.object(M, "CalcWorker") as worker:
            self.dlg.generate_cubes(["1"], force=True)
            worker.assert_called_once()

    def test_display_cube_forwards_when_not_silent(self):
        self.dlg.generation_silent = False
        self.dlg._display_cube("x.cube")
        self.assertEqual(self.shown, ["x.cube"])


# ---------------------------------------------------------------------------
# Launcher / teardown handshake
# ---------------------------------------------------------------------------


class TestCompareLauncher(_CompareCase):
    def test_the_button_opens_the_comparison_window(self):
        with patch.object(M, "MOCompareDialog") as cls:
            self.dlg.show_compare_dialog()
            cls.assert_called_once_with(self.dlg)

    def test_a_second_click_raises_the_existing_window(self):
        with patch.object(M, "MOCompareDialog") as cls:
            self.dlg.show_compare_dialog()
            self.dlg.show_compare_dialog()
            self.assertEqual(cls.call_count, 1)
        self.dlg.compare_dlg.raise_.assert_called_once()

    def test_a_destroyed_window_is_rebuilt_rather_than_raised(self):
        stale = MagicMock()
        stale.raise_.side_effect = RuntimeError("wrapped C/C++ object deleted")
        self.dlg.compare_dlg = stale
        with patch.object(M, "MOCompareDialog") as cls:
            self.dlg.show_compare_dialog()
            cls.assert_called_once()

    def test_closing_the_window_clears_the_reference(self):
        self.dlg.compare_dlg = MagicMock()
        self.dlg.on_compare_closed()
        self.assertIsNone(self.dlg.compare_dlg)

    def test_the_missing_module_is_reported_not_ignored(self):
        with patch.object(M, "MOCompareDialog", None):
            with patch.object(M, "QMessageBox") as mb:
                self.dlg.show_compare_dialog()
                mb.warning.assert_called_once()

    def test_closing_the_mo_dialog_closes_the_comparison(self):
        child = MagicMock()
        self.dlg.compare_dlg = child
        self.dlg.closeEvent(MagicMock())
        child.close.assert_called_once()
        self.assertIsNone(self.dlg.compare_dlg)


# ---------------------------------------------------------------------------
# Dialog construction
# ---------------------------------------------------------------------------


class TestCompareConstruction(_CompareCase):
    def test_there_are_four_slots(self):
        self.assertEqual(len(self._make_compare().slots), CMP.SLOT_COUNT)

    def test_every_slot_lists_every_orbital(self):
        cmp_dlg = self._make_compare()
        for slot in cmp_dlg.slots:
            self.assertEqual(slot.combo_mo.count(), 4)

    def test_orbital_entries_are_labelled_with_their_frontier_tag(self):
        labels = [lbl for lbl, _k, _d in self._make_compare().collect_orbitals()]
        self.assertTrue(any("HOMO" in lbl for lbl in labels))

    def test_orbital_entries_carry_the_lookup_key_and_display_id(self):
        cmp_dlg = self._make_compare()
        key, display_id = cmp_dlg.slots[0].combo_mo.itemData(0)
        # The top row is the highest-energy orbital, as in the MO table.
        self.assertEqual(display_id, "3")
        self.assertEqual(key, self._row("3").data(0, M.Qt.ItemDataRole.UserRole))

    def test_rows_without_a_lookup_key_are_skipped(self):
        self._row("1").setData(0, M.Qt.ItemDataRole.UserRole, None)
        self.assertEqual(len(self._make_compare().collect_orbitals()), 3)

    def test_no_orbitals_without_a_table(self):
        self.dlg.tree = None
        self.assertEqual(self._make_compare().collect_orbitals(), [])

    def test_each_slot_gets_its_own_actor_namespace(self):
        prefixes = [s.prefix for s in self._make_compare().slots]
        self.assertEqual(len(set(prefixes)), CMP.SLOT_COUNT)

    def test_only_the_first_slot_starts_enabled(self):
        slots = self._make_compare().slots
        self.assertTrue(slots[0].is_on())
        self.assertFalse(any(s.is_on() for s in slots[1:]))


class TestNoOrbitals(_CompareCase):
    """An output whose MO table is empty must still open a usable window."""

    def setUp(self):
        super().setUp()
        self.dlg.tree.clear()

    def test_the_window_still_builds(self):
        self.assertEqual(len(self._make_compare().slots), CMP.SLOT_COUNT)

    def test_a_slot_with_nothing_to_show_reports_no_selection(self):
        self.assertEqual(self._make_compare().slots[0].selection(), (None, None))

    def test_update_generates_nothing(self):
        cmp_dlg = self._make_compare()
        with patch.object(self.dlg, "generate_cubes") as gen:
            cmp_dlg.update_view()
            gen.assert_not_called()

    def test_render_draws_nothing(self):
        self._make_compare().render_all()
        self.assertEqual(_FakeVisualizer.calls, [])

    def test_a_table_row_that_vanished_is_skipped(self):
        self.dlg.tree.topLevelItem = lambda i: None
        self.dlg.tree.topLevelItemCount = lambda: 2
        self.assertEqual(self._make_compare().collect_orbitals(), [])


class TestFirstSlotMirrorsTheMODialog(_CompareCase):
    def test_the_positive_lobe_colour_matches(self):
        self.dlg.set_btn_color(self.dlg.btn_color_p, "#123456")
        self.assertEqual(self._make_compare().slots[0].color("p"), "#123456")

    def test_the_negative_lobe_colour_matches(self):
        self.dlg.set_btn_color(self.dlg.btn_color_n, "#654321")
        self.assertEqual(self._make_compare().slots[0].color("n"), "#654321")

    def test_the_other_slots_get_distinct_colours(self):
        cmp_dlg = self._make_compare()
        pos = [s.color("p") for s in cmp_dlg.slots]
        self.assertEqual(len(set(pos)), CMP.SLOT_COUNT)

    def test_the_isovalue_matches(self):
        self.dlg.spin_iso.setValue(0.077)
        self.assertAlmostEqual(self._make_compare().slots[0].spin_iso.value(), 0.077)

    def test_the_opacity_matches(self):
        self.dlg.spin_opacity.setValue(0.3)
        self.assertAlmostEqual(self._make_compare().slots[0].spin_opacity.value(), 0.3)

    def test_the_style_matches(self):
        self.dlg.combo_style.setCurrentText("Wireframe")
        self.assertEqual(
            self._make_compare().slots[0].combo_style.currentText(), "Wireframe"
        )

    def test_the_smooth_shading_flag_matches(self):
        self.dlg.check_smooth.setChecked(False)
        self.assertFalse(self._make_compare().slots[0].check_smooth.isChecked())

    def test_the_selected_orbital_matches(self):
        self.dlg.tree.setCurrentItem(self._row("2"))
        slot = self._make_compare().slots[0]
        self.assertEqual(slot.selection()[1], "2")

    def test_an_unreadable_parent_colour_falls_back_to_the_default(self):
        self.dlg.get_color_hex = MagicMock(side_effect=AttributeError("gone"))
        self.assertEqual(
            self._make_compare().slots[0].color("p"), CMP.DEFAULT_COLORS[0][0]
        )


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


class TestRendering(_CompareCase):
    def _enable(self, cmp_dlg, index, display_id):
        slot = cmp_dlg.slots[index]
        slot.check_on.setChecked(True)
        for i in range(slot.combo_mo.count()):
            if slot.combo_mo.itemData(i)[1] == display_id:
                slot.combo_mo.setCurrentIndex(i)
                break
        return slot

    def test_an_enabled_slot_with_a_cube_is_drawn(self):
        self._write_cube("1")
        cmp_dlg = self._make_compare()
        self._enable(cmp_dlg, 0, "1")
        cmp_dlg.render_all()
        self.assertEqual(len(_FakeVisualizer.calls), 1)

    def test_four_orbitals_are_drawn_at_once(self):
        cmp_dlg = self._make_compare()
        for i, mo in enumerate(["0", "1", "2", "3"]):
            self._write_cube(mo)
            self._enable(cmp_dlg, i, mo)
        cmp_dlg.render_all()
        self.assertEqual(len(_FakeVisualizer.calls), 4)

    def test_each_orbital_gets_its_own_actor_names(self):
        cmp_dlg = self._make_compare()
        for i, mo in enumerate(["0", "1", "2", "3"]):
            self._write_cube(mo)
            self._enable(cmp_dlg, i, mo)
        cmp_dlg.render_all()
        names = [c["name_prefix"] for c in _FakeVisualizer.calls]
        self.assertEqual(len(set(names)), 4)

    def test_each_slot_renders_with_its_own_colours(self):
        cmp_dlg = self._make_compare()
        for i, mo in enumerate(["0", "1"]):
            self._write_cube(mo)
            self._enable(cmp_dlg, i, mo)
        cmp_dlg.slots[1].set_color("p", "#00ff00")
        cmp_dlg.render_all()
        self.assertEqual(_FakeVisualizer.calls[1]["color_p"], "#00ff00")

    def test_per_slot_isovalue_and_opacity_reach_the_renderer(self):
        self._write_cube("1")
        cmp_dlg = self._make_compare()
        slot = self._enable(cmp_dlg, 0, "1")
        slot.spin_iso.setValue(0.055)
        slot.spin_opacity.setValue(0.8)
        cmp_dlg.render_all()
        call = _FakeVisualizer.calls[0]
        self.assertAlmostEqual(call["iso"], 0.055)
        self.assertAlmostEqual(call["opacity"], 0.8)

    def test_the_style_is_lower_cased_for_pyvista(self):
        self._write_cube("1")
        cmp_dlg = self._make_compare()
        slot = self._enable(cmp_dlg, 0, "1")
        slot.combo_style.setCurrentText("Wireframe")
        cmp_dlg.render_all()
        self.assertEqual(_FakeVisualizer.calls[0]["style"], "wireframe")

    def test_a_disabled_slot_is_not_drawn(self):
        self._write_cube("1")
        cmp_dlg = self._make_compare()
        self._enable(cmp_dlg, 0, "1").check_on.setChecked(False)
        cmp_dlg.render_all()
        self.assertEqual(_FakeVisualizer.calls, [])

    def test_a_disabled_slot_has_its_actors_removed(self):
        cmp_dlg = self._make_compare()
        cmp_dlg.slots[1].check_on.setChecked(False)
        cmp_dlg.render_all()
        removed = [c[0][0] for c in self.dlg.mw.plotter.remove_actor.call_args_list]
        self.assertIn("mo_cmp1_p", removed)

    def test_the_single_orbital_actors_are_cleared_first(self):
        """Otherwise the MO dialog's own lobes double-draw over slot 1."""
        self._write_cube("1")
        cmp_dlg = self._make_compare()
        self._enable(cmp_dlg, 0, "1")
        cmp_dlg.render_all()
        removed = [c[0][0] for c in self.dlg.mw.plotter.remove_actor.call_args_list]
        self.assertIn("mo_iso_p", removed)
        self.assertIn("mo_iso_n", removed)

    def test_a_missing_cube_is_reported_not_drawn(self):
        cmp_dlg = self._make_compare()
        self._enable(cmp_dlg, 0, "1")
        cmp_dlg.render_all()
        self.assertEqual(_FakeVisualizer.calls, [])
        self.assertIn("missing", cmp_dlg.lbl_status.text())

    def test_an_unreadable_cube_is_skipped(self):
        self._write_cube("1")
        _FakeVisualizer.loadable = False
        cmp_dlg = self._make_compare()
        self._enable(cmp_dlg, 0, "1")
        cmp_dlg.render_all()
        self.assertEqual(_FakeVisualizer.calls, [])

    def test_the_status_line_counts_what_was_drawn(self):
        self._write_cube("1")
        cmp_dlg = self._make_compare()
        self._enable(cmp_dlg, 0, "1")
        cmp_dlg.render_all()
        self.assertIn("1 orbital", cmp_dlg.lbl_status.text())

    def test_a_missing_visualiser_is_reported(self):
        cmp_dlg = self._make_compare()
        with patch.object(CMP, "CubeVisualizer", None):
            with patch.object(CMP, "QMessageBox") as mb:
                cmp_dlg.render_all()
                mb.warning.assert_called_once()

    def test_rendering_without_a_main_window_is_a_no_op(self):
        cmp_dlg = self._make_compare()
        cmp_dlg.mw = None
        cmp_dlg.render_all()
        self.assertEqual(_FakeVisualizer.calls, [])


# ---------------------------------------------------------------------------
# Update / generation flow
# ---------------------------------------------------------------------------


class TestUpdateView(_CompareCase):
    def _enable(self, cmp_dlg, index, display_id):
        slot = cmp_dlg.slots[index]
        slot.check_on.setChecked(True)
        for i in range(slot.combo_mo.count()):
            if slot.combo_mo.itemData(i)[1] == display_id:
                slot.combo_mo.setCurrentIndex(i)
                break
        return slot

    def test_cached_orbitals_are_drawn_without_recomputing(self):
        self._write_cube("1")
        cmp_dlg = self._make_compare()
        self._enable(cmp_dlg, 0, "1")
        with patch.object(self.dlg, "generate_cubes") as gen:
            cmp_dlg.update_view()
            gen.assert_not_called()
        self.assertEqual(len(_FakeVisualizer.calls), 1)

    def test_a_missing_cube_is_requested_from_the_mo_dialog(self):
        cmp_dlg = self._make_compare()
        self._enable(cmp_dlg, 0, "1")
        with patch.object(self.dlg, "generate_cubes") as gen:
            cmp_dlg.update_view()
            self.assertEqual(gen.call_args[0][0], ["1"])

    def test_the_render_runs_once_generation_finishes(self):
        cmp_dlg = self._make_compare()
        self._enable(cmp_dlg, 0, "1")
        with patch.object(self.dlg, "generate_cubes") as gen:
            cmp_dlg.update_view()
        self.assertEqual(gen.call_args[1]["on_done"], cmp_dlg.render_all)

    def test_only_missing_cubes_are_generated(self):
        self._write_cube("1")
        cmp_dlg = self._make_compare()
        self._enable(cmp_dlg, 0, "1")
        self._enable(cmp_dlg, 1, "2")
        with patch.object(self.dlg, "generate_cubes") as gen:
            cmp_dlg.update_view()
            self.assertEqual(gen.call_args[0][0], ["2"])

    def test_the_same_orbital_in_two_slots_is_generated_once(self):
        cmp_dlg = self._make_compare()
        self._enable(cmp_dlg, 0, "1")
        self._enable(cmp_dlg, 1, "1")
        with patch.object(self.dlg, "generate_cubes") as gen:
            cmp_dlg.update_view()
            self.assertEqual(gen.call_args[0][0], ["1"])

    def test_disabled_slots_do_not_trigger_generation(self):
        cmp_dlg = self._make_compare()
        with patch.object(self.dlg, "generate_cubes") as gen:
            cmp_dlg.slots[0].check_on.setChecked(False)
            cmp_dlg.update_view()
            gen.assert_not_called()

    def test_a_generation_failure_is_reported_to_the_user(self):
        cmp_dlg = self._make_compare()
        self._enable(cmp_dlg, 0, "1")
        with patch.object(self.dlg, "generate_cubes", side_effect=RuntimeError("boom")):
            with patch.object(CMP, "QMessageBox") as mb:
                cmp_dlg.update_view()
                mb.warning.assert_called_once()


# ---------------------------------------------------------------------------
# Colour picking, clearing and teardown
# ---------------------------------------------------------------------------


class TestSlotColorPicking(_CompareCase):
    def _pick(self, cmp_dlg, slot, which, hex_c, valid=True):
        chosen = MagicMock()
        chosen.isValid.return_value = valid
        chosen.name.return_value = hex_c
        dialog = MagicMock()
        dialog.getColor.return_value = chosen
        with gui_harness.qt_available():
            with patch.object(CMP, "QColorDialog", dialog):
                cmp_dlg.pick_slot_color(slot, which)

    def test_a_chosen_colour_is_stored_on_the_slot(self):
        cmp_dlg = self._make_compare()
        self._pick(cmp_dlg, cmp_dlg.slots[2], "p", "#abcdef")
        self.assertEqual(cmp_dlg.slots[2].color("p"), "#abcdef")

    def test_choosing_a_colour_redraws(self):
        self._write_cube("1")
        cmp_dlg = self._make_compare()
        self._pick(cmp_dlg, cmp_dlg.slots[0], "n", "#abcdef")
        self.assertTrue(_FakeVisualizer.calls)

    def test_cancelling_leaves_the_colour_alone(self):
        cmp_dlg = self._make_compare()
        before = cmp_dlg.slots[3].color("p")
        self._pick(cmp_dlg, cmp_dlg.slots[3], "p", "#abcdef", valid=False)
        self.assertEqual(cmp_dlg.slots[3].color("p"), before)

    def test_a_colour_only_affects_its_own_slot(self):
        cmp_dlg = self._make_compare()
        other = cmp_dlg.slots[1].color("p")
        self._pick(cmp_dlg, cmp_dlg.slots[2], "p", "#abcdef")
        self.assertEqual(cmp_dlg.slots[1].color("p"), other)

    def test_a_slot_with_no_stylesheet_reports_a_default(self):
        cmp_dlg = self._make_compare()
        cmp_dlg.slots[0].btn_p.setStyleSheet("")
        self.assertEqual(cmp_dlg.slots[0].color("p"), "#ff0000")


class TestContrastText(unittest.TestCase):
    def test_dark_backgrounds_get_white_text(self):
        self.assertEqual(CMP.contrast_text("#000080"), "white")

    def test_light_backgrounds_get_black_text(self):
        self.assertEqual(CMP.contrast_text("#ffff00"), "black")

    def test_a_malformed_colour_does_not_raise(self):
        self.assertEqual(CMP.contrast_text("nonsense"), "black")

    def test_a_non_hex_triplet_does_not_raise(self):
        self.assertEqual(CMP.contrast_text("#zzzzzz"), "black")


class TestClearAndClose(_CompareCase):
    def test_clear_all_switches_every_slot_off(self):
        cmp_dlg = self._make_compare()
        cmp_dlg.clear_all()
        self.assertFalse(any(s.is_on() for s in cmp_dlg.slots))

    def test_clear_all_removes_every_actor(self):
        cmp_dlg = self._make_compare()
        cmp_dlg.clear_all()
        removed = [c[0][0] for c in self.dlg.mw.plotter.remove_actor.call_args_list]
        for i in range(CMP.SLOT_COUNT):
            self.assertIn(f"mo_cmp{i}_p", removed)
            self.assertIn(f"mo_cmp{i}_n", removed)

    def test_closing_removes_every_actor(self):
        cmp_dlg = self._make_compare()
        cmp_dlg.closeEvent(MagicMock())
        removed = [c[0][0] for c in self.dlg.mw.plotter.remove_actor.call_args_list]
        self.assertIn("mo_cmp3_n", removed)

    def test_closing_tells_the_mo_dialog_to_forget_it(self):
        self.dlg.compare_dlg = object()
        cmp_dlg = self._make_compare()
        cmp_dlg.closeEvent(MagicMock())
        self.assertIsNone(self.dlg.compare_dlg)

    def test_closing_accepts_the_event_rather_than_recursing(self):
        cmp_dlg = self._make_compare()
        event = MagicMock()
        cmp_dlg.closeEvent(event)
        event.accept.assert_called_once()

    def test_escape_routes_through_the_cleanup(self):
        cmp_dlg = self._make_compare()
        with patch.object(cmp_dlg, "close") as close:
            cmp_dlg.reject()
            close.assert_called_once()


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------


class TestSettingsPersistence(_CompareCase):
    def _settings(self):
        import json

        with open(os.path.join(self.tmp, "settings.json"), encoding="utf-8") as fh:
            return json.load(fh)

    def test_closing_writes_the_slot_settings(self):
        cmp_dlg = self._make_compare()
        cmp_dlg.slots[1].set_color("p", "#010203")
        cmp_dlg.closeEvent(MagicMock())
        self.assertEqual(
            self._settings()["mo_compare"]["slots"][1]["color_p"], "#010203"
        )

    def test_every_slot_is_written(self):
        self._make_compare().closeEvent(MagicMock())
        self.assertEqual(len(self._settings()["mo_compare"]["slots"]), CMP.SLOT_COUNT)

    def test_rendering_saves_too(self):
        self._write_cube("1")
        cmp_dlg = self._make_compare()
        cmp_dlg.slots[0].spin_iso.setValue(0.031)
        cmp_dlg.render_all()
        self.assertAlmostEqual(self._settings()["mo_compare"]["slots"][0]["iso"], 0.031)

    def test_saved_colours_come_back_on_reopen(self):
        first = self._make_compare()
        first.slots[2].set_color("n", "#0a0b0c")
        first.closeEvent(MagicMock())
        self.assertEqual(self._make_compare().slots[2].color("n"), "#0a0b0c")

    def test_saved_isovalue_opacity_and_style_come_back(self):
        first = self._make_compare()
        first.slots[1].spin_iso.setValue(0.044)
        first.slots[1].spin_opacity.setValue(0.9)
        first.slots[1].combo_style.setCurrentText("Points")
        first.slots[1].check_smooth.setChecked(False)
        first.closeEvent(MagicMock())

        slot = self._make_compare().slots[1]
        self.assertAlmostEqual(slot.spin_iso.value(), 0.044)
        self.assertAlmostEqual(slot.spin_opacity.value(), 0.9)
        self.assertEqual(slot.combo_style.currentText(), "Points")
        self.assertFalse(slot.check_smooth.isChecked())

    def test_saved_enabled_state_comes_back(self):
        first = self._make_compare()
        first.slots[3].check_on.setChecked(True)
        first.closeEvent(MagicMock())
        self.assertTrue(self._make_compare().slots[3].is_on())

    def test_a_saved_slot_one_colour_wins_over_the_mo_dialog_default(self):
        first = self._make_compare()
        first.slots[0].set_color("p", "#111111")
        first.closeEvent(MagicMock())
        self.dlg.set_btn_color(self.dlg.btn_color_p, "#999999")
        self.assertEqual(self._make_compare().slots[0].color("p"), "#111111")

    def test_the_mo_dialog_settings_are_not_clobbered(self):
        self.dlg.presets["Mine"] = {"iso": 0.05}
        self.dlg.settings_file = os.path.join(self.tmp, "settings.json")
        self.dlg.save_settings()
        self._make_compare().closeEvent(MagicMock())
        self.assertIn("Mine", self._settings()["mo_settings"]["presets"])

    def test_nothing_saved_yet_leaves_the_defaults_alone(self):
        cmp_dlg = self._make_compare()
        self.assertEqual(cmp_dlg.slots[1].color("p"), CMP.DEFAULT_COLORS[1][0])

    def test_a_corrupt_settings_file_is_survivable(self):
        with open(os.path.join(self.tmp, "settings.json"), "w", encoding="utf-8") as fh:
            fh.write("{not json")
        self.assertEqual(len(self._make_compare().slots), CMP.SLOT_COUNT)

    def test_a_settings_file_without_a_compare_section_is_survivable(self):
        with open(os.path.join(self.tmp, "settings.json"), "w", encoding="utf-8") as fh:
            fh.write('{"mo_settings": {}}')
        self.assertEqual(len(self._make_compare().slots), CMP.SLOT_COUNT)

    def test_malformed_slot_entries_are_ignored(self):
        with open(os.path.join(self.tmp, "settings.json"), "w", encoding="utf-8") as fh:
            fh.write('{"mo_compare": {"slots": ["nonsense", {"iso": "big"}]}}')
        cmp_dlg = self._make_compare()
        self.assertAlmostEqual(cmp_dlg.slots[1].spin_iso.value(), 0.02)

    def test_an_unwritable_settings_path_does_not_break_closing(self):
        cmp_dlg = self._make_compare()
        CMP.__file__ = os.path.join(self.tmp, "job.out", "mo_compare.py")
        with open(self.out, "w", encoding="utf-8") as fh:
            fh.write("x")
        event = MagicMock()
        cmp_dlg.closeEvent(event)
        event.accept.assert_called_once()


if __name__ == "__main__":
    unittest.main()
