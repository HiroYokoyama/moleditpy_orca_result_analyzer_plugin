"""
tests/test_cli.py
Covers orca_result_analyzer/cli.py: ``moleditpy -r orca_result_analyzer parse``,
which reads an output file without a window and prints a summary or JSON.
"""

import contextlib
import io
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(__file__))
import gui_harness

M = gui_harness.load_isolated("cli")

SAMPLES = os.path.join(os.path.dirname(__file__), "sample_outputs")


def _run(*argv):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        try:
            code = M.main(list(argv))
        except SystemExit as exc:
            code = exc.code
    return code, out.getvalue(), err.getvalue()


def _sample(name):
    return os.path.join(SAMPLES, name)


class Summary(unittest.TestCase):
    def test_an_optimisation_with_frequencies(self):
        code, out, err = _run("parse", _sample("acetone-opt.out"))
        self.assertEqual(code, 0, err)
        self.assertIn("Program:     ORCA 6.1.1", out)
        self.assertIn("Status:      Terminated normally", out)
        self.assertIn("Molecule:    C3H6O, 10 atoms, charge 0, multiplicity 1", out)
        self.assertIn("SCF energy:  -192.8958723644 Eh\n", out)
        self.assertIn("Optimised:   converged", out)
        self.assertIn("Frequencies: 24 modes, 0 imaginary", out)
        self.assertIn("Gibbs:       -192.84077533 Eh at 298.15 K", out)
        self.assertIn("Steps:       5 opt_cycle, 1 opt_final", out)

    def test_an_open_shell_reports_spin_contamination(self):
        code, out, _ = _run("parse", _sample("o2-opt_orb.out"))
        self.assertEqual(code, 0)
        self.assertIn("multiplicity 3", out)
        self.assertIn("<S**2>:      2.0070 (ideal 2.0)", out)

    def test_nbo_and_the_sections_present(self):
        code, out, _ = _run("parse", _sample("c33cl-nbo7.out"))
        self.assertEqual(code, 0)
        self.assertIn("Molecule:    C33H67Cl, 101 atoms", out)
        self.assertIn("NBO:         815 orbitals, 474 E(2) entries", out)
        sections = out.rsplit("Sections:", 1)[1]
        self.assertIn("nbo_population", sections)
        self.assertNotIn("tddft", sections)

    def test_a_scan_is_not_called_an_optimisation(self):
        code, out, _ = _run("parse", _sample("ethane-scan.out"))
        self.assertEqual(code, 0)
        self.assertIn("5 scan_step", out)
        self.assertNotIn("Optimised:", out)

    def test_nmr(self):
        code, out, _ = _run("parse", _sample("benzene-opt-nmr.out"))
        self.assertEqual(code, 0)
        self.assertIn("NMR:         12 shieldings", out)


class Json(unittest.TestCase):
    def test_the_whole_data_dict(self):
        code, out, _ = _run("parse", _sample("acetone-opt.out"), "--json")
        self.assertEqual(code, 0)
        data = json.loads(out)
        self.assertEqual(len(data["atoms"]), 10)
        self.assertEqual(data["thermal"]["imaginary_freq_count"], 0)

    def test_selected_sections_only(self):
        code, out, _ = _run(
            "parse",
            _sample("c33cl-nbo7.out"),
            "--json",
            "-s",
            "scf_energy",
            "-s",
            "mult",
        )
        self.assertEqual(code, 0)
        self.assertEqual(sorted(json.loads(out)), ["mult", "scf_energy"])

    def test_an_unknown_section_lists_the_available_ones(self):
        code, out, err = _run(
            "parse", _sample("acetone-opt.out"), "--json", "-s", "nope"
        )
        self.assertEqual(code, 2)
        self.assertEqual(out, "")
        self.assertIn("no section nope", err)
        self.assertIn("scf_energy", err)

    def test_section_without_json_is_refused(self):
        code, _, err = _run("parse", _sample("acetone-opt.out"), "-s", "atoms")
        self.assertEqual(code, 2)
        self.assertIn("--section needs --json", err)

    def test_output_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "data.json")
            code, out, _ = _run(
                "parse", _sample("acetone-opt.out"), "--json", "-o", path
            )
            self.assertEqual((code, out), (0, ""))
            with open(path, encoding="utf-8") as handle:
                self.assertEqual(json.load(handle)["version"], "6.1.1")


def _frames(text):
    """Split multi-frame XYZ text into (comment, rows) pairs."""
    lines, frames = text.splitlines(), []
    while lines:
        count = int(lines[0])
        frames.append((lines[1], [row.split() for row in lines[2 : 2 + count]]))
        lines = lines[2 + count :]
    return frames


class Xyz(unittest.TestCase):
    def test_the_final_structure(self):
        code, out, _ = _run("xyz", _sample("acetone-opt.out"))
        self.assertEqual(code, 0)
        ((comment, rows),) = _frames(out)
        self.assertEqual(comment, "acetone-opt.out  E = -192.8958723644 Eh")
        self.assertEqual(len(rows), 10)
        self.assertEqual(rows[2][0], "O")
        self.assertAlmostEqual(float(rows[2][3]), -1.786414)

    def test_a_scan_writes_one_frame_per_point(self):
        code, out, _ = _run("xyz", _sample("ethane-scan.out"), "--all")
        self.assertEqual(code, 0)
        frames = _frames(out)
        self.assertEqual(len(frames), 5)
        self.assertEqual(
            frames[0][0],
            "ethane-scan.out  point 1  E = -79.7918467962 Eh  coordinate = -60.0",
        )
        self.assertTrue(all(len(rows) == 8 for _, rows in frames))

    def test_an_optimisation_writes_every_cycle(self):
        code, out, _ = _run("xyz", _sample("acetone-opt.out"), "--all")
        self.assertEqual(code, 0)
        comments = [comment for comment, _ in _frames(out)]
        self.assertEqual(len(comments), 6)
        self.assertIn("opt_cycle 1", comments[0])
        self.assertIn("opt_final 6", comments[-1])

    def test_all_without_steps_is_an_error(self):
        code, out, err = _run("xyz", _sample("benzene-opt-nmr.out"), "--all")
        self.assertEqual((code, out), (1, ""))
        self.assertIn("no optimisation, scan or IRC steps", err)

    def test_output_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "final.xyz")
            code, out, _ = _run("xyz", _sample("o2-opt_orb.out"), "-o", path)
            self.assertEqual((code, out), (0, ""))
            with open(path, encoding="utf-8") as handle:
                self.assertEqual(handle.readline(), "2\n")

    def test_a_comment_stays_on_one_line(self):
        text = M.xyz_frame(["H"], [(0.0, 0.0, 0.0)], "a\nb")
        self.assertEqual(text.splitlines()[1], "a b")


class Errors(unittest.TestCase):
    def test_a_missing_file(self):
        code, out, err = _run("parse", os.path.join(SAMPLES, "missing.out"))
        self.assertEqual((code, out), (1, ""))
        self.assertIn("cannot read", err)

    def test_a_file_that_is_not_orca_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "notes.txt")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write("just some text\n")
            code, out, err = _run("parse", path)
        self.assertEqual((code, out), (1, ""))
        self.assertIn("does not look like ORCA output", err)

    def test_no_command_prints_help(self):
        code, out, _ = _run()
        self.assertEqual(code, 2)
        self.assertIn("parse", out)


class Frontier(unittest.TestCase):
    def test_beta_orbitals_are_ignored(self):
        orbitals = [
            {"spin": "alpha", "occupation": 1.0, "energy_eh": -0.3},
            {"spin": "alpha", "occupation": 0.0, "energy_eh": 0.1},
            {"spin": "beta", "occupation": 1.0, "energy_eh": -0.1},
            {"spin": "beta", "occupation": 0.0, "energy_eh": 0.05},
        ]
        self.assertEqual(M._frontier(orbitals), (-0.3, 0.1))

    def test_an_optimisation_without_its_final_step_is_not_converged(self):
        data = {"scan_steps": [{"type": "opt_cycle"}, {"type": "opt_cycle"}]}
        self.assertFalse(M._optimisation_converged(data))
        lines = M.summary_lines("x.out", dict(data, version="6.1.1"))
        self.assertIn("Optimised:   not converged", lines)

    def test_formula_is_hill_ordered(self):
        self.assertEqual(M._formula(["Cl", "H", "C", "H", "O"]), "CH2ClO")
        self.assertEqual(M._formula(["O", "H", "H"]), "H2O")


class EntryPoint(unittest.TestCase):
    def test_main_module_runs_main(self):
        path = os.path.join(os.path.dirname(M.__file__), "__main__.py")
        with open(path, encoding="utf-8") as handle:
            source = handle.read()
        self.assertIn("from .cli import main", source)
        self.assertIn("sys.exit(main())", source)


if __name__ == "__main__":
    unittest.main()
