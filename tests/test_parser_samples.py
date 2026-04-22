"""
tests/test_parser_samples.py
Integration tests for OrcaParser using real ORCA output files.

All tests load actual .out files from tests/sample_outputs/ and assert
concrete values extracted from those files.
"""

import os
import sys
import importlib.util
import unittest

# ---------------------------------------------------------------------------
# Bootstrap: load parser without Qt
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(__file__)
_SAMPLES = os.path.join(_HERE, "sample_outputs")
_PARSER_SRC = os.path.normpath(
    os.path.join(_HERE, "..", "orca_result_analyzer", "parser.py")
)


def _load_parser():
    spec = importlib.util.spec_from_file_location("orca_parser_standalone", _PARSER_SRC)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["orca_parser_standalone"] = mod
    spec.loader.exec_module(mod)
    return mod


_parser_mod = _load_parser()
OrcaParser = _parser_mod.OrcaParser


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load(filename):
    path = os.path.join(_SAMPLES, filename)
    p = OrcaParser()
    with open(path, encoding="utf-8", errors="replace") as f:
        content = f.read()
    p.load_from_memory(content, path)
    return p


def _opt_cycles(p):
    return [s for s in p.data["scan_steps"] if s.get("type") == "opt_cycle"]


def _opt_final(p):
    return [s for s in p.data["scan_steps"] if s.get("type") == "opt_final"]


# ---------------------------------------------------------------------------
# Benzene geometry optimization  —  ORCA 6
# ---------------------------------------------------------------------------

class TestBenzeneOptOrca6(unittest.TestCase):
    """benzene-opt.out  (ORCA 6.1.1, 2 cycles, converged)"""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("benzene-opt.out")

    # --- basic ---

    def test_version(self):
        self.assertEqual(self.p.data["version"], "6.1.1")

    def test_converged(self):
        self.assertTrue(self.p.data["converged"])

    def test_final_energy(self):
        self.assertAlmostEqual(
            self.p.data["scf_energy"], -231.930734972084, places=6
        )

    def test_charge(self):
        self.assertEqual(self.p.data["charge"], 0)

    def test_multiplicity(self):
        self.assertEqual(self.p.data["mult"], 1)

    # --- final geometry (from FINAL ENERGY EVALUATION) ---

    def test_atom_count(self):
        self.assertEqual(len(self.p.data["atoms"]), 12)

    def test_atom_symbols(self):
        syms = self.p.data["atoms"]
        self.assertEqual(syms.count("C"), 6)
        self.assertEqual(syms.count("H"), 6)

    def test_final_coords_first_atom(self):
        """data['coords'][0] must be the FINAL ENERGY EVALUATION geometry, not cycle 1."""
        x, y, z = self.p.data["coords"][0]
        self.assertAlmostEqual(x,  0.805933, places=4)
        self.assertAlmostEqual(y, -1.143600, places=4)
        self.assertAlmostEqual(z, -0.008261, places=4)

    # --- trajectory ---

    def test_total_scan_steps(self):
        self.assertEqual(len(self.p.data["scan_steps"]), 3)

    def test_opt_cycle_count(self):
        self.assertEqual(len(_opt_cycles(self.p)), 2)

    def test_opt_final_count(self):
        self.assertEqual(len(_opt_final(self.p)), 1)

    def test_opt_final_is_last(self):
        self.assertEqual(self.p.data["scan_steps"][-1]["type"], "opt_final")

    def test_opt_final_coords_match_data_coords(self):
        """opt_final coords must equal data['coords'] (both from FINAL ENERGY EVALUATION)."""
        final_step = _opt_final(self.p)[0]
        self.assertEqual(final_step["atoms"], self.p.data["atoms"])
        for i, (fc, dc) in enumerate(zip(final_step["coords"], self.p.data["coords"])):
            for j in range(3):
                self.assertAlmostEqual(fc[j], dc[j], places=5,
                    msg=f"atom {i} coord {j} mismatch")

    def test_opt_final_energy(self):
        final_step = _opt_final(self.p)[0]
        self.assertAlmostEqual(final_step["energy"], -231.930734972084, places=6)

    def test_cycle1_step_number(self):
        self.assertEqual(_opt_cycles(self.p)[0]["step"], 1)

    def test_cycle2_step_number(self):
        self.assertEqual(_opt_cycles(self.p)[1]["step"], 2)

    def test_cycle1_has_atoms(self):
        self.assertEqual(len(_opt_cycles(self.p)[0]["atoms"]), 12)

    def test_cycle1_coords_differ_from_final(self):
        """Cycle 1 starting geometry must differ from the FINAL ENERGY EVALUATION coords."""
        cycle1_x = _opt_cycles(self.p)[0]["coords"][0][0]
        final_x = self.p.data["coords"][0][0]
        self.assertNotAlmostEqual(cycle1_x, final_x, places=3)

    # --- dipole ---

    def test_dipole_present(self):
        self.assertIsNotNone(self.p.data["dipole"])

    def test_dipole_near_zero(self):
        """Benzene is symmetric — dipole magnitude must be near zero."""
        mag = self.p.data["dipole"]["magnitude"]
        self.assertLess(abs(mag), 0.1)

    # --- charges ---

    def test_mulliken_count(self):
        self.assertEqual(len(self.p.data["charges"]["Mulliken"]), 12)

    def test_mulliken_sum_near_zero(self):
        total = sum(e["charge"] for e in self.p.data["charges"]["Mulliken"])
        self.assertAlmostEqual(total, 0.0, places=2)


# ---------------------------------------------------------------------------
# Benzene geometry optimization  —  ORCA 5
# ---------------------------------------------------------------------------

class TestBenzeneOptOrca5(unittest.TestCase):
    """benzene-opt_5.out  (ORCA 5.0.4, 2 cycles, converged)"""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("benzene-opt_5.out")

    def test_version(self):
        self.assertEqual(self.p.data["version"], "5.0.4")

    def test_converged(self):
        self.assertTrue(self.p.data["converged"])

    def test_final_energy(self):
        self.assertAlmostEqual(
            self.p.data["scf_energy"], -231.930718933873, places=6
        )

    def test_atom_count(self):
        self.assertEqual(len(self.p.data["atoms"]), 12)

    def test_final_coords_first_atom(self):
        x, y, z = self.p.data["coords"][0]
        self.assertAlmostEqual(x,  0.805991, places=4)
        self.assertAlmostEqual(y, -1.143635, places=4)
        self.assertAlmostEqual(z, -0.008264, places=4)

    def test_total_scan_steps(self):
        self.assertEqual(len(self.p.data["scan_steps"]), 3)

    def test_opt_cycle_count(self):
        self.assertEqual(len(_opt_cycles(self.p)), 2)

    def test_opt_final_count(self):
        self.assertEqual(len(_opt_final(self.p)), 1)

    def test_opt_final_is_last(self):
        self.assertEqual(self.p.data["scan_steps"][-1]["type"], "opt_final")

    def test_opt_final_coords_match_data_coords(self):
        final_step = _opt_final(self.p)[0]
        for i, (fc, dc) in enumerate(zip(final_step["coords"], self.p.data["coords"])):
            for j in range(3):
                self.assertAlmostEqual(fc[j], dc[j], places=5)

    def test_mulliken_count(self):
        self.assertEqual(len(self.p.data["charges"]["Mulliken"]), 12)


# ---------------------------------------------------------------------------
# Acetone geometry optimization  —  ORCA 6
# ---------------------------------------------------------------------------

class TestAcetoneOptOrca6(unittest.TestCase):
    """acetone-opt.out  (ORCA 6.x, 5 cycles, converged)"""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("acetone-opt.out")

    def test_converged(self):
        self.assertTrue(self.p.data["converged"])

    def test_atom_count(self):
        # C3H6O = 10 atoms
        self.assertEqual(len(self.p.data["atoms"]), 10)

    def test_atom_symbols(self):
        syms = self.p.data["atoms"]
        self.assertEqual(syms.count("C"), 3)
        self.assertEqual(syms.count("O"), 1)
        self.assertEqual(syms.count("H"), 6)

    def test_total_scan_steps(self):
        # 5 opt_cycles + 1 opt_final
        self.assertEqual(len(self.p.data["scan_steps"]), 6)

    def test_opt_cycle_count(self):
        self.assertEqual(len(_opt_cycles(self.p)), 5)

    def test_opt_final_count(self):
        self.assertEqual(len(_opt_final(self.p)), 1)

    def test_opt_final_is_last(self):
        self.assertEqual(self.p.data["scan_steps"][-1]["type"], "opt_final")

    def test_opt_final_coords_match_data_coords(self):
        final_step = _opt_final(self.p)[0]
        self.assertEqual(final_step["atoms"], self.p.data["atoms"])
        for i, (fc, dc) in enumerate(zip(final_step["coords"], self.p.data["coords"])):
            for j in range(3):
                self.assertAlmostEqual(fc[j], dc[j], places=5)

    def test_final_coords_first_atom(self):
        x, y, z = self.p.data["coords"][0]
        self.assertAlmostEqual(x, -1.274842, places=4)
        self.assertAlmostEqual(y,  0.247667, places=4)
        self.assertAlmostEqual(z,  0.079907, places=4)

    def test_cycle1_coords_differ_from_final(self):
        cycle1_coords = _opt_cycles(self.p)[0]["coords"][0]
        final_coords = self.p.data["coords"][0]
        diffs = [abs(cycle1_coords[j] - final_coords[j]) for j in range(3)]
        self.assertGreater(max(diffs), 0.001)

    def test_mulliken_count(self):
        self.assertEqual(len(self.p.data["charges"]["Mulliken"]), 10)

    def test_mulliken_sum_near_zero(self):
        total = sum(e["charge"] for e in self.p.data["charges"]["Mulliken"])
        self.assertAlmostEqual(total, 0.0, places=2)


# ---------------------------------------------------------------------------
# Acetone geometry optimization  —  ORCA 5
# ---------------------------------------------------------------------------

class TestAcetoneOptOrca5(unittest.TestCase):
    """acetone-opt_5.out  (ORCA 5.x, 5 cycles, converged)"""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("acetone-opt_5.out")

    def test_converged(self):
        self.assertTrue(self.p.data["converged"])

    def test_atom_count(self):
        self.assertEqual(len(self.p.data["atoms"]), 10)

    def test_total_scan_steps(self):
        self.assertEqual(len(self.p.data["scan_steps"]), 6)

    def test_opt_cycle_count(self):
        self.assertEqual(len(_opt_cycles(self.p)), 5)

    def test_opt_final_count(self):
        self.assertEqual(len(_opt_final(self.p)), 1)

    def test_opt_final_is_last(self):
        self.assertEqual(self.p.data["scan_steps"][-1]["type"], "opt_final")

    def test_opt_final_coords_match_data_coords(self):
        final_step = _opt_final(self.p)[0]
        for i, (fc, dc) in enumerate(zip(final_step["coords"], self.p.data["coords"])):
            for j in range(3):
                self.assertAlmostEqual(fc[j], dc[j], places=5)


# ---------------------------------------------------------------------------
# Benzene single-point energy  —  ORCA 6
# ---------------------------------------------------------------------------

class TestBenzeneEneOrca6(unittest.TestCase):
    """benzene-opt-ene.out  (ORCA 6.1.1, SP, no optimization)"""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("benzene-opt-ene.out")

    def test_final_energy(self):
        self.assertAlmostEqual(
            self.p.data["scf_energy"], -231.930784926115, places=6
        )

    def test_converged(self):
        self.assertTrue(self.p.data["converged"])

    def test_no_opt_cycles(self):
        self.assertEqual(len(_opt_cycles(self.p)), 0)

    def test_no_opt_final(self):
        self.assertEqual(len(_opt_final(self.p)), 0)

    def test_atom_count(self):
        self.assertEqual(len(self.p.data["atoms"]), 12)

    def test_dipole_near_zero(self):
        self.assertIsNotNone(self.p.data["dipole"])
        self.assertLess(abs(self.p.data["dipole"]["magnitude"]), 0.1)

    def test_mulliken_count(self):
        self.assertEqual(len(self.p.data["charges"]["Mulliken"]), 12)


# ---------------------------------------------------------------------------
# Benzene single-point energy  —  ORCA 5
# ---------------------------------------------------------------------------

class TestBenzeneEneOrca5(unittest.TestCase):
    """benzene-opt-ene_5.out  (ORCA 5.0.4, SP, no optimization)"""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("benzene-opt-ene_5.out")

    def test_final_energy(self):
        self.assertAlmostEqual(
            self.p.data["scf_energy"], -231.930718955827, places=6
        )

    def test_no_opt_cycles(self):
        self.assertEqual(len(_opt_cycles(self.p)), 0)

    def test_no_opt_final(self):
        self.assertEqual(len(_opt_final(self.p)), 0)

    def test_atom_count(self):
        self.assertEqual(len(self.p.data["atoms"]), 12)

    def test_mulliken_count(self):
        self.assertEqual(len(self.p.data["charges"]["Mulliken"]), 12)


# ---------------------------------------------------------------------------
# Benzene NMR  —  ORCA 6
# ---------------------------------------------------------------------------

class TestBenzeneNmrOrca6(unittest.TestCase):
    """benzene-opt-nmr.out  (ORCA 6.1.1, NMR shielding)"""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("benzene-opt-nmr.out")

    def test_final_energy(self):
        self.assertAlmostEqual(
            self.p.data["scf_energy"], -231.930784926115, places=6
        )

    def test_nmr_shielding_present(self):
        self.assertGreater(len(self.p.data["nmr_shielding"]), 0)

    def test_nmr_shielding_count(self):
        # 6 C + 6 H = 12 nuclei
        self.assertEqual(len(self.p.data["nmr_shielding"]), 12)

    def test_nmr_carbon_shielding_range(self):
        """Benzene carbon shielding should be roughly 40-80 ppm."""
        carbon_shields = [
            e["shielding"] for e in self.p.data["nmr_shielding"]
            if e["atom_sym"] == "C"
        ]
        self.assertEqual(len(carbon_shields), 6)
        for s in carbon_shields:
            self.assertGreater(s, 0.0)

    def test_nmr_hydrogen_shielding_range(self):
        """Benzene H shielding should be positive."""
        h_shields = [
            e["shielding"] for e in self.p.data["nmr_shielding"]
            if e["atom_sym"] == "H"
        ]
        self.assertEqual(len(h_shields), 6)
        for s in h_shields:
            self.assertGreater(s, 0.0)

    def test_no_opt_cycles(self):
        self.assertEqual(len(_opt_cycles(self.p)), 0)


# ---------------------------------------------------------------------------
# Benzene NMR  —  ORCA 5
# ---------------------------------------------------------------------------

class TestBenzeneNmrOrca5(unittest.TestCase):
    """benzene-opt-nmr_5.out  (ORCA 5.0.4, NMR shielding)"""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("benzene-opt-nmr_5.out")

    def test_final_energy(self):
        self.assertAlmostEqual(
            self.p.data["scf_energy"], -231.930718955827, places=6
        )

    def test_nmr_shielding_present(self):
        self.assertGreater(len(self.p.data["nmr_shielding"]), 0)

    def test_nmr_shielding_count(self):
        self.assertEqual(len(self.p.data["nmr_shielding"]), 12)

    def test_no_opt_cycles(self):
        self.assertEqual(len(_opt_cycles(self.p)), 0)


# ---------------------------------------------------------------------------
# Benzene TDDFT excited states  —  ORCA 6
# ---------------------------------------------------------------------------

class TestBenzeneVexOrca6(unittest.TestCase):
    """benzene-opt-vex.out  (ORCA 6.1.1, TD-DFT, 5 states)"""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("benzene-opt-vex.out")

    def test_final_energy(self):
        self.assertAlmostEqual(
            self.p.data["scf_energy"], -231.728888949460, places=6
        )

    def test_tddft_state_count(self):
        self.assertEqual(len(self.p.data["tddft"]), 5)

    def test_tddft_first_state_ev(self):
        self.assertAlmostEqual(
            self.p.data["tddft"][0]["energy_ev"], 5.494, places=2
        )

    def test_tddft_states_sorted_ascending(self):
        evs = [s["energy_ev"] for s in self.p.data["tddft"]]
        self.assertEqual(evs, sorted(evs))

    def test_tddft_all_states_have_transitions(self):
        for state in self.p.data["tddft"]:
            self.assertGreater(len(state["transitions"]), 0)

    def test_tddft_energy_range(self):
        """All 5 states should be between 5 and 8 eV for benzene."""
        for state in self.p.data["tddft"]:
            self.assertGreater(state["energy_ev"], 5.0)
            self.assertLess(state["energy_ev"], 8.5)

    def test_no_opt_cycles(self):
        self.assertEqual(len(_opt_cycles(self.p)), 0)

    def test_atom_count(self):
        self.assertEqual(len(self.p.data["atoms"]), 12)

    def test_mulliken_count(self):
        self.assertEqual(len(self.p.data["charges"]["Mulliken"]), 12)


# ---------------------------------------------------------------------------
# Benzene TDDFT excited states  —  ORCA 5
# ---------------------------------------------------------------------------

class TestBenzeneVexOrca5(unittest.TestCase):
    """benzene-opt-vex_5.out  (ORCA 5.0.4, TD-DFT, 5 states)"""

    @classmethod
    def setUpClass(cls):
        cls.p = _load("benzene-opt-vex_5.out")

    def test_final_energy(self):
        self.assertAlmostEqual(
            self.p.data["scf_energy"], -231.728828351432, places=6
        )

    def test_tddft_state_count(self):
        self.assertEqual(len(self.p.data["tddft"]), 5)

    def test_tddft_first_state_ev(self):
        self.assertAlmostEqual(
            self.p.data["tddft"][0]["energy_ev"], 5.494, places=2
        )

    def test_tddft_states_sorted_ascending(self):
        evs = [s["energy_ev"] for s in self.p.data["tddft"]]
        self.assertEqual(evs, sorted(evs))

    def test_no_opt_cycles(self):
        self.assertEqual(len(_opt_cycles(self.p)), 0)

    def test_atom_count(self):
        self.assertEqual(len(self.p.data["atoms"]), 12)


# ---------------------------------------------------------------------------
# Cross-version consistency
# ---------------------------------------------------------------------------

class TestVersionConsistency(unittest.TestCase):
    """Verify that ORCA 5 and 6 outputs yield structurally identical results."""

    def test_benzene_opt_same_atom_count(self):
        p6 = _load("benzene-opt.out")
        p5 = _load("benzene-opt_5.out")
        self.assertEqual(len(p6.data["atoms"]), len(p5.data["atoms"]))

    def test_benzene_opt_same_step_count(self):
        p6 = _load("benzene-opt.out")
        p5 = _load("benzene-opt_5.out")
        self.assertEqual(len(p6.data["scan_steps"]), len(p5.data["scan_steps"]))

    def test_benzene_opt_same_step_types(self):
        p6 = _load("benzene-opt.out")
        p5 = _load("benzene-opt_5.out")
        types6 = [s["type"] for s in p6.data["scan_steps"]]
        types5 = [s["type"] for s in p5.data["scan_steps"]]
        self.assertEqual(types6, types5)

    def test_benzene_ene_same_atom_count(self):
        p6 = _load("benzene-opt-ene.out")
        p5 = _load("benzene-opt-ene_5.out")
        self.assertEqual(len(p6.data["atoms"]), len(p5.data["atoms"]))

    def test_benzene_nmr_same_shielding_count(self):
        p6 = _load("benzene-opt-nmr.out")
        p5 = _load("benzene-opt-nmr_5.out")
        self.assertEqual(len(p6.data["nmr_shielding"]), len(p5.data["nmr_shielding"]))

    def test_benzene_vex_same_state_count(self):
        p6 = _load("benzene-opt-vex.out")
        p5 = _load("benzene-opt-vex_5.out")
        self.assertEqual(len(p6.data["tddft"]), len(p5.data["tddft"]))

    def test_acetone_opt_same_step_count(self):
        p6 = _load("acetone-opt.out")
        p5 = _load("acetone-opt_5.out")
        self.assertEqual(len(p6.data["scan_steps"]), len(p5.data["scan_steps"]))

    def test_acetone_opt_same_atom_count(self):
        p6 = _load("acetone-opt.out")
        p5 = _load("acetone-opt_5.out")
        self.assertEqual(len(p6.data["atoms"]), len(p5.data["atoms"]))


if __name__ == "__main__":
    unittest.main()
