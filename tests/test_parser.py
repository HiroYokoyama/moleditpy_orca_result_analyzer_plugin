"""
tests/test_parser.py
Unit tests for OrcaParser (orca_result_analyzer/parser.py).

parser.py imports only `re` and `logging` — no Qt, no third-party deps —
so OrcaParser is loaded directly without any stubs.
"""

import os
import sys
import importlib.util
import unittest

_PARSER_SRC = os.path.normpath(
    os.path.join(os.path.dirname(__file__), "..", "orca_result_analyzer", "parser.py")
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


def _parse(content):
    """Load content via load_from_memory (runs parse_all)."""
    p = OrcaParser()
    p.load_from_memory(content, "test.out")
    return p


def _parse_method(content, method_name):
    """Run a single parse method on content without calling parse_all."""
    p = OrcaParser()
    p.raw_content = content
    p.lines = content.splitlines()
    p.filename = "test.out"
    getattr(p, method_name)()
    return p


# ---------------------------------------------------------------------------
# TestParseBasic
# ---------------------------------------------------------------------------


class TestParseBasic(unittest.TestCase):
    def test_scf_energy(self):
        p = _parse_method(
            "  FINAL SINGLE POINT ENERGY    -76.23456789\n", "parse_basic"
        )
        self.assertAlmostEqual(p.data["scf_energy"], -76.23456789, places=6)

    def test_scf_energy_multiple_occurrences_last_wins(self):
        content = (
            "  FINAL SINGLE POINT ENERGY    -76.10000000\n"
            "  FINAL SINGLE POINT ENERGY    -76.50000000\n"
        )
        p = _parse_method(content, "parse_basic")
        self.assertAlmostEqual(p.data["scf_energy"], -76.50000000, places=6)

    def test_convergence_scf_converged(self):
        p = _parse_method("  SCF CONVERGED after 5 cycles\n", "parse_basic")
        self.assertTrue(p.data["converged"])

    def test_convergence_optimization_converged(self):
        p = _parse_method("  OPTIMIZATION CONVERGED\n", "parse_basic")
        self.assertTrue(p.data["converged"])

    def test_convergence_hurray(self):
        p = _parse_method(
            "  HURRAY  - the optimization has converged!\n", "parse_basic"
        )
        self.assertTrue(p.data["converged"])

    def test_not_converged_by_default(self):
        p = _parse_method("  FINAL SINGLE POINT ENERGY    -76.0\n", "parse_basic")
        self.assertFalse(p.data["converged"])

    def test_charge_zero(self):
        p = _parse_method("  Total Charge                 0\n", "parse_basic")
        self.assertEqual(p.data["charge"], 0)

    def test_charge_nonzero(self):
        p = _parse_method("  Total Charge                 -1\n", "parse_basic")
        self.assertEqual(p.data["charge"], -1)

    def test_multiplicity(self):
        p = _parse_method("  Multiplicity                 3\n", "parse_basic")
        self.assertEqual(p.data["mult"], 3)

    def test_version(self):
        p = _parse_method("   Program Version 5.0.3 -  RELEASE  -\n", "parse_basic")
        self.assertEqual(p.data["version"], "5.0.3")

    def test_geometry_parsing(self):
        content = (
            "CARTESIAN COORDINATES (ANGSTROEM)\n"
            "---------------------------------\n"
            "  O      0.00000    0.00000    0.11779\n"
            "  H      0.75545    0.00000   -0.47116\n"
            "  H     -0.75545    0.00000   -0.47116\n"
            "\n"
        )
        p = _parse_method(content, "parse_basic")
        self.assertEqual(p.data["atoms"], ["O", "H", "H"])
        self.assertAlmostEqual(p.data["coords"][0][2], 0.11779, places=4)
        self.assertAlmostEqual(p.data["coords"][1][0], 0.75545, places=4)
        self.assertAlmostEqual(p.data["coords"][2][0], -0.75545, places=4)

    def test_geometry_au_not_parsed(self):
        """CARTESIAN COORDINATES (A.U.) must be skipped."""
        content = (
            "CARTESIAN COORDINATES (A.U.)\n"
            "----------------------------\n"
            "  O      0.00000    0.00000    0.22257\n"
            "\n"
        )
        p = _parse_method(content, "parse_basic")
        self.assertEqual(p.data["atoms"], [])

    def test_empty_content(self):
        p = _parse_method("", "parse_basic")
        self.assertIsNone(p.data["scf_energy"])
        self.assertFalse(p.data["converged"])
        self.assertEqual(p.data["atoms"], [])

    def test_is_scan_flag(self):
        p = _parse_method("  RELAXED SURFACE SCAN\n", "parse_basic")
        self.assertTrue(p.data.get("is_scan", False))

    def test_is_neb_flag(self):
        p = _parse_method(
            "   NUDGED ELASTIC BAND calculation starting\n", "parse_basic"
        )
        self.assertTrue(p.data.get("is_neb", False))


# ---------------------------------------------------------------------------
# TestParseScfTrace
# ---------------------------------------------------------------------------

# Minimal SCF ITERATIONS block:
#   - "SCF ITERATIONS" triggers the search
#   - a line with "ITER" AND "ENERGY" is the header
#   - optional "---" separator after header is skipped
#   - data rows: int float ...
#   - empty line ends the block

_SCF_SINGLE = """\
SCF ITERATIONS
 ITER      Energy       Delta-E
 ---
    0   -75.81234000     0.00000
    1   -76.14567000     0.33333
    2   -76.23456000     0.08889

"""

_SCF_WITH_OPT_CYCLE = """\
                         **** OPTIMIZATION CYCLE    2 ****
SCF ITERATIONS
 ITER      Energy       Delta-E
 ---
    0   -76.20000000     0.00000
    1   -76.30000000     0.10000

"""

_SCF_TWO_CYCLES = """\
SCF ITERATIONS
 ITER      Energy       Delta-E
 ---
    0   -75.00000000     0.00000
    1   -75.50000000     0.50000

                         **** OPTIMIZATION CYCLE    1 ****
SCF ITERATIONS
 ITER      Energy       Delta-E
 ---
    0   -76.00000000     0.00000
    1   -76.23456000     0.23456

"""


class TestParseScfTrace(unittest.TestCase):
    def test_single_block_count(self):
        p = _parse_method(_SCF_SINGLE, "parse_scf_trace")
        self.assertEqual(len(p.data["scf_traces"]), 1)

    def test_single_block_iterations(self):
        p = _parse_method(_SCF_SINGLE, "parse_scf_trace")
        iters = p.data["scf_traces"][0]["iterations"]
        self.assertEqual(len(iters), 3)
        self.assertEqual(iters[0]["iter"], 0)
        self.assertAlmostEqual(iters[0]["energy"], -75.81234, places=4)
        self.assertEqual(iters[2]["iter"], 2)
        self.assertAlmostEqual(iters[2]["energy"], -76.23456, places=4)

    def test_single_block_default_label(self):
        p = _parse_method(_SCF_SINGLE, "parse_scf_trace")
        self.assertEqual(p.data["scf_traces"][0]["step"], "Initial")

    def test_opt_cycle_label(self):
        p = _parse_method(_SCF_WITH_OPT_CYCLE, "parse_scf_trace")
        self.assertEqual(len(p.data["scf_traces"]), 1)
        self.assertIn("Cycle", p.data["scf_traces"][0]["step"])
        self.assertIn("2", p.data["scf_traces"][0]["step"])

    def test_two_cycles_produce_two_traces(self):
        p = _parse_method(_SCF_TWO_CYCLES, "parse_scf_trace")
        self.assertEqual(len(p.data["scf_traces"]), 2)

    def test_duplicate_label_gets_suffix(self):
        # Two "Initial" blocks → second should be "Initial (2)"
        content = """\
SCF ITERATIONS
 ITER      Energy       Delta-E
 ---
    0   -75.00000000     0.00000

SCF ITERATIONS
 ITER      Energy       Delta-E
 ---
    0   -75.50000000     0.00000

"""
        p = _parse_method(content, "parse_scf_trace")
        labels = [t["step"] for t in p.data["scf_traces"]]
        self.assertEqual(labels[0], "Initial")
        self.assertIn("2", labels[1])

    def test_empty_content_produces_no_traces(self):
        p = _parse_method("", "parse_scf_trace")
        self.assertEqual(p.data["scf_traces"], [])


# ---------------------------------------------------------------------------
# TestParseDipole
# ---------------------------------------------------------------------------

_DIPOLE_WITH_MAG = """\
Total Dipole Moment    :    0.12345   0.00000   1.23456
Magnitude (Debye)      :    1.24074
"""

_DIPOLE_NO_MAG_LINE = """\
Total Dipole Moment    :    0.00000   0.00000   1.85468
  (not a magnitude line)
"""


class TestParseDipole(unittest.TestCase):
    def test_vector_parsed(self):
        p = _parse_method(_DIPOLE_WITH_MAG, "parse_dipole")
        self.assertIsNotNone(p.data["dipole"])
        x, y, z = p.data["dipole"]["vector"]
        self.assertAlmostEqual(x, 0.12345, places=4)
        self.assertAlmostEqual(y, 0.00000, places=4)
        self.assertAlmostEqual(z, 1.23456, places=4)

    def test_magnitude_from_line(self):
        p = _parse_method(_DIPOLE_WITH_MAG, "parse_dipole")
        self.assertAlmostEqual(p.data["dipole"]["magnitude"], 1.24074, places=4)

    def test_magnitude_calculated_if_no_line(self):
        import math

        p = _parse_method(_DIPOLE_NO_MAG_LINE, "parse_dipole")
        expected = math.sqrt(0.0 + 0.0 + 1.85468**2)
        self.assertAlmostEqual(p.data["dipole"]["magnitude"], expected, places=4)

    def test_no_dipole(self):
        p = _parse_method("  FINAL SINGLE POINT ENERGY   -76.0\n", "parse_dipole")
        self.assertIsNone(p.data["dipole"])

    def test_multiple_occurrences_last_wins(self):
        content = (
            "Total Dipole Moment    :    0.00000   0.00000   0.50000\n"
            "Total Dipole Moment    :    0.00000   0.00000   1.23456\n"
        )
        p = _parse_method(content, "parse_dipole")
        self.assertAlmostEqual(p.data["dipole"]["vector"][2], 1.23456, places=4)


# ---------------------------------------------------------------------------
# TestParseCharges
# ---------------------------------------------------------------------------

_MULLIKEN = """\
MULLIKEN ATOMIC CHARGES
-----------------------
   0   O    :   -0.76543
   1   H    :    0.38271
   2   H    :    0.38272
Sum of atomic charges:         0.00000
"""

_LOEWDIN = """\
LOEWDIN ATOMIC CHARGES
-----------------------
   0   C    :   -0.12345
   1   H    :    0.04115
"""

_HIRSHFELD = """\
HIRSHFELD ANALYSIS
------------------
   0   O    -0.32100  0.12345
   1   H     0.16050  0.00000
"""


class TestParseCharges(unittest.TestCase):
    def test_mulliken_count(self):
        p = _parse_method(_MULLIKEN, "parse_charges")
        self.assertEqual(len(p.data["charges"]["Mulliken"]), 3)

    def test_mulliken_values(self):
        p = _parse_method(_MULLIKEN, "parse_charges")
        entry = p.data["charges"]["Mulliken"][0]
        self.assertEqual(entry["atom_sym"], "O")
        self.assertAlmostEqual(entry["charge"], -0.76543, places=4)

    def test_mulliken_last_atom(self):
        p = _parse_method(_MULLIKEN, "parse_charges")
        entry = p.data["charges"]["Mulliken"][2]
        self.assertEqual(entry["atom_sym"], "H")
        self.assertAlmostEqual(entry["charge"], 0.38272, places=4)

    def test_loewdin_parsed(self):
        p = _parse_method(_LOEWDIN, "parse_charges")
        self.assertIn("Loewdin", p.data["charges"])
        self.assertAlmostEqual(
            p.data["charges"]["Loewdin"][0]["charge"], -0.12345, places=4
        )

    def test_hirshfeld_has_spin(self):
        p = _parse_method(_HIRSHFELD, "parse_charges")
        self.assertIn("Hirshfeld", p.data["charges"])
        entry = p.data["charges"]["Hirshfeld"][0]
        self.assertAlmostEqual(entry["charge"], -0.32100, places=4)
        self.assertAlmostEqual(entry["spin"], 0.12345, places=4)

    def test_empty_content_no_charges(self):
        p = _parse_method("", "parse_charges")
        self.assertEqual(p.data["charges"], {})


# ---------------------------------------------------------------------------
# TestParseFrequencies
# ---------------------------------------------------------------------------

_VIBRATIONAL = """\
VIBRATIONAL FREQUENCIES
-----------------------
   0:         0.00 cm**-1
   1:         0.00 cm**-1
   2:         0.00 cm**-1
   3:      1623.45 cm**-1
   4:      3652.84 cm**-1
   5:      3754.93 cm**-1

"""

_IMAGINARY_FREQ = """\
VIBRATIONAL FREQUENCIES
-----------------------
   0:      -234.56 cm**-1
   1:         0.00 cm**-1
   2:      1200.00 cm**-1

"""


class TestParseFrequencies(unittest.TestCase):
    def test_count(self):
        p = _parse_method(_VIBRATIONAL, "parse_frequencies")
        self.assertEqual(len(p.data["frequencies"]), 6)

    def test_zero_modes(self):
        p = _parse_method(_VIBRATIONAL, "parse_frequencies")
        self.assertAlmostEqual(p.data["frequencies"][0]["freq"], 0.0)
        self.assertAlmostEqual(p.data["frequencies"][1]["freq"], 0.0)

    def test_real_frequencies(self):
        p = _parse_method(_VIBRATIONAL, "parse_frequencies")
        self.assertAlmostEqual(p.data["frequencies"][3]["freq"], 1623.45, places=2)
        self.assertAlmostEqual(p.data["frequencies"][5]["freq"], 3754.93, places=2)

    def test_imaginary_frequency_negative(self):
        p = _parse_method(_IMAGINARY_FREQ, "parse_frequencies")
        self.assertEqual(len(p.data["frequencies"]), 3)
        self.assertLess(p.data["frequencies"][0]["freq"], 0.0)

    def test_empty_content(self):
        p = _parse_method("", "parse_frequencies")
        self.assertEqual(p.data["frequencies"], [])

    def test_cm_minus_1_alternative_format(self):
        content = (
            "VIBRATIONAL FREQUENCIES\n"
            "-----------------------\n"
            "   0:      1234.56 cm-1\n"
        )
        p = _parse_method(content, "parse_frequencies")
        self.assertEqual(len(p.data["frequencies"]), 1)
        self.assertAlmostEqual(p.data["frequencies"][0]["freq"], 1234.56, places=2)


# ---------------------------------------------------------------------------
# TestParseNmr
# ---------------------------------------------------------------------------

_NMR_SHIELDING = """\
CHEMICAL SHIELDING SUMMARY (PPM)
---------------------------------
  N    Nucleus    Shielding
  ---
  0      C       123.456
  1      H        31.234
  2      H        32.100

"""

_NMR_WITH_COUPLINGS = """\
CHEMICAL SHIELDING SUMMARY (PPM)
---------------------------------
  N    Nucleus    Shielding
  ---
  0      C       100.000

SUMMARY OF ISOTROPIC COUPLING CONSTANTS (Hz)
----------------
   0 C   1 H
   0     C     0.000   12.345
   1     H    12.345    0.000
"""


class TestParseNmr(unittest.TestCase):
    def test_shielding_count(self):
        p = _parse_method(_NMR_SHIELDING, "parse_nmr")
        self.assertEqual(len(p.data["nmr_shielding"]), 3)

    def test_shielding_values(self):
        p = _parse_method(_NMR_SHIELDING, "parse_nmr")
        self.assertEqual(p.data["nmr_shielding"][0]["atom_sym"], "C")
        self.assertAlmostEqual(
            p.data["nmr_shielding"][0]["shielding"], 123.456, places=2
        )
        self.assertAlmostEqual(
            p.data["nmr_shielding"][1]["shielding"], 31.234, places=2
        )

    def test_empty_content(self):
        p = _parse_method("", "parse_nmr")
        self.assertEqual(p.data["nmr_shielding"], [])
        self.assertEqual(p.data["nmr_couplings"], [])

    def test_couplings_parsed(self):
        p = _parse_method(_NMR_WITH_COUPLINGS, "parse_nmr")
        self.assertEqual(len(p.data["nmr_couplings"]), 1)
        coupling = p.data["nmr_couplings"][0]
        self.assertEqual(coupling["atom_idx1"], 0)
        self.assertEqual(coupling["atom_idx2"], 1)
        self.assertAlmostEqual(coupling["coupling"], 12.345, places=3)


# ---------------------------------------------------------------------------
# TestParseTddft
# ---------------------------------------------------------------------------

_TDDFT_CONTENT = """\
STATE  1:  E=   0.12345 au      3.360 eV    367.7 nm
    41a ->  42a  :     0.987 (c= 0.993)
STATE  2:  E=   0.15678 au      4.267 eV    290.6 nm
    40a ->  42a  :     0.982 (c= 0.991)
"""


class TestParseTddft(unittest.TestCase):
    def test_states_count(self):
        p = _parse_method(_TDDFT_CONTENT, "parse_tddft")
        self.assertEqual(len(p.data["tddft"]), 2)

    def test_state_energies(self):
        p = _parse_method(_TDDFT_CONTENT, "parse_tddft")
        states = p.data["tddft"]
        self.assertAlmostEqual(states[0]["energy_ev"], 3.360, places=2)
        self.assertAlmostEqual(states[0]["energy_nm"], 367.7, places=0)

    def test_states_sorted_by_energy(self):
        p = _parse_method(_TDDFT_CONTENT, "parse_tddft")
        evs = [s["energy_ev"] for s in p.data["tddft"]]
        self.assertEqual(evs, sorted(evs))

    def test_transitions_captured(self):
        p = _parse_method(_TDDFT_CONTENT, "parse_tddft")
        self.assertTrue(len(p.data["tddft"][0]["transitions"]) > 0)

    def test_empty_content(self):
        p = _parse_method("", "parse_tddft")
        self.assertEqual(p.data["tddft"], [])


# ---------------------------------------------------------------------------
# TestParseTrajectory — Optimization Cycles
# ---------------------------------------------------------------------------

_OPT_CYCLE_CONTENT = """\
                         **** OPTIMIZATION CYCLE    1 ****
  FINAL SINGLE POINT ENERGY    -76.23456789
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  O      0.00000    0.00000    0.11779
  H      0.75545    0.00000   -0.47116
  H     -0.75545    0.00000   -0.47116

"""

_TWO_OPT_CYCLES = """\
                         **** OPTIMIZATION CYCLE    1 ****
  FINAL SINGLE POINT ENERGY    -76.10000000
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  O      0.10000    0.00000    0.11779
  H      0.75545    0.00000   -0.47116
  H     -0.75545    0.00000   -0.47116

                         **** OPTIMIZATION CYCLE    2 ****
  FINAL SINGLE POINT ENERGY    -76.23456789
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  O      0.00000    0.00000    0.11779
  H      0.75545    0.00000   -0.47116
  H     -0.75545    0.00000   -0.47116

"""

_SCAN_STEP_CONTENT = """\
RELAXED SURFACE SCAN STEP   3
  Actual scan coordinate      ...   1.500000
  FINAL SINGLE POINT ENERGY    -76.23456789
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  C      0.00000    0.00000    0.00000
  H      1.50000    0.00000    0.00000

"""


class TestParseTrajectory(unittest.TestCase):
    def test_single_opt_cycle(self):
        p = _parse_method(_OPT_CYCLE_CONTENT, "parse_trajectory")
        steps = p.data["scan_steps"]
        self.assertEqual(len(steps), 1)
        self.assertAlmostEqual(steps[0]["energy"], -76.23456789, places=5)
        self.assertEqual(steps[0]["atoms"], ["O", "H", "H"])

    def test_single_opt_cycle_type(self):
        p = _parse_method(_OPT_CYCLE_CONTENT, "parse_trajectory")
        self.assertEqual(p.data["scan_steps"][0]["type"], "opt_cycle")

    def test_two_opt_cycles(self):
        p = _parse_method(_TWO_OPT_CYCLES, "parse_trajectory")
        steps = p.data["scan_steps"]
        self.assertEqual(len(steps), 2)
        self.assertAlmostEqual(steps[0]["energy"], -76.10000000, places=6)
        self.assertAlmostEqual(steps[1]["energy"], -76.23456789, places=5)

    def test_scan_step(self):
        p = _parse_method(_SCAN_STEP_CONTENT, "parse_trajectory")
        steps = p.data["scan_steps"]
        self.assertEqual(len(steps), 1)
        self.assertEqual(steps[0]["type"], "scan_step")
        self.assertEqual(steps[0]["step"], 3)
        self.assertAlmostEqual(steps[0]["scan_coord"], 1.5, places=4)
        self.assertAlmostEqual(steps[0]["energy"], -76.23456789, places=5)

    def test_empty_content(self):
        p = _parse_method("", "parse_trajectory")
        self.assertEqual(p.data["scan_steps"], [])


# ---------------------------------------------------------------------------
# TestParseGradients
# ---------------------------------------------------------------------------

_GRADIENT_COLON = """\
CARTESIAN GRADIENT
------------------
   1   O   :    0.00012   0.00000  -0.00034
   2   H   :    0.00003   0.00000   0.00017
   3   H   :   -0.00003   0.00000   0.00017
-------
"""

_GRADIENT_NO_COLON = """\
CARTESIAN GRADIENT
------------------
   1   O    0.00012   0.00000  -0.00034
   2   H    0.00003   0.00000   0.00017
-------
"""

_GRADIENT_NORM_SKIPPED = """\
CARTESIAN GRADIENT NORM: 0.00123
  Some other content
"""


class TestParseGradients(unittest.TestCase):
    def test_colon_format_count(self):
        p = _parse_method(_GRADIENT_COLON, "parse_gradients")
        self.assertEqual(len(p.data["gradients"]), 3)

    def test_colon_format_values(self):
        p = _parse_method(_GRADIENT_COLON, "parse_gradients")
        g = p.data["gradients"][0]
        self.assertEqual(g["atom_sym"], "O")
        self.assertAlmostEqual(g["vector"][0], 0.00012, places=5)
        self.assertAlmostEqual(g["vector"][2], -0.00034, places=5)

    def test_no_colon_format(self):
        p = _parse_method(_GRADIENT_NO_COLON, "parse_gradients")
        self.assertEqual(len(p.data["gradients"]), 2)

    def test_norm_line_not_parsed(self):
        """Lines with NORM must be ignored."""
        p = _parse_method(_GRADIENT_NORM_SKIPPED, "parse_gradients")
        self.assertEqual(len(p.data["gradients"]), 0)

    def test_all_gradients_list(self):
        p = _parse_method(_GRADIENT_COLON, "parse_gradients")
        self.assertEqual(len(p.data["all_gradients"]), 1)
        self.assertEqual(len(p.data["all_gradients"][0]["grads"]), 3)

    def test_empty_content(self):
        p = _parse_method("", "parse_gradients")
        self.assertEqual(p.data["gradients"], [])


# ---------------------------------------------------------------------------
# TestParseXyzContent
# ---------------------------------------------------------------------------

_XYZ_MULTI = """\
3
Energy: -76.234
O      0.00000    0.00000    0.11779
H      0.75545    0.00000   -0.47116
H     -0.75545    0.00000   -0.47116
3
Energy: -76.345
O      0.00000    0.00000    0.12000
H      0.75000    0.00000   -0.47000
H     -0.75000    0.00000   -0.47000
"""

_XYZ_WITH_TS = """\
3
TS Structure Energy: -76.100
O      0.00000    0.00000    0.20000
H      0.70000    0.00000   -0.40000
H     -0.70000    0.00000   -0.40000
3
Energy: -76.300
O      0.00000    0.00000    0.11779
H      0.75545    0.00000   -0.47116
H     -0.75545    0.00000   -0.47116
"""

_XYZ_WITH_CI = """\
3
Image 3 (CI) Energy: -76.100
O      0.00000    0.00000    0.20000
H      0.70000    0.00000   -0.40000
H     -0.70000    0.00000   -0.40000
3
Energy: -76.500
O      0.00000    0.00000    0.11779
H      0.75545    0.00000   -0.47116
H     -0.75545    0.00000   -0.47116
"""

_XYZ_CI_NEB_NOT_FILTERED = """\
3
CI-NEB calculation Energy: -76.500
O      0.00000    0.00000    0.11779
H      0.75545    0.00000   -0.47116
H     -0.75545    0.00000   -0.47116
"""


class TestParseXyzContent(unittest.TestCase):
    def test_multi_frame_count(self):
        p = OrcaParser()
        steps = p.parse_xyz_content(_XYZ_MULTI)
        self.assertEqual(len(steps), 2)

    def test_multi_frame_energies(self):
        p = OrcaParser()
        steps = p.parse_xyz_content(_XYZ_MULTI)
        self.assertAlmostEqual(steps[0]["energy"], -76.234, places=3)
        self.assertAlmostEqual(steps[1]["energy"], -76.345, places=3)

    def test_multi_frame_atoms(self):
        p = OrcaParser()
        steps = p.parse_xyz_content(_XYZ_MULTI)
        self.assertEqual(steps[0]["atoms"], ["O", "H", "H"])
        self.assertAlmostEqual(steps[0]["coords"][0][2], 0.11779, places=4)

    def test_ts_frame_filtered(self):
        p = OrcaParser()
        steps = p.parse_xyz_content(_XYZ_WITH_TS)
        self.assertEqual(len(steps), 1)
        self.assertAlmostEqual(steps[0]["energy"], -76.300, places=3)

    def test_ci_frame_filtered(self):
        p = OrcaParser()
        steps = p.parse_xyz_content(_XYZ_WITH_CI)
        self.assertEqual(len(steps), 1)

    def test_ci_neb_not_filtered(self):
        """CI-NEB should NOT be excluded — only plain 'CI' labels are."""
        p = OrcaParser()
        steps = p.parse_xyz_content(_XYZ_CI_NEB_NOT_FILTERED)
        self.assertEqual(len(steps), 1)

    def test_empty_content(self):
        p = OrcaParser()
        steps = p.parse_xyz_content("")
        self.assertEqual(steps, [])


# ---------------------------------------------------------------------------
# TestParseSpinContamination
# ---------------------------------------------------------------------------


class TestParseSpinContamination(unittest.TestCase):
    def test_open_shell_s2(self):
        content = (
            "Expectation value of <S**2>     :     2.007028\n"
            "Ideal value S*(S+1) for S=1.0   :     2.000000\n"
        )
        p = _parse_method(content, "parse_spin_contamination")
        self.assertIsNotNone(p.data["spin_s2"])
        self.assertAlmostEqual(p.data["spin_s2"]["actual"], 2.007028, places=5)
        self.assertAlmostEqual(p.data["spin_s2"]["ideal"], 2.0, places=5)
        self.assertAlmostEqual(p.data["spin_s2"]["contamination"], 0.007028, places=5)

    def test_closed_shell_none(self):
        p = _parse_method(
            "FINAL SINGLE POINT ENERGY -76.0\n", "parse_spin_contamination"
        )
        self.assertIsNone(p.data["spin_s2"])

    def test_last_occurrence_wins(self):
        content = (
            "Expectation value of <S**2>     :     0.760000\n"
            "Ideal value S*(S+1) for S=0.5   :     0.750000\n"
            "Expectation value of <S**2>     :     2.007028\n"
            "Ideal value S*(S+1) for S=1.0   :     2.000000\n"
        )
        p = _parse_method(content, "parse_spin_contamination")
        self.assertAlmostEqual(p.data["spin_s2"]["actual"], 2.007028, places=5)


# ---------------------------------------------------------------------------
# TestParseDispersion
# ---------------------------------------------------------------------------


class TestParseDispersion(unittest.TestCase):
    def test_dispersion_value(self):
        p = _parse_method(
            "Dispersion correction           -0.000882087\n", "parse_dispersion"
        )
        self.assertAlmostEqual(p.data["dispersion"], -0.000882087, places=7)

    def test_no_dispersion_none(self):
        p = _parse_method("FINAL SINGLE POINT ENERGY -76.0\n", "parse_dispersion")
        self.assertIsNone(p.data["dispersion"])

    def test_prose_line_not_matched(self):
        # The descriptive "... London dispersion correction" line carries no
        # number and must not be parsed as a value.
        p = _parse_method(
            "   A generally applicable atomic-charge dependent "
            "London dispersion correction\n",
            "parse_dispersion",
        )
        self.assertIsNone(p.data["dispersion"])


# ---------------------------------------------------------------------------
# TestParseMayerBondOrders
# ---------------------------------------------------------------------------

_MAYER_BONDS = """\
  Mayer bond orders larger than 0.100000
B(  0-C ,  1-C ) :   1.3925 B(  0-C ,  5-C ) :   1.3924 B(  0-C ,  6-H ) :   0.9840
B(  1-C ,  2-C ) :   1.3927 B(  1-C ,  7-H ) :   0.9838

"""


class TestParseMayerBondOrders(unittest.TestCase):
    def test_count_multi_per_line(self):
        p = _parse_method(_MAYER_BONDS, "parse_mayer_bond_orders")
        self.assertEqual(len(p.data["mayer_bond_orders"]), 5)

    def test_values_and_symbols(self):
        p = _parse_method(_MAYER_BONDS, "parse_mayer_bond_orders")
        b0 = p.data["mayer_bond_orders"][0]
        self.assertEqual(b0["atom_idx1"], 0)
        self.assertEqual(b0["atom_idx2"], 1)
        self.assertEqual(b0["atom_sym1"], "C")
        self.assertAlmostEqual(b0["order"], 1.3925, places=4)

    def test_indices_ordered(self):
        p = _parse_method(_MAYER_BONDS, "parse_mayer_bond_orders")
        for b in p.data["mayer_bond_orders"]:
            self.assertLess(b["atom_idx1"], b["atom_idx2"])

    def test_blank_line_terminates(self):
        content = _MAYER_BONDS + "B(  9-C , 10-C ) :   1.0000\n"
        p = _parse_method(content, "parse_mayer_bond_orders")
        # the trailing bond is past the blank line and must NOT be included
        self.assertEqual(len(p.data["mayer_bond_orders"]), 5)

    def test_empty_content(self):
        p = _parse_method("", "parse_mayer_bond_orders")
        self.assertEqual(p.data["mayer_bond_orders"], [])


# ---------------------------------------------------------------------------
# TestParseTerminationStatus
# ---------------------------------------------------------------------------


class TestParseTerminationStatus(unittest.TestCase):
    def test_running_by_default(self):
        p = _parse_method("Some running calculation lines\n", "parse_basic")
        self.assertEqual(p.data["termination_status"], "Running")

    def test_terminated_normally(self):
        content = (
            "Some lines...\n"
            "****ORCA TERMINATED NORMALLY****\n"
            "Some trailing blank lines\n"
        )
        p = _parse_method(content, "parse_basic")
        self.assertEqual(p.data["termination_status"], "Terminated normally")

    def test_module_error_termination(self):
        content = "Some lines...\nORCA finished by error termination in ORCA_GSTEP\n"
        p = _parse_method(content, "parse_basic")
        self.assertEqual(p.data["termination_status"], "ERROR")

    def test_input_error_with_file_line(self):
        content = (
            "            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
            "                                 INPUT ERROR\n"
            "            UNRECOGNIZED OR DUPLICATED KEYWORD(S) IN SIMPLE INPUT LINE\n"
            "                       MK  \n"
            "            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
            "[file orca_main/main_input_keywordline.cpp, line 13872]: \n"
        )
        p = _parse_method(content, "parse_basic")
        self.assertEqual(p.data["termination_status"], "ERROR")

    def test_error_bang_with_file_line(self):
        content = (
            "    ----------------------------------------------------------------------------\n"
            "                                   ERROR !!!\n"
            "       The optimization did not converge but reached the maximum \n"
            "    ----------------------------------------------------------------------------\n"
            "\n"
            "[file orca_main/run.cpp, line 12331]: ORCA finished with error return - aborting the run\n"
        )
        p = _parse_method(content, "parse_basic")
        self.assertEqual(p.data["termination_status"], "ERROR")

    def test_error_return_with_file_line_no_highlevel(self):
        content = (
            "Random lines...\n"
            "[file orca_main/run.cpp, line 12331]: ORCA finished with error return - aborting the run\n"
        )
        p = _parse_method(content, "parse_basic")
        self.assertEqual(p.data["termination_status"], "ERROR")


if __name__ == "__main__":
    unittest.main()
