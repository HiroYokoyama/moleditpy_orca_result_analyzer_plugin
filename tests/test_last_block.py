"""Repeated ORCA sections: the last block wins, and blocks do not mix.

Compound jobs, optimizations and excited-state runs print the same section
several times. Every parser must read the final one, and must not borrow data
from an earlier block when the final one lacks it.
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(__file__))
from _parser_loader import load_standalone_parser  # noqa: E402

OrcaParser = load_standalone_parser("orca_parser_last_block").OrcaParser


def _parse(content, *methods):
    p = OrcaParser()
    p.raw_content = content
    p.lines = content.splitlines()
    p.filename = "test.out"
    p.data["atoms"] = ["O", "H", "H"]
    for m in methods:
        getattr(p, m)()
    return p.data


def _freq_block(values):
    rows = "\n".join(f"   {i}:      {v:.2f} cm**-1" for i, v in enumerate(values))
    return f"""
-----------------------
VIBRATIONAL FREQUENCIES
-----------------------

{rows}
"""


_IR_BLOCK = """
-----------
IR SPECTRUM
-----------

 Mode   freq       eps      Int      T**2         TX        TY        TZ
       cm**-1   L/(mol*cm) km/mol    a.u.
----------------------------------------------------------------------------
  6:   1600.00   0.012345   62.39  0.002367  ( 0.040181 -0.009108  0.002223)
  7:   3700.00   0.012345   11.00  0.002367  ( 0.040181 -0.009108  0.002223)
  8:   3800.00   0.012345   22.00  0.002367  ( 0.040181 -0.009108  0.002223)

"""


class TestFrequenciesTakeTheLastBlock(unittest.TestCase):
    def test_the_last_frequency_block_is_read(self):
        text = _freq_block([0.0] * 6 + [1600, 3700, 3800]) + _freq_block(
            [0.0] * 6 + [1650, 3750, 3850]
        )
        freqs = _parse(text, "parse_frequencies")["frequencies"]
        self.assertEqual([f["freq"] for f in freqs[6:]], [1650, 3750, 3850])

    def test_an_earlier_ir_table_is_not_borrowed(self):
        # Run 1 has IR, run 2 (the one we read) does not.
        text = (
            _freq_block([0.0] * 6 + [1600, 3700, 3800])
            + _IR_BLOCK
            + _freq_block([0.0] * 6 + [1650, 3750, 3850])
        )
        freqs = _parse(text, "parse_frequencies")["frequencies"]
        self.assertEqual([f["ir"] for f in freqs], [0.0] * 9)

    def test_the_ir_table_of_the_last_block_is_used(self):
        text = _freq_block([0.0] * 6 + [1600, 3700, 3800]) + _IR_BLOCK
        freqs = _parse(text, "parse_frequencies")["frequencies"]
        self.assertEqual(freqs[6]["ir"], 62.39)


def _nmr_summary(values):
    rows = "\n".join(f"  {i}       H          {v:.3f}         10.000" for i, v in values)
    return f"""
--------------------------
CHEMICAL SHIELDING SUMMARY (ppm)
--------------------------

  Nucleus  Element    Isotropic     Anisotropy
  -------  -------  ------------   ------------
{rows}

"""


class TestNMRTakesTheLastBlock(unittest.TestCase):
    def test_the_last_shielding_summary_is_read(self):
        text = _nmr_summary([(1, 31.0), (2, 31.1)]) + _nmr_summary(
            [(1, 30.0), (2, 30.1)]
        )
        shield = _parse(text, "parse_nmr")["nmr_shielding"]
        self.assertEqual([s["shielding"] for s in shield], [30.0, 30.1])


def _tddft_run(e1, e2, header="TD-DFT/TDA EXCITED STATES (SINGLETS)"):
    return f"""
{header}

STATE  1:  E=   0.100000 au      {e1:.3f} eV    24000.0 cm**-1
    10a ->  11a  :     0.980000 (c= -0.98994949)

STATE  2:  E=   0.200000 au      {e2:.3f} eV    48000.0 cm**-1
    10a ->  12a  :     0.970000 (c= -0.98488578)

"""


class TestTDDFTTakesTheLastRun(unittest.TestCase):
    def test_only_the_last_run_is_kept(self):
        # Run 1 had three roots, run 2 only two: state 3 must not linger.
        run1 = _tddft_run(3.0, 4.0) + (
            "STATE  3:  E=   0.300000 au      5.000 eV    72000.0 cm**-1\n"
        )
        data = _parse(run1 + _tddft_run(3.5, 4.5), "parse_tddft")
        self.assertEqual([s["state"] for s in data["tddft"]], [1, 2])
        self.assertEqual([s["energy_ev"] for s in data["tddft"]], [3.5, 4.5])

    def test_triplets_do_not_overwrite_singlets(self):
        text = _tddft_run(3.0, 4.0) + _tddft_run(
            1.5, 2.5, header="TD-DFT/TDA EXCITED STATES (TRIPLETS)"
        )
        data = _parse(text, "parse_tddft")
        self.assertEqual([s["energy_ev"] for s in data["tddft"]], [3.0, 4.0])


if __name__ == "__main__":
    unittest.main()
