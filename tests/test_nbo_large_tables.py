"""
tests/test_nbo_large_tables.py
The extended NBO tables of a large molecule are read to their last row, and
the Lewis / non-Lewis summary line survives a partial set of totals.

The tables used to be read through a fixed window of lines, so a molecule
past 54 atoms lost its later NPA rows, past 78 atoms its later electron
configurations, and a bending table past 199 rows its later rows.
"""

import os
import sys
import unittest
from typing import ClassVar

sys.path.insert(0, os.path.dirname(__file__))
import gui_harness

P = gui_harness.load_isolated("parser")
B = gui_harness.load_isolated("bond_analysis")

ATOMS = 120
BENDS = 250


def _population():
    lines = [
        " Summary of Natural Population Analysis:",
        "",
        "                                     Natural Population",
        "             Natural    ---------------------------------------------",
        "  Atom No    Charge        Core      Valence    Rydberg      Total",
        " " + "-" * 68,
    ]
    lines += [
        f"    C{i:3d}   -0.10000      1.99900     4.09000    0.01100     6.10000"
        for i in range(1, ATOMS + 1)
    ]
    return lines + [
        " " + "=" * 68,
        " * Total *  0.00000      1.99996     7.97158    0.02845    10.00000",
        "",
    ]


def _configurations():
    lines = [
        "    Atom No         Natural Electron Configuration",
        " " + "-" * 76,
    ]
    lines += [
        f"      C{i:3d}      [core]2s( 1.00)2p( 3.10)" for i in range(1, ATOMS + 1)
    ]
    return lines + ["", ""]


def _nao():
    lines = [
        " NATURAL POPULATIONS:  Natural atomic orbital occupancies",
        "",
        "  NAO Atom No lang   Type(AO)    Occupancy      Energy",
        " " + "-" * 55,
    ]
    index = 0
    for atom in range(1, ATOMS + 1):
        for _ in range(50):
            index += 1
            lines.append(
                f"{index:5d}    C{atom:3d}  s      Val( 2s)     1.00000      -0.50000"
            )
        lines.append("")
    return lines + ["", " Summary of something else", ""]


def _bending():
    lines = [
        " NHO DIRECTIONALITY AND BOND BENDING (deviation from line of nuclear centers at",
        "",
        (
            "           NBO            Theta   Phi    Theta   Phi    Dev"
            "   Theta   Phi    Dev"
        ),
        " " + "=" * 79,
    ]
    lines += [
        (
            f"  {k:3d}. BD (   1) C{k:3d} - C{k + 1:3d}     90.0   10.0   88.0   12.0"
            "   2.1    92.0  190.0   2.3"
        )
        for k in range(1, BENDS + 1)
    ]
    return lines + ["", ""]


class LargeMoleculeTables(unittest.TestCase):
    def setUp(self):
        self.parser = P.OrcaParser()
        self.parser.lines = _population() + _configurations() + _nao() + _bending()

    def test_every_atom_of_the_population_summary_is_read(self):
        self.parser.parse_nbo_population()
        rows = self.parser.data["nbo_population"]
        self.assertEqual(len(rows), ATOMS)
        self.assertEqual(rows[-1]["index"], ATOMS)

    def test_every_atom_configuration_is_read(self):
        self.parser.parse_nbo_configurations()
        rows = self.parser.data["nbo_configurations"]
        self.assertEqual(len(rows), ATOMS)
        self.assertEqual(rows[-1]["orbitals"], {"2s": 1.0, "2p": 3.1})

    def test_every_nao_row_is_read(self):
        self.parser.parse_nbo_nao()
        rows = self.parser.data["nbo_nao"]
        self.assertEqual(len(rows), ATOMS * 50)
        self.assertEqual(rows[-1]["atom"], ATOMS)

    def test_every_bending_row_is_read(self):
        self.parser.parse_nbo_bending()
        rows = self.parser.data["nbo_bending"]
        self.assertEqual(len(rows), BENDS)
        self.assertEqual(rows[-1]["atoms"], [("C", BENDS), ("C", BENDS + 1)])

    def test_the_population_table_stops_at_its_rule(self):
        self.parser.lines = _population() + _population()
        self.parser.parse_nbo_population()
        self.assertEqual(len(self.parser.data["nbo_population"]), ATOMS)


class TotalsSummary(unittest.TestCase):
    FULL: ClassVar[dict] = {
        "lewis": 9.98992,
        "lewis_pct": 99.899,
        "valence_nonlewis": 0.00027,
        "valence_nonlewis_pct": 0.003,
        "rydberg_nonlewis": 0.00982,
        "rydberg_nonlewis_pct": 0.098,
        "total_nonlewis": 0.01008,
        "total_nonlewis_pct": 0.101,
    }

    def test_full_totals(self):
        self.assertEqual(
            B._nbo_totals_text(self.FULL),
            "NBO Lewis structure: 9.98992 e (99.90%) · non-Lewis: 0.01008 e "
            "(valence 0.00027 + Rydberg 0.00982)",
        )

    def test_a_missing_non_lewis_total_is_summed(self):
        totals = dict(self.FULL)
        del totals["total_nonlewis"], totals["total_nonlewis_pct"]
        self.assertIn("non-Lewis: 0.01009 e", B._nbo_totals_text(totals))

    def test_only_the_lewis_total(self):
        totals = {"lewis": 9.98992, "lewis_pct": 99.899}
        self.assertEqual(
            B._nbo_totals_text(totals), "NBO Lewis structure: 9.98992 e (99.90%)"
        )

    def test_one_part_without_a_total(self):
        totals = {"lewis": 9.98992, "lewis_pct": 99.899, "valence_nonlewis": 0.00027}
        self.assertTrue(
            B._nbo_totals_text(totals).endswith("· non-Lewis: valence 0.00027")
        )

    def test_nothing_without_a_lewis_total(self):
        self.assertEqual(B._nbo_totals_text(None), "")
        self.assertEqual(B._nbo_totals_text({}), "")


# Rows copied from a real NBO 7.0.10 run on C33H67Cl (101 atoms, Cl last):
# from atom 100 on, NBO prints the number fused to the element symbol.
REAL_POPULATION = [
    " Summary of Natural Population Analysis:",
    "",
    "  Atom No    Charge        Core      Valence    Rydberg      Total",
    " " + "-" * 68,
    "    H 99    0.22924      0.00000     0.76783    0.00294     0.77076",
    "    H100    0.22358      0.00000     0.77300    0.00342     0.77642",
    "   Cl101   -0.10391      9.99997     7.08912    0.01483    17.10391",
    " " + "=" * 68,
]
REAL_CONFIGURATIONS = [
    "    Atom No         Natural Electron Configuration",
    " " + "-" * 76,
    "      H100            1s( 0.77)",
    "     Cl101      [core]3s( 1.85)3p( 5.24)3d( 0.01)",
    "",
]
REAL_NAO = [
    " NATURAL POPULATIONS:  Natural atomic orbital occupancies",
    "",
    "  NAO Atom No lang   Type(AO)    Occupancy      Energy",
    " " + "-" * 55,
    " 793    H100  s      Val( 1s)     0.77300",
    " 815   Cl101  dz2    Ryd( 3d)     0.00112",
    "",
]
REAL_BENDING = [
    " NHO DIRECTIONALITY AND BOND BENDING (deviation from line of nuclear centers at",
    "           NBO            Theta   Phi    Theta   Phi    Dev   Theta   Phi    Dev",
    " " + "=" * 79,
    "  40. LP ( 2)Cl101          --     --    137.8  218.8   --      --     --    --",
    " 140. BD ( 1) C 33- H100   50.2   80.6    48.4   80.2   1.8     --     --    --",
    "",
]


class FusedAtomLabels(unittest.TestCase):
    def setUp(self):
        self.parser = P.OrcaParser()

    def _parse(self, lines, method, key):
        self.parser.lines = lines
        getattr(self.parser, method)()
        return self.parser.data[key]

    def test_population(self):
        rows = self._parse(REAL_POPULATION, "parse_nbo_population", "nbo_population")
        self.assertEqual(
            [(r["element"], r["index"]) for r in rows],
            [("H", 99), ("H", 100), ("Cl", 101)],
        )
        self.assertEqual(rows[-1]["total"], 17.10391)

    def test_configurations(self):
        rows = self._parse(
            REAL_CONFIGURATIONS, "parse_nbo_configurations", "nbo_configurations"
        )
        self.assertEqual(
            [(r["element"], r["index"]) for r in rows], [("H", 100), ("Cl", 101)]
        )
        self.assertEqual(rows[1]["orbitals"], {"3s": 1.85, "3p": 5.24, "3d": 0.01})

    def test_nao(self):
        rows = self._parse(REAL_NAO, "parse_nbo_nao", "nbo_nao")
        self.assertEqual(
            [(r["element"], r["atom"]) for r in rows], [("H", 100), ("Cl", 101)]
        )
        self.assertIsNone(rows[0]["energy"])

    def test_bending(self):
        rows = self._parse(REAL_BENDING, "parse_nbo_bending", "nbo_bending")
        self.assertEqual(rows[0]["atoms"], [("Cl", 101)])
        self.assertEqual(rows[1]["atoms"], [("C", 33), ("H", 100)])
        self.assertEqual(rows[1]["h1_dev"], 1.8)


class RealOrcaSample(unittest.TestCase):
    """c33cl-nbo7.out: ORCA 6.1.1 + NBO 7.0.10, C33H67Cl with Cl as atom 101."""

    @classmethod
    def setUpClass(cls):
        path = os.path.join(
            os.path.dirname(__file__), "sample_outputs", "c33cl-nbo7.out"
        )
        with open(path, encoding="utf-8") as handle:
            text = handle.read()
        cls.parser = P.OrcaParser()
        cls.parser.load_from_memory(text, path)
        cls.data = cls.parser.data

    def test_every_table_is_read_to_its_last_row(self):
        self.assertEqual(len(self.data["nbo_population"]), 101)
        self.assertEqual(len(self.data["nbo_configurations"]), 101)
        self.assertEqual(len(self.data["nbo_nao"]), 815)
        self.assertEqual(len(self.data["nbo_orbitals"]), 815)
        self.assertEqual(len(self.data["nbo_perturbation"]), 474)
        self.assertEqual(len(self.data["nbo_bending"]), 8)

    def test_the_nbo_charges_reach_atom_101(self):
        last = self.data["charges"]["NBO"][-1]
        self.assertEqual((last["atom_sym"], last["atom_idx"]), ("Cl", 100))
        self.assertEqual(len(self.data["charges"]["NBO"]), 101)

    def test_a_bond_to_atom_101_has_both_hybrids(self):
        bond = next(o for o in self.data["nbo_orbitals"] if o["index"] == 141)
        self.assertEqual(bond["atom_indices"], [32, 100])
        self.assertEqual(
            [(h["atom_sym"], h["atom_idx"]) for h in bond["hybrids"]],
            [("C", 32), ("Cl", 100)],
        )

    def test_a_chlorine_lone_pair_keeps_its_atom_in_the_bending_table(self):
        self.assertEqual(self.data["nbo_bending"][0]["atoms"], [("Cl", 101)])


if __name__ == "__main__":
    unittest.main()
