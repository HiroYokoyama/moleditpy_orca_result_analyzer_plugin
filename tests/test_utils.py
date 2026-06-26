"""
tests/test_utils.py
Unit tests for orca_result_analyzer/utils.py (pure Python, no stubs required).
"""

import os
import sys
import importlib.util
import unittest

_SRC = os.path.normpath(
    os.path.join(os.path.dirname(__file__), "..", "orca_result_analyzer", "utils.py")
)


def _load_utils():
    spec = importlib.util.spec_from_file_location("orca_utils_mod", _SRC)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["orca_utils_mod"] = mod
    spec.loader.exec_module(mod)
    return mod


_utils = _load_utils()
get_default_export_path = _utils.get_default_export_path


class TestGetDefaultExportPath(unittest.TestCase):
    def test_csv_extension(self):
        result = get_default_export_path(
            "/some/dir/job.out", suffix="_scf_trace", extension=".csv"
        )
        self.assertEqual(result, os.path.join("/some/dir", "job_scf_trace.csv"))

    def test_default_suffix(self):
        result = get_default_export_path("/data/mol.out", extension=".csv")
        self.assertEqual(result, os.path.join("/data", "mol_analyzed.csv"))

    def test_empty_base_path_returns_empty_string(self):
        self.assertEqual(get_default_export_path(""), "")
        self.assertEqual(get_default_export_path(None), "")

    def test_no_directory(self):
        result = get_default_export_path("job.out", suffix="_result", extension=".txt")
        self.assertEqual(result, "job_result.txt")

    def test_no_extension_arg(self):
        result = get_default_export_path("/dir/calc.out", suffix="_data")
        self.assertEqual(result, os.path.join("/dir", "calc_data"))

    def test_preserves_directory(self):
        base = os.path.join("a", "b", "c", "mol.out")
        result = get_default_export_path(base, suffix="_x", extension=".png")
        expected_dir = os.path.join("a", "b", "c")
        self.assertIn(expected_dir, result)
        self.assertTrue(result.endswith("mol_x.png"))


normalize_atom_symbol = _utils.normalize_atom_symbol
determine_bonds_without_dummies = _utils.determine_bonds_without_dummies


class TestNormalizeAtomSymbol(unittest.TestCase):
    """Tests for normalize_atom_symbol — no RDKit required for dummy cases."""

    def _norm(self, raw):
        return normalize_atom_symbol(raw)

    def test_known_dummy_da(self):
        self.assertEqual(self._norm("DA"), "*")

    def test_known_dummy_bq(self):
        self.assertEqual(self._norm("BQ"), "*")

    def test_known_dummy_du(self):
        self.assertEqual(self._norm("DU"), "*")

    def test_known_dummy_asterisk(self):
        self.assertEqual(self._norm("*"), "*")

    def test_colon_suffix_dummy(self):
        # "X:1" should strip to "X" → dummy
        self.assertEqual(self._norm("X:1"), "*")

    def test_colon_suffix_real(self):
        # "C:2" should strip to "C" → Carbon (real element)
        result = self._norm("C:2")
        self.assertEqual(result, "C")

    def test_whitespace_stripped(self):
        self.assertEqual(self._norm("  DA  "), "*")

    def test_real_element_carbon(self):
        self.assertEqual(self._norm("C"), "C")

    def test_real_element_iron(self):
        self.assertEqual(self._norm("Fe"), "Fe")

    def test_real_element_lowercase(self):
        # lowercase input should be capitalised
        self.assertEqual(self._norm("fe"), "Fe")


# ---------------------------------------------------------------------------
# RDKit-dependent tests — skipped when RDKit is not installed
# ---------------------------------------------------------------------------
try:
    from rdkit import Chem
    from rdkit.Geometry import Point3D

    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False


@unittest.skipUnless(_RDKIT_AVAILABLE, "RDKit not installed")
class TestDetermineBondsWithoutDummies(unittest.TestCase):
    """Tests for determine_bonds_without_dummies()."""

    def _make_mol(self, symbols, coords_angstrom):
        """Build an RWMol with a conformer from lists of symbols and (x,y,z) tuples."""
        mol = Chem.RWMol()
        conf = Chem.Conformer(len(symbols))
        for i, (sym, (x, y, z)) in enumerate(zip(symbols, coords_angstrom)):
            mol.AddAtom(Chem.Atom(sym))
            conf.SetAtomPosition(i, Point3D(x, y, z))
        mol.AddConformer(conf)
        return mol

    def test_pure_real_atoms_get_bonds(self):
        """H2O: both O-H bonds should be found."""
        # Simple water geometry (approximate)
        symbols = ["O", "H", "H"]
        coords = [(0.0, 0.0, 0.0), (0.96, 0.0, 0.0), (-0.24, 0.93, 0.0)]
        mol = self._make_mol(symbols, coords)
        determine_bonds_without_dummies(mol, charge=0, bond_orders=False)
        bond_pairs = {(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()}
        # O(0)-H(1) and O(0)-H(2) must exist
        self.assertTrue(
            (0, 1) in bond_pairs or (1, 0) in bond_pairs,
            f"O-H bond missing; got {bond_pairs}",
        )
        self.assertTrue(
            (0, 2) in bond_pairs or (2, 0) in bond_pairs,
            f"O-H bond missing; got {bond_pairs}",
        )

    def test_dummy_atom_has_no_bonds(self):
        """Molecule with a dummy atom: no bond should touch the dummy index."""
        # H2O (indices 0,1,2) + dummy at index 3 far away
        symbols = ["O", "H", "H", "*"]
        coords = [
            (0.0, 0.0, 0.0),
            (0.96, 0.0, 0.0),
            (-0.24, 0.93, 0.0),
            (99.0, 99.0, 99.0),  # far dummy
        ]
        mol = self._make_mol(symbols, coords)
        determine_bonds_without_dummies(mol, charge=0, bond_orders=False)
        dummy_idx = 3
        for bond in mol.GetBonds():
            self.assertNotIn(
                dummy_idx,
                (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
                "Dummy atom must not participate in any bond",
            )

    def test_all_dummy_atoms_no_crash(self):
        """All-dummy molecule must not raise."""
        symbols = ["*", "*"]
        coords = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
        mol = self._make_mol(symbols, coords)
        try:
            determine_bonds_without_dummies(mol, charge=0, bond_orders=False)
        except Exception as exc:
            self.fail(f"determine_bonds_without_dummies raised unexpectedly: {exc}")
        self.assertEqual(mol.GetNumBonds(), 0)

    def test_non_fatal_on_bad_charge(self):
        """Passing an impossible charge must not crash (exception is swallowed)."""
        symbols = ["O", "H", "H"]
        coords = [(0.0, 0.0, 0.0), (0.96, 0.0, 0.0), (-0.24, 0.93, 0.0)]
        mol = self._make_mol(symbols, coords)
        try:
            determine_bonds_without_dummies(mol, charge=999, bond_orders=True)
        except Exception as exc:
            self.fail(f"Should be non-fatal, but raised: {exc}")


if __name__ == "__main__":
    unittest.main()
