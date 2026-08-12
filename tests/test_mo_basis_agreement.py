"""
tests/test_mo_basis_agreement.py

Regenerates MO 21 of benzene through BasisSetEngine from two ORCA outputs of
the same optimized geometry -- def2-SVP and def2-QZVP -- and checks that the
two orbitals describe the same physical state.

What this covers: the def2-QZVP job carries f and g shells that def2-SVP does
not, so regenerating from it drives the high-angular-momentum branches of
_prepare_definitions/_precompute_shells end to end, from ORCA text through the
parser to grid values, and pins the result against a stored reference. A basis
mismatch, a dropped shell or a coefficient misalignment shows up here.

What this does NOT cover: the *accuracy* of the f/g normalization. MO 21 is a
benzene pi* orbital, ~88% p by coefficient weight, with g carrying 0.27%; a
wrong g prefactor moves the grid by ~4e-5 relative, against a stored-cube text
precision of ~1e-5. There is no headroom to assert on. Deliberately mangling a
g polynomial leaves every assertion here green. The real guards for that are
the unit-norm and nodal-cone tests in test_mo_engine.py, which probe the
angular functions directly; an output whose frontier orbitals carry real f/g
weight would be needed to add a system-level one.

Nothing is written to disk -- the grids are evaluated in memory. The stored
cubes are read only for their grid definition and as regression references.
"""

import importlib.util
import os
import sys
import types
import unittest
from unittest.mock import MagicMock

import numpy as np

_HERE = os.path.dirname(__file__)
_SAMPLES = os.path.join(_HERE, "sample_outputs")
_PKG = os.path.normpath(os.path.join(_HERE, "..", "orca_result_analyzer"))

# The engine works in Bohr; the parser reports shell centers in Angstrom.
BOHR_TO_ANG = 0.529177249

# Both jobs are single points on the same optimized geometry, so the plugin
# writes both cubes on an identical grid -- no interpolation is needed to
# compare them pointwise.
_CASES = {
    "svp": (
        "benzene-opt-ene.out",
        os.path.join("benzene-opt-ene_cubes", "benzene-opt-ene_MO_21.cube"),
    ),
    "qzvp": (
        "benzene-opt-eneQZ.out",
        os.path.join("benzene-opt-eneQZ_cubes", "benzene-opt-eneQZ_MO_21.cube"),
    ),
}

_MO_KEY = "21_restricted"
_MO_IDX = 21


def _load_module(name, filename):
    """Load one plugin module against a minimal QtCore stub.

    mo_engine needs QThread to be a real, subclassable class for CalcWorker's
    class body to execute. The stub is torn down afterwards so it cannot leak
    into test modules that install their own, richer PyQt6 stand-ins.
    """

    class _QThread:
        def __init__(self, *a, **k):
            pass

    core = types.ModuleType("PyQt6.QtCore")
    core.QThread = _QThread
    core.pyqtSignal = lambda *a, **k: MagicMock()
    pyqt6 = types.ModuleType("PyQt6")
    pyqt6.QtCore = core

    saved = {k: sys.modules.get(k) for k in ("PyQt6", "PyQt6.QtCore")}
    sys.modules["PyQt6"] = pyqt6
    sys.modules["PyQt6.QtCore"] = core
    try:
        spec = importlib.util.spec_from_file_location(
            name, os.path.join(_PKG, filename)
        )
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        spec.loader.exec_module(mod)
        return mod
    finally:
        for key, value in saved.items():
            if value is None:
                sys.modules.pop(key, None)
            else:
                sys.modules[key] = value


_parser_mod = _load_module("orca_parser_basis_agreement", "parser.py")
_engine_mod = _load_module("orca_mo_engine_basis_agreement", "mo_engine.py")
OrcaParser = _parser_mod.OrcaParser
BasisSetEngine = _engine_mod.BasisSetEngine


def _read_cube(path):
    """Grid definition and volumetric data of a cube file.

    Deliberately independent of vis._parse_cube: this module is checking the
    numbers, so it should not inherit a bug from the reader under test.
    """
    with open(path, encoding="utf-8") as fh:
        lines = fh.readlines()

    tokens = lines[2].split()
    n_atoms = abs(int(tokens[0]))
    origin = np.array([float(v) for v in tokens[1:4]])

    header = [lines[i].split() for i in (3, 4, 5)]
    dims = tuple(int(h[0]) for h in header)
    vectors = np.array([[float(v) for v in h[1:4]] for h in header])

    values = np.array(
        [float(v) for v in " ".join(lines[6 + n_atoms :]).split()], dtype=float
    )
    return {"origin": origin, "dims": dims, "vectors": vectors, "values": values}


def _grid_points(cube):
    """The cube's sample points as an (N, 3) array in Bohr, in file order.

    Cube files run z fastest, then y, then x -- which is exactly the order
    np.indices produces for the (nx, ny, nz) shape, so the flattened result
    lines up with the file's value sequence.
    """
    idx = np.indices(cube["dims"]).reshape(3, -1).T
    return cube["origin"] + idx @ cube["vectors"]


def _regenerate(out_name, grid):
    """Evaluate MO 21 on `grid`, driving the same path the dialog drives."""
    path = os.path.join(_SAMPLES, out_name)
    parser = OrcaParser()
    with open(path, encoding="utf-8", errors="replace") as fh:
        parser.load_from_memory(fh.read(), path)

    shells = []
    for shell in parser.data.get("basis_set_shells", []):
        center_ang = np.array(shell.get("origin", shell.get("center", [0, 0, 0])))
        shells.append(
            {
                "type": shell.get("l", shell.get("type", 0)),
                "center": center_ang / BOHR_TO_ANG,
                "exps": np.array(shell["exps"]),
                "coeffs": np.array(shell["coeffs"]),
            }
        )

    engine = BasisSetEngine(shells)
    coeffs = np.array(
        [c["coeff"] for c in parser.data["mo_coeffs"][_MO_KEY]["coeffs"]], dtype=float
    )
    return engine, coeffs, engine.evaluate_mo_on_grid(_MO_IDX, grid, coeffs)


def _cosine(a, b):
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))


class TestBasisSetAgreement(unittest.TestCase):
    """MO 21 regenerated from def2-SVP and def2-QZVP outputs."""

    @classmethod
    def setUpClass(cls):
        cls.cubes = {
            tag: _read_cube(os.path.join(_SAMPLES, cube))
            for tag, (_, cube) in _CASES.items()
        }
        cls.grid = _grid_points(cls.cubes["svp"])
        cls.engines, cls.coeffs, cls.values = {}, {}, {}
        for tag, (out_name, _) in _CASES.items():
            engine, coeffs, values = _regenerate(out_name, cls.grid)
            cls.engines[tag], cls.coeffs[tag], cls.values[tag] = engine, coeffs, values

    def test_both_jobs_share_one_grid(self):
        """The pointwise comparisons below are only meaningful if they do."""
        svp, qzvp = self.cubes["svp"], self.cubes["qzvp"]
        self.assertEqual(svp["dims"], qzvp["dims"])
        np.testing.assert_allclose(svp["origin"], qzvp["origin"])
        np.testing.assert_allclose(svp["vectors"], qzvp["vectors"])

    def test_the_two_bases_span_different_angular_momenta(self):
        """Guards the premise: the QZVP path must actually evaluate f and g.

        If a future sample swap quietly dropped the high-l shells, the
        comparison below would still pass while exercising strictly less.
        """
        types_of = {
            tag: {int(s["type"]) for s in engine.shells}
            for tag, engine in self.engines.items()
        }
        self.assertEqual(types_of["svp"], {0, 1, 2})
        self.assertEqual(types_of["qzvp"], {0, 1, 2, 3, 4})

    def test_every_basis_function_gets_a_coefficient(self):
        for tag in _CASES:
            with self.subTest(basis=tag):
                self.assertEqual(len(self.coeffs[tag]), self.engines[tag].n_basis)

    def test_regeneration_reproduces_the_stored_cube(self):
        """Regression anchor, per basis, before the two are compared."""
        for tag in _CASES:
            with self.subTest(basis=tag):
                self.assertAlmostEqual(
                    _cosine(self.values[tag], self.cubes[tag]["values"]), 1.0, places=6
                )

    def test_svp_and_qzvp_describe_the_same_orbital(self):
        """Same physical orbital, two basis sets: the shapes must match.

        Compared on magnitude because an MO's overall phase is arbitrary --
        ORCA hands these two out with opposite signs. What must hold is that
        the orbitals are parallel, not that they point the same way.

        Threshold sits below the measured 0.979: the residual is real basis-set
        difference, not error, so this catches a gross break (wrong orbital,
        shifted coefficients, mis-mapped shells) rather than small numerics.
        """
        cosine = _cosine(self.values["svp"], self.values["qzvp"])
        self.assertGreater(
            abs(cosine),
            0.97,
            f"SVP and QZVP MO {_MO_IDX} disagree (|cos| = {abs(cosine):.4f}); "
            "the two jobs should describe the same orbital",
        )

    def test_the_orbitals_carry_comparable_weight_on_the_grid(self):
        """A prefactor slip can leave shape intact but change the scale."""
        norms = {tag: float(np.linalg.norm(v)) for tag, v in self.values.items()}
        self.assertAlmostEqual(norms["svp"], norms["qzvp"], delta=0.05 * norms["svp"])

    def test_no_cube_is_written_while_regenerating(self):
        """The engine evaluates in memory; only CubeWriter should touch disk."""
        before = set(os.listdir(os.path.join(_SAMPLES, "benzene-opt-eneQZ_cubes")))
        _regenerate(_CASES["qzvp"][0], self.grid[:64])
        after = set(os.listdir(os.path.join(_SAMPLES, "benzene-opt-eneQZ_cubes")))
        self.assertEqual(before, after)


if __name__ == "__main__":
    unittest.main()
