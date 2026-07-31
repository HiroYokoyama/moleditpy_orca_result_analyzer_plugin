"""
tests/test_mo_engine.py
Unit tests for mo_engine.py: CubeWriter file format and BasisSetEngine
GTO evaluation. Previously this module had no direct tests.

mo_engine imports numpy (real) and PyQt6.QtCore (stubbed: QThread must be an
inheritable class for CalcWorker's class definition to succeed).
"""

import importlib.util
import math
import os
import sys
import types
import unittest
from unittest.mock import MagicMock

import numpy as np

_SRC = os.path.normpath(
    os.path.join(
        os.path.dirname(__file__), "..", "orca_result_analyzer", "mo_engine.py"
    )
)


def _load_mo_engine():
    """Load mo_engine with a temporary PyQt6.QtCore stub.

    The stub is removed from sys.modules afterwards (originals restored) so
    the minimal QtCore — it only has QThread/pyqtSignal — cannot leak into
    other test modules that install their own, richer PyQt6 stubs.
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
        spec = importlib.util.spec_from_file_location("mo_engine_standalone", _SRC)
        mod = importlib.util.module_from_spec(spec)
        sys.modules["mo_engine_standalone"] = mod
        spec.loader.exec_module(mod)
        return mod
    finally:
        for k, v in saved.items():
            if v is None:
                sys.modules.pop(k, None)
            else:
                sys.modules[k] = v


_mod = _load_mo_engine()
CubeWriter = _mod.CubeWriter
BasisSetEngine = _mod.BasisSetEngine

ANG_TO_BOHR = 1.0 / 0.529177249


# ---------------------------------------------------------------------------
# CubeWriter
# ---------------------------------------------------------------------------


class TestCubeWriter(unittest.TestCase):
    def _write(self, tmpdir, data=None, atoms_sym=(8, 1, 1)):
        # Atomic numbers (ints) avoid any rdkit dependency in to_z().
        if data is None:
            data = np.arange(2 * 2 * 2, dtype=float).reshape(2, 2, 2)
        path = os.path.join(tmpdir, "out.cube")
        CubeWriter.write(
            path,
            list(atoms_sym),
            [(0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0)][: len(atoms_sym)],
            origin=(-1.0, -2.0, -3.0),
            vectors=np.diag([0.5, 0.6, 0.7]),
            data=data,
            comment="MO 5",
        )
        with open(path, encoding="utf-8") as f:
            return f.read().splitlines()

    def test_comment_lines(self):
        import tempfile

        with tempfile.TemporaryDirectory() as td:
            lines = self._write(td)
        self.assertIn("MO 5", lines[0])
        self.assertTrue(lines[1])

    def test_atom_count_and_origin(self):
        import tempfile

        with tempfile.TemporaryDirectory() as td:
            lines = self._write(td)
        parts = lines[2].split()
        self.assertEqual(int(parts[0]), 3)
        self.assertAlmostEqual(float(parts[1]), -1.0, places=5)
        self.assertAlmostEqual(float(parts[3]), -3.0, places=5)

    def test_grid_vector_lines(self):
        import tempfile

        with tempfile.TemporaryDirectory() as td:
            lines = self._write(td)
        nx = lines[3].split()
        ny = lines[4].split()
        nz = lines[5].split()
        self.assertEqual([int(nx[0]), int(ny[0]), int(nz[0])], [2, 2, 2])
        self.assertAlmostEqual(float(nx[1]), 0.5, places=5)
        self.assertAlmostEqual(float(ny[2]), 0.6, places=5)
        self.assertAlmostEqual(float(nz[3]), 0.7, places=5)

    def test_atom_lines_z_charge_and_bohr_coords(self):
        import tempfile

        with tempfile.TemporaryDirectory() as td:
            lines = self._write(td)
        first = lines[6].split()
        self.assertEqual(int(first[0]), 8)
        self.assertAlmostEqual(float(first[1]), 8.0, places=5)
        # Second atom at (0, 0, 1.0 Angstrom) -> z in Bohr
        second = lines[7].split()
        self.assertAlmostEqual(float(second[4]), ANG_TO_BOHR, places=5)

    def test_data_six_values_per_line_and_total(self):
        import tempfile

        data = np.arange(2 * 2 * 2, dtype=float).reshape(2, 2, 2)
        with tempfile.TemporaryDirectory() as td:
            lines = self._write(td, data=data)
        data_lines = lines[9:]  # 2 comments + 4 header + 3 atoms
        vals = [float(v) for ln in data_lines for v in ln.split()]
        self.assertEqual(len(vals), 8)
        self.assertEqual(vals, sorted(vals))  # arange order preserved (C-flatten)
        self.assertTrue(all(len(ln.split()) <= 6 for ln in data_lines))

    def test_creates_missing_directory(self):
        import tempfile

        with tempfile.TemporaryDirectory() as td:
            nested = os.path.join(td, "a", "b")
            path = os.path.join(nested, "out.cube")
            CubeWriter.write(
                path,
                [1],
                [(0.0, 0.0, 0.0)],
                origin=(0, 0, 0),
                vectors=np.eye(3),
                data=np.zeros((2, 2, 2)),
            )
            self.assertTrue(os.path.exists(path))


# ---------------------------------------------------------------------------
# BasisSetEngine
# ---------------------------------------------------------------------------


def _s_shell(alpha=1.0, coeff=1.0, center=(0.0, 0.0, 0.0)):
    return {
        "type": 0,
        "center": np.array(center, dtype=float),
        "exps": np.array([alpha]),
        "coeffs": np.array([coeff]),
    }


class TestBasisSetEngine(unittest.TestCase):
    def test_n_basis_counts_components(self):
        shells = [
            _s_shell(),
            {**_s_shell(), "type": 1},  # P -> 3
            {**_s_shell(), "type": 2},  # D (spherical) -> 5
            {**_s_shell(), "type": 3},  # F (spherical) -> 7
            {**_s_shell(), "type": 4},  # G (spherical) -> 9
        ]
        eng = BasisSetEngine(shells)
        self.assertEqual(eng.n_basis, 1 + 3 + 5 + 7 + 9)

    def test_unsupported_shell_type_skipped(self):
        eng = BasisSetEngine([_s_shell(), {**_s_shell(), "type": 99}])
        self.assertEqual(eng.n_basis, 1)
        val = eng.evaluate_mo_on_grid(0, np.zeros((1, 3)), np.array([1.0]))
        self.assertEqual(val.shape, (1,))

    def test_s_orbital_value_at_center(self):
        eng = BasisSetEngine([_s_shell(alpha=1.0, coeff=1.0)])
        val = eng.evaluate_mo_on_grid(0, np.zeros((1, 3)), np.array([1.0]))
        expected = (2.0 / math.pi) ** 0.75  # GTO normalization prefactor at r=0
        self.assertAlmostEqual(val[0], expected, places=8)

    def test_s_orbital_gaussian_decay(self):
        alpha = 0.8
        eng = BasisSetEngine([_s_shell(alpha=alpha)])
        pts = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.5]])
        vals = eng.evaluate_mo_on_grid(0, pts, np.array([1.0]))
        self.assertAlmostEqual(vals[1] / vals[0], math.exp(-alpha * 1.5**2), places=8)

    def test_p_shell_order_pz_px_py(self):
        eng = BasisSetEngine([{**_s_shell(), "type": 1}])
        pt_z = np.array([[0.0, 0.0, 0.7]])
        pt_x = np.array([[0.7, 0.0, 0.0]])
        # Component 0 = pz: nonzero along z, zero along x
        c_pz = np.array([1.0, 0.0, 0.0])
        self.assertGreater(eng.evaluate_mo_on_grid(0, pt_z, c_pz)[0], 0.0)
        self.assertAlmostEqual(eng.evaluate_mo_on_grid(0, pt_x, c_pz)[0], 0.0)
        # Component 1 = px: nonzero along x, zero along z
        c_px = np.array([0.0, 1.0, 0.0])
        self.assertGreater(eng.evaluate_mo_on_grid(0, pt_x, c_px)[0], 0.0)
        self.assertAlmostEqual(eng.evaluate_mo_on_grid(0, pt_z, c_px)[0], 0.0)

    def test_p_orbital_antisymmetric(self):
        eng = BasisSetEngine([{**_s_shell(), "type": 1}])
        pts = np.array([[0.0, 0.0, 0.9], [0.0, 0.0, -0.9]])
        vals = eng.evaluate_mo_on_grid(0, pts, np.array([1.0, 0.0, 0.0]))
        self.assertAlmostEqual(vals[0], -vals[1], places=10)

    def test_shell_centered_off_origin(self):
        center = (1.0, -2.0, 0.5)
        eng = BasisSetEngine([_s_shell(alpha=1.0, center=center)])
        vals = eng.evaluate_mo_on_grid(
            0, np.array([center, (0.0, 0.0, 0.0)]), np.array([1.0])
        )
        self.assertGreater(vals[0], vals[1])
        self.assertAlmostEqual(vals[0], (2.0 / math.pi) ** 0.75, places=8)

    def test_too_few_coeffs_raises(self):
        eng = BasisSetEngine([{**_s_shell(), "type": 1}])  # needs 3
        with self.assertRaises(ValueError):
            eng.evaluate_mo_on_grid(0, np.zeros((1, 3)), np.array([1.0]))

    def test_near_zero_mo_coeff_skipped(self):
        eng = BasisSetEngine([_s_shell()])
        vals = eng.evaluate_mo_on_grid(0, np.zeros((1, 3)), np.array([1e-12]))
        self.assertEqual(vals[0], 0.0)

    def test_normalization_prefactor_s_and_p(self):
        eng = BasisSetEngine([_s_shell()])
        alpha = 1.3
        n_s = eng._normalization_prefactor(alpha, 0, 0, 0)
        self.assertAlmostEqual(n_s, (2 * alpha / math.pi) ** 0.75, places=10)
        # P: extra sqrt(8*alpha * 1!/2!) = sqrt(4*alpha)
        n_p = eng._normalization_prefactor(alpha, 1, 0, 0)
        self.assertAlmostEqual(
            n_p, (2 * alpha / math.pi) ** 0.75 * math.sqrt(4 * alpha), places=10
        )


class TestSphericalShellNormalization(unittest.TestCase):
    """d/f/g components must be unit normalized and have the right shape.

    The hand-typed angular prefactors these replaced were wrong two ways:
    the f-shell fed raw polynomial coefficients to components with different
    Cartesian normalizations (putting the f(z^3) nodal cone at 29 deg instead
    of 39.2), and the g-shell prefactors disagreed between components (g2 came
    out ~2x the norm of g1), so g lobes had the wrong relative sizes.
    """

    @staticmethod
    def _eng(shell_type):
        return BasisSetEngine([{**_s_shell(alpha=1.0, coeff=1.0), "type": shell_type}])

    def _component_norms(self, shell_type):
        """<phi|phi> per component, on a spherical quadrature grid.

        Gauss-Legendre in cos(theta), uniform in phi, and a fine trapezoid in
        r — the Gaussian decays fast enough that a plain radial grid is well
        inside the tolerance asserted below.
        """
        eng = self._eng(shell_type)
        n_r, n_t, n_phi = 400, 80, 80
        r = np.linspace(1e-6, 6.0, n_r)
        ct, wt = np.polynomial.legendre.leggauss(n_t)
        phi = np.linspace(0.0, 2 * np.pi, n_phi, endpoint=False)
        dr = r[1] - r[0]
        dphi = 2 * np.pi / n_phi

        R, CT, PH = np.meshgrid(r, ct, phi, indexing="ij")
        ST = np.sqrt(1.0 - CT**2)
        pts = np.stack(
            [
                (R * ST * np.cos(PH)).ravel(),
                (R * ST * np.sin(PH)).ravel(),
                (R * CT).ravel(),
            ],
            axis=1,
        )
        weights = (R**2 * np.broadcast_to(wt[None, :, None], R.shape)).ravel()
        weights = weights * dr * dphi

        norms = []
        for i in range(eng.n_basis):
            c = np.zeros(eng.n_basis)
            c[i] = 1.0
            v = eng.evaluate_mo_on_grid(0, pts, c)
            norms.append(float(np.sum(v * v * weights)))
        return norms

    def test_d_components_are_unit_normalized(self):
        for n in self._component_norms(2):
            self.assertAlmostEqual(n, 1.0, places=4)

    def test_f_components_are_unit_normalized(self):
        for n in self._component_norms(3):
            self.assertAlmostEqual(n, 1.0, places=4)

    def test_g_components_are_unit_normalized(self):
        for n in self._component_norms(4):
            self.assertAlmostEqual(n, 1.0, places=4)

    def _first_node_deg(self, shell_type, component, lo=0.5, hi=89.5):
        eng = self._eng(shell_type)
        c = np.zeros(eng.n_basis)
        c[component] = 1.0
        thetas = np.linspace(lo, hi, 20000)
        pts = np.stack(
            [
                np.sin(np.radians(thetas)),
                np.zeros_like(thetas),
                np.cos(np.radians(thetas)),
            ],
            axis=1,
        )
        v = eng.evaluate_mo_on_grid(0, pts, c)
        crossings = thetas[:-1][np.sign(v[:-1]) != np.sign(v[1:])]
        return list(crossings)

    def test_f_z3_nodal_cone(self):
        """f(z^3) changes sign at cos^2(theta) = 3/5."""
        exact = math.degrees(math.acos(math.sqrt(0.6)))
        nodes = self._first_node_deg(3, 0)
        self.assertEqual(len(nodes), 1)
        self.assertAlmostEqual(nodes[0], exact, places=1)

    def test_g_z4_nodal_cones(self):
        """g(z^4) has two cones, at cos^2(theta) = (3 +/- 2*sqrt(6/5))/7."""
        exact = sorted(
            math.degrees(math.acos(math.sqrt(x)))
            for x in (
                (3 + 2 * math.sqrt(1.2)) / 7,
                (3 - 2 * math.sqrt(1.2)) / 7,
            )
        )
        nodes = sorted(self._first_node_deg(4, 0))
        self.assertEqual(len(nodes), 2)
        for got, want in zip(nodes, exact):
            self.assertAlmostEqual(got, want, places=1)

    def test_f_xyz_is_sign_alternating_across_octants(self):
        eng = self._eng(3)
        c = np.zeros(eng.n_basis)
        c[4] = 1.0  # f(-2) = xyz
        pts = np.array([[0.7, 0.7, 0.7], [-0.7, 0.7, 0.7]])
        v = eng.evaluate_mo_on_grid(0, pts, c)
        self.assertAlmostEqual(v[0], -v[1], places=10)


if __name__ == "__main__":
    unittest.main()


class TestCalcWorkerGridGuard(unittest.TestCase):
    """v3.9.1: n_points < 2 used to divide span by zero — numpy silently
    produced inf grid vectors and a corrupt cube file was written."""

    def _make_worker(self, n_points, output_path=None):
        worker = _mod.CalcWorker(
            engine=MagicMock(),
            mo_idx=0,
            n_points=n_points,
            margin=3.0,
            atoms_sym=["H"],
            atoms_coords=[[0.0, 0.0, 0.0]],
            mo_coeffs=[1.0],
            # Guard tests never reach the writer; the success-path test
            # passes a real temp path instead.
            output_path=output_path or os.path.join("unused_dir", "x.cube"),
        )
        worker.finished_sig = MagicMock()
        worker.progress_sig = MagicMock()
        return worker

    def test_single_point_grid_rejected_not_silently_written(self):
        worker = self._make_worker(n_points=1)
        worker.run()
        success, msg = worker.finished_sig.emit.call_args[0]
        self.assertFalse(success)
        self.assertIn("at least 2 points", msg)
        worker.engine.evaluate_mo_on_grid.assert_not_called()

    def test_zero_points_also_rejected(self):
        worker = self._make_worker(n_points=0)
        worker.run()
        success, _ = worker.finished_sig.emit.call_args[0]
        self.assertFalse(success)

    def test_two_points_passes_the_guard(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            worker = self._make_worker(
                n_points=2, output_path=os.path.join(tmp, "x.cube")
            )
            worker.engine.evaluate_mo_on_grid.return_value = np.zeros(8)
            worker.run()
        success, _ = worker.finished_sig.emit.call_args[0]
        # Guard passed: either full success (cube written) or a failure from a
        # later stage — but never the grid-resolution message.
        args = worker.finished_sig.emit.call_args[0]
        self.assertNotIn("at least 2 points", str(args[1]))
