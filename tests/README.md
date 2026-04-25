# ORCA Result Analyzer — Test Suite

382 tests across 7 files. All run headlessly — no display server, no PyQt6
installation required. Qt and matplotlib Qt backends are stubbed at module level
in files that need them. `parser.py` and `utils.py` are loaded directly without
any stubs (they have no Qt or third-party dependencies).

---

## Running the tests

```bash
# Full suite
python -m pytest tests/ -v

# With coverage
python -m pytest tests/ --cov=orca_result_analyzer --cov-report=term-missing

# Single file
python -m pytest tests/test_parser_samples.py -v

# Single test
python -m pytest tests/ -k "TestParseBasic"
```

---

## Test files

| File | Tests | Area |
|---|---|---|
| `test_parser_samples.py` | 186 | Real ORCA output files — parser regression tests |
| `test_parser.py` | ~55 | OrcaParser unit tests (isolated method snippets) |
| `test_parser_extended.py` | ~51 | Extended parser coverage: thermal, orbital energies, scan, NEB |
| `test_cube_and_mo.py` | 48 | Cube file parsing, MO coefficients, orbital energies |
| `test_traj_analysis.py` | 20 | TrajectoryResultDialog logic methods |
| `test_init.py` | ~15 | Plugin initialization contract |
| `test_utils.py` | ~7 | Pure utility functions |

---

## Test files — detailed

### `test_parser.py` — OrcaParser unit tests

Each class tests one `parse_*` method in isolation using minimal ORCA output
snippets. No stubs — `parser.py` imports only `re` and `logging`.

| Class | Methods tested |
|---|---|
| `TestParseBasic` | SCF energy (last wins), convergence detection, charge, multiplicity, version, Cartesian geometry (Å vs a.u.), `is_scan`/`is_neb` flags |
| `TestParseScfTrace` | Single block, default label, opt-cycle label, two cycles, duplicate label → suffixed |
| `TestParseDipole` | Components, magnitude (explicit / calculated as √(x²+y²+z²)), multiple occurrences (last wins) |
| `TestParseCharges` | Mulliken, Loewdin, Hirshfeld with spin |
| `TestParseFrequencies` | Count, zero/real/imaginary values, `cm-1` unit |
| `TestParseNmr` | Shielding count/values, J-coupling table (unique pairs) |
| `TestParseTddft` | State count, eV/nm energies, sorted by energy, transitions |
| `TestParseTrajectory` | Single/two opt cycles, scan step type, scan_coord, step index |
| `TestParseGradients` | Colon/no-colon formats, NORM exclusion, `all_gradients` list |
| `TestParseXyzContent` | Multi-frame count, energies, TS/CI filtering, CI-NEB not filtered |

---

### `test_parser_extended.py` — Extended parser coverage

| Class | What is tested |
|---|---|
| `TestParseThermal` | Electronic energy, ZPE, enthalpy, Gibbs, entropy, temperature |
| `TestParseOrbitalEnergies` | Restricted RHF (3 MOs), UHF alpha+beta spins, energy_eh/ev, occupation, type; backward-compat `mos` list |
| `TestParseFrequenciesIR` | IR intensities by mode index (T² column); Raman activity by mode index |
| `TestParseChargesMayer` | Mayer QA/VA/BVA/FA; fallback into Mulliken when no separate block |
| `TestParseScanResultsTable` | `scan_steps` from "Actual Energy" summary; coord/energy mapping |
| `TestParseTrajectoryNeb` | NEB PATH SUMMARY → `neb_image` steps, energies, distances |
| `TestParseOptCycleConvergence` | Convergence dict from GEOMETRY CONVERGENCE block; YES flags |
| `TestParseGradientsMultipleBlocks` | Two gradient blocks; default = last; first accessible |

---

### `test_parser_samples.py` — Real ORCA output files (186 tests)

Tests against real `.out` files in `tests/sample_outputs/`. Catches parser
regressions against production data, not just synthetic snippets.

**Sample files used:**

| File | ORCA | Type | Atoms |
|---|---|---|---|
| `benzene-opt.out` / `benzene-opt_5.out` | 6.1.1 / 5.0.4 | Geometry opt (2 cycles) | C₆H₆ (12) |
| `benzene-opt-ene.out` / `benzene-opt-ene_5.out` | 6.1.1 / 5.0.4 | Single-point | C₆H₆ (12) |
| `benzene-opt-nmr.out` / `benzene-opt-nmr_5.out` | 6.1.1 / 5.0.4 | NMR | C₆H₆ (12) |
| `benzene-opt-vex.out` / `benzene-opt-vex_5.out` | 6.1.1 / 5.0.4 | TD-DFT (5 states) | C₆H₆ (12) |
| `acetone-opt.out` / `acetone-opt_5.out` | 6.x / 5.x | Geometry opt (5 cycles) | C₃H₆O (10) |

Key assertions per file type: SCF energy (Eh), convergence, atom count/symbols,
coordinates, dipole (Debye), Mulliken charges, NMR shieldings (ppm), TDDFT
states (eV/nm), `opt_cycle` / `opt_final` step structure.

`TestVersionConsistency` cross-compares ORCA 5 vs 6 outputs for the same
calculation to ensure version-independent parsing.

---

### `test_cube_and_mo.py` — Cube files and MO coefficients (48 tests)

Uses `benzene-opt-ene_MO_21.cube` (40×40×40 grid, def2-SVP basis, C₆H₆).

| Class | What is tested |
|---|---|
| `TestParseCube` | `_parse_cube`: n_atoms=12, dims=(40,40,40), origin in Bohr, step vectors, 64000 data points, +/− phases present |
| `TestParseOrbitalEnergies` | 32 MOs (21 occupied + 11 virtual); HOMO/LUMO indices; gap positive; `energy` alias = `energy_eh`; `mos` length matches |
| `TestParseMoCoeffs` | 114 entries (full def2-SVP basis); all keys end in `_restricted`; HOMO occupied / LUMO virtual; every MO has ≥1 coefficient; correct atom/orbital/value fields |

**Note:** `_parse_cube` file-reading logic is tested without PyVista by
reimplementing it in the test. `_build_grid` (PyVista-dependent) is not
exercised.

---

### `test_traj_analysis.py` — TrajectoryResultDialog logic (20 tests)

Logic methods are called as unbound methods on minimal fake objects so
`TrajectoryResultDialog.__init__` is never invoked.

`QDialog` and `FigureCanvasQTAgg` are stubbed as real inheritable Python
classes (not `MagicMock()` instances) — required because
`TrajectoryResultDialog` and `MplCanvas` inherit from them.

| Class | What is tested |
|---|---|
| `TestComputeScanPoints` | No scan IDs → all steps; groups by `scan_step_id` (last wins); sorted; `opt_final` preferred; `opt_cycle` preferred over other types; single step per group |
| `TestUpdateDisplayValues` | kJ/mol relative; kcal/mol relative; eV relative; absolute kJ/mol; unknown unit → factor 1.0; minimum energy = zero in relative mode |

---

### `test_init.py` — Plugin initialization contract

| Class | What is tested |
|---|---|
| `TestMetadata` | `PLUGIN_NAME`, `PLUGIN_VERSION`, `PLUGIN_AUTHOR`, `PLUGIN_DESCRIPTION` present and non-empty |
| `TestInitialize` | `initialize(ctx)` registers `.out` file opener (priority 100) and drop handler (priority 100) |
| `TestDropHandler` | Non-.out / .log / non-existent / empty / no-ORCA-header files all return `False`; case-insensitive extension |
| `TestInitializeIdempotent` | Fresh context: exactly 1 opener, 1 handler |

---

### `test_utils.py` — Pure utility functions

Covers `get_default_export_path`: CSV/PNG extensions, custom suffix, empty
input, no-directory, directory preservation.

---

## Mocking strategy

| Module | Approach |
|---|---|
| `parser.py` | Loaded directly — no stubs (stdlib only) |
| `utils.py` | Loaded directly — no stubs |
| `vis.py` (`_parse_cube` only) | Cube-reading logic reimplemented in test — no PyVista |
| `traj_analysis.py` | Qt + matplotlib Qt backend stubbed; `QDialog`/`FigureCanvasQTAgg` as real stub classes |
| `__init__.py` | Minimal PyQt6 (`QMessageBox`, `QFileDialog`) stubbed as `MagicMock()` |

---

## Coverage summary

| Module | Coverage |
|---|---|
| `utils.py` | **100%** |
| `parser.py` | **79%** |
| `traj_analysis.py` | ~12% (Qt UI paths excluded) |
| `__init__.py` | ~20% (`open_orca_file()` dialog excluded) |
| `vis.py` | partial (`_parse_cube` only; `_build_grid` needs PyVista) |
| Qt UI modules | 0% (require live display) |

### Known parser gaps (need new sample files)

| Gap | Reason |
|---|---|
| NEB `.xyz` multi-frame | No NEB sample file |
| NBO / FMO charges | No NBO analysis in samples |
| TDDFT short-table format B | Benzene uses arrow format only |
| Scan results table | No scan output in samples |
