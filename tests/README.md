# ORCA Result Analyzer — Test Suite

480 tests across 11 files. All run headlessly — no display server, no PyQt6
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

# Single test class
python -m pytest tests/ -k "TestParseBasic"

# Real-context integration tests (main app required)
CI_MAIN_APP_SRC=/path/to/python_molecular_editor/moleditpy/src \
  python -m pytest tests/test_plugin_integration.py tests/test_api.py -v
```

---

## Test files

| File | Tests | Area |
|---|---|---|
| `test_parser_samples.py` | 189 | Real ORCA output files — parser regression tests |
| `test_parser.py` | 66 | `OrcaParser` unit tests (isolated method snippets) |
| `test_cube_and_mo.py` | 48 | Cube file parsing, MO coefficients, orbital energies |
| `test_api.py` | 47 | `plugin_api_checker.py` infrastructure + API contract against main app |
| `test_parser_extended.py` | 40 | Extended parser coverage: thermal, orbital energies, scan, NEB |
| `test_init.py` | 22 | Plugin initialization contract |
| `test_traj_analysis.py` | 20 | `TrajectoryResultDialog` logic methods |
| `test_plugin_integration.py` | 17 | `PluginContext` contract — stub mode (always) + real-context mode (optional) |
| `test_freq_analysis.py` | 13 | `FrequencyDialog.populate_list()` and `update_data()` |
| `test_nmr_mo_fixes.py` | 12 | NMR timer safety and MO `apply_preset()` signal blocking |
| `test_utils.py` | 6 | Pure utility functions |

---

## Test files — detailed

### `test_parser.py` — OrcaParser unit tests (66 tests)

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

### `test_parser_extended.py` — Extended parser coverage (40 tests)

| Class | What is tested |
|---|---|
| `TestParseThermal` | Electronic energy, ZPE, enthalpy, Gibbs, entropy, temperature |
| `TestParseOrbitalEnergies` | Restricted RHF (3 MOs), UHF alpha+beta spins, `energy_eh`/`energy_ev`, occupation, type; backward-compat `mos` list |
| `TestParseFrequenciesIR` | IR intensities by mode index (T² column); Raman activity by mode index |
| `TestParseChargesMayer` | Mayer QA/VA/BVA/FA; fallback into Mulliken when no separate block |
| `TestParseScanResultsTable` | `scan_steps` from "Actual Energy" summary; coord/energy mapping |
| `TestParseTrajectoryNeb` | NEB PATH SUMMARY → `neb_image` steps, energies, distances |
| `TestParseOptCycleConvergence` | Convergence dict from GEOMETRY CONVERGENCE block; YES flags |
| `TestParseGradientsMultipleBlocks` | Two gradient blocks; default = last; first accessible |

---

### `test_parser_samples.py` — Real ORCA output files (189 tests)

Tests against real `.out` files in `tests/sample_outputs/`. Catches parser
regressions against production data, not just synthetic snippets. Covers both
ORCA 5 and ORCA 6 output for every sample type.

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

Uses `benzene-opt-ene_cubes/benzene-opt-ene_MO_21.cube` (40×40×40 grid, def2-SVP basis, C₆H₆).

| Class | What is tested |
|---|---|
| `TestParseCube` | `_parse_cube`: n_atoms=12, dims=(40,40,40), origin in Bohr, step vectors, 64 000 data points, +/− phases present |
| `TestParseOrbitalEnergies` | 32 MOs (21 occupied + 11 virtual); HOMO/LUMO indices; gap positive; `energy` alias = `energy_eh`; `mos` length matches |
| `TestParseMoCoeffs` | 114 entries (full def2-SVP basis); all keys end in `_restricted`; HOMO occupied / LUMO virtual; every MO has ≥1 coefficient; correct atom/orbital/value fields |

`_parse_cube` file-reading logic is tested without PyVista by reimplementing it
in the test. `_build_grid` (PyVista-dependent) is not exercised.

---

### `test_traj_analysis.py` — TrajectoryResultDialog logic (20 tests)

Logic methods called as unbound methods on minimal fake objects so
`TrajectoryResultDialog.__init__` is never invoked.

`QDialog` and `FigureCanvasQTAgg` are stubbed as real inheritable Python classes
(not `MagicMock()` instances) — required because `TrajectoryResultDialog` and
`MplCanvas` inherit from them.

| Class | What is tested |
|---|---|
| `TestComputeScanPoints` | No scan IDs → all steps; groups by `scan_step_id` (last wins); sorted; `opt_final` preferred; `opt_cycle` preferred over other types; single step per group |
| `TestUpdateDisplayValues` | kJ/mol relative; kcal/mol relative; eV relative; absolute kJ/mol; unknown unit → factor 1.0; minimum energy = zero in relative mode |

---

### `test_init.py` — Plugin initialization contract (22 tests)

`__init__.py` imports `QMessageBox` and `QFileDialog` at module level so minimal
stubs are installed before loading. Qt-dependent paths inside `open_orca_file()`
are not exercised here.

| Class | What is tested |
|---|---|
| `TestMetadata` | `PLUGIN_NAME`, `PLUGIN_VERSION`, `PLUGIN_AUTHOR`, `PLUGIN_DESCRIPTION` — present and non-empty |
| `TestInitialize` | `initialize(ctx)` registers `.out` file opener (priority 100) and drop handler (priority 100) |
| `TestDropHandler` | Non-.out / .log / non-existent / empty / no-ORCA-header files all return `False`; case-insensitive extension |
| `TestInitializeIdempotent` | Fresh context: exactly 1 opener, 1 handler |
| `TestContextRegistryAPI` | `StubContext` implements `register_window`, `get_window`, `mark_project_modified`; registry roundtrip; `initialize()` doesn't call absent methods |

---

### `test_freq_analysis.py` — FrequencyDialog logic (13 tests)

PyQt6, pyvista, rdkit, numpy, and PIL are fully stubbed.

| Class | What is tested |
|---|---|
| `TestPopulateList` | Imaginary modes (raw freq < 0) appear red; real modes use default colour; Unscaled column always shows the raw value |
| `TestUpdateData` | Re-applies (or removes) red colour based on unscaled value when scaling coefficients change; only the Scaled column (col 1) is updated; Unscaled column (col 2) stays fixed |

---

### `test_nmr_mo_fixes.py` — NMR timer safety and MO signal blocking (12 tests)

Tests for bug-fixes applied in v2.6.0. PyQt6, pyvista, rdkit, numpy, matplotlib,
and PIL are fully stubbed.

| Class | What is tested |
|---|---|
| `TestNmrCloseEvent` | `sel_timer` is stopped in `closeEvent` so it cannot fire on a dead widget |
| `TestNmrResetSelection` | `reset_selection()` exists and stops/restarts the timer while clearing selection state |
| `TestMoApplyPreset` | `apply_preset()` blocks signals on all vis widgets so only ONE `show_cube()` call is made per preset application (no phantom intermediate renders) |

---

### `test_plugin_integration.py` — PluginContext contract (17 tests)

Two execution modes:

**Stub mode** (always runs): `_StubContext` mirrors the real `PluginContext` API.

**Real-context mode** (runs when `python_molecular_editor` is present locally or
via `CI_MAIN_APP_SRC`): exercises `initialize()` against the actual `PluginContext`.

| Class | Tests | What is verified |
|---|---|---|
| `TestMetadata` | 2 | `PLUGIN_NAME` contains "ORCA"; `PLUGIN_VERSION` is semver X.Y.Z |
| `TestInitialize` | 5 | `.out` file opener registered (priority ≥ 100); drop handler registered (priority ≥ 100); both are callable |
| `TestDropHandler` | 5 | Non-.out and wrong-content files return `False`; ORCA banner file returns `True`; "Program Version" header returns `True` |
| `TestWithRealPluginContext` | 3 | `initialize(real_ctx)` doesn't raise; real ctx is `PluginContext` instance; `register_file_opener`, `register_drop_handler`, `get_main_window` exist on real class |

**CI:** the `test-integration` job in `.github/workflows/tests.yml` clones
`python_molecular_editor` and runs this file alongside `test_api.py` with the
real `PluginContext`.

---

### `test_api.py` — API checker infrastructure + integration (47 tests)

Tests `plugin_api_checker.py` itself using **synthetic code and temp directories** — no
main app required for the unit tests. All classes except `TestAPIChecker` run
unconditionally.

| Class | Tests | What is tested |
|---|---|---|
| `TestIssue` | 4 | `Issue` str format (`[try]` tag), `key()` 4-tuple, `in_try` excluded from key |
| `TestMergeAllowlists` | 5 | `mw`/`manager`/`context` set unions; empty merge; single dict preserved |
| `TestLoadSiteAllowlist` | 6 | List form, dict form, manager form, context form; missing file → `{}`; invalid JSON → `{}` |
| `TestPluginFileChecker` | 23 | Unknown MW attr flagged; known attr OK; private/Qt-inherited skipped; `hasattr` not flagged; try-block `in_try` flag; `get_main_window()` alias tracked; `self.main_window` alias tracked; MW/manager allowlist suppression; unknown/known manager attr; unknown/known context attr; context check off by default; plugin-specific attrs (`register_file_opener`, `register_drop_handler`, `register_window`, `get_window`) OK; syntax error; dedup |
| `TestAppAPIExtractor` | 8 | Extracts MW methods, properties, class attrs, manager attr, manager members, context members from synthetic app tree; scans `self.host.X` assignments |
| `TestAPIChecker` | 1 | Scans all plugin source files against the real MainWindow/PluginContext API (skipped unless main app present) |

---

### `test_utils.py` — Pure utility functions (6 tests)

Covers `get_default_export_path`: CSV/PNG extensions, custom suffix, empty
input, no-directory, directory preservation.

---

## Integration test strategy

```
Stub mode (always runs)
  _StubContext mirrors the real PluginContext API.
  No main app required. Catches interface mismatches immediately.

Real-context mode (runs when main app is available)
  Uses the actual PluginContext from python_molecular_editor.
  Catches drift invisible to the stub.
```

Real-context mode activates when repos are siblings:

```
<parent>/
    moleditpy_orca_result_analyzer_plugin/   ← this plugin
    python_molecular_editor/                 ← main app
```

### CI

| Job | Python | Main app | Tests |
|---|---|---|---|
| `test` | 3.11, 3.12, 3.13 | Cloned with `\|\| true` | Full suite; `TestWithRealPluginContext` + `TestAPIChecker` run if clone succeeded |
| `test-integration` | 3.11 | Cloned (hard failure if clone fails) | `test_plugin_integration.py` + `test_api.py` including real-context tier |

---

## Mocking strategy

| Module | Approach |
|---|---|
| `parser.py` | Loaded directly — no stubs (stdlib only) |
| `utils.py` | Loaded directly — no stubs |
| `vis.py` (`_parse_cube` only) | Cube-reading logic reimplemented in test — no PyVista |
| `traj_analysis.py` | Qt + matplotlib Qt backend stubbed; `QDialog`/`FigureCanvasQTAgg` as real stub classes |
| `freq_analysis.py` | Qt/numpy/pyvista/rdkit/PIL fully stubbed; `QDialog`/`QWidget` as real stub classes |
| `nmr_analysis.py` / `mo_analysis.py` | Qt/pyvista/numpy/matplotlib/PIL fully stubbed |
| `__init__.py` | Minimal PyQt6 (`QMessageBox`, `QFileDialog`) stubbed as `MagicMock()` |
| `plugin_api_checker.py` | No stubs — pure stdlib (`ast`, `json`, `pathlib`) |

---

## Coverage summary

| Module | Coverage | Notes |
|---|---|---|
| `utils.py` | **100%** | All branches hit |
| `parser.py` | **~79%** | All `parse_*` methods exercised; remaining gaps require new sample files |
| `traj_analysis.py` | **~12%** | `compute_scan_points`, `update_display_values`; Qt UI paths excluded |
| `__init__.py` | **~20%** | `initialize()`, `handle_drop()` False paths; `open_orca_file()` excluded |
| `vis.py` | partial | `_parse_cube` only; `_build_grid`/`show_iso` need PyVista + display |
| Qt UI modules | **0%** | `gui.py`, `*_analysis.py`, etc. require a live display |

### Known parser gaps (need new sample files)

| Gap | Reason |
|---|---|
| NEB `.xyz` multi-frame | No NEB sample file |
| NBO / FMO charges | No NBO analysis in samples |
| TDDFT short-table format B | Benzene uses arrow format only |
| Scan results table | No scan output in samples |
