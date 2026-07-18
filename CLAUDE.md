# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What this is

A MoleditPy plugin that parses ORCA quantum-chemistry output (`.out`) and
presents it as analysis dialogs — MO visualisation, vibrational modes, NMR,
TD-DFT, charges, forces, trajectories. PyQt6 for the UI, pyvista for the 3D
viewer, matplotlib for plots, rdkit for molecule handling.

`orca_result_analyzer/` is the whole plugin; there is no `src/` layout. It ships
as a zip of that directory.

## Layout

- `__init__.py` — plugin entry points (`initialize`, `run`) and the
  `PLUGIN_VERSION` constant. Registers windows through the host's context
  registry rather than holding references itself.
- `parser.py` — all ORCA text parsing. `load_from_memory(content, filename)`
  runs `parse_all()`, filling a single `self.data` dict; every dialog reads from
  that dict. Adding support for a new ORCA block means a new `parse_*` method
  called from `parse_all()`, writing one more key. Ordering matters there —
  gradients are parsed before the trajectory so they can be linked to it.
- `gui.py` — `OrcaResultAnalyzerDialog`, the main window. Owns the file loading,
  the 3D structure, atom picking, and one `show_*` launcher per analysis dialog.
- One module per analysis type (`nmr_analysis.py`, `mo_analysis.py`, …), each a
  self-contained `QDialog`.
- `utils.py` — shared helpers, notably `save_json_atomic` and
  `get_default_export_path`.

## Version bumping

`PLUGIN_VERSION` in `orca_result_analyzer/__init__.py` is the single source of
truth. The release workflow verifies it matches the pushed tag and fails
otherwise, so bump it in the same commit as the fix it describes.

## Testing

```bash
python run_tests.py                 # whole suite, exactly as CI runs it
python run_tests.py -k parser       # pytest args are forwarded
```

**Always use `run_tests.py`, not `pytest` directly.** It sets
`PYTEST_DISABLE_PLUGIN_AUTOLOAD=1` so local runs match CI, which installs no Qt
binding. An installed `pytest-qt` would otherwise import a real PyQt6 during
collection and defeat the stubs the tests rely on — producing failures and
segfaults that do not reproduce in CI.

Coverage (plugin autoload is off, so `-p pytest_cov` is required):

```bash
python run_tests.py -p pytest_cov --cov=orca_result_analyzer --cov-report=term-missing
```

Target is above 80%. CI reports it but does not gate on it.

### Python version

CI runs 3.11–3.13. The source uses 3.10+ syntax (`str | None`, `list[str]`), so
on an older interpreter `tests/test_api.py` fails at collection. If you are on
3.9, `--ignore=tests/test_api.py` lets the rest of the suite run, but treat CI
as authoritative.

## Writing dialog tests

The dialogs need Qt and a live 3D viewer at import time. `tests/gui_harness.py`
provides `load_isolated("module_name")`, which loads one source module against
subclassable Qt stand-ins that keep real state (spin boxes, combos, buttons,
trees, tables, lists remember what you set), then restores the shared
`sys.modules`. `tests/test_nmr_dialog.py` is a representative example.

Traps worth knowing, each of which has cost real debugging time here:

**The widget stubs are permissive.** An unknown attribute returns a MagicMock
instead of raising, so a mistyped method name in a test *passes silently*.
Assert on observable state, not merely that a call happened.

**Settings paths come from the module directory**, sometimes recomputed inside
each method rather than stored on the dialog. Redirect the loaded module's
`__file__` at a temp dir or the suite writes a `settings.json` into the package
source tree (it is gitignored, and was deliberately untracked, so it is easy to
miss):

```python
saved = M.__file__
M.__file__ = os.path.join(self.tmp, "mod.py")
self.addCleanup(lambda: setattr(M, "__file__", saved))
```

**Never install a crippled third-party module into `sys.modules` permanently.**
Several test modules stub `numpy`, `matplotlib` and `PIL`; those stubs must be
*fallbacks used only when the real package is absent*. A stub that overwrites a
working package breaks every later test doing real numerics — matplotlib
resolves `numpy` and its Agg backend from `sys.modules` at call time, long after
import.

**Qt idioms that hang against a MagicMock.** `while layout.count(): takeAt(0)`
and `while it.value(): ...; it += 1` never terminate when the counter is a
MagicMock. The harness provides stateful layouts and a terminating tree iterator
for this reason; if you add a stub, make sure loop conditions can go falsy.

**Patch where the name is bound.** Some methods import Qt lazily inside the
function body (`from PyQt6.QtWidgets import QInputDialog`) — wrap those calls in
`gui_harness.qt_available(QInputDialog=fake)`. Others import at module level —
patch the loaded module's attribute (`patch.object(M, "QInputDialog", fake)`).
Sibling dialogs imported lazily (`from .bond_analysis import ...`) are patched
through `sys.modules[f"{G.__package__}.bond_analysis"]`.

**Avoid patching `builtins.open`** — it breaks pytest's own I/O and hangs the
run. Point at a genuinely unwritable path instead (e.g. a path under a file).

**Keep exporter tests off the rasteriser.** Assert that `savefig` was driven
rather than writing real bytes; a real render pulls in Agg and Pillow, which
other test modules stub out.

## Fixture shapes

Match what `parser.py` actually produces, or tests pass against fiction:

- MO dicts are keyed `"0_alpha"` but carry a numeric `"id"`.
- Convergence entries store the raw ORCA verdict string (`"YES"`/`"NO"`), not a
  bool.
- `spin_s2` is a dict (`actual`/`ideal`/`contamination`), not a float.
- NBO hybrids are dicts with `atom_sym`, `atom_idx`, `s_pct`, `p_pct`, `d_pct`,
  `label`.
- FMO rows label atoms with a single `"0-C"` token.

## Conventions

Error handling policy is in `CONTRIBUTING.md` and is enforced in review: never
hide errors, never crash. UI slots and callbacks catch, log, and tell the user;
internal helpers propagate. Empty `except` blocks are not acceptable.
