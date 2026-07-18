"""
tests/gui_harness.py
Isolated headless-Qt loader for instantiation-style coverage tests.

Every other test module in this suite installs its own partial PyQt6 stub into
the shared ``sys.modules`` and relies on those exact objects at run time. So we
must NOT mutate the shared stubs in a way that outlives our own imports.

``load_isolated(modname)`` therefore:

  1. temporarily installs a complete set of real, subclassable Qt base classes
     (plus a matplotlib Qt-agg backend and a functional QColor) onto the shared
     stub modules, remembering whatever was there before;
  2. loads ``orca_result_analyzer/<modname>.py`` under a private throwaway
     package whose ``__path__`` points at the real source dir, so the module's
     ``from .sibling import x`` relative imports auto-resolve to fresh copies of
     the siblings (also bound to our complete classes) without executing the
     heavy package ``__init__``;
  3. restores the shared stub modules to exactly what they were.

The loaded module keeps the classes it captured at import time, so it works at
run time while the shared stubs the other tests depend on are left untouched.
"""

import os
import sys
import types
import importlib.util
from unittest.mock import MagicMock

_PKG_DIR = os.path.normpath(
    os.path.join(os.path.dirname(__file__), "..", "orca_result_analyzer")
)


class _PermissiveMeta(type):
    def __getattr__(cls, name):
        if name.startswith("__") and name.endswith("__"):
            raise AttributeError(name)
        return MagicMock()


_NOOP_METHODS = (
    "closeEvent", "showEvent", "resizeEvent", "keyPressEvent", "paintEvent",
    "moveEvent", "accept", "reject", "done", "close", "setWindowTitle",
    "resize", "setLayout", "setAcceptDrops",
)


def _noop(self, *a, **k):
    return None


def _make_widget(name):
    ns = {m: _noop for m in _NOOP_METHODS}

    class DialogCode:
        Accepted = 1
        Rejected = 0

    ns["DialogCode"] = DialogCode

    def __init__(self, *a, **k):
        pass

    def __getattr__(self, attr):
        if attr.startswith("_"):
            raise AttributeError(attr)
        return MagicMock()

    ns["__init__"] = __init__
    ns["__getattr__"] = __getattr__
    return _PermissiveMeta(name, (object,), ns)


# Widget base classes subclassed anywhere in the source.
_QTW_BASES = ["QDialog", "QWidget", "QTableWidget", "QSlider", "QLabel", "QMainWindow"]
_QTC_BASES = ["QThread", "QObject"]
_BACKEND_BASES = ["FigureCanvasQTAgg", "NavigationToolbar2QT"]


def _ensure_module(name):
    mod = sys.modules.get(name)
    if mod is None:
        mod = MagicMock()
        mod.__name__ = name
        sys.modules[name] = mod
    return mod


_SENTINEL = object()


def _fresh_module(name, bases, base_factory, extra=None):
    """A MagicMock module (so `from mod import Anything` auto-resolves) carrying
    real, subclassable base classes plus any extras."""
    mod = MagicMock()
    mod.__name__ = name
    for b in bases:
        setattr(mod, b, base_factory(b))
    for k, v in (extra or {}).items():
        setattr(mod, k, v)
    return mod


def _pyvista_module():
    needed = ("Arrow", "Box", "Sphere", "Line", "PolyData")
    pv = sys.modules.get("pyvista")
    if pv is not None and all(hasattr(pv, n) for n in needed):
        return pv
    try:
        import pyvista as real  # noqa: F401

        return real
    except Exception:
        return MagicMock()


_load_counter = [0]

# sys.modules keys we swap for the duration of an isolated load.
_SWAP_KEYS = [
    "PyQt6.QtWidgets",
    "PyQt6.QtCore",
    "PyQt6.QtGui",
    "matplotlib.backends.backend_qtagg",
    "pyvista",
]


def load_isolated(modname):
    _ensure_module("PyQt6")

    swapped = {
        "PyQt6.QtWidgets": _fresh_module("PyQt6.QtWidgets", _QTW_BASES, _make_widget),
        "PyQt6.QtCore": _fresh_module("PyQt6.QtCore", _QTC_BASES, _make_widget),
        "PyQt6.QtGui": _fresh_module("PyQt6.QtGui", [], _make_widget,
                                     extra={"QColor": QColorStub}),
        "matplotlib.backends.backend_qtagg": _fresh_module(
            "matplotlib.backends.backend_qtagg", _BACKEND_BASES, _make_widget),
        "pyvista": _pyvista_module(),
    }
    saved = {k: sys.modules.get(k, _SENTINEL) for k in _SWAP_KEYS}
    for k, v in swapped.items():
        sys.modules[k] = v

    _load_counter[0] += 1
    pkg_name = f"_ora_iso_{_load_counter[0]}"
    pkg = types.ModuleType(pkg_name)
    pkg.__path__ = [_PKG_DIR]
    sys.modules[pkg_name] = pkg
    try:
        spec = importlib.util.spec_from_file_location(
            f"{pkg_name}.{modname}", os.path.join(_PKG_DIR, f"{modname}.py")
        )
        mod = importlib.util.module_from_spec(spec)
        mod.__package__ = pkg_name
        sys.modules[f"{pkg_name}.{modname}"] = mod
        spec.loader.exec_module(mod)
        return mod
    finally:
        for k in _SWAP_KEYS:
            if saved[k] is _SENTINEL:
                sys.modules.pop(k, None)
            else:
                sys.modules[k] = saved[k]


# ---------------------------------------------------------------------------
# Functional QColor for color-interpolation math.
# ---------------------------------------------------------------------------

_NAMED = {
    "red": (255, 0, 0), "green": (0, 128, 0), "blue": (0, 0, 255),
    "white": (255, 255, 255), "black": (0, 0, 0), "purple": (128, 0, 128),
    "cyan": (0, 255, 255), "yellow": (255, 255, 0), "gray": (128, 128, 128),
}


class QColorStub:
    def __init__(self, *a):
        r = g = b = 0
        if len(a) == 1 and isinstance(a[0], str):
            s = a[0].strip().lower()
            if s.startswith("#") and len(s) == 7:
                r, g, b = int(s[1:3], 16), int(s[3:5], 16), int(s[5:7], 16)
            else:
                r, g, b = _NAMED.get(s, (0, 0, 0))
        elif len(a) >= 3:
            r, g, b = int(a[0]), int(a[1]), int(a[2])
        self._r, self._g, self._b = r, g, b

    def red(self):
        return self._r

    def green(self):
        return self._g

    def blue(self):
        return self._b

    def name(self):
        return f"#{self._r:02x}{self._g:02x}{self._b:02x}"

    def isValid(self):
        return True

    def getRgb(self):
        return (self._r, self._g, self._b, 255)
