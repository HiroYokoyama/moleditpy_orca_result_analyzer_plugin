"""
tests/gui_harness.py
Shared headless Qt/matplotlib harness for instantiation-style coverage tests.

Installs a permissive PyQt6 stub: QtWidgets/QtCore/QtGui are MagicMock modules
so *any* widget import auto-resolves, while the specific classes used as base
classes in the source (QDialog, QWidget, QTableWidget, QThread, QSlider,
QObject, QLabel) are real, subclassable Python classes whose instances answer
arbitrary attribute access with a MagicMock. The matplotlib Qt-agg backend is
likewise stubbed with a real FigureCanvasQTAgg base.

This mirrors the stub that tests/test_about_menu.py sets up, generalised so the
whole plugin package imports cleanly and its dialogs can be constructed with the
Qt constructor running (which is what drives the bulk of the coverage).

Import this module *before* importing any orca_result_analyzer submodule.
"""

import sys
import types
from unittest.mock import MagicMock


class _PermissiveMeta(type):
    """Answer missing *class-level* attribute access (enums like
    QTableWidget.EditTrigger, QHeaderView.ResizeMode) with a MagicMock."""

    def __getattr__(cls, name):
        if name.startswith("__") and name.endswith("__"):
            raise AttributeError(name)
        return MagicMock()


class _Permissive(metaclass=_PermissiveMeta):
    """Base for stubbed Qt widgets: swallow constructor args, answer any
    attribute (methods, signals) with a fresh MagicMock."""

    class DialogCode:
        Accepted = 1
        Rejected = 0

    def __init__(self, *a, **k):
        pass

    def __getattr__(self, name):
        if name.startswith("__") and name.endswith("__"):
            raise AttributeError(name)
        return MagicMock()


def _ensure_module(name):
    mod = sys.modules.get(name)
    if mod is None or not isinstance(mod, types.ModuleType):
        # A MagicMock module gives us auto-resolving `from mod import Anything`.
        mod = MagicMock()
        mod.__name__ = name
        sys.modules[name] = mod
    return mod


def _set_real_class(mod, name, bases=(_Permissive,), members=None):
    """Force *mod.name* to be a real subclassable class (overriding any
    previously-installed MagicMock or partial stub)."""
    existing = getattr(mod, name, None)
    if isinstance(existing, type) and issubclass(existing, _Permissive):
        return existing
    cls = type(name, bases, dict(members or {}))
    setattr(mod, name, cls)
    return cls


def install():
    _ensure_module("PyQt6")
    qtw = _ensure_module("PyQt6.QtWidgets")
    qtc = _ensure_module("PyQt6.QtCore")
    _ensure_module("PyQt6.QtGui")
    sys.modules["PyQt6"].QtWidgets = qtw
    sys.modules["PyQt6"].QtCore = qtc

    # Real, subclassable base classes for every class subclassed in the source.
    for name in ["QDialog", "QWidget", "QTableWidget", "QSlider", "QLabel"]:
        _set_real_class(qtw, name)
    for name in ["QThread", "QObject"]:
        _set_real_class(qtc, name)

    # matplotlib Qt-agg backend: FigureCanvasQTAgg / NavigationToolbar2QT are
    # subclassed by the spectrum / trajectory / frequency widgets.
    backend = _ensure_module("matplotlib.backends.backend_qtagg")
    _set_real_class(backend, "FigureCanvasQTAgg")
    _set_real_class(backend, "NavigationToolbar2QT")

    # A functional QColor so color-interpolation math in the charge / dipole /
    # energy-diagram dialogs runs for real instead of tripping over MagicMock
    # arithmetic and diving straight into exception handlers.
    _install_qcolor(sys.modules["PyQt6.QtGui"])

    # pyvista is real if installed; only stub when absent.
    if "pyvista" not in sys.modules:
        try:
            import pyvista  # noqa: F401
        except Exception:
            sys.modules["pyvista"] = MagicMock()


_NAMED = {
    "red": (255, 0, 0),
    "green": (0, 128, 0),
    "blue": (0, 0, 255),
    "white": (255, 255, 255),
    "black": (0, 0, 0),
    "purple": (128, 0, 128),
    "cyan": (0, 255, 255),
    "yellow": (255, 255, 0),
    "gray": (128, 128, 128),
}


class _QColor:
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


def _install_qcolor(qtg):
    if not isinstance(getattr(qtg, "QColor", None), type):
        qtg.QColor = _QColor


install()
