"""Reading and parsing ORCA output behind a progress dialog.

A large ``.out`` is read and parsed on the GUI thread, so without feedback the
application looks hung.
"""

import logging
import os

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QApplication, QProgressDialog

from .parser import OrcaParser, ParseCancelled

ENCODINGS = ("utf-8", "utf-16", "latin-1", "cp1252")

#: Reading counts as this fraction of the whole load in the progress bar.
_READ_FRACTION = 0.25

#: Don't flash a dialog for a small file that loads instantly (ms).
_MIN_DURATION_MS = 300


def read_orca_text(path):
    """Return the text of *path*, trying several encodings.

    Undecodable bytes are replaced rather than failing the load; I/O errors
    propagate.
    """
    for enc in ENCODINGS:
        try:
            with open(path, "r", encoding=enc) as f:
                return f.read()
        except UnicodeError:
            continue
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        return f.read()


class LoadProgress:
    """Modal progress dialog for one file load."""

    def __init__(self, parent, path):
        self._cancelled = False
        self._dlg = QProgressDialog(
            f"Reading {os.path.basename(path)}...", "Cancel", 0, 100, parent
        )
        self._dlg.setWindowTitle("Loading ORCA Output")
        self._dlg.setWindowModality(Qt.WindowModality.WindowModal)
        self._dlg.setMinimumDuration(_MIN_DURATION_MS)
        self._dlg.setAutoClose(False)
        self._dlg.setAutoReset(False)
        # A signal, not wasCanceled(): polling a permissive test stub would
        # read as cancelled and abort every load.
        try:
            self._dlg.canceled.connect(self._on_cancel)
        except (AttributeError, TypeError) as exc:
            logging.debug("LoadProgress: cancel not wired — %s", exc)
        self._dlg.setValue(0)
        QApplication.processEvents()

    @property
    def cancelled(self):
        return self._cancelled

    def _on_cancel(self):
        self._cancelled = True

    def _pump(self, value, label):
        self._dlg.setLabelText(label)
        self._dlg.setValue(int(value))
        QApplication.processEvents()

    def read_done(self, path):
        self._pump(_READ_FRACTION * 100, f"Parsing {os.path.basename(path)}...")

    def on_parse_step(self, done, total, label):
        if self._cancelled:
            raise ParseCancelled(label)
        frac = _READ_FRACTION + (1.0 - _READ_FRACTION) * (done / max(total, 1))
        self._pump(frac * 100, label if done >= total else f"{label}...")

    def close(self):
        try:
            self._dlg.hide()
            self._dlg.deleteLater()
        except (AttributeError, RuntimeError) as exc:
            import logging

            logging.debug("LoadProgress: hide/delete failed — %s", exc)
        self._dlg.close()


def load_orca_parser(path, parent):
    """Read and parse *path* behind a progress dialog.

    Returns a loaded :class:`OrcaParser`, or ``None`` if the user cancelled.
    I/O errors propagate to the caller, which owns the error reporting.
    """
    progress = LoadProgress(parent, path)
    try:
        content = read_orca_text(path)
        if progress.cancelled:
            return None
        progress.read_done(path)

        parser = OrcaParser()
        try:
            parser.load_from_memory(content, path, progress=progress.on_parse_step)
        except ParseCancelled:
            return None
        return parser
    finally:
        progress.close()
