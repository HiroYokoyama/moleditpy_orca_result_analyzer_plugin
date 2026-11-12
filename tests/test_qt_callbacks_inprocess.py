"""The queued-callback helper's own logic, in-process so coverage sees it.

test_qt_callbacks.py proves the Qt behaviour with real PySide6 in child
processes, which coverage cannot follow. These tests load the same file
against a small stand-in QtCore and check what the helper itself decides.
"""

import importlib.util
import sys
import types
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

HELPER = Path(__file__).resolve().parents[1] / "orca_result_analyzer/_qt_callbacks.py"


def load_helper(application):
    core = types.ModuleType("PySide6.QtCore")

    class QObject:
        def __init__(self, parent=None):
            self.parent = parent

        def deleteLater(self):
            pass

    core.QObject = QObject
    core.QCoreApplication = MagicMock()
    core.QCoreApplication.instance.return_value = application
    core.Qt = MagicMock()
    core.Slot = lambda *types_: (lambda function: function)
    package = types.ModuleType("PySide6")
    package.QtCore = core
    with patch.dict(sys.modules, {"PySide6": package, "PySide6.QtCore": core}):
        spec = importlib.util.spec_from_file_location("qt_callbacks_inprocess", HELPER)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    return module, core


class InProcessCallbackTests(unittest.TestCase):
    def test_dispatch_hands_the_signal_arguments_to_the_callback(self):
        callbacks, core = load_helper(application=object())
        received, signal, terminal = [], MagicMock(), MagicMock()
        receiver = callbacks.connect_queued(
            signal, lambda *args: received.append(args), signature=(int, str),
            cleanup=(terminal,),
        )
        signal.connect.assert_called_once_with(
            receiver.dispatch, core.Qt.ConnectionType.QueuedConnection
        )
        terminal.connect.assert_called_once_with(receiver.deleteLater)
        receiver.dispatch(7, "seven")
        self.assertEqual(received, [(7, "seven")])

    def test_without_an_owner_the_application_owns_the_receiver(self):
        application = object()
        callbacks, _ = load_helper(application)
        receiver = callbacks.connect_queued(MagicMock(), lambda: None)
        self.assertIs(receiver.parent, application)

    def test_no_application_is_reported_instead_of_an_orphan_receiver(self):
        callbacks, _ = load_helper(application=None)
        signal = MagicMock()
        with self.assertRaisesRegex(RuntimeError, "Qt application"):
            callbacks.connect_queued(signal, lambda: None)
        signal.connect.assert_not_called()

    def test_receiver_classes_are_shared_per_signature(self):
        callbacks, _ = load_helper(application=object())
        first = callbacks.connect_queued(MagicMock(), lambda value: None, signature=(int,))
        second = callbacks.connect_queued(MagicMock(), lambda value: None, signature=(int,))
        third = callbacks.connect_queued(MagicMock(), lambda value: None, signature=(str,))
        self.assertIs(type(first), type(second))
        self.assertIsNot(type(first), type(third))


if __name__ == "__main__":
    unittest.main()
