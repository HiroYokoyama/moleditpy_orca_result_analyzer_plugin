"""
Tests for the Help > About menu action.
"""

import sys
import types
import unittest
from unittest.mock import MagicMock

# Stub Qt
if "PyQt6" not in sys.modules:
    pyqt6 = MagicMock()
    qtw = MagicMock()
    
    # We only need QDialog to be a class that we can instantiate
    class _QDialog:
        def __init__(self, *a, **k):
            self.menu_bar_mock = MagicMock()
        def menuBar(self):
            return self.menu_bar_mock
        def setWindowTitle(self, *a): pass
        def resize(self, *a): pass
        def setAcceptDrops(self, *a): pass
        def close(self): pass
        
    qtw.QDialog = _QDialog
        
    sys.modules["PyQt6"] = pyqt6
    sys.modules["PyQt6.QtWidgets"] = qtw
    sys.modules["PyQt6.QtCore"] = MagicMock()
    sys.modules["PyQt6.QtGui"] = MagicMock()

    # Stub matplotlib backends
    sys.modules["matplotlib.backends.backend_qtagg"] = MagicMock()
    sys.modules["pyvista"] = MagicMock()

from orca_result_analyzer.gui import OrcaResultAnalyzerDialog
from orca_result_analyzer.parser import OrcaParser


class TestAboutMenu(unittest.TestCase):
    @unittest.mock.patch("PyQt6.QtWidgets.QMessageBox.about", create=True)
    def test_show_about(self, mock_about):
        """Test that show_about calls QMessageBox.about."""
        parser = OrcaParser()
        # Mock Context
        ctx = MagicMock()

        # Instantiate Dialog
        dlg = OrcaResultAnalyzerDialog(None, parser, "", ctx)

        # Call show_about directly
        dlg.show_about()

        # Verify it called QMessageBox.about
        mock_about.assert_called_once()
        args = mock_about.call_args[0]

        self.assertIs(args[0], dlg)
        self.assertEqual(args[1], "About ORCA Result Analyzer")
        self.assertIn("Author:</b> Hiromichi Yokoyama", args[2])
        self.assertIn("GitHub:", args[2])


if __name__ == "__main__":
    unittest.main()
