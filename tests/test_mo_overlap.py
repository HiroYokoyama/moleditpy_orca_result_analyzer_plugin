"""
Tests for orbital energy diagram overlap shifting logic.
"""

import sys
from unittest.mock import MagicMock

# Stub PyQt6 to prevent DLL loading in test environment
if "PyQt6" not in sys.modules:
    pyqt6 = MagicMock()
    sys.modules["PyQt6"] = pyqt6
    sys.modules["PyQt6.QtWidgets"] = pyqt6.QtWidgets
    sys.modules["PyQt6.QtCore"] = pyqt6.QtCore
    sys.modules["PyQt6.QtGui"] = pyqt6.QtGui

from orca_result_analyzer.energy_diag import calculate_arrow_shifts


def test_calculate_arrow_shifts_no_overlap():
    def mock_val_to_y(energy):
        # Higher energy -> smaller Y coordinate (higher up on screen)
        return 1000 - energy * 100

    items = [(0, 3.0, 1.0), (1, 2.0, 1.0), (2, 1.0, 1.0)]

    shifts = calculate_arrow_shifts(items, mock_val_to_y)
    assert shifts == {0: 0, 1: 0, 2: 0}


def test_calculate_arrow_shifts_two_overlap():
    def mock_val_to_y(energy):
        return 1000 - energy * 100

    items = [(0, 3.0, 1.0), (1, 2.9, 1.0), (2, 1.0, 1.0)]

    shifts = calculate_arrow_shifts(items, mock_val_to_y)
    assert shifts == {0: 10, 1: 0, 2: 0}


def test_calculate_arrow_shifts_multiple_overlap():
    def mock_val_to_y(energy):
        return 1000 - energy * 100

    items = [(0, 3.0, 1.0), (1, 2.95, 1.0), (2, 2.90, 1.0)]

    shifts = calculate_arrow_shifts(items, mock_val_to_y)
    assert shifts == {0: 20, 1: 10, 2: 0}
