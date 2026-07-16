"""
tests/test_reject_close.py
Static rule: every QDialog subclass that does cleanup in closeEvent must
also override reject() — QDialog.reject (Esc key) only hides the dialog,
so without the override Esc skips timer stops, 3D actor/label removal,
settings saves, and the NMR unsaved-merge prompt.
"""

import ast
import glob
import os
import unittest

_SRC_DIR = os.path.normpath(
    os.path.join(os.path.dirname(__file__), "..", "orca_result_analyzer")
)


def _dialog_classes(tree):
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef):
            base_names = {
                b.id if isinstance(b, ast.Name) else getattr(b, "attr", None)
                for b in node.bases
            }
            if "QDialog" in base_names:
                yield node


class TestRejectRoutesThroughClose(unittest.TestCase):
    def test_every_qdialog_with_closeevent_overrides_reject(self):
        offenders = []
        for path in glob.glob(os.path.join(_SRC_DIR, "*.py")):
            with open(path, encoding="utf-8") as f:
                tree = ast.parse(f.read(), filename=path)
            for cls in _dialog_classes(tree):
                methods = {
                    n.name for n in cls.body if isinstance(n, ast.FunctionDef)
                }
                if "closeEvent" in methods and "reject" not in methods:
                    offenders.append(f"{os.path.basename(path)}:{cls.name}")
        self.assertEqual(
            offenders,
            [],
            "QDialog subclasses with closeEvent cleanup but no reject() "
            f"override (Esc would skip the cleanup): {offenders}",
        )


if __name__ == "__main__":
    unittest.main()
