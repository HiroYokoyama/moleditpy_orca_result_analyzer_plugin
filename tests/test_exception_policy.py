"""
tests/test_exception_policy.py
Enforces CONTRIBUTING.md section 4 B/C on the package source.

Section 4B requires UI slots and callbacks to catch broadly so a slot can
never crash the app, and to log or report what they caught. Section 4C
requires internal helper code to name specific exception types instead.

These are structural checks over the AST, so they cost nothing at runtime and
fail the moment a new blind handler is added. The broad-handler count is a
ratchet: it may fall freely, but raising it is a deliberate act that has to
change this file too.
"""

import ast
import os
import unittest

_PKG = os.path.normpath(
    os.path.join(os.path.dirname(__file__), "..", "orca_result_analyzer")
)

# Ratchet. Lower it whenever handlers are narrowed; never raise it without a
# reason in the commit message. Was 212 before the 3.12.2 narrowing pass.
_MAX_BROAD_HANDLERS = 49

# parser.py is pure text parsing with no UI slots, so every handler in it can
# name its types. Verified bit-identical over the 19 sample outputs.
_FULLY_NARROWED = {"parser.py"}


def _handlers():
    """Yield (filename, lineno, caught_names, handler_node, try_node)."""
    for name in sorted(os.listdir(_PKG)):
        if not name.endswith(".py"):
            continue
        path = os.path.join(_PKG, name)
        with open(path, encoding="utf-8") as fh:
            src = fh.read()
        for node in ast.walk(ast.parse(src)):
            if not isinstance(node, ast.Try):
                continue
            for h in node.handlers:
                t = h.type
                if t is None:
                    caught = None  # bare except:
                elif isinstance(t, ast.Name):
                    caught = {t.id}
                elif isinstance(t, ast.Tuple):
                    caught = {e.id for e in t.elts if isinstance(e, ast.Name)}
                else:
                    caught = {ast.unparse(t)}
                yield name, h.lineno, caught, h, node


def _is_broad(caught):
    return caught is None or bool(caught & {"Exception", "BaseException"})


class TestExceptionPolicy(unittest.TestCase):
    def test_no_bare_except_anywhere(self):
        """`except:` also swallows KeyboardInterrupt and SystemExit."""
        bare = [
            f"{name}:{ln}" for name, ln, caught, _h, _t in _handlers() if caught is None
        ]
        self.assertEqual(bare, [], f"bare `except:` found at {bare}")

    def test_parser_names_every_exception_it_catches(self):
        offenders = [
            f"{name}:{ln}"
            for name, ln, caught, _h, _t in _handlers()
            if name in _FULLY_NARROWED and _is_broad(caught)
        ]
        self.assertEqual(
            offenders,
            [],
            "parser.py is internal helper code with no UI slots; every handler "
            f"must name its types. Broad handlers at {offenders}",
        )

    def test_broad_handler_count_does_not_grow(self):
        broad = [
            f"{name}:{ln}" for name, ln, caught, _h, _t in _handlers() if _is_broad(caught)
        ]
        self.assertLessEqual(
            len(broad),
            _MAX_BROAD_HANDLERS,
            f"{len(broad)} broad handlers, ceiling is {_MAX_BROAD_HANDLERS}. "
            "Name the exception types, or raise the ratchet deliberately.",
        )

    def test_the_ratchet_is_tight(self):
        """A stale ceiling stops catching regressions, so keep it exact."""
        broad = [1 for _n, _l, caught, _h, _t in _handlers() if _is_broad(caught)]
        self.assertEqual(
            len(broad),
            _MAX_BROAD_HANDLERS,
            "handlers were narrowed without lowering _MAX_BROAD_HANDLERS",
        )

    def test_no_handler_silently_discards_its_exception(self):
        """CONTRIBUTING 4B: empty except blocks are not acceptable."""
        silent = []
        for name, ln, _caught, h, _t in _handlers():
            body = [n for n in h.body if not isinstance(n, ast.Pass)]
            if not body:
                silent.append(f"{name}:{ln}")
        self.assertEqual(
            silent, [], f"handlers that only `pass` (no log, no report): {silent}"
        )

    def test_broad_handlers_log_or_report(self):
        """What survives as broad must at least say so."""
        mute = []
        for name, ln, caught, h, _t in _handlers():
            if not _is_broad(caught):
                continue
            text = " ".join(ast.unparse(n) for n in h.body)
            reports = any(
                k in text
                for k in (
                    "logging.",
                    "logger.",
                    "QMessageBox",
                    "show_status_message",
                    "finished_sig",
                    "raise",
                    "self.fail",
                )
            )
            # A handler that only substitutes a fallback value is fine.
            if not reports and not any(
                isinstance(n, (ast.Assign, ast.Return, ast.FunctionDef)) for n in h.body
            ):
                mute.append(f"{name}:{ln}")
        self.assertEqual(mute, [], f"broad handlers that neither log nor report: {mute}")


if __name__ == "__main__":
    unittest.main()
