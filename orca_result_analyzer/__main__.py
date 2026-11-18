"""Entry point for ``moleditpy -r orca_result_analyzer`` and ``python -m``."""

import sys

from .cli import main

sys.exit(main())
