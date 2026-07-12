"""Backward-compatibility shim for the old top-level `bse_solver` extension name.

The compiled module is now `chiq._bse_solver`. This shim re-exports it so the
legacy `import bse_solver; bse_solver.BSESolver` keeps working. Deprecated;
removed in ChiQ 2.0 -- use `chiq.solver.get_solver` instead.
"""
import warnings

from chiq._bse_solver import *  # noqa: F401,F403

warnings.warn(
    "Importing the top-level 'bse_solver' module is deprecated and will be "
    "removed in ChiQ 2.0; use chiq.solver.get_solver(...) instead.",
    FutureWarning,
    stacklevel=2,
)
