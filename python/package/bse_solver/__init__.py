"""Backward-compatibility shim for the old top-level `bse_solver` extension name.

The compiled module is now `chiq._bse_solver`. This shim re-exports it so the
legacy `import bse_solver; bse_solver.BSESolver` keeps working. Deprecated;
removed in ChiQ 2.0 -- use `chiq.solver.get_solver` instead.
"""
import importlib
import warnings

try:
    _impl = importlib.import_module("chiq._bse_solver")
except ModuleNotFoundError as exc:
    if exc.name != "chiq._bse_solver":
        raise
    _impl = importlib.import_module("_bse_solver")

BSESolver = _impl.BSESolver

__all__ = ["BSESolver"]

try:
    _warning_emitted
except NameError:
    _warning_emitted = False

if not _warning_emitted:
    warnings.warn(
        "Importing the top-level 'bse_solver' module is deprecated and will be "
        "removed in ChiQ 2.0; use chiq.solver.get_solver(...) instead.",
        FutureWarning,
        stacklevel=2,
    )
    _warning_emitted = True
