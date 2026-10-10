"""Backward-compatibility shim for the old `bse` package name.

`bse` was renamed to `chiq`. This shim forwards a fixed allowlist of public
submodules to `chiq.*` (deterministic, IDE/type-checker friendly). It is
deprecated and will be REMOVED in ChiQ 2.0 -- import from `chiq` instead.
"""
import warnings

try:
    _warning_emitted
except NameError:
    _warning_emitted = False

if not _warning_emitted:
    warnings.warn(
        "The 'bse' package is deprecated and will be removed in ChiQ 2.0; "
        "import from 'chiq' instead (e.g. `from chiq import h5bse`).",
        FutureWarning,
        stacklevel=2,
    )
    _warning_emitted = True
