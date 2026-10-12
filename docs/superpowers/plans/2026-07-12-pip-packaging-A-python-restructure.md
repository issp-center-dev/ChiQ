# pip packaging — Plan A: Python restructure Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Restructure the ChiQ Python layer so it is `pip`-packageable — move the CLI into an importable `chiq.cli` subpackage, replace the `bse` symlink with a static forwarding-shim package, and make `matplotlib` a lazy (optional) import — without any build-system change and without breaking the existing CMake workflow or tests.

**Architecture:** Pure-Python reorganization. CLI `main()` implementations move from standalone `python/scripts/*.py` into `chiq/cli/<name>.py` (so PEP 621 `console_scripts` can target them in Plan B); the old script files become thin wrappers. `python/package/bse` becomes a real package of thin `from chiq.<mod> import *` forwarding modules. `chiq.g2scl_core` and the plotting commands import `matplotlib` lazily so a core install needs no matplotlib.

**Tech Stack:** Python ≥3.8, pytest. No C++/CMake/pyproject changes in this plan (those are Plan B).

## Global Constraints

- **Python floor 3.8** — no 3.9+ syntax (no PEP 585 generics, no `str.removeprefix`, etc.).
- **No behavior change** to what commands compute; only their location/import structure changes. The existing integration tests (`tests/python/non-mpi/bsetool_*`, `gen_qpath`, `bsepost_linear_combination`) must stay green — they invoke `chiq_main.py`, `chiq_post.py`, `gen_qpath.py` by name via `os.system`, so those script names must keep working (as wrappers).
- **`bse` shim = static forwarding modules**, allowlisted, lazy (importing `bse` must not import `mpi4py`/`dcore`/`matplotlib`). Allowlist: `h5bse`, `bse_toml`, `matrix_dict`, `point_group`, `sumk_dft_chi`, `g2scl_core`, `tools`, `index_pair`, `mpi`, and the `point_group_data` subpackage (`base, C1, C2, D3, D4, D6, O, Oh`).
- **`import *` forwards only public names** (respects each module's `__all__`); this is the documented compatibility surface.
- **Deprecation is visible:** `import bse` emits a `FutureWarning` (shown by default, unlike `DeprecationWarning`) naming the replacement and removal version **2.0**.
- **`chiq._bse_solver` does NOT exist yet** (that rename is Plan B). This plan does not touch the extension or `chiq/solver/cpp.py`; the top-level `bse_solver` shim is Plan B.
- **Do not touch** `python/CMakeLists.txt` install rules in this plan (Plan B). The legacy CMake build still globs `python/scripts/*.py` — the wrappers keep that working.
- Spec: `docs/superpowers/specs/2026-07-12-pip-packaging-design.md`.
- **Test env:** `PYTHONPATH="$REPO/python/package:$REPO/build/src" python3 -m pytest tests/python/non-mpi/... -q` (the `build/src` entry provides the built `bse_solver` extension for the integration tests; source `chiqvars.sh`-equivalent PATH for the `os.system` script calls: `export PATH="$REPO/python/scripts:$PATH"`).

## File Structure

- New: `python/package/chiq/cli/__init__.py` and one module per command:
  `chiq_main.py`, `chiq_post.py`, `chiq_fft.py`, `gen_qpath.py`, `gen_allq.py`,
  `calc_Iq.py`, `calc_Iq_scl.py`, `plot_chiq_path.py`, `plot_Ir.py`,
  `eigenvec_viewer.py`, `dcore_chiq.py` — each exposing `def main()`. Plus
  `chiq/cli/_deprecated.py` (`.py`-command alias wrappers).
- Replace: `python/scripts/*.py` — each becomes a thin wrapper delegating to `chiq.cli.<name>:main`.
- Replace: `python/package/bse` symlink → real `python/package/bse/` package (forwarding modules + `point_group_data/` mirror).
- Modify: `python/package/chiq/g2scl_core.py` — lazy `matplotlib`.
- Tests: `tests/python/non-mpi/packaging/` (new dir) — `test_bse_shim.py`, `test_cli_package.py`, `test_matplotlib_lazy.py`.

---

### Task 1: Make matplotlib a lazy import in `chiq.g2scl_core`

**Files:**
- Modify: `python/package/chiq/g2scl_core.py:9-12` (module-scope matplotlib imports)
- Test: `tests/python/non-mpi/packaging/test_matplotlib_lazy.py`

**Interfaces:**
- Consumes: nothing.
- Produces: `chiq.g2scl_core` importable with no `matplotlib` installed.

- [ ] **Step 1: Write the failing test**

```python
# tests/python/non-mpi/packaging/test_matplotlib_lazy.py
import builtins
import importlib
import sys
import pytest

def test_g2scl_core_imports_without_matplotlib(monkeypatch):
    # Simulate matplotlib being absent: block its import, drop cached modules.
    for name in list(sys.modules):
        if name == "matplotlib" or name.startswith("matplotlib.") or name == "chiq.g2scl_core":
            monkeypatch.delitem(sys.modules, name, raising=False)
    real_import = builtins.__import__
    def blocked_import(name, *args, **kwargs):
        if name == "matplotlib" or name.startswith("matplotlib"):
            raise ImportError("matplotlib is blocked for this test")
        return real_import(name, *args, **kwargs)
    monkeypatch.setattr(builtins, "__import__", blocked_import)
    # Importing the module must NOT import matplotlib at module scope.
    importlib.import_module("chiq.g2scl_core")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH="$(pwd)/python/package:$(pwd)/build/src" python3 -m pytest tests/python/non-mpi/packaging/test_matplotlib_lazy.py -v`
Expected: FAIL with `ImportError: matplotlib is blocked` (module-scope import triggers it).

- [ ] **Step 3: Move the matplotlib imports into the functions that use them**

In `python/package/chiq/g2scl_core.py`, delete the module-scope lines (currently near lines 9-12):
```python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
```
Then, in each function that references `plt` or `cm`, add at the top of that function:
```python
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import cm
```
(Find the referencing functions with `grep -n "plt\.\|cm\." python/package/chiq/g2scl_core.py`; add the four lines inside each. If several plotting helpers share a call chain, importing inside the outermost plotting entry point is enough — but importing inside each function that names `plt`/`cm` is safe and simplest.)

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH="$(pwd)/python/package:$(pwd)/build/src" python3 -m pytest tests/python/non-mpi/packaging/test_matplotlib_lazy.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add python/package/chiq/g2scl_core.py tests/python/non-mpi/packaging/test_matplotlib_lazy.py
git commit -m "refactor(g2scl_core): import matplotlib lazily so it stays optional"
```

---

### Task 2: Replace the `bse` symlink with a static forwarding-shim package

**Files:**
- Delete: `python/package/bse` (symlink)
- Create: `python/package/bse/__init__.py`, `python/package/bse/<mod>.py` (×9 allowlisted modules), `python/package/bse/point_group_data/__init__.py` + `python/package/bse/point_group_data/<X>.py` (×8)
- Test: `tests/python/non-mpi/packaging/test_bse_shim.py`

**Interfaces:**
- Consumes: `chiq.*` modules (unchanged).
- Produces: `import bse`, `from bse import <mod>`, `import bse.point_group_data.<X>` work; `import bse` pulls no optional deps; a `FutureWarning` fires once on `import bse`.

- [ ] **Step 1: Write the failing test**

```python
# tests/python/non-mpi/packaging/test_bse_shim.py
import builtins
import importlib
import sys
import warnings
import pytest

def _fresh(monkeypatch, *drop_prefixes):
    for name in list(sys.modules):
        if any(name == p or name.startswith(p + ".") or name == p for p in drop_prefixes):
            monkeypatch.delitem(sys.modules, name, raising=False)

def test_import_bse_warns_and_pulls_no_optional_deps(monkeypatch):
    _fresh(monkeypatch, "bse")
    real_import = builtins.__import__
    def blocked(name, *a, **k):
        if name.split(".")[0] in ("mpi4py", "dcore", "matplotlib"):
            raise ImportError(f"{name} blocked")
        return real_import(name, *a, **k)
    monkeypatch.setattr(builtins, "__import__", blocked)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        importlib.import_module("bse")            # must not import optional deps
    assert any(issubclass(x.category, FutureWarning) for x in w)

def test_from_bse_import_public_submodule():
    from bse import h5bse
    import chiq.h5bse
    assert h5bse.h5BSE is chiq.h5bse.h5BSE       # class identity preserved

def test_nested_point_group_data():
    import bse.point_group_data.C1
    import chiq.point_group_data.C1
    # a public symbol is the same object
    assert dir(bse.point_group_data.C1)

def test_all_allowlisted_modules_importable():
    for mod in ["h5bse", "bse_toml", "matrix_dict", "point_group", "tools", "index_pair", "mpi"]:
        importlib.import_module(f"bse.{mod}")
    for x in ["base", "C1", "C2", "D3", "D4", "D6", "O", "Oh"]:
        importlib.import_module(f"bse.point_group_data.{x}")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH="$(pwd)/python/package:$(pwd)/build/src" python3 -m pytest tests/python/non-mpi/packaging/test_bse_shim.py -v`
Expected: FAIL — currently `bse` is a symlink to `chiq`, so `import bse` gives no warning and `bse.h5bse is chiq.h5bse`'s module differs / the warning assertion fails.

- [ ] **Step 3: Remove the symlink and create the real shim package**

```bash
git rm python/package/bse            # removes the symlink
mkdir -p python/package/bse/point_group_data
```

Create `python/package/bse/__init__.py`:
```python
"""Backward-compatibility shim for the old `bse` package name.

`bse` was renamed to `chiq`. This shim forwards a fixed allowlist of public
submodules to `chiq.*` (deterministic, IDE/type-checker friendly). It is
deprecated and will be REMOVED in ChiQ 2.0 -- import from `chiq` instead.
"""
import warnings

warnings.warn(
    "The 'bse' package is deprecated and will be removed in ChiQ 2.0; "
    "import from 'chiq' instead (e.g. `from chiq import h5bse`).",
    FutureWarning,
    stacklevel=2,
)
```

For each of the 9 allowlisted top-level modules, create `python/package/bse/<mod>.py` with exactly:
```python
from chiq.<mod> import *  # noqa: F401,F403  (backward-compat forwarding shim)
```
i.e. `bse/h5bse.py` → `from chiq.h5bse import *`, and likewise for `bse_toml`, `matrix_dict`, `point_group`, `sumk_dft_chi`, `g2scl_core`, `tools`, `index_pair`, `mpi`.

Create `python/package/bse/point_group_data/__init__.py`:
```python
from chiq.point_group_data import *  # noqa: F401,F403
```
and for each of `base, C1, C2, D3, D4, D6, O, Oh` create
`python/package/bse/point_group_data/<X>.py`:
```python
from chiq.point_group_data.<X> import *  # noqa: F401,F403
```

Helper to generate the forwarding files (run once, then delete the helper — or create by hand):
```bash
cd python/package/bse
for m in h5bse bse_toml matrix_dict point_group sumk_dft_chi g2scl_core tools index_pair mpi; do
  printf 'from chiq.%s import *  # noqa: F401,F403  (backward-compat forwarding shim)\n' "$m" > "$m.py"
done
printf 'from chiq.point_group_data import *  # noqa: F401,F403\n' > point_group_data/__init__.py
for x in base C1 C2 D3 D4 D6 O Oh; do
  printf 'from chiq.point_group_data.%s import *  # noqa: F401,F403\n' "$x" > "point_group_data/$x.py"
done
cd -
```

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH="$(pwd)/python/package:$(pwd)/build/src" python3 -m pytest tests/python/non-mpi/packaging/test_bse_shim.py -v`
Expected: PASS (4 tests). Note `test_from_bse_import_public_submodule` asserts *class* identity (`h5BSE`), which `import *` preserves.

- [ ] **Step 5: Commit**

```bash
git add python/package/bse tests/python/non-mpi/packaging/test_bse_shim.py
git commit -m "feat(bse): replace symlink with static forwarding-shim package (FutureWarning, allowlisted, lazy)"
```

---

### Task 3: Create the `chiq.cli` package (move CLI implementations)

**Files:**
- Create: `python/package/chiq/cli/__init__.py` (empty: `# chiq CLI entry-point modules`)
- Move: `python/scripts/<name>.py` → `python/package/chiq/cli/<name>.py` for the 11 commands
- Modify (in the moved files): sibling imports; `dcore_chiq` `run`→`main`; lazy `dcore`/`matplotlib` imports
- Test: `tests/python/non-mpi/packaging/test_cli_package.py`

**Interfaces:**
- Consumes: `chiq.*` (unchanged absolute imports inside the scripts).
- Produces: `chiq.cli.<name>:main` for each command (`chiq_main, chiq_post, chiq_fft, gen_qpath, gen_allq, calc_Iq, calc_Iq_scl, plot_chiq_path, plot_Ir, eigenvec_viewer, dcore_chiq`); `chiq.cli.gen_qpath.GenQPath` for intra-package use.

- [ ] **Step 1: Write the failing test**

```python
# tests/python/non-mpi/packaging/test_cli_package.py
import importlib
import subprocess
import sys
import pytest

CORE_COMMANDS = [
    "chiq_main", "chiq_post", "chiq_fft", "gen_qpath", "gen_allq",
    "calc_Iq", "calc_Iq_scl", "plot_chiq_path", "plot_Ir", "eigenvec_viewer",
]

@pytest.mark.parametrize("name", CORE_COMMANDS + ["dcore_chiq"])
def test_cli_module_has_main(name):
    mod = importlib.import_module(f"chiq.cli.{name}")
    assert callable(mod.main)

def test_gen_qpath_exports_GenQPath():
    from chiq.cli.gen_qpath import GenQPath
    assert GenQPath is not None

def test_chiq_main_version_runs_on_core():
    # --version must work without optional deps / without the built extension.
    # `python -c "<code>" --version` puts "--version" in sys.argv[1], which argparse reads.
    import os
    env = dict(os.environ)
    r = subprocess.run([sys.executable, "-c",
                        "from chiq.cli.chiq_main import main; main()", "--version"],
                       env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "ChiQ version" in (r.stdout + r.stderr)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH="$(pwd)/python/package:$(pwd)/build/src" python3 -m pytest tests/python/non-mpi/packaging/test_cli_package.py -v`
Expected: FAIL — `chiq.cli` does not exist.

- [ ] **Step 3: Move the command modules and fix imports**

```bash
mkdir -p python/package/chiq/cli
printf '# chiq CLI entry-point modules\n' > python/package/chiq/cli/__init__.py
for n in chiq_main chiq_post chiq_fft gen_qpath gen_allq calc_Iq calc_Iq_scl plot_chiq_path plot_Ir eigenvec_viewer dcore_chiq; do
  git mv "python/scripts/$n.py" "python/package/chiq/cli/$n.py"
done
```

Then edit the moved files:

1. **`chiq/cli/chiq_post.py`** — sibling import. Change:
   ```python
   import chiq_main # in the same directory
   ```
   to:
   ```python
   from chiq.cli import chiq_main
   ```

2. **`chiq/cli/chiq_fft.py`** and **`chiq/cli/gen_allq.py`** — sibling import. Change:
   ```python
   from gen_qpath import GenQPath # in the same directory
   ```
   to:
   ```python
   from chiq.cli.gen_qpath import GenQPath
   ```

3. **`chiq/cli/dcore_chiq.py`** — rename the entry function `run` → `main` (its `def run():` becomes `def main():`; update the `if __name__ == '__main__': run()` to `main()`), and make the DCore import lazy: move the module-scope
   ```python
   from dcore._dispatcher import HDFArchive, dyson
   from dcore.dmft_core import DMFTCoreSolver
   ```
   (and any other `dcore`/`mpi4py` module-scope imports) to the top of `main()` (or the worker function `main()` calls), wrapped so a missing extra is actionable:
   ```python
   def main():
       try:
           from dcore._dispatcher import HDFArchive, dyson
           from dcore.dmft_core import DMFTCoreSolver
       except ImportError as e:
           import sys
           sys.exit("dcore_chiq requires the 'dcore' extra: pip install chiq[dcore]  (%s)" % e)
       # ... rest of main, using the now-local imports ...
   ```
   Ensure `--help`/`--version` (argparse) are parsed **before** the try/except import, so `dcore_chiq --version` works on a core install.

4. **`chiq/cli/plot_chiq_path.py`** and **`chiq/cli/plot_Ir.py`** — move the module-scope
   ```python
   import matplotlib
   matplotlib.use("Agg")
   import matplotlib.pyplot as plt
   ```
   into `main()` (after arg parsing), with the same try/except pattern giving
   `"... requires the 'plot' extra: pip install chiq[plot]"` on `ImportError`.

5. Any remaining absolute `from chiq import ...` imports in the moved files stay unchanged.

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH="$(pwd)/python/package:$(pwd)/build/src" python3 -m pytest tests/python/non-mpi/packaging/test_cli_package.py -v`
Expected: PASS. If a `--version` test fails because a command lacks a `--version` flag, adjust that command's test parametrization to one that has it (`chiq_main` has `--version`); the key assertion is that the core command modules import and expose `main`.

- [ ] **Step 5: Commit**

```bash
git add python/package/chiq/cli tests/python/non-mpi/packaging/test_cli_package.py
git commit -m "feat(cli): move CLI implementations into importable chiq.cli package"
```

---

### Task 4: Legacy `python/scripts/*.py` wrappers + deprecated `.py` aliases

**Files:**
- Create: `python/scripts/<name>.py` (×11 thin wrappers — same names as before)
- Create: `python/package/chiq/cli/_deprecated.py` (`.py`-alias entry callables for Plan B)
- Test: `tests/python/non-mpi/packaging/test_cli_package.py` (extend); the existing integration tests are the real regression gate.

**Interfaces:**
- Consumes: `chiq.cli.<name>:main`.
- Produces: runnable `python/scripts/<name>.py` (for the legacy CMake `bin` install); `chiq.cli._deprecated.<name>_py` callables (Plan B entry points).

- [ ] **Step 1: Write the failing test**

```python
# append to tests/python/non-mpi/packaging/test_cli_package.py
import os

def test_scripts_are_thin_wrappers():
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
    p = os.path.join(root, "python", "scripts", "chiq_main.py")
    text = open(p).read()
    assert "from chiq.cli.chiq_main import main" in text

def test_deprecated_alias_prints_stderr(capfd):
    from chiq.cli import _deprecated
    assert hasattr(_deprecated, "chiq_main_py")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH="$(pwd)/python/package:$(pwd)/build/src" python3 -m pytest tests/python/non-mpi/packaging/test_cli_package.py -k "wrapper or deprecated" -v`
Expected: FAIL — the moved scripts no longer exist at `python/scripts/`, and `_deprecated` is absent.

- [ ] **Step 3: Create the wrappers and the `_deprecated` module**

For each of the 11 commands, create `python/scripts/<name>.py`:
```python
#!/usr/bin/env python3
from chiq.cli.<name> import main

if __name__ == "__main__":
    main()
```
Generate them:
```bash
for n in chiq_main chiq_post chiq_fft gen_qpath gen_allq calc_Iq calc_Iq_scl plot_chiq_path plot_Ir eigenvec_viewer dcore_chiq; do
  printf '#!/usr/bin/env python3\nfrom chiq.cli.%s import main\n\nif __name__ == "__main__":\n    main()\n' "$n" > "python/scripts/$n.py"
  chmod +x "python/scripts/$n.py"
done
```

Create `python/package/chiq/cli/_deprecated.py`:
```python
"""Deprecated `<name>.py` console-script aliases.

Each callable prints a stderr deprecation notice (so shell users see it even
though DeprecationWarning is filtered) then delegates to the real entry point.
The `.py` command names are removed in ChiQ 2.0.
"""
import sys
import importlib


def _make(name):
    def _entry():
        sys.stderr.write(
            "warning: the '%s.py' command is deprecated and will be removed in "
            "ChiQ 2.0; use '%s' instead.\n" % (name, name)
        )
        importlib.import_module("chiq.cli." + name).main()
    _entry.__name__ = name + "_py"
    return _entry


chiq_main_py = _make("chiq_main")
chiq_post_py = _make("chiq_post")
chiq_fft_py = _make("chiq_fft")
gen_qpath_py = _make("gen_qpath")
gen_allq_py = _make("gen_allq")
calc_Iq_py = _make("calc_Iq")
calc_Iq_scl_py = _make("calc_Iq_scl")
plot_chiq_path_py = _make("plot_chiq_path")
plot_Ir_py = _make("plot_Ir")
eigenvec_viewer_py = _make("eigenvec_viewer")
dcore_chiq_py = _make("dcore_chiq")
```

- [ ] **Step 4: Run tests to verify they pass**

Run the packaging tests:
`PYTHONPATH="$(pwd)/python/package:$(pwd)/build/src" python3 -m pytest tests/python/non-mpi/packaging/ -v`
Expected: PASS.
Run the **full non-MPI suite** to prove no regression (scripts still run by name):
```bash
chmod +x python/scripts/*.py
export PYTHONPATH="$(pwd)/python/package:$(pwd)/build/src"
export PATH="$(pwd)/python/scripts:$PATH"
cd build && python3 -m pytest ../tests/python/non-mpi -q
```
Expected: all previously-passing tests still pass (the `os.system("chiq_main.py ...")` calls now run the wrappers, which import `chiq.cli.chiq_main`).

- [ ] **Step 5: Commit**

```bash
git add python/scripts python/package/chiq/cli/_deprecated.py tests/python/non-mpi/packaging/test_cli_package.py
git commit -m "feat(cli): legacy script wrappers + deprecated .py-alias entry callables"
```

---

## Final verification

- [ ] Full non-MPI suite green (integration + new packaging tests):
  `export PYTHONPATH="$(pwd)/python/package:$(pwd)/build/src"; export PATH="$(pwd)/python/scripts:$PATH"; cd build && python3 -m pytest ../tests/python/non-mpi -q`
- [ ] `import bse` emits a `FutureWarning` and pulls no optional deps; `from bse import h5bse` works.
- [ ] `python3 -c "import chiq.cli.chiq_main as m; m.main"` works; `chiq_main.py --version` (wrapper) prints the version.
- [ ] `chiq.g2scl_core` imports with matplotlib blocked.

## Out of scope (Plan B — build & packaging)

`pyproject.toml` (scikit-build-core, deps, `[project.scripts]` entry points incl. the
`chiq.cli._deprecated.*_py` aliases, extras), the `_bse_solver` pybind/CMake rename +
`chiq/solver/cpp.py` import change + top-level `bse_solver` shim package, `CHIQ_WHEEL_BUILD`
CMake option + `find_package(pybind11)`, the legacy `python/CMakeLists.txt` install-rule
updates (install `bse`/`bse_solver` shims; `file(REMOVE)` the old `bse` symlink), wheel/sdist/
editable clean-venv install tests, CI job, and `README.md`/`doc/install.rst` updates.
