#!/usr/bin/env bash
set -euo pipefail

REPO="$(cd "$(dirname "$0")/../../../.." && pwd)"
PYTHON=${PYTHON:-$(command -v python3)}
if [[ "$(uname -s)" == "Darwin" ]]; then
  SDK=${SDK:-/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk}
  if brew --prefix eigen@3 >/dev/null 2>&1; then
    EIGEN3_INCLUDE_DIR="$(brew --prefix eigen@3)/include/eigen3"
  else
    EIGEN3_INCLUDE_DIR="$(brew --prefix eigen)/include/eigen3"
  fi
else
  EIGEN3_INCLUDE_DIR=${EIGEN3_INCLUDE_DIR:-/usr/include/eigen3}
fi

VENV=$(mktemp -d)/venv
"$PYTHON" -m venv "$VENV"
# shellcheck disable=SC1091
source "$VENV/bin/activate"
pip install -U pip build

PYTHON_INCLUDE_DIR=$(python -c "import sysconfig; print(sysconfig.get_path('include'))")
PYTHON_LIBRARY=$(python -c "import os, sysconfig; print(os.path.join(sysconfig.get_config_var('LIBDIR'), sysconfig.get_config_var('LIBRARY') or sysconfig.get_config_var('LDLIBRARY')))")
export CMAKE_ARGS="-DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DEIGEN3_INCLUDE_DIR=$EIGEN3_INCLUDE_DIR \
  -DPYTHON_INCLUDE_DIR=$PYTHON_INCLUDE_DIR \
  -DPYTHON_LIBRARY=$PYTHON_LIBRARY"
if [[ "$(uname -s)" == "Darwin" ]]; then
  CMAKE_ARGS="$CMAKE_ARGS -DCMAKE_CXX_FLAGS=-isysroot\\ $SDK\\ -isystem\\ $SDK/usr/include/c++/v1"
  export CMAKE_ARGS
fi

# --- wheel: build, install into the clean venv OUTSIDE the repo, import-check ---
WORK=$(mktemp -d)
cp -R "$REPO" "$WORK/ChiQ"
cd "$WORK/ChiQ"
git submodule update --init >/dev/null 2>&1 || true
pip install .
cd /tmp
python -c "import chiq, chiq._bse_solver; print('wheel import chiq OK')"
python -W ignore::FutureWarning -c "import bse, bse_solver; print('wheel shims OK')"
python -c "import sys; import chiq.g2scl_core; assert 'matplotlib' not in sys.modules; print('matplotlib optional OK')"
chiq_main --version
python - <<'PY'
import subprocess

r = subprocess.run(["chiq_main.py", "--version"], capture_output=True, text=True)
assert "deprecated" in r.stderr.lower(), r.stderr
print("deprecated alias notice OK")
PY
SITE=$(python -c "import chiq, os; print(os.path.dirname(os.path.dirname(chiq.__file__)))")
! find "$SITE/chiq" -name '*.h5' | grep -q . && echo "no test .h5 in wheel OK"
pip uninstall -y chiq

# --- sdist: build, inspect archive contents, install from it ---
cd "$WORK/ChiQ"
python -m build --sdist --outdir dist_sdist
SDIST=$(ls dist_sdist/*.tar.gz | head -1)
tar tzf "$SDIST" | grep -q 'src/bse_solver_pybind.cpp' && echo "sdist has C++ sources OK"
tar tzf "$SDIST" | grep -q 'tests/.*\.h5' && echo "sdist has test .h5 OK"
pip install "$SDIST"
cd /tmp
python -c "import chiq, chiq._bse_solver; print('sdist install OK')"
pip uninstall -y chiq

# --- editable ---
cd "$WORK/ChiQ"
pip install -e .
cd /tmp
python -c "import chiq, chiq._bse_solver; print('editable install OK')"
deactivate
echo "ALL PIP INSTALL PATHS VERIFIED"
