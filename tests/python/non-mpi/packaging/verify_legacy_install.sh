#!/usr/bin/env bash
# Verifies the legacy (CHIQ_WHEEL_BUILD=OFF) install layout into a temp prefix.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/../../../.." && pwd)"
PREFIX=$(mktemp -d)
SDK=${SDK:-/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk}
PYTHON=${PYTHON:-$(command -v python3)}
PYTHON_INCLUDE_DIR=${PYTHON_INCLUDE_DIR:-$("$PYTHON" -c "import sysconfig; print(sysconfig.get_path('include'))")}
PYTHON_LIBRARY=${PYTHON_LIBRARY:-$("$PYTHON" -c "import os, sysconfig; print(os.path.join(sysconfig.get_config_var('LIBDIR'), sysconfig.get_config_var('LIBRARY') or sysconfig.get_config_var('LDLIBRARY')))")}

if [[ -z "${EIGEN3_INCLUDE_DIR:-}" ]]; then
  if brew --prefix eigen@3 >/dev/null 2>&1; then
    EIGEN3_INCLUDE_DIR="$(brew --prefix eigen@3)/include/eigen3"
  else
    EIGEN3_INCLUDE_DIR="$(brew --prefix eigen)/include/eigen3"
  fi
fi

cd "$REPO"
rm -rf build_legacy
mkdir build_legacy
cd build_legacy
cmake .. -DTesting=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DEIGEN3_INCLUDE_DIR="$EIGEN3_INCLUDE_DIR" \
  -DPYTHON_EXECUTABLE="$PYTHON" \
  -DPYTHON_INCLUDE_DIR="$PYTHON_INCLUDE_DIR" \
  -DPYTHON_LIBRARY="$PYTHON_LIBRARY" \
  -DCMAKE_CXX_FLAGS="-isysroot $SDK -isystem $SDK/usr/include/c++/v1"
make _bse_solver
make install
cd "$REPO"

LIBDIR=$(find "$PREFIX" -type d -name bse-python)
test -f "$LIBDIR/chiq/__init__.py"
test -f "$LIBDIR"/chiq/_bse_solver*.so
test -f "$LIBDIR/bse/__init__.py" && ! test -L "$LIBDIR/bse"
test -f "$LIBDIR/bse_solver/__init__.py"
PYTHONPATH="$LIBDIR" "$PYTHON" -c "from chiq import _bse_solver; import bse, bse_solver; print('legacy layout OK')"
echo "LEGACY INSTALL VERIFIED at $PREFIX"
