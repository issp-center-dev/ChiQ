#!/usr/bin/env bash
# Verifies the legacy (CHIQ_WHEEL_BUILD=OFF) install layout into a temp prefix.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/../../../.." && pwd)"
PREFIX=$(mktemp -d)
PYTHON=${PYTHON:-$(command -v python3)}
PYTHON_INCLUDE_DIR=${PYTHON_INCLUDE_DIR:-$("$PYTHON" -c "import sysconfig; print(sysconfig.get_path('include'))")}
PYTHON_LIBRARY=${PYTHON_LIBRARY:-$("$PYTHON" -c "import os, sysconfig; print(os.path.join(sysconfig.get_config_var('LIBDIR'), sysconfig.get_config_var('LIBRARY') or sysconfig.get_config_var('LDLIBRARY')))")}
CMAKE_PLATFORM_ARGS=()

if [[ -z "${EIGEN3_INCLUDE_DIR:-}" ]]; then
  if [[ "$(uname -s)" == "Darwin" ]]; then
    if brew --prefix eigen@3 >/dev/null 2>&1; then
      EIGEN3_INCLUDE_DIR="$(brew --prefix eigen@3)/include/eigen3"
    else
      EIGEN3_INCLUDE_DIR="$(brew --prefix eigen)/include/eigen3"
    fi
  else
    EIGEN3_INCLUDE_DIR=/usr/include/eigen3
  fi
fi

if [[ "$(uname -s)" == "Darwin" ]]; then
  SDK=${SDK:-/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk}
  CMAKE_PLATFORM_ARGS+=(
    "-DCMAKE_CXX_FLAGS=-isysroot $SDK -isystem $SDK/usr/include/c++/v1"
  )
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
  "${CMAKE_PLATFORM_ARGS[@]}"
make _bse_solver

# The uninstalled build tree must work with the environment exported by
# chiqvars.sh. Use a clean package copy so a developer's local .so cannot mask
# an incorrect build-tree import path.
SOURCE_PACKAGE=$(mktemp -d)
cp -R "$REPO/python/package" "$SOURCE_PACKAGE/package"
find "$SOURCE_PACKAGE/package" -name '_bse_solver*.so' -delete
PYTHONPATH="$SOURCE_PACKAGE/package:$REPO/build_legacy/src" "$PYTHON" -c \
  "from chiq.solver import get_solver; get_solver('cpp', 1.0, [1], [1], [1]); print('legacy build tree OK')"

make install
cd "$REPO"

LIBDIR=$(find "$PREFIX" -type d -name bse-python)
test -f "$LIBDIR/chiq/__init__.py"
test -f "$LIBDIR"/chiq/_bse_solver*.so
test -f "$LIBDIR/bse/__init__.py" && ! test -L "$LIBDIR/bse"
test -f "$LIBDIR/bse_solver/__init__.py"
test -x "$PREFIX/bin/chiq_main"
PYTHONPATH="$LIBDIR" "$PREFIX/bin/chiq_main" --version
PYTHONPATH="$LIBDIR" "$PYTHON" -c "from chiq import _bse_solver; import bse, bse_solver; print('legacy layout OK')"
echo "LEGACY INSTALL VERIFIED at $PREFIX"
