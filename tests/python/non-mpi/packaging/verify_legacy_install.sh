#!/usr/bin/env bash
# Authoritative verification for the legacy CHIQ_WHEEL_BUILD=OFF layout.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/../../../.." && pwd)"
WORK_ROOT=$(mktemp -d "${TMPDIR:-/tmp}/chiq-legacy-verification.XXXXXX")
WORK_ROOT=$(cd "$WORK_ROOT" && pwd -P)
BUILD="$WORK_ROOT/build"
PREFIX="$WORK_ROOT/prefix"
EXTERNAL_CWD="$WORK_ROOT/external-cwd"
EMPTY_PYTEST_CONFIG="$WORK_ROOT/pytest.ini"
PYTHON=${PYTHON:-$(command -v python3)}
JOBS=${JOBS:-2}
PYTHON_INCLUDE_DIR=${PYTHON_INCLUDE_DIR:-$("$PYTHON" -c "import sysconfig; print(sysconfig.get_path('include'))")}
PYTHON_LIBRARY=${PYTHON_LIBRARY:-$("$PYTHON" -c "import os, sysconfig; print(os.path.join(sysconfig.get_config_var('LIBDIR'), sysconfig.get_config_var('LIBRARY') or sysconfig.get_config_var('LDLIBRARY')))" )}
CMAKE_PLATFORM_ARGS=()

test -f "$REPO/extern/pybind11/include/pybind11/pybind11.h"
PYBIND11_TAG=$(git -C "$REPO/extern/pybind11" describe --tags --exact-match HEAD)
if [[ "$PYBIND11_TAG" != "v2.13.6" ]]; then
  echo "Expected populated pybind11 v2.13.6, found $PYBIND11_TAG" >&2
  exit 1
fi

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

mkdir -p "$BUILD" "$PREFIX" "$EXTERNAL_CWD"
printf '[pytest]\n' >"$EMPTY_PYTEST_CONFIG"
cmake -S "$REPO" -B "$BUILD" \
  -DCHIQ_WHEEL_BUILD=OFF \
  -DTesting=ON \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DCMAKE_INSTALL_LIBDIR=lib \
  -DEIGEN3_INCLUDE_DIR="$EIGEN3_INCLUDE_DIR" \
  -DPYTHON_EXECUTABLE="$PYTHON" \
  -DPYTHON_INCLUDE_DIR="$PYTHON_INCLUDE_DIR" \
  -DPYTHON_LIBRARY="$PYTHON_LIBRARY" \
  "${CMAKE_PLATFORM_ARGS[@]}"

(
  cd "$BUILD"
  make -j"$JOBS"
  ctest --output-on-failure
)

# Exercise both narrow native cleanup and exact stale-bse-symlink replacement.
LEGACY_ROOT="$PREFIX/lib/bse-python"
mkdir -p "$LEGACY_ROOT/chiq" "$WORK_ROOT/old-bse-target"
printf 'obsolete\n' >"$LEGACY_ROOT/bse_solver-old.so"
printf 'obsolete\n' >"$LEGACY_ROOT/chiq/_bse_solver-old.so"
printf 'keep\n' >"$LEGACY_ROOT/native-neighbor.so"
printf 'keep\n' >"$WORK_ROOT/old-bse-target/sentinel.txt"
ln -s "$WORK_ROOT/old-bse-target" "$LEGACY_ROOT/bse"

(
  cd "$BUILD"
  make install
)

test ! -e "$LEGACY_ROOT/bse_solver-old.so"
test ! -e "$LEGACY_ROOT/chiq/_bse_solver-old.so"
test "$(cat "$LEGACY_ROOT/native-neighbor.so")" = "keep"
test -d "$LEGACY_ROOT/bse" && ! test -L "$LEGACY_ROOT/bse"
test "$(cat "$WORK_ROOT/old-bse-target/sentinel.txt")" = "keep"
test -f "$LEGACY_ROOT/chiq/__init__.py"
test -f "$LEGACY_ROOT/bse/__init__.py"
test -f "$LEGACY_ROOT/bse_solver/__init__.py"

NATIVE_COUNT=$(find "$LEGACY_ROOT" -type f \( \
  -name '_bse_solver*.so' -o -name '_bse_solver*.pyd' -o \
  -name '_bse_solver*.dylib' \) | wc -l | tr -d ' ')
test "$NATIVE_COUNT" = 1

(
  cd "$EXTERNAL_CWD"
  export PYTHONPATH="$LEGACY_ROOT"
  export PATH="$PREFIX/bin:$PATH"
  "$PYTHON" "$REPO/tests/python/non-mpi/packaging/runtime_smoke.py" \
    --mode wheel --expected-root "$PREFIX"

  "$PYTHON" -m pytest -c "$EMPTY_PYTEST_CONFIG" "$BUILD/tests/python/non-mpi" -q \
    --ignore="$BUILD/tests/python/non-mpi/packaging"
  "$PYTHON" -m pytest -c "$EMPTY_PYTEST_CONFIG" \
    "$REPO/tests/python/non-mpi/packaging" -q
)

echo "LEGACY INSTALL VERIFIED"
echo "Diagnostics retained at $WORK_ROOT"
