#!/bin/sh
set -eu

# POSIX bootstrap: /bin/sh does not read Bash's BASH_ENV startup hook. No
# command lookup or Python execution is allowed before the empty environment.
if [ "${CHIQ_LEGACY_VERIFIER_INTERNAL:-0}" != 1 ]; then
  exec /usr/bin/env -i \
    CHIQ_LEGACY_VERIFIER_INTERNAL=1 \
    "HOME=${HOME:-}" \
    "PATH=/usr/bin:/bin:/usr/sbin:/sbin:/opt/homebrew/bin:/usr/local/bin" \
    "TMPDIR=${TMPDIR:-/tmp}" \
    "PYTHON_REQUEST=${PYTHON:-python3}" \
    "LANG=${LANG:-}" \
    "LC_ALL=${LC_ALL:-}" \
    "LC_CTYPE=${LC_CTYPE:-}" \
    "CC=${CC:-}" \
    "CXX=${CXX:-}" \
    "SDKROOT=${SDKROOT:-${SDK:-}}" \
    "EIGEN3_INCLUDE_DIR=${EIGEN3_INCLUDE_DIR:-}" \
    "PYTHON_INCLUDE_DIR=${PYTHON_INCLUDE_DIR:-}" \
    "PYTHON_LIBRARY=${PYTHON_LIBRARY:-}" \
    "JOBS=${JOBS:-}" \
    /bin/bash "$0" "$@"
  exit 127
fi

if [ -z "${BASH_VERSION:-}" ]; then
  echo "internal legacy verifier requires Bash" >&2
  exit 2
fi
set -euo pipefail
SAFE_PATH="/usr/bin:/bin:/usr/sbin:/sbin:/opt/homebrew/bin:/usr/local/bin"

# Fail closed if anyone bypasses the POSIX bootstrap. Bash itself exports
# PWD/SHLVL/_, while every other exported name must be explicitly allowlisted.
audit_error=""
while IFS= read -r name; do
  case "$name" in
    CHIQ_LEGACY_VERIFIER_INTERNAL|HOME|PATH|TMPDIR|PYTHON_REQUEST|LANG|LC_ALL|LC_CTYPE|CC|CXX|SDKROOT|EIGEN3_INCLUDE_DIR|PYTHON_INCLUDE_DIR|PYTHON_LIBRARY|JOBS|PWD|SHLVL|_)
      ;;
    PYTHON*|PYTEST_*|PYTHONPATH|PYTHONHOME|VIRTUAL_ENV|CONDA_PREFIX|PIP_*|GIT_*|BASH_ENV|ENV|CDPATH|SHELLOPTS|BASHOPTS|GLOBIGNORE)
      audit_error="forbidden exported variable $name"
      break
      ;;
    *)
      audit_error="unexpected exported variable $name"
      break
      ;;
  esac
done < <(compgen -e)
if [[ ${CHIQ_LEGACY_VERIFIER_INTERNAL:-0} != 1 || "$PATH" != "$SAFE_PATH" || -n "$audit_error" ]]; then
  echo "unsafe internal legacy verifier environment: ${audit_error:-invalid marker or PATH}" >&2
  exit 2
fi

case "$PYTHON_REQUEST" in
  /*) PYTHON_CANDIDATE=$PYTHON_REQUEST ;;
  */*)
    echo "unsafe internal legacy verifier environment: PYTHON must be absolute or a bare command" >&2
    exit 2
    ;;
  *)
    PYTHON_CANDIDATE=$(command -v -- "$PYTHON_REQUEST") || {
      echo "selected Python is unavailable on the safe PATH: $PYTHON_REQUEST" >&2
      exit 2
    }
    ;;
esac
[[ -x "$PYTHON_CANDIDATE" ]] || {
  echo "selected Python is not executable: $PYTHON_CANDIDATE" >&2
  exit 2
}
PYTHON="$($PYTHON_CANDIDATE -I -c 'import os, sys; print(os.path.realpath(sys.executable))')"
PYTHON_BIN=${PYTHON%/*}
PATH="$PYTHON_BIN:$SAFE_PATH"
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1
export PATH PYTEST_DISABLE_PLUGIN_AUTOLOAD

REPO="$(cd "$(/usr/bin/dirname "$0")/../../../.." && pwd)"
TMPDIR="$($PYTHON -I - "$TMPDIR" "$REPO" <<'PY'
import os
import sys

requested = os.path.realpath(sys.argv[1])
repository = os.path.realpath(sys.argv[2])
try:
    inside = os.path.commonpath((requested, repository)) == repository
except ValueError:
    inside = True
print(os.path.realpath("/tmp") if inside else requested)
PY
)"
export TMPDIR

print_environment_policy() {
  echo "Environment policy: env -i explicit allowlist"
  echo "Environment audit OK"
  echo "Restored names: HOME LANG LC_ALL LC_CTYPE PATH TMPDIR PYTHON CC CXX SDKROOT EIGEN3_INCLUDE_DIR PYTHON_INCLUDE_DIR PYTHON_LIBRARY JOBS"
  echo "Controlled names: PYTEST_DISABLE_PLUGIN_AUTOLOAD=1"
  echo "Scrubbed names: PYTHON* PYTEST_* VIRTUAL_ENV CONDA_PREFIX PIP_* GIT_* BASH_ENV ENV CDPATH SHELLOPTS BASHOPTS GLOBIGNORE"
  for name in HOME LANG LC_ALL LC_CTYPE PATH TMPDIR PYTHON CC CXX SDKROOT EIGEN3_INCLUDE_DIR PYTHON_INCLUDE_DIR PYTHON_LIBRARY JOBS PYTEST_DISABLE_PLUGIN_AUTOLOAD; do
    if [[ -n ${!name:-} ]]; then
      printf '%s=%s\n' "$name" "${!name}"
    fi
  done
}

if [[ ${1:-} == "--audit-environment-only" ]]; then
  print_environment_policy
  exit 0
fi

# Authoritative verification for the legacy CHIQ_WHEEL_BUILD=OFF layout.
WORK_ROOT=$(mktemp -d "${TMPDIR:-/tmp}/chiq-legacy-verification.XXXXXX")
WORK_ROOT=$(cd "$WORK_ROOT" && pwd -P)
BUILD="$WORK_ROOT/build"
PREFIX="$WORK_ROOT/prefix"
EXTERNAL_CWD="$WORK_ROOT/external-cwd"
EMPTY_PYTEST_CONFIG="$WORK_ROOT/pytest.ini"
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
  SDK=${SDKROOT:-/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk}
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

  PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 \
    "$PYTHON" -m pytest -c "$EMPTY_PYTEST_CONFIG" "$BUILD/tests/python/non-mpi" -q \
    --ignore="$BUILD/tests/python/non-mpi/packaging"
  PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 \
    "$PYTHON" -m pytest -c "$EMPTY_PYTEST_CONFIG" \
    "$REPO/tests/python/non-mpi/packaging" -q
)

echo "LEGACY INSTALL VERIFIED"
echo "Diagnostics retained at $WORK_ROOT"
