#!/bin/sh
set -eu

# POSIX bootstrap: /bin/sh does not read Bash's BASH_ENV startup hook. This
# branch performs no command lookup and starts no Python before env -i.
if [ "${CHIQ_VERIFIER_INTERNAL:-0}" != 1 ]; then
  exec /usr/bin/env -i \
    CHIQ_VERIFIER_INTERNAL=1 \
    "PATH=/usr/bin:/bin:/usr/sbin:/sbin:/opt/homebrew/bin:/usr/local/bin" \
    "TMPDIR=${TMPDIR:-/tmp}" \
    "PYTHON_REQUEST=${PYTHON:-python3}" \
    "LANG=${LANG:-}" \
    "LC_ALL=${LC_ALL:-}" \
    "LC_CTYPE=${LC_CTYPE:-}" \
    "CC=${CC:-}" \
    "CXX=${CXX:-}" \
    "SDKROOT=${SDKROOT:-}" \
    "CMAKE_ARGS=${CMAKE_ARGS:-}" \
    "EIGEN3_INCLUDE_DIR=${EIGEN3_INCLUDE_DIR:-}" \
    /bin/bash "$0" "$@"
  exit 127
fi

# Everything below is parsed and executed only by the clean Bash child.
if [ -z "${BASH_VERSION:-}" ]; then
  echo "internal verifier requires Bash" >&2
  exit 2
fi
set -euo pipefail
SAFE_PATH="/usr/bin:/bin:/usr/sbin:/sbin:/opt/homebrew/bin:/usr/local/bin"

# Direct access to the internal branch fails closed under any ambient state.
# Bash itself may export PWD/SHLVL/_, but every other name is explicit.
audit_error=""
while IFS= read -r name; do
  case "$name" in
    CHIQ_VERIFIER_INTERNAL|PATH|TMPDIR|PYTHON_REQUEST|LANG|LC_ALL|LC_CTYPE|CC|CXX|SDKROOT|CMAKE_ARGS|EIGEN3_INCLUDE_DIR|PWD|SHLVL|_)
      ;;
    PYTHONPATH|PYTHONHOME|VIRTUAL_ENV|CONDA_PREFIX|CMAKE_PREFIX_PATH|CMAKE_BUILD_PARALLEL_LEVEL|SKBUILD_*|PIP_*|GIT_*|BASH_ENV|ENV|CDPATH|SHELLOPTS|BASHOPTS|GLOBIGNORE)
      audit_error="forbidden exported variable $name"
      break
      ;;
    *)
      audit_error="unexpected exported variable $name"
      break
      ;;
  esac
done < <(compgen -e)
if [[ ${CHIQ_VERIFIER_INTERNAL:-0} != 1 || "$PATH" != "$SAFE_PATH" || -n "$audit_error" ]]; then
  echo "unsafe internal verifier environment: ${audit_error:-invalid marker or PATH}" >&2
  exit 2
fi

case "$PYTHON_REQUEST" in
  /*) PYTHON_CANDIDATE=$PYTHON_REQUEST ;;
  */*)
    echo "unsafe internal verifier environment: PYTHON must be absolute or a bare command" >&2
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

REPO="$(cd "$(/usr/bin/dirname "$0")/../../../.." && pwd)"
PACKAGING="$REPO/tests/python/non-mpi/packaging"
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

POLICY_ROOT=$(mktemp -d "$TMPDIR/chiq-verifier-policy.XXXXXX")
HOME="$POLICY_ROOT/home"
mkdir "$HOME"
chmod 700 "$POLICY_ROOT" "$HOME"
PIP_CONFIG_FILE=/dev/null
GIT_CONFIG_GLOBAL=/dev/null
GIT_CONFIG_NOSYSTEM=1
CMAKE_FIND_USE_PACKAGE_REGISTRY=FALSE
export HOME PIP_CONFIG_FILE GIT_CONFIG_GLOBAL GIT_CONFIG_NOSYSTEM
export CMAKE_FIND_USE_PACKAGE_REGISTRY

print_environment_policy() {
  echo "Environment policy: env -i explicit allowlist"
  echo "Environment audit OK"
  echo "Restored names: LANG LC_ALL LC_CTYPE PATH TMPDIR PYTHON CC CXX SDKROOT CMAKE_ARGS EIGEN3_INCLUDE_DIR"
  echo "Controlled names: HOME(private) PIP_CONFIG_FILE GIT_CONFIG_GLOBAL GIT_CONFIG_NOSYSTEM CMAKE_FIND_USE_PACKAGE_REGISTRY"
  echo "Scrubbed names: PYTHONPATH PYTHONHOME VIRTUAL_ENV CONDA_PREFIX CMAKE_PREFIX_PATH CMAKE_BUILD_PARALLEL_LEVEL SKBUILD_* PIP_* GIT_* BASH_ENV ENV CDPATH SHELLOPTS BASHOPTS GLOBIGNORE"
  for name in HOME PIP_CONFIG_FILE GIT_CONFIG_GLOBAL GIT_CONFIG_NOSYSTEM CMAKE_FIND_USE_PACKAGE_REGISTRY LANG LC_ALL LC_CTYPE PATH TMPDIR PYTHON CC CXX SDKROOT CMAKE_ARGS EIGEN3_INCLUDE_DIR; do
    if [[ -n ${!name:-} ]]; then
      printf '%s=%q\n' "$name" "${!name}"
    fi
  done
}

if [[ ${1:-} == "--audit-environment-only" ]]; then
  print_environment_policy
  exit 0
fi

if [[ -z ${EIGEN3_INCLUDE_DIR:-} ]]; then
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
if [[ ${CMAKE_ARGS:-} != *EIGEN3_INCLUDE_DIR* ]]; then
  CMAKE_ARGS="${CMAKE_ARGS:-} -DEIGEN3_INCLUDE_DIR=$EIGEN3_INCLUDE_DIR"
fi
if [[ ${CMAKE_ARGS:-} != *CMAKE_FIND_USE_PACKAGE_REGISTRY* ]]; then
  CMAKE_ARGS="$CMAKE_ARGS -DCMAKE_FIND_USE_PACKAGE_REGISTRY=FALSE -DCMAKE_FIND_PACKAGE_NO_PACKAGE_REGISTRY=TRUE"
fi
export CMAKE_ARGS EIGEN3_INCLUDE_DIR

DIAGNOSTICS="$(mktemp -d "${TMPDIR:-/tmp}/chiq-pip-verification.XXXXXX")"
trap 'status=$?; echo "Diagnostics retained at $DIAGNOSTICS"; exit "$status"' EXIT

SNAPSHOT="$DIAGNOSTICS/source-snapshot"
SNAPSHOT_MANIFEST="$DIAGNOSTICS/snapshot-manifest.json"
ARTIFACTS="$DIAGNOSTICS/artifacts"
VALIDATED_ARTIFACTS="$DIAGNOSTICS/validated-artifacts"
LOGS="$DIAGNOSTICS/logs"
EXTRACT="$DIAGNOSTICS/extracted-sdist"
EXTERNAL_CWD="$DIAGNOSTICS/external-runtime-cwd"
FRONTEND_VENV="$DIAGNOSTICS/frontend-venv"
WHEEL_VENV="$DIAGNOSTICS/wheel-venv-A"
SDIST_VENV="$DIAGNOSTICS/sdist-venv-B"
EDITABLE_VENV="$DIAGNOSTICS/editable-venv-C"

mkdir -m 700 "$VALIDATED_ARTIFACTS"
mkdir -p "$ARTIFACTS/sdist" "$ARTIFACTS/wheel" "$LOGS" "$EXTERNAL_CWD"
for mode in frontend sdist-build wheel-build wheel sdist editable; do
  mkdir -p "$DIAGNOSTICS/tmp/$mode" "$DIAGNOSTICS/cache/$mode"
done

print_environment_policy >"$DIAGNOSTICS/sanitized-environment.txt"
cat "$DIAGNOSTICS/sanitized-environment.txt"

"$PYTHON" "$PACKAGING/source_snapshot.py" \
  "$REPO" "$SNAPSHOT" --manifest "$SNAPSHOT_MANIFEST" \
  >"$LOGS/source-snapshot.stdout.log" 2>"$LOGS/source-snapshot.stderr.log"

"$PYTHON" -m venv "$FRONTEND_VENV"
if "$FRONTEND_VENV/bin/python" -c 'import sys; raise SystemExit(0 if sys.version_info[:2] == (3, 8) else 1)'; then
  FRONTEND_REQUIREMENTS=("pip<25.1" "build<1.3" packaging toml)
  BACKEND_REQUIREMENT="scikit-build-core>=0.10,<0.11"
else
  FRONTEND_REQUIREMENTS=(pip build packaging toml)
  BACKEND_REQUIREMENT="scikit-build-core>=0.10,<1"
fi
TMPDIR="$DIAGNOSTICS/tmp/frontend" PIP_CACHE_DIR="$DIAGNOSTICS/cache/frontend" \
  "$FRONTEND_VENV/bin/python" -m pip install --no-cache-dir --upgrade \
  "${FRONTEND_REQUIREMENTS[@]}" \
  >"$LOGS/frontend-install.log" 2>&1
TMPDIR="$DIAGNOSTICS/tmp/frontend" PIP_CACHE_DIR="$DIAGNOSTICS/cache/frontend" \
  "$FRONTEND_VENV/bin/python" -m pip install --no-cache-dir \
  "$BACKEND_REQUIREMENT" 'pybind11>=2.12,<3' 'cmake>=3.15' \
  >"$LOGS/backend-probe-install.log" 2>&1
{
  "$FRONTEND_VENV/bin/python" --version
  "$FRONTEND_VENV/bin/python" -m pip --version
  "$FRONTEND_VENV/bin/python" -m build --version
} >"$DIAGNOSTICS/frontend-versions.txt" 2>&1
"$FRONTEND_VENV/bin/python" -m pip list --format=json \
  >"$DIAGNOSTICS/frontend-pip-list.json"
"$FRONTEND_VENV/bin/python" - "$DIAGNOSTICS/frontend-pip-list.json" \
  "$DIAGNOSTICS/backend-distributions.json" <<'PY'
# BACKEND_PROBE_PARSER_BEGIN
import json
import sys
from packaging.utils import canonicalize_name

with open(sys.argv[1], encoding="utf-8") as stream:
    installed = json.load(stream)
wanted = {canonicalize_name(name) for name in ("scikit-build-core", "pybind11", "cmake")}
selected = [
    row for row in installed if canonicalize_name(row["name"]) in wanted
]
if {canonicalize_name(row["name"]) for row in selected} != wanted:
    raise SystemExit("build-system distribution probe is incomplete")
with open(sys.argv[2], "w", encoding="utf-8") as stream:
    json.dump(selected, stream, indent=2, sort_keys=True)
    stream.write("\n")
# BACKEND_PROBE_PARSER_END
PY

TMPDIR="$DIAGNOSTICS/tmp/sdist-build" PIP_CACHE_DIR="$DIAGNOSTICS/cache/sdist-build" \
  "$FRONTEND_VENV/bin/python" -m build --sdist --no-isolation --verbose \
  --outdir "$ARTIFACTS/sdist" "$SNAPSHOT" \
  >"$LOGS/sdist-build.log" 2>&1

shopt -s nullglob
sdists=("$ARTIFACTS/sdist"/*.tar.gz)
[[ ${#sdists[@]} -eq 1 ]] || { echo "expected exactly one sdist" >&2; exit 1; }
SDIST_BUILD=${sdists[0]}

hash_file() {
  "$FRONTEND_VENV/bin/python" - "$1" <<'PY'
import hashlib
import sys

digest = hashlib.sha256()
with open(sys.argv[1], "rb") as stream:
    for block in iter(lambda: stream.read(1024 * 1024), b""):
        digest.update(block)
print(digest.hexdigest())
PY
}

SDIST_BUILD_HASH=$(hash_file "$SDIST_BUILD")
"$FRONTEND_VENV/bin/python" "$PACKAGING/artifact_inspector.py" sdist "$SDIST_BUILD" \
  --extract "$EXTRACT" --publish-directory "$VALIDATED_ARTIFACTS/sdist" \
  --manifest "$DIAGNOSTICS/sdist-manifest.json" \
  >"$LOGS/sdist-inspection.log" 2>&1

manifest_value() {
  "$FRONTEND_VENV/bin/python" -I - "$1" "$2" <<'PY'
import json
import sys
with open(sys.argv[1], encoding="utf-8") as stream:
    print(json.load(stream)[sys.argv[2]])
PY
}
SDIST=$(manifest_value "$DIAGNOSTICS/sdist-manifest.json" published_path)
SDIST_HASH=$(manifest_value "$DIAGNOSTICS/sdist-manifest.json" archive_sha256)
[[ "$SDIST_HASH" == "$SDIST_BUILD_HASH" ]] || { echo "sdist inspection digest mismatch" >&2; exit 1; }
printf '%s  %s\n' "$SDIST_HASH" "$SDIST" >"$DIAGNOSTICS/sdist.sha256"

extracted_roots=("$EXTRACT"/*)
[[ ${#extracted_roots[@]} -eq 1 && -d ${extracted_roots[0]} ]] || {
  echo "sdist extraction did not produce one root" >&2
  exit 1
}
EXTRACTED_ROOT=${extracted_roots[0]}
if find "$EXTRACTED_ROOT" -name .git -print -quit | grep -q .; then
  echo "extracted sdist contains Git metadata" >&2
  exit 1
fi
if [[ -e "$EXTRACTED_ROOT/extern/pybind11" ]] && \
   find "$EXTRACTED_ROOT/extern/pybind11" -mindepth 1 -print -quit | grep -q .; then
  echo "extracted sdist contains extern/pybind11 submodule contents" >&2
  exit 1
fi

WHEEL_SOURCE="$EXTRACTED_ROOT"
(
  cd "$WHEEL_SOURCE"
  TMPDIR="$DIAGNOSTICS/tmp/wheel-build" PIP_CACHE_DIR="$DIAGNOSTICS/cache/wheel-build" \
    "$FRONTEND_VENV/bin/python" -m build --wheel --verbose \
    --outdir "$ARTIFACTS/wheel"
) >"$LOGS/wheel-build.log" 2>&1

wheels=("$ARTIFACTS/wheel"/*.whl)
[[ ${#wheels[@]} -eq 1 ]] || { echo "expected exactly one wheel" >&2; exit 1; }
WHEEL_BUILD=${wheels[0]}
WHEEL_BUILD_HASH=$(hash_file "$WHEEL_BUILD")
"$FRONTEND_VENV/bin/python" "$PACKAGING/artifact_inspector.py" wheel "$WHEEL_BUILD" \
  --publish-directory "$VALIDATED_ARTIFACTS/wheel" \
  --manifest "$DIAGNOSTICS/wheel-manifest.json" \
  >"$LOGS/wheel-inspection.log" 2>&1
WHEEL=$(manifest_value "$DIAGNOSTICS/wheel-manifest.json" published_path)
WHEEL_HASH=$(manifest_value "$DIAGNOSTICS/wheel-manifest.json" archive_sha256)
[[ "$WHEEL_HASH" == "$WHEEL_BUILD_HASH" ]] || { echo "wheel inspection digest mismatch" >&2; exit 1; }
printf '%s  %s\n' "$WHEEL_HASH" "$WHEEL" >"$DIAGNOSTICS/wheel.sha256"

verify_hash() {
  local artifact=$1 expected=$2
  "$FRONTEND_VENV/bin/python" -I - "$artifact" "$expected" "$VALIDATED_ARTIFACTS" <<'PY'
import hashlib
import os
import stat
import sys

path = os.path.abspath(sys.argv[1])
expected = sys.argv[2]
trusted_root = os.path.realpath(sys.argv[3])
trusted_parent = os.path.realpath(os.path.dirname(path))
# pip accepts only a pathname, not an inherited descriptor.  The remaining
# boundary is therefore the verifier's private 0700 diagnostics tree.  Bind
# and hash the 0600 leaf immediately before handing that pathname to pip.
if os.path.commonpath((trusted_parent, trusted_root)) != trusted_root:
    raise SystemExit("published artifact escaped the trusted private directory")
descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
try:
    before = os.fstat(descriptor)
    if not stat.S_ISREG(before.st_mode) or stat.S_IMODE(before.st_mode) != 0o600:
        raise SystemExit("published artifact is not a private regular file")
    digest = hashlib.sha256()
    for block in iter(lambda: os.read(descriptor, 1024 * 1024), b""):
        digest.update(block)
    after = os.fstat(descriptor)
    current = os.stat(path, follow_symlinks=False)
    identity = lambda value: (
        value.st_dev, value.st_ino, stat.S_IFMT(value.st_mode), value.st_size,
        value.st_mtime_ns, value.st_ctime_ns,
    )
    if identity(before) != identity(after) or identity(after) != identity(current):
        raise SystemExit("published artifact identity changed before install")
    if digest.hexdigest() != expected:
        raise SystemExit("published artifact digest changed before install")
finally:
    os.close(descriptor)
PY
}

"$PYTHON" -m venv "$WHEEL_VENV"
verify_hash "$WHEEL" "$WHEEL_HASH"
TMPDIR="$DIAGNOSTICS/tmp/wheel" PIP_CACHE_DIR="$DIAGNOSTICS/cache/wheel" \
  "$WHEEL_VENV/bin/python" -m pip install --no-cache-dir "$WHEEL[mpi]" \
  >"$LOGS/wheel-install.log" 2>&1
"$WHEEL_VENV/bin/python" -m pip list --format=json \
  >"$DIAGNOSTICS/wheel-pip-list.json"
(
  cd "$EXTERNAL_CWD"
  TMPDIR="$DIAGNOSTICS/tmp/wheel" PIP_CACHE_DIR="$DIAGNOSTICS/cache/wheel" \
    PATH="$WHEEL_VENV/bin:$PATH" \
    "$WHEEL_VENV/bin/python" "$PACKAGING/runtime_smoke.py" --mode wheel \
    --expected-root "$WHEEL_VENV" --snapshot-root "$SNAPSHOT"
) >"$LOGS/wheel-runtime.log" 2>&1

"$PYTHON" -m venv "$SDIST_VENV"
verify_hash "$SDIST" "$SDIST_HASH"
TMPDIR="$DIAGNOSTICS/tmp/sdist" PIP_CACHE_DIR="$DIAGNOSTICS/cache/sdist" \
  "$SDIST_VENV/bin/python" -m pip install --no-cache-dir "$SDIST[mpi]" \
  >"$LOGS/sdist-install.log" 2>&1
"$SDIST_VENV/bin/python" -m pip list --format=json \
  >"$DIAGNOSTICS/sdist-pip-list.json"
(
  cd "$EXTERNAL_CWD"
  TMPDIR="$DIAGNOSTICS/tmp/sdist" PIP_CACHE_DIR="$DIAGNOSTICS/cache/sdist" \
    PATH="$SDIST_VENV/bin:$PATH" \
    "$SDIST_VENV/bin/python" "$PACKAGING/runtime_smoke.py" --mode sdist \
    --expected-root "$SDIST_VENV" --snapshot-root "$SNAPSHOT"
) >"$LOGS/sdist-runtime.log" 2>&1

"$PYTHON" -m venv "$EDITABLE_VENV"
TMPDIR="$DIAGNOSTICS/tmp/editable" PIP_CACHE_DIR="$DIAGNOSTICS/cache/editable" \
  "$EDITABLE_VENV/bin/python" -m pip install --no-cache-dir -e "$SNAPSHOT[mpi]" \
  >"$LOGS/editable-install.log" 2>&1
"$EDITABLE_VENV/bin/python" -m pip list --format=json \
  >"$DIAGNOSTICS/editable-pip-list.json"
(
  cd "$EXTERNAL_CWD"
  TMPDIR="$DIAGNOSTICS/tmp/editable" PIP_CACHE_DIR="$DIAGNOSTICS/cache/editable" \
    PATH="$EDITABLE_VENV/bin:$PATH" \
    "$EDITABLE_VENV/bin/python" "$PACKAGING/runtime_smoke.py" --mode editable \
    --expected-root "$EDITABLE_VENV" --snapshot-root "$SNAPSHOT"
) >"$LOGS/editable-runtime.log" 2>&1

echo "ALL PIP INSTALL PATHS VERIFIED"
