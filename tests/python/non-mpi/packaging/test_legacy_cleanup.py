import os
from pathlib import Path
import subprocess
import sys

import pytest


REPOSITORY = Path(__file__).resolve().parents[4]
TEMPLATE = REPOSITORY / "cmake" / "RemoveObsoleteLegacySolver.cmake.in"


def _write_project(source):
    template = TEMPLATE.as_posix()
    (source / "CMakeLists.txt").write_text(
        """cmake_minimum_required(VERSION 3.15)
project(cleanup_test NONE)
include(GNUInstallDirs)
configure_file(\"%s\" \"${CMAKE_CURRENT_BINARY_DIR}/cleanup.cmake\" @ONLY)
install(SCRIPT \"${CMAKE_CURRENT_BINARY_DIR}/cleanup.cmake\")
""" % template,
        encoding="utf-8",
    )


def _configure(tmp_path, libdir="lib", prefix=None):
    source = tmp_path / "source"
    build = tmp_path / "build"
    source.mkdir()
    _write_project(source)
    if prefix is None:
        prefix = tmp_path / "prefix"
    result = subprocess.run(
        [
            "cmake", "-S", os.fspath(source), "-B", os.fspath(build),
            "-DCMAKE_INSTALL_PREFIX=%s" % prefix,
            "-DCMAKE_INSTALL_LIBDIR=%s" % libdir,
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return build, Path(prefix)


def _install(build, destdir=None):
    environment = os.environ.copy()
    if destdir is not None:
        environment["DESTDIR"] = os.fspath(destdir)
    return subprocess.run(
        ["cmake", "--install", os.fspath(build)],
        env=environment,
        capture_output=True,
        text=True,
    )


def _legacy_root(prefix, libdir="lib", destdir=None):
    libdir_path = Path(libdir)
    configured = libdir_path if libdir_path.is_absolute() else prefix / libdir_path
    if destdir is not None:
        configured = Path(os.fspath(destdir) + os.fspath(configured))
    return configured / "bse-python"


def _populate_obsolete(root):
    (root / "chiq").mkdir(parents=True)
    obsolete = (
        root / "bse_solver.old.so",
        root / "bse_solver.pyd",
        root / "bse_solver-debug.dylib",
        root / "chiq" / "_bse_solver.old.so",
        root / "chiq" / "_bse_solver.pyd",
        root / "chiq" / "_bse_solver-debug.dylib",
    )
    for path in obsolete:
        path.write_bytes(b"obsolete")
    (root / "bse_solver-link.so").symlink_to("missing-target")
    (root / "chiq" / "_bse_solver-link.so").symlink_to("missing-target")
    return obsolete


@pytest.mark.parametrize("libdir", ["lib/custom", "lib64"])
def test_cleanup_uses_configure_time_prefix_for_relative_libdir(tmp_path, libdir):
    configured_prefix = tmp_path / "configured-prefix"
    build, prefix = _configure(tmp_path, libdir=libdir, prefix=configured_prefix)
    root = _legacy_root(prefix, libdir)
    obsolete = _populate_obsolete(root)
    sentinels = (
        root / "bse_solver.txt",
        root / "unrelated.so",
        root / "chiq" / "_bse_solver.py",
        root / "chiq" / "neighbor.dylib",
    )
    for path in sentinels:
        path.write_bytes(b"sentinel")

    result = _install(build)

    assert result.returncode == 0, result.stdout + result.stderr
    assert not any(path.exists() for path in obsolete)
    assert not (root / "bse_solver-link.so").is_symlink()
    assert not (root / "chiq" / "_bse_solver-link.so").is_symlink()
    assert all(path.read_bytes() == b"sentinel" for path in sentinels)


def test_cleanup_combines_destdir_with_absolute_libdir(tmp_path):
    absolute_libdir = tmp_path / "absolute-lib"
    destdir = tmp_path / "stage"
    build, prefix = _configure(
        tmp_path, libdir=os.fspath(absolute_libdir), prefix=tmp_path / "ignored-prefix"
    )
    root = _legacy_root(prefix, absolute_libdir, destdir)
    obsolete = _populate_obsolete(root)
    sentinel = root / "keep.so"
    sentinel.write_bytes(b"keep")

    result = _install(build, destdir=destdir)

    assert result.returncode == 0, result.stdout + result.stderr
    assert not any(path.exists() for path in obsolete)
    assert sentinel.read_bytes() == b"keep"


def test_matching_symlink_to_directory_is_removed_without_touching_target(tmp_path):
    build, prefix = _configure(tmp_path)
    root = _legacy_root(prefix)
    root.mkdir(parents=True)
    target = tmp_path / "directory-target"
    target.mkdir()
    sentinel = target / "sentinel.txt"
    sentinel.write_bytes(b"unchanged")
    candidate = root / "bse_solver-link.so"
    candidate.symlink_to(target, target_is_directory=True)

    result = _install(build)

    assert result.returncode == 0, result.stdout + result.stderr
    assert not candidate.is_symlink()
    assert sentinel.read_bytes() == b"unchanged"


def test_missing_legacy_directory_is_successful_noop(tmp_path):
    build, _prefix = _configure(tmp_path)

    result = _install(build)

    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("relative", ["bse_solver-bad.so", "chiq/_bse_solver-bad.pyd"])
def test_matching_directory_fails_without_touching_sentinel(tmp_path, relative):
    build, prefix = _configure(tmp_path)
    root = _legacy_root(prefix)
    candidate = root / relative
    candidate.mkdir(parents=True)
    sentinel = root / "sentinel.txt"
    sentinel.write_bytes(b"unchanged")

    result = _install(build)

    assert result.returncode != 0
    assert candidate.is_dir()
    assert sentinel.read_bytes() == b"unchanged"


@pytest.mark.skipif(not hasattr(os, "mkfifo"), reason="FIFO is not available")
def test_matching_non_regular_file_fails_without_touching_sentinel(tmp_path):
    build, prefix = _configure(tmp_path)
    root = _legacy_root(prefix)
    root.mkdir(parents=True)
    candidate = root / "bse_solver-fifo.so"
    os.mkfifo(candidate)
    sentinel = root / "sentinel.txt"
    sentinel.write_bytes(b"unchanged")

    result = _install(build)

    assert result.returncode != 0
    assert candidate.exists()
    assert sentinel.read_bytes() == b"unchanged"


def test_symlink_ancestor_fails_closed(tmp_path):
    build, prefix = _configure(tmp_path)
    outside = tmp_path / "outside"
    outside.mkdir()
    (outside / "bse-python" / "chiq").mkdir(parents=True)
    candidate = outside / "bse-python" / "bse_solver-old.so"
    candidate.write_bytes(b"obsolete")
    sentinel = outside / "bse-python" / "sentinel.txt"
    sentinel.write_bytes(b"unchanged")
    prefix.mkdir()
    (prefix / "lib").symlink_to(outside)

    result = _install(build)

    assert result.returncode != 0
    assert candidate.read_bytes() == b"obsolete"
    assert sentinel.read_bytes() == b"unchanged"


def test_project_registers_cleanup_immediately_before_legacy_target_install():
    text = (REPOSITORY / "src" / "CMakeLists.txt").read_text(encoding="utf-8")
    cleanup = text.index("install(SCRIPT ${CHIQ_LEGACY_CLEANUP_SCRIPT})")
    target = text.index("install(TARGETS _bse_solver", cleanup)
    between = text[cleanup:target]
    assert between.count("install(") == 1
    assert cleanup < target


def test_python_directory_installs_exclude_caches_bytecode_and_native_files():
    text = (REPOSITORY / "python" / "CMakeLists.txt").read_text(encoding="utf-8")
    exclusions_start = text.index("set(_chiq_install_excludes")
    exclusions_end = text.index(")", exclusions_start)
    exclusions = text[exclusions_start:exclusions_end]
    assert 'PATTERN "__pycache__" EXCLUDE' in exclusions
    for suffix in ("*.pyc", "*.so", "*.pyd", "*.dylib"):
        assert 'PATTERN "%s" EXCLUDE' % suffix in exclusions
    for package in ("chiq", "bse", "bse_solver"):
        start = text.index("install(DIRECTORY package/%s" % package)
        end = text.index(")", start)
        block = text[start:end]
        assert "${_chiq_install_excludes}" in block
    assert "if(IS_SYMLINK" in text and "\\${_obsolete_bse}" in text
    assert "file(REMOVE" in text


def test_legacy_dot_py_commands_are_generated_deprecated_wrappers():
    text = (REPOSITORY / "python" / "CMakeLists.txt").read_text(encoding="utf-8")
    assert "from chiq.cli._deprecated import ${command_name}_py" in text
    assert "${command_name}_py()" in text
    assert "install(FILES ${scripts}" not in text
    assert "install(PROGRAMS \"${deprecated_wrapper}\"" in text


def test_legacy_verifier_covers_isolated_build_runtime_and_full_suites():
    text = (
        REPOSITORY / "tests" / "python" / "non-mpi" / "packaging"
        / "verify_legacy_install.sh"
    ).read_text(encoding="utf-8")
    assert "build_legacy" not in text
    assert 'CHIQ_WHEEL_BUILD=OFF' in text
    assert '-DTesting=ON' in text
    assert 'v2.13.6' in text
    assert 'ctest --output-on-failure' in text
    assert 'runtime_smoke.py' in text
    assert '--mode wheel' in text
    assert 'cd "$EXTERNAL_CWD"' in text
    assert '--ignore="$BUILD/tests/python/non-mpi/packaging"' in text
    assert '"$REPO/tests/python/non-mpi/packaging"' in text
    assert text.count('-c "$EMPTY_PYTEST_CONFIG"') == 2
    assert "PYTEST_DISABLE_PLUGIN_AUTOLOAD=1" in text


def test_legacy_verifier_executes_no_commands_before_empty_environment_boundary():
    script = (
        REPOSITORY / "tests" / "python" / "non-mpi" / "packaging"
        / "verify_legacy_install.sh"
    ).read_text(encoding="utf-8")
    assert script.startswith("#!/bin/sh\n")
    wrapper = script.split("exec /usr/bin/env -i", 1)[0]
    assert " -c " not in wrapper
    assert "command -v" not in wrapper
    assert "dirname" not in wrapper
    assert "$(" not in wrapper
    assert "[[" not in wrapper


def test_legacy_verifier_scrubs_ambient_code_and_tool_configuration(tmp_path):
    script = (
        REPOSITORY / "tests" / "python" / "non-mpi" / "packaging"
        / "verify_legacy_install.sh"
    )
    poison = tmp_path / "poison"
    poison.mkdir()
    site_marker = tmp_path / "sitecustomize-ran"
    bash_marker = tmp_path / "bash-env-ran"
    bash_env = tmp_path / "bash-env.sh"
    bash_env.write_text("/usr/bin/touch %s\n" % bash_marker, encoding="utf-8")
    (poison / "sitecustomize.py").write_text(
        "from pathlib import Path\nPath(%r).write_text('ran')\n" % str(site_marker),
        encoding="utf-8",
    )
    environment = {
        "HOME": os.environ.get("HOME", str(tmp_path)),
        "PATH": "/usr/bin:/bin",
        "TMPDIR": str(tmp_path),
        "PYTHON": os.fspath(Path(sys.executable).resolve()),
        "PYTHONPATH": str(poison),
        "PYTHONOPTIMIZE": "1",
        "PYTEST_ADDOPTS": "--collect-only",
        "PYTEST_PLUGINS": "poison_plugin",
        "PIP_CONFIG_FILE": str(tmp_path / "pip.conf"),
        "GIT_DIR": str(tmp_path / "invalid.git"),
        "BASH_ENV": str(bash_env),
    }

    result = subprocess.run(
        [os.fspath(script), "--audit-environment-only"],
        cwd=os.fspath(REPOSITORY),
        env=environment,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert "Environment policy: env -i explicit allowlist" in result.stdout
    assert "Environment audit OK" in result.stdout
    assert "PYTEST_DISABLE_PLUGIN_AUTOLOAD=1" in result.stdout
    assert not site_marker.exists()
    assert not bash_marker.exists()


def test_legacy_verifier_runtime_path_starts_with_selected_python(tmp_path):
    script = (
        REPOSITORY / "tests" / "python" / "non-mpi" / "packaging"
        / "verify_legacy_install.sh"
    )
    selected = Path(sys.executable).resolve()
    environment = {
        "HOME": os.environ.get("HOME", str(tmp_path)),
        "PATH": "/usr/bin:/bin",
        "TMPDIR": str(tmp_path),
        "PYTHON": os.fspath(selected),
        "GIT_DIR": str(tmp_path / "invalid.git"),
    }

    result = subprocess.run(
        [os.fspath(script), "--audit-environment-only"],
        cwd=os.fspath(REPOSITORY),
        env=environment,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    path_line = next(
        line for line in result.stdout.splitlines() if line.startswith("PATH=")
    )
    assert path_line.split("=", 1)[1].split(":", 1)[0] == os.fspath(selected.parent)
    assert "PYTHON=%s\n" % selected in result.stdout


def test_legacy_verifier_rejects_forged_internal_environment(tmp_path):
    script = (
        REPOSITORY / "tests" / "python" / "non-mpi" / "packaging"
        / "verify_legacy_install.sh"
    )
    environment = {
        "CHIQ_LEGACY_VERIFIER_INTERNAL": "1",
        "PATH": "/usr/bin:/bin",
        "PYTHON_REQUEST": os.fspath(Path(sys.executable).resolve()),
        "PYTHONPATH": os.fspath(tmp_path),
        "GIT_DIR": os.fspath(tmp_path / "invalid.git"),
    }

    result = subprocess.run(
        ["/bin/bash", os.fspath(script), "--audit-environment-only"],
        cwd=os.fspath(REPOSITORY),
        env=environment,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "unsafe internal legacy verifier environment" in result.stderr
    assert "Environment policy" not in result.stdout
