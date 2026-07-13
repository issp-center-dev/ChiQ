import hashlib
import io
import json
import errno
import os
from pathlib import Path
import stat
import subprocess
import sys

import pytest

import source_snapshot


def git(repo, *args, input_data=None):
    return subprocess.run(
        ["git", "-C", str(repo)] + list(args),
        input=input_data,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    ).stdout


@pytest.fixture
def repo(tmp_path):
    root = tmp_path / "repo"
    root.mkdir()
    git(root, "init", "-q")
    git(root, "config", "user.email", "snapshot@example.invalid")
    git(root, "config", "user.name", "Snapshot Test")
    (root / "README").write_text("initial\n")
    git(root, "add", "README")
    git(root, "commit", "-qm", "initial")
    return root


def snapshot(repo, tmp_path, name="snapshot"):
    stderr = io.StringIO()
    destination = tmp_path / name
    manifest = source_snapshot.create_snapshot(repo, destination, stderr)
    return destination, manifest, stderr.getvalue()


def assert_fd_closed(descriptor):
    with pytest.raises(OSError) as caught:
        os.fstat(descriptor)
    assert caught.value.errno == errno.EBADF


def test_copies_current_tracked_bytes_modes_and_sorted_manifest(repo, tmp_path):
    script = repo / "bin" / "run"
    script.parent.mkdir()
    script.write_bytes(b"#!/bin/sh\necho staged\n")
    script.chmod(0o755)
    git(repo, "add", "bin/run")
    (repo / "README").write_bytes(b"unstaged bytes\n")

    destination, manifest, warnings = snapshot(repo, tmp_path)

    assert (destination / "README").read_bytes() == b"unstaged bytes\n"
    assert (destination / "bin/run").read_bytes() == b"#!/bin/sh\necho staged\n"
    assert stat.S_IMODE((destination / "README").stat().st_mode) == 0o644
    assert stat.S_IMODE((destination / "bin/run").stat().st_mode) == 0o755
    assert stat.S_IMODE(destination.stat().st_mode) == 0o700
    assert [row["path"] for row in manifest["files"]] == ["README", "bin/run"]
    assert manifest["files"] == [
        {
            "path": "README",
            "index_mode": "100644",
            "sha256": hashlib.sha256(b"unstaged bytes\n").hexdigest(),
        },
        {
            "path": "bin/run",
            "index_mode": "100755",
            "sha256": hashlib.sha256(b"#!/bin/sh\necho staged\n").hexdigest(),
        },
    ]
    assert manifest["head"] == git(repo, "rev-parse", "HEAD").decode().strip()
    assert manifest["diff_stat"] == git(repo, "diff", "HEAD", "--stat").decode()
    assert warnings == ""


def test_warns_and_excludes_nonignored_untracked_files(repo, tmp_path):
    (repo / "notes.txt").write_text("do not copy")

    destination, manifest, warnings = snapshot(repo, tmp_path)

    assert not (destination / "notes.txt").exists()
    assert manifest["diagnostics"]["excluded_untracked"] == ["notes.txt"]
    assert "notes.txt" in warnings


def test_diagnoses_ignored_and_untracked_package_artifacts(repo, tmp_path):
    (repo / ".gitignore").write_text("*.pyc\n*.so\n")
    git(repo, "add", ".gitignore")
    git(repo, "commit", "-qm", "ignore caches")
    package = repo / "python" / "package" / "chiq"
    (package / "__pycache__").mkdir(parents=True)
    (package / "__pycache__" / "module.pyc").write_bytes(b"cache")
    (package / "_bse_solver.so").write_bytes(b"native")
    (package / "local.pyd").write_bytes(b"native")

    destination, manifest, warnings = snapshot(repo, tmp_path)

    artifacts = manifest["diagnostics"]["excluded_package_artifacts"]
    assert artifacts == [
        "python/package/chiq/__pycache__/module.pyc",
        "python/package/chiq/_bse_solver.so",
        "python/package/chiq/local.pyd",
    ]
    assert all(not (destination / path).exists() for path in artifacts)
    assert "package cache/native artifact" in warnings


@pytest.mark.parametrize(
    "path",
    [
        "python/package/chiq/__pycache__/module.py",
        "python/package/chiq/module.pyc",
        "python/package/chiq/_bse_solver.so",
        "python/package/chiq/_bse_solver.pyd",
        "python/package/chiq/_bse_solver.dylib",
    ],
)
def test_rejects_tracked_package_cache_or_native_entries(repo, tmp_path, path):
    tracked = repo / path
    tracked.parent.mkdir(parents=True, exist_ok=True)
    tracked.write_bytes(b"tracked contamination")
    git(repo, "add", path)

    with pytest.raises(source_snapshot.SnapshotError, match="tracked package cache/native"):
        source_snapshot.create_snapshot(repo, tmp_path / "snapshot", io.StringIO())


def test_skips_gitlinks_without_copying_contents(repo, tmp_path):
    head = git(repo, "rev-parse", "HEAD").decode().strip()
    git(repo, "update-index", "--add", "--cacheinfo", "160000,%s,vendor/dependency" % head)
    contents = repo / "vendor" / "dependency"
    contents.mkdir(parents=True)
    (contents / "CMakeLists.txt").write_text("attacker bytes")

    destination, manifest, warnings = snapshot(repo, tmp_path)

    assert not (destination / "vendor").exists()
    assert [row["path"] for row in manifest["files"]] == ["README"]
    assert "vendor/dependency" in warnings


def test_rejects_unmerged_index_entries(repo, tmp_path):
    blob = git(repo, "hash-object", "README").decode().strip()
    git(repo, "rm", "--cached", "-q", "README")
    records = (
        "100644 %s 1\tREADME\n100644 %s 2\tREADME\n" % (blob, blob)
    ).encode()
    git(repo, "update-index", "--index-info", input_data=records)

    with pytest.raises(source_snapshot.SnapshotError, match="stage 0|unmerged"):
        source_snapshot.create_snapshot(repo, tmp_path / "snapshot", io.StringIO())


def test_rejects_tracked_symlink(repo, tmp_path):
    link = repo / "link"
    link.symlink_to("README")
    git(repo, "add", "link")

    with pytest.raises(source_snapshot.SnapshotError, match="symlink"):
        source_snapshot.create_snapshot(repo, tmp_path / "snapshot", io.StringIO())


def test_index_parser_rejects_unknown_modes():
    raw = b"100600 0000000000000000000000000000000000000000 0\todd\0"

    with pytest.raises(source_snapshot.SnapshotError, match="unsupported index mode"):
        source_snapshot._parse_index(raw)


def test_rejects_nonregular_leaf(repo, tmp_path):
    os.unlink(str(repo / "README"))
    os.mkfifo(str(repo / "README"))

    with pytest.raises(source_snapshot.SnapshotError, match="regular file"):
        source_snapshot.create_snapshot(repo, tmp_path / "snapshot", io.StringIO())


def test_rejects_executable_mode_mismatch(repo, tmp_path):
    (repo / "README").chmod(0o755)

    with pytest.raises(source_snapshot.SnapshotError, match="executable mode"):
        source_snapshot.create_snapshot(repo, tmp_path / "snapshot", io.StringIO())


def test_rejects_symlink_ancestor(repo, tmp_path):
    tracked = repo / "safe" / "file.txt"
    tracked.parent.mkdir()
    tracked.write_text("safe")
    git(repo, "add", "safe/file.txt")
    git(repo, "commit", "-qm", "nested")
    outside = tmp_path / "outside"
    outside.mkdir()
    (outside / "file.txt").write_text("attacker")
    tracked.unlink()
    tracked.parent.rmdir()
    (repo / "safe").symlink_to(outside, target_is_directory=True)

    with pytest.raises(source_snapshot.SnapshotError, match="ancestor|directory|symlink"):
        source_snapshot.create_snapshot(repo, tmp_path / "snapshot", io.StringIO())


def test_leaf_replacement_between_stat_and_open_fails_closed(repo, tmp_path, monkeypatch):
    real_open = source_snapshot._open
    replaced = {"done": False}

    def replacing_open(path, flags, mode=0o777, *, dir_fd=None):
        if path == "README" and flags & os.O_RDONLY == os.O_RDONLY and not replaced["done"]:
            replaced["done"] = True
            (repo / "README").unlink()
            (repo / "README").write_bytes(b"attacker replacement\n")
        return real_open(path, flags, mode, dir_fd=dir_fd)

    monkeypatch.setattr(source_snapshot, "_open", replacing_open)

    with pytest.raises(source_snapshot.SnapshotError, match="replaced|race"):
        source_snapshot.create_snapshot(repo, tmp_path / "snapshot", io.StringIO())

    copied = tmp_path / "snapshot" / "README"
    assert not copied.exists() or copied.read_bytes() != b"attacker replacement\n"


def test_leaf_mode_change_between_stat_and_open_fails_closed(repo, tmp_path, monkeypatch):
    script = repo / "script"
    script.write_bytes(b"#!/bin/sh\n")
    script.chmod(0o755)
    git(repo, "add", "script")
    script.chmod(0o644)
    real_open = source_snapshot._open
    changed = {"done": False}

    def changing_open(path, flags, mode=0o777, *, dir_fd=None):
        if path == "script" and not changed["done"]:
            changed["done"] = True
            script.chmod(0o755)
        return real_open(path, flags, mode, dir_fd=dir_fd)

    monkeypatch.setattr(source_snapshot, "_open", changing_open)

    with pytest.raises(source_snapshot.SnapshotError, match="mode.*race|changed.*mode"):
        source_snapshot.create_snapshot(repo, tmp_path / "snapshot", io.StringIO())


def test_ancestor_replacement_after_open_uses_safe_descriptor(repo, tmp_path, monkeypatch):
    tracked = repo / "safe" / "file.txt"
    tracked.parent.mkdir()
    tracked.write_bytes(b"safe bytes\n")
    git(repo, "add", "safe/file.txt")
    real_open = source_snapshot._open
    replaced = {"done": False}

    def replacing_open(path, flags, mode=0o777, *, dir_fd=None):
        fd = real_open(path, flags, mode, dir_fd=dir_fd)
        if path == "safe" and flags & os.O_DIRECTORY and not replaced["done"]:
            replaced["done"] = True
            (repo / "safe").rename(repo / "opened-safe")
            (repo / "safe").mkdir()
            (repo / "safe" / "file.txt").write_bytes(b"attacker bytes\n")
        return fd

    monkeypatch.setattr(source_snapshot, "_open", replacing_open)

    destination, manifest, warnings = snapshot(repo, tmp_path)

    assert (destination / "safe/file.txt").read_bytes() == b"safe bytes\n"


def test_checkout_root_replacement_between_stat_and_open_fails_closed(
    repo, tmp_path, monkeypatch
):
    original = tmp_path / "original-repo"
    real_open = source_snapshot._open
    replaced = {"done": False}

    def replacing_open(path, flags, mode=0o777, *, dir_fd=None):
        if (
            os.fspath(path) == os.fspath(repo)
            and flags & os.O_DIRECTORY
            and not replaced["done"]
        ):
            replaced["done"] = True
            repo.rename(original)
            repo.mkdir()
            git(repo, "init", "-q")
            git(repo, "config", "user.email", "attacker@example.invalid")
            git(repo, "config", "user.name", "Attacker")
            (repo / "README").write_bytes(b"attacker bytes\n")
            git(repo, "add", "README")
            git(repo, "commit", "-qm", "attacker")
        return real_open(path, flags, mode, dir_fd=dir_fd)

    monkeypatch.setattr(source_snapshot, "_open", replacing_open)

    with pytest.raises(source_snapshot.SnapshotError, match="checkout root.*race"):
        source_snapshot.create_snapshot(repo, tmp_path / "snapshot", io.StringIO())

    copied = tmp_path / "snapshot" / "README"
    assert not copied.exists() or copied.read_bytes() != b"attacker bytes\n"


def test_git_queries_and_copy_stay_bound_to_open_checkout(repo, tmp_path, monkeypatch):
    original_head = git(repo, "rev-parse", "HEAD").decode().strip()
    original_diff_stat = git(repo, "diff", "HEAD", "--stat").decode()
    original = tmp_path / "opened-repo"
    real_open = source_snapshot._open
    replaced = {"done": False}

    def replacing_open(path, flags, mode=0o777, *, dir_fd=None):
        fd = real_open(path, flags, mode, dir_fd=dir_fd)
        if (
            os.fspath(path) == os.fspath(repo)
            and flags & os.O_DIRECTORY
            and not replaced["done"]
        ):
            replaced["done"] = True
            repo.rename(original)
            repo.mkdir()
            git(repo, "init", "-q")
            git(repo, "config", "user.email", "attacker@example.invalid")
            git(repo, "config", "user.name", "Attacker")
            (repo / "README").write_bytes(b"attacker bytes\n")
            git(repo, "add", "README")
            git(repo, "commit", "-qm", "attacker")
            (repo / "README").write_bytes(b"attacker modified bytes\n")
        return fd

    monkeypatch.setattr(source_snapshot, "_open", replacing_open)

    destination, manifest, warnings = snapshot(repo, tmp_path)

    assert (destination / "README").read_bytes() == b"initial\n"
    assert manifest["head"] == original_head
    assert manifest["diff_stat"] == original_diff_stat


def test_destination_directory_replacement_before_open_fails_closed(
    repo, tmp_path, monkeypatch
):
    destination = tmp_path / "snapshot"
    detached = tmp_path / "detached-snapshot"
    attacker = tmp_path / "attacker-directory"
    attacker.mkdir()
    (attacker / "marker").write_bytes(b"unchanged\n")
    real_open = source_snapshot._open
    replaced = {"done": False}

    def replacing_open(path, flags, mode=0o777, *, dir_fd=None):
        is_destination = os.fspath(path) == os.fspath(destination) or (
            path == destination.name and dir_fd is not None
        )
        if is_destination and flags & os.O_DIRECTORY and not replaced["done"]:
            replaced["done"] = True
            destination.rename(detached)
            attacker.rename(destination)
        return real_open(path, flags, mode, dir_fd=dir_fd)

    monkeypatch.setattr(source_snapshot, "_open", replacing_open)

    with pytest.raises(source_snapshot.SnapshotError, match="destination root.*race"):
        source_snapshot.create_snapshot(repo, destination, io.StringIO())

    assert (destination / "marker").read_bytes() == b"unchanged\n"
    assert not (destination / "README").exists()


def test_destination_symlink_swap_after_creation_fails_without_touching_target(
    repo, tmp_path, monkeypatch
):
    destination = tmp_path / "snapshot"
    external = tmp_path / "external"
    external.mkdir(mode=0o711)
    original_mode = stat.S_IMODE(external.stat().st_mode)
    real_mkdir = os.mkdir

    def swapping_mkdir(path, mode=0o777, *, dir_fd=None):
        result = real_mkdir(path, mode, dir_fd=dir_fd)
        if path == destination.name and dir_fd is not None:
            destination.rmdir()
            destination.symlink_to(external, target_is_directory=True)
        return result

    monkeypatch.setattr(source_snapshot, "_mkdir", swapping_mkdir, raising=False)

    with pytest.raises(source_snapshot.SnapshotError, match="destination root|symlink"):
        source_snapshot.create_snapshot(repo, destination, io.StringIO())

    assert stat.S_IMODE(external.stat().st_mode) == original_mode
    assert list(external.iterdir()) == []


def test_nested_destination_replacement_before_open_fails_closed(
    repo, tmp_path, monkeypatch
):
    script = repo / "bin" / "run"
    script.parent.mkdir()
    script.write_bytes(b"#!/bin/sh\n")
    script.chmod(0o755)
    git(repo, "add", "bin/run")
    destination = tmp_path / "snapshot"
    detached = tmp_path / "detached-bin"
    attacker = tmp_path / "attacker-bin"
    attacker.mkdir(mode=0o711)
    (attacker / "marker").write_bytes(b"unchanged\n")
    original_mode = stat.S_IMODE(attacker.stat().st_mode)
    real_open = source_snapshot._open
    destination_root = {"fd": None}
    replaced = {"done": False}

    def replacing_open(path, flags, mode=0o777, *, dir_fd=None):
        if path == destination.name and flags & os.O_DIRECTORY and dir_fd is not None:
            fd = real_open(path, flags, mode, dir_fd=dir_fd)
            destination_root["fd"] = fd
            return fd
        if (
            path == "bin"
            and flags & os.O_DIRECTORY
            and dir_fd == destination_root["fd"]
            and not replaced["done"]
        ):
            replaced["done"] = True
            (destination / "bin").rename(detached)
            attacker.rename(destination / "bin")
        return real_open(path, flags, mode, dir_fd=dir_fd)

    monkeypatch.setattr(source_snapshot, "_open", replacing_open)

    with pytest.raises(source_snapshot.SnapshotError, match="destination directory.*race"):
        source_snapshot.create_snapshot(repo, destination, io.StringIO())

    assert (destination / "bin/marker").read_bytes() == b"unchanged\n"
    assert stat.S_IMODE((destination / "bin").stat().st_mode) == original_mode
    assert not (destination / "bin/run").exists()


def test_checkout_root_fd_closes_when_fstat_raises(repo, tmp_path, monkeypatch):
    real_open = source_snapshot._open
    real_fstat = os.fstat
    target = {"fd": None, "raised": False}

    def recording_open(path, flags, mode=0o777, *, dir_fd=None):
        fd = real_open(path, flags, mode, dir_fd=dir_fd)
        if os.fspath(path) == os.fspath(repo) and flags & os.O_DIRECTORY:
            target["fd"] = fd
        return fd

    def failing_fstat(fd):
        if fd == target["fd"] and not target["raised"]:
            target["raised"] = True
            raise OSError(errno.EIO, "injected checkout-root fstat failure")
        return real_fstat(fd)

    monkeypatch.setattr(source_snapshot, "_open", recording_open)
    monkeypatch.setattr(source_snapshot.os, "fstat", failing_fstat)

    with pytest.raises(OSError, match="injected checkout-root"):
        source_snapshot.create_snapshot(repo, tmp_path / "snapshot", io.StringIO())

    assert_fd_closed(target["fd"])


@pytest.mark.parametrize("failure", ["fstat", "fchmod"])
def test_destination_root_fd_closes_after_open_failure(
    repo, tmp_path, monkeypatch, failure
):
    destination = tmp_path / "snapshot"
    real_open = source_snapshot._open
    real_fstat = os.fstat
    real_fchmod = os.fchmod
    target = {"fd": None, "raised": False}

    def recording_open(path, flags, mode=0o777, *, dir_fd=None):
        fd = real_open(path, flags, mode, dir_fd=dir_fd)
        if path == destination.name and flags & os.O_DIRECTORY and dir_fd is not None:
            target["fd"] = fd
        return fd

    def failing_fstat(fd):
        if failure == "fstat" and fd == target["fd"] and not target["raised"]:
            target["raised"] = True
            raise OSError(errno.EIO, "injected destination-root fstat failure")
        return real_fstat(fd)

    def failing_fchmod(fd, mode):
        if failure == "fchmod" and fd == target["fd"] and not target["raised"]:
            target["raised"] = True
            raise OSError(errno.EIO, "injected destination-root fchmod failure")
        return real_fchmod(fd, mode)

    monkeypatch.setattr(source_snapshot, "_open", recording_open)
    monkeypatch.setattr(source_snapshot.os, "fstat", failing_fstat)
    monkeypatch.setattr(source_snapshot.os, "fchmod", failing_fchmod)

    with pytest.raises(OSError, match="injected destination-root"):
        source_snapshot.create_snapshot(repo, destination, io.StringIO())

    assert_fd_closed(target["fd"])


@pytest.mark.parametrize("kind", ["ancestor", "leaf"])
def test_source_fd_closes_when_fstat_raises(repo, tmp_path, monkeypatch, kind):
    nested = repo / "safe" / "file.txt"
    nested.parent.mkdir()
    nested.write_text("safe")
    git(repo, "add", "safe/file.txt")
    real_open = source_snapshot._open
    real_fstat = os.fstat
    target = {"fd": None, "raised": False}

    def recording_open(path, flags, mode=0o777, *, dir_fd=None):
        fd = real_open(path, flags, mode, dir_fd=dir_fd)
        is_source_ancestor = kind == "ancestor" and path == "safe" and flags & os.O_DIRECTORY
        is_source_leaf = (
            kind == "leaf"
            and path == "README"
            and not flags & os.O_DIRECTORY
            and not flags & (os.O_WRONLY | os.O_RDWR)
        )
        if is_source_ancestor or is_source_leaf:
            target["fd"] = fd
        return fd

    def failing_fstat(fd):
        if fd == target["fd"] and not target["raised"]:
            target["raised"] = True
            raise OSError(errno.EIO, "injected source-%s fstat failure" % kind)
        return real_fstat(fd)

    monkeypatch.setattr(source_snapshot, "_open", recording_open)
    monkeypatch.setattr(source_snapshot.os, "fstat", failing_fstat)

    with pytest.raises(OSError, match="injected source-%s" % kind):
        source_snapshot.create_snapshot(repo, tmp_path / "snapshot", io.StringIO())

    assert_fd_closed(target["fd"])


@pytest.mark.parametrize("failure", ["fstat", "fchmod"])
def test_nested_destination_fd_closes_after_open_failure(
    repo, tmp_path, monkeypatch, failure
):
    script = repo / "bin" / "run"
    script.parent.mkdir()
    script.write_bytes(b"#!/bin/sh\n")
    script.chmod(0o755)
    git(repo, "add", "bin/run")
    destination = tmp_path / "snapshot"
    real_open = source_snapshot._open
    real_fstat = os.fstat
    real_fchmod = os.fchmod
    destination_root = {"fd": None}
    target = {"fd": None, "raised": False}

    def recording_open(path, flags, mode=0o777, *, dir_fd=None):
        fd = real_open(path, flags, mode, dir_fd=dir_fd)
        if path == destination.name and flags & os.O_DIRECTORY and dir_fd is not None:
            destination_root["fd"] = fd
        elif path == "bin" and flags & os.O_DIRECTORY and dir_fd == destination_root["fd"]:
            target["fd"] = fd
        return fd

    def failing_fstat(fd):
        if failure == "fstat" and fd == target["fd"] and not target["raised"]:
            target["raised"] = True
            raise OSError(errno.EIO, "injected nested-destination fstat failure")
        return real_fstat(fd)

    def failing_fchmod(fd, mode):
        if failure == "fchmod" and fd == target["fd"] and not target["raised"]:
            target["raised"] = True
            raise OSError(errno.EIO, "injected nested-destination fchmod failure")
        return real_fchmod(fd, mode)

    monkeypatch.setattr(source_snapshot, "_open", recording_open)
    monkeypatch.setattr(source_snapshot.os, "fstat", failing_fstat)
    monkeypatch.setattr(source_snapshot.os, "fchmod", failing_fchmod)

    with pytest.raises(OSError, match="injected nested-destination"):
        source_snapshot.create_snapshot(repo, destination, io.StringIO())

    assert_fd_closed(target["fd"])


@pytest.mark.parametrize("missing", ["nofollow", "dir_fd", "stat_nofollow"])
def test_fails_closed_when_descriptor_safety_is_unsupported(repo, tmp_path, monkeypatch, missing):
    monkeypatch.setattr(source_snapshot, "_capability_override", missing)

    with pytest.raises(source_snapshot.SnapshotError, match="unsupported"):
        source_snapshot.create_snapshot(repo, tmp_path / "snapshot", io.StringIO())


def test_fails_closed_without_descriptor_relative_mkdir(repo, tmp_path, monkeypatch):
    supported = set(os.supports_dir_fd)
    supported.discard(os.mkdir)
    monkeypatch.setattr(source_snapshot.os, "supports_dir_fd", supported)

    with pytest.raises(source_snapshot.SnapshotError, match="unsupported.*mkdir"):
        source_snapshot.create_snapshot(repo, tmp_path / "snapshot", io.StringIO())


@pytest.mark.parametrize(
    "names",
    [
        ("caf\N{LATIN SMALL LETTER E WITH ACUTE}.txt", "cafe\N{COMBINING ACUTE ACCENT}.txt"),
        ("README.txt", "readme.TXT"),
    ],
)
def test_rejects_unicode_nfc_casefold_path_aliases(names):
    records = []
    for name in names:
        records.append(
            b"100644 0000000000000000000000000000000000000000 0\t"
            + name.encode("utf-8")
            + b"\0"
        )
    with pytest.raises(source_snapshot.SnapshotError, match="collision"):
        source_snapshot._parse_index(b"".join(records))


def test_destination_must_be_new(repo, tmp_path):
    destination = tmp_path / "snapshot"
    destination.mkdir()

    with pytest.raises((FileExistsError, source_snapshot.SnapshotError)):
        source_snapshot.create_snapshot(repo, destination, io.StringIO())


def test_cli_writes_manifest_outside_snapshot(repo, tmp_path):
    destination = tmp_path / "snapshot"
    manifest_path = tmp_path / "snapshot-manifest.json"
    script = Path(source_snapshot.__file__)

    subprocess.run(
        [sys.executable, str(script), str(repo), str(destination), "--manifest", str(manifest_path)],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    manifest = json.loads(manifest_path.read_text())
    assert manifest["files"][0]["path"] == "README"
    assert not (destination / ".git").exists()
    assert not (destination / "snapshot-manifest.json").exists()
