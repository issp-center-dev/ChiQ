#!/usr/bin/env python3
"""Create a descriptor-safe snapshot of the tracked Git worktree.

This is verification infrastructure, not part of ChiQ's build backend.  It
deliberately requires both Git and a checkout so that unstaged tracked bytes
can be verified without trusting path-based copies.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import unicodedata


REGULAR_MODES = {"100644": 0o644, "100755": 0o755}
GITLINK_MODE = "160000"
SYMLINK_MODE = "120000"
NATIVE_SUFFIXES = (".so", ".pyd", ".dylib")
PACKAGE_PREFIX = "python/package/"


class SnapshotError(RuntimeError):
    """The checkout cannot be copied without weakening snapshot safety."""


# Small seams used by race and platform-capability tests.
_open = os.open
_stat = os.stat
_capability_override = None


def _git(repo, *args):
    command = ["git", "-C", os.fspath(repo)] + list(args)
    try:
        return subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        ).stdout
    except subprocess.CalledProcessError as error:
        detail = error.stderr.decode("utf-8", "replace").strip()
        raise SnapshotError("Git command failed: %s" % detail) from error


def _decode_path(value):
    try:
        return value.decode("utf-8")
    except UnicodeDecodeError as error:
        raise SnapshotError("tracked paths must be valid UTF-8") from error


def _validate_relative_path(path):
    if not path or path.startswith("/"):
        raise SnapshotError("invalid tracked path: %r" % path)
    components = path.split("/")
    if any(component in ("", ".", "..") for component in components):
        raise SnapshotError("invalid tracked path: %r" % path)
    return components


def _parse_index(raw):
    """Parse the NUL-delimited output of ``git ls-files --stage``."""
    entries = []
    collision_keys = {}
    for record in raw.split(b"\0"):
        if not record:
            continue
        try:
            metadata, raw_path = record.split(b"\t", 1)
            raw_mode, _object_id, raw_stage = metadata.split(b" ", 2)
            mode = raw_mode.decode("ascii")
            stage = int(raw_stage)
        except (UnicodeDecodeError, ValueError) as error:
            raise SnapshotError("malformed git index record") from error
        path = _decode_path(raw_path)
        _validate_relative_path(path)
        if stage != 0:
            raise SnapshotError("unmerged/non-stage-0 index entry for %s" % path)
        if mode == SYMLINK_MODE:
            raise SnapshotError("tracked symlink is unsupported: %s" % path)
        if mode not in REGULAR_MODES and mode != GITLINK_MODE:
            raise SnapshotError("unsupported index mode %s for %s" % (mode, path))
        key = unicodedata.normalize("NFC", path).casefold()
        previous = collision_keys.get(key)
        if previous is not None and previous != path:
            raise SnapshotError("normalized path collision: %s and %s" % (previous, path))
        collision_keys[key] = path
        entries.append((path, mode))
    return entries


def _is_package_artifact(path):
    if not path.startswith(PACKAGE_PREFIX):
        return False
    components = path.split("/")
    return "__pycache__" in components or path.lower().endswith(NATIVE_SUFFIXES + (".pyc",))


def _nul_paths(raw):
    return [_decode_path(item) for item in raw.split(b"\0") if item]


def _untracked_paths(repo):
    raw = _git(repo, "status", "--porcelain=v1", "-z", "--untracked-files=all")
    paths = []
    records = raw.split(b"\0")
    index = 0
    while index < len(records):
        record = records[index]
        index += 1
        if not record:
            continue
        if len(record) < 4:
            raise SnapshotError("malformed git status record")
        status_code = record[:2]
        path = _decode_path(record[3:])
        if status_code == b"??":
            paths.append(path)
        # Porcelain -z rename/copy entries carry an additional source path.
        if b"R" in status_code or b"C" in status_code:
            index += 1
    return sorted(paths)


def _ignored_paths(repo):
    return sorted(
        _nul_paths(_git(repo, "ls-files", "-z", "--others", "--ignored", "--exclude-standard"))
    )


def _require_capabilities():
    if _capability_override is not None:
        raise SnapshotError("unsupported descriptor safety: %s" % _capability_override)
    if not hasattr(os, "O_NOFOLLOW") or not hasattr(os, "O_DIRECTORY"):
        raise SnapshotError("unsupported descriptor safety: O_NOFOLLOW/O_DIRECTORY")
    if os.open not in getattr(os, "supports_dir_fd", set()):
        raise SnapshotError("unsupported descriptor safety: open(dir_fd=...)")
    if os.stat not in getattr(os, "supports_dir_fd", set()):
        raise SnapshotError("unsupported descriptor safety: stat(dir_fd=...)")
    if os.stat not in getattr(os, "supports_follow_symlinks", set()):
        raise SnapshotError("unsupported descriptor safety: stat(follow_symlinks=False)")


def _same_identity(left, right):
    return (
        left.st_dev == right.st_dev
        and left.st_ino == right.st_ino
        and stat.S_IFMT(left.st_mode) == stat.S_IFMT(right.st_mode)
    )


def _safe_stat(component, parent_fd, description):
    try:
        return _stat(component, dir_fd=parent_fd, follow_symlinks=False)
    except (OSError, TypeError, NotImplementedError) as error:
        raise SnapshotError("cannot no-follow stat %s" % description) from error


def _safe_open(component, flags, parent_fd, description):
    try:
        return _open(component, flags, 0o777, dir_fd=parent_fd)
    except (OSError, TypeError, NotImplementedError) as error:
        raise SnapshotError("cannot safely open %s" % description) from error


def _open_source_leaf(root_fd, path, index_mode):
    components = _validate_relative_path(path)
    parent_fd = root_fd
    owned = []
    try:
        for position, component in enumerate(components[:-1]):
            description = "ancestor %s of %s" % ("/".join(components[: position + 1]), path)
            before = _safe_stat(component, parent_fd, description)
            if not stat.S_ISDIR(before.st_mode):
                raise SnapshotError("non-directory or symlink ancestor for %s" % path)
            child_fd = _safe_open(
                component,
                os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                parent_fd,
                description,
            )
            after = os.fstat(child_fd)
            if not stat.S_ISDIR(after.st_mode) or not _same_identity(before, after):
                os.close(child_fd)
                raise SnapshotError("ancestor replacement race for %s" % path)
            owned.append(child_fd)
            parent_fd = child_fd

        leaf = components[-1]
        before = _safe_stat(leaf, parent_fd, "leaf %s" % path)
        if not stat.S_ISREG(before.st_mode):
            raise SnapshotError("tracked leaf is not a regular file: %s" % path)
        leaf_fd = _safe_open(leaf, os.O_RDONLY | os.O_NOFOLLOW, parent_fd, "leaf %s" % path)
        after = os.fstat(leaf_fd)
        if not stat.S_ISREG(after.st_mode) or not _same_identity(before, after):
            os.close(leaf_fd)
            raise SnapshotError("leaf replacement race for %s" % path)
        executable = bool(after.st_mode & 0o111)
        expected_executable = index_mode == "100755"
        if executable != expected_executable:
            os.close(leaf_fd)
            raise SnapshotError("executable mode differs from index for %s" % path)
        return leaf_fd
    finally:
        for descriptor in reversed(owned):
            os.close(descriptor)


def _open_destination_parent(root_fd, components):
    parent_fd = root_fd
    owned = []
    try:
        for component in components:
            try:
                os.mkdir(component, 0o755, dir_fd=parent_fd)
            except FileExistsError:
                pass
            child_fd = _safe_open(
                component,
                os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                parent_fd,
                "destination directory %s" % component,
            )
            if not stat.S_ISDIR(os.fstat(child_fd).st_mode):
                os.close(child_fd)
                raise SnapshotError("destination ancestor is not a directory")
            os.fchmod(child_fd, 0o755)
            owned.append(child_fd)
            parent_fd = child_fd
        if not owned:
            return os.dup(root_fd)
        result = os.dup(parent_fd)
        return result
    finally:
        for descriptor in reversed(owned):
            os.close(descriptor)


def _copy_open_file(source_fd, destination_root_fd, path, output_mode):
    components = _validate_relative_path(path)
    parent_fd = _open_destination_parent(destination_root_fd, components[:-1])
    output_fd = None
    digest = hashlib.sha256()
    try:
        try:
            output_fd = _open(
                components[-1],
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                output_mode,
                dir_fd=parent_fd,
            )
        except (OSError, TypeError, NotImplementedError) as error:
            raise SnapshotError("cannot create destination file %s exclusively" % path) from error
        os.fchmod(output_fd, output_mode)
        while True:
            block = os.read(source_fd, 1024 * 1024)
            if not block:
                break
            digest.update(block)
            offset = 0
            while offset < len(block):
                offset += os.write(output_fd, block[offset:])
        return digest.hexdigest()
    finally:
        if output_fd is not None:
            os.close(output_fd)
        os.close(parent_fd)


def _warn(stderr, message):
    stderr.write("source snapshot: warning: %s\n" % message)


def create_snapshot(repo, destination, stderr):
    """Copy validated tracked worktree bytes and return the manifest dict."""
    _require_capabilities()
    repo = Path(repo)
    destination = Path(destination)
    entries = _parse_index(_git(repo, "ls-files", "--stage", "-z"))

    tracked_artifacts = sorted(path for path, mode in entries if mode != GITLINK_MODE and _is_package_artifact(path))
    if tracked_artifacts:
        raise SnapshotError(
            "tracked package cache/native artifact: %s" % ", ".join(tracked_artifacts)
        )

    untracked = _untracked_paths(repo)
    ignored = _ignored_paths(repo)
    package_artifacts = sorted(set(path for path in untracked + ignored if _is_package_artifact(path)))
    gitlinks = sorted(path for path, mode in entries if mode == GITLINK_MODE)
    for path in untracked:
        _warn(stderr, "excluded untracked path %s" % path)
    for path in package_artifacts:
        _warn(stderr, "excluded package cache/native artifact %s" % path)
    for path in gitlinks:
        _warn(stderr, "excluded gitlink %s" % path)

    try:
        os.mkdir(os.fspath(destination), 0o700)
    except OSError:
        raise
    os.chmod(os.fspath(destination), 0o700)

    source_root_fd = None
    destination_root_fd = None
    rows = []
    try:
        try:
            source_root_fd = _open(
                os.fspath(repo), os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, 0o777
            )
        except (OSError, TypeError, NotImplementedError) as error:
            raise SnapshotError("cannot safely open checkout root") from error
        try:
            destination_root_fd = _open(
                os.fspath(destination), os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, 0o777
            )
        except (OSError, TypeError, NotImplementedError) as error:
            raise SnapshotError("cannot safely open destination root") from error

        for path, index_mode in sorted(entries):
            if index_mode == GITLINK_MODE:
                continue
            source_fd = _open_source_leaf(source_root_fd, path, index_mode)
            try:
                digest = _copy_open_file(
                    source_fd, destination_root_fd, path, REGULAR_MODES[index_mode]
                )
            finally:
                os.close(source_fd)
            rows.append({"path": path, "index_mode": index_mode, "sha256": digest})
    finally:
        if destination_root_fd is not None:
            os.close(destination_root_fd)
        if source_root_fd is not None:
            os.close(source_root_fd)

    return {
        "head": _git(repo, "rev-parse", "HEAD").decode("ascii").strip(),
        "diff_stat": _git(repo, "diff", "HEAD", "--stat").decode("utf-8"),
        "files": rows,
        "diagnostics": {
            "excluded_untracked": untracked,
            "excluded_package_artifacts": package_artifacts,
            "excluded_gitlinks": gitlinks,
        },
    }


def _is_within(path, root):
    try:
        return os.path.commonpath([os.fspath(path), os.fspath(root)]) == os.fspath(root)
    except ValueError:
        return False


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repo")
    parser.add_argument("destination")
    parser.add_argument("--manifest", required=True)
    arguments = parser.parse_args(argv)
    destination = Path(arguments.destination).resolve()
    manifest_path = Path(arguments.manifest).resolve()
    if _is_within(manifest_path, destination):
        parser.error("--manifest must be outside DESTINATION")
    manifest = create_snapshot(arguments.repo, arguments.destination, sys.stderr)
    with open(os.fspath(manifest_path), "w", encoding="utf-8") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
