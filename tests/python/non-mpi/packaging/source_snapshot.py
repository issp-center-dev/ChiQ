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
_mkdir = os.mkdir
_read = os.read
_write = os.write
_capability_override = None


_GIT_EXEC_HELPER = """\
import os
import sys

repo_fd = int(sys.argv[1])
os.fchdir(repo_fd)
environment = {
    key: value
    for key, value in os.environ.items()
    if not key.upper().startswith("GIT_")
}
environment["GIT_CONFIG_GLOBAL"] = os.devnull
environment["GIT_CONFIG_NOSYSTEM"] = "1"
os.execvpe("git", ["git"] + sys.argv[2:], environment)
"""


def _sanitized_git_environment():
    """Remove every inherited Git override at the checkout trust boundary."""
    environment = {
        key: value
        for key, value in os.environ.items()
        if not key.upper().startswith("GIT_")
    }
    environment["GIT_CONFIG_GLOBAL"] = os.devnull
    environment["GIT_CONFIG_NOSYSTEM"] = "1"
    return environment


def _git(repo_fd, *args):
    command = [
        sys.executable,
        "-I",
        "-c",
        _GIT_EXEC_HELPER,
        str(repo_fd),
    ] + list(args)

    try:
        return subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            pass_fds=(repo_fd,),
            env=_sanitized_git_environment(),
        ).stdout
    except subprocess.CalledProcessError as error:
        detail = error.stderr.decode("utf-8", "replace").strip()
        raise SnapshotError("Git command failed: %s" % detail) from error
    except (OSError, subprocess.SubprocessError) as error:
        raise SnapshotError("cannot bind Git to the open checkout") from error


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
        components = _validate_relative_path(path)
        if stage != 0:
            raise SnapshotError("unmerged/non-stage-0 index entry for %s" % path)
        if mode == SYMLINK_MODE:
            raise SnapshotError("tracked symlink is unsupported: %s" % path)
        if mode not in REGULAR_MODES and mode != GITLINK_MODE:
            raise SnapshotError("unsupported index mode %s for %s" % (mode, path))
        normalized_prefix = []
        for position, component in enumerate(components):
            normalized_prefix.append(unicodedata.normalize("NFC", component).casefold())
            key = tuple(normalized_prefix)
            display = "/".join(components[: position + 1])
            previous = collision_keys.get(key)
            if previous is not None and previous != display:
                raise SnapshotError(
                    "normalized path collision: %s and %s" % (previous, display)
                )
            collision_keys[key] = display
        entries.append((path, mode))
    return entries


def _is_package_artifact(path):
    components = [
        unicodedata.normalize("NFC", component).casefold()
        for component in path.split("/")
    ]
    package_components = PACKAGE_PREFIX.rstrip("/").split("/")
    if components[: len(package_components)] != package_components:
        return False
    normalized_path = "/".join(components)
    return "__pycache__" in components or normalized_path.endswith(
        NATIVE_SUFFIXES + (".pyc",)
    )


def _nul_paths(raw):
    return [_decode_path(item) for item in raw.split(b"\0") if item]


def _untracked_paths(repo_fd):
    raw = _git(repo_fd, "status", "--porcelain=v1", "-z", "--untracked-files=all")
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


def _ignored_paths(repo_fd):
    return sorted(
        _nul_paths(
            _git(
                repo_fd,
                "ls-files",
                "-z",
                "--others",
                "--ignored",
                "--exclude-standard",
            )
        )
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
    if os.mkdir not in getattr(os, "supports_dir_fd", set()):
        raise SnapshotError("unsupported descriptor safety: mkdir(dir_fd=...)")
    if os.stat not in getattr(os, "supports_follow_symlinks", set()):
        raise SnapshotError("unsupported descriptor safety: stat(follow_symlinks=False)")
    if os.name != "posix" or not hasattr(os, "fchdir"):
        raise SnapshotError("unsupported descriptor safety: fchdir-bound Git")


def _same_identity(left, right):
    return (
        left.st_dev == right.st_dev
        and left.st_ino == right.st_ino
        and stat.S_IFMT(left.st_mode) == stat.S_IFMT(right.st_mode)
    )


def _source_state(value):
    """Return all source attributes whose change invalidates copied bytes."""
    return (
        value.st_dev,
        value.st_ino,
        stat.S_IFMT(value.st_mode),
        bool(value.st_mode & 0o111),
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
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


def _close_preserving_error(descriptor):
    try:
        os.close(descriptor)
    except OSError:
        pass


def _close_all(descriptors):
    primary_error = sys.exc_info()[1]
    close_error = None
    for descriptor in descriptors:
        if descriptor is None:
            continue
        try:
            os.close(descriptor)
        except BaseException as error:
            if close_error is None:
                close_error = error
    if primary_error is None and close_error is not None:
        raise close_error


def _open_validated_directory(path, description):
    try:
        before = _stat(os.fspath(path), follow_symlinks=False)
    except (OSError, TypeError, NotImplementedError) as error:
        raise SnapshotError("cannot no-follow stat %s" % description) from error
    if not stat.S_ISDIR(before.st_mode):
        raise SnapshotError("%s is not a directory" % description)
    try:
        descriptor = _open(
            os.fspath(path), os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, 0o777
        )
    except (OSError, TypeError, NotImplementedError) as error:
        raise SnapshotError("cannot safely open %s" % description) from error
    try:
        after = os.fstat(descriptor)
        if not stat.S_ISDIR(after.st_mode) or not _same_identity(before, after):
            raise SnapshotError("%s replacement race" % description)
    except BaseException:
        _close_preserving_error(descriptor)
        raise
    return descriptor


def _create_destination_root(destination):
    name = destination.name
    if name in ("", ".", "..") or "/" in name:
        raise SnapshotError("invalid destination path")
    parent_fd = _open_validated_directory(destination.parent, "destination parent")
    descriptor = None
    try:
        _mkdir(name, 0o700, dir_fd=parent_fd)
        before = _safe_stat(name, parent_fd, "destination root")
        if not stat.S_ISDIR(before.st_mode):
            raise SnapshotError("destination root is not a directory or is a symlink")
        descriptor = _safe_open(
            name,
            os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
            parent_fd,
            "destination root",
        )
        try:
            after = os.fstat(descriptor)
            if not stat.S_ISDIR(after.st_mode) or not _same_identity(before, after):
                raise SnapshotError("destination root replacement race")
            current = _safe_stat(name, parent_fd, "destination root")
            if not _same_identity(after, current):
                raise SnapshotError("destination root replacement race")
            os.fchmod(descriptor, 0o700)
        except BaseException:
            _close_preserving_error(descriptor)
            raise
        return descriptor
    finally:
        try:
            _close_all([parent_fd])
        except BaseException:
            _close_all([descriptor])
            raise


def _open_source_leaf(root_fd, path, index_mode):
    components = _validate_relative_path(path)
    parent_fd = root_fd
    owned = []
    leaf_fd = None
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
            try:
                after = os.fstat(child_fd)
                if not stat.S_ISDIR(after.st_mode) or not _same_identity(before, after):
                    raise SnapshotError("ancestor replacement race for %s" % path)
            except BaseException:
                _close_preserving_error(child_fd)
                raise
            owned.append(child_fd)
            parent_fd = child_fd

        leaf = components[-1]
        before = _safe_stat(leaf, parent_fd, "leaf %s" % path)
        if not stat.S_ISREG(before.st_mode):
            raise SnapshotError("tracked leaf is not a regular file: %s" % path)
        leaf_fd = _safe_open(leaf, os.O_RDONLY | os.O_NOFOLLOW, parent_fd, "leaf %s" % path)
        try:
            after = os.fstat(leaf_fd)
            if not stat.S_ISREG(after.st_mode) or not _same_identity(before, after):
                raise SnapshotError("leaf replacement race for %s" % path)
            before_executable = bool(before.st_mode & 0o111)
            executable = bool(after.st_mode & 0o111)
            if before_executable != executable:
                raise SnapshotError("executable mode race for %s" % path)
            expected_executable = index_mode == "100755"
            if executable != expected_executable:
                raise SnapshotError("executable mode differs from index for %s" % path)
        except BaseException:
            _close_preserving_error(leaf_fd)
            raise
        return leaf_fd
    finally:
        try:
            _close_all(reversed(owned))
        except BaseException:
            _close_all([leaf_fd])
            raise


def _open_destination_parent(root_fd, components):
    parent_fd = root_fd
    owned = []
    result = None
    try:
        for component in components:
            try:
                _mkdir(component, 0o755, dir_fd=parent_fd)
            except FileExistsError:
                pass
            before = _safe_stat(
                component, parent_fd, "destination directory %s" % component
            )
            if not stat.S_ISDIR(before.st_mode):
                raise SnapshotError("destination ancestor is not a directory")
            child_fd = _safe_open(
                component,
                os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                parent_fd,
                "destination directory %s" % component,
            )
            try:
                after = os.fstat(child_fd)
                if not stat.S_ISDIR(after.st_mode) or not _same_identity(before, after):
                    raise SnapshotError(
                        "destination directory %s replacement race" % component
                    )
                current = _safe_stat(
                    component, parent_fd, "destination directory %s" % component
                )
                if not _same_identity(after, current):
                    raise SnapshotError(
                        "destination directory %s replacement race" % component
                    )
                os.fchmod(child_fd, 0o755)
            except BaseException:
                _close_preserving_error(child_fd)
                raise
            owned.append(child_fd)
            parent_fd = child_fd
        if not owned:
            result = os.dup(root_fd)
            return result
        result = os.dup(parent_fd)
        return result
    finally:
        try:
            _close_all(reversed(owned))
        except BaseException:
            _close_all([result])
            raise


def _copy_open_file(source_fd, destination_root_fd, path, output_mode):
    components = _validate_relative_path(path)
    parent_fd = _open_destination_parent(destination_root_fd, components[:-1])
    output_fd = None
    digest = hashlib.sha256()
    try:
        before = os.fstat(source_fd)
        if not stat.S_ISREG(before.st_mode):
            raise SnapshotError("tracked source is not regular during copy: %s" % path)
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
            block = _read(source_fd, 1024 * 1024)
            if not block:
                break
            digest.update(block)
            offset = 0
            while offset < len(block):
                written = _write(output_fd, block[offset:])
                if written <= 0:
                    raise SnapshotError("destination write made no progress for %s" % path)
                offset += written
        after = os.fstat(source_fd)
        if _source_state(before) != _source_state(after):
            raise SnapshotError("tracked source changed during copy: %s" % path)
        return digest.hexdigest()
    finally:
        _close_all([output_fd, parent_fd])


def _warn(stderr, message):
    stderr.write("source snapshot: warning: %s\n" % message)


def create_snapshot(repo, destination, stderr):
    """Copy validated tracked worktree bytes and return the manifest dict."""
    _require_capabilities()
    repo = Path(repo)
    destination = Path(destination)
    source_root_fd = _open_validated_directory(repo, "checkout root")
    destination_root_fd = None
    rows = []
    try:
        initial_index = _git(source_root_fd, "ls-files", "--stage", "-z")
        entries = _parse_index(initial_index)
        initial_head = _git(source_root_fd, "rev-parse", "HEAD")
        head = initial_head.decode("ascii").strip()

        tracked_artifacts = sorted(
            path
            for path, mode in entries
            if mode != GITLINK_MODE and _is_package_artifact(path)
        )
        if tracked_artifacts:
            raise SnapshotError(
                "tracked package cache/native artifact: %s" % ", ".join(tracked_artifacts)
            )

        untracked = _untracked_paths(source_root_fd)
        ignored = _ignored_paths(source_root_fd)
        package_artifacts = sorted(
            set(path for path in untracked + ignored if _is_package_artifact(path))
        )
        gitlinks = sorted(path for path, mode in entries if mode == GITLINK_MODE)
        for path in untracked:
            _warn(stderr, "excluded untracked path %s" % path)
        for path in package_artifacts:
            _warn(stderr, "excluded package cache/native artifact %s" % path)
        for path in gitlinks:
            _warn(stderr, "excluded gitlink %s" % path)

        destination_root_fd = _create_destination_root(destination)

        for path, index_mode in sorted(entries):
            if index_mode == GITLINK_MODE:
                continue
            source_fd = _open_source_leaf(source_root_fd, path, index_mode)
            try:
                digest = _copy_open_file(
                    source_fd, destination_root_fd, path, REGULAR_MODES[index_mode]
                )
            finally:
                _close_all([source_fd])
            rows.append({"path": path, "index_mode": index_mode, "sha256": digest})
        final_index = _git(source_root_fd, "ls-files", "--stage", "-z")
        final_head = _git(source_root_fd, "rev-parse", "HEAD")
        if final_index != initial_index or final_head != initial_head:
            raise SnapshotError("Git index/HEAD generation changed during snapshot")
        diff_stat = _git(source_root_fd, "diff", "HEAD", "--stat").decode("utf-8")
    finally:
        _close_all([destination_root_fd, source_root_fd])

    return {
        "head": head,
        "diff_stat": diff_stat,
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
