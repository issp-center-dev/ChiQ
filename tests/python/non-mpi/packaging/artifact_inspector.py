#!/usr/bin/env python3
"""Validate packaging artifacts before using them in installation tests.

This module is deliberately independent of the ChiQ package.  It is security
infrastructure for hostile build artifacts, and therefore validates the full
archive graph before creating an extraction destination.
"""

import argparse
from email import policy
from email.parser import BytesParser
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import sys
import tarfile
import unicodedata

from packaging.utils import canonicalize_name, parse_sdist_filename
from packaging.version import InvalidVersion, Version
import toml


MAX_ARCHIVE = 512 * 1024 * 1024
MAX_MEMBERS = 50000
MAX_FILE = 256 * 1024 * 1024
MAX_TOTAL = 2 * 1024 * 1024 * 1024
MAX_RATIO = 200

_NATIVE_SUFFIXES = (".so", ".pyd", ".dylib")
_REGULAR_TYPES = (tarfile.REGTYPE, tarfile.AREGTYPE)
_VERSION_PROVIDER = "scikit_build_core.metadata.regex"
_VERSION_INPUT = "python/package/chiq/__init__.py"
_DRIVE = re.compile(r"^[A-Za-z]:")


class ArtifactError(RuntimeError):
    """An artifact violates a structural, resource, or project contract."""


# Small seams used by limit, platform-resolution, and race tests.
_physical_size = os.path.getsize
_resolve_path = lambda path: Path(path).resolve(strict=False)
_open = os.open
_mkdir = os.mkdir


def _identity_component(component):
    return unicodedata.normalize("NFC", component).casefold()


def canonical_member_path(path):
    """Return strict POSIX components for one canonical archive member path."""
    if not isinstance(path, str) or not path:
        raise ArtifactError("invalid member path: %r" % path)
    if "\x00" in path or "\\" in path:
        raise ArtifactError("invalid member path: %r" % path)
    if path.startswith("/") or path.startswith("//") or _DRIVE.match(path):
        raise ArtifactError("absolute/drive/UNC member path: %r" % path)
    components = path.split("/")
    if any(component in ("", ".", "..") for component in components):
        raise ArtifactError("non-canonical member path: %r" % path)
    canonical = PurePosixPath(*components).as_posix()
    if canonical != path:
        raise ArtifactError("non-canonical member path: %r" % path)
    return tuple(components)


def _require_contained(root, candidate):
    """Resolve both paths and reject a candidate outside root."""
    resolved_root = _resolve_path(root)
    resolved_candidate = _resolve_path(candidate)
    try:
        resolved_candidate.relative_to(resolved_root)
    except ValueError as error:
        raise ArtifactError("resolved extraction path escapes destination") from error
    return resolved_candidate


def _checked_archive_size(path):
    try:
        size = _physical_size(os.fspath(path))
    except OSError as error:
        raise ArtifactError("cannot stat archive") from error
    if size < 0 or size > MAX_ARCHIVE:
        raise ArtifactError("archive physical size exceeds 512 MiB limit")
    return size


def _check_member_kind(member):
    if member.type in _REGULAR_TYPES:
        if any(
            key.startswith("GNU.sparse") or key == "SCHILY.realsize"
            for key in member.pax_headers
        ):
            raise ArtifactError("sparse tar members are forbidden")
        return "file"
    if member.type == tarfile.DIRTYPE:
        return "directory"
    raise ArtifactError("tar members must be regular files or directories; special member found")


def _validated_graph(members, archive_size):
    if len(members) > MAX_MEMBERS:
        raise ArtifactError("archive member count exceeds limit")

    entries = set()
    nodes = {}
    spellings = {}
    roots = set()
    result = []
    total = 0
    for member in members:
        components = canonical_member_path(member.name)
        kind = _check_member_kind(member)
        path = "/".join(components)
        if path in entries:
            raise ArtifactError("duplicate archive entry: %s" % path)
        entries.add(path)
        roots.add(components[0])

        normalized = []
        for index, component in enumerate(components):
            normalized.append(_identity_component(component))
            key = tuple(normalized)
            display = "/".join(components[: index + 1])
            old_spelling = spellings.get(key)
            if old_spelling is not None and old_spelling != display:
                raise ArtifactError("normalized path collision/alias: %s and %s" % (old_spelling, display))
            spellings[key] = display

        for index in range(1, len(components)):
            ancestor = "/".join(components[:index])
            if nodes.get(ancestor) == "file":
                raise ArtifactError("ancestor-file conflict at %s" % ancestor)
            nodes.setdefault(ancestor, "directory")
        previous_kind = nodes.get(path)
        if previous_kind is not None and previous_kind != kind:
            raise ArtifactError("file/directory conflict at %s" % path)
        if kind == "file" and previous_kind == "directory":
            raise ArtifactError("file conflicts with directory graph at %s" % path)
        nodes[path] = kind

        if member.size < 0:
            raise ArtifactError("negative file size")
        if kind == "file":
            if member.size > MAX_FILE:
                raise ArtifactError("regular file exceeds size limit")
            total += member.size
            if total > MAX_TOTAL:
                raise ArtifactError("total uncompressed size exceeds limit")
        elif member.size != 0:
            raise ArtifactError("directory has nonzero file size")
        result.append(
            {
                "path": path,
                "type": kind,
                "size": member.size,
                "mode": member.mode,
            }
        )

    if len(roots) != 1:
        raise ArtifactError("archive must have one shared top-level root")
    if archive_size == 0:
        if total:
            raise ArtifactError("aggregate compression ratio exceeds limit")
    elif total > MAX_RATIO * archive_size:
        raise ArtifactError("aggregate compression ratio exceeds limit")
    return roots.pop(), sorted(result, key=lambda row: row["path"]), total


def _read_exact(archive, member):
    stream = archive.extractfile(member)
    if stream is None:
        raise ArtifactError("cannot read regular member: %s" % member.name)
    remaining = member.size
    chunks = []
    while remaining:
        chunk = stream.read(min(1024 * 1024, remaining))
        if not chunk:
            raise ArtifactError("archive member is shorter than declared")
        chunks.append(chunk)
        remaining -= len(chunk)
    if stream.read(1):
        raise ArtifactError("archive member is longer than declared")
    return b"".join(chunks)


def _member_map(members):
    return {member.name: member for member in members if member.type in _REGULAR_TYPES}


def _require_file(paths, relative):
    if relative not in paths:
        raise ArtifactError("missing required sdist content: %s" % relative)


def _validate_content(relative_files):
    paths = set(relative_files)
    for required in (
        "pyproject.toml",
        "PKG-INFO",
        "LICENSE",
        "README.md",
        "CMakeLists.txt",
        "python/CMakeLists.txt",
        "src/CMakeLists.txt",
        "python/package/chiq/__init__.py",
        "python/package/bse/__init__.py",
        "python/package/bse_solver/__init__.py",
    ):
        _require_file(paths, required)

    requirements = (
        (lambda path: path.startswith("src/") and path.endswith((".cpp", ".cc", ".cxx")), "C++ source"),
        (lambda path: path.startswith("src/") and path.endswith((".h", ".hpp", ".hxx")), "C++ header"),
        (lambda path: path.startswith("cmake/") and not path.endswith("/"), "CMake helper"),
        (lambda path: path.startswith("python/package/chiq/cli/") and path.endswith(".py") and not path.endswith("/__init__.py"), "chiq/cli module"),
        (lambda path: path.startswith("python/package/chiq/point_group_data/") and path.endswith(".py") and not path.endswith("/__init__.py"), "point-group module"),
        (lambda path: path.startswith("tests/") and path.casefold().endswith((".h5", ".hdf5")), "test HDF5 fixture"),
    )
    for predicate, description in requirements:
        if not any(predicate(path) for path in paths):
            raise ArtifactError("missing required %s content" % description)

    for path in paths:
        components = [_identity_component(item) for item in path.split("/")]
        basename = components[-1]
        if basename.endswith(_NATIVE_SUFFIXES):
            raise ArtifactError("forbidden native binary in sdist: %s" % path)
        if "__pycache__" in components or basename.endswith(".pyc"):
            raise ArtifactError("forbidden cache artifact in sdist: %s" % path)
        build_directories = components[:-1]
        if any(
            component == "build"
            or component.startswith("build-")
            or component.startswith("build_")
            or component in ("dist", "_skbuild")
            for component in build_directories
        ):
            raise ArtifactError("forbidden build tree in sdist: %s" % path)
        installed_layout = any(
            components[index:index + 2] == ["lib", "bse-python"]
            for index in range(len(components) - 1)
        )
        if installed_layout or basename in ("chiqvars.sh", "install_manifest.txt"):
            raise ArtifactError("forbidden installed layout in sdist: %s" % path)


def _one_header(message, name):
    values = message.get_all(name, [])
    if len(values) != 1 or not str(values[0]).strip():
        raise ArtifactError("malformed or duplicate PKG-INFO %s" % name)
    return str(values[0]).strip()


def _parse_pkg_info(data):
    try:
        message = BytesParser(policy=policy.default).parsebytes(data)
    except (TypeError, ValueError) as error:
        raise ArtifactError("malformed PKG-INFO metadata") from error
    if message.defects:
        raise ArtifactError("malformed PKG-INFO metadata")
    return _one_header(message, "Name"), _one_header(message, "Version")


def _validated_version(value, source):
    try:
        return Version(value)
    except InvalidVersion as error:
        raise ArtifactError("invalid version in %s" % source) from error


def _parse_project_metadata(pyproject_data, init_data):
    try:
        document = toml.loads(pyproject_data.decode("utf-8"))
        project = document["project"]
        provider = document["tool"]["scikit-build"]["metadata"]["version"]
    except (UnicodeDecodeError, KeyError, TypeError, toml.TomlDecodeError) as error:
        raise ArtifactError("malformed pyproject metadata") from error
    if canonicalize_name(project.get("name", "")) != "chiq":
        raise ArtifactError("pyproject project name must remain chiq")
    dynamic = project.get("dynamic")
    if not isinstance(dynamic, list) or dynamic.count("version") != 1:
        raise ArtifactError("pyproject must declare one dynamic version")
    if "version" in project:
        raise ArtifactError("pyproject version must be dynamic")
    if provider.get("provider") != _VERSION_PROVIDER or provider.get("input") != _VERSION_INPUT:
        raise ArtifactError("dynamic regex provider must read chiq.__version__")
    expression = provider.get("regex")
    if not isinstance(expression, str) or "?P<value>" not in expression or "__version__" not in expression:
        raise ArtifactError("malformed dynamic version regex provider")
    try:
        text = init_data.decode("utf-8")
        matches = list(re.finditer(expression, text))
    except (UnicodeDecodeError, re.error) as error:
        raise ArtifactError("cannot extract version from chiq.__version__") from error
    if len(matches) != 1 or matches[0].groupdict().get("value") is None:
        raise ArtifactError("dynamic version provider did not find exactly one version")
    return project["name"], matches[0].group("value")


def _sdist_filename(path):
    try:
        return parse_sdist_filename(Path(path).name)
    except Exception as error:
        raise ArtifactError("malformed sdist filename") from error


def _validate_metadata(path, root, archive, members):
    files = _member_map(members)
    prefix = root + "/"
    relative_files = {
        name[len(prefix):]: member
        for name, member in files.items()
        if name.startswith(prefix)
    }
    _validate_content(relative_files)
    pkg_name, pkg_version = _parse_pkg_info(_read_exact(archive, relative_files["PKG-INFO"]))
    project_name, source_version = _parse_project_metadata(
        _read_exact(archive, relative_files["pyproject.toml"]),
        _read_exact(archive, relative_files[_VERSION_INPUT]),
    )
    filename_name, filename_version = _sdist_filename(path)
    root_name, root_version = _sdist_filename(root + ".tar.gz")

    names = (filename_name, root_name, pkg_name, project_name)
    if any(canonicalize_name(name) != "chiq" for name in names):
        raise ArtifactError("sdist filename/root/metadata project name disagreement")
    versions = (
        filename_version,
        root_version,
        _validated_version(pkg_version, "PKG-INFO"),
        _validated_version(source_version, "chiq.__version__"),
    )
    if any(version != versions[0] for version in versions[1:]):
        raise ArtifactError("sdist filename/root/metadata version disagreement")
    return str(versions[0])


def _inspect_open(path):
    archive_size = _checked_archive_size(path)
    try:
        archive = tarfile.open(os.fspath(path), "r:*")
        members = archive.getmembers()
        root, manifest_members, total = _validated_graph(members, archive_size)
        version = _validate_metadata(path, root, archive, members)
    except ArtifactError:
        if "archive" in locals():
            archive.close()
        raise
    except (OSError, EOFError, tarfile.TarError, UnicodeError, ValueError) as error:
        if "archive" in locals():
            archive.close()
        raise ArtifactError("malformed tar archive") from error
    manifest = {
        "kind": "sdist",
        "name": "chiq",
        "version": version,
        "root": root,
        "archive_size": archive_size,
        "uncompressed_size": total,
        "members": manifest_members,
    }
    return archive, members, manifest


def validate_sdist(path):
    """Return a validated manifest without extracting any member."""
    archive, _members, manifest = _inspect_open(path)
    archive.close()
    return manifest


def _directory_flags():
    return os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW


def _same_identity(left, right):
    return (
        left.st_dev == right.st_dev
        and left.st_ino == right.st_ino
        and stat.S_IFMT(left.st_mode) == stat.S_IFMT(right.st_mode)
    )


def _require_root_identity(destination, root_fd):
    opened = os.fstat(root_fd)
    current = os.stat(os.fspath(destination), follow_symlinks=False)
    if not _same_identity(opened, current) or not stat.S_ISDIR(current.st_mode):
        raise ArtifactError("extraction destination changed during extraction")


def _open_directory_chain(root_fd, components):
    descriptor = os.dup(root_fd)
    try:
        for component in components:
            next_descriptor = _open(component, _directory_flags(), dir_fd=descriptor)
            os.close(descriptor)
            descriptor = next_descriptor
        return descriptor
    except (OSError, TypeError, NotImplementedError) as error:
        try:
            os.close(descriptor)
        except OSError:
            pass
        raise ArtifactError("unsafe directory race during extraction") from error


def _make_directory(root_fd, components):
    parent_fd = _open_directory_chain(root_fd, components[:-1])
    try:
        try:
            _mkdir(components[-1], 0o755, dir_fd=parent_fd)
        except FileExistsError:
            existing = os.stat(components[-1], dir_fd=parent_fd, follow_symlinks=False)
            if not stat.S_ISDIR(existing.st_mode):
                raise ArtifactError("file/directory collision during extraction")
        descriptor = _open(components[-1], _directory_flags(), dir_fd=parent_fd)
        try:
            os.fchmod(descriptor, 0o755)
        finally:
            os.close(descriptor)
    except ArtifactError:
        raise
    except (OSError, TypeError, NotImplementedError) as error:
        raise ArtifactError("unsafe directory creation during extraction") from error
    finally:
        os.close(parent_fd)


def _copy_member(archive, member, root_fd, components):
    parent_fd = _open_directory_chain(root_fd, components[:-1])
    descriptor = None
    try:
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW
        descriptor = _open(components[-1], flags, 0o600, dir_fd=parent_fd)
        source = archive.extractfile(member)
        if source is None:
            raise ArtifactError("cannot read regular member during extraction")
        remaining = member.size
        while remaining:
            chunk = source.read(min(1024 * 1024, remaining))
            if not chunk:
                raise ArtifactError("extraction byte-count mismatch")
            view = memoryview(chunk)
            while view:
                written = os.write(descriptor, view)
                if written <= 0:
                    raise ArtifactError("short write during extraction")
                view = view[written:]
            remaining -= len(chunk)
        if source.read(1):
            raise ArtifactError("extraction byte-count mismatch")
        os.fchmod(descriptor, 0o755 if member.mode & stat.S_IXUSR else 0o644)
    except ArtifactError:
        raise
    except (OSError, TypeError, NotImplementedError) as error:
        raise ArtifactError("unsafe exclusive file creation during extraction") from error
    finally:
        if descriptor is not None:
            os.close(descriptor)
        os.close(parent_fd)


def extract_sdist(path, destination):
    """Validate the complete graph, then manually extract regular files."""
    archive, members, manifest = _inspect_open(path)
    destination = Path(destination)
    try:
        _mkdir(os.fspath(destination), 0o700)
    except (FileExistsError, OSError, TypeError, NotImplementedError) as error:
        archive.close()
        raise ArtifactError("extraction destination must be new and non-symlink") from error
    root_fd = None
    try:
        root_fd = _open(os.fspath(destination), _directory_flags())
        os.fchmod(root_fd, 0o700)
        _require_root_identity(destination, root_fd)
        resolved_root = _resolve_path(destination)
        directories = set()
        for member in members:
            components = canonical_member_path(member.name)
            _require_contained(resolved_root, destination.joinpath(*components))
            if member.type == tarfile.DIRTYPE:
                directories.add(components)
            for index in range(1, len(components)):
                directories.add(components[:index])
        for components in sorted(directories, key=lambda item: (len(item), item)):
            _make_directory(root_fd, components)
        for member in members:
            if member.type in _REGULAR_TYPES:
                _copy_member(archive, member, root_fd, canonical_member_path(member.name))
        _require_root_identity(destination, root_fd)
    except ArtifactError:
        raise
    except (OSError, RuntimeError, TypeError, NotImplementedError) as error:
        raise ArtifactError("unsafe extraction destination or race") from error
    finally:
        if root_fd is not None:
            os.close(root_fd)
        archive.close()
    return manifest


def _write_manifest(path, manifest):
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("x", encoding="utf-8") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True)
        stream.write("\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="kind", required=True)
    sdist = subparsers.add_parser("sdist")
    sdist.add_argument("archive")
    sdist.add_argument("--extract")
    sdist.add_argument("--manifest")
    args = parser.parse_args(argv)

    if args.manifest and args.extract:
        manifest_path = _resolve_path(args.manifest)
        extraction_path = _resolve_path(args.extract)
        try:
            manifest_path.relative_to(extraction_path)
        except ValueError:
            pass
        else:
            parser.error("--manifest must be outside --extract")
    if args.extract:
        manifest = extract_sdist(args.archive, args.extract)
    else:
        manifest = validate_sdist(args.archive)
    if args.manifest:
        _write_manifest(args.manifest, manifest)
    else:
        json.dump(manifest, sys.stdout, indent=2, sort_keys=True)
        sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
