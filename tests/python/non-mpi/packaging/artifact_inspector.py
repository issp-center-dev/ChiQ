#!/usr/bin/env python3
"""Validate packaging artifacts before using them in installation tests.

This module is deliberately independent of the ChiQ package.  It is security
infrastructure for hostile build artifacts, and therefore validates the full
archive graph before creating an extraction destination.
"""

import argparse
import base64
import bz2
import configparser
import ctypes
import csv
from email import policy
from email.parser import BytesParser
import errno
import gzip
import hashlib
import io
import json
import lzma
import os
from pathlib import Path, PurePosixPath
import re
import secrets
import stat
import struct
import sys
import tarfile
import tempfile
import unicodedata
import warnings
import zipfile
import zlib

from packaging.requirements import InvalidRequirement, Requirement
from packaging.specifiers import InvalidSpecifier, SpecifierSet
from packaging.utils import canonicalize_name, parse_sdist_filename, parse_wheel_filename
from packaging.version import InvalidVersion, Version
import toml

from contracts import CANONICAL_COMMANDS, NATIVE_SUFFIXES, POINT_GROUP_MODULES, SHIM_MODULES


MAX_ARCHIVE = 512 * 1024 * 1024
MAX_MEMBERS = 50000
MAX_FILE = 256 * 1024 * 1024
MAX_TOTAL = 2 * 1024 * 1024 * 1024
MAX_RATIO = 200
MAX_CONTROL_FILE = 16 * 1024 * 1024

_REGULAR_TYPES = (tarfile.REGTYPE, tarfile.AREGTYPE)
_VERSION_PROVIDER = "scikit_build_core.metadata.regex"
_VERSION_INPUT = "python/package/chiq/__init__.py"
_DRIVE = re.compile(r"^[A-Za-z]:")
_ZERO_BLOCK = b"\0" * tarfile.BLOCKSIZE
_PAX_TYPES = (tarfile.XHDTYPE, tarfile.XGLTYPE)
_MAX_TAR_PADDING = tarfile.RECORDSIZE
_EOCD = struct.Struct("<4s4H2LH")
_CENTRAL = struct.Struct("<4s6H3L5H2L")
_LOCAL = struct.Struct("<4s5H3L2H")
_ZIP_METHODS = (zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED)
_ZIP_UTF8 = 0x800
_ZIP_DEFLATE_OPTIONS = 0x6
_ZIP_OUTPUT_CHUNK = 1024 * 1024
_NATIVE_MODULE_NAME = re.compile(
    r"^_bse_solver(?:"
    r"\.(?:so|pyd|dylib)"
    r"|\.abi3\.(?:so|pyd|dylib)"
    r"|\.cpython-[0-9]+[dt]?-(?:darwin|"
    r"(?:(?:x86_64|aarch64|powerpc64le|s390x|i386|riscv64|loongarch64)-"
    r"(?:linux-gnu|linux-musl)|arm-linux-gnueabihf))\.so"
    r"|\.cp[0-9]+[dt]?-(?:win32|win_amd64|win_arm64)\.pyd"
    r")$"
)
_CORE_REQUIREMENTS = ("numpy>=1.23", "scipy", "more-itertools", "h5py", "toml")
_EXTRA_REQUIREMENTS = {
    "plot": ("matplotlib",),
    "mpi": ("mpi4py",),
    "dcore": ("dcore==4.2.0", "mpi4py"),
    "gpu": ("cupy",),
    "test": ("pytest",),
}
_SUPPORTED_CORE_METADATA = frozenset(("2.1", "2.2", "2.3", "2.4", "2.5"))
_CONSOLE_SCRIPTS = {
    name: "chiq.cli.%s:main" % name for name in CANONICAL_COMMANDS
}
_CONSOLE_SCRIPTS.update({
    name + ".py": "chiq.cli._deprecated:%s_py" % name
    for name in tuple(_CONSOLE_SCRIPTS)
})


class ArtifactError(RuntimeError):
    """An artifact violates a structural, resource, or project contract."""


# Small seams used by limit, platform-resolution, and race tests.
_physical_size = os.path.getsize
_resolve_path = lambda path: Path(path).resolve(strict=False)
_open = os.open
_mkdir = os.mkdir
_stat = os.stat
_fstat = os.fstat
_fchmod = os.fchmod
_write = os.write
_close = os.close
_capability_override = None
_after_directory_mkdir = lambda parent_fd, name: None
_before_publish_identity = lambda parent_fd, staging, destination: None


def _zip_entry_ratio_exceeded(size, compressed_size):
    return bool(size and (not compressed_size or size > MAX_RATIO * compressed_size))


def _require_capabilities():
    if _capability_override is not None:
        raise ArtifactError("unsupported POSIX descriptor capability: %s" % _capability_override)
    if os.name != "posix" or not hasattr(os, "O_DIRECTORY") or not hasattr(os, "O_NOFOLLOW"):
        raise ArtifactError("unsupported POSIX descriptor capabilities")
    if os.open not in getattr(os, "supports_dir_fd", set()):
        raise ArtifactError("unsupported open(dir_fd) capability")
    if os.stat not in getattr(os, "supports_dir_fd", set()):
        raise ArtifactError("unsupported stat(dir_fd) capability")
    if os.mkdir not in getattr(os, "supports_dir_fd", set()):
        raise ArtifactError("unsupported mkdir(dir_fd) capability")
    if os.stat not in getattr(os, "supports_follow_symlinks", set()):
        raise ArtifactError("unsupported no-follow stat capability")


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


def _copy_private_artifact(path, expected_size):
    source = None
    private = None
    primary_error = None
    try:
        source = open(os.fspath(path), "rb")
        opened_size = _fstat(source.fileno()).st_size
        if opened_size != expected_size or opened_size > MAX_ARCHIVE:
            raise ArtifactError("archive changed before private copy")
        private = tempfile.TemporaryFile(mode="w+b")
        remaining = expected_size
        while remaining:
            chunk = source.read(min(1024 * 1024, remaining))
            if not chunk:
                raise ArtifactError("archive is shorter than physical size")
            written = private.write(chunk)
            if written != len(chunk):
                raise ArtifactError("short write to private artifact copy")
            remaining -= len(chunk)
        if source.read(1):
            raise ArtifactError("archive grew during private copy")
        private.flush()
        private.seek(0)
        return private
    except BaseException as error:
        primary_error = error
        if private is not None:
            try:
                private.close()
            except BaseException:
                pass
        raise
    finally:
        if source is not None:
            try:
                source.close()
            except BaseException as error:
                if primary_error is None:
                    if private is not None:
                        try:
                            private.close()
                        except BaseException:
                            pass
                    raise ArtifactError("cannot close source artifact") from error


def _read_stream_exact(stream, size, description):
    chunks = []
    remaining = size
    while remaining:
        chunk = stream.read(min(1024 * 1024, remaining))
        if not chunk:
            raise ArtifactError("truncated %s" % description)
        chunks.append(chunk)
        remaining -= len(chunk)
    return b"".join(chunks)


def _discard_exact(stream, size, description="tar payload"):
    remaining = size
    while remaining:
        chunk = stream.read(min(1024 * 1024, remaining))
        if not chunk:
            raise ArtifactError("truncated %s" % description)
        remaining -= len(chunk)


def _strict_raw_string(field, description):
    terminator = field.find(b"\0")
    value = field if terminator < 0 else field[:terminator]
    if terminator >= 0 and any(field[terminator + 1:]):
        raise ArtifactError("ambiguous raw %s has non-NUL padding" % description)
    try:
        return value.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ArtifactError("raw tar %s must be UTF-8" % description) from error


def _tar_number(field, description):
    try:
        value = tarfile.nti(field)
    except (ValueError, tarfile.InvalidHeaderError) as error:
        raise ArtifactError("malformed tar %s" % description) from error
    if value is None or value < 0:
        raise ArtifactError("malformed tar %s" % description)
    return value


def _validate_header_checksum(header):
    expected = _tar_number(header[148:156], "checksum")
    check_header = header[:148] + b" " * 8 + header[156:]
    unsigned = sum(check_header)
    signed = sum(value if value < 128 else value - 256 for value in check_header)
    if expected not in (unsigned, signed):
        raise ArtifactError("invalid tar header checksum")


def _raw_header(header):
    _validate_header_checksum(header)
    type_ = header[156:157]
    name = _strict_raw_string(header[0:100], "name")
    linkname = _strict_raw_string(header[157:257], "linkname")
    prefix = _strict_raw_string(header[345:500], "prefix")
    path = prefix + "/" + name if prefix else name
    if type_ == tarfile.DIRTYPE and path.endswith("/") and not path.endswith("//"):
        path = path[:-1]
    if type_ in _PAX_TYPES:
        if path not in ("PaxHeader", "././@PaxHeader", "pax_global_header"):
            canonical_member_path(path)
    else:
        canonical_member_path(path)
    if linkname:
        canonical_member_path(linkname)
    return type_, _tar_number(header[124:136], "size")


def _parse_pax(payload):
    records = {}
    position = 0
    while position < len(payload):
        separator = payload.find(b" ", position)
        if separator < 0:
            raise ArtifactError("malformed PAX record length")
        raw_length = payload[position:separator]
        if not raw_length.isdigit() or raw_length.startswith(b"0"):
            raise ArtifactError("malformed PAX record length")
        length = int(raw_length)
        end = position + length
        if end > len(payload) or end <= separator + 2:
            raise ArtifactError("malformed PAX record length")
        record = payload[separator + 1:end]
        if not record.endswith(b"\n"):
            raise ArtifactError("malformed PAX record framing")
        raw_key, marker, raw_value = record[:-1].partition(b"=")
        if not marker or not raw_key:
            raise ArtifactError("malformed PAX record")
        try:
            key = raw_key.decode("utf-8")
            value = raw_value.decode("utf-8")
        except UnicodeDecodeError as error:
            raise ArtifactError("PAX records must be UTF-8") from error
        if key in records:
            raise ArtifactError("duplicate PAX key: %s" % key)
        if any(character in key for character in "\x00\r\n="):
            raise ArtifactError("malformed PAX record key")
        records[key] = value
        position = end
    for key in ("path", "linkpath"):
        if key in records:
            canonical_member_path(records[key])
    if "size" in records:
        value = records["size"]
        if not value.isascii() or not value.isdecimal():
            raise ArtifactError("PAX size must be a nonnegative canonical decimal")
        if value != "0" and value.startswith("0"):
            raise ArtifactError("PAX size must be a nonnegative canonical decimal")
        parsed_size = int(value)
        if parsed_size > sys.maxsize:
            raise ArtifactError("PAX size exceeds representable bounds")
    return records


def _compression_stream(raw_file):
    raw_file.seek(0)
    magic = raw_file.read(6)
    raw_file.seek(0)
    if magic.startswith(b"\x1f\x8b"):
        return gzip.GzipFile(fileobj=raw_file, mode="rb")
    if magic.startswith(b"BZh"):
        return bz2.BZ2File(raw_file, mode="rb")
    if magic.startswith(b"\xfd7zXZ\x00"):
        return lzma.LZMAFile(raw_file, mode="rb")
    return raw_file


def _merge_pax(target, incoming, other, description):
    for key, value in incoming.items():
        if key in target:
            raise ArtifactError("duplicate PAX key across %s headers: %s" % (description, key))
        if key in other and other[key] != value:
            raise ArtifactError("conflicting global/per-file PAX metadata: %s" % key)
        target[key] = value


def _effective_size(raw_size, global_pax, local_pax):
    value = local_pax.get("size", global_pax.get("size"))
    return raw_size if value is None else int(value)


def _validate_trailing_padding(stream):
    total = 0
    while True:
        chunk = stream.read(min(4096, _MAX_TAR_PADDING - total + 1))
        if not chunk:
            return
        if any(chunk):
            raise ArtifactError("trailing nonzero bytes after tar end marker")
        total += len(chunk)
        if total > _MAX_TAR_PADDING:
            raise ArtifactError("excessive tar padding exceeds limit")


def _preflight_tar(raw_file, archive_size):
    """Bound and disambiguate raw tar headers before semantic parsing."""
    stream = _compression_stream(raw_file)
    global_pax = {}
    local_pax = {}
    count = 0
    total = 0
    zero_blocks = 0
    try:
        while zero_blocks < 2:
            header = _read_stream_exact(stream, tarfile.BLOCKSIZE, "tar header")
            if header == _ZERO_BLOCK:
                zero_blocks += 1
                continue
            if zero_blocks:
                raise ArtifactError("nonzero tar header after end marker")
            type_, raw_size = _raw_header(header)
            count += 1
            if count > MAX_MEMBERS:
                raise ArtifactError("archive member count exceeds limit")
            if type_ not in _REGULAR_TYPES + (tarfile.DIRTYPE,) + _PAX_TYPES:
                raise ArtifactError("tar members must be regular files or directories; special member found")
            if type_ in _PAX_TYPES:
                size = raw_size
                if size > MAX_FILE:
                    raise ArtifactError("PAX metadata size exceeds limit")
            else:
                size = _effective_size(raw_size, global_pax, local_pax)
                if type_ in _REGULAR_TYPES and size > MAX_FILE:
                    raise ArtifactError("regular file size exceeds limit")
                if type_ == tarfile.DIRTYPE and size:
                    raise ArtifactError("directory has nonzero file size")
            total += size
            if total > MAX_TOTAL:
                raise ArtifactError("total uncompressed size exceeds limit")
            if archive_size == 0 or total > MAX_RATIO * archive_size:
                raise ArtifactError("aggregate compression ratio exceeds limit")

            padded_size = (size + tarfile.BLOCKSIZE - 1) // tarfile.BLOCKSIZE * tarfile.BLOCKSIZE
            if type_ in _PAX_TYPES:
                payload = _read_stream_exact(stream, size, "PAX payload")
                records = _parse_pax(payload)
                _discard_exact(stream, padded_size - size, "PAX padding")
                if type_ == tarfile.XGLTYPE:
                    _merge_pax(global_pax, records, {}, "global")
                else:
                    _merge_pax(local_pax, records, global_pax, "per-file")
                continue

            if local_pax:
                for key, value in local_pax.items():
                    if key in global_pax and global_pax[key] != value:
                        raise ArtifactError("conflicting global/per-file PAX metadata: %s" % key)
                local_pax = {}
            _discard_exact(stream, padded_size)
        if local_pax:
            raise ArtifactError("orphan per-file PAX metadata")
        _validate_trailing_padding(stream)
    except (OSError, EOFError, gzip.BadGzipFile, lzma.LZMAError) as error:
        raise ArtifactError("malformed or truncated compressed tar archive") from error
    finally:
        if stream is not raw_file:
            stream.close()


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
    primary_error = None
    try:
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
    except BaseException as error:
        primary_error = error
        raise
    finally:
        try:
            stream.close()
        except BaseException as error:
            if primary_error is None:
                raise ArtifactError("cannot close archive member stream") from error


def _member_map(members):
    return {member.name: member for member in members if member.type in _REGULAR_TYPES}


def _require_file(paths, relative):
    if relative not in paths:
        raise ArtifactError("missing required sdist content: %s" % relative)


def _validate_content(relative_files, relative_entries):
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

    for path, kind in relative_entries.items():
        components = [_identity_component(item) for item in path.split("/")]
        basename = components[-1]
        if basename.endswith(NATIVE_SUFFIXES):
            raise ArtifactError("forbidden native binary in sdist: %s" % path)
        if "__pycache__" in components or basename.endswith(".pyc"):
            raise ArtifactError("forbidden cache artifact in sdist: %s" % path)
        build_directories = components if kind == "directory" else components[:-1]
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
    relative_entries = {
        member.name[len(prefix):]: (
            "directory" if member.type == tarfile.DIRTYPE else "file"
        )
        for member in members
        if member.name.startswith(prefix)
    }
    _validate_content(relative_files, relative_entries)
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
    raw_file = None
    archive = None
    try:
        raw_file = _copy_private_artifact(path, archive_size)
        _preflight_tar(raw_file, archive_size)
        raw_file.seek(0)
        archive = tarfile.open(fileobj=raw_file, mode="r:*")
        archive._artifact_raw_file = raw_file
        members = archive.getmembers()
        root, manifest_members, total = _validated_graph(members, archive_size)
        version = _validate_metadata(path, root, archive, members)
    except ArtifactError:
        if archive is not None:
            _close_archive(archive)
        elif raw_file is not None:
            raw_file.close()
        raise
    except (OSError, EOFError, tarfile.TarError, UnicodeError, ValueError) as error:
        if archive is not None:
            _close_archive(archive)
        elif raw_file is not None:
            raw_file.close()
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


def _close_archive(archive):
    raw_file = getattr(archive, "_artifact_raw_file", None)
    primary_error = sys.exc_info()[1]
    close_error = None
    try:
        archive.close()
    except BaseException as error:
        close_error = error
    if raw_file is not None and not raw_file.closed:
        try:
            raw_file.close()
        except BaseException as error:
            if close_error is None:
                close_error = error
    if primary_error is None and close_error is not None:
        raise ArtifactError("cannot close archive resources") from close_error


def validate_sdist(path):
    """Return a validated manifest without extracting any member."""
    archive, _members, manifest = _inspect_open(path)
    _close_archive(archive)
    return manifest


def _directory_flags():
    return os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW


def _same_identity(left, right):
    return (
        left.st_dev == right.st_dev
        and left.st_ino == right.st_ino
        and stat.S_IFMT(left.st_mode) == stat.S_IFMT(right.st_mode)
    )


def _close_descriptors(descriptors):
    primary_error = sys.exc_info()[1]
    close_error = None
    for descriptor in descriptors:
        if descriptor is None:
            continue
        try:
            _close(descriptor)
        except BaseException as error:
            if close_error is None:
                close_error = error
    if primary_error is None and close_error is not None:
        raise ArtifactError("cannot close extraction descriptor") from close_error


def _safe_stat(name, parent_fd, description):
    try:
        return _stat(name, dir_fd=parent_fd, follow_symlinks=False)
    except (OSError, TypeError, NotImplementedError) as error:
        raise ArtifactError("cannot no-follow stat %s" % description) from error


def _open_bound_directory(name, parent_fd, description):
    before = _safe_stat(name, parent_fd, description)
    if not stat.S_ISDIR(before.st_mode):
        raise ArtifactError("%s is not a directory" % description)
    descriptor = None
    try:
        descriptor = _open(name, _directory_flags(), dir_fd=parent_fd)
        opened = _fstat(descriptor)
        current = _safe_stat(name, parent_fd, description)
        if not _same_identity(before, opened) or not _same_identity(opened, current):
            raise ArtifactError("directory identity race for %s" % description)
        if not stat.S_ISDIR(opened.st_mode):
            raise ArtifactError("%s is not a directory" % description)
        return descriptor
    except ArtifactError:
        if descriptor is not None:
            _close_descriptors([descriptor])
        raise
    except (OSError, TypeError, NotImplementedError) as error:
        if descriptor is not None:
            _close_descriptors([descriptor])
        raise ArtifactError("cannot safely open %s" % description) from error


def _open_validated_parent(path):
    descriptor = None
    try:
        before = _stat(os.fspath(path), follow_symlinks=False)
        if not stat.S_ISDIR(before.st_mode):
            raise ArtifactError("destination parent is not a directory")
        descriptor = _open(os.fspath(path), _directory_flags())
        opened = _fstat(descriptor)
        current = _stat(os.fspath(path), follow_symlinks=False)
        if not _same_identity(before, opened) or not _same_identity(opened, current):
            raise ArtifactError("destination parent identity race")
        return descriptor
    except ArtifactError:
        if descriptor is not None:
            _close_descriptors([descriptor])
        raise
    except (OSError, TypeError, NotImplementedError) as error:
        if descriptor is not None:
            _close_descriptors([descriptor])
        raise ArtifactError("cannot safely open destination parent") from error


def _create_bound_directory(parent_fd, name, mode, description):
    descriptor = None
    try:
        _mkdir(name, mode, dir_fd=parent_fd)
        _after_directory_mkdir(parent_fd, name)
        descriptor = _open_bound_directory(name, parent_fd, description)
        _fchmod(descriptor, mode)
        current = _safe_stat(name, parent_fd, description)
        opened = _fstat(descriptor)
        if not _same_identity(opened, current):
            raise ArtifactError("directory identity race for %s" % description)
        return descriptor
    except ArtifactError:
        if descriptor is not None:
            _close_descriptors([descriptor])
        raise
    except (OSError, TypeError, NotImplementedError) as error:
        if descriptor is not None:
            _close_descriptors([descriptor])
        raise ArtifactError("cannot safely create %s" % description) from error


def _open_directory_chain(root_fd, components):
    descriptor = os.dup(root_fd)
    try:
        for component in components:
            next_descriptor = _open_bound_directory(component, descriptor, component)
            previous_descriptor = descriptor
            descriptor = None
            try:
                _close_descriptors([previous_descriptor])
            except BaseException:
                _close_descriptors([next_descriptor])
                raise
            descriptor = next_descriptor
        return descriptor
    except BaseException:
        _close_descriptors([descriptor])
        raise


def _make_directory(root_fd, components):
    parent_fd = _open_directory_chain(root_fd, components[:-1])
    descriptor = None
    try:
        try:
            descriptor = _create_bound_directory(
                parent_fd,
                components[-1],
                0o755,
                "/".join(components),
            )
        except FileExistsError:
            descriptor = _open_bound_directory(
                components[-1], parent_fd, "/".join(components)
            )
            _fchmod(descriptor, 0o755)
    except ArtifactError:
        raise
    except (OSError, TypeError, NotImplementedError) as error:
        raise ArtifactError("unsafe directory creation during extraction") from error
    finally:
        _close_descriptors([descriptor, parent_fd])


def _copy_member(archive, member, root_fd, components):
    parent_fd = _open_directory_chain(root_fd, components[:-1])
    descriptor = None
    source = None
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
                written = _write(descriptor, view)
                if written <= 0:
                    raise ArtifactError("short write during extraction")
                view = view[written:]
            remaining -= len(chunk)
        if source.read(1):
            raise ArtifactError("extraction byte-count mismatch")
        _fchmod(descriptor, 0o755 if member.mode & stat.S_IXUSR else 0o644)
    except ArtifactError:
        raise
    except (OSError, TypeError, NotImplementedError) as error:
        raise ArtifactError("unsafe exclusive file creation during extraction") from error
    finally:
        primary_error = sys.exc_info()[1]
        close_error = None
        if source is not None:
            try:
                source.close()
            except BaseException as error:
                close_error = error
        try:
            _close_descriptors([descriptor, parent_fd])
        except BaseException as error:
            if close_error is None:
                close_error = error
        if primary_error is None and close_error is not None:
            raise ArtifactError("cannot close extracted member") from close_error


def _require_entry_identity(parent_fd, name, descriptor, description):
    try:
        current = _safe_stat(name, parent_fd, description)
    except ArtifactError as error:
        raise ArtifactError(
            "%s identity is unavailable; manual cleanup may be required" % description
        ) from error
    try:
        opened = _fstat(descriptor)
    except (OSError, TypeError, NotImplementedError) as error:
        raise ArtifactError("cannot fstat %s" % description) from error
    if not stat.S_ISDIR(opened.st_mode) or not _same_identity(opened, current):
        raise ArtifactError(
            "%s identity mismatch; manual cleanup may be required" % description
        )
    return opened


def _cleanup_staging(parent_fd, name, staging_fd, populated):
    if staging_fd is None:
        raise ArtifactError(
            "staging ownership is unavailable; manual cleanup may be required"
        )
    _require_entry_identity(parent_fd, name, staging_fd, "staging cleanup entry")
    if populated:
        raise ArtifactError(
            "populated staging is retained fail-closed; manual cleanup is required"
        )
    raise ArtifactError(
        "empty staging is retained fail-closed; manual cleanup is required"
    )


def _diagnose_cleanup(error):
    message = str(error)
    try:
        warnings.warn(message, RuntimeWarning)
        return
    except BaseException:
        pass
    try:
        sys.stderr.write("artifact inspector cleanup diagnostic: %s\n" % message)
    except BaseException:
        pass


def _entry_exists(parent_fd, name):
    try:
        _stat(name, dir_fd=parent_fd, follow_symlinks=False)
        return True
    except FileNotFoundError:
        return False
    except (OSError, TypeError, NotImplementedError) as error:
        raise ArtifactError("cannot check extraction destination") from error


def _publish_noreplace(parent_fd, source, destination):
    libc = ctypes.CDLL(None, use_errno=True)
    if sys.platform == "darwin":
        function = getattr(libc, "renameatx_np", None)
        flags = 0x00000004
    elif sys.platform.startswith("linux"):
        function = getattr(libc, "renameat2", None)
        flags = 1
    else:
        function = None
        flags = 0
    if function is None:
        raise ArtifactError("atomic no-replace directory publication is unsupported")
    function.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint]
    function.restype = ctypes.c_int
    result = function(
        parent_fd,
        os.fsencode(source),
        parent_fd,
        os.fsencode(destination),
        flags,
    )
    if result != 0:
        number = ctypes.get_errno()
        if number in (errno.EEXIST, errno.ENOTEMPTY):
            raise ArtifactError("extraction destination already exists")
        raise ArtifactError("cannot atomically publish extraction destination") from OSError(number, os.strerror(number))


def extract_sdist(path, destination):
    """Validate the complete graph, then manually extract regular files."""
    _require_capabilities()
    archive, members, manifest = _inspect_open(path)
    destination = Path(destination)
    parent_fd = None
    staging_fd = None
    renamed = False
    committed = False
    staging_name = None
    staging_populated = False
    try:
        parent = _resolve_path(destination.parent)
        if not destination.name or destination.name in (".", ".."):
            raise ArtifactError("invalid extraction destination")
        _require_contained(parent, parent / destination.name)
        parent_fd = _open_validated_parent(parent)
        if _entry_exists(parent_fd, destination.name):
            raise ArtifactError("extraction destination already exists")
        staging_name = ".%s.tmp-%s" % (
            destination.name,
            secrets.token_hex(8),
        )
        staging_fd = _create_bound_directory(
            parent_fd, staging_name, 0o700, "private extraction staging root"
        )
        directories = set()
        for member in members:
            components = canonical_member_path(member.name)
            if member.type == tarfile.DIRTYPE:
                directories.add(components)
            for index in range(1, len(components)):
                directories.add(components[:index])
        for components in sorted(directories, key=lambda item: (len(item), item)):
            staging_populated = True
            _make_directory(staging_fd, components)
        for member in members:
            if member.type in _REGULAR_TYPES:
                staging_populated = True
                _copy_member(archive, member, staging_fd, canonical_member_path(member.name))
        # Every archive/member/private-copy resource is closed before the
        # directory publication commit point.
        _close_archive(archive)
        archive = None
        if _entry_exists(parent_fd, destination.name):
            raise ArtifactError("extraction destination already exists")
        _before_publish_identity(parent_fd, staging_name, destination.name)
        _require_entry_identity(
            parent_fd, staging_name, staging_fd, "staging publication source"
        )
        # renameatx_np/renameat2 prevent replacement of the destination, but
        # pathname rename APIs cannot prevent a same-UID attacker from swapping
        # the source name.  The pre/post descriptor identity checks ensure such
        # a swap can never be reported as successful publication.
        _publish_noreplace(parent_fd, staging_name, destination.name)
        renamed = True
        _require_entry_identity(
            parent_fd, destination.name, staging_fd, "published destination"
        )
        committed = True
    except ArtifactError:
        raise
    except (OSError, RuntimeError, TypeError, NotImplementedError) as error:
        raise ArtifactError("unsafe extraction destination or race") from error
    finally:
        primary_error = sys.exc_info()[1]
        cleanup_error = None
        if not renamed and parent_fd is not None and staging_name is not None:
            try:
                _cleanup_staging(
                    parent_fd, staging_name, staging_fd, staging_populated
                )
            except BaseException as error:
                cleanup_error = error
                if primary_error is not None:
                    _diagnose_cleanup(error)
        try:
            _close_descriptors([staging_fd, parent_fd])
        except BaseException as error:
            if cleanup_error is None:
                cleanup_error = error
        if archive is not None:
            try:
                _close_archive(archive)
            except BaseException as error:
                if cleanup_error is None:
                    cleanup_error = error
        # Once destination identity is verified, publication is committed.
        # Cleanup-only close failures cannot truthfully turn that into failure.
        if not committed and primary_error is None and cleanup_error is not None:
            raise cleanup_error
    return manifest


def _zip_read_at(stream, offset, size, description):
    if offset < 0 or size < 0:
        raise ArtifactError("negative ZIP range for %s" % description)
    try:
        stream.seek(offset)
        data = stream.read(size)
    except OSError as error:
        raise ArtifactError("cannot read ZIP %s" % description) from error
    if len(data) != size:
        raise ArtifactError("out-of-range ZIP %s" % description)
    return data


def _zip_extra_fields(data):
    offset = 0
    while offset < len(data):
        if len(data) - offset < 4:
            raise ArtifactError("malformed ZIP extra field")
        identifier, size = struct.unpack_from("<HH", data, offset)
        offset += 4
        if size > len(data) - offset:
            raise ArtifactError("malformed ZIP extra field length")
        if identifier == 1:
            raise ArtifactError("ZIP64 is unsupported")
        offset += size


def _zip_name(raw_name, flags):
    if not raw_name or b"\0" in raw_name:
        raise ArtifactError("invalid ZIP member name")
    if not flags & _ZIP_UTF8 and any(byte >= 128 for byte in raw_name):
        raise ArtifactError("non-ASCII ZIP name lacks UTF-8 flag")
    try:
        return raw_name.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ArtifactError("ZIP member name is not valid UTF-8") from error


def _zip_components(name, is_directory):
    if is_directory:
        if not name.endswith("/") or name.endswith("//"):
            raise ArtifactError("non-canonical ZIP directory path")
        name = name[:-1]
    elif name.endswith("/"):
        raise ArtifactError("ambiguous ZIP file path")
    return canonical_member_path(name)


def _zip_graph(entries):
    exact = set()
    nodes = {}
    spellings = {}
    result = []
    for entry in entries:
        components = entry["components"]
        path = "/".join(components)
        kind = "directory" if entry["directory"] else "file"
        if path in exact:
            raise ArtifactError("duplicate ZIP archive entry: %s" % path)
        exact.add(path)
        normalized = []
        for index, component in enumerate(components):
            normalized.append(_identity_component(component))
            key = tuple(normalized)
            spelling = "/".join(components[:index + 1])
            previous = spellings.get(key)
            if previous is not None and previous != spelling:
                raise ArtifactError("normalized ZIP path collision/alias")
            spellings[key] = spelling
        for index in range(1, len(components)):
            ancestor = "/".join(components[:index])
            if nodes.get(ancestor) == "file":
                raise ArtifactError("ZIP ancestor-file conflict")
            nodes.setdefault(ancestor, "directory")
        previous_kind = nodes.get(path)
        if previous_kind is not None and previous_kind != kind:
            raise ArtifactError("ZIP file/directory conflict")
        if kind == "file" and previous_kind == "directory":
            raise ArtifactError("ZIP file conflicts with directory graph")
        nodes[path] = kind
        entry["path"] = path
        result.append(entry)
    return result


def _preflight_zip(stream, archive_size):
    tail_size = min(archive_size, 65557)
    tail_offset = archive_size - tail_size
    tail = _zip_read_at(stream, tail_offset, tail_size, "EOCD search window")
    relative = tail.rfind(b"PK\x05\x06")
    if relative < 0:
        raise ArtifactError("missing ZIP EOCD")
    eocd_offset = tail_offset + relative
    if archive_size - eocd_offset < _EOCD.size:
        raise ArtifactError("truncated ZIP EOCD")
    values = _EOCD.unpack(_zip_read_at(stream, eocd_offset, _EOCD.size, "EOCD"))
    _signature, disk, central_disk, disk_count, count, central_size, central_offset, comment_size = values
    if any(value == 0xFFFF for value in (disk, central_disk, disk_count, count)) or any(
        value == 0xFFFFFFFF for value in (central_size, central_offset)
    ):
        raise ArtifactError("ZIP64 is unsupported")
    if disk or central_disk:
        raise ArtifactError("multi-disk ZIP archives are forbidden")
    if disk_count != count:
        raise ArtifactError("inconsistent ZIP EOCD member count")
    if count > MAX_MEMBERS:
        raise ArtifactError("ZIP member count exceeds limit")
    if eocd_offset + _EOCD.size + comment_size != archive_size:
        raise ArtifactError("inconsistent ZIP EOCD comment size")
    if comment_size:
        raise ArtifactError("ZIP archive comments are unsupported")
    if central_offset + central_size != eocd_offset:
        raise ArtifactError("inconsistent ZIP central directory offset/size")
    if central_offset < 0 or central_offset > eocd_offset:
        raise ArtifactError("out-of-range ZIP central directory")

    entries = []
    cursor = central_offset
    total = 0
    compressed_total = 0
    local_offsets = set()
    for _index in range(count):
        if cursor + _CENTRAL.size > eocd_offset:
            raise ArtifactError("truncated ZIP central directory")
        header = _CENTRAL.unpack(_zip_read_at(stream, cursor, _CENTRAL.size, "central header"))
        if header[0] != b"PK\x01\x02":
            raise ArtifactError("malformed ZIP central directory signature")
        flags, method = header[3], header[4]
        crc, compressed_size, size = header[7], header[8], header[9]
        name_size, extra_size, member_comment = header[10], header[11], header[12]
        disk_start, local_offset = header[13], header[16]
        external_attributes = header[15]
        variable_size = name_size + extra_size + member_comment
        variable = _zip_read_at(
            stream, cursor + _CENTRAL.size, variable_size, "central variable fields"
        )
        raw_name = variable[:name_size]
        extra = variable[name_size:name_size + extra_size]
        if disk_start:
            raise ArtifactError("multi-disk ZIP member is forbidden")
        if member_comment:
            raise ArtifactError("ZIP member comments are unsupported")
        if flags & 1:
            raise ArtifactError("encrypted ZIP entries are forbidden")
        if flags & 8:
            raise ArtifactError("ZIP data descriptors are unsupported")
        allowed_flags = _ZIP_UTF8 | (_ZIP_DEFLATE_OPTIONS if method == zipfile.ZIP_DEFLATED else 0)
        if flags & ~allowed_flags:
            raise ArtifactError("unsafe or unsupported ZIP flags")
        if method not in _ZIP_METHODS:
            raise ArtifactError("unsupported wheel ZIP compression method")
        if method == zipfile.ZIP_STORED and compressed_size != size:
            raise ArtifactError("stored ZIP member has inconsistent sizes")
        _zip_extra_fields(extra)
        if local_offset in local_offsets:
            raise ArtifactError("duplicate ZIP local header offset")
        local_offsets.add(local_offset)
        name = _zip_name(raw_name, flags)
        directory = name.endswith("/")
        components = _zip_components(name, directory)
        unix_mode = (external_attributes >> 16) & 0xFFFF
        file_type = stat.S_IFMT(unix_mode)
        allowed_type = stat.S_IFDIR if directory else stat.S_IFREG
        if file_type not in (0, allowed_type):
            raise ArtifactError("ZIP symlink or special file type is forbidden")
        if bool(external_attributes & 0x10) != directory:
            raise ArtifactError("ZIP directory external attributes disagree with name")
        if directory and (size or crc):
            raise ArtifactError("ZIP directory has file data")
        if size > MAX_FILE:
            raise ArtifactError("ZIP member file size exceeds limit")
        total += size
        compressed_total += compressed_size
        if total > MAX_TOTAL:
            raise ArtifactError("ZIP total uncompressed size exceeds limit")
        if _zip_entry_ratio_exceeded(size, compressed_size):
            raise ArtifactError("ZIP per-entry compression ratio exceeds limit")

        local_raw = _zip_read_at(stream, local_offset, _LOCAL.size, "local header")
        local = _LOCAL.unpack(local_raw)
        if local[0] != b"PK\x03\x04":
            raise ArtifactError("malformed or out-of-range ZIP local header")
        local_name_size, local_extra_size = local[9], local[10]
        local_variable = _zip_read_at(
            stream,
            local_offset + _LOCAL.size,
            local_name_size + local_extra_size,
            "local variable fields",
        )
        local_name = local_variable[:local_name_size]
        local_extra = local_variable[local_name_size:]
        _zip_extra_fields(local_extra)
        if (
            local[2] != flags
            or local[3] != method
            or local[6] != crc
            or local[7] != compressed_size
            or local[8] != size
            or local_name != raw_name
        ):
            raise ArtifactError("ZIP local/central header disagreement")
        payload_offset = local_offset + _LOCAL.size + local_name_size + local_extra_size
        end = payload_offset + compressed_size
        if local_offset < 0 or end > central_offset:
            raise ArtifactError("ZIP local header or payload range overlaps central directory")
        entries.append({
            "raw_name": raw_name,
            "zip_name": name,
            "components": components,
            "directory": directory,
            "size": size,
            "compressed_size": compressed_size,
            "crc": crc,
            "method": method,
            "flags": flags,
            "local_offset": local_offset,
            "payload_offset": payload_offset,
            "range": (local_offset, end),
        })
        cursor += _CENTRAL.size + variable_size
    if cursor != eocd_offset or cursor - central_offset != central_size:
        raise ArtifactError("inconsistent ZIP central directory size/count")
    if total and (not compressed_total or total > MAX_RATIO * compressed_total):
        raise ArtifactError("ZIP aggregate compression ratio exceeds limit")
    ordered = sorted((entry["range"], entry["zip_name"]) for entry in entries)
    expected = 0
    for (start, end), _name in ordered:
        if start != expected:
            raise ArtifactError("overlapping or ambiguous ZIP local/payload regions")
        if end < start:
            raise ArtifactError("invalid ZIP payload range")
        expected = end
    if expected != central_offset:
        raise ArtifactError("ZIP local regions overlap or leave ambiguous data")
    return _zip_graph(entries), total


def _read_zip_members(stream, entries, retained_paths):
    records = {}
    retained = {}
    actual_total = 0
    compressed_total = 0
    try:
        for entry in entries:
            keep = entry["path"] in retained_paths
            if keep and entry["size"] > MAX_CONTROL_FILE:
                raise ArtifactError("wheel control file exceeds in-memory validation limit")
            chunks = [] if keep else None
            digest = hashlib.sha256()
            crc = 0
            actual_size = 0

            def consume(output):
                nonlocal actual_size, crc
                if not output:
                    return
                actual_size += len(output)
                if actual_size > entry["size"]:
                    raise ArtifactError("ZIP stream output exceeds declared size")
                if actual_size > MAX_FILE:
                    raise ArtifactError("ZIP stream output exceeds file size limit")
                if actual_total + actual_size > MAX_TOTAL:
                    raise ArtifactError("ZIP actual total output exceeds size limit")
                if entry["compressed_size"] and actual_size > MAX_RATIO * entry["compressed_size"]:
                    raise ArtifactError("ZIP actual stream compression ratio exceeds limit")
                digest.update(output)
                crc = zlib.crc32(output, crc)
                if chunks is not None:
                    chunks.append(output)

            stream.seek(entry["payload_offset"])
            compressed_remaining = entry["compressed_size"]
            if entry["method"] == zipfile.ZIP_STORED:
                if compressed_remaining != entry["size"]:
                    raise ArtifactError("stored ZIP stream size disagreement")
                while compressed_remaining:
                    chunk = stream.read(min(1024 * 1024, compressed_remaining))
                    if not chunk:
                        raise ArtifactError("truncated stored ZIP stream")
                    compressed_remaining -= len(chunk)
                    consume(chunk)
            elif entry["method"] == zipfile.ZIP_DEFLATED:
                decompressor = zlib.decompressobj(-15)
                while compressed_remaining:
                    compressed = stream.read(min(1024 * 1024, compressed_remaining))
                    if not compressed:
                        raise ArtifactError("truncated deflate ZIP payload")
                    compressed_remaining -= len(compressed)
                    pending = compressed
                    while pending:
                        output_limit = min(
                            _ZIP_OUTPUT_CHUNK,
                            entry["size"] - actual_size + 1,
                            MAX_FILE - actual_size + 1,
                            MAX_TOTAL - actual_total - actual_size + 1,
                        )
                        if entry["compressed_size"]:
                            output_limit = min(
                                output_limit,
                                MAX_RATIO * entry["compressed_size"] - actual_size + 1,
                            )
                        if output_limit <= 0:
                            raise ArtifactError("deflate ZIP output exceeds resource limit")
                        output = decompressor.decompress(pending, output_limit)
                        consume(output)
                        if decompressor.unused_data:
                            raise ArtifactError("trailing or concatenated deflate ZIP stream")
                        pending = decompressor.unconsumed_tail
                        if decompressor.eof and (pending or compressed_remaining):
                            raise ArtifactError("trailing or concatenated deflate ZIP stream")
                if not decompressor.eof:
                    raise ArtifactError("incomplete or truncated deflate ZIP stream")
                if decompressor.unconsumed_tail or decompressor.unused_data:
                    raise ArtifactError("trailing or incomplete deflate ZIP stream")
                flush_limit = min(
                    _ZIP_OUTPUT_CHUNK,
                    entry["size"] - actual_size + 1,
                    MAX_FILE - actual_size + 1,
                    MAX_TOTAL - actual_total - actual_size + 1,
                )
                if flush_limit <= 0:
                    raise ArtifactError("deflate ZIP flush exceeds resource limit")
                consume(decompressor.flush(flush_limit))
            else:
                raise ArtifactError("unsupported ZIP stream compression method")

            if actual_size != entry["size"]:
                raise ArtifactError("ZIP stream output size disagrees with declaration")
            if crc & 0xFFFFFFFF != entry["crc"]:
                raise ArtifactError("ZIP stream CRC disagrees with declaration")
            actual_total += actual_size
            compressed_total += entry["compressed_size"]
            if entry["directory"]:
                if actual_size:
                    raise ArtifactError("ZIP directory decompressed to file data")
            else:
                records[entry["path"]] = {
                    "size": actual_size,
                    "sha256": digest.digest(),
                }
                if chunks is not None:
                    retained[entry["path"]] = b"".join(chunks)
        if actual_total and (
            not compressed_total or actual_total > MAX_RATIO * compressed_total
        ):
            raise ArtifactError("ZIP actual aggregate compression ratio exceeds limit")
    except ArtifactError:
        raise
    except (OSError, EOFError, RuntimeError, zlib.error) as error:
        raise ArtifactError("malformed or invalid ZIP member stream") from error
    return records, retained


def _wheel_identity(path, entries):
    try:
        filename_name, filename_version, _build, tags = parse_wheel_filename(Path(path).name)
    except Exception as error:
        raise ArtifactError("malformed wheel filename") from error
    dist_infos = set()
    for entry in entries:
        components = entry["components"]
        if any(
            component.casefold().endswith(".dist-info")
            for component in components[1:]
        ):
            raise ArtifactError("nested .dist-info directory is forbidden")
        first = components[0]
        if first.casefold().endswith(".dist-info"):
            dist_infos.add(first)
    if len(dist_infos) != 1:
        raise ArtifactError("wheel must contain exactly one .dist-info directory")
    dist_info = dist_infos.pop()
    stem = dist_info[:-10]
    if "-" not in stem:
        raise ArtifactError("malformed wheel .dist-info identity")
    dist_name, dist_version_text = stem.rsplit("-", 1)
    try:
        dist_version = Version(dist_version_text)
    except InvalidVersion as error:
        raise ArtifactError("invalid wheel .dist-info version") from error
    if canonicalize_name(filename_name) != "chiq" or canonicalize_name(dist_name) != "chiq":
        raise ArtifactError("wheel filename/dist-info project name disagreement")
    if filename_version != dist_version:
        raise ArtifactError("wheel filename/dist-info version disagreement")
    return dist_info, dist_version, {str(tag) for tag in tags}


def _strict_metadata(data):
    try:
        message = BytesParser(policy=policy.default).parsebytes(data)
    except (TypeError, ValueError) as error:
        raise ArtifactError("malformed wheel METADATA") from error
    if message.defects:
        raise ArtifactError("malformed wheel METADATA")
    values = {}
    for name in ("Metadata-Version", "Name", "Version", "Requires-Python"):
        all_values = message.get_all(name, [])
        if len(all_values) != 1 or not str(all_values[0]).strip():
            raise ArtifactError("malformed or duplicate wheel metadata %s" % name)
        values[name] = str(all_values[0]).strip()
    if values["Metadata-Version"] not in _SUPPORTED_CORE_METADATA:
        raise ArtifactError("unexpected wheel metadata version")
    try:
        python = SpecifierSet(values["Requires-Python"])
    except InvalidSpecifier as error:
        raise ArtifactError("malformed Requires-Python metadata") from error
    if str(python) != ">=3.8":
        raise ArtifactError("wheel Requires-Python must be >=3.8")
    extras = message.get_all("Provides-Extra", [])
    if len(extras) != len(set(extras)) or set(map(str, extras)) != set(_EXTRA_REQUIREMENTS):
        raise ArtifactError("wheel dependency extra declarations disagree with pyproject")
    expected_text = list(_CORE_REQUIREMENTS)
    expected_text.extend(
        '%s; extra == "%s"' % (requirement, extra)
        for extra, requirements in _EXTRA_REQUIREMENTS.items()
        for requirement in requirements
    )
    try:
        actual = [Requirement(str(value)) for value in message.get_all("Requires-Dist", [])]
        expected = [Requirement(value) for value in expected_text]
    except InvalidRequirement as error:
        raise ArtifactError("malformed wheel Requires-Dist metadata") from error
    if len(actual) != len(set(actual)) or set(actual) != set(expected):
        raise ArtifactError("wheel dependency declarations disagree with pyproject")
    try:
        version = Version(values["Version"])
    except InvalidVersion as error:
        raise ArtifactError("invalid wheel METADATA version") from error
    return values["Name"], version


def _validate_entry_points(data):
    try:
        text = data.decode("utf-8", errors="strict")
        if re.search(r"(?mi)^\s*\[DEFAULT\]\s*$", text):
            raise ArtifactError("wheel entry points forbid [DEFAULT]")
        parser = configparser.ConfigParser(
            interpolation=None, strict=True, delimiters=("=",), comment_prefixes=()
        )
        parser.optionxform = str
        parser.read_string(text)
    except (UnicodeDecodeError, configparser.Error, ValueError) as error:
        raise ArtifactError("malformed or duplicate wheel entry points") from error
    if parser.sections() != ["console_scripts"]:
        raise ArtifactError("wheel entry points must contain exactly [console_scripts]")
    actual = {name: value.strip() for name, value in parser.items("console_scripts")}
    if actual != _CONSOLE_SCRIPTS:
        raise ArtifactError("wheel console script entry points disagree with pyproject")


def _valid_native_module_name(basename):
    return _NATIVE_MODULE_NAME.fullmatch(basename) is not None


def _validate_wheel_control(data, filename_tags):
    try:
        message = BytesParser(policy=policy.default).parsebytes(data)
    except (TypeError, ValueError) as error:
        raise ArtifactError("malformed WHEEL metadata") from error
    if message.defects:
        raise ArtifactError("malformed WHEEL metadata")
    values = {}
    for name in ("Wheel-Version", "Root-Is-Purelib"):
        all_values = message.get_all(name, [])
        if len(all_values) != 1 or not str(all_values[0]).strip():
            raise ArtifactError("malformed or duplicate WHEEL %s" % name)
        values[name] = str(all_values[0]).strip()
    if values["Wheel-Version"] != "1.0":
        raise ArtifactError("unsupported WHEEL version")
    if values["Root-Is-Purelib"].casefold() != "false":
        raise ArtifactError("native ChiQ WHEEL Root-Is-Purelib must be false")
    tags = [str(value).strip() for value in message.get_all("Tag", [])]
    if not tags or len(tags) != len(set(tags)) or set(tags) != filename_tags:
        raise ArtifactError("WHEEL tags disagree with wheel filename")


def _validate_record(records, data, dist_info):
    record_path = dist_info + "/RECORD"
    for signature in (dist_info + "/RECORD.jws", dist_info + "/RECORD.p7s"):
        if signature in records:
            raise ArtifactError("wheel RECORD signatures are forbidden")
    try:
        text = data.decode("utf-8", errors="strict")
        rows = list(csv.reader(io.StringIO(text, newline=""), strict=True))
    except (UnicodeDecodeError, csv.Error) as error:
        raise ArtifactError("wheel RECORD must be strict UTF-8 CSV") from error
    seen = set()
    aliases = {}
    for row in rows:
        if len(row) != 3:
            raise ArtifactError("wheel RECORD rows must have exactly three columns")
        name, encoded_hash, encoded_size = row
        components = canonical_member_path(name)
        alias = tuple(_identity_component(item) for item in components)
        previous = aliases.get(alias)
        if previous is not None and previous != name:
            raise ArtifactError("wheel RECORD path alias/collision")
        aliases[alias] = name
        if name in seen:
            raise ArtifactError("duplicate wheel RECORD row")
        seen.add(name)
        if name == record_path:
            if encoded_hash or encoded_size:
                raise ArtifactError("wheel RECORD row must have empty hash and size")
            continue
        if name not in records:
            raise ArtifactError("extra wheel RECORD row")
        if not re.fullmatch(r"sha256=[A-Za-z0-9_-]{43}", encoded_hash):
            raise ArtifactError("wheel RECORD hash must be canonical unpadded URL-safe sha256")
        digest_text = encoded_hash[len("sha256="):]
        try:
            digest = base64.b64decode(digest_text + "=", altchars=b"-_", validate=True)
        except (ValueError, TypeError) as error:
            raise ArtifactError("malformed wheel RECORD base64 hash") from error
        expected_digest = records[name]["sha256"]
        canonical_digest = base64.urlsafe_b64encode(digest).rstrip(b"=").decode("ascii")
        if canonical_digest != digest_text:
            raise ArtifactError("wheel RECORD base64 hash is not canonical")
        if digest != expected_digest:
            raise ArtifactError("wheel RECORD hash mismatch")
        if not re.fullmatch(r"0|[1-9][0-9]*", encoded_size):
            raise ArtifactError("wheel RECORD size is not canonical decimal")
        if int(encoded_size) != records[name]["size"]:
            raise ArtifactError("wheel RECORD size mismatch")
    if seen != set(records):
        raise ArtifactError("wheel RECORD has missing or extra file rows")


def _validate_wheel_content(records, dist_info, entries):
    required = (
        "chiq/__init__.py",
        "chiq/cli/__init__.py",
        "chiq/cli/_deprecated.py",
        "chiq/point_group_data/__init__.py",
        "chiq/solver/__init__.py",
        "chiq/solver/cpp.py",
        "bse/__init__.py",
        "bse/point_group_data/__init__.py",
        "bse_solver/__init__.py",
        dist_info + "/METADATA",
        dist_info + "/WHEEL",
        dist_info + "/entry_points.txt",
        dist_info + "/RECORD",
    )
    for path in required:
        if path not in records:
            raise ArtifactError("missing required wheel package/content: %s" % path)
    contract_modules = set(
        ["chiq/cli/%s.py" % name for name in CANONICAL_COMMANDS]
        + ["chiq/point_group_data/%s.py" % name for name in POINT_GROUP_MODULES]
        + ["chiq/%s.py" % name for name in SHIM_MODULES]
        + ["bse/%s.py" % name for name in SHIM_MODULES]
        + ["bse/point_group_data/%s.py" % name for name in POINT_GROUP_MODULES]
        + ["chiq/cli/_common.py", "chiq/cli/_deprecated.py"]
        + [
            "chiq/solver/%s.py" % name
            for name in ("base", "cpp", "numpy", "cupy", "kernels", "layout")
        ]
    )
    missing_modules = contract_modules - set(records)
    if missing_modules:
        raise ArtifactError(
            "missing required wheel contract module: %s" % sorted(missing_modules)[0]
        )

    native = []
    for entry in entries:
        path = entry["path"]
        components = path.split("/")
        folded = [_identity_component(component) for component in components]
        basename = folded[-1]
        original_basename = components[-1]
        if not entry["directory"] and basename.endswith(NATIVE_SUFFIXES):
            if (
                len(components) == 2
                and components[0] == "chiq"
                and _valid_native_module_name(original_basename)
            ):
                native.append(path)
            else:
                raise ArtifactError("forbidden native/solver extension in wheel: %s" % path)
        if "__pycache__" in folded or basename.endswith((".pyc", ".pyo")):
            raise ArtifactError("forbidden cache content in wheel: %s" % path)
        if not entry["directory"] and basename.endswith((".h5", ".hdf5")):
            raise ArtifactError("forbidden HDF5 fixture/content in wheel: %s" % path)
        if components[0] not in ("chiq", "bse", "bse_solver", dist_info):
            raise ArtifactError("forbidden source/build/legacy wheel content: %s" % path)
        if any(component in ("tests", "test", "src", "source", "build", "cmake", "lib", "python", "scripts") for component in folded[:-1]):
            raise ArtifactError("forbidden source/build/legacy wheel layout: %s" % path)
        if basename in ("cmakelists.txt", "chiqvars.sh", "install_manifest.txt") or basename.endswith((".cpp", ".cc", ".cxx", ".h", ".hpp", ".o", ".a")):
            raise ArtifactError("forbidden source/build material in wheel: %s" % path)
    if len(native) != 1:
        raise ArtifactError("wheel must contain exactly one native chiq/_bse_solver extension")


def _close_private_wheel(private):
    primary_error = sys.exc_info()[1]
    try:
        private.close()
    except BaseException as error:
        if primary_error is None:
            raise ArtifactError("cannot close private wheel artifact resource") from error


def validate_wheel(path):
    """Validate ZIP safety, wheel metadata, RECORD, and ChiQ content."""
    archive_size = _checked_archive_size(path)
    private = None
    try:
        private = _copy_private_artifact(path, archive_size)
        entries, total = _preflight_zip(private, archive_size)
        dist_info, version, filename_tags = _wheel_identity(path, entries)
        metadata_path = dist_info + "/METADATA"
        wheel_path = dist_info + "/WHEEL"
        entry_points_path = dist_info + "/entry_points.txt"
        record_path = dist_info + "/RECORD"
        controls = (metadata_path, wheel_path, entry_points_path, record_path)
        records, retained = _read_zip_members(private, entries, set(controls))
        for required in controls:
            if required not in records:
                raise ArtifactError("missing required wheel dist-info file: %s" % required)
        metadata_name, metadata_version = _strict_metadata(retained[metadata_path])
        if canonicalize_name(metadata_name) != "chiq" or metadata_version != version:
            raise ArtifactError("wheel filename/dist-info/METADATA name or version disagreement")
        _validate_wheel_control(retained[wheel_path], filename_tags)
        _validate_entry_points(retained[entry_points_path])
        _validate_record(records, retained[record_path], dist_info)
        _validate_wheel_content(records, dist_info, entries)
    except ArtifactError:
        raise
    except (OSError, EOFError, UnicodeError, ValueError, struct.error) as error:
        raise ArtifactError("malformed wheel ZIP archive") from error
    finally:
        if private is not None:
            _close_private_wheel(private)
    return {
        "kind": "wheel",
        "name": "chiq",
        "version": str(version),
        "archive_size": archive_size,
        "uncompressed_size": total,
        "members": sorted(
            [
                {
                    "path": entry["path"],
                    "type": "directory" if entry["directory"] else "file",
                    "size": entry["size"],
                }
                for entry in entries
            ],
            key=lambda row: row["path"],
        ),
    }


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
    wheel = subparsers.add_parser("wheel")
    wheel.add_argument("archive")
    wheel.add_argument("--manifest")
    args = parser.parse_args(argv)

    extract = getattr(args, "extract", None)
    if args.manifest and extract:
        manifest_path = _resolve_path(args.manifest)
        extraction_path = _resolve_path(extract)
        try:
            manifest_path.relative_to(extraction_path)
        except ValueError:
            pass
        else:
            parser.error("--manifest must be outside --extract")
    if args.kind == "wheel":
        manifest = validate_wheel(args.archive)
    elif extract:
        manifest = extract_sdist(args.archive, extract)
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
