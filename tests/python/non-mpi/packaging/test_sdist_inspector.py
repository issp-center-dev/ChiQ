import io
import bz2
import gzip
import json
import lzma
import os
from pathlib import Path
import stat
import subprocess
import sys
import tarfile

import pytest

import artifact_inspector


VERSION = "1.2.3"
ROOT = "chiq-" + VERSION


def required_files(version=VERSION):
    pyproject = """\
[project]
name = "chiq"
dynamic = ["version"]

[tool.scikit-build.metadata.version]
provider = "scikit_build_core.metadata.regex"
input = "python/package/chiq/__init__.py"
regex = '''__version__\\s*=\\s*["'](?P<value>[^"']+)["']'''
"""
    return {
        "pyproject.toml": pyproject.encode(),
        "PKG-INFO": (
            "Metadata-Version: 2.3\nName: chiq\nVersion: %s\n\n" % version
        ).encode(),
        "LICENSE": b"license\n",
        "README.md": b"readme\n",
        "CMakeLists.txt": b"project(ChiQ)\n",
        "python/CMakeLists.txt": b"# python\n",
        "src/CMakeLists.txt": b"# sources\n",
        "src/solver.cpp": b"int solve() { return 0; }\n",
        "src/solver.hpp": b"int solve();\n",
        "cmake/FindEigen3.cmake": b"# helper\n",
        "python/package/chiq/__init__.py": (
            '__version__ = "%s"\n' % version
        ).encode(),
        "python/package/chiq/cli/__init__.py": b"",
        "python/package/chiq/cli/chiq_main.py": b"",
        "python/package/chiq/point_group_data/__init__.py": b"",
        "python/package/chiq/point_group_data/C1.py": b"",
        "python/package/bse/__init__.py": b"",
        "python/package/bse_solver/__init__.py": b"",
        "tests/python/non-mpi/fixture/ref/data.h5": b"HDF5 fixture",
    }


def make_sdist(tmp_path, files=None, entries=(), filename=None, root=ROOT):
    if files is None:
        files = required_files()
    path = tmp_path / (filename or (root + ".tar.gz"))
    with tarfile.open(str(path), "w:gz") as archive:
        for relative, data in files.items():
            info = tarfile.TarInfo(root + "/" + relative)
            info.size = len(data)
            info.mode = 0o644
            archive.addfile(info, io.BytesIO(data))
        for info, data in entries:
            archive.addfile(info, io.BytesIO(data) if data is not None else None)
    return path


def tar_info(name, data=b"", type_=tarfile.REGTYPE, mode=0o644):
    info = tarfile.TarInfo(name)
    info.type = type_
    info.mode = mode
    info.size = len(data) if type_ in (tarfile.REGTYPE, tarfile.AREGTYPE) else 0
    return info, data


def checksum_header(header):
    header[148:156] = b"        "
    header[148:156] = ("%06o\0 " % sum(header)).encode("ascii")
    return bytes(header)


def raw_archive(path):
    return bytearray(gzip.decompress(path.read_bytes()))


def write_raw_archive(path, raw):
    path.write_bytes(gzip.compress(bytes(raw)))
    return path


def pax_record(key, value):
    body = (key + "=" + value + "\n").encode("utf-8")
    length = len(body) + 2
    while True:
        record = str(length).encode("ascii") + b" " + body
        if len(record) == length:
            return record
        length = len(record)


def extension_header(type_, payload):
    info = tarfile.TarInfo("PaxHeader")
    info.type = type_
    info.size = len(payload)
    header = info.tobuf(format=tarfile.USTAR_FORMAT)
    padding = b"\0" * ((-len(payload)) % tarfile.BLOCKSIZE)
    return header + payload + padding


def zero_raw_header_size(header):
    changed = bytearray(header)
    changed[124:136] = b"00000000000\0"
    return checksum_header(changed)


def replace_file(files, path, data):
    changed = dict(files)
    changed[path] = data
    return changed


def test_validates_without_extracting_and_accepts_implicit_root(tmp_path):
    archive = make_sdist(tmp_path)

    manifest = artifact_inspector.validate_sdist(archive)

    assert manifest["kind"] == "sdist"
    assert manifest["name"] == "chiq"
    assert manifest["version"] == VERSION
    assert manifest["root"] == ROOT
    assert manifest["members"] == sorted(manifest["members"], key=lambda row: row["path"])
    assert not (tmp_path / ROOT).exists()


def test_accepts_explicit_root_directory(tmp_path):
    archive = make_sdist(
        tmp_path,
        entries=[tar_info(ROOT, type_=tarfile.DIRTYPE)],
    )

    assert artifact_inspector.validate_sdist(archive)["root"] == ROOT


@pytest.mark.parametrize(
    "path",
    [
        "../escape",
        "/absolute",
        "C:/drive",
        "//server/share",
        "\\\\server\\share",
        "root\\child",
        "root//child",
        "root/./child",
        "root/../child",
        "root/",
        "root\x00child",
    ],
)
def test_canonical_path_rejects_unsafe_spelling(path):
    with pytest.raises(artifact_inspector.ArtifactError, match="path"):
        artifact_inspector.canonical_member_path(path)


@pytest.mark.parametrize(
    "start,end,description",
    [(0, 100, "name"), (157, 257, "linkname"), (345, 500, "prefix")],
)
def test_rejects_nonpadding_bytes_after_raw_path_field_nul(
    tmp_path, start, end, description
):
    archive = make_sdist(tmp_path)
    raw = raw_archive(archive)
    terminator = raw[start:end].index(0) + start
    raw[terminator + 1] = ord("X")
    raw[:512] = checksum_header(raw[:512])
    write_raw_archive(archive, raw)

    with pytest.raises(
        artifact_inspector.ArtifactError,
        match="padding|raw.*%s|ambiguous" % description,
    ):
        artifact_inspector.validate_sdist(archive)


def test_preflight_rejects_invalid_raw_header_checksum(tmp_path):
    archive = make_sdist(tmp_path)
    raw = raw_archive(archive)
    raw[100] ^= 1
    write_raw_archive(archive, raw)

    with pytest.raises(artifact_inspector.ArtifactError, match="checksum"):
        artifact_inspector.validate_sdist(archive)


def test_accepts_canonical_raw_header_nul_padding(tmp_path):
    archive = make_sdist(tmp_path)

    assert artifact_inspector.validate_sdist(archive)["root"] == ROOT


@pytest.mark.parametrize("key", ["path", "linkpath"])
def test_rejects_duplicate_pax_path_keys_before_tarfile_collapse(tmp_path, key):
    archive = make_sdist(tmp_path)
    raw = raw_archive(archive)
    safe = ROOT + "/pyproject.toml"
    payload = pax_record(key, "../unsafe") + pax_record(key, safe)
    write_raw_archive(archive, extension_header(tarfile.XHDTYPE, payload) + raw)

    with pytest.raises(artifact_inspector.ArtifactError, match="duplicate.*PAX|PAX.*duplicate"):
        artifact_inspector.validate_sdist(archive)


@pytest.mark.parametrize(
    "payload",
    [
        b"99 path=safe\n",
        b"4 x\n",
        b"x path=safe\n",
        b"12 path=\xffxx\n",
    ],
)
def test_rejects_malformed_pax_record_framing_length_or_utf8(tmp_path, payload):
    archive = make_sdist(tmp_path)
    raw = raw_archive(archive)
    write_raw_archive(archive, extension_header(tarfile.XHDTYPE, payload) + raw)

    with pytest.raises(artifact_inspector.ArtifactError, match="PAX.*(record|length|UTF-8)|malformed PAX"):
        artifact_inspector.validate_sdist(archive)


@pytest.mark.parametrize("key", ["path", "linkpath"])
def test_rejects_unsafe_value_in_every_pax_path_occurrence(tmp_path, key):
    archive = make_sdist(tmp_path)
    raw = raw_archive(archive)
    payload = pax_record(key, "../unsafe")
    write_raw_archive(archive, extension_header(tarfile.XHDTYPE, payload) + raw)

    with pytest.raises(artifact_inspector.ArtifactError, match="path"):
        artifact_inspector.validate_sdist(archive)


def test_rejects_conflicting_global_and_per_file_pax_metadata(tmp_path):
    archive = make_sdist(tmp_path)
    raw = raw_archive(archive)
    global_payload = pax_record("comment", "global")
    local_payload = pax_record("comment", "local")
    prefixed = (
        extension_header(tarfile.XGLTYPE, global_payload)
        + extension_header(tarfile.XHDTYPE, local_payload)
        + raw
    )
    write_raw_archive(archive, prefixed)

    with pytest.raises(artifact_inspector.ArtifactError, match="conflicting.*PAX|PAX.*conflict"):
        artifact_inspector.validate_sdist(archive)


@pytest.mark.parametrize("type_", [tarfile.XHDTYPE, tarfile.XGLTYPE])
def test_pax_effective_size_limit_precedes_raw_zero_payload(tmp_path, monkeypatch, type_):
    archive = tmp_path / (ROOT + ".tar.gz")
    payload = pax_record("size", str(artifact_inspector.MAX_FILE + 1))
    member = tarfile.TarInfo(ROOT + "/huge")
    member.size = 0
    raw = extension_header(type_, payload) + member.tobuf(format=tarfile.USTAR_FORMAT)
    raw += b"\0" * (2 * tarfile.BLOCKSIZE)
    write_raw_archive(archive, raw)
    monkeypatch.setattr(
        artifact_inspector.tarfile,
        "open",
        lambda *args, **kwargs: pytest.fail("over-limit effective size must preflight"),
    )

    with pytest.raises(artifact_inspector.ArtifactError, match="file.*size|size.*limit"):
        artifact_inspector.validate_sdist(archive)


@pytest.mark.parametrize("value", ["-1", "01", "+1", "1.0", "x", str(1 << 100)])
def test_rejects_noncanonical_or_unrepresentable_pax_size(tmp_path, value):
    archive = make_sdist(tmp_path)
    raw = raw_archive(archive)
    payload = pax_record("size", value)
    write_raw_archive(archive, extension_header(tarfile.XHDTYPE, payload) + raw)

    with pytest.raises(artifact_inspector.ArtifactError, match="PAX.*size|size.*PAX"):
        artifact_inspector.validate_sdist(archive)


def test_rejects_duplicate_pax_size_key(tmp_path):
    archive = make_sdist(tmp_path)
    raw = raw_archive(archive)
    payload = pax_record("size", "1") + pax_record("size", "2")
    write_raw_archive(archive, extension_header(tarfile.XHDTYPE, payload) + raw)

    with pytest.raises(artifact_inspector.ArtifactError, match="duplicate.*PAX"):
        artifact_inspector.validate_sdist(archive)


def test_rejects_conflicting_global_and_local_pax_size(tmp_path):
    archive = make_sdist(tmp_path)
    raw = raw_archive(archive)
    prefixed = (
        extension_header(tarfile.XGLTYPE, pax_record("size", "1"))
        + extension_header(tarfile.XHDTYPE, pax_record("size", "2"))
        + raw
    )
    write_raw_archive(archive, prefixed)

    with pytest.raises(artifact_inspector.ArtifactError, match="conflicting.*PAX"):
        artifact_inspector.validate_sdist(archive)


def test_accepts_safe_local_pax_size_override_with_raw_zero(tmp_path):
    archive = make_sdist(tmp_path)
    raw = raw_archive(archive)
    first_size = len(required_files()["pyproject.toml"])
    raw[:512] = zero_raw_header_size(raw[:512])
    prefixed = extension_header(
        tarfile.XHDTYPE,
        pax_record("size", str(first_size)),
    ) + raw
    write_raw_archive(archive, prefixed)

    assert artifact_inspector.validate_sdist(archive)["root"] == ROOT


def test_pending_size_override_does_not_resize_following_pax_carrier(tmp_path):
    archive = make_sdist(tmp_path)
    raw = raw_archive(archive)
    first_size = len(required_files()["pyproject.toml"])
    raw[:512] = zero_raw_header_size(raw[:512])
    prefixed = (
        extension_header(tarfile.XHDTYPE, pax_record("size", str(first_size)))
        + extension_header(tarfile.XHDTYPE, pax_record("comment", "safe"))
        + raw
    )
    write_raw_archive(archive, prefixed)

    assert artifact_inspector.validate_sdist(archive)["root"] == ROOT


@pytest.mark.parametrize(
    "first,second",
    [
        ("chiq-1.2.3/A.py", "chiq-1.2.3/a.py"),
        ("chiq-1.2.3/caf\N{LATIN SMALL LETTER E WITH ACUTE}/x", "chiq-1.2.3/cafe\N{COMBINING ACUTE ACCENT}/y"),
        ("chiq-1.2.3/Dir/x", "chiq-1.2.3/dir/y"),
    ],
)
def test_rejects_normalized_full_and_prefix_aliases(tmp_path, first, second):
    archive = make_sdist(
        tmp_path,
        entries=[tar_info(first, b"one"), tar_info(second, b"two")],
    )

    with pytest.raises(artifact_inspector.ArtifactError, match="collision|alias"):
        artifact_inspector.validate_sdist(archive)


def test_rejects_duplicate_entries(tmp_path):
    duplicate = ROOT + "/README.md"
    archive = make_sdist(tmp_path, entries=[tar_info(duplicate, b"again")])

    with pytest.raises(artifact_inspector.ArtifactError, match="duplicate"):
        artifact_inspector.validate_sdist(archive)


@pytest.mark.parametrize(
    "entries",
    [
        [tar_info(ROOT + "/conflict", b"file"), tar_info(ROOT + "/conflict", type_=tarfile.DIRTYPE)],
        [tar_info(ROOT + "/ancestor", b"file"), tar_info(ROOT + "/ancestor/child", b"child")],
        [tar_info(ROOT + "/ancestor/child", b"child"), tar_info(ROOT + "/ancestor", b"file")],
    ],
)
def test_rejects_file_directory_and_ancestor_file_conflicts(tmp_path, entries):
    archive = make_sdist(tmp_path, entries=entries)

    with pytest.raises(artifact_inspector.ArtifactError, match="conflict|ancestor|duplicate"):
        artifact_inspector.validate_sdist(archive)


@pytest.mark.parametrize(
    "type_",
    [
        tarfile.SYMTYPE,
        tarfile.LNKTYPE,
        tarfile.CHRTYPE,
        tarfile.BLKTYPE,
        tarfile.FIFOTYPE,
        tarfile.GNUTYPE_SPARSE,
        tarfile.CONTTYPE,
    ],
)
def test_rejects_links_sparse_and_special_members(tmp_path, type_):
    info, data = tar_info(ROOT + "/special", type_=type_)
    if type_ in (tarfile.SYMTYPE, tarfile.LNKTYPE):
        info.linkname = "target"
    archive = make_sdist(tmp_path, entries=[(info, data)])

    with pytest.raises(artifact_inspector.ArtifactError, match="regular|directory|special|sparse"):
        artifact_inspector.validate_sdist(archive)


def test_rejects_physical_archive_limit_before_tar_parsing(tmp_path, monkeypatch):
    archive = make_sdist(tmp_path)
    monkeypatch.setattr(artifact_inspector, "_physical_size", lambda path: artifact_inspector.MAX_ARCHIVE + 1)
    monkeypatch.setattr(
        artifact_inspector.tarfile,
        "open",
        lambda *args, **kwargs: pytest.fail("tar parser must not run"),
    )

    with pytest.raises(artifact_inspector.ArtifactError, match="archive.*size|512"):
        artifact_inspector.validate_sdist(archive)


def test_preflight_rejects_oversized_truncated_member_at_header(tmp_path, monkeypatch):
    archive = tmp_path / (ROOT + ".tar.gz")
    info = tarfile.TarInfo(ROOT + "/huge")
    info.size = artifact_inspector.MAX_FILE + 1
    archive.write_bytes(gzip.compress(info.tobuf(format=tarfile.USTAR_FORMAT)))
    monkeypatch.setattr(
        artifact_inspector,
        "_discard_exact",
        lambda *args, **kwargs: pytest.fail("oversized payload must not be consumed"),
        raising=False,
    )

    with pytest.raises(artifact_inspector.ArtifactError, match="file.*limit|size.*limit"):
        artifact_inspector.validate_sdist(archive)


def test_preflight_enforces_member_count_before_tarfile_getmembers(tmp_path, monkeypatch):
    archive = make_sdist(tmp_path)
    monkeypatch.setattr(artifact_inspector, "MAX_MEMBERS", 1)
    monkeypatch.setattr(
        artifact_inspector.tarfile,
        "open",
        lambda *args, **kwargs: pytest.fail("tarfile must not enumerate over-limit archive"),
    )

    with pytest.raises(artifact_inspector.ArtifactError, match="member"):
        artifact_inspector.validate_sdist(archive)


def test_preflight_enforces_ratio_before_payload_consumption(tmp_path, monkeypatch):
    archive = make_sdist(tmp_path)
    monkeypatch.setattr(artifact_inspector, "MAX_RATIO", 0)
    monkeypatch.setattr(
        artifact_inspector,
        "_discard_exact",
        lambda *args, **kwargs: pytest.fail("over-ratio payload must not be consumed"),
        raising=False,
    )
    monkeypatch.setattr(
        artifact_inspector.tarfile,
        "open",
        lambda *args, **kwargs: pytest.fail("tarfile must not parse over-ratio archive"),
    )

    with pytest.raises(artifact_inspector.ArtifactError, match="ratio"):
        artifact_inspector.validate_sdist(archive)


def test_rejects_checksum_valid_header_after_end_marker(tmp_path):
    archive = make_sdist(tmp_path)
    raw = raw_archive(archive)
    hidden = tarfile.TarInfo(ROOT + "/hidden")
    hidden.size = 0
    raw.extend(hidden.tobuf(format=tarfile.USTAR_FORMAT))
    write_raw_archive(archive, raw)

    with pytest.raises(artifact_inspector.ArtifactError, match="after.*end|trailing.*nonzero"):
        artifact_inspector.validate_sdist(archive)


def test_rejects_excessive_compressed_zero_padding(tmp_path):
    archive = make_sdist(tmp_path)
    raw = raw_archive(archive)
    raw.extend(b"\0" * (2 * tarfile.RECORDSIZE))
    write_raw_archive(archive, raw)

    with pytest.raises(artifact_inspector.ArtifactError, match="padding.*limit|excessive.*padding"):
        artifact_inspector.validate_sdist(archive)


@pytest.mark.parametrize(
    "compress",
    [lambda data: data, gzip.compress, bz2.compress, lzma.compress],
    ids=["plain", "gzip", "bzip2", "xz"],
)
def test_accepts_normal_record_padding_for_supported_tar_streams(tmp_path, compress):
    archive = make_sdist(tmp_path)
    raw = raw_archive(archive)
    archive.write_bytes(compress(raw))

    assert artifact_inspector.validate_sdist(archive)["root"] == ROOT


def test_accepts_resource_boundaries_and_rejects_overflow(tmp_path, monkeypatch):
    archive = make_sdist(tmp_path)
    count = len(required_files())
    logical = sum(len(value) for value in required_files().values())
    monkeypatch.setattr(artifact_inspector, "MAX_MEMBERS", count)
    monkeypatch.setattr(artifact_inspector, "MAX_FILE", max(map(len, required_files().values())))
    monkeypatch.setattr(artifact_inspector, "MAX_TOTAL", logical)
    monkeypatch.setattr(artifact_inspector, "MAX_RATIO", logical / archive.stat().st_size)
    artifact_inspector.validate_sdist(archive)

    for constant, value, message in [
        ("MAX_MEMBERS", count - 1, "member"),
        ("MAX_FILE", 1, "file"),
        ("MAX_TOTAL", logical - 1, "total"),
        ("MAX_RATIO", (logical / archive.stat().st_size) - 0.000001, "ratio"),
    ]:
        monkeypatch.setattr(artifact_inspector, constant, value)
        with pytest.raises(artifact_inspector.ArtifactError, match=message):
            artifact_inspector.validate_sdist(archive)
        monkeypatch.undo()


def test_rejects_logical_huge_member_without_allocating(tmp_path, monkeypatch):
    archive = make_sdist(tmp_path)
    real_open = artifact_inspector.tarfile.open

    class LogicalArchive:
        def __init__(self):
            self.archive = real_open(str(archive), "r:*")

        def close(self):
            self.archive.close()

        def getmembers(self):
            members = self.archive.getmembers()
            members[0].size = artifact_inspector.MAX_FILE + 1
            return members

    monkeypatch.setattr(artifact_inspector.tarfile, "open", lambda *args, **kwargs: LogicalArchive())
    with pytest.raises(artifact_inspector.ArtifactError, match="file"):
        artifact_inspector.validate_sdist(archive)


def test_invalid_complete_graph_writes_nothing(tmp_path):
    invalid = tar_info(ROOT + "/z-link", type_=tarfile.SYMTYPE)[0]
    invalid.linkname = "elsewhere"
    archive = make_sdist(tmp_path, entries=[(invalid, None)])
    destination = tmp_path / "extract"

    with pytest.raises(artifact_inspector.ArtifactError):
        artifact_inspector.extract_sdist(archive, destination)

    assert not destination.exists()


def test_source_mutation_during_destination_creation_cannot_change_extracted_bytes(
    tmp_path, monkeypatch
):
    archive = make_sdist(tmp_path)
    original = bytes(raw_archive(archive))
    archive.write_bytes(original)
    mutated = original.replace(b"readme\n", b"mutate\n", 1)
    assert len(mutated) == len(original)
    real_mkdir = artifact_inspector._mkdir
    changed = []

    def mutating_mkdir(path, *args, **kwargs):
        if not changed:
            archive.write_bytes(mutated)
            changed.append(True)
        return real_mkdir(path, *args, **kwargs)

    monkeypatch.setattr(artifact_inspector, "_mkdir", mutating_mkdir)
    destination = tmp_path / "extract"

    manifest = artifact_inspector.extract_sdist(archive, destination)

    assert manifest["version"] == VERSION
    assert (destination / ROOT / "README.md").read_bytes() == b"readme\n"


def test_source_mutation_between_raw_and_semantic_passes_is_ignored(tmp_path, monkeypatch):
    archive = make_sdist(tmp_path)
    original = bytes(raw_archive(archive))
    archive.write_bytes(original)
    mutated = original.replace(b"Version: 1.2.3", b"Version: 9.9.9", 1)
    assert len(mutated) == len(original)
    real_preflight = artifact_inspector._preflight_tar

    def mutating_preflight(raw_file, archive_size):
        result = real_preflight(raw_file, archive_size)
        archive.write_bytes(mutated)
        return result

    monkeypatch.setattr(artifact_inspector, "_preflight_tar", mutating_preflight)

    assert artifact_inspector.validate_sdist(archive)["version"] == VERSION


def test_extracts_with_sanitized_modes(tmp_path):
    files = required_files()
    executable = tar_info(ROOT + "/python/tool.py", b"#!/usr/bin/env python\n", mode=0o6771)
    nonexecutable = tar_info(ROOT + "/python/data.txt", b"data", mode=0o666)
    archive = make_sdist(tmp_path, files=files, entries=[executable, nonexecutable])
    destination = tmp_path / "extract"

    manifest = artifact_inspector.extract_sdist(archive, destination)

    assert manifest["root"] == ROOT
    assert stat.S_IMODE(destination.stat().st_mode) == 0o700
    assert stat.S_IMODE((destination / ROOT / "python").stat().st_mode) == 0o755
    assert stat.S_IMODE((destination / ROOT / "python" / "tool.py").stat().st_mode) == 0o755
    assert stat.S_IMODE((destination / ROOT / "python" / "data.txt").stat().st_mode) == 0o644


def test_root_permissions_are_applied_through_open_descriptor(tmp_path, monkeypatch):
    archive = make_sdist(tmp_path)
    monkeypatch.setattr(
        artifact_inspector.os,
        "chmod",
        lambda *args, **kwargs: pytest.fail("path-based chmod is race-prone"),
    )

    artifact_inspector.extract_sdist(archive, tmp_path / "extract")


def test_rejects_missing_posix_descriptor_capabilities(tmp_path, monkeypatch):
    archive = make_sdist(tmp_path)
    monkeypatch.setattr(
        artifact_inspector,
        "_capability_override",
        "O_NOFOLLOW unavailable",
        raising=False,
    )

    with pytest.raises(artifact_inspector.ArtifactError, match="capabil|O_NOFOLLOW"):
        artifact_inspector.extract_sdist(archive, tmp_path / "extract")

    assert not (tmp_path / "extract").exists()


def test_staging_root_replacement_does_not_chmod_or_write_attacker(tmp_path, monkeypatch):
    archive = make_sdist(tmp_path)
    attacker = tmp_path / "attacker"
    attacker.mkdir(mode=0o711)
    original_mode = stat.S_IMODE(attacker.stat().st_mode)

    def replace_root(parent_fd, name):
        if name.startswith(".extract.tmp-"):
            os.rename(name, name + ".moved", src_dir_fd=parent_fd, dst_dir_fd=parent_fd)
            os.symlink(str(attacker), name, dir_fd=parent_fd)

    monkeypatch.setattr(
        artifact_inspector,
        "_after_directory_mkdir",
        replace_root,
        raising=False,
    )

    with pytest.warns(RuntimeWarning, match="manual cleanup"):
        with pytest.raises(artifact_inspector.ArtifactError, match="identity|symlink|race|directory"):
            artifact_inspector.extract_sdist(archive, tmp_path / "extract")

    assert stat.S_IMODE(attacker.stat().st_mode) == original_mode
    assert not any(attacker.iterdir())
    assert not (tmp_path / "extract").exists()


def test_nested_directory_replacement_does_not_write_attacker(tmp_path, monkeypatch):
    archive = make_sdist(tmp_path)
    attacker = tmp_path / "attacker"
    attacker.mkdir()

    def replace_nested(parent_fd, name):
        if name == "python":
            os.rename(name, name + ".moved", src_dir_fd=parent_fd, dst_dir_fd=parent_fd)
            os.symlink(str(attacker), name, dir_fd=parent_fd)

    monkeypatch.setattr(
        artifact_inspector,
        "_after_directory_mkdir",
        replace_nested,
        raising=False,
    )

    with pytest.warns(RuntimeWarning, match="manual cleanup"):
        with pytest.raises(artifact_inspector.ArtifactError, match="identity|symlink|race|directory"):
            artifact_inspector.extract_sdist(archive, tmp_path / "extract")

    assert not any(attacker.iterdir())
    assert not (tmp_path / "extract").exists()


def test_write_failure_leaves_destination_absent_and_retryable(tmp_path, monkeypatch):
    archive = make_sdist(tmp_path)
    destination = tmp_path / "extract"
    real_write = os.write
    calls = []

    def failing_write(descriptor, data):
        calls.append(True)
        if len(calls) == 2:
            raise OSError("injected write failure")
        return real_write(descriptor, data)

    monkeypatch.setattr(artifact_inspector, "_write", failing_write, raising=False)
    with pytest.warns(RuntimeWarning, match="manual cleanup"):
        with pytest.raises(artifact_inspector.ArtifactError, match="write|extraction"):
            artifact_inspector.extract_sdist(archive, destination)
    assert not destination.exists()
    assert list(tmp_path.glob(".extract.tmp-*"))

    monkeypatch.setattr(artifact_inspector, "_write", real_write)
    artifact_inspector.extract_sdist(archive, destination)
    assert (destination / ROOT / "README.md").is_file()


def test_short_write_leaves_destination_absent(tmp_path, monkeypatch):
    archive = make_sdist(tmp_path)
    destination = tmp_path / "extract"
    monkeypatch.setattr(artifact_inspector, "_write", lambda descriptor, data: 0)

    with pytest.warns(RuntimeWarning, match="manual cleanup"):
        with pytest.raises(artifact_inspector.ArtifactError, match="short write"):
            artifact_inspector.extract_sdist(archive, destination)

    assert not destination.exists()
    assert list(tmp_path.glob(".extract.tmp-*"))


def test_nested_file_replacement_is_retained_on_failed_cleanup(tmp_path, monkeypatch):
    archive = make_sdist(tmp_path)
    destination = tmp_path / "extract"
    replacement = []

    def replace_then_fail(descriptor, data):
        staging = next(tmp_path.glob(".extract.tmp-*"))
        target = next(path for path in staging.rglob("*") if path.is_file())
        owned = target.with_name(target.name + ".owned")
        target.rename(owned)
        target.write_bytes(b"attacker replacement")
        replacement.extend([owned, target])
        raise artifact_inspector.ArtifactError("primary write failure")

    monkeypatch.setattr(artifact_inspector, "_write", replace_then_fail)
    with pytest.warns(RuntimeWarning, match="manual cleanup"):
        with pytest.raises(artifact_inspector.ArtifactError, match="primary write failure"):
            artifact_inspector.extract_sdist(archive, destination)

    assert not destination.exists()
    assert replacement[0].exists()
    assert replacement[1].read_bytes() == b"attacker replacement"


def test_swap_during_publish_cannot_return_false_success_or_touch_replacement(
    tmp_path, monkeypatch
):
    archive = make_sdist(tmp_path)
    destination = tmp_path / "extract"
    real_publish = artifact_inspector._publish_noreplace
    replacement_mode = 0o711

    def swapping_publish(parent_fd, source, target):
        os.rename(source, source + ".owned", src_dir_fd=parent_fd, dst_dir_fd=parent_fd)
        os.mkdir(source, replacement_mode, dir_fd=parent_fd)
        replacement_fd = os.open(source, os.O_RDONLY | os.O_DIRECTORY, dir_fd=parent_fd)
        try:
            marker_fd = os.open(
                "marker",
                os.O_WRONLY | os.O_CREAT | os.O_EXCL,
                0o600,
                dir_fd=replacement_fd,
            )
            os.write(marker_fd, b"attacker")
            os.close(marker_fd)
        finally:
            os.close(replacement_fd)
        return real_publish(parent_fd, source, target)

    monkeypatch.setattr(artifact_inspector, "_publish_noreplace", swapping_publish)

    with pytest.raises(artifact_inspector.ArtifactError, match="identity|false success|published"):
        artifact_inspector.extract_sdist(archive, destination)

    assert (destination / "marker").read_bytes() == b"attacker"
    assert stat.S_IMODE(destination.stat().st_mode) == replacement_mode


def test_cleanup_leaves_mismatched_staging_entry_untouched(tmp_path, monkeypatch):
    archive = make_sdist(tmp_path)
    destination = tmp_path / "extract"
    replacement = []

    def swap_before_identity(parent_fd, staging, target):
        os.rename(staging, staging + ".owned", src_dir_fd=parent_fd, dst_dir_fd=parent_fd)
        os.mkdir(staging, 0o711, dir_fd=parent_fd)
        replacement.append(staging)

    monkeypatch.setattr(
        artifact_inspector,
        "_before_publish_identity",
        swap_before_identity,
        raising=False,
    )

    with pytest.warns(RuntimeWarning, match="manual cleanup"):
        with pytest.raises(artifact_inspector.ArtifactError, match="identity|manual cleanup"):
            artifact_inspector.extract_sdist(archive, destination)

    replacement_path = tmp_path / replacement[0]
    assert replacement_path.is_dir()
    assert stat.S_IMODE(replacement_path.stat().st_mode) == 0o711
    assert not destination.exists()


def test_post_publish_close_failure_does_not_reverse_committed_success(
    tmp_path, monkeypatch
):
    archive = make_sdist(tmp_path)
    destination = tmp_path / "extract"
    real_publish = artifact_inspector._publish_noreplace
    real_close = os.close
    committed = []

    def tracking_publish(parent_fd, source, target):
        result = real_publish(parent_fd, source, target)
        committed.append(True)
        return result

    def failing_close(descriptor):
        real_close(descriptor)
        if committed:
            raise OSError("post-publish close failure")

    monkeypatch.setattr(artifact_inspector, "_publish_noreplace", tracking_publish)
    monkeypatch.setattr(artifact_inspector, "_close", failing_close)

    manifest = artifact_inspector.extract_sdist(archive, destination)

    assert manifest["root"] == ROOT
    assert (destination / ROOT / "README.md").is_file()


def test_directory_chain_close_failure_closes_new_descriptor(tmp_path, monkeypatch):
    root = tmp_path / "root"
    (root / "child").mkdir(parents=True)
    root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
    real_open = artifact_inspector._open
    real_close = os.close
    opened = []
    close_attempts = []

    def tracking_open(path, flags, *args, **kwargs):
        descriptor = real_open(path, flags, *args, **kwargs)
        if path == "child":
            opened.append(descriptor)
        return descriptor

    def failing_close(descriptor):
        close_attempts.append(descriptor)
        real_close(descriptor)
        raise OSError("injected close failure")

    monkeypatch.setattr(artifact_inspector, "_open", tracking_open)
    monkeypatch.setattr(artifact_inspector, "_close", failing_close)
    try:
        with pytest.raises(artifact_inspector.ArtifactError, match="close"):
            artifact_inspector._open_directory_chain(root_fd, ("child",))
    finally:
        os.close(root_fd)

    assert opened
    for descriptor in opened:
        with pytest.raises(OSError) as caught:
            os.fstat(descriptor)
        assert caught.value.errno == artifact_inspector.errno.EBADF
    assert all(close_attempts.count(descriptor) == 1 for descriptor in set(close_attempts))


def test_fchmod_failure_leaves_destination_absent_and_retryable(tmp_path, monkeypatch):
    archive = make_sdist(tmp_path)
    destination = tmp_path / "extract"
    real_fchmod = os.fchmod
    failed = []

    def failing_fchmod(descriptor, mode):
        if mode == 0o644 and not failed:
            failed.append(True)
            raise OSError("injected fchmod failure")
        return real_fchmod(descriptor, mode)

    monkeypatch.setattr(artifact_inspector, "_fchmod", failing_fchmod)
    with pytest.warns(RuntimeWarning, match="manual cleanup"):
        with pytest.raises(artifact_inspector.ArtifactError, match="extraction|file|create"):
            artifact_inspector.extract_sdist(archive, destination)
    assert not destination.exists()
    assert list(tmp_path.glob(".extract.tmp-*"))

    monkeypatch.setattr(artifact_inspector, "_fchmod", real_fchmod)
    artifact_inspector.extract_sdist(archive, destination)


def test_fstat_failure_leaves_destination_absent(tmp_path, monkeypatch):
    archive = make_sdist(tmp_path)
    destination = tmp_path / "extract"
    real_fstat = os.fstat
    calls = []

    def failing_fstat(descriptor):
        calls.append(True)
        if len(calls) == 2:
            raise OSError("injected fstat failure")
        return real_fstat(descriptor)

    monkeypatch.setattr(artifact_inspector, "_fstat", failing_fstat)
    with pytest.raises(artifact_inspector.ArtifactError, match="destination|open|safe|malformed"):
        artifact_inspector.extract_sdist(archive, destination)
    assert not destination.exists()


def test_cleanup_diagnostic_cannot_mask_primary_when_warnings_are_errors(
    tmp_path, monkeypatch
):
    archive = make_sdist(tmp_path)
    destination = tmp_path / "extract"

    def primary_failure(descriptor, data):
        raise artifact_inspector.ArtifactError("primary sentinel")

    monkeypatch.setattr(artifact_inspector, "_write", primary_failure)
    with artifact_inspector.warnings.catch_warnings():
        artifact_inspector.warnings.simplefilter("error", RuntimeWarning)
        with pytest.raises(artifact_inspector.ArtifactError, match="primary sentinel"):
            artifact_inspector.extract_sdist(archive, destination)

    assert not destination.exists()


def test_close_helper_attempts_every_descriptor_without_masking_primary(monkeypatch):
    descriptors = os.pipe()
    real_close = os.close
    closed = []

    def failing_close(descriptor):
        closed.append(descriptor)
        real_close(descriptor)
        raise OSError("injected close failure")

    monkeypatch.setattr(artifact_inspector, "_close", failing_close)
    with pytest.raises(artifact_inspector.ArtifactError, match="close"):
        artifact_inspector._close_descriptors(descriptors)
    assert closed == list(descriptors)

    second = os.pipe()
    closed[:] = []
    with pytest.raises(ValueError, match="primary"):
        try:
            raise ValueError("primary")
        finally:
            artifact_inspector._close_descriptors(second)
    assert closed == list(second)


def test_validation_explicitly_closes_tar_member_streams(tmp_path, monkeypatch):
    archive = make_sdist(tmp_path)
    real_extractfile = tarfile.TarFile.extractfile
    opened = []
    closed = []

    class Proxy:
        def __init__(self, stream):
            self.stream = stream
            opened.append(self)

        def read(self, size=-1):
            return self.stream.read(size)

        def close(self):
            closed.append(self)
            self.stream.close()

    def tracking_extractfile(instance, member):
        stream = real_extractfile(instance, member)
        return None if stream is None else Proxy(stream)

    monkeypatch.setattr(tarfile.TarFile, "extractfile", tracking_extractfile)
    artifact_inspector.validate_sdist(archive)

    assert opened
    assert closed == opened


def test_requires_new_destination_and_rejects_symlink_collision(tmp_path):
    archive = make_sdist(tmp_path)
    destination = tmp_path / "extract"
    destination.symlink_to(tmp_path / "elsewhere")

    with pytest.raises(artifact_inspector.ArtifactError, match="destination|symlink|exist"):
        artifact_inspector.extract_sdist(archive, destination)

    assert not (tmp_path / "elsewhere").exists()


def test_does_not_overwrite_preexisting_empty_destination(tmp_path):
    archive = make_sdist(tmp_path)
    destination = tmp_path / "extract"
    destination.mkdir()

    with pytest.raises(artifact_inspector.ArtifactError, match="destination|exist"):
        artifact_inspector.extract_sdist(archive, destination)

    assert destination.is_dir()
    assert not any(destination.iterdir())


def test_file_creation_is_exclusive_and_nofollow(tmp_path, monkeypatch):
    archive = make_sdist(tmp_path)
    real_open = artifact_inspector._open
    checked = []

    def checking_open(path, flags, *args, **kwargs):
        if flags & artifact_inspector.os.O_CREAT:
            checked.append(flags)
            assert flags & artifact_inspector.os.O_EXCL
            assert flags & artifact_inspector.os.O_NOFOLLOW
        return real_open(path, flags, *args, **kwargs)

    monkeypatch.setattr(artifact_inspector, "_open", checking_open)
    artifact_inspector.extract_sdist(archive, tmp_path / "extract")

    assert checked


def test_resolved_containment_rejects_nested_symlink(tmp_path):
    root = tmp_path / "root"
    outside = tmp_path / "outside"
    root.mkdir()
    outside.mkdir()
    (root / "link").symlink_to(outside, target_is_directory=True)

    with pytest.raises(artifact_inspector.ArtifactError, match="escapes"):
        artifact_inspector._require_contained(root, root / "link" / "file")


def test_resolved_containment_handles_tmp_aliases(tmp_path, monkeypatch):
    root = Path("/private/tmp/extract")
    candidate = Path("/tmp/extract") / ROOT / "README.md"
    real_resolve = artifact_inspector._resolve_path

    def macos_resolve(path):
        value = Path(path)
        if str(value).startswith("/tmp/"):
            return Path("/private" + str(value))
        return real_resolve(value)

    monkeypatch.setattr(artifact_inspector, "_resolve_path", macos_resolve)
    artifact_inspector._require_contained(root, candidate)


@pytest.mark.parametrize(
    "missing",
    [
        "pyproject.toml",
        "LICENSE",
        "README.md",
        "CMakeLists.txt",
        "python/CMakeLists.txt",
        "src/CMakeLists.txt",
        "src/solver.cpp",
        "src/solver.hpp",
        "cmake/FindEigen3.cmake",
        "python/package/chiq/cli/chiq_main.py",
        "python/package/chiq/point_group_data/C1.py",
        "python/package/bse/__init__.py",
        "python/package/bse_solver/__init__.py",
        "tests/python/non-mpi/fixture/ref/data.h5",
    ],
)
def test_rejects_missing_required_sdist_content(tmp_path, missing):
    files = required_files()
    del files[missing]
    archive = make_sdist(tmp_path, files=files)

    with pytest.raises(artifact_inspector.ArtifactError, match="required|missing|HDF5|content"):
        artifact_inspector.validate_sdist(archive)


@pytest.mark.parametrize(
    "forbidden",
    [
        "python/package/chiq/native.so",
        "python/package/chiq/native.pyd",
        "python/package/chiq/native.dylib",
        "python/package/chiq/cache.pyc",
        "python/package/chiq/__pycache__/cache",
        "build/output.txt",
        "build-macos/output.txt",
        "nested/build/output.txt",
        "lib/bse-python/chiq/__init__.py",
        "prefix/lib/bse-python/chiq/__init__.py",
        "chiqvars.sh",
        "install_manifest.txt",
    ],
)
def test_rejects_forbidden_built_or_installed_artifacts(tmp_path, forbidden):
    archive = make_sdist(tmp_path, entries=[tar_info(ROOT + "/" + forbidden, b"bad")])

    with pytest.raises(artifact_inspector.ArtifactError, match="forbidden|native|cache|build|installed"):
        artifact_inspector.validate_sdist(archive)


@pytest.mark.parametrize(
    "forbidden",
    [
        "build",
        "python/package/chiq/__pycache__",
        "prefix/lib/bse-python",
    ],
)
def test_rejects_empty_forbidden_directories(tmp_path, forbidden):
    archive = make_sdist(
        tmp_path,
        entries=[tar_info(ROOT + "/" + forbidden, type_=tarfile.DIRTYPE)],
    )

    with pytest.raises(artifact_inspector.ArtifactError, match="forbidden|cache|build|installed"):
        artifact_inspector.validate_sdist(archive)


@pytest.mark.parametrize(
    "filename,root,files",
    [
        ("other-1.2.3.tar.gz", ROOT, required_files()),
        (ROOT + ".tar.gz", "chiq-9.9", required_files()),
        (ROOT + ".tar.gz", ROOT, replace_file(required_files(), "PKG-INFO", b"Metadata-Version: 2.3\nName: other\nVersion: 1.2.3\n\n")),
        (ROOT + ".tar.gz", ROOT, replace_file(required_files(), "PKG-INFO", b"Metadata-Version: 2.3\nName: chiq\nVersion: 9.9\n\n")),
        (ROOT + ".tar.gz", ROOT, replace_file(required_files(), "python/package/chiq/__init__.py", b'__version__ = "9.9"\n')),
        (ROOT + ".tar.gz", ROOT, replace_file(required_files(), "pyproject.toml", required_files()["pyproject.toml"].replace(b'name = "chiq"', b'name = "other"'))),
        (ROOT + ".tar.gz", ROOT, replace_file(required_files(), "pyproject.toml", required_files()["pyproject.toml"].replace(b'dynamic = ["version"]', b'dynamic = ["version"]\nversion = "1.2.3"'))),
        (ROOT + ".tar.gz", ROOT, replace_file(required_files(), "pyproject.toml", required_files()["pyproject.toml"].replace(b"scikit_build_core.metadata.regex", b"other.provider"))),
        (ROOT + ".tar.gz", ROOT, replace_file(required_files(), "pyproject.toml", required_files()["pyproject.toml"].replace(b"python/package/chiq/__init__.py", b"VERSION"))),
    ],
)
def test_rejects_metadata_disagreement(tmp_path, filename, root, files):
    archive = make_sdist(tmp_path, filename=filename, root=root, files=files)

    with pytest.raises(artifact_inspector.ArtifactError, match="name|version|metadata|provider|filename|root"):
        artifact_inspector.validate_sdist(archive)


@pytest.mark.parametrize(
    "pkg_info",
    [
        b"Metadata-Version: 2.3\nName: chiq\nName: chiq\nVersion: 1.2.3\n\n",
        b"Metadata-Version: 2.3\nName: chiq\nVersion: 1.2.3\nVersion: 1.2.3\n\n",
        b"not valid metadata",
    ],
)
def test_rejects_malformed_or_duplicate_pkg_info(tmp_path, pkg_info):
    archive = make_sdist(
        tmp_path,
        files=replace_file(required_files(), "PKG-INFO", pkg_info),
    )

    with pytest.raises(artifact_inspector.ArtifactError, match="PKG-INFO|metadata"):
        artifact_inspector.validate_sdist(archive)


def test_requires_one_shared_top_level_component(tmp_path):
    archive = make_sdist(tmp_path, entries=[tar_info("other-root/file", b"bad")])

    with pytest.raises(artifact_inspector.ArtifactError, match="top-level|root"):
        artifact_inspector.validate_sdist(archive)


def test_cli_writes_manifest_outside_extraction_root(tmp_path):
    archive = make_sdist(tmp_path)
    destination = tmp_path / "extract"
    manifest = tmp_path / "reports" / "sdist.json"
    script = Path(artifact_inspector.__file__)

    subprocess.run(
        [sys.executable, str(script), "sdist", str(archive), "--extract", str(destination), "--manifest", str(manifest)],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    result = json.loads(manifest.read_text())
    assert result["root"] == ROOT
    assert (destination / ROOT / "README.md").is_file()
    assert not (destination / ROOT / "reports" / "sdist.json").exists()
