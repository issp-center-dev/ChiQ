import base64
import csv
import hashlib
import io
from pathlib import Path
import struct
import stat
import subprocess
import sys
import zipfile
import zlib

import pytest
import toml

import artifact_inspector
from contracts import CANONICAL_COMMANDS, POINT_GROUP_MODULES, SHIM_MODULES


VERSION = "1.2.3"
DIST_INFO = "chiq-%s.dist-info" % VERSION
EOCD = struct.Struct("<4s4H2LH")
CENTRAL = struct.Struct("<4s6H3L5H2L")
LOCAL = struct.Struct("<4s5H3L2H")


def _project_contract():
    root = Path(__file__).resolve().parents[4]
    project = toml.load(str(root / "pyproject.toml"))["project"]
    scripts = project["scripts"]
    entry_points = "[console_scripts]\n" + "".join(
        "%s = %s\n" % item for item in scripts.items()
    )
    dependencies = list(project["dependencies"])
    extras = project["optional-dependencies"]
    requires_dist = dependencies + [
        '%s; extra == "%s"' % (requirement, extra)
        for extra, requirements in extras.items()
        for requirement in requirements
    ]
    metadata = [
        "Metadata-Version: 2.3",
        "Name: chiq",
        "Version: %s" % VERSION,
        "Requires-Python: >=3.8",
    ]
    metadata.extend("Provides-Extra: %s" % extra for extra in extras)
    metadata.extend("Requires-Dist: %s" % requirement for requirement in requires_dist)
    return ("\n".join(metadata) + "\n\n").encode(), entry_points.encode()


def wheel_files():
    metadata, entry_points = _project_contract()
    files = {
        "chiq/__init__.py": b'__version__ = "1.2.3"\n',
        "chiq/cli/__init__.py": b"",
        "chiq/cli/_common.py": b"",
        "chiq/cli/_deprecated.py": b"",
        "chiq/point_group_data/__init__.py": b"",
        "chiq/solver/__init__.py": b"",
        "chiq/_bse_solver.cpython-38-x86_64-linux-gnu.so": b"native",
        "bse/__init__.py": b"",
        "bse/point_group_data/__init__.py": b"",
        "bse_solver/__init__.py": b"",
        DIST_INFO + "/METADATA": metadata,
        DIST_INFO + "/WHEEL": b"Wheel-Version: 1.0\nGenerator: test\nRoot-Is-Purelib: false\nTag: cp38-abi3-manylinux_2_17_x86_64\n",
        DIST_INFO + "/entry_points.txt": entry_points,
    }
    files.update(
        ("chiq/cli/%s.py" % name, b"def main(): pass\n")
        for name in CANONICAL_COMMANDS
    )
    files.update(
        ("chiq/point_group_data/%s.py" % name, b"DATA = 1\n")
        for name in POINT_GROUP_MODULES
    )
    files.update(
        ("chiq/solver/%s.py" % name, b"")
        for name in ("base", "cpp", "numpy", "cupy", "kernels", "layout")
    )
    for module in SHIM_MODULES:
        files["chiq/%s.py" % module] = b""
        files["bse/%s.py" % module] = b""
    files.update(
        ("bse/point_group_data/%s.py" % name, b"")
        for name in POINT_GROUP_MODULES
    )
    return files


def record_bytes(files, changes=None, dist_info=DIST_INFO):
    rows = []
    for name, data in files.items():
        digest = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=").decode()
        rows.append([name, "sha256=" + digest, str(len(data))])
    rows.append([dist_info + "/RECORD", "", ""])
    if changes is not None:
        rows = changes(rows)
    stream = io.StringIO(newline="")
    csv.writer(stream, lineterminator="\n").writerows(rows)
    return stream.getvalue().encode()


def make_wheel(tmp_path, files=None, filename=None, record_changes=None,
               compression=zipfile.ZIP_DEFLATED, directories=(),
               directory_compression=zipfile.ZIP_STORED):
    files = dict(wheel_files() if files is None else files)
    dist_infos = {
        name.split("/", 1)[0] for name in files
        if "/" in name and name.split("/", 1)[0].casefold().endswith(".dist-info")
    }
    record_dist_info = next(iter(dist_infos)) if len(dist_infos) == 1 else DIST_INFO
    record = record_bytes(files, record_changes, record_dist_info)
    files[record_dist_info + "/RECORD"] = record
    path = tmp_path / (filename or "chiq-%s-cp38-abi3-manylinux_2_17_x86_64.whl" % VERSION)
    with zipfile.ZipFile(str(path), "w", compression=compression) as archive:
        for name in directories:
            info = zipfile.ZipInfo(name.rstrip("/") + "/")
            info.external_attr = (0o40775 << 16) | 0x10
            info.compress_type = directory_compression
            archive.writestr(info, b"")
        for name, data in files.items():
            archive.writestr(name, data)
    return path


def _eocd(raw):
    offset = raw.rfind(b"PK\x05\x06")
    return offset, list(EOCD.unpack_from(raw, offset))


def _central_entries(raw):
    _offset, eocd = _eocd(raw)
    offset = eocd[6]
    result = []
    for _index in range(eocd[4]):
        values = list(CENTRAL.unpack_from(raw, offset))
        name_start = offset + CENTRAL.size
        name_end = name_start + values[10]
        result.append((offset, values, bytes(raw[name_start:name_end])))
        offset = name_end + values[11] + values[12]
    return result


def _rewrite(path, mutate):
    raw = bytearray(path.read_bytes())
    mutate(raw)
    path.write_bytes(raw)
    return path


def _raw_deflate(data):
    compressor = zlib.compressobj(level=9, wbits=-15)
    return compressor.compress(data) + compressor.flush()


def _replace_deflate_stream(path, name, compressed, declared):
    raw = bytearray(path.read_bytes())
    entries = _central_entries(raw)
    target = next(item for item in entries if item[2].decode("utf-8") == name)
    _target_central, target_values, _target_name = target
    local_offset = target_values[16]
    local = list(LOCAL.unpack_from(raw, local_offset))
    payload_offset = local_offset + LOCAL.size + local[9] + local[10]
    old_compressed_size = target_values[8]
    delta = len(compressed) - old_compressed_size
    old_eocd_offset, eocd = _eocd(raw)
    old_central_offset = eocd[6]
    changed = bytearray(
        raw[:payload_offset]
        + compressed
        + raw[payload_offset + old_compressed_size:]
    )
    crc = zlib.crc32(declared) & 0xFFFFFFFF
    local[6] = crc
    local[7] = len(compressed)
    local[8] = len(declared)
    LOCAL.pack_into(changed, local_offset, *local)
    for central_offset, values, raw_name in entries:
        values = list(values)
        if raw_name.decode("utf-8") == name:
            values[7] = crc
            values[8] = len(compressed)
            values[9] = len(declared)
        if values[16] > local_offset:
            values[16] += delta
        CENTRAL.pack_into(changed, central_offset + delta, *values)
    eocd[6] = old_central_offset + delta
    EOCD.pack_into(changed, old_eocd_offset + delta, *eocd)
    path.write_bytes(changed)
    return path


def test_valid_wheel_returns_manifest_without_extracting(tmp_path):
    archive = make_wheel(tmp_path)
    manifest = artifact_inspector.validate_wheel(archive)
    assert manifest["kind"] == "wheel"
    assert manifest["name"] == "chiq"
    assert manifest["version"] == VERSION
    assert manifest["members"] == sorted(manifest["members"], key=lambda row: row["path"])
    assert not (tmp_path / "chiq").exists()


def test_accepts_and_crc_reads_deflated_explicit_package_directory(tmp_path):
    archive = make_wheel(
        tmp_path,
        directories=("chiq/",),
        directory_compression=zipfile.ZIP_DEFLATED,
    )
    assert artifact_inspector.validate_wheel(archive)["kind"] == "wheel"


@pytest.mark.parametrize("compression", [zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED])
def test_independent_stream_verifier_accepts_valid_wheels(tmp_path, compression):
    archive = make_wheel(tmp_path, compression=compression)
    assert artifact_inspector.validate_wheel(archive)["kind"] == "wheel"


def test_stream_verification_does_not_use_zipfile_reader(tmp_path, monkeypatch):
    archive = make_wheel(tmp_path)
    monkeypatch.setattr(
        artifact_inspector.zipfile,
        "ZipFile",
        lambda *args, **kwargs: pytest.fail("trusted zipfile reader"),
    )
    assert artifact_inspector.validate_wheel(archive)["kind"] == "wheel"


def test_rejects_deflate_output_beyond_declared_size(tmp_path):
    files = wheel_files()
    files["chiq/payload.dat"] = b"x"
    archive = make_wheel(tmp_path, files)
    attack = _raw_deflate(b"x" * (10 * 1024 * 1024))
    _replace_deflate_stream(archive, "chiq/payload.dat", attack, b"x")
    with pytest.raises(artifact_inspector.ArtifactError, match="output|size|declared|stream"):
        artifact_inspector.validate_wheel(archive)


def test_rejects_deflate_output_smaller_than_declared_size(tmp_path):
    files = wheel_files()
    files["chiq/payload.dat"] = b"xx"
    archive = make_wheel(tmp_path, files)
    _replace_deflate_stream(archive, "chiq/payload.dat", _raw_deflate(b"x"), b"xx")
    with pytest.raises(artifact_inspector.ArtifactError, match="output|size|declared|stream"):
        artifact_inspector.validate_wheel(archive)


def test_raw_deflate_verifier_rejects_incorrect_crc(tmp_path):
    archive = make_wheel(tmp_path)
    def mutate(raw):
        central_offset, values, _name = _central_entries(raw)[0]
        values[7] ^= 1
        CENTRAL.pack_into(raw, central_offset, *values)
        local = list(LOCAL.unpack_from(raw, values[16]))
        local[6] = values[7]
        LOCAL.pack_into(raw, values[16], *local)
    _rewrite(archive, mutate)
    with pytest.raises(artifact_inspector.ArtifactError, match="CRC"):
        artifact_inspector.validate_wheel(archive)


@pytest.mark.parametrize(
    "stream",
    [
        _raw_deflate(b"x")[:-1],
        _raw_deflate(b"x") + _raw_deflate(b"second"),
        _raw_deflate(b"x") + b"trailing",
    ],
    ids=("truncated", "concatenated", "trailing"),
)
def test_rejects_incomplete_or_trailing_deflate_stream(tmp_path, stream):
    files = wheel_files()
    files["chiq/payload.dat"] = b"x"
    archive = make_wheel(tmp_path, files)
    _replace_deflate_stream(archive, "chiq/payload.dat", stream, b"x")
    with pytest.raises(artifact_inspector.ArtifactError, match="deflate|stream|trailing|complete"):
        artifact_inspector.validate_wheel(archive)


def test_physical_size_is_checked_before_zip_parsing(tmp_path, monkeypatch):
    archive = make_wheel(tmp_path)
    monkeypatch.setattr(artifact_inspector, "MAX_ARCHIVE", archive.stat().st_size - 1)
    monkeypatch.setattr(artifact_inspector.zipfile, "ZipFile", lambda *args, **kwargs: pytest.fail("parsed ZIP"))
    with pytest.raises(artifact_inspector.ArtifactError, match="size|512"):
        artifact_inspector.validate_wheel(archive)


def test_validated_wheel_is_published_from_private_bytes_with_digest(tmp_path):
    archive = make_wheel(tmp_path)
    expected = hashlib.sha256(archive.read_bytes()).hexdigest()
    published_directory = tmp_path / "validated-wheel"

    manifest = artifact_inspector.validate_wheel(
        archive, publish_directory=published_directory
    )
    published = Path(manifest["published_path"])
    archive.write_bytes(b"attacker replacement")

    assert manifest["archive_sha256"] == expected
    assert hashlib.sha256(published.read_bytes()).hexdigest() == expected
    assert stat.S_IMODE(published_directory.stat().st_mode) == 0o700
    assert stat.S_IMODE(published.stat().st_mode) == 0o600


@pytest.mark.parametrize("field", [1, 2, 3, 4])
def test_rejects_multidisk_or_inconsistent_eocd_fields(tmp_path, field):
    archive = make_wheel(tmp_path)
    def mutate(raw):
        offset, values = _eocd(raw)
        values[field] += 1
        EOCD.pack_into(raw, offset, *values)
    _rewrite(archive, mutate)
    with pytest.raises(artifact_inspector.ArtifactError, match="disk|EOCD|count"):
        artifact_inspector.validate_wheel(archive)


@pytest.mark.parametrize("field", [5, 6])
def test_rejects_inconsistent_central_directory_size_or_offset(tmp_path, field):
    archive = make_wheel(tmp_path)
    def mutate(raw):
        offset, values = _eocd(raw)
        values[field] += 1
        EOCD.pack_into(raw, offset, *values)
    _rewrite(archive, mutate)
    with pytest.raises(artifact_inspector.ArtifactError, match="central|offset|size|EOCD"):
        artifact_inspector.validate_wheel(archive)


def test_rejects_encrypted_entry_using_real_flag_bits(tmp_path):
    archive = make_wheel(tmp_path)
    def mutate(raw):
        central_offset, values, _name = _central_entries(raw)[0]
        values[3] |= 1
        CENTRAL.pack_into(raw, central_offset, *values)
        local = list(LOCAL.unpack_from(raw, values[16]))
        local[2] |= 1
        LOCAL.pack_into(raw, values[16], *local)
    _rewrite(archive, mutate)
    with pytest.raises(artifact_inspector.ArtifactError, match="encrypt"):
        artifact_inspector.validate_wheel(archive)


def test_rejects_data_descriptors_and_zip64(tmp_path):
    archive = make_wheel(tmp_path)
    def mutate(raw):
        central_offset, values, _name = _central_entries(raw)[0]
        values[3] |= 8
        CENTRAL.pack_into(raw, central_offset, *values)
        local = list(LOCAL.unpack_from(raw, values[16]))
        local[2] |= 8
        LOCAL.pack_into(raw, values[16], *local)
    _rewrite(archive, mutate)
    with pytest.raises(artifact_inspector.ArtifactError, match="descriptor|flag"):
        artifact_inspector.validate_wheel(archive)


def test_rejects_zip64_extra_fields(tmp_path):
    archive = make_wheel(tmp_path)
    def mutate(raw):
        offset, values, _name = _central_entries(raw)[0]
        name_end = offset + CENTRAL.size + values[10]
        raw[name_end:name_end] = b"\x01\x00\x00\x00"
        values[11] += 4
        CENTRAL.pack_into(raw, offset, *values)
        eocd_offset, eocd = _eocd(raw)
        eocd[5] += 4
        EOCD.pack_into(raw, eocd_offset, *eocd)
    _rewrite(archive, mutate)
    with pytest.raises(artifact_inspector.ArtifactError, match="ZIP64"):
        artifact_inspector.validate_wheel(archive)


def test_rejects_duplicate_or_overlapping_local_regions(tmp_path):
    archive = make_wheel(tmp_path)
    def mutate(raw):
        entries = _central_entries(raw)
        second_offset, second, _name = entries[1]
        second[16] = entries[0][1][16]
        CENTRAL.pack_into(raw, second_offset, *second)
    _rewrite(archive, mutate)
    with pytest.raises(artifact_inspector.ArtifactError, match="duplicate|overlap|local"):
        artifact_inspector.validate_wheel(archive)


def test_rejects_out_of_range_local_header(tmp_path):
    archive = make_wheel(tmp_path)
    def mutate(raw):
        offset, values, _name = _central_entries(raw)[0]
        values[16] = len(raw) + 1
        CENTRAL.pack_into(raw, offset, *values)
    _rewrite(archive, mutate)
    with pytest.raises(artifact_inspector.ArtifactError, match="range|offset|local"):
        artifact_inspector.validate_wheel(archive)


@pytest.mark.parametrize("central_field,local_field", [(4, 3), (7, 6), (8, 7), (9, 8)])
def test_rejects_local_central_method_crc_or_size_disagreement(tmp_path, central_field, local_field):
    archive = make_wheel(tmp_path)
    def mutate(raw):
        _offset, values, _name = _central_entries(raw)[0]
        local = list(LOCAL.unpack_from(raw, values[16]))
        local[local_field] ^= 1
        LOCAL.pack_into(raw, values[16], *local)
    _rewrite(archive, mutate)
    with pytest.raises(artifact_inspector.ArtifactError, match="local|central|disagree"):
        artifact_inspector.validate_wheel(archive)


def test_rejects_bad_crc_after_reading_member(tmp_path):
    archive = make_wheel(tmp_path, compression=zipfile.ZIP_STORED)
    def mutate(raw):
        _offset, values, _name = _central_entries(raw)[0]
        local_offset = values[16]
        local = LOCAL.unpack_from(raw, local_offset)
        payload = local_offset + LOCAL.size + local[9] + local[10]
        raw[payload] ^= 1
    _rewrite(archive, mutate)
    with pytest.raises(artifact_inspector.ArtifactError, match="CRC|corrupt|malformed"):
        artifact_inspector.validate_wheel(archive)


def test_rejects_symlink_external_attributes(tmp_path):
    archive = make_wheel(tmp_path)
    def mutate(raw):
        offset, values, _name = _central_entries(raw)[0]
        values[15] = 0o120777 << 16
        CENTRAL.pack_into(raw, offset, *values)
    _rewrite(archive, mutate)
    with pytest.raises(artifact_inspector.ArtifactError, match="symlink|file type|special"):
        artifact_inspector.validate_wheel(archive)


@pytest.mark.parametrize("unsafe", [
    "../escape.py", "/absolute.py", "C:/drive.py", "//server/share.py",
    "chiq\\bad.py", "chiq//bad.py", "chiq/./bad.py", "chiq/../bad.py",
])
def test_rejects_unsafe_member_paths(tmp_path, unsafe):
    files = wheel_files()
    files[unsafe] = b"bad"
    archive = make_wheel(tmp_path, files)
    with pytest.raises(artifact_inspector.ArtifactError, match="path|absolute|canonical"):
        artifact_inspector.validate_wheel(archive)


@pytest.mark.parametrize("left,right", [
    ("chiq/A.py", "chiq/a.py"),
    ("chiq/caf\u00e9.py", "chiq/cafe\u0301.py"),
    ("chiq/node", "chiq/node/child.py"),
])
def test_rejects_alias_and_ancestor_file_conflicts(tmp_path, left, right):
    files = wheel_files()
    files[left] = b"a"
    files[right] = b"b"
    archive = make_wheel(tmp_path, files)
    with pytest.raises(artifact_inspector.ArtifactError, match="alias|collision|ancestor|conflict"):
        artifact_inspector.validate_wheel(archive)


def test_resource_limits_and_compression_ratios_use_real_zip_metadata(tmp_path, monkeypatch):
    files = wheel_files()
    files["chiq/data.txt"] = b"x" * 4096
    archive = make_wheel(tmp_path, files)
    monkeypatch.setattr(artifact_inspector, "MAX_FILE", 2048)
    with pytest.raises(artifact_inspector.ArtifactError, match="file|member|size"):
        artifact_inspector.validate_wheel(archive)

    monkeypatch.setattr(artifact_inspector, "MAX_FILE", 8192)
    monkeypatch.setattr(artifact_inspector, "MAX_RATIO", 2)
    with pytest.raises(artifact_inspector.ArtifactError, match="ratio"):
        artifact_inspector.validate_wheel(archive)


def test_aggregate_compression_ratio_has_an_independent_seam(tmp_path, monkeypatch):
    files = wheel_files()
    files["chiq/data.txt"] = b"x" * 4096
    archive = make_wheel(tmp_path, files)
    monkeypatch.setattr(
        artifact_inspector,
        "_zip_entry_ratio_exceeded",
        lambda size, compressed_size: False,
    )
    monkeypatch.setattr(artifact_inspector, "MAX_RATIO", 2)
    with pytest.raises(artifact_inspector.ArtifactError, match="aggregate.*ratio|ratio.*aggregate"):
        artifact_inspector.validate_wheel(archive)


def test_member_count_and_total_limits_have_lowerable_seams(tmp_path, monkeypatch):
    archive = make_wheel(tmp_path)
    monkeypatch.setattr(artifact_inspector, "MAX_MEMBERS", 2)
    with pytest.raises(artifact_inspector.ArtifactError, match="member"):
        artifact_inspector.validate_wheel(archive)


def test_retained_control_files_have_a_lowerable_memory_bound(tmp_path, monkeypatch):
    archive = make_wheel(tmp_path)
    monkeypatch.setattr(artifact_inspector, "MAX_CONTROL_FILE", 1)
    with pytest.raises(artifact_inspector.ArtifactError, match="control.*memory|control.*limit"):
        artifact_inspector.validate_wheel(archive)
    monkeypatch.setattr(artifact_inspector, "MAX_MEMBERS", 50000)
    monkeypatch.setattr(artifact_inspector, "MAX_TOTAL", 2)
    with pytest.raises(artifact_inspector.ArtifactError, match="total|size"):
        artifact_inspector.validate_wheel(archive)


@pytest.mark.parametrize("filename", [
    "other-1.2.3-cp38-abi3-manylinux_2_17_x86_64.whl",
    "chiq-9.9-cp38-abi3-manylinux_2_17_x86_64.whl",
    "not-a-wheel.whl",
])
def test_filename_identity_must_match_project(tmp_path, filename):
    archive = make_wheel(tmp_path, filename=filename)
    with pytest.raises(artifact_inspector.ArtifactError, match="filename|name|version|wheel"):
        artifact_inspector.validate_wheel(archive)


def test_normalized_filename_and_dist_info_name_are_accepted(tmp_path):
    files = {name.replace(DIST_INFO, "ChiQ-1.2.3.dist-info"): data for name, data in wheel_files().items()}
    archive = make_wheel(tmp_path, files, filename="ChiQ-1.2.3-cp38-abi3-manylinux_2_17_x86_64.whl")
    assert artifact_inspector.validate_wheel(archive)["version"] == VERSION


@pytest.mark.parametrize("change", [
    lambda files: files.pop(DIST_INFO + "/METADATA"),
    lambda files: files.pop(DIST_INFO + "/WHEEL"),
    lambda files: files.pop(DIST_INFO + "/entry_points.txt"),
    lambda files: files.update({"other-1.2.3.dist-info/METADATA": b"x"}),
])
def test_requires_exactly_one_complete_dist_info(tmp_path, change):
    files = wheel_files()
    change(files)
    archive = make_wheel(tmp_path, files)
    with pytest.raises(artifact_inspector.ArtifactError, match="dist-info|METADATA|WHEEL|entry"):
        artifact_inspector.validate_wheel(archive)


def test_explicit_empty_second_dist_info_directory_is_rejected(tmp_path):
    archive = make_wheel(tmp_path, directories=("other-1.0.dist-info/",))
    with pytest.raises(artifact_inspector.ArtifactError, match="one|dist-info"):
        artifact_inspector.validate_wheel(archive)


@pytest.mark.parametrize("header,value", [
    ("Name", "other"), ("Version", "9.9"), ("Requires-Python", ">=3.9"),
])
def test_metadata_singletons_and_identity_are_exact(tmp_path, header, value):
    files = wheel_files()
    metadata = files[DIST_INFO + "/METADATA"].decode()
    metadata = metadata.replace("%s: " % header, "%s: %s\n%s: " % (header, value, header), 1)
    files[DIST_INFO + "/METADATA"] = metadata.encode()
    archive = make_wheel(tmp_path, files)
    with pytest.raises(artifact_inspector.ArtifactError, match="metadata|duplicate|name|version|Python"):
        artifact_inspector.validate_wheel(archive)


@pytest.mark.parametrize("metadata_version", ["2.1", "2.2", "2.3", "2.4", "2.5"])
def test_accepts_supported_core_metadata_versions(tmp_path, metadata_version):
    files = wheel_files()
    metadata_path = DIST_INFO + "/METADATA"
    files[metadata_path] = files[metadata_path].replace(
        b"Metadata-Version: 2.3", ("Metadata-Version: %s" % metadata_version).encode()
    )
    archive = make_wheel(tmp_path, files)

    assert artifact_inspector.validate_wheel(archive)["version"] == VERSION


@pytest.mark.parametrize("metadata_version", ["1.2", "2.0", "2.6", "3.0", "garbage"])
def test_rejects_unsupported_core_metadata_versions(tmp_path, metadata_version):
    files = wheel_files()
    metadata_path = DIST_INFO + "/METADATA"
    files[metadata_path] = files[metadata_path].replace(
        b"Metadata-Version: 2.3", ("Metadata-Version: %s" % metadata_version).encode()
    )
    archive = make_wheel(tmp_path, files)

    with pytest.raises(artifact_inspector.ArtifactError, match="metadata version"):
        artifact_inspector.validate_wheel(archive)


def test_rejects_duplicate_core_metadata_version(tmp_path):
    files = wheel_files()
    metadata_path = DIST_INFO + "/METADATA"
    files[metadata_path] = files[metadata_path].replace(
        b"Metadata-Version: 2.3",
        b"Metadata-Version: 2.2\nMetadata-Version: 2.3",
    )
    archive = make_wheel(tmp_path, files)

    with pytest.raises(artifact_inspector.ArtifactError, match="duplicate.*Metadata-Version"):
        artifact_inspector.validate_wheel(archive)


@pytest.mark.parametrize("mutation", ["missing-core", "extra-core", "missing-extra", "wrong-marker"])
def test_dependency_and_extra_contract_is_exact(tmp_path, mutation):
    files = wheel_files()
    text = files[DIST_INFO + "/METADATA"].decode()
    if mutation == "missing-core":
        text = text.replace("Requires-Dist: scipy\n", "")
    elif mutation == "extra-core":
        text = text.replace("\n\n", "\nRequires-Dist: requests\n\n")
    elif mutation == "missing-extra":
        text = text.replace("Provides-Extra: gpu\n", "")
    else:
        text = text.replace('extra == "plot"', 'extra == "gpu"')
    files[DIST_INFO + "/METADATA"] = text.encode()
    archive = make_wheel(tmp_path, files)
    with pytest.raises(artifact_inspector.ArtifactError, match="depend|extra|Requires-Dist|metadata"):
        artifact_inspector.validate_wheel(archive)


@pytest.mark.parametrize("mutation", ["version", "purelib", "tag", "duplicate"])
def test_wheel_control_metadata_is_bound_to_filename(tmp_path, mutation):
    files = wheel_files()
    text = files[DIST_INFO + "/WHEEL"].decode()
    if mutation == "version":
        text = text.replace("Wheel-Version: 1.0", "Wheel-Version: 2.0")
    elif mutation == "purelib":
        text = text.replace("Root-Is-Purelib: false", "Root-Is-Purelib: true")
    elif mutation == "tag":
        text = text.replace("cp38-abi3-manylinux_2_17_x86_64", "py3-none-any")
    else:
        text = text.replace("Wheel-Version: 1.0\n", "Wheel-Version: 1.0\nWheel-Version: 1.0\n")
    files[DIST_INFO + "/WHEEL"] = text.encode()
    archive = make_wheel(tmp_path, files)
    with pytest.raises(artifact_inspector.ArtifactError, match="WHEEL|wheel|tag|Purelib"):
        artifact_inspector.validate_wheel(archive)


@pytest.mark.parametrize("mutation", ["missing", "extra", "target", "section", "duplicate"])
def test_entry_points_are_exact_and_strict(tmp_path, mutation):
    files = wheel_files()
    text = files[DIST_INFO + "/entry_points.txt"].decode()
    first = next(line for line in text.splitlines() if " = " in line)
    if mutation == "missing":
        text = text.replace(first + "\n", "")
    elif mutation == "extra":
        text += "surprise = chiq.cli.chiq_main:main\n"
    elif mutation == "target":
        text = text.replace(first.split(" = ")[1], "wrong.module:main", 1)
    elif mutation == "section":
        text += "[gui_scripts]\nx = y:z\n"
    else:
        text += first + "\n"
    files[DIST_INFO + "/entry_points.txt"] = text.encode()
    archive = make_wheel(tmp_path, files)
    with pytest.raises(artifact_inspector.ArtifactError, match="entry|console|duplicate|script"):
        artifact_inspector.validate_wheel(archive)


@pytest.mark.parametrize("defaults", ["", "shadow = wrong.module:main\n"])
def test_entry_points_forbid_default_section(tmp_path, defaults):
    files = wheel_files()
    entry_points = files[DIST_INFO + "/entry_points.txt"]
    files[DIST_INFO + "/entry_points.txt"] = (
        "[DEFAULT]\n" + defaults
    ).encode() + entry_points
    archive = make_wheel(tmp_path, files)
    with pytest.raises(artifact_inspector.ArtifactError, match="entry|DEFAULT|section"):
        artifact_inspector.validate_wheel(archive)


@pytest.mark.parametrize("signature", ["RECORD.jws", "RECORD.p7s"])
def test_forbids_record_signatures(tmp_path, signature):
    files = wheel_files()
    files[DIST_INFO + "/" + signature] = b"signature"
    archive = make_wheel(tmp_path, files)
    with pytest.raises(artifact_inspector.ArtifactError, match="signature|RECORD"):
        artifact_inspector.validate_wheel(archive)


def _record_mutation(kind):
    def change(rows):
        if kind == "duplicate": rows.append(list(rows[0]))
        elif kind == "missing": rows.pop(0)
        elif kind == "extra": rows.append(["ghost.py", "sha256=" + "A" * 43, "1"])
        elif kind == "columns": rows[0].append("fourth")
        elif kind == "algorithm": rows[0][1] = rows[0][1].replace("sha256=", "md5=")
        elif kind == "padded": rows[0][1] += "="
        elif kind == "unsafe64": rows[0][1] = "sha256=" + "+" * 43
        elif kind == "digest": rows[0][1] = "sha256=" + "A" * 43
        elif kind == "size": rows[0][2] = str(int(rows[0][2]) + 1)
        elif kind == "record-hash": rows[-1][1:] = ["sha256=" + "A" * 43, "0"]
        return rows
    return change


@pytest.mark.parametrize("kind", [
    "duplicate", "missing", "extra", "columns", "algorithm", "padded",
    "unsafe64", "digest", "size", "record-hash",
])
def test_record_contract_rejects_malformed_rows(tmp_path, kind):
    archive = make_wheel(tmp_path, record_changes=_record_mutation(kind))
    with pytest.raises(artifact_inspector.ArtifactError, match="RECORD|hash|size|row|base64"):
        artifact_inspector.validate_wheel(archive)


def test_record_rejects_noncanonical_base64_pad_bits(tmp_path):
    def change(rows):
        alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_"
        digest = rows[0][1][len("sha256="):]
        index = alphabet.index(digest[-1])
        rows[0][1] = "sha256=" + digest[:-1] + alphabet[(index & 0x30) | 1]
        return rows
    archive = make_wheel(tmp_path, record_changes=change)
    with pytest.raises(artifact_inspector.ArtifactError, match="canonical|base64|hash"):
        artifact_inspector.validate_wheel(archive)


def test_record_must_be_strict_utf8_csv(tmp_path):
    files = wheel_files()
    files[DIST_INFO + "/RECORD"] = b"\xff"
    path = tmp_path / "chiq-1.2.3-cp38-abi3-manylinux_2_17_x86_64.whl"
    with zipfile.ZipFile(str(path), "w") as archive:
        for name, data in files.items(): archive.writestr(name, data)
    with pytest.raises(artifact_inspector.ArtifactError, match="UTF-8|CSV|RECORD"):
        artifact_inspector.validate_wheel(path)


@pytest.mark.parametrize("missing", [
    "chiq/__init__.py", "chiq/cli/__init__.py", "chiq/point_group_data/C1.py",
    "bse/__init__.py", "bse_solver/__init__.py",
])
def test_requires_package_and_data_content(tmp_path, missing):
    files = wheel_files(); files.pop(missing)
    archive = make_wheel(tmp_path, files)
    with pytest.raises(artifact_inspector.ArtifactError, match="content|package|point|cli|missing"):
        artifact_inspector.validate_wheel(archive)


@pytest.mark.parametrize(
    "missing",
    ["chiq/cli/%s.py" % name for name in CANONICAL_COMMANDS]
    + ["chiq/point_group_data/%s.py" % name for name in POINT_GROUP_MODULES]
    + ["chiq/%s.py" % name for name in SHIM_MODULES]
    + ["bse/%s.py" % name for name in SHIM_MODULES]
    + ["bse/point_group_data/%s.py" % name for name in POINT_GROUP_MODULES]
    + ["chiq/cli/_common.py", "chiq/cli/_deprecated.py"]
    + ["chiq/solver/%s.py" % name for name in ("base", "cpp", "numpy", "cupy", "kernels", "layout")],
)
def test_requires_every_contract_module(tmp_path, missing):
    files = wheel_files(); files.pop(missing)
    archive = make_wheel(tmp_path, files)
    with pytest.raises(artifact_inspector.ArtifactError, match="missing|required|content|module"):
        artifact_inspector.validate_wheel(archive)


def test_rejects_nested_dist_info_component(tmp_path):
    files = wheel_files(); files["chiq/evil.dist-info/payload"] = b"bad"
    archive = make_wheel(tmp_path, files)
    with pytest.raises(artifact_inspector.ArtifactError, match="dist-info|nested|directory"):
        artifact_inspector.validate_wheel(archive)


@pytest.mark.parametrize("forbidden", [
    "python/scripts/chiq_main.py", "lib/bse-python/chiq/__init__.py", "chiqvars.sh",
    "tests/data.h5", "chiq/__pycache__/x.pyc", "src/solver.cpp",
    "build/temp.o", "CMakeLists.txt", "chiq/data.h5",
])
def test_rejects_source_fixture_cache_and_legacy_layouts(tmp_path, forbidden):
    files = wheel_files(); files[forbidden] = b"bad"
    archive = make_wheel(tmp_path, files)
    with pytest.raises(artifact_inspector.ArtifactError, match="forbidden|content|layout|cache|fixture"):
        artifact_inspector.validate_wheel(archive)


def test_explicit_forbidden_directory_is_rejected(tmp_path):
    archive = make_wheel(tmp_path, directories=("tests/",))
    with pytest.raises(artifact_inspector.ArtifactError, match="forbidden|content|layout"):
        artifact_inspector.validate_wheel(archive)


@pytest.mark.parametrize("native", [
    "chiq/other.so", "chiq/_bse_solver.so.backup.dylib", "bse_solver/_bse_solver.so",
    "_bse_solver.pyd", "chiq/sub/_bse_solver.so",
])
def test_rejects_native_files_outside_exact_solver_allowlist(tmp_path, native):
    files = wheel_files(); files[native] = b"native"
    archive = make_wheel(tmp_path, files)
    with pytest.raises(artifact_inspector.ArtifactError, match="native|solver|extension"):
        artifact_inspector.validate_wheel(archive)


def _replace_native(files, basename):
    native = next(
        name for name in files
        if name.startswith("chiq/_bse_solver") and name.endswith((".so", ".pyd", ".dylib"))
    )
    files.pop(native)
    files["chiq/" + basename] = b"native"


@pytest.mark.parametrize("basename", [
    "_bse_solverevil.so",
    "_bse_solver_evil.so",
    "_bse_solver-evil.so",
    "_bse_solver..so",
    "_bse_solver.-evil.so",
    "_bse_solver._evil.so",
    "_bse_solver.evil.so",
    "_bse_solver.backup.dylib",
    "_bse_solver.foo.bar.pyd",
    "_bse_solver.so.backup.dylib",
    "_bse_solver.cpython-313-evil.so",
    "_bse_solver.cpython-313-darwin.pyd",
    "_bse_solver.cpython-313-darwin.dylib",
    "_bse_solver.cpython-313-arm64-linux-gnu.so",
    "_bse_solver.cpython-313-i686-linux-gnu.so",
    "_bse_solver.cpython-313-ppc64le-linux-gnu.so",
    "_bse_solver.cp313-win_evil.pyd",
    "_bse_solver.cp313-win_amd64.so",
    "_bse_solver.cp310.abi3.win_amd64.pyd",
])
def test_replacement_native_rejects_invalid_stem_or_tags(tmp_path, basename):
    files = wheel_files()
    _replace_native(files, basename)
    archive = make_wheel(tmp_path, files)
    with pytest.raises(artifact_inspector.ArtifactError, match="native|solver|extension"):
        artifact_inspector.validate_wheel(archive)


@pytest.mark.parametrize("basename", [
    "_bse_solver.so",
    "_bse_solver.cpython-310-x86_64-linux-gnu.so",
    "_bse_solver.cpython-313-darwin.so",
    "_bse_solver.cpython-38-darwin.so",
    "_bse_solver.cpython-313t-darwin.so",
    "_bse_solver.cpython-313d-aarch64-linux-gnu.so",
    "_bse_solver.cpython-313t-x86_64-linux-gnu.so",
    "_bse_solver.cpython-313-powerpc64le-linux-gnu.so",
    "_bse_solver.cpython-313-loongarch64-linux-gnu.so",
    "_bse_solver.cpython-38-arm-linux-gnueabihf.so",
    "_bse_solver.cp38-win32.pyd",
    "_bse_solver.cp310-win_amd64.pyd",
    "_bse_solver.cp313t-win_arm64.pyd",
    "_bse_solver.abi3.so",
    "_bse_solver.abi3.pyd",
    "_bse_solver.abi3.dylib",
])
def test_replacement_native_accepts_portable_extension_names(tmp_path, basename):
    files = wheel_files()
    _replace_native(files, basename)
    archive = make_wheel(tmp_path, files)
    assert artifact_inspector.validate_wheel(archive)["kind"] == "wheel"


def test_requires_exactly_one_native_solver(tmp_path):
    files = wheel_files()
    native = next(name for name in files if name.endswith(".so"))
    files.pop(native)
    archive = make_wheel(tmp_path, files)
    with pytest.raises(artifact_inspector.ArtifactError, match="native|solver"):
        artifact_inspector.validate_wheel(archive)


class _CloseFailure:
    def __init__(self, wrapped):
        self._wrapped = wrapped

    def __getattr__(self, name):
        return getattr(self._wrapped, name)

    def close(self):
        self._wrapped.close()
        raise OSError("injected private wheel close failure")


def _inject_private_close_failure(monkeypatch):
    real_copy = artifact_inspector._copy_private_artifact
    monkeypatch.setattr(
        artifact_inspector,
        "_copy_private_artifact",
        lambda path: _CloseFailure(real_copy(path)),
    )


def test_private_wheel_close_failure_is_wrapped_after_success(tmp_path, monkeypatch):
    archive = make_wheel(tmp_path)
    _inject_private_close_failure(monkeypatch)
    with pytest.raises(artifact_inspector.ArtifactError, match="close|private|resource"):
        artifact_inspector.validate_wheel(archive)


def test_private_wheel_close_failure_preserves_active_artifact_error(tmp_path, monkeypatch):
    archive = make_wheel(
        tmp_path,
        filename="other-1.2.3-cp38-abi3-manylinux_2_17_x86_64.whl",
    )
    _inject_private_close_failure(monkeypatch)
    with pytest.raises(artifact_inspector.ArtifactError, match="filename|name|disagreement"):
        artifact_inspector.validate_wheel(archive)


def test_cli_wheel_manifest_is_external_and_does_not_extract(tmp_path):
    archive = make_wheel(tmp_path)
    manifest = tmp_path / "manifest.json"
    script = Path(artifact_inspector.__file__)
    subprocess.run(
        [sys.executable, str(script), "wheel", str(archive), "--manifest", str(manifest)],
        check=True,
    )
    assert '"kind": "wheel"' in manifest.read_text()
    assert not (tmp_path / "chiq").exists()
