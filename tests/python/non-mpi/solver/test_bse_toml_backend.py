"""Tests for the `[chiq_main] backend` option parsed by bse_toml.

`load_params_from_toml` is a CLI-only helper (called by chiq_main.py /
chiq_post.py); on an invalid `backend` it exits the process via ``sys.exit``
(matching the sibling ``_check_if_dict_empty`` convention), so an invalid value
raises ``SystemExit`` rather than ``ValueError``.
"""

import textwrap

import pytest

from chiq import bse_toml


def _write_toml(tmp_path, backend_line):
    text = textwrap.dedent(
        f"""
        [chiq_common]
        input = "dmft_bse.h5"
        type = ["chi0"]

        [chiq_main]
        {backend_line}

        [chiq_post]
        """
    )
    p = tmp_path / "bse.toml"
    p.write_text(text)
    return str(p)


def test_default_backend_is_cpp(tmp_path):
    toml = _write_toml(tmp_path, "")  # no backend key
    _, dict_tool, _ = bse_toml.load_params_from_toml(toml, print_summary=False)
    assert dict_tool["backend"] == "cpp"


@pytest.mark.parametrize("backend", ["cpp", "numpy", "cupy"])
def test_valid_backends_parse(tmp_path, backend):
    toml = _write_toml(tmp_path, f'backend = "{backend}"')
    _, dict_tool, _ = bse_toml.load_params_from_toml(toml, print_summary=False)
    assert dict_tool["backend"] == backend


def test_invalid_backend_exits(tmp_path):
    toml = _write_toml(tmp_path, 'backend = "gpu999"')
    with pytest.raises(SystemExit) as exc:
        bse_toml.load_params_from_toml(toml, print_summary=False)
    assert "backend must be" in str(exc.value)
