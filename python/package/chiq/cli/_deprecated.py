"""Deprecated `<name>.py` console-script aliases.

Each callable prints a stderr deprecation notice (so shell users see it even
though DeprecationWarning is filtered) then delegates to the real entry point.
The `.py` command names are removed in ChiQ 2.0.
"""
import sys
import importlib


def _make(name):
    def _entry():
        sys.stderr.write(
            "warning: the '%s.py' command is deprecated and will be removed in "
            "ChiQ 2.0; use '%s' instead.\n" % (name, name)
        )
        return importlib.import_module("chiq.cli." + name).main()
    _entry.__name__ = name + "_py"
    return _entry


chiq_main_py = _make("chiq_main")
chiq_post_py = _make("chiq_post")
chiq_fft_py = _make("chiq_fft")
gen_qpath_py = _make("gen_qpath")
gen_allq_py = _make("gen_allq")
calc_Iq_py = _make("calc_Iq")
calc_Iq_scl_py = _make("calc_Iq_scl")
plot_chiq_path_py = _make("plot_chiq_path")
plot_Ir_py = _make("plot_Ir")
eigenvec_viewer_py = _make("eigenvec_viewer")
dcore_chiq_py = _make("dcore_chiq")
