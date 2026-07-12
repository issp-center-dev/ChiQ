"""The five susceptibility calculations on block matrices (layout ops).

Each returns a dict of named output block matrices. Operation order mirrors
src/bse.hpp verbatim (explicit block_inverse then matmul) so the results
agree with the C++ backend within tolerance.
"""

from . import layout as L


def calc_chi0(state, beta, nb, nw, n_in):
    info_B = [n_in] * (nb * nw)  # noqa: F841 (documents B dims)
    diff = L.sub(state["X0_q"], state["X0_loc"])
    chi0_q = L.add(state["chi0_loc"], L.scale(1.0 / beta, L.sumfreq_b(diff, nb, nw, n_in)))
    return {"chi0_q": chi0_q}
