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


def calc_bse(state, beta, nb, nw, n_in):
    info_A = [n_in * nw] * nb
    info_B = [n_in] * (nb * nw)  # noqa: F841
    info_C = [n_in] * nb  # noqa: F841
    X0_loc = state["X0_loc"]
    X0_q = state["X0_q"]
    X_loc = state["X_loc"]
    chi_loc = state["chi_loc"]

    Pq = L.sub(L.block_inverse(X0_loc, [n_in] * (nb * nw)),
               L.block_inverse(X0_q, [n_in] * (nb * nw)))            # B-type
    Xloc_Pq = L.prod_ab(X_loc, Pq, nb, nw, n_in)                     # A-type
    ident_A = L.identity(info_A)
    inner = L.sub(L.block_inverse(L.sub(ident_A, Xloc_Pq), info_A), ident_A)  # A-type
    Zq = L.matmul(inner, X_loc, info_A)                             # A-type
    chi_q = L.add(chi_loc, L.scale(1.0 / beta, L.sumfreq_a(Zq, nb, nw, n_in)))  # C-type
    Iq = L.sub(L.block_inverse(chi_loc, [n_in] * nb),
               L.block_inverse(chi_q, [n_in] * nb))                  # C-type
    return {"chi_q": chi_q, "I_q": Iq}
