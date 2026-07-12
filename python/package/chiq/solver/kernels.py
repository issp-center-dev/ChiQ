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


def calc_rpa(state, beta, nb, nw, n_in):
    info_C = [n_in] * nb
    chi0_q = state["chi0_q"]
    gamma0 = state["gamma0"]
    m = L.sub(L.identity(info_C), L.matmul(chi0_q, gamma0, info_C))
    chi_q_rpa = L.matmul(L.block_inverse(m, info_C), chi0_q, info_C)
    return {"chi_q_rpa": chi_q_rpa}


def calc_rrpa(state, beta, nb, nw, n_in):
    info_C = [n_in] * nb
    chi0_loc = state["chi0_loc"]
    chi_loc = state["chi_loc"]
    chi0_q = state["chi0_q"]
    Ueff = L.sub(L.block_inverse(chi0_loc, info_C), L.block_inverse(chi_loc, info_C))
    m = L.sub(L.identity(info_C), L.matmul(chi0_q, Ueff, info_C))
    chi_q_rrpa = L.matmul(L.block_inverse(m, info_C), chi0_q, info_C)
    return {"chi_q_rrpa": chi_q_rrpa}
