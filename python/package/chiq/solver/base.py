"""Solver facade base class and shared name tables."""

import abc

# set-name -> layout type
SET_LAYOUT = {
    "X_loc": "A", "X0_loc": "B", "X0_q": "B", "Phi": "B",
    "chi_loc": "C", "chi0_loc": "C", "chi0_q": "C", "gamma0": "C", "Phi_sum": "C",
}
CALC_TYPES = {"chi0", "bse", "rpa", "rrpa", "scl"}
# valid get() output names (from bse_solver_pybind.cpp getMatrix)
GET_NAMES = {"chi0_q", "chi_q", "chi_q_rpa", "chi_q_rrpa", "chi_q_scl", "I_q", "I_q_scl"}


class SolverBase(abc.ABC):
    def __init__(self, beta, info_A, info_B, info_C):
        self.beta = float(beta)
        self.info_A = list(int(x) for x in info_A)
        self.info_B = list(int(x) for x in info_B)
        self.info_C = list(int(x) for x in info_C)

    @abc.abstractmethod
    def set(self, bm, name):
        ...

    @abc.abstractmethod
    def calc(self, calc_type):
        ...

    @abc.abstractmethod
    def get(self, name):
        ...
