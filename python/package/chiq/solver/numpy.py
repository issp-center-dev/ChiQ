"""NumPy backend: runs the kernels on complex128 block matrices."""

import numpy as np

from . import kernels
from .base import SolverBase, SET_LAYOUT, CALC_TYPES, GET_NAMES
from .layout import parse_matrix_info

# internal require-names used in the C++ error strings (bse.hpp)
_REQUIRE_NAME = {"X0_loc": "X0_Loc"}  # others equal their set name

# calc_type -> list of required set names (from bse.hpp require() calls)
_REQUIRED = {
    "chi0": ["X0_q", "X0_loc", "chi0_loc"],
    "bse": ["X0_q", "X0_loc", "X_loc", "chi_loc"],
    "rpa": ["gamma0", "chi0_q"],
    "rrpa": ["chi0_loc", "chi_loc", "chi0_q"],
    "scl": ["X0_q", "X0_loc", "Phi", "Phi_sum"],
}


class NumpySolver(SolverBase):
    def __init__(self, beta, info_A, info_B, info_C):
        super().__init__(beta, info_A, info_B, info_C)
        self.nb, self.nw, self.n_in = parse_matrix_info(info_A, info_B, info_C)
        self._in = {}
        self._out = {}

    def set(self, bm, name):
        if name not in SET_LAYOUT:
            raise RuntimeError(f"Invalid type '{name}'")
        self._in[name] = {k: np.asarray(v, dtype=complex) for k, v in bm.items()}

    def _require(self, calc_type):
        for name in _REQUIRED[calc_type]:
            if name not in self._in:
                internal = _REQUIRE_NAME.get(name, name)
                raise RuntimeError(f"'{internal}' must be set before calling calc.")

    def calc(self, calc_type):
        ct = str(calc_type).lower()
        if ct not in CALC_TYPES:
            raise RuntimeError(f"Invalid type '{calc_type}'")
        self._require(ct)
        fn = {
            "chi0": kernels.calc_chi0,
            # bse/rpa/rrpa/scl added in later tasks
        }[ct]
        result = fn(self._in, self.beta, self.nb, self.nw, self.n_in)
        self._out.update(result)
        # I_q and I_q_scl are one slot in C++ (both GetIq()): whichever
        # irreducible vertex was computed most recently answers to both names.
        if "I_q" in result:
            self._out["I_q_scl"] = result["I_q"]
        if "I_q_scl" in result:
            self._out["I_q"] = result["I_q_scl"]

    def get(self, name):
        if name not in GET_NAMES:
            raise RuntimeError(f"Invalid type '{name}'")
        return self._out.get(name, {})
