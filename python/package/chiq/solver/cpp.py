"""C++ backend: thin wrapper over the pybind module `bse_solver`."""

import numpy as np
import bse_solver  # top-level pybind module (packaging phase renames to chiq._bse_solver)

from .base import SolverBase


class CppSolver(SolverBase):
    def __init__(self, beta, info_A, info_B, info_C):
        super().__init__(beta, info_A, info_B, info_C)
        self._impl = bse_solver.BSESolver(
            self.beta,
            np.array(self.info_A),
            np.array(self.info_B),
            np.array(self.info_C),
        )

    def set(self, bm, name):
        self._impl.set(bm, name)

    def calc(self, calc_type):
        self._impl.calc(calc_type)

    def get(self, name):
        return self._impl.get(name)
