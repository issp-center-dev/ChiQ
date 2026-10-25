"""CuPy/GPU backend skeleton. Enabled in a later GPU phase."""

from .base import SolverBase


class CupySolver(SolverBase):
    def __init__(self, beta, info_A, info_B, info_C):
        super().__init__(beta, info_A, info_B, info_C)
        raise NotImplementedError(
            "cupy backend not yet enabled (planned for the GPU phase)."
        )

    def set(self, bm, name):
        raise NotImplementedError(
            "cupy backend not yet enabled (planned for the GPU phase)."
        )

    def calc(self, calc_type):
        raise NotImplementedError(
            "cupy backend not yet enabled (planned for the GPU phase)."
        )

    def get(self, name):
        raise NotImplementedError(
            "cupy backend not yet enabled (planned for the GPU phase)."
        )
