"""chiq solver subpackage: backend factory."""

from .base import SolverBase


def get_solver(backend, beta, info_A, info_B, info_C):
    """Return a solver facade for the named backend.

    backend: "cpp" (default in the driver), "numpy", or "cupy".
    """
    name = str(backend).lower()
    if name == "cpp":
        from .cpp import CppSolver
        return CppSolver(beta, info_A, info_B, info_C)
    if name == "numpy":
        from .numpy import NumpySolver
        return NumpySolver(beta, info_A, info_B, info_C)
    if name == "cupy":
        from .cupy import CupySolver
        return CupySolver(beta, info_A, info_B, info_C)
    raise ValueError(
        f"Unknown backend {backend!r}; expected one of 'cpp', 'numpy', 'cupy'."
    )
