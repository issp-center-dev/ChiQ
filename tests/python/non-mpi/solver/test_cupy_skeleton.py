import pytest
from chiq.solver import get_solver

def test_cupy_backend_not_enabled():
    with pytest.raises((NotImplementedError, ImportError, RuntimeError)) as exc:
        get_solver("cupy", 10.0, [2], [1, 1], [1])
    assert "cupy" in str(exc.value).lower()
