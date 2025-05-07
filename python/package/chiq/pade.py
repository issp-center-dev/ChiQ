# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2024- SpM-lab
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

from typing import Union

import copy

import numpy as np
from numpy.typing import NDArray

from pprint import pprint

class Pade:
    """
    Pade approximation implementation with numerical stabilization
    F(z) = a0/[1 +a1(z-z1)/[1 +a2(z-z2)/[1 +...

    Parameters
    ----------
    z : NDArray[np.complex128]
        Input points
    u : NDArray[np.complex128]
        Function values at input points
    """

    N: int
    "datasize"
    a: NDArray[np.complex128]
    "coefficient"
    z: NDArray[np.complex128]
    "position"
    _debug: bool = False
    "debug mode flag"
    _eps: float = 1e-15
    "numerical zero threshold"

    def __init__(self, z: NDArray[np.complex128], u: NDArray[np.complex128], debug: bool = False):
        assert z.size == u.size
        self.N = z.size
        self._debug = debug
        
        self.z = copy.copy(z)
        
        # Calculate coefficients
        a = copy.copy(u)
        self.re_shift = np.min(np.real(a)) + 1.0
        a += self.re_shift

        if self._debug:
            print(f"Initial values: {a}")
        
        # Calculate Pade approximation coefficients with numerical stabilization
        for p in range(1, self.N):
            for q in range(p, self.N):
                numerator = a[p - 1] - a[q]
                dz = z[q] - z[p - 1]
                
                # Check for too close points
                if abs(dz) < self._eps:
                    raise ValueError(f"Too close points at p={p}, q={q}: z[q]={z[q]}, z[p-1]={z[p-1]}")
                
                # Numerical stabilization for coefficient calculation
                if abs(a[q]) < self._eps:
                    # For zero or near-zero values, use a regularized value
                    a_reg = self._eps * (1.0 if abs(a[p-1]) < self._eps else np.sign(a[p-1]))
                    a[q] = numerator / (dz * a_reg)
                else:
                    # For normal values, use scaled calculation
                    scale = max(abs(a[p-1]), abs(a[q]))
                    if scale > self._eps:
                        a[q] = (numerator / scale) / (dz * (a[q] / scale))
                    else:
                        a[q] = 0.0
                        
            if self._debug:
                print(f"Step {p}: a[{p}] = {a[p]}")
        
        self.a = a

    def evaluate(
        self, w: Union[np.complex128, NDArray[np.complex128]]
    ) -> NDArray[np.complex128]:
        """
        Evaluate Pade approximation at given points
        
        Parameters
        ----------
        w : Union[np.complex128, NDArray[np.complex128]]
            Points to evaluate at
            
        Returns
        -------
        NDArray[np.complex128]
            Approximated values
        """
        w = np.array(w).flatten()
        n = w.size

        P0 = np.zeros(n, dtype=np.complex128)
        P1 = np.zeros(n, dtype=np.complex128)
        P2 = np.zeros(n, dtype=np.complex128)
        Q0 = np.ones(n, dtype=np.complex128)
        Q1 = np.ones(n, dtype=np.complex128)
        Q2 = np.ones(n, dtype=np.complex128)

        P1[:] = self.a[0]
        
        for p in range(1, self.N):
            coeffs = self.a[p] * (w - self.z[p - 1])
            P2[:] = P1 + coeffs * P0
            Q2[:] = Q1 + coeffs * Q0

            if self._debug:
                print(f"Step {p}:")
                print(f"  coeffs = {coeffs}")
                print(f"  P2 = {P2}")
                print(f"  Q2 = {Q2}")

            # Scale each point independently
            scale = np.maximum(np.abs(Q2), self._eps)
            P0[:] = P1[:] / scale
            P1[:] = P2[:] / scale
            P2[:] = P2[:] / scale
            Q0[:] = Q1[:] / scale
            Q1[:] = Q2[:] / scale
            # Q2[:] = np.ones(n, dtype=np.complex128)
            Q2[:] = Q2[:] / scale

        # Final scaling and regularization
        # mask = np.abs(Q1) < self._eps
        # P1[mask] = self.a[0]  # Revert to initial value for unstable points

        return P1/Q2 - self.re_shift


def test_pade():
    # Test with various input patterns
    test_cases = [
        # Normal case
        (np.array([1.0, 2.0, 3.0], dtype=np.complex128),
         np.array([1.0, 2.0, 3.0], dtype=np.complex128)),
        # Zero at first element
        (np.array([1.0, 2.0, 3.0], dtype=np.complex128),
         np.array([0.0, 2.0, 3.0], dtype=np.complex128)),
        # Zero at second element
        (np.array([1.0, 2.0, 3.0], dtype=np.complex128),
         np.array([1.0, 0.0, 3.0], dtype=np.complex128)),
        # All zeros
        (np.array([1.0, 2.0, 3.0], dtype=np.complex128),
         np.array([0.0, 0.0, 0.0], dtype=np.complex128)),
        # Small values
        (np.array([1.0, 2.0, 3.0], dtype=np.complex128),
         np.array([1e-12, 1e-11, 1e-10], dtype=np.complex128)),
        # Matsubara frequency test
        (np.array([1j, 2j, 3j], dtype=np.complex128),
         np.array([1.0, 0.5, 0.33], dtype=np.complex128))
    ]
    
    Failed_cases = []
    for i, (z, u) in enumerate(test_cases):
        print(f"\nTest case {i}:")
        print(f"z = {z}")
        print(f"u = {u}")
        
        pade = Pade(z, u, debug=True)
        u_pade = pade.evaluate(z)
        
        print(f"u_pade = {u_pade}")
        error = np.max(np.abs(u - u_pade))
        print(f"Max error = {error}")
        
        try:
            # Adjust tolerance based on input magnitude
            rtol = max(1e-10, np.max(np.abs(u)) * 1e-10)
            atol = 1e-15
            assert np.allclose(u, u_pade, rtol=rtol, atol=atol)
            print("Test passed!")
        except AssertionError:
            print("Test failed!")
            Failed_cases.append(i)

    if len(Failed_cases) == 0:
        print("All tests passed!")
    else:
        print("Some tests failed!")
        print(f"Failed cases: {Failed_cases}")


if __name__ == "__main__":
    test_pade()
