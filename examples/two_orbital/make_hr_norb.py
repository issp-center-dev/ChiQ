"""
Make Wannier90 format file on a square lattice

Usage:
    $ python3 make_hr_norb.py 2

2 is the number of orbitals
"""

import numpy as np
from itertools import product
import sys

# Check if norb is provided as a command line argument
if len(sys.argv) != 2:
    print("Usage: python3 make_hr_norb.py <norb>")
    sys.exit(1)

# Get norb from command line argument
try:
    norb = int(sys.argv[1])
except ValueError:
    print("Error: <norb> should be an integer")
    sys.exit(1)

print("square lattice")
print(norb)
print(9)
print(" ".join(["1"] * 9))

hop = np.zeros((3, 3, 1, norb, norb))

# nearest neighbor
t = -1.0
hop[-1,  0, 0, :, :] = t * np.identity(norb)
hop[ 1,  0, 0, :, :] = t * np.identity(norb)
hop[ 0, -1, 0, :, :] = t * np.identity(norb)
hop[ 0,  1, 0, :, :] = t * np.identity(norb)

z = 0  # 2D
for x, y in product([-1, 0, 1], repeat=2):
    for i, j in product(range(norb), repeat=2):
        t = hop[x, y, z, i, j]
        print(f"{x:2d} {y:2d} {z}  {i+1} {j+1}  {t:.1f} 0.0")
