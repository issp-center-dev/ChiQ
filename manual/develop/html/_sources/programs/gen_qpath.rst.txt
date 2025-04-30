.. _program_gen_qpath:

gen_qpath
=========

.. code-block:: bash

    gen_qpath.py [-h] [-o OUTFILE] file_param file_qpoints

Description
-----------

``gen_qpath.py`` is a script for generating a list of q-points on a specified q-path.

Positional Arguments
---------------------

**file_param**
    Config file or HDF5 file. The config file should include:
        - ``(nk0, nk1, nk2)`` or ``nk`` in the ``[model]`` section (**DCore** input),
        - ``(Lx, Ly, Lz)`` in the ``[H0]`` section, or
        - ``dft_h5file`` in the ``[DMFT]`` section.
    and
        - ``bvec`` in ``[model]`` section (**DCore** input). If not given, unit vectors are used.

**file_qpoints**
    File containing a desired q-path.

    Example

    .. code-block:: none

        # qx qy qz name  (in unit of reciprocal lattice vectors)
        1/2 0 0  X
        0 0 0  Gamma
        1/2 1/2 0  M
        1/2 0 0  X
        1/4 1/4 0  M'

Optional Arguments
------------------

**-h, --help**
    Show the help message and exit.

**-o OUTFILE, --outfile OUTFILE**
    Specify the output file name. Default is ``"q_path.dat"``.

Example
-------

Basic usage:

.. code-block:: console

    $ gen_qpath.py dcore.in qpath.in

Example of the standard output:

.. code-block:: console

    $ gen_qpath.py dmft_square.in qpath.in

    Read k-grid from file 'dmft_square.in'.
    [model]
    (nk0, nk1, nk2) = (32, 32, 1)

    Read bvec (reciprocal vectors) from file 'dmft_square.in'.
    bvec not found. Assume the identity matrix.

    Read q-points from 'qpath.in' and convert to grid points.
    q-path runs through (w.r.t. reciprocal primitive vectors)
    [16  0  0] X
    [0 0 0] Gamma
    [16 16  0] M
    [16  0  0] X
    [8 8 0] M'

    File 'q_path.dat' generated

    x-coordinate for plot
    0.00000 X
    0.50000 Gamma
    1.20711 M
    1.70711 X
    2.06066 M'

Example for a generated file ``q_path.dat``:

.. code-block:: console

    $ cat q_path.dat
    0 16.00.00   0.00000 X
    0 15.00.00   0.03125
    0 14.00.00   0.06250
    0 13.00.00   0.09375
    0 12.00.00   0.12500
    0 11.00.00   0.15625
    0 10.00.00   0.18750
    0 09.00.00   0.21875
    0 08.00.00   0.25000
    0 07.00.00   0.28125
    0 06.00.00   0.31250
    0 05.00.00   0.34375
    0 04.00.00   0.37500
    0 03.00.00   0.40625
    0 02.00.00   0.43750
    0 01.00.00   0.46875
    0 00.00.00   0.50000 Gamma
    0 01.01.00   0.54419
    0 02.02.00   0.58839
    ...
