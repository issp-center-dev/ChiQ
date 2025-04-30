gen_allq
========

.. code-block:: bash

    gen_allq.py [-h] [-o OUTFILE] file_param

Description
-----------

``gen_allq.py`` is a script for generating a list of all q-points in the Brillouin zone.

Positional Arguments
---------------------

**file_param**
    Config file or HDF5 file. The config file should include:
        - ``(nk0, nk1, nk2)`` or ``nk`` in the ``[model]`` section (**DCore** input),
        - ``(Lx, Ly, Lz)`` in the ``[H0]`` section, or
        - ``dft_h5file`` in the ``[DMFT]`` section.

Options
-------

**-h, --help**
    Show this help message and exit.

**-o OUTFILE, --outfile OUTFILE**
    Specify the output file name. Default is ``"q_fbz.dat"``.

Example
-------

Basic usage:

.. code-block:: console

    $ gen_allq.py dcore.in

Example of the standard output:

.. code-block:: console

    $ gen_allq.py dmft_square.in

    Read k-grid from file 'dmft_square.in'.
    [model]
    (nk0, nk1, nk2) = (32, 32, 1)

    File 'q_fbz.dat' generated.

Example of the output file ``q_fbz.dat``:

.. code-block:: console

    $ cat q_fbz.dat
    0 00.00.00
    0 00.01.00
    0 00.02.00
    0 00.03.00
    0 00.04.00
    0 00.05.00
    0 00.06.00
    0 00.07.00
    0 00.08.00
    0 00.09.00
    0 00.10.00
    0 00.11.00
    0 00.12.00
    0 00.13.00
    0 00.14.00
    0 00.15.00
    0 00.16.00
    0 00.17.00
    0 00.18.00
    0 00.19.00
    0 00.20.00
    0 00.21.00
    0 00.22.00
    0 00.23.00
    0 00.24.00
    0 00.25.00
    0 00.26.00
    0 00.27.00
    0 00.28.00
    0 00.29.00
    0 00.30.00
    0 00.31.00
    0 01.00.00
    0 01.01.00
    0 01.02.00
    ...
