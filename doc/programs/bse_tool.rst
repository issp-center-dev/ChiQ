.. _program_bse_tool:

bse_tool
========

.. code-block:: bash

    bse_tool.py [-h] [--version] toml

Description
-----------

``bse_tool.py`` is a main script for solving the Bethe-Salpeter equation. MPI-based parallel computation is supported. A shared library implemented in C++ is called internally.

Positional Arguments
---------------------

**toml**
    Parameter file in TOML format. See :ref:`input` for details.

Options
-------

**-h, --help**
    Show this help message and exit.

**--version**
    Display the program's version number and exit.

Example
-------

Basic usage:

.. code-block:: console

    $ mpirun -np 4 bse_tool.py parameters.toml

Print version:

.. code-block:: console

    $ bse_tool.py --version
    BSE version 0.17.0
