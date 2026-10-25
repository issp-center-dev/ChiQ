.. _program_chiq_main:

chiq_main
=========

.. code-block:: bash

    chiq_main.py [-h] [--version] toml

Description
-----------

``chiq_main.py`` is a main script for solving the Bethe-Salpeter equation. MPI-based parallel computation is supported. The solver backend is selected by the ``backend`` option in the ``[chiq_main]`` section (see :ref:`input`): ``"cpp"`` (default) calls the compiled C++ shared library, ``"numpy"`` uses an equivalent pure-Python implementation, and ``"cupy"`` is reserved for a future GPU backend.

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

    $ mpirun -np 4 chiq_main.py parameters.toml

Print version:

.. code-block:: console

    $ chiq_main.py --version
    ChiQ version 1.0-beta
