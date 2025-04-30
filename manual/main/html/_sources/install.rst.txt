
.. highlight:: bash

.. _installation:

Installation
============

Prerequisites
-------------

Before beginning the installation, ensure that you have the following software installed on your system:

- **Git**: Used for cloning the repositories.
- **Python**: Required for using pybind11, which interfaces C++ with Python. Ensure Python is compatible with the required version of pybind11.
- **CMake** (version 3 or higher): Required for building the software.
- **GNU Make** or equivalent: Used for running make commands.
- **C++ Compiler** (supporting C++11 standard): Necessary for compiling the source code.
- **Eigen3** (version 3.1 or higher): C++ header-only library for linear algebra.

Additionally, make sure that you have write access to the installation directory (`$HOME/local` in the example below) and that Python's bin directory is in your system's PATH.

Installation
------------------

First, clone the repository from GitHub:

.. code-block:: bash

    git clone https://github.com/issp-center-dev/ChiQ.git

You also need to clone pybind11, which is linked from extern/pybind11 directory as a submodule. Type

.. code-block:: bash

    git submodule init
    git submodule update

in the ChiQ directory. pybind11 is cloned into ``extern/pybind11``.

For build, make a new directory and go into it

.. code-block:: bash

    mkdir ChiQ.build; cd ChiQ.build

Then, execute the commands below

.. code-block:: bash

    cmake -DCMAKE_BUILD_TYPE=Release \
          -DTesting=ON \
          -DCMAKE_INSTALL_PREFIX=$HOME/local \
          ../ChiQ

to configure the build.
If ``cmake`` is version 2, use ``cmake3`` instead.
If Eigen3 is not found, tell CMake the path to Eigen3 by ``-DEIGEN3_DIR=/path/to/eigen3``.
For example, ``$HOME/local/include/eigen3`` is the path to Eigen3 headers, then add ``-DEIGEN3_DIR=$HOME/local`` to the above cmake command.

Pybind11 may fail to detect the right python version.
In this case, you can pass ``-DPYBIND11_FINDPYTHON=ON`` to the cmake command to use the new CMake FindPython instead of pybind11's custom search.
Or, you can pass ``-DPYTHON_EXECUTABLE=/path/to/python`` to the cmake command to specify the python executable.
See `pybind11 FAQ <https://pybind11.readthedocs.io/en/stable/faq.html#cmake-doesn-t-detect-the-right-python-version>`_ for more details.

After configuration, type the following to build, test, and install

.. code-block:: bash

    make
    make test  # when -DTesting=ON is activated in cmake
    make install

Python scripts such as ``chiq_main.py`` are installed in ``$HOME/local/bin``.
A python package ``chiq`` and a shared library ``bse_solver.cpython-XXX-YYY.so`` (``XXX`` is the python version, and ``YYY`` is the os info) is installed in ``Home/local/lib/bse-python`` (or ``lib64/bse-python``).
A configurations file ``chiqvars.sh`` is installed in ``$HOME/local/share``, see the next section for details.

Environment variables
---------------------

Before running any tests or Python scripts, you need to set up environment variables.
Once you install the package, the configuration file ``chiqvars.sh`` is installed in ``$HOME/local/share``.
To set up the environment variables, type the following command in your terminal:

.. code-block:: bash

    source $HOME/local/share/chiqvars.sh


Dependency on external python packages
--------------------------------------

Provided Python package and scripts depend on external Python packages. Before running the script and tests, install required packages using ``pip`` command:

.. code-block:: bash

    python3 -m pip install numpy matplotlib scipy more-itertools

Test of python package/scripts
------------------------------

The standard ``make test`` command only verifies the C++ components of the software.
To test Python codes, which includes various scripts and packages, you should follow these additional steps after setting up the environment variables:

First, install ``pytest`` python package via ``pip``.

.. code-block:: bash

    python3 -m pip install pytest

Then, in the **build** directory, use the following command to run all Python tests. This includes any script whose name starts with ``test_*.py``.

.. code-block:: bash

    python3 -m pytest

If you prefer to run specific Python tests rather than all tests, you can specify the path to the test directory or file. For example, to run only the non-MPI tests or a specific test:

.. code-block:: bash

    python3 -m pytest tests/python/non-mpi/
    python3 -m pytest tests/python/non-mpi/bsetool_BSE/test_bsetool_BSE.py

This targeted approach helps in debugging specific components without running the entire test suite.


Build documentation
-------------------

To build the documentation, type the following command in the build directory:

.. code-block:: bash

    python3 -m pip install sphinx wild_sphinx_theme
    python3 -m sphinx -b html ../ChiQ/doc html

The documentation is built in the ``html`` directory.
``index.html`` is the main page of the documentation.

If you want to build the PDF documentation, type the following command:

.. code-block:: bash

    python3 -m sphinx -b latex ../ChiQ/doc latex

and then move into ``latex`` directory and use ``make`` to compile the PDF file.

.. code-block:: bash

    cd latex
    make

The PDF file is ``chiq.pdf``.
