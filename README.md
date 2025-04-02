# ChiQ

## What is ChiQ?

This package solves the Bethe-Salpeter equation (BSE) in the dynamical mean-field theory (DMFT) and calculates the momentum-dependent static and dynamical susceptibilities $\chi(q,\nu_m)$ in correlated electron systems.
Necessary input to this package can be generated using [DCore](https://github.com/issp-center-dev/DCore). Combined with DCore, this package offers investigations of two-particle responses in models and materials.

Features:

- Computes the momentum-dependent spin/charge/orbital susceptibility $\chi_{ijkl}(\boldsymbol{q}, i\nu_m)$ in Matsubara-frequency domain (see [algorightm in doc](doc/algorithms.rst) for the explicit definition).

- Provides several approximation schemes, including BSE, SCL, RPA, and RRPA (see [algorithms in doc](doc/altorithms.rst) for details).

- Runs as a post-processing step of DCore.

  - Allows input from DFT packages via tight-binding bands in the Wannier90 format.

  - Supports the inclusion of spin-orbit coupling (SOC).

## Setup

### Requirements

- CMake (>= 3.4)
- C++ compiler compatible with C++11
- [Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page) (>= 3.1)
- Python
  - more-itertools package
    - `python3 -m pip install more-itertools`

### Download

Clone the repository

``` bash
git clone https://github.com/issp-center-dev/ChiQ
```

You also need to clone pybind11, which is linked from extern/pybind11 directory as a submodule. Type

``` bash
git submodule init
git submodule update
```

in the ChiQ directory. pybind11 is cloned into /extern/pybind11.

### Configure to build

Make a new directory and go into it

``` bash
mkdir ChiQ.build; cd ChiQ.build

``` bash
cmake -DCMAKE_BUILD_TYPE=Release \
      -DTesting=ON \
      -DCMAKE_INSTALL_PREFIX=$HOME/local \
      ../ChiQ
```

#### Troubleshooting

- Q: `cmake` is cmake version 2.
  - A: Your system may have `cmake3` instead of `cmake`. Use `cmake3`.
- Q: Eigen3 is not found although it is installed.
  - A: Tell CMake the path to Eigen3 by `-DEIGEN3_DIR=/path/to/eigen3`.
  - For example, `$HOME/local/include/eigen3` is the path to Eigen3 headers, then add `-DEIGEN3_DIR=$HOME/local` to the above cmake command.

### Build and install

After configuration, type the following to build, test, and install

``` bash
make
make test  # when -DTesting=ON is activated in cmake
make install
```

Python scripts such as `bse_tool.py` are installed in **$HOME/local/bin**.
A python package **bse** and a shared library **bse_solver.cpython-XXX-YYY.so** (XXX is the python version, and YYY is the os info) is installed in **$Home/local/lib/bse-python** (or **lib64/bse-python**).
A configurations file **chiqvars.sh** is installed in **$HOME/local/share**, see the next section.

You can build the documentation as follows.

``` bash
pip3 install sphinx wild_sphinx_theme matplotlib
sphinx-build -b html ../bse/doc html
```

### Environment variables

To set up the environment variables **PATH** and **PYTHONPATH** properly, use the configurations file **chiqvars.sh** as follows.

``` bash
source $HOME/local/share/chiqvars.sh

## Test of python package/scripts

The command `make test` tests only C++ codes.
To test Python codes, **complete the above setup first**.
Then, in the build directory,

``` bash
pytest
```

All python tests (script with name `test_*.py`) will run.
You can instead run some tests selectively as

``` bash
pytest tests/python/non-mpi/
pytest tests/python/non-mpi/bsetool_BSE/test_bsetool_BSE.py
```

## How to use

``` bash
mpiexec -np 4 chiq_main.py chiq.toml
chiq_post.py chiq.toml
```

### Input parameters

Parameters are given by a file with TOML format.
An example is shown below.

``` toml
[chiq_common]
input = "dmft_bse.h5"
output = "dmft_bse.out.h5"
type = ['chi0', 'bse', 'rpa', 'rrpa']  # 'chi0', 'bse', 'scl', 'rpa', 'rrpa'
omega_q = "q_path.dat"

[chiq_main]
work_dir = "work/chiq_main"
# num_wf = 20  # If not specified, the value is determined from X_loc

[chiq_post]
output_dir = ""
mode = ['eigen']  # 'matrix_element', 'eigen', 'linear_combination'

# for mode='eigen'
#vector = true  # output eigenvectors
#order = 'file'  # 'descend' (default), 'overlap', 'file'
#order_file = "/path/to/eigenvec.in"  # for order='file'

# for mode='linear_combination'
#coefs_file = "/path/to/coefs.in"
```

### For detailed usage

See the `manual <https://issp-center-dev.github.io/ChiQ/manual/main/html/index.html>`_ for the detailed usage.

## Release History

See the `release page on GitHub <https://github.com/issp-center-dev/ChiQ/releases>`_ for the release history.
