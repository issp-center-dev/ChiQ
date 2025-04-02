.. _program_calc_iq_scl:

calc_Iq_scl
===========

.. code-block:: bash

    calc_Iq_scl.py [-h] [--verbose] input_file

Description
-----------

``calc_Iq_scl.py`` is a script for calculating the q-dependent interactions and susceptibility based on the SCL formula. For the detailed method, see :ref:`Algorithm_scl`.

Positional Arguments
---------------------

**input_file**
    Input config file. The parameters below should be given in  ``[SCL]`` block.

    **Mandatory parameters**:

    - **E_plus** : ``(float)`` Excitation energy from n-particle state to (n+1)-particle states.
    - **E_minus**: ``(float)`` Excitation energy from n-particle state to (n-1)-particle states.

    **Optional parameters**:

    - **input**: (string) HDF5 file containing input data. Default is "dmft_bse.h5".
    - **output**: (string) HDF5 file in which output data will be stored. Default is "dmft_bse.out.h5".
    - **iw_cutoff**: (int) Number of Matsubara frequencies used in the summation. Default is 0 (no cutoff).
    - **wb**: (int) Bosonic Matsubara frequency. Default is 0 (static component).
    - **verbose**: (bool) Verbose output. Default is ``False``.

Options
-------

**-h, --help**
    Show this help message and exit.

**--verbose**
    Enable verbose output for detailed logs during execution.

Example
-------

Basic usage:

.. code-block:: console

    $ calc_Iq_scl.py scl_2pole.in

Example of the input file ``scl_2pole.in``:

.. code-block:: none

    [SCL]
    E_plus = 4.0
    E_minus = 4.0

    # input = dmft_bse.h5
    # output = dmft_bse.out.h5
    # iw_cutoff = 20
    # wb = 0
    # verbose = True

The input file ``dmft_bse.h5`` should contain the following data:

.. code-block:: console

    $ h5ls dmft_bse.h5/bse/input
    X0_loc                   Group
    X0_q                     Group
    chi_loc_in               Group

..
    $ h5ls dmft_bse.h5/bse/input
    X0_loc                   Group
    X0_q                     Group
    chi0_loc_in              Group
    chi_loc_in               Group
    gamma0                   Group

Example of the standard output:

.. code-block:: console

    $ calc_Iq_scl.py scl_2pole.in
    Namespace(input_file="scl_2pole.in", verbose=False)

    Read file 'scl_2pole.in'
    input      = dmft_bse.h5
    output     = dmft_bse.out.h5
    E_plus     = 4.0
    E_minus    = 4.0
    iw_cutoff  = 0
    wb         = 0
    verbose    = False

    Load data from 'dmft_bse.h5'
    beta = 2.0
    inner_name: [b'0-0']
    block_name: [b'0-up-0-up', b'0-up-0-down', b'0-down-0-up', b'0-down-0-down']

    Load 'chi_loc'
    'chi_loc_in' found

    Eigenvalues of chi_loc
    0: 0.9996646498695337
    1: 0.9996646498695336
    2: 0.9996646498695334
    3: 0.0003353501304669271

    Load 'X0_loc'
    n_iw = 400
    iw_max = 626.7477343911637

    55 data of 'X0_q' found
    Compute 'I_q_scl' and 'chi_q_scl' (existing data are overwritten)

    Results are saved to 'dmft_bse.out.h5'

    Start q-loop
    End q-loop

The generated output file ``dmft_bse.out.h5`` contain data as follows:

.. code-block:: console

    $ h5ls dmft_bse.out.h5/bse/output
    I_q_scl                  Group
    chi_loc                  Group
    chi_q_scl                Group
