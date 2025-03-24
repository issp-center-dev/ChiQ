.. _calc_Iq:

calc_Iq
=======

.. code-block:: bash

    calc_Iq.py [-h] [-f FILE] [--remove REMOVE | --retain RETAIN] [-w W]
               [--verbose]

Description
-----------

``calc_Iq.py`` is a script for calculating the momentum-dependent interactions from the susceptibility, assuming a localized model. This script should be run after the BSE calculation by ``bse_tool.py``. For the explicit equation, see :ref:`Algorithm_Iq`.

Options
-------

**-h, --help**
    Show this help message and exit.

**-f FILE, --file FILE**
    Specify the HDF5 filename. Default is ``"dmft_bse.out.h5"``. The file should contain ``chi_loc`` and ``chi_q`` data. The results are stored as ``I_q`` (overwritten if data exist).

**--remove REMOVE**
    Specify the number of eigenmodes to remove. Default is ``1``, typically used to exclude the charge susceptibility.

**--retain RETAIN**
    Specify the number of eigenmodes to retain. Cannot be used simultaneously with ``--remove``.

**-w W**
    Specify the bosonic Matsubara frequency. Default is ``0`` (static susceptibility).

**--verbose**
    Enable verbose output for detailed execution logs.

.. note::

    Options ``--remove`` and ``--retain`` are mutually exclusive.

Example
-------

Example for using ``calc_Iq.py`` in the single-orbital Hubbard model:

.. code-block:: console

    $ calc_Iq.py
    Namespace(file='dmft_bse.out.h5', remove=1, retain=None, w=0, verbose=False)

    Load data from 'dmft_bse.out.h5'
    inner_name: [b'0-0']
    block_name: [b'0-up-0-up', b'0-up-0-down', b'0-down-0-up', b'0-down-0-down']

    Load 'chi_loc'
    'chi_loc' found

    Diagonalize chi_loc
    Number of eigenmodes retained = 3
    0: 0.9996646498695337
    1: 0.9996646498695336
    2: 0.9996646498695334
    ------------------------- Remove eigenmodes below
    3: 0.0003353501304669271

    55 data of 'chi_q' found
    Compute 'I_q' (existing data are overwritten)

    Start q-loop
    End q-loop

The description "------------------------- Remove eigenmodes below" indicates that the charge mode with the smallest eigenvalue of ``chi_loc`` is removed by default (``--remove 1``). This treatment ensures a stable inversion of the ``chi_loc`` matrix.

The results are saved as ``I_q`` in the HDF5 file.

.. code-block:: none

    $ h5ls dmft_bse.out.h5/bse/output
    I_q                      Group
    chi0_loc                 Group
    chi0_q                   Group
    chi_loc                  Group
    chi_q                    Group

.. note::

    The file ``dmft_bse.out.h5`` originally contains ``I_q`` computed in ``bse_tool.py``. It corresponds to a result with no remove option ``--remove 0``. However, removing the charge mode is crucial especially in the strong-coupling regime. We recommend using ``calc_Iq.py`` and removing small eigenvalues manually.
