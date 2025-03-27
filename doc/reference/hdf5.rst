.. _Reference_hdf5:

Structure of HDF5 file
----------------------

Basics
~~~~~~

All input and output data are stored in HDF5 format.
The file names are specified by parameters

.. code-block:: bash

    [bse_common]
    input = "dmft_bse.h5"
    output = "dmft_bse.out.h5"

One may specify the same file for both input and output.

.. note::

    One can print the stored data using ``h5ls`` command as shown below. However, it is recommended to use the class `h5BSE` to save/load data. See :ref:`Reference_h5bse` for details.

.. _reference_hdf5_input:

**dmft_bse.h5**
^^^^^^^^^^^^^^^

One can print the structure of the stored data using ``h5ls`` command.

.. code-block:: bash

    $ h5ls dmft_bse.h5
    bse                      Group

All data are stored in the **bse** group.
One can print subgroups in **bse** group by

.. code-block:: bash

    $ h5ls dmft_bse.h5/bse
    info                     Group
    input                    Group

**info** subgroup contains basic information

.. code-block:: bash

    $ h5ls dmft_bse.h5/bse/info
    beta                     Dataset {SCALAR}
    block_name               Group
    inner_name               Group

**input** subgroup contains input data necessary for BSE calculations

.. code-block:: bash

    $ h5ls dmft_bse.h5/bse/input
    X0_loc                   Group
    X0_q                     Group
    X_loc                    Group
    chi0_loc_in              Group
    chi_loc_in               Group
    gamma0                   Group

.. _reference_hdf5_output:

**dmft_bse.out.h5**
^^^^^^^^^^^^^^^^^^^

The output file **dmft_bse.out.h5** contain **output** subgroup as follows.

.. code-block:: bash

    $ h5ls dmft_bse.out.h5/bse
    info                     Group
    output                   Group

    $ h5ls dmft_bse.out.h5/bse/output
    I_q                      Group
    chi0_loc                 Group
    chi0_q                   Group
    chi_loc                  Group
    chi_q                    Group
    chi_q_rpa                Group
    chi_q_rrpa               Group



Structure of data
~~~~~~~~~~~~~~~~~

The structure of the data is classified into the following three types.

**dict A**

.. code-block:: python

    {(block1, block2): ndarray.complex[N_in1, N_in2, N_w, N_w] }

**dict B**

.. code-block:: python

    {(block1, block2): ndarray.complex[N_in1, N_in2, N_w] }

**dict C**

.. code-block:: python

    {(block1, block2): ndarray.complex[N_in1, N_in2] }

- ``block1``, ``block2``: block index
- ``N_in1``: number of inner_indices in ``block1``
- ``N_in2``: number of inner_indices in ``block2``
- ``N_w``: number of fermionic Matsubara frequencies

dict A is used for two-particle Green's function such as :math:`[X_\text{loc}(i\omega_n, i\omega_{n'})]_{12,34}`, where :math:`1 \equiv (a_1, \sigma_1, m_1)`.
The bosonic frequency has been omitted.
This quantity is reformatted as
:math:`X_\text{loc}(a_1 \sigma_1 \sigma_2, m_1 m_2, n; a_3 \sigma_3 \sigma_4, m_3 m_4, n')`.
Here, we used :math:`a_1=a_2` and :math:`a_3=a_4`.
The indices before and after the semicolon represent row and column, respectively.
Then, this is expressed as a block matrix as ``X_loc[(block1, block2)][i1, i2, w1, w2]``, where :math:`\mathtt{block1} = (a_1 \sigma_1 \sigma_2)`,
:math:`\mathtt{block2} = (a_3 \sigma_3 \sigma_4)`,
:math:`\mathtt{i1} = (m_1 m_2)`,
:math:`\mathtt{i2} = (m_3 m_4)`,
:math:`\mathtt{w1}=n`.
and :math:`\mathtt{w2}=n'`.

dict B is used for bare two-particle Green's function :math:`[X_\text{0,loc}(i\omega_n)]_{12,34}`.
The corresponding block matrix is given by ``X0_loc[(block1, block2)][i1, i2, w]``, where :math:`\mathtt{w}=n`.

dict C represents the structure of the susceptibility such as :math:`[\chi(\boldsymbol{q})]_{12,34}`. There is no fermionic Matsubara frequency.
Therefore, the block matrix is given by ``chi_q[(block1, block2)][i1, i2]``.

List of Data
~~~~~~~~~~~~

**info** subgroup
^^^^^^^^^^^^^^^^^^^

The table below shows the data in **info** subgroup.

.. csv-table::
    :header: type, dtype

    ``beta``, float
    ``block_name``, list[string]
    ``inner_name``, list[string]

**input** subgroup
^^^^^^^^^^^^^^^^^^^^

The table below shows the data in **input** subgroup.

The columns from 'bse' to 'scl' indicate whether the data is referred to by ``chiq_main.py`` and ``chiq_post.py`` when the corresponding approximation scheme (selected by 'type' parameter) is specified. The last two columns indicate data access by ``calc_Iq.py`` and ``calc_Iq_scl.py``.

.. csv-table::
    :header: type, w, q, dtype, 'bse', 'chi0', 'rpa', 'rrpa', 'scl', ``calc_Iq``, ``calc_Iq_scl``

    ``X_loc``, ✓, --, dict A, ◇, --, --, (△*1), (△*1), (△*1), (△*1)
    ``X0_loc``, ✓, --, dict B, ◇, ◇, --, --, ◇, --, ◇
    ``X0_q``, ✓, ✓, dict B, ◇, ◇, --, --, ◇, --, ◇
    ``chi_loc_in``, ✓, --, dict C, △*1, --, --, △*1, △*1, △*1, △*1
    ``chi0_loc_in``, ✓, --, dict C, △*2, △*2, --, --, △*2, --, --
    ``gamma0``, --, --, dict C, --, --, ◇, --, --, --, --
    ``Phi``, ✓, --, dict B, --, --, --, --, ◇, --, --

- ✓: Indicates required elements to get/save the data by ``h5BSE`` class (see :ref:`Reference_h5bse`). For example, ``X_loc`` can be accessed with ``key=('X_loc', w)``, whereas ``X0_q`` is accessed with ``key=('X0_q', w, q)``.
- ◇: Input (mandatory)
- △\*1: Input (recommended). Generated from ``X_loc`` if data is not provided.
- (△\*1): Used if ``chi_loc_in`` is not provided.
- △\*2: Input (recommended). Generated from ``X0_loc`` if data is not provided.

**'output'** subgroup
^^^^^^^^^^^^^^^^^^^^^

The table below shows the data in **output** subgroup and which mode output them.

.. csv-table::
    :header: type, w, q, dtype, 'bse', 'chi0', 'rpa', 'rrpa', 'scl', ``calc_Iq``, ``calc_Iq_scl``

    ``chi0_loc``, ✓, --, dict C, ◆, ◆, --, ◇*3, ◆, --, --
    ``chi0_q``, ✓, ✓, dict C, ◆, ◆, ◇*3, ◇*3, ◆, --, --
    ``chi_loc``, ✓, --, dict C, ◆, --, --, ◆, ◆, --, ◆
    ``chi_q``, ✓, ✓, dict C, ◆, --, --, --, --, ◇*4, --
    ``chi_q_rpa``, ✓, ✓, dict C, --, --, ◆, --, --, --, --
    ``chi_q_rrpa``, ✓, ✓, dict C, --, --, --, ◆, --, --, --
    ``chi_q_scl``, ✓, ✓, dict C, --, --, --, --, ◆, --, ◆
    ``I_q``, ✓, ✓, dict C, ◆, --, --, --, --, ◆, --
    ``I_r``, ✓, ✓*6, dict C, ◆*5, --, --, --, --, --, --
    ``I_q_scl``, ✓, ✓, dict C, --, --, --, --, ◆, --, ◆
    ``I_r_scl``, ✓, ✓*6, dict C, --, --, --, --, ◆*5, --, --


- ◆: Output
- ◇\*3: Used as input. ``chiq_main.py`` with mode='chi0' should be run in advance.
- ◇\*4: Used as input. ``chiq_main.py`` with mode='bse' should be run in advance.
- ◆\*5: Output by ``chiq_fft``.
- ✓*6: ``q`` is used for the real space coordinate ``r``.
