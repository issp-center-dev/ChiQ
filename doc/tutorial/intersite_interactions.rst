Intersite interactions
===============================

In this tutorial, we calculate the intersite interactions of the single band Hubbard model on the square lattice.
For details of the intersite interactions, see :ref:`Intersite interactions in Algorithms <Algorithm_Iq>`.
Sample files for this tutorial are available at ``examples/square_Ir``.

Model
-----

The model parameters such as :math:`U` are the same as in the previous tutorial.
The input file of DCore is ``dmft_square.in``:

.. literalinclude:: ../../examples/square_Ir/dmft_square.in
  :language: ini

Compared to the previous example, the ``X0q_qpoints_saved`` parameter in the ``[bse]`` section is set to ``q_fbz.dat`` from ``q_path.dat``.
This is because we need to calculate the intersite interactions for all the points in the BZ (``fbz`` means full Brillouin zone) in order to perform FFT from the momentum space to the real space.

Workflow
--------

Calculation of susceptibility :math:`\chi(q)`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, we perform the DMFT calculation as usual.

.. code-block:: bash

  dcore_pre dmft_square.in
  dcore --np 4 dmft_square.in

Before running ``dcore_chiq.py``, we need to prepare the ``q_fbz.dat`` file by using the ``gen_allq.py`` script.

.. code-block:: bash

  gen_allq.py dmft_square.in

Run ``dcore_chiq.py`` and ``chiq_main.py`` to calculate the susceptibility in the momentum space, :math:`\chi(q)` and :math:`\chi_\text{loc}(q)`.

.. code-block:: bash

  dcore_chiq.py --np 4 dmft_square.in
  mpiexec -np 4 chiq_main.py bse.in

The input file of ChiQ tools, ``bse.in``, is also the same as in the previous tutorial:

.. literalinclude:: ../../examples/square_Ir/bse.in
  :language: toml

The only difference is that the ``q_fbz.dat`` file is used instead of ``q_path.dat`` in ``omega_q`` parameter in the ``[bse_common]`` section.


Calculation of intersite interactions :math:`I(q)` and :math:`I(r)`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once these susceptibilities are calculated, we can run ``calc_Iq.py`` to calculate the intersite interactions in the momentum space, :math:`I(q)`.

.. code-block:: bash

  calc_Iq.py -f dmft_bse.out.h5 --remove 1

``-f`` is used to specify the HDF5 file, which is the output file of ``chiq_main.py`` (in this case, ``dmft_bse.out.h5``).
``--remove`` specifies the number of the modes to be removed to stabilize the calculation.
In this case, we remove one mode with the smallest eigenvalue.
The obtained :math:`I(q)` is saved in the HDF5 file as ``dmft_bse.out.h5/bse/output/I_q``.

Next, we perform the FFT from the momentum space :math:`I(q)` to the real space :math:`I(r)` by using ``bse_fft.py``.

.. code-block:: bash

  bse_fft.py -f dmft_bse.out.h5 --input_dname I_q --output_dname I_r dmft_square.in

``-f`` is used to specify the HDF5 file, which is used in ``calc_Iq.py`` (in this case, ``dmft_bse.out.h5``).
``--input_dname`` is used to specify the name of the dataset in the HDF5 file of :math:`I(q)` (in this case, ``I_q``).
``--output_dname`` is used to specify the name of the dataset in the HDF5 file of :math:`I(r)` (in this case, ``I_r``).
In this case, ``bse_fft.py`` reads the :math:`I(q)` from ``dmft_bse.out.h5/bse/output/I_q`` and writes the :math:`I(r)` to ``dmft_bse.out.h5/bse/output/I_r``.

Finally, as in the case of :math:`\chi(q)`, we need to transform (or diagonalize) :math:`I(q)` to :math:`I(r)` by using ``chiq_post.py``.

.. code-block:: bash

  mpiexec --np 4 chiq_post.py bse.in

This also transforms :math:`I(q)` and :math:`I(r)` in the same way as :math:`\chi(q)`.
The results, ``I_q_eigen.dat`` and ``I_r_eigen.dat``, are saved in the directory specified by ``output_dir`` parameter in the ``[bse_post]`` section of ``bse.in`` (in this case, ``bse`` directory).
The file format is also the same as in the :math:`\chi(q)` case.

Plotting the results
~~~~~~~~~~~~~~~~~~~~

:math:`I(q)` can be plotted in the same way as :math:`\chi(q)` in the previous tutorial.

First, we prepare a q-path file, ``q_path.dat``:

.. code-block:: bash

  gen_qpath.py dmft_square.in qpath.in

Then, we plot :math:`I(q)` as follows in the ``bse`` directory:

.. code-block:: bash

  plot_chiq_path.py --label ../label2.in --mode=Iq ../q_path.in I_q_eigen.dat

The following figure ``I_q_eigen_path.pdf`` is saved.

.. figure:: intersite_interactions/I_q_eigen_path.*
  :align: center

  :math:`I(q)` calculated by the BSE.


On the other hand, :math:`I(r)` can be plotted by using ``plot_Ir.py``.
This script plots :math:`I(r)` as a function of the distance from the origin.
In the ``bse`` directory, run the following command:

.. code-block:: bash

  plot_Ir.py --label ../label2.in I_r_eigen.dat

The following figure ``I_r_eigen_distance.pdf`` is saved.

.. figure:: intersite_interactions/I_r_eigen_distance.*
  :align: center

  :math:`I(r)` calculated by the BSE.
