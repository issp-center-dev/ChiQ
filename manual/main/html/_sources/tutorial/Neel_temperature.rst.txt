Search for the Neel temperature
===============================

In this tutorial, we search for the Neel temperature of the single band Hubbard model on the square lattice.
To this end, we calculate the susceptibility at different temperatures and find the temperature where the inverse susceptibility becomes zero.
Sample files for this tutorial are available at ``examples/square_bse_TNeel``.

Model
-----

The model parameters such as :math:`U` are the same as in the previous tutorial.
The common part of the input file of DCore is ``dmft_square.in``:

.. literalinclude:: ../../examples/square_bse_TNeel/dmft_square.in
  :language: ini

Compared to the previous example, we move the ``[system]`` section to the bottom of the file and remove the ``T`` parameter from it. This is because we will vary the temperature by adding different ``T`` parameters when running the script.

Workflow and results
--------------------

A script file, ``run.sh``, which is shown below:

.. literalinclude:: ../../examples/square_bse_TNeel/run.sh
  :language: bash

To run the script, you need to set the ``NPROCS`` (number of MPI processes). For example, ``NPROCS=1 sh run.sh``.

1. Perform preprocess

  * ``square.h5`` and ``q_path.dat`` are common files for all temperatures.
  * ``chi.dat`` is an empty file to store the susceptibility at Q_M for each temperature.

2. Get the momentum of M point as ``Q_M``.

  * In ``q_path.dat``, the second field of the line with ``M`` indicates the momentum of M point:

    .. code-block:: text

      0 16.16.00   1.20711 M

3. Create a directory for each temperature, copy necessary files and add T to the input file.

  * Since we moved the ``[system]`` section to the bottom, we can simply append ``T = $T`` to the file.

4. Calculate the susceptibility by ChiQ.

  * The procedure is the same as in the previous tutorial.
  * ``bse.in`` is a little bit different from the previous example. We will describe it in the next section.

5. Get the susceptibility at Q_M and save to ``chi.dat``.

  * Find the line with the second field equal to ``Q_M`` in ``bse/chi_q_eigen.dat`` and get the third field.

After running the script, we obtain ``chi.dat`` storing pairs of :math:`T` and :math:`\chi(q_{\boldsymbol{M}})`.

By fitting the inverse susceptibility :math:`1/\chi(q_{\boldsymbol{M}})` by a function :math:`1/\chi = a(T-T_\text{N})`, we estimate the Neel temperature :math:`T_\text{N} = 0.467(1)`.

.. figure:: Neel_temperature/chiinv.*
  :width: 70%
  :align: center

Sorting the eigenvalues
-----------------------

In this example, :math:`\chi(\boldsymbol{q})` has four modes; three of them are degenerate over the whole momentum space and the other is far away from them.
The degenerated modes diverge at M point at the Neel temperature as follows.
As temperature approaches the Neel temperature from above, the susceptibility :math:`\chi(q_{\boldsymbol{M}})` grows and finally diverges.
Below the Neel temperature, :math:`\chi(q_{\boldsymbol{M}})` obtained by solving BSE becomes negative.
In other words, the inverse susceptibility :math:`1/\chi(q_{\boldsymbol{M}})` becomes zero from above at the Neel temperature, and becomes negative below the Neel temperature.
The following figures show :math:`1/\chi(\boldsymbol{q})` for temperatures :math:`T = 0.5 > T_\text{Neel}` and :math:`T = 0.4 < T_\text{Neel}`.

.. figure:: Neel_temperature/chi_q_eigen_path_inv_0.5.*
  :align: center

  :math:`1/\chi(\boldsymbol{q})` for :math:`T = 0.5 > T_\text{Neel}`.

.. figure:: Neel_temperature/chi_q_eigen_path_inv_0.4.*
  :align: center

  :math:`1/\chi(\boldsymbol{q})` for :math:`T = 0.4 < T_\text{Neel}`.

We need to sort the eigenvalues properly to pursue the diverging modes.
By setting ``order = overlap`` in the ``[chiq_post]`` section, the eigenvalues are sorted based on the overlap with the eigenvector at the previous momentum (at the first momentum, X in this case, the eigenvalues are sorted in descending order).
With this change, the ``bse.in`` file is as follows:

.. literalinclude:: ../../examples/square_bse_TNeel/bse.in
  :language: toml
