Multiple band Hubbard model
=============================

In this tutorial, we will calculate the susceptibility of the multiple band Hubbard model.
For details, see the section VII of the reference paper, `J. Otsuki, et al., Phys. Rev. B 99, 165134 (2019) <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.99.165134>`_.
Sample files for this tutorial are available at ``examples/two_orbital``.

Model
-----

.. math::

  \mathcal{H}
  = \sum_{\boldsymbol{k}, \alpha, \sigma} \varepsilon_{\boldsymbol{k}} c_{\alpha\sigma}^\dagger(\boldsymbol{k}) c_{\alpha\sigma}(\boldsymbol{k})
  + \frac{1}{2} \sum_{i, \alpha, \sigma} \Delta_\alpha n_{i\alpha\sigma}
  + U \sum_{i,\alpha} n_{i\alpha\uparrow} n_{i\alpha\downarrow} \\
  + \frac{U'}{2} \sum_{i,\alpha\ne\beta,\sigma} n_{i\alpha\sigma}n_{i\beta\bar{\sigma}}
  + \frac{U'-J}{2} \sum_{i,\alpha\ne\beta,\sigma} n_{i\alpha\sigma}n_{i\beta\sigma} \\
  + J \sum_{i,\alpha \ne \beta} \left[
  c_{i\alpha\downarrow}^{\dagger}c_{i\beta\uparrow}^{\dagger} c_{i\beta\uparrow}c_{i\alpha\downarrow}
  +
  c_{i\alpha\uparrow}^{\dagger}c_{i\alpha\downarrow}^{\dagger} c_{i\beta\downarrow}c_{i\beta\uparrow}
  \right],

where :math:`i` is the site index and :math:`\alpha`, :math:`\beta` are the orbital indices.
:math:`\sigma` is the spin index and :math:`\bar{\sigma}` means the opposite spin of :math:`\sigma`.
:math:`\varepsilon_{\boldsymbol{k}} = -2t(\cos k_x + \cos k_y)` is the dispersion relation, :math:`\Delta_\alpha` is the crystal field splitting, :math:`U` is the intra-orbital Coulomb interaction, :math:`U'` is the inter-orbital Coulomb interaction, and :math:`J` is the Hund's coupling.
In this tutorial, we consider the two-band Hubbard model with the following parameters:
:math:`t = 1` (set as the unit of energy), :math:`U = 12.0`, :math:`U' = 6.0`, :math:`J = 3.0`, :math:`\Delta_0 = 8.8`, and :math:`\Delta_1 = -8.8`.

DMFT calculation with DCore
----------------------------

To implement this setup in DCore, we need to prepare the input file as follows:

.. literalinclude:: ../../examples/two_orbital/dmft_2orb.in

As for the previous example for the single-band Hubbard model, the ``[model]`` section specifies the model parameters.
``norb`` is set to ``2`` to specify the number of orbitals. ``nelec = 2`` means that the number of electrons per site is 2 (half-filling).
The hopping constants are specified by the Wannier90 format file ``square_2orb_hr.dat``.
When ``norb`` increases, the number of lines in this file will increase as :math:`O(n_{\text{orb}}^2)`.
Hence, a python script ``make_hr_norb.py`` is offered to generate this file as follows:

.. code-block:: bash

  python3 make_hr_norb.py 2 > square_2orb_hr.dat

The ``kanamori`` section specifies `the Kanamori interaction parameters <https://issp-center-dev.github.io/DCore/develop/reference/interaction.html>`_ :math:`(U, U', J) = (12.0, 6.0, 3.0)`.
``local_potential_matrix`` and ``local_potential_factor`` are used to specify `the local_potential parameters <https://issp-center-dev.github.io/DCore/develop/reference/local_potential.html>`_.
In this case, ``local_potential_matrix`` points to the file containing the local potential matrix as follows:

.. literalinclude:: ../../examples/two_orbital/local_pot.in

and ``local_potential_factor`` is set to ``8.8``.

Once the input file is prepared, we can run the DMFT calculation as the previous example:

.. code-block:: bash

  dcore_pre dmft_2orb.in
  dcore --np 4 dmft_2orb.in
  dcore_check dmft_2orb.in

BSE calculation
---------------

First, we need to generate the q-path file from the input file

.. literalinclude:: ../../examples/two_orbital/qpath.in

as follows:

.. code-block:: bash

  gen_qpath.py dmft_2orb.in qpath.in

Next, we need to prepare the input file for the BSE calculation.

.. literalinclude:: ../../examples/two_orbital/bse.in

The ``[bse_common]`` and ``[bse_tool]`` sections are the same as the previous example.
BSE calculation can be performed using the following commands:

.. code-block:: bash

  dcore_chiq.py --np 4 bse.in
  mpiexec --np 4 chiq_main.py bse.in

The obtained susceptibility is a 16x16 complex-valued matrix.
In this example, we transform the susceptibility to the basis described in the reference paper.
The transformation matrix is given by ``eigenvec.in`` as follows:

.. literalinclude:: ../../examples/two_orbital/eigenvec.in

This has the same format as the eigenvector file outputted by ``chiq_post.py``.
To check these basis, use ``eigenvec_viewer.py`` as follows:

.. code-block:: bash

  eigenvec_viewer.py eigenvec.in -a 1 --spin-charge

For example, the 11th eigenvector (``#10``) is as follows:

.. code-block:: bash

  #10
  atom 0
  ch
    0               0
    0               0
  sp+
    0               0
    0               0
  sp-
    0               0
    0               0
  spz
    0               0.707
    0.707           0

``ch``, ``sp+``, ``sp-``, and ``spz`` are the charge :math:`n`, spin-up :math:`\sigma^+`, spin-down :math:`\sigma^-`, and spin-z :math:`\sigma^z`, respectively.
:math:`2 \times 2` matrices are the orbital components :math:`\tau` coupled to each spin-charge component.
Therefore, this says that the 11th eigenvector is :math:`\sigma^z \times \tau^x` (one of the :math:`T_{1u}`).

Basis transformation is performed by ``chiq_post.py`` as follows:

.. code-block:: bash

  chiq_post.py bse.in

The transformed susceptibility is stored in ``bse/chi_q_eigen.dat``.

In the plotting tool, ``plot_chiq_path.py``, you can select the susceptibility to plot by specifying the eigenvector index by using the following file, ``label.in``:

.. literalinclude:: ../../examples/two_orbital/label.in

This file has a python dictionary, where the key is the eigenvector index and the value is the label in the plot.
Note that the line starting with ``#`` is a comment, and then the four bands, 1, 4, 6, and 10, are specified in this case.
``label.in`` is passed to ``plot_chiq_path.py`` as follows:

.. code-block:: bash

  plot_chiq_path.py q_path.dat bse/chi_q_eigen.dat --label=label.in

The obtained figure shows the susceptibility along the path in the Brillouin zone as follows:

.. figure:: multiple_orb/chi_q_eigen_path.*
  :width: 70%
  :align: center

  The susceptibility :math:`\chi(\boldsymbol{q})` calculated by the BSE.
