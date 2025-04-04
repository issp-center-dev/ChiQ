Output files
=============

The prefix of the output files
--------------------------------

The prefix of the output files shows the type of the quantity and the used method.

.. csv-table::
   :header-rows: 1
   :widths: 50 50

   type part, description
   ``chi_q``, :math:`\chi(q)`
   ``chi0_q``, :math:`\chi_0(q)`
   ``chi_loc``, :math:`\chi_\text{loc}`
   ``I_q``, :math:`I(q)`
   ``I_r``, :math:`I(r)`

.. csv-table::
   :header-rows: 1
   :widths: 50 50

   method part, description
   (nothing), BSE
   ``_rpa``, RPA
   ``_rrpa``, RRPA
   ``_scl``, SCL

For example, ``chi_q_eigen.dat`` means that the susceptibility :math:`\chi(q)` is calculated by the BSE method,
and ``I_q_rpa_eigen.dat`` means that the intersite interaction :math:`I(q)` is calculated by the RPA method.

This scheme is applied to all the output files descibed below.

``chi_q_eigen.dat``
-------------------

This file is generated by ``chiq_post.py``, which contains `the eigenvalues <Algorithm_Eigen_values>`_ or `the linear combinations <Algorithm_Eigen_linear_combination>`_ of the two-particle susceptibilities :math:`\chi_{\xi}(\boldsymbol{q}, i\Omega_m)` or :math:`\chi_{\gamma\gamma}(\boldsymbol{q}, i\Omega_m)`, respectively.

This file has three or more columns.
The first column is the index :math:`m` of the Matsubara frequency :math:`\Omega_m`.
The second column is the three integers ``X.Y.Z`` separated by a dot without any spaces, which specify the momentum :math:`\boldsymbol{q} = 2\pi/a \times [X, Y, Z]` where :math:`a` is the lattice constant.
The third and subsequent columns are the susceptibilities, :math:`\chi_{0}(\boldsymbol{q}, i\Omega_m), \chi_{1}(\boldsymbol{q}, i\Omega_m), \cdots`.
The following is the first few lines of an example.

.. code-block:: text

  #Temperature: 0.5
  0 00.00.00 0.5851252932088625 0.5851252932088622 0.5851252932088622 -0.010185534590226983
  0 01.00.00 0.5871784814635085 0.587178481463508 0.5871784814635079 -0.010043078925395688
  0 01.01.00 0.5892514011643883 0.5892514011643873 0.5892514011643873 -0.009901741925596026

.. _output_eigenvec_dat:

``chi_q_eigenvector.XX.YY.ZZ.dat``
-------------------------------------

These files are generated by ``chiq_post.py`` with ``vector = true`` in ``[chiq_post]``, which contain `the eigenvectors <Algorithm_Eigen>`_ of the two-particle susceptibilities, :math:`U_{kl}^{(\xi)}(\boldsymbol{q}, i\Omega_m)`.
``XX.YY.ZZ`` is the momentum :math:`\boldsymbol{q} = 2\pi/a \times [X, Y, Z]` where :math:`a` is the lattice constant.
By combining :math:`k` and :math:`l` into one index, :math:`U_{kl}^{(\xi)}` is represented as a :math:`N \times N` complex-valued matrix with :math:`\xi` as row and :math:`(k,l)` as column.

These files have :math:`N` rows and :math:`2N` columns.
One row corresponds to one eigenvector.
In one row, the real part and the imaginary part are alternated; real, imaginary, real, imaginary, and so on.
These files can be analyzed by `eigenvec_viewer.py <program_eigenvec_viewer>`_.

``chi_q_eigen_path.pdf``
------------------------

This file is generated by ``plot_chiq_path.py``, which shows the plot of the susceptibilities :math:`\chi_{\xi}(q)` along the path.

``chi_q_eigen_path_inv.pdf``
----------------------------

This file is generated by ``plot_chiq_path.py``, which shows the plot of the inverse susceptibilities :math:`1/\chi_{\xi}(q)` along the path.

``I_r_eigen_distance.pdf``
----------------------------

This file is generated by ``plot_Ir.py``, which shows the plot of the intersite interaction :math:`I(r)` as a function of the distance from the origin.

``chi_q_eigen_path.dat``
------------------------

This file is generated by ``plot_chiq_path.py`` or ``plot_Ir.py`` by option ``--data_out FILENAME``.
This file stores the numerical data of the susceptibilities along the path, which helps to draw the plot of the susceptibilities along the path manually.
The first column is distance from the origin along the path (``plot_chiq_path.py``) or the distance from the origin (``plot_Ir.py``).
The remaining columns are the values :math:`\chi_{0}(\boldsymbol{q}, i\Omega_m), \chi_{1}(\boldsymbol{q}, i\Omega_m), \cdots`.
