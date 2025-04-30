.. _input:

Input files
===========

.. _reference_bse_in:

``bse.in``
----------

``bse.in`` is a file that contains parameters for ``chiq_main.py`` and ``chiq_post.py`` in the TOML format.
Below is a template of the input file.

.. code-block:: toml

    [chiq_common]
    input = "dmft_bse.h5"
    output = "dmft_bse.out.h5"
    type = ["chi0", "bse", "rpa", "rrpa"]  # "chi0", "bse", "scl", "rpa", "rrpa"
    omega_q = "q_path.dat"

    [chiq_main]
    work_dir = "work/chiq"
    # num_wf = 20  # If not specified, the value is determined from X_loc

    [chiq_post]
    output_dir = ""
    mode = ["eigen"]  # "matrix_element", "eigen", "linear_combination"

    # for mode="eigen"
    #vector = true  # output eigenvectors
    #order = "file"  # "descend" (default), "overlap", "file"
    #order_file = "/path/to/eigenvec.in"  # for order="file"

    # for mode="linear_combination"
    # coefs_file = "/path/to/coefs.in"


The ``[chiq_common]`` section contains general settings that apply to both programs ``chiq_main.py`` and ``chiq_post.py``. In contrast, the ``[chiq_main]`` and ``[chiq_post]`` sections provide configurations that are specific to each program.

[chiq_common] section
~~~~~~~~~~~~~~~~~~~~~~

.. csv-table::
   :widths: 10, 10, 20, 60
   :header: "Name", "Type", "Default", "Description"

   "input", "string", ``"dmft_bse.h5"``, "Filename that contain input data in HDF5 format."
   "output", "string", ``"dmft_bse.out.h5"``, "Filename to store output data in HDF5 format. This can be the same as parameter ``input``."
   "type", "list", "``[""bse""]``", "List of calculation schemes. Specify one or more from ``""chi0"", ""bse"", ""scl"", ""rpa"", ""rrpa""``."
   "omega_q", "string", "Not specified", "Filename storing the omega-q list to be computed. If not specified, omega-q list are generated from X0(q) in the input HDF5 file."


[chiq_main] section
~~~~~~~~~~~~~~~~~~~~


.. csv-table::
   :widths: 10, 10, 20, 60
   :header: "Name", "Type", "Default", "Description"

   "work_dir", "string", ``""``, "Working directory where intermediate files are stored. The empty string (default value) means the current directory."
   "num_wf", "int", "Not specified", "Number of fermionic Matsubara frequencies. If not specified, num_wf is determined from X_loc in the input HDF5 file. num_wf can be smaller than the actual number of Matsubara frequencies in the input file."

[chiq_post] section
~~~~~~~~~~~~~~~~~~~~

.. csv-table::
   :widths: 10, 10, 20, 60
   :header: "Name", "Type", "Default", "Description"

   "output_dir", "string", ``""``, "Output directory for the results. The empty string (default value) means the current directory."
   "mode", "list", ``["eigen"]``, "Types of the output. Specify one or more from ``""matrix_element""``, ``""eigen""``, ``""linear_combination""``."
   "vector", "bool", ``false``, "[mode='eigen'] If set to ture, eigenvectors are output."
   "order", "string", ``"descend"``, "[When ``mode=""eigen""``] Sorting order of eigenvalues: ``""descend""``, ``""overlap""``, ``""file""``. For ``""descend""``, eigenvalues are sorted in descending order for eqch q. For ``""overlap""``, eigenvalues are sorted based on overlap with the eigenvectors for the previous q. For ``""file""``, the order is specified by an external file."
   "order_file", "string", ``"eigenvec.in"``, "[When ``mode=""eigen""`` and ``order=""file""``] Filename that contains sequence of vectors that specify the order of the eigenvalues."
   "coefs_file", "string", ``"coefs.in"``, "[When ``mode=""linear_combination""``] Filename that contains coefficients for the linear combination of eigenvectors."

.. _reference_qpath_in:

``qpath.in``
------------

``qpath.in`` is an input file for ``gen_qpath.py`` for generating a q-path.
The file specifies the symmetry points in the Brillouin zone and the path connecting them.
Below is an example of the input file.

.. code-block:: text

   # qx qy qz name  (in unit of reciprocal lattice vectors)
   1/2 0 0  X
   0 0 0  Gamma
   1/2 1/2 0  M
   1/2 0 0  X
   1/4 1/4 0  M'

Lines starting with ``#`` are comments.
Each line consists of four values: the components of the wave vector and the name of the point.
The wave vector is specified in units of the reciprocal lattice vectors, and users can use fractions like ``1/2`` and ``1/4``.
The name of the point will be used in the plotting.

.. _reference_q_path_dat:

``q_path.dat``
----------------

``q_path.dat`` specifies the :math:`(\boldsymbol{q}, i\Omega_m)` points for the calculation of the two-particle susceptibilities.
This has four columns.
The first column is the index :math:`m` of the Matsubara frequency :math:`\Omega_m`.
The second column is the three integers ``X.Y.Z`` separated by a dot without any spaces, which specify the momentum :math:`\boldsymbol{q} = 2\pi/a \times [X, Y, Z]` where :math:`a` is the lattice constant.
The third column is the distance from the starting point along the path.
The optional fourth column is the name of the point, which is used for labeling the plot.
An example of ``q_path.dat`` is as follows

.. code-block:: text

  0 16.00.00   0.00000 X
  0 15.00.00   0.03125
  0 14.00.00   0.06250
  0 13.00.00   0.09375
  0 12.00.00   0.12500
  0 11.00.00   0.15625
     ... continued ...

``eigenvec.in``
----------------

``eigenvec.in`` is an input file for ``chiq_post.py`` to specify the order of the eigenvectors when ``mode = "eigen"`` and ``order = "file"``.
This files specifies the transform matrix :math:`U` in the equation :eq:`chi_eigen` in :ref:`the Algorithm section <Algorithm_Eigen>`.
The format is the same as that of the eigenvectors files generated by ``chiq_main.py``, so see :ref:`the reference of the output files <output_eigenvec_dat>` for details.

``coeff.in``
------------

``coeff.in`` is an input file for ``chiq_post.py`` to specify the coefficients for the linear combination of the density operators when ``mode = "linear_combination"``.
This files specifies the transform matrix :math:`C` in the equation :eq:`chi_linear_combination` in :ref:`the Algorithm section <Algorithm_Eigen>`.
The format is the same as that of the eigenvectors files generated by ``chiq_main.py``, so see :ref:`the reference of the output files <output_eigenvec_dat>` for details.
