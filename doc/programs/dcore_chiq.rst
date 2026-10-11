dcore_chiq
==========

.. code-block:: bash

    dcore_chiq [-h] --np NP [--version] [path_input_files ...]

Description
-----------

``dcore_chiq`` is an auxiliary program for **DCore**. It generates necessary data for running **ChiQ**.
``dcore_chiq`` should be executed after completing **DCore** calculations using ``dcore_pre`` and ``dcore``.

Positional Arguments
---------------------

**path_input_files**
    Input filename(s).

**--np NP**
    Number of processors.

Options
-------

**-h, --help**
    Show this help message and exit.

**--version**
    Display the program's version number and exit.

Example
-------

The first usage involves adding the ChiQ input parameters directly to the DCore input file (``square.ini``).

.. code-block:: console

    $ dcore_chiq --np 4 square.ini

The second usage keeps the DCore input file (``dcore.ini``) unchanged while specifying the ChiQ input parameters in a separate file (``bse.in``).

.. code-block:: console

    $ dcore_chiq --np 4 dcore.ini bse.ini


Input parameters
----------------

``dcore_chiq`` requires parameters in [bse] block in addition to other parameters in **DCore**.
See `DCore document <https://issp-center-dev.github.io/DCore/master/reference/input.html>`_ for **DCore** parameters.

[bse] block
^^^^^^^^^^^

.. csv-table::
   :header: "Name", "Type", "Default", "Description"
   :widths: auto

   "num_wb", "Integer", "1", "Number of bosonic frequencies (>0)"
   "num_wf", "Integer", "10", "Number of fermionic frequencies (>0)"
   "h5_output_file", "String", "dmft_bse.h5", "Output HDF5 file for bse data"
   "skip_X0q", "Bool", "False", "Skip X0_q calc"
   "skip_Xloc", "Bool", "False", "Skip X_loc calc (for RPA)"
   "calc_only_chiloc", "Bool", "False", "Calculate only chi_loc but no X_loc (for SCL, rRPA). Do not activate skip_Xloc when using this option."
   "use_temp_file", "Bool", "False", "Whether or not temporary file is used in computing X0_q. This option will reduce the memory footprints."
   "X0q_qpoints_saved", "String", "quadrant", "Specifies for which q points X0q are saved in a HDF file. quadrant or path to a q_path.dat file."
   "X0q_mode", "String", "save", "How X0(q) is provided to chiq_main: 'save' stores X0(q); 'gk' stores G(k) instead (see below)."
   "h5_compression", "String", "gzip", "Compression algorithm for HDF5 output"
   "h5_compression_opts", "Integer", "4", "Compression level of gzip for HDF5 output"

Reducing the size of the HDF5 file (``X0q_mode``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default (``X0q_mode = "save"``), X0(q) is stored for every q point selected by ``X0q_qpoints_saved``.
Its size grows as the fourth power of the number of orbitals and linearly in the number of q points, and it usually dominates the file size.

With ``X0q_mode = "gk"``, ``dcore_chiq`` stores the k-resolved Green's function G(k) instead, which grows only as the square of the number of orbitals.
``chiq_main`` then computes X0(q) for exactly the (omega, q) listed in ``omega_q``; any q point on the k mesh can be requested, and ``X0q_qpoints_saved`` is ignored.
The results agree with ``X0q_mode = "save"`` up to floating-point rounding.

Notes:

- ``omega_q`` is required in ``chiq_main`` for files written with ``X0q_mode = "gk"``.
- ``calc_Iq_scl`` does not support ``X0q_mode = "gk"`` yet.
- Rerunning ``dcore_chiq`` in an existing HDF5 file replaces the previous X0 data, but HDF5 does not return the freed space to the file system. Use a new ``h5_output_file`` or ``h5repack`` to obtain a smaller file.
- After rerunning ``dcore_chiq``, recompute ``chi0``/``bse`` into a new ``chiq_main`` output file before running ``rpa``/``rrpa``, which read ``chi0_q`` from that output file.
