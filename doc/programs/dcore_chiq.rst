dcore_chiq
==========

.. code-block:: bash

    dcore_chiq.py [-h] --np NP [--version] [path_input_files ...]

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

    $ dcore_chiq.py --np 4 square.ini

The second usage keeps the DCore input file (``dcore.ini``) unchanged while specifying the ChiQ input parameters in a separate file (``bse.in``).

.. code-block:: console

    $ dcore_chiq.py --np 4 dcore.ini bse.ini


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
   "skip_X0q", "Bool", "False", "Skip X_0(q) calc"
   "skip_Xloc", "Bool", "False", "Skip X_loc calc (for RPA)"
   "calc_only_chiloc", "Bool", "False", "Calculate only chi_loc but no X_loc (for SCL, rRPA). Do not activate skip_Xloc when using this option."
   "use_temp_file", "Bool", "False", "Whether or not temporary file is used in computing X0_q. This option will reduce the memory footprints."
   "X0q_qpoints_saved", "String", "quadrant", "Specifies for which q points X0q are saved in a HDF file. quadrant or path to a q_path.dat file."
   "h5_compression", "String", "gzip", "Compression algorithm for HDF5 output"
   "h5_compression_opts", "Integer", "4", "Compression level of gzip for HDF5 output"