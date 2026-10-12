Frequently asked questions
==========================

HDF5 files
----------

ChiQ exchanges all data through HDF5 files (see :ref:`Reference_hdf5`).
This section collects practical points about handling these files.

Why is ``dmft_bse.h5`` so large, and how can I make it smaller?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The file is usually dominated by ``X0_q``, the bare susceptibility stored for every q point.
Its size grows as the fourth power of the number of orbitals,
and linearly in the number of q points (``X0q_qpoints_saved``),
in the number of fermionic frequencies (``num_wf``) and in the number of bosonic frequencies (``num_wb``).

To reduce it:

- Save only the q points you need: set ``X0q_qpoints_saved`` in ``dcore_chiq`` to a ``q_path.dat`` file instead of ``quadrant``.
- Use the smallest ``num_wf`` and ``num_wb`` that give converged results.
- Keep compression enabled (``h5_compression = "gzip"``, the default). A higher ``h5_compression_opts`` makes the file somewhat smaller at the cost of slower writing.

See :doc:`../programs/dcore_chiq` for these parameters.

I deleted or overwrote data, but the file did not shrink
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

HDF5 does not return the space of deleted datasets to the file system.
When ChiQ saves a dataset that already exists, it deletes the old dataset and writes a new one,
so running ``dcore_chiq`` or ``chiq_main`` again on the same file can make it grow.
To reclaim the space, copy the file with ``h5repack`` (part of the HDF5 tools):

.. code-block:: bash

    $ h5repack dmft_bse.h5 dmft_bse_packed.h5
    $ mv dmft_bse_packed.h5 dmft_bse.h5

Writing the output of ``chiq_main`` to a separate file (``output`` in ``[chiq_common]``)
keeps the input file unchanged, so it can be reused without growing.

A file named ``G_k_iw_temp.h5`` is left in the working directory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``dcore_chiq`` writes the k-resolved Green's function to this temporary file when ``use_temp_file = True``.
It is not needed after ``dcore_chiq`` has finished and can be deleted.

``unable to lock file`` errors on a shared file system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

HDF5 locks a file while it is open.
On some network or parallel file systems (for example NFS or Lustre) locking is not supported,
and opening a file fails with an error such as ``unable to lock file``.
In that case, disable file locking before running ChiQ:

.. code-block:: bash

    $ export HDF5_USE_FILE_LOCKING=FALSE

With locking disabled, make sure that no two programs write to the same file at the same time.
For the same reason, do not run ``chiq_main`` or ``chiq_post`` on a file while ``dcore_chiq`` is still writing it.

A run was interrupted. Can I trust the file?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Not necessarily. If a program stops while writing (for example, killed by the job scheduler at the time limit),
the file may be incomplete or, in rare cases, unreadable.
Rerun the step that was interrupted.
For long calculations, keep a copy of ``dmft_bse.h5`` after ``dcore_chiq`` has finished,
so that ``chiq_main`` can be rerun without recomputing ``X0_q``.

How can I look inside a file?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The HDF5 command line tools show the structure and contents:

.. code-block:: bash

    $ h5ls -r dmft_bse.h5          # list all groups and datasets
    $ h5dump -d /bse/info/beta dmft_bse.h5

To read data from Python, use the ``h5BSE`` class described in :ref:`Reference_h5bse`.
