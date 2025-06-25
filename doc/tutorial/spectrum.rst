Spectrum
~~~~~~~~~

We present how to calculate the spectrum. As the simplest example, we use the single-orbital Hubbard model.
The static susceptibility calculation is presented in :doc:`/tutorial/single_orb`. In the following, we show how to extend it to spectrum calculations.


Sample files
------------

Sample files for this tutorial are available at

  `examples/spectrum/ <https://github.com/issp-center-dev/ChiQ/tree/main/examples/spectrum/>`_

in the `GitHub repository <https://github.com/issp-center-dev/ChiQ/tree/main/>`_.


DMFT calculation with DCore
---------------------------

The input file for the DMFT calculation is provided as the file name **dmft_square.in**.

.. literalinclude:: ../../examples/spectrum/dmft_square.in
  :language: ini


The difference with the previous example for the static susceptibility is the value of ``num_wb`` parameter in the ``[bse]`` section. This parameter specifies the number of the bosonic Matsubara frequencies. For spectrum calculation, we need to set some finite value to perform analytical continuation to real frequencies. ``num_wb`` should be large enough, but the computational cost increases proportional to ``num_wb``.

The calculation process is the same as the static susceptibility calculation.

.. code-block:: bash

  dcore_pre dmft_square.in
  dcore --np 4 dmft_square.in
  gen_qpath.py dmft_square.in qpath.in
  dcore_chiq.py --np 4 dmft_square.in

The generated HDF5 file contatins the two-particle Green's function for finite bosonic frequencies. We can confirm the data using ``h5ls`` command as follows.

.. code-block:: bash

  $ h5ls dmft_bse.h5/bse/input/X_loc
  w0                       Group
  w1                       Group
  w2                       Group
  ...

w0, w1, w2, ... indicate the bosonic Matsubara frequencies.

BSE calculation with **ChiQ**
------------------------------

The input file for the BSE calculation is provided as the file name **bse.in**.

.. literalinclude:: ../../examples/spectrum/bse.in
  :language: ini

Two points should be noted in this input file.
``num_wb`` parameter in ``[chiq_common]`` section specifies the number of bosonic frequencies to be computed in **ChiQ**. The default is `-1`, which means all bosonic frequencies in the HDF5 file are computed. If one wants to decrease the number, one can specify a value that is smaller than ``num_wb`` in the **DCore** input, **dmft_square.in**.

The other important point is ``mode`` parameter in ``[chiq_post]`` section. For the spectrum calculations, ``mode='linear_comnination'`` is recommended.
When the magnetic excitation is of interest, one should specify the definition of the magnetic moment by an external file that is specified by ``coefs_file`` parameter. The default file name is **coefs.in**.

In this example, :math:`S_z` operator can be specified by **coefs.in** as follows:

.. literalinclude:: ../../examples/spectrum/coefs.in
  :language: ini

The longitudinal spin excitations are the obtained. For the transverse spin excitations, one specify :math:`S_+` operator indicated in the comment.


The **ChiQ** calculation is performed by the same command as in the static susceptibility calculation.

.. code-block:: bash

  mpiexec -np 4 chiq_main.py bse.in
  chiq_post.py bse.in

The generated output file **chi_q_linear_combination.dat** contains the results for finite bosonic frequencies.
In this file, the first column (integer) is the index :math:`m` of the bosonic frequency :math:`\Omega_m = 2\pi m T`.

We perform an analytical continuation from the Matsubara frequency to real frequencies.
This can be done by the script **chiq_anacont.py** as follows:

.. code-block:: bash

  mpiexec -np 4 chiq_anacont.py bse.in

The output files are saved in the directory specified by ``output_dir`` in the ``[chiq_anacont]`` section of the input file.
If omitted, the output directory is ``anacont`` in the output directory of ``chiq_post.py`` specified by ``output_dir`` in the ``[chiq_post]`` section of the input file.
In this example, therefore, the output directory is ``bse/anacont``.

The method for the analytical continuation is specified by ``method`` parameter in the ``[chiq_anacont]`` section of the input file.
The default method is the Padé approximant, ``method='pade'``.
Alternatively, one can use the method of the sparse modeling (SpM), ``method='spm'`` (see `J. Otsuki et al., Phys. Rev. E 95, 061302(R) (2017) <https://doi.org/10.1103/PhysRevE.95.061302>`).
The SpM method is more robust against the noise than the Padé approximant, but it is more time-consuming and has more hyperparameters to be tuned.

``wmax``, ``wnum``, and ``wmin`` parameters in the ``[chiq_anacont]`` section specify the maximum, minimum, and number of frequencies to be computed, respectively.
The default values are ``wmax=10.0``, ``wmin=0.0``, and ``wnum=101``.

The other parameters are described in :ref:`reference manual of the input file <reference_bse_in>`.

The results are saved in the output directory.
Each q-point data is saved in a separate file as **chi_q_w_XX.YY.ZZ.dat**, where **XX.YY.ZZ** is the q-point index (the second column of the ``chi_q_linear_combination.dat``).
The first column is the real frequency :math:`\omega`.
The second and subsequent columns are the susceptibilities, :math:`\chi_{0}(\boldsymbol{q}, \omega), \chi_{1}(\boldsymbol{q}, \omega), \cdots`.
The susceptibilities are complex-valued, and the real part and the imaginary part are separated by a space as ``Re[chi_0] Im[chi_0] Re[chi_1] Im[chi_1] ...``.

The results can be plotted by the script **plot_chiq_spectrum.py** as

.. code-block:: bash

  plot_chiq_spectrum.py q_path.dat bse/anacont

The first argument is the path to the q-path file. The second argument is the directory containing the output files of ``chiq_anacont.py``.
By default, the script plots the real part of the susceptibility.
To plot the imaginary part, one can specify ``--part imag`` option.

.. figure:: spectrum/chi_q_w.*
  :align: center

  The spectrum of the two-body susceptibility calculated by the BSE.
