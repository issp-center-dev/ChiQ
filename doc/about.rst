About ChiQ
===========

Overview
--------

This package solves the Bethe-Salpeter equation (BSE) in the dynamical mean-field theory (DMFT) and calculates the momentum-dependent static and dynamical susceptibilities in correlated electron systems.
Necessary input to this package can be generated using **DCore**. Combined with **DCore**, this package offers investigations of two-particle responses in models and materials.

Features:

- Computes the momentum-dependent spin/charge/orbital susceptibility :math:`\chi_{ijkl}(\boldsymbol{q}, i\nu_m)` in Matsubara-frequency domain (see :doc:`Algorithms <algorithms>` for the explicit definition).

- Provides several approximation schemes, including BSE, SCL, RPA, and RRPA (see :doc:`Algorithms <algorithms>` for details).

- Runs as a post-processing step of **DCore**.

  - Allows input from DFT packages via tight-binding bands in the Wannier90 format.

  - Supports the inclusion of spin-orbit coupling (SOC).

License
-------

The **ChiQ** package is published under the
`GNU General Public License, version 3 <http://www.gnu.org/licenses/gpl.html>`_.

This package of ver.1.0 was developed under the support of "Project for advancement of software usability in materials science" by The Institute for Solid State Physics, The University of Tokyo. The copyright of DCore ver.1.0 belongs to The University of Tokyo.

.. _chiq_paper:

Citing ChiQ
-------------------

We have not yet published a paper of ChiQ.
We will post the information here as soon as a paper is published.

Furthermore, please cite DCore and related papers if you use them. Please refer to `DCore documentation <https://issp-center-dev.github.io/DCore/master/index.html>`_ for citing DCore.


Developers
-------------------

- ver. 1.0

  J. Otsuki, K. Yoshimi, H. Shinaoka, Y. Motoyama, T. Aoyama


GitHub repository
-----------------

You can download the source code of **ChiQ** from https://github.com/issp-center-dev/ChiQ.

See the `release page on GitHub <https://github.com/issp-center-dev/ChiQ/releases>`_ for the release history.


Disclaimer
----------

The program is provided as is, i.e. WITHOUT ANY WARRANTY of any kind, as
stated in the license.  In particular, its authors and contributors will take
no responsibility for any possible bugs or any improper use of these programs,
including those resulting in incorrect scientific publications.
