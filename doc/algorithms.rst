Algorithms
==========

Definition of the susceptibility
--------------------------------

The momentum-dependent susceptibility computed by **ChiQ** is defined by

.. math::
   :label: chi_q

   \chi_{ij,kl}(\boldsymbol{q}, i\Omega_m) = \frac{1}{N} \sum_{\boldsymbol{r}} \int_{0}^{\beta} d\tau \langle O_{ij}(\boldsymbol{r}, \tau) O_{lk}(\boldsymbol{0}, 0) \rangle  e^{i\Omega_m \tau} e^{-i\boldsymbol{q} \cdot \boldsymbol{r}},

where :math:`\boldsymbol{q}` represents the momentum, :math:`\Omega_m` represents the bosonic Matsubara frequency, :math:`\beta` is the inverse temperature, and
:math:`\boldsymbol{r}` denotes the lattice positions (the origin of each unit cell) with :math:`N` being their total number.
The superscript :math:`i, j, k, l` indicates the combined spin and orbital indices. Here, the term "orbital" includes the site within a unit cell, as defined in the tight-binding model in the Wannier90 format.
The angle brackets denote the thermal average.
The operator :math:`O_{ij}(\boldsymbol{r}, \tau)` is defined by

.. math::

   O_{ij}(\boldsymbol{r}) = c_{i}^{\dagger}(\boldsymbol{r}) c_{j}(\boldsymbol{r}),

and the argument :math:`\tau` indicates the time evolution in the Heisenberg picture.


Calculation of the susceptibility
---------------------------------

.. _Algorithm_BSE:

Bethe-Salpeter equation (BSE)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This part is based on Section II.A in Ref. :ref:`[Otsuki et al. 2019] <references>`.

In DMFT, we compute the two-particle Green's function :math:`X_{ij,kl}(i\omega_n, i\omega_{n'}; i\Omega_m)`, which has two additional fermionic Matsubara frequencies :math:`i\omega_n` and :math:`i\omega_{n'}`. The susceptibility :math:`\chi_{ij,kl}(\boldsymbol{q}, i\Omega_m)` in Eq. :eq:`chi_q` is obtained by summing over the fermionic frequencies as

.. math::
    :label: X_to_chi

    \chi_{ij,kl}(\boldsymbol{q}, i\Omega_m) = T \sum_{nn'} X_{ij,kl}(i\omega_n, i\omega_{n'}; \boldsymbol{q}, i\Omega_m).

We introduce a matrix notation by combining spin-orbital indices and the frequency index as

.. math::

    [\hat{X}(\boldsymbol{q}, i\Omega_{m})]_{(ijn),(kln')} \equiv X_{ij,kl}(i\omega_{n}, i\omega_{n'}; \boldsymbol{q}, i\Omega_{m}).

Then, the Bethe-Salpeter (BS) equation that :math:`\hat{X}(\boldsymbol{q}, i\Omega_{m})` follows is given by

.. math::

    \hat{X}^{-1}(\boldsymbol{q}, i\Omega_{m}) = \hat{X}_0^{-1}(\boldsymbol{q}, i\Omega_{m}) - \hat{\Gamma}_{\text{loc}}(i\Omega_{m}),

where :math:`\hat{X}_0(\boldsymbol{q}, i\Omega_{m})` is the bare two-particle Green's function

.. math::
    :label: X0q

    [\hat{X}_0(\boldsymbol{q}, i\Omega_{m})]_{(ijn),(kln')} = -\delta_{nn'} \frac{1}{N} \sum_k G_{ki}(\boldsymbol{k}, i\omega_{n}) G_{jl}(\boldsymbol{k} + \boldsymbol{q}, i\omega_n + i\Omega_{m}).

Here, :math:`N` denotes the number of lattice sites, and :math:`G_{ki}(\boldsymbol{k}, i\omega_n)` is the single-particle Green function within DMFT.
:math:`\Gamma_{\text{loc}}(i\Omega_m)` is the vertex part. The vertex part is local (no dependence on :math:`\boldsymbol{q}`) in DMFT. The vertex part can be obtained from the local two-particle Green's function :math:`\hat{X}_{\text{loc}}(i\Omega_m)` defined in the effective impurity problem in the DMFT:

.. math::

    \hat{X}_\text{loc}^{-1}(i\Omega_{m}) = \hat{X}_\text{0,loc}^{-1}(i\Omega_{m}) - \hat{\Gamma}_{\text{loc}}(i\Omega_{m}).

Eliminating :math:`\Gamma_{\text{loc}}(i\Omega_{m})` from the above two BS equations, we obtain the final equation

.. math::
    :label: Xq

    \hat{X}(\boldsymbol{q}, i\Omega_{m}) = \left[ \hat{X}_\text{loc}(i\Omega_{m})^{-1} - \hat{X}_\text{0,loc}(i\Omega_{m})^{-1} + \hat{X}_{0}(\boldsymbol{q}, i\Omega_{m})^{-1} \right]^{-1}.

Implementation notes:

- :math:`\hat{X}_0(\boldsymbol{q}, i\Omega_{m})` [Eq. :eq:`X0q`] is computed in ``dcore_chiq.py``.
- :math:`\hat{X}_\text{loc}(i\Omega_{m})` is computed in the impurity solver, which is called from ``dcore_chiq.py``.
- :math:`\hat{X}(\boldsymbol{q}, i\Omega_{m})` [Eq. :eq:`Xq`] and :math:`\hat{\chi}(\boldsymbol{q}, i\Omega_{m})` [Eq. :eq:`X_to_chi`] are computed in C++ part, which is called from ``chiq_main.py``.
- The fermionic Matsubara frequencies are truncated. The number of fermionic Matsubara frequency is specified by the parameter ``num_wf``. The maximum frequency is :math:`(2 \mathtt{num\_wf} + 1) \pi T`.
- The size of the matrices is :math:`(2N_{\text{orb}})^2 \times \mathtt{num\_wf}`, where :math:`N_{\text{orb}}` is the number of the total orbitals of the correlated shells in the unit cell.


.. _Algorithm_SCL:

Strong-coupling limit formula (SCL)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The desciption below is based on Section III.C in Ref. :ref:`[Otsuki et al. 2024] <references>`.

The strong-coupling limit formula (SCL) is a simplified approach to compute the momentum-dependent static susceptibilities in DMFT. The derivation is based on the strong-coupling feature of the two-particle Green's function. The validity of the SCL formula has been demonstrated for the single-orbital and two-orbital Hubbard models :ref:`[Otsuki et al. 2019] <references>` and for the multipolar ordering in CeB\ :sub:`6` :ref:`[Otsuki et al. 2024] <references>`.

In the original paper :ref:`[Otsuki et al. 2019] <references>`, three levels of approximations termed SCL1, SCL2, and SCL3 were proposed. **ChiQ** package provides SCL3, which is the simplest method among the three and does not require the vertex part.

The momentum-dependent static (:math:`\Omega_m=0`) susceptibility in the SCL formula is given by

.. math::
    :label: chiq_SCL

    \hat{\chi}_\text{SCL}(\boldsymbol{q}) = \left[ \hat{\chi}_\text{loc}^{-1} - \hat{I}_\text{SCL}(\boldsymbol{q}) \right]^{-1},

where hat indicates the matrix form defined as

.. math::

    [\hat{\chi}(\boldsymbol{q})]_{(ij),(kl)} \equiv \chi_{ij,kl}(\boldsymbol{q}).

:math:`\hat{\chi}_\text{loc}` is the local susceptibility computed in the effective impurity model.
:math:`\hat{I}_\text{SCL}(\boldsymbol{q})` is the effective intersite interaction defined by

.. math::
    :label: Iq_SCL

    \hat{I}_\text{SCL}(\boldsymbol{q}) = T \sum_{n} \phi(i\omega_n)^2 \hat{\Lambda}(i\omega_n, \boldsymbol{q}).

The function :math:`\hat{\Lambda}(i\omega_n, \boldsymbol{q})` is computed from the bare two-particle Green's function in the effective impurity and lattice models:

.. math::

    \hat{\Lambda}(i\omega_n, \boldsymbol{q}) = \hat{X}_\text{0,loc}(i\omega_n)^{-1} - \hat{X}_{0}(i\omega_n, \boldsymbol{q})^{-1}.

:math:`\phi(i\omega_n)` is given by

.. math::
    :label: phi

    \phi(i\omega_n) = \frac{1}{i\omega_n + \Delta_{-}} - \frac{1}{i\omega_n - \Delta_{+}},

where :math:`\Delta_{-}` and :math:`\Delta_{+}` are the excitation energy from :math:`n`-electron state to  :math:`(n-1)`-electron states and :math:`(n+1)`-electron states, respectively.
In the single-band Hubbard model, :math:`\Delta_{-}` and :math:`\Delta_{+}` correspond to the position of the lower and upper Hubbard bands, respectively.

Implementation notes:

- :math:`\hat{X}_0(\boldsymbol{q}, i\Omega_{m})` [Eq. :eq:`X0q`] is computed in ``dcore_chiq.py``.
- :math:`\hat{\chi}_\text{SCL}(\boldsymbol{q})` and :math:`\hat{I}_\text{SCL}(\boldsymbol{q})` [Eqs. :eq:`chiq_SCL` and :eq:`Iq_SCL`] are computed in the script ``calc_Iq_scl.py`` (``chiq_main.py`` is not used).
- The value of :math:`\Delta_{-}` and :math:`\Delta_{+}` are specified by the parameter ``delta_minus`` and ``delta_plus``, respectively.
- The size of the matrices is :math:`(2N_{\text{orb}})^2`.


.. _Algorithm_RPA:

Random phase approximation (RPA)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The random-phase approximation (RPA) is a simple weak-coupling approximation to compute the momentum-dependent susceptibility. The RPA susceptibility is given by

.. math::
    :label: chiq_RPA

    \hat{\chi}_\text{RPA}(\boldsymbol{q}, i\Omega_m) = \left[ \hat{\chi}_0(\boldsymbol{q}, i\Omega_{m})^{-1} - \hat{\gamma}_0 \right]^{-1},

where :math:`\hat{\chi}_0(\boldsymbol{q}, i\Omega_{m})` is the bare susceptibility, which is calculated from :math:`\hat{X}_0` in Eq. :eq:`X0q` by

.. math::
    :label: chi0q

    [\hat{\chi}_0(\boldsymbol{q}, i\Omega_{m})]_{(ij),(kl)} = T \sum_{n} [\hat{X}_0(\boldsymbol{q}, i\Omega_{m})]_{(ijn),(kln)}.

:math:`\hat{\gamma}` is the symmetrized vertex part given by

.. math::

    [ \hat{\gamma}_0 ]_{(ij),(kl)} = U_{ijkl} - U_{ijlk}.

Implementation notes:

- :math:`\hat{\gamma}_0` is generated in ``dcore_chiq.py``.
- :math:`\hat{X}_0(\boldsymbol{q}, i\Omega_{m})` [Eq. :eq:`X0q`] is computed in ``dcore_chiq.py``.
- :math:`\hat{\chi}_\text{RPA}(\boldsymbol{q}, i\Omega_{m})` and :math:`\hat{\chi}_0(\boldsymbol{q}, i\Omega_{m})` [Eqs. :eq:`chiq_RPA` and :eq:`chi0q`] are computed in C++ part, which is called from ``chiq_main.py``.
- The size of the matrices is :math:`(2N_{\text{orb}})^2`.


.. _Algorithm_RRPA:

Renormalized random-phase approximation (RRPA)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This part is based on Section II.B in Ref. :ref:`[Otsuki et al. 2019] <references>`.

The RPA above overestimates the fluctuations. Improvement to the RPA is made by including correlations from DMFT. By taking the Matsubara summation independently for each :math:`X` in Eq. :eq:`Xq`, we obtain the renormalized RPA (RRPA) susceptibility as

.. math::
    :label: chi_RRPA

    \hat{\chi}_\text{RRPA}(\boldsymbol{q}, i\Omega_{m}) = \left[ \hat{\chi}_{0}(\boldsymbol{q}, i\Omega_{m})^{-1} - \hat{\chi}_\text{0,loc}(i\Omega_{m})^{-1} + \hat{\chi}_\text{loc}(i\Omega_{m})^{-1} \right]^{-1}.

Compared with the RPA formula in Eq. :eq:`chiq_RPA`, the bare interaction :math:`\hat{\gamma}` is replaced by :math:`\hat{\chi}_\text{0,loc}(i\Omega_{m})^{-1} - \hat{\chi}_\text{loc}(i\Omega_{m})^{-1}`, which takes the local correlations into account.

Implementation notes:

- :math:`\hat{X}_0(\boldsymbol{q}, i\Omega_{m})` [Eq. :eq:`X0q`] is computed in ``dcore_chiq.py``.
- :math:`\hat{\chi}_\text{RRPA}(\boldsymbol{q}, i\Omega_{m})` and :math:`\hat{\chi}_0(\boldsymbol{q}, i\Omega_{m})` [Eqs. :eq:`chi_RRPA` and :eq:`chi0q`] are computed in C++ part, which is called from ``chiq_main.py``.
- The size of the matrices is :math:`(2N_{\text{orb}})^2`.


.. _Algorithm_Iq:

Intersite interactions
----------------------

This part is based on Section III.B in Ref. :ref:`[Otsuki et al. 2024] <references>`.

We can estimate the intersite interactions inversly from :math:`\hat{\chi}(\boldsymbol{q})` computed by solving the BSE. Assuming that the final results for :math:`\hat{\chi}(\boldsymbol{q})` was obtained from a localized model within the mean-field approximation, we can define the momentum-dependent interactions :math:`\hat{I}(\boldsymbol{q})` by

.. math::
    :label: Iq

    \hat{I}(\boldsymbol{q}) = \hat{\chi}_\text{loc}^{-1} - \hat{\chi}(\boldsymbol{q})^{-1}.

We note that, in the SCL, we directly compute :math:`\hat{I}(\boldsymbol{q})` by Eq. :eq:`Iq_SCL`.

Performing the Fourier transform, we obtain the interaction coefficients in the real space

.. math::
    :label: Ir

    \hat{I}(\boldsymbol{r}) = \frac{1}{N} \sum_{\boldsymbol{q}} \hat{I}(\boldsymbol{q}) e^{i\boldsymbol{q} \cdot \boldsymbol{r}},

where the summation of :math:`\boldsymbol{q}` is taken over the Brillouin zone.
This quantity defines the effective localized Hamiltonian

.. math::

    \mathcal{H}_\text{eff} = -\sum_{\boldsymbol{r} \boldsymbol{r}'} \sum_{ijkl} O_{ij}(\boldsymbol{r}) I_{ij,kl}(\boldsymbol{r} - \boldsymbol{r}') O_{lk}(\boldsymbol{r}').

Implementation notes:

- :math:`\hat{I}(\boldsymbol{q})` [Eq. :eq:`Iq`] is computed in ``calc_Iq.py``.
- Calculation of Eq. :eq:`Iq` is unstable when the charge fluctuation is tiny. The script ``calc_Iq.py`` takes a special care of the charge fluctuation to avoid the instability. See Appendix B in Ref. :ref:`[Otsuki et al. 2024] <references>` for details.
- ``chiq_main.py`` also output :math:`\hat{I}(\boldsymbol{q})`, but using ``calc_Iq.py`` is recommended because the special treatment of the charge fluctuation is done only in ``calc_Iq.py``.
- :math:`\hat{I}(\boldsymbol{r})` [Eq. :eq:`Ir`] is computed in ``bse_fft.py``.


.. _Algorithm_Eigen:

Analyzing the susceptibility
----------------------------

:math:`\chi_{ij,kl}(\boldsymbol{q}, i\Omega_m)` is a complicated object with many components. This section describes how to analyze it to extract the physical information. For simplicity, we consider only the static component, :math:`\Omega_{m}=0`.

.. _Algorithm_Eigen_values:

Eigenvalues
~~~~~~~~~~~

A simple way is to diagonalize :math:`\chi_{ij,kl}(\boldsymbol{q})` by taking the combined index :math:`(ij)` as row and :math:`(kl)` as column.

.. math::
    :label: chi_eigen

    \chi_{\xi}(\boldsymbol{q}) = \sum_{ijkl} U^{(\xi)}_{ij}(\boldsymbol{q})^{\ast} \chi_{ij,kl}(\boldsymbol{q}) U^{(\xi)}_{kl}(\boldsymbol{q}).


The eigenvalues :math:`\chi_{\xi}(\boldsymbol{q})` represent physical susceptibilities. The eigenvectors represent the physical quantitiy corresponding to the fluctuation.

Impelementation notes:

- ``chiq_post.py`` outputs eigenvalues by default (``mode = 'eigen'``).
- The eigenvalues are sorted in descending order by default. This behavior can be changed by parameter ``order``.
- The eigenvectors are output by setting ``vector = true`` (default is ``false``).

.. _Algorithm_Eigen_linear_combination:

Linear combinations
~~~~~~~~~~~~~~~~~~~

We can define the spin, charge, orbital, or multipole operator by taking linear comination of the density operator :math:`O_{\gamma} (\boldsymbol{r})` as

.. math::

   O_{\gamma}(\boldsymbol{r}) = \sum_{ij} C^{(\gamma)}_{ij} O_{ij}(\boldsymbol{r}).

With this coefficient :math:`C^{(\gamma)}_{ij}`, we can transform :math:`\chi_{ij,kl}(\boldsymbol{q})` as follows:

.. math::
    :label: chi_linear_combination

    \chi_{\gamma\gamma'}(\boldsymbol{q}) = \sum_{ijkl} C^{(\gamma)}_{ij} \chi_{ij,kl}(\boldsymbol{q}) C^{(\gamma') \ast}_{kl}.

Implementation notes:

- ``chiq_post.py`` outputs the diagonal componenents :math:`\chi_{\gamma\gamma}(\boldsymbol{q})` when ``mode = 'linear_combination'``.
- The coefficients :math:`C^{(\gamma)}_{ij}` are input by an external file.
- The number of coefficients can be smaller than the total number of operators.


.. _references:

References
----------

- \J. Otsuki, K. Yoshimi, H. Shinaoka, and Y. Nomura, "Strong-coupling formula for momentum-dependent susceptibilities in dynamical mean-field theory", `Phys. Rev. B, 99, 165134 (2019) <https://doi.org/10.1103/PhysRevB.99.165134>`_.
- \J. Otsuki, K. Yoshimi, H. Shinaoka, and H. Jeschke, "Multipolar ordering from dynamical mean field theory with application to CeB\ :sub:`6`\ ", `Phys. Rev. B, 110, 035104 (2024) <https://doi.org/10.1103/PhysRevB.110.035104>`_.


