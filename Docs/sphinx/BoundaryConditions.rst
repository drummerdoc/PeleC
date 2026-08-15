
 .. role:: cpp(code)
    :language: c++


.. _BCs:

Boundary Conditions
-------------------

PeleC manages boundary conditions in a form consistent with many AMReX codes. Ghost cell data are updated over an AMR level during a ``FillPatch`` operation and fluxes are then computed over the entire box without specifically recognizing boundary cells. A generic boundary filler function fills standard boundary condition types that do not require user input, including:

* *Interior* - Copy-in-intersect in index space (same as periodic boundary conditions). Periodic boundaries are set in the PeleC inputs file
* *Symmetry* - All conserved quantities and the tangential momentum component are reflected from interior cells without
  sign change (REFLECT_EVEN) while the normal component is reflected with a sign change (REFLECT_ODD)
* *NoSlipWall* - REFLECT_EVEN is applied to all conserved quantities except for both tangential and normal momentum components which are updated
  using REFLECT_ODD
* *SlipWall*  - SlipWall is identical to Symmetry
* *FOExtrap* - First-order extrapolation: the value in the ghost-cells are a copy of the last interior cell.

More complex boundary conditions require user input that is prescribed explicitly. Boundaries identified as ``UserBC`` or ``Hard`` in the inputs will be tagged as ``EXT_DIR`` in ``pc_hypfill``.  Users will then fill the boundary values, by calling the helper function, ``bcnormal``. The ``bcnormal`` function fills an exterior (ghost) cell based on the value of the outermost interior cell. Its arguments include a problem-specific data structure, the location, direction, and orientation of the boundary being filled, and potentially fluctuating turbulent velocities from the `TurbInflow <https://amrex-combustion.github.io/PelePhysics/Utility.html#turbulent-inflows>`_ utility in PelePhysics. See the ``Exec/RegTests/EB-ConvergingNozzle`` case for an example of how to use the TurbInflow utility. This gives the user flexibility to specify a variety of boundary conditions, including faces that contain both walls and inflow regions. Note that the external state ``s_ext`` is prepopulated with the same values as are used for the ``NoSlipWall`` condition, so the default if the ``bcnormal`` function does nothing is to specify a ``NoSlipWall``.

.. note::
   To ensure conservation, when Godunov schemes are used the order of accuracy is reduced at boundaries specified using ``bcnormal``; PLM is used and the predictor step is omitted when computing fluxes through these boundaries. This does not affect any other boundary types or simulations using MOL.

Special care should be taken when prescribing subsonic ``Inflow`` or an ``Outflow`` boundary conditions. It might be tempting to directly impose target values in the boundary filler function (for ``Inflow``), or to perform a simple extrapolation (for ``Outflow``).  However, this approach would fail to correctly respect the flow of information along solution characteristics - the system would be ill-posed and would lead to unphysical behavior. In particular, at a subsonic inflow boundary, at a subsonic inlet there is one outgoing characteristic, so one flow variable must be specified using information from inside the domain. Similarly, there is one incoming characteristic at outflow boundaries. The :ref:`NSCBC method <NSCBC>` described below is the preferred way to account for this. A simpler alternative, adequate when the boundary is far from any acoustic source of interest, is:

* Subsonic Inflows: Specify the desired temperature, velocity, and composition (if relevant) in the ghost cells. Take the pressure from the domain interior. Based on these values, compute the density, internal energy, and total energy for the ghost cells.

* Subsonic Outflows: Specify the desired outlet pressure and extrapolate the other flow quantities. In particular, we recommend following the simple characteristic-based extrapolation proposed by Whitfield and Janus (`Three-Dimensional Unsteady Euler Equations Solution Using Flux Vector Splitting. AIAA Paper 84-1552, 1984. <https://arc.aiaa.org/doi/abs/10.2514/6.1984-1552>`_) and described in Ch. 8 of Blazek's textbook (`Computational Fluid Dunamics - Principles and Applications <https://www.sciencedirect.com/book/9780080445069/computational-fluid-dynamics-principles-and-applications>`_). Implementations of this method can be found in the ``bcnormal`` function in ``prob.H`` for various test cases, including EB-C10 and EB-ConvergingNozzle.

A detailed analysis comparing various boundary condition strategies and demonstrating their implementation is available for the :ref:`EB-ConvergingNozzle` case.

Isothermal Walls
~~~~~~~~~~~~~~~~

By default, the boundaries specified as ``NoSlipWall`` and ``SlipWall`` are adiabatic. For isothermal wall boundaries, energy fluxes through the isothermal wall are computed separately, rather than being based on values populated in the ghost cells. To activate computation of isothermal wallfluxes, use the input file option ``pelec.do_isothermal_walls = 1`` and then specify the desired wall temperatures using, for example, ``pelec.domlo_wall_temp = -1 -1 300.0`` and ``pelec.domhi_isothermal_temp = -1 -1 400.0``, which would leave the x and y boundaries as adiabatic, make the lower z boundary isothermal at 300 K, and make the upper z boundary isothermal at 400 K. Any boundary with a negative (or zero) value for the specified temperature is treated as adiabatic; boundaries that are not ``NoSlipWall``, ``SlipWall``, ``UserBC``, or ``Hard`` must always have a negative value specified.

.. _NSCBC:

Navier-Stokes Characteristic Boundary Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Directly imposing a target state at a subsonic inflow, or extrapolating at a subsonic outflow, does not respect the
flow of information along the solution characteristics: the resulting problem is ill-posed and acoustic waves reflect
off the boundary back into the domain. The remedy is to decompose the boundary state into characteristic waves, take
the outgoing ones from the interior, and model only the incoming ones. This is the Navier-Stokes Characteristic
Boundary Condition (NSCBC) strategy of `Poinsot and Lele (1992) JCP
<https://www.sciencedirect.com/science/article/pii/0021999192900462>`_, in the ghost-cell form of `Motheau et al.
(2017) AIAA Journal <https://ccse.lbl.gov/people/motheau/Manuscripts_website/2017_AIAA_CFD_Motheau.pdf>`_.

Enable it with::

    pelec.bc_nscbc = 1

and provide a ``bcnormal_nscbc`` hook in the problem's ``prob.H`` (see below). The switch alone does nothing: the
treatment is applied only at boundary *points* for which that hook returns a live target, so a problem that does not
provide it is unaffected. A face must be ``Hard`` or ``UserBC`` for the treatment to reach it.

.. note::

   By default (``bc_nscbc_beta = bc_nscbc_beta_s = 1``) the implementation is the **1-D-normal LODI limit**: the
   modelled incoming waves carry neither transverse nor reaction terms. Setting :math:`\beta < 1` enables the
   transverse terms and is worth a factor of a few whenever the boundary is met obliquely — see *Transverse terms*
   below. :math:`\beta_s < 1` enables the reaction source. The **viscous** contribution to the incoming wave is
   genuinely absent. Place an outflow away from flames regardless: see *What a front crossing the outflow actually
   costs*, which makes that advice quantitative rather than folkloric.

Transverse terms
""""""""""""""""""""

The 1-D characteristic decomposition assumes the wave meets the boundary head-on. When it does not, the
decomposition attributes part of a purely *outgoing* wave to the incoming characteristic, and the boundary then
reflects it. The transverse terms correct for that:

.. math::

   L_{\rm in} = K (\phi - \phi_{\rm target}) - (1 - \beta) \, T_{\rm transverse}

For a plane wave at angle :math:`\theta` to the boundary normal the correction that exactly cancels the
obliqueness error works out to :math:`(1-\beta) = \cos\theta / (1 + \cos\theta)`, so

.. math::

   \beta_{\rm opt} = 1 - \frac{\cos\theta}{1 + \cos\theta}

which is **0.5 at normal incidence** and rises toward 1 as the wave becomes grazing (0.67 at 60°, 0.79 at 75°).
Both measured optima in ``Exec/RegTests/NSCBC-COVO`` land on that curve.

.. table::

   ==========================  ==========================  =============================
   :math:`\beta`               convected vortex, rms        circular pulse, front
                               residual vs no-boundary      amplitude spread
   ==========================  ==========================  =============================
   hard boundary               14.0×                        21.4%
   1.0 (transverse off)        10.0×                        0.90%
   0.8                         7.3×                         **0.066%**
   0.5                         **3.8×**                     1.08%
   0.2 (= local Mach)          16.1×                        —
   0.0 (full)                  132×                         2.95%
   ==========================  ==========================  =============================

Note the shape of that curve: too little correction costs a factor of a few, **too much is catastrophic**.
:math:`\beta` below about 0.3 is worse than having no transverse terms at all, and :math:`\beta = 0` is close to
unstable. That asymmetry is why the default is 1 rather than the optimum — a wrong :math:`\beta` is far more
dangerous than an absent one. Start at 0.5 and raise it toward 1 if the boundary is mostly met at grazing
incidence.

A negative :math:`\beta` selects the local Mach number, which is the legacy Fortran's convention. The measurements
do not support it: at M = 0.2 it gives :math:`\beta = 0.2`, on the wrong side of the minimum.

The reaction source term
""""""""""""""""""""""""

A flame near an outflow raises the pressure in the boundary cell at a rate the 1-D relaxation has no model for, so
:math:`\sigma` has to absorb it and the mean pressure sits off target. ``pelec.bc_nscbc_beta_s`` adds the missing
term, with the same convention as :math:`\beta`.

It is worth knowing what this quantity is *not*. Chemistry conserves mass, and PeleC's constant-volume reactor
conserves total internal energy because the formation enthalpies are carried inside :math:`e` — so the reaction
contributes neither a density source nor a total-energy source, and ``srcq(QPRES)`` (which is built only from those
two) is essentially zero for pure chemistry. The entire effect is compositional:

.. math::

   \left.\frac{\partial p}{\partial t}\right|_{\rm react}
   = \sum_k \left.\frac{\partial p}{\partial Y_k}\right|_{\rho,e} \frac{\dot\omega_k}{\rho}
   = R_u T \sum_k \frac{\dot\omega_k}{W_k} \;-\; \frac{p}{\rho T c_v}\sum_k e_k \dot\omega_k

a mole-change term plus the familiar heat-release term. That closed form is exact for an ideal-gas mixture. For a
real-gas EOS there is no such one-liner, so a directional finite difference along the reaction path is used instead;
because :math:`\sum_k \dot\omega_k = 0` the perturbation already preserves :math:`\sum_k Y_k = 1`, so this costs a
single extra :math:`(\rho, e, Y) \rightarrow T` solve rather than one per species. The two agree to eight
significant figures (check C7 in ``Verification/NSCBC1D``).

The correct first response to heat release at an outflow is still to **move the boundary**; this term is for when
that is not possible. How little it buys you when that is not possible is the subject of the next section.

What a front crossing the outflow actually costs
""""""""""""""""""""""""""""""""""""""""""""""""

``Exec/RegTests/NSCBC-FlameOutflow`` parks a wrinkled H\ :sub:`2`/air flame sheet on the outflow plane so that half
the boundary is burnt, half is fresh, and the reaction zone lies in the boundary cells at the two crossings. Errors
below are against a reference whose outflow is 2.4 cm further downstream, outside the domain of dependence of the
comparison region, so the difference is the outflow's own error. :math:`R` is the acoustic reflection at the same
:math:`\sigma` from ``Verification/NSCBC1D``.

.. table::

   ==============================  ==========  ==========  =============  ==========
   Outflow                         :math:`\sigma`  :math:`\beta_s`  mean :math:`\Delta p`  :math:`R` [%]
   ==============================  ==========  ==========  =============  ==========
   hard ``p = p_amb``              --          --          --527          --
   characteristic                  0.25        0           +2185          0.76
   characteristic                  1           1           +2065          2.56
   characteristic                  1           0           +2074          2.56
   characteristic                  4           0           +1450          7.20
   characteristic                  16          0           **+218**       28.1
   ``pin_farfield``                --          0           +748           0.015
   ==============================  ==========  ==========  =============  ==========

in dyn/cm\ :sup:`2`, against :math:`p_{\rm amb} = 1.013\times10^6`. Three conclusions, none of them comfortable.

**The reaction source term is not what matters here.** :math:`\beta_s` moves the error by under 1%, and the whole
error is still present with ``pelec.do_react = 0``. The dominant mechanism is the ghost pressure. Per ghost layer
:math:`\ell`,

.. math::

   p_g - p_N = \underbrace{\tfrac{1}{2}\rho c\, \ell\, \delta R_+}_{\rm extrapolation}
             - \underbrace{\frac{\ell\,\sigma}{2 n_x}\,(p_N - p_\infty)}_{\rm anchoring} ,

and the anchoring increment carries a factor :math:`1/n_x` that the extrapolation term does not. Balancing them,

.. math::

   \Delta p \;\approx\; \rho\, c\, L_{\rm ref} \left. \frac{\partial u_{\rm out}}{\partial n}\right|_{b} \Big/ \sigma .

A **normal velocity gradient** at the boundary biases the ghost pressure, and :math:`\sigma` is the only thing
fighting it. A flame crossing the outflow is a large such gradient — the gas accelerates from 45.6 to 123.8 cm/s
across 0.07 cm — and it is dilatational, not acoustic, so the invariant algebra is not at fault: LODI has no way to
tell the two apart. The formula predicts the measured :math:`\sigma^{-1}` trend, and predicts that the offset is
grid-*converged* rather than refining away, which it is (:math:`n_x` = 48, 96, 192 give +531, +409, +371).

**The inert default of** :math:`\sigma = 0.25` **is the worst available choice** in this situation, and only around
:math:`\sigma \approx 10` does the characteristic outflow beat a hard pressure outflow — at roughly 20% acoustic
reflection. You are choosing which error to have, not removing one.

**``order = 2`` is load-bearing, not a refinement.** First order flips the sign and gives an error seven times
larger than :math:`\sigma = 16`: extrapolating the outgoing invariant is what lets the front's structure leave.

What is specified, and what is not
""""""""""""""""""""""""""""""""""

In an outward-normal frame the wave speeds are :math:`u_n - c`, :math:`u_n` (with multiplicity
:math:`D - 1 + N_{\rm spec} + \ldots`) and :math:`u_n + c`. At a subsonic boundary :math:`u_n + c` always leaves and
:math:`u_n - c` always enters; the :math:`u_n` family leaves at an outflow and enters at an inflow. Hence:

* **Subsonic outflow** has exactly one incoming characteristic, so exactly one quantity may be specified: the
  pressure. Density, temperature, velocity, mass fractions and passive scalars are all *extrapolated*. In particular
  **composition must not be imposed at an outflow** — doing so over-specifies the problem and plants a spurious
  composition front on the boundary.
* **Subsonic inflow** has :math:`D + N_{\rm spec} + 1` incoming characteristics: velocity, temperature and
  composition are specified, and only the outgoing acoustic comes from the interior.
* **Supersonic outflow** needs nothing; a zero-gradient copy is exact.
* **Supersonic inflow** is a full Dirichlet condition.

The implementation dispatches on the local Mach number and handles all of these, including transient flow reversal
through a face, without any user intervention.

Controls
""""""""

There are three physical dials, and for the commonest case — one non-reflecting subsonic outflow — none of them
needs to be changed.

.. table::

   ==================================  =========  ==========================================================
   Parameter                           Default    Meaning
   ==================================  =========  ==========================================================
   ``pelec.bc_nscbc``                  0          Master switch.
   ``pelec.bc_nscbc_sigma``            0.25       Outflow pressure relaxation (Poinsot-Lele :math:`\sigma`).
   ``pelec.bc_nscbc_relax_u``          2.0        Inflow normal-velocity relaxation.
   ``pelec.bc_nscbc_relax_t``          0.2        Inflow temperature and tangential-velocity relaxation.
   ``pelec.bc_nscbc_beta``             1.0        Transverse-term weight. 1 = off; **try 0.5**.
   ``pelec.bc_nscbc_beta_s``           1.0        Reaction-source weight. 1 = off.
   ``pelec.bc_nscbc_order``            2          Extrapolation order, 1 or 2. Verification knob only.
   ``pelec.bc_nscbc_pin_farfield``     0          Pin the incoming acoustic instead of relaxing it.
   ``pelec.bc_nscbc_extrap_temperature`` 0        Give the outflow face a controlled conductive flux.
   ``pelec.bc_nscbc_extrap_material``  0          Continue material structure through the outflow.
   ==================================  =========  ==========================================================

**All of these are positive.** The legacy Fortran implementation required ``relax_T`` to be negative and its
documentation and its code disagreed about the sign of ``sigma_out``; every internal sign is now handled by the
kernel, and a negative coefficient is a fatal error.

There is deliberately **no relaxation coefficient for composition**. At an inflow every species characteristic is
incoming, so hard imposition is well posed and needs no dial; at an outflow every species characteristic is
outgoing, so composition is extrapolated and there is nothing to relax.

There is also no separate reference-length parameter. The relaxation *rate* is

.. math::

   K = \sigma \, (1 - M^2) \, c / L , \qquad \tau_{\rm relax} = 1/K ,

with :math:`L` the domain extent along the boundary normal. Only the ratio :math:`\sigma/L` is physical, so exposing
both would be two dials for one degree of freedom. Dividing through by the acoustic transit time :math:`L/c` gives
the number that actually matters:

.. math::

   \frac{\tau_{\rm relax}}{\tau_{\rm acoustic}} = \frac{1}{\sigma (1 - M^2)} ,

so :math:`\sigma = 0.25` at low Mach means *the boundary pressure is pulled back to the target over about four
acoustic transit times*. That is the whole content of "0.25 is a good default", and it gives two conditions to
check: :math:`\tau_{\rm relax}` must be much larger than the acoustic transit time, so that outgoing waves leave
before the pressure is pulled back, and much smaller than the run time, so that the mean pressure does not drift.
PeleC prints both at startup. Rudy and Strikwerda (*JCP* **36**, 55, 1980) find :math:`\sigma \approx 0.27` optimal
for a 1-D pulse; the usual band is 0.15-0.3. :math:`\sigma = 0` is *perfectly* non-reflecting and leaves the mean
pressure entirely unanchored.

``bc_nscbc_pin_farfield`` replaces the relaxation by a hard pin of the incoming invariant to its far-field value.

``bc_nscbc_extrap_temperature`` addresses something the characteristic algebra does not: PeleC's diffusion operator forms
the conductive and species fluxes at a physical boundary face from the NSCBC ghost cells, so whatever normal temperature
gradient the ghost carries *is* the heat flux leaving the domain. By default the outflow closes its :math:`\lambda_0`
family on the linearised entropy invariant, and the ghost temperature is then whatever the EOS returns from the
extrapolated density and the acoustically-set pressure — chosen by nothing that has a heat flux in mind. Check C8 in
``Verification/NSCBC1D`` measures the consequence on a flame-like ramp: the face temperature gradient is overstated by
**40%**. Setting this to 1 extrapolates the temperature on the same limited slope as everything else and derives the
density from the EOS, so the face gradient is the interior one exactly, at no cost on the hyperbolic side (a uniform
state is still reproduced to round-off).

It is off by default because it changes every outflow result. Turn it on when a thermal or compositional structure is
near the boundary. In ``Exec/RegTests/NSCBC-FlameOutflow`` it is worth 9% at :math:`\sigma = 1`, where the ghost-pressure
bias still dominates, and 69% at :math:`\sigma = 16`, where that bias is suppressed and the diffusive error is most of
what remains — taking the mean-pressure error to 67 dyn/cm², eight times *better* than a hard pressure outflow.

``bc_nscbc_extrap_material`` attacks the ghost-pressure bias itself rather than out-relaxing it. The default outflow
extrapolates the outgoing invariant :math:`R_+` with its full interior slope while the incoming :math:`R_-` carries only
the relaxation increment, so across a dilatational velocity gradient the ghost picks up the
:math:`\tfrac{1}{2}\rho c\,\ell\, \delta u` bias of check C9(a) — the dominant error at a front-crossing outflow. The
material part of the :math:`R_-` slope cannot be measured from :math:`R_-` directly (its gradient also contains every
incoming wave, including the one the relaxation itself launches — extrapolating that is positive feedback), so it is
bounded through the *entropy* family, which is outgoing at an outflow and carries no acoustic content: steady continuity
turns :math:`dS = d\rho - dp/c^2` into :math:`du_{\rm mat} = -u\, dS / (\rho\,(1 - M^2))`, and the applied slope is
``minmod`` of the measured :math:`R_-` slope and that bound. Across a flame the two agree and the ghost continues the
interior's :math:`u` and :math:`p` slopes exactly, killing the bias at *any* :math:`\sigma`; for pure acoustics the
entropy bound vanishes and nothing changes (reflection, anchoring and the relaxation rate are unmoved in
``Verification/NSCBC1D``). Check C11 gates it on a manufactured flame — a sustained velocity/density ramp with an exact
steady solution straddling the outflow: the entropy closure drifts by :math:`15707\ {\rm dyn/cm^2}` while the boundary
equilibrates and distorts :math:`\partial u/\partial n` at the face to 175% of its exact value; with this flag the drift
is :math:`3175` and the face gradient stays within the band a ghost fill of the *exact* continuation also occupies.

Measured in ``Exec/RegTests/NSCBC-FlameOutflow`` (quasi-frozen front on the outflow, :math:`\sigma = 1`): 5% alone —
which settles the open question from the C9/C10 commits: the extrapolation bias is real and exactly reproducible, but it
is **not** what dominates a reacting front-crossing outflow, whose residual points at the unmodelled diffusive and
reactive enthalpy deposition in the boundary cells. Combined with ``bc_nscbc_extrap_temperature``, however, the two
closures compound to **42%** — together they make the ghost a fully consistent material continuation.

Its sharp edge is a front that *leaves*. The continuation is a quasi-steady model, and during a fast transit at small
:math:`\sigma` it turns the crossing into a runaway (see *A structure leaving through the outflow* below): use it for
fronts that sit near the boundary, keep :math:`\sigma > 0` always — the relaxation increment is still the only
anchoring — and prefer it off, or paired with a large :math:`\sigma`, when structures actually cross.

A structure leaving through the outflow
"""""""""""""""""""""""""""""""""""""""

``Exec/RegTests/NSCBC-FlameOutflow`` (``nscbc-flameexit.inp``) measures a wrinkled flame sheet *actually leaving* at
U = 40 :math:`S_L`, with the passage complete inside the run and the post-exit truth known exactly (a uniform fresh
stream). The result inverts the acoustic intuition this page is otherwise built on:

* A fast dilatational transit wants **anchoring, not transparency**. What crosses is mass and enthalpy at low Mach,
  with essentially no acoustic content, and every treatment ranks by anchoring strength: the hard Dirichlet is
  near-perfect; :math:`\sigma = 16` + ``extrap_temperature`` exits the flame with front kinematics and wrinkle
  amplitude indistinguishable from it (transient 0.18% of ambient, final state exact to 8 ppm) — and the ~28% acoustic
  reflection that :math:`\sigma = 16` costs elsewhere never appears, because nothing acoustic is present.
* The inert default :math:`\sigma = 0.25` is **unstable** for a transit: the anchoring time is comparable to the
  crossing, the flux imbalance integrates into a pressure ramp that pushes the front back upstream, pumps the wrinkle,
  and the run dies. Raise :math:`\sigma` to O(10) *before* a front reaches the outflow.
* ``pin_farfield`` is stable but anchors a through-flow duct to :math:`p + \rho c u` — the operating point shifts and
  stays shifted. Reservoir boundaries only.
That is simultaneously non-reflecting and anchored, which a rate-based relaxation cannot be, but it anchors to
:math:`p_{\rm target} + \rho c u_n` rather than to :math:`p_{\rm target}`, and because it constrains a *value*
rather than a rate its effective relaxation rate is :math:`c/\Delta x` and it does not converge under mesh
refinement. Use it for an open boundary onto a large quiescent reservoir; do not use it for a duct exhausting into
a plenum whose true mean pressure is not the target.

Providing a target
""""""""""""""""""

Override ``bcnormal_nscbc`` in the problem's ``ProblemSpecificFunctions``. The decision is made per boundary
*point*, so a single face may mix an inflow, an outflow and a wall::

    struct MyProblemSpecificFunctions : public DefaultProblemSpecificFunctions
    {
      AMREX_GPU_DEVICE
      AMREX_FORCE_INLINE
      static pc_nscbc::Target bcnormal_nscbc(
        const amrex::Real* x,
        const amrex::Real* /*s_int*/,
        const int idir,
        const int sgn,
        const amrex::Real /*time*/,
        amrex::GeometryData const& /*geomdata*/,
        ProbParmDevice const& prob_parm)
      {
        pc_nscbc::Target t;
        if (idir == 0 && sgn == -1) {          // x-hi: non-reflecting outflow
          t.type = pc_nscbc::Type::outflow;
          t.p = prob_parm.p_amb;               // the ONLY thing specified
        } else if (idir == 0 && sgn == +1) {   // x-lo: jet inlet in a wall
          if (std::abs(x[1]) < prob_parm.jet_radius) {
            t.type = pc_nscbc::Type::inflow;
            t.u[0] = prob_parm.u_jet;
            t.T = prob_parm.T_jet;
            for (int n = 0; n < NUM_SPECIES; n++) {
              t.Y[n] = prob_parm.Y_jet[n];
            }
          }
          // outside the jet, t.type stays `off` and the point falls through
          // to the ordinary bcnormal(), which can make it a wall
        }
        return t;
      }
    };
    using ProblemSpecificFunctions = MyProblemSpecificFunctions;

``x`` is the location on the boundary *plane*, not the ghost-cell centre: the target is a property of the boundary
point and must not vary with the ghost layer. Returning the default-constructed ``Target`` (type ``off``) leaves the
point to the ordinary ``bcnormal`` path.

Where to put the boundary
"""""""""""""""""""""""""

The quality of a characteristic boundary condition is set more by where the boundary is than by the dials.

* Put a subsonic outflow at least one acoustic wavelength of the lowest frequency of interest, and roughly ten flame
  thicknesses, downstream of any reaction zone. The modelled incoming wave carries no reaction-source term, so heat
  release in the boundary cell shows up as a mean pressure offset that :math:`\sigma` must absorb.
* Avoid placing an outflow across a strong shear layer, a vortex core or a composition front. The acoustic impedance
  :math:`\rho c` is frozen at the boundary cell, and that linearisation is weakest exactly there. More to the point,
  the mean-pressure error scales as :math:`\rho c L_{\rm ref} (\partial u_{\rm out}/\partial n) / \sigma` and does
  not refine away — see *What a front crossing the outflow actually costs* above. Anything that puts a normal
  velocity gradient on the boundary is expensive, and a flame is the most expensive of them.
* Do not refine an AMR level at an outflow face. The residual reflection is :math:`O(\Delta x^p)` and the
  extrapolation stencil is level-local, so a fine patch on the boundary introduces a level-dependent artefact.
* If a sponge or artificial-viscosity ramp is currently needed at an outflow, removing it is the acceptance test for
  the boundary condition, not something to keep for safety.

Diagnosing
""""""""""

Every fallback path in the kernel is counted, and the counts are reported alongside the other integrated quantities
when ``pelec.v > 0``. A silent fallback is a bug that will not be found.

.. table::

   ===============================================  ==========================================================
   Symptom                                          Action
   ===============================================  ==========================================================
   Outgoing pulse visibly reflects                  Reduce :math:`\sigma` toward 0.15; confirm ``order = 2``.
   Mean pressure drifts over a long run             Increase :math:`\sigma` so :math:`\tau_{\rm relax}` is
                                                    well below the run time; or move the outflow away from
                                                    the flame.
   Reflects *and* drifts                            The boundary is sitting on a flame or shear layer. No
                                                    dial fixes this; move it. If you cannot, raise
                                                    :math:`\sigma` to O(10) and accept the reflection, or
                                                    use ``pin_farfield`` and accept a fixed offset of order
                                                    :math:`\rho c\, u_{\rm out}`.
   A flame or front must pass OUT through the       Raise :math:`\sigma` to O(10) and set
   outflow                                          ``extrap_temperature = 1`` before it arrives. The
                                                    default :math:`\sigma = 0.25` can push the front back
                                                    and crash the run; ``extrap_material`` off during the
                                                    transit. See *A structure leaving through the outflow*.
   Mean pressure off target at a reacting outflow,  Expected. :math:`\beta_s` is a second-order correction
   and ``beta_s`` does not help                     on top of a first-order problem; the error is driven by
                                                    :math:`\partial u_{\rm out}/\partial n`, not by the heat
                                                    release. Raise :math:`\sigma`.
   Reflection persists at any :math:`\sigma`, 2-D    Missing transverse terms. Move the boundary further out.
   Inflow never reaches the target velocity         Increase ``relax_u``.
   Inflow oscillates or goes unstable               Reduce ``relax_u``; a value much above 10 is a hard
                                                    Dirichlet condition in disguise, which is ill-posed for
                                                    the number of incoming characteristics.
   Inlet temperature sags below target              Increase ``relax_t``.
   Turbulent inlet :math:`u'` below the prescribed   ``relax_u`` low-pass filters the injected spectrum.
   intensity                                        Increase it, or rescale the prescribed intensity to hit
                                                    the realised one.
   ``flow reversal`` count sustained above zero     The outflow is misplaced, or the face should be an
                                                    inflow. Transient counts are benign.
   ``EB body state`` count above zero               Covered cells are adjacent to the boundary face; the
                                                    stencil order was degraded to avoid them.
   ===============================================  ==========================================================

Verification
""""""""""""

``Verification/NSCBC1D`` compiles the production kernel unmodified against a minimal 1-D solver and asserts uniform-
state consistency, relaxation direction, species handling, acoustic reflection, grid-independence of the relaxation
rate, and every fallback path. It also produces the reference :math:`\sigma` sweep. Run it before trusting a change
to ``Source/NSCBC.H``.

