.. _TutorialNSCBC:

Characteristic boundary conditions, hands on
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``Exec/Tutorials/NSCBC`` is a ladder of small, preconfigured runs that shows what the characteristic (NSCBC)
boundary condition of :ref:`NSCBC` does, how each of its settings changes the answer, and where it stops being the
right tool. Every rung builds an existing case, runs a handful of variants, and plots one figure; the plot script
prints the measured numbers next to the ones below, so a user knows when the figure has been reproduced. The
figures here are those reference figures.

By the end a user should be able to choose :math:`\sigma`, ``relax_u``, ``relax_t`` and :math:`\beta` for a new
problem, recognise a boundary that is in the wrong place, know which structures may cross an outflow and which
may not, and know when the right answer is a longer domain rather than a better setting.

Each rung has a ``README.md`` with the full walk-through; what follows is the summary.

.. code-block:: sh

   cd Exec/Tutorials/NSCBC/01-outflow-pulse
   ./run.sh            # builds the case it needs, runs every variant into runs/
   python plot.py      # writes the figure and prints measured against reference

Prerequisites: a PeleC build environment, MPI, and Python 3 with ``numpy`` and ``matplotlib``. ``run.sh`` skips
variants that already finished; ``NPROCS`` (default 4) sets the ranks. Pressures are in dyn/cm²; the target is
1 atm = 1.01325 × 10⁶.

.. list-table::
   :header-rows: 1
   :widths: 8 62 30

   * - Rung
     - Question
     - Case
   * - 1
     - Why a characteristic outflow?
     - NSCBC-Acoustic, 2-D
   * - 2
     - :math:`\sigma`: letting waves out against holding the mean
     - NSCBC-Acoustic, 2-D
   * - 3
     - Subsonic inflow: ``relax_u``, forcing (``dudt``), the NRI register
     - NSCBC-Acoustic duct
   * - 4
     - A turbulent inflow through the characteristic inlet
     - NSCBC-Acoustic 3-D duct, TurbInflow
   * - 5
     - Oblique waves and vortices: the transverse weight :math:`\beta`
     - NSCBC-COVO
   * - 6
     - What may cross an outflow: hot fronts and flames
     - NSCBC-FlameOutflow
   * - 7
     - Flow entering through an outflow
     - NSCBC-Acoustic
   * - 8
     - When the answer is a longer domain
     - NSCBC-Chamber

Rung 1 — why a characteristic outflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: /images/nscbc_tutorial/01_outflow_pulse.png
   :align: center
   :figwidth: 100%

A 0.1% acoustic pulse leaves a 10 cm duct. Imposing the ambient pressure in the ghost cells reflects it inverted
and almost whole (R = 97%); the characteristic outflow at the default :math:`\sigma = 0.25` reflects 0.8% and
leaves a small uniform offset that the relaxation then removes.

Rung 2 — :math:`\sigma`
^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: /images/nscbc_tutorial/02_outflow_sigma.png
   :align: center
   :figwidth: 100%

The same duct, for :math:`\sigma` from 0 to 16: the reflection of a pulse (left) and the return of a uniform
overpressure to the target (right). :math:`\sigma = 0` is perfectly transparent and never anchors; raising
:math:`\sigma` anchors faster and reflects more (27.5% at 16), and from :math:`\sigma \approx 4` the duct rings.
The default 0.25 reflects under 1% and returns the mean without overshoot. The NDNR register breaks the trade
for sound: 0.9% at :math:`\sigma = 16`, and the mean back within 1% in 9 transits.

Rung 3 — the subsonic inflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: /images/nscbc_tutorial/03_inflow.png
   :align: center
   :figwidth: 100%

A step in the inlet velocity (left) settles in 5–8 acoustic round trips for ``relax_u`` between 0.3 and 1; too
soft never arrives, too stiff rings like a wall. A forced inlet (right) injects the requested wave only with the
analytic derivative ``Target::dudt`` *and* the NRI register (``bc_nscbc_nri = 1``): relaxation alone lags or
overdrives at every stiffness.

Rung 4 — a turbulent inflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: /images/nscbc_tutorial/04_turbulent_inflow.png
   :align: center
   :figwidth: 100%

TurbInflow fluctuations fed through the characteristic inlet (the hook adds ``turb_fluc`` to its target). The
defaults pass about half the streamwise and under a third of the tangential intensity; ``relax_u = relax_t = 20``
pass 0.93 and 0.90. A stiff ``relax_u`` reflects sound, so a turbulent inlet trades turbulence fidelity against
acoustic transparency. Do not combine NRI with a turbulent inflow: it amplifies the streamwise fluctuations (1.4×).

Rung 5 — the transverse weight :math:`\beta`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: /images/nscbc_tutorial/05_transverse.png
   :align: center
   :figwidth: 100%

A circular pulse leaving a box at every angle (left) and a vortex leaving through an outflow (right). The
transverse terms, weighted by :math:`1 - \beta`, remove most of what the 1-D-normal boundary does wrong:
:math:`\beta = 0.8` for sound arriving at shallow angles, :math:`\beta = 0.5` for a vortex (3.9× the
no-boundary floor, against 10× without transverse terms). Too much correction is catastrophic: :math:`\beta = 0`
goes unstable on the vortex.

Rung 6 — what may cross an outflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: /images/nscbc_tutorial/06_crossing.png
   :align: center
   :figwidth: 100%

A hot front, inert (A) and burning (B), carried out through the outflow. A front is dilatation, not a wave, and
the relaxation has no term for it. With all defaults a flame never leaves: the overpressure (5.8% of ambient)
reverses the outflow and pushes the flame back upstream, and the boundary cells heat to 3,400 K (lower right).
The flame closures (``extrap_temperature = 1``, ``beta_s = 0``) let it out on schedule at :math:`\sigma = 0.25`
with a 0.9% dip; :math:`\sigma = 16` lets it out alone (0.15%), and with the closures has about twice a hard
outflow's excursion. An inert hot
front needs stiff anchoring instead (:math:`\sigma \approx 16`): at :math:`\sigma = 0.25` the defaults leave a 2%
excursion and ``extrap_temperature = 1`` alone over-vents (−7.7%). NDNR does no harm while a front crosses
(−1,062 at :math:`\sigma = 16` with the closures, against −432 without).

.. list-table::
   :header-rows: 1
   :widths: 35 45 10

   * - Crossing the outflow
     - Settings
     - Rung
   * - sound, head-on
     - defaults; ``order = 2``
     - 1, 2
   * - sound at an angle, vortices, shear
     - ``beta = 0.5`` (0.8 for sound at shallow angles)
     - 5
   * - an inert hot or compositional front
     - :math:`\sigma \approx 16`, or a hard outflow if no sound must leave
     - 6
   * - a flame
     - ``extrap_temperature = 1``, ``beta_s = 0``
     - 6
   * - flow coming *in*
     - ``backflow_material = 1`` with ``t.T``, ``t.Y`` in the target
     - 7

Rung 7 — flow entering through an outflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: /images/nscbc_tutorial/07_backflow.png
   :align: center
   :figwidth: 100%

A 600 K duct draws 300 K air in through its outflow. Brief reversals need nothing; sustained backflow by default
hands back the interior's own material, so the duct draws its own exhaust and never cools.
``backflow_material = 1``, with the reservoir's ``T`` and ``Y`` in the outflow target, flushes it in 17 ms. A
``flow reversal`` count that keeps growing in the fallback report means a face in sustained backflow.

Rung 8 — when the answer is a longer domain
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: /images/nscbc_tutorial/08_longer_domain.png
   :align: center
   :figwidth: 100%

A closed-end H₂/air chamber whose flame pushes fresh gas out through a vent, with the characteristic boundary at
the vent and with the duct continued 2.4 cm beyond it (air-filled, ``prob.charge_end``), each at
:math:`\sigma = 0.25` and 16. While the flame is inside the chamber the vent placement is fine: the same 73 cm²/s
leaves either way, for a third of the cells. When the flame reaches the vent it is not: at :math:`\sigma = 0.25`
the vent flow stalls and reverses and the closed-end pressure spikes to 3.7% of ambient, 13× anything the extended
duct shows; at :math:`\sigma = 16` the chamber rings at its quarter-wave. With the boundary downstream the flame
crosses the vent plane as interior physics, and goes on to burn the charge it pushed out ahead of itself — an
external burn no boundary at the vent can represent.

The diagnostic is the two-:math:`\sigma` test: run at :math:`\sigma = 0.25` and 16 and compare. Here the vent-plane
flux differs by 67 cm²/s rms between the two with the boundary at the vent, and by 2.5 with it downstream. What a
longer domain does not settle is what the far boundary stands for: a soft one is an endless duct, a stiff one an
opening into a large space, and the closed-end pressure during the burn still differs between them. Choose
:math:`\sigma` for what really lies beyond the domain; remember that the relaxation rate
:math:`K = \sigma (1 - M^2) c / L` scales with the domain length.

Make the domain longer when something other than sound and vortices reaches the boundary during the measurement
and rung 6's settings cost too much, when the answer fails the two-:math:`\sigma` test, when what leaves keeps doing
physics outside, or when a face is in sustained backflow whose incoming state cannot be specified.
