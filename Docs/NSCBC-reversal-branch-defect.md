# NSCBC reversal-branch defect — investigation notes (2026-08-24)

**Status: FIXED and gated, same day — and the fix itself took two layers,
because the first layer's production rerun exposed a second, older defect.**

**Layer 1 — unify the acoustics.** The dedicated reversal closure is removed:
a reversal is counted (`Diag::reversed`) and falls through to the standard
outflow closure, whose acoustic side is well defined and restoring for
`u_out < 0` (`lambda_in = c − u_out` grows, `1 − M²` is even in M, and
`K (p − p_tgt)` both pulls the ghost pressure toward the target and pushes
the ghost normal velocity outward — the restoring feedback the
frozen-velocity branch lacked). `pin_farfield` now keeps operating during
reversal instead of being preempted, which is correct for a reservoir
boundary.

**Layer 2 — upwind the material family.** The unified closure alone NaN'd
BOTH vent variants at t ≈ 4.05 ms — reproducing the *original pin's* crash
signature, and the fields identified the mechanism the pin-era forensics
("a 138 K cell") had already photographed: with backflow, an EXTRAPOLATED
λ₀ ghost is advected back in and feeds on itself. The boundary column went
385 → 241 → 89 K in 42 µs while the inflow accelerated −52 → −1425 →
−11,773 cm/s: the outward-cooling T ramp refrigerates the boundary cell,
the steepened ramp extrapolates colder still, and `1/(ρc)` amplifies the
relaxation increment as T falls. The old soft branch survived the crossing
*because* it froze material content — its defect was only what it dropped
from the acoustics. So the final closure keeps the unified acoustics and
applies an upwinding factor `w_mat` to the material slopes (`dS`, `dY`,
`dut`, `dT`, and `extrap_material`'s `dRm`): 1 for `u_out ≥ 0` (every
forward-flow result is untouched to the bit), 0 under firm reversal (the
ghost keeps the interior's material VALUES), ramped across a narrow band
(10⁻³ c) below zero — a continuity regularisation of the fill, not a
physics knob. Treating *sustained* recirculation as a local inflow with a
target state remains queue item 2.

Gate: driver **C13** (`Verification/NSCBC1D/nscbc1d.cpp`,
`check_reversal_continuity()`), static, all four builds, three checks:
(a) ghost p, u AND T are continuous through `u_out = 0` on a swept
flame-like stencil (10 cm/s steps resolve the `w_mat` band); (b) an
over-pressured firm reversal still relaxes — ghost p toward target, ghost
`u_out` pushed outward; (c) a firm reversal does not extrapolate material —
with T falling toward the face, ghost T is the interior's own. Measured
red/green: the b764cfb branch fails (a) with a 4.220e3 dyn/cm² crossing
jump vs 1.180e-2 sweep variation (3.6e5×) and (b) with `du_out = +0.0000`
exactly; the layer-1-only kernel passes (a)/(b) but fails (c) with ghost
T = 341.5 K vs interior 400 K; the shipped kernel passes all three —
(a) crossing jump equal to the sweep's own variation, (b) `du_out = +1.49`
cm/s with `dp = −25.1`, (c) ghost T = interior T.

The sections below are the original investigation record, kept as written
(they are the reproduction protocol and the defect's evidence base); the
"suggested next steps" at the bottom are superseded by the fix above, and
the post-fix production A/B lives in the NSCBC-Chamber README.

---

This is the write-up of a Sonnet-5
session's investigation into a pressure/flow anomaly found in the first
post-fix `NSCBC-Chamber` production runs (commit `b764cfb`, "the reversal pin
becomes a soft closure"). Read this before touching `Source/NSCBC.H`'s
outflow-reversal branch or before regenerating the Lesson-9 comparison table
(work queue item 3 in `RESTART-PROMPT.md`); it also bears directly on queue
item 2 (the sustained-recirculation gate).

## Summary

The reversal branch added by `b764cfb` (`Source/NSCBC.H:604-638`) is
**physically incomplete relative to the forward outflow branch**: it drops
the reaction-source term (`S_p`, the `beta_s` mechanism), the transverse term
(`T_in`), and the material-slope continuation (`dRm`) that the forward branch
carries. This asymmetry across `u_out_N = 0` is invisible in every case the
suite has run to date because none of them puts a flame's *end of transit*
right at the boundary while the flow is dithering near zero. The
`NSCBC-Chamber` vent variant does exactly that, and production data shows the
closure switching back and forth across the branch boundary **with growing
backflow amplitude each crossing** — the signature of a feedback loop, not
passive weak damping — culminating in a bounded but large (~30% chamber
overpressure) spurious pressure/reversal event that does not appear when the
same physical location has no boundary condition on it (the plenum control).

This is not a sign error, not an EOS failure, not a NaN/corruption event —
the field stays finite and physical throughout. It is a **missing-physics
gap in one branch of a two-branch closure**, reachable specifically by the
flame-crossing recipe (`bc_nscbc_beta_s = 0`) this case is built to exercise.

## How to reproduce

```sh
cd Exec/RegTests/NSCBC-Chamber
mpiexec -n 6 ./PeleC2d.llvm.MPI.ex chamber-vent.inp \
  prob.pmf_datafile=$PWD/LiDryer_H2_p1_phi0_4000tu0300.dat
```
Runs to `stop_time = 8e-3` in ~1 CPU-hour (see case README). The event is at
t = 4.0-4.7 ms. The comparison data referenced below is the completed
2026-08-24 production set:
* `vent/` — boundary at the vent plane, `bc_nscbc_sigma = 0.25` (case default)
* `plenum/` — identical settings, boundary two chamber-lengths downstream
* `vent-sigma16/` — same as `vent/` but `pelec.bc_nscbc_sigma=16`

Extracted per-plotfile traces cached at `<rundir>-traces.csv` (columns:
`t, dp_closed, dp_mean, burnt_frac, Vdot_vent, sigma_exp`) by
`chamber_qoi.py`/`chamber_metrics.py` (needs `Verification/NSCBCFields/fielddump`
on `PATH` or via `FIELDDUMP=`).

## The event, as measured

Closed-end overpressure (`dp_closed`, dyn/cm² above ambient) and vent
volumetric flux (`Vdot_vent`, cm²/s, positive = outflow), same physical
location (x ≤ 1.2 cm) in all three runs:

| t [ms] | vent dp_closed | vent Vdot | σ16 dp_closed | σ16 Vdot | plenum dp_closed | plenum Vdot |
|---|---|---|---|---|---|---|
| 3.775 | 607 | 74.3 | — | — | 1871 | 74.0 |
| 4.075 | 16,166 | 9.4 | — | — | 1819 | 72.2 |
| 4.375 | **143,540** | −114.0 | — | — | ~1815 | ~30 |
| 4.525 | **321,236 (peak)** | — | — | — | — | — |
| 4.647 | — | — | **12,684 (peak)** | — | — | — |
| 4.675 | 32,869 | 188.1 | — | — | 1905 | 7.5 |
| 4.975 | 404 | −9.6 | — | — | 2051 | 1.2 |

σ = 16 (16× the relaxation rate) damps the peak 26× (321,236 → 12,684) but
does not eliminate it — consistent with the closure being *weak-but-present*
and scaling with σ as designed, not with a sign error. Only removing the
boundary from that location (plenum) eliminates the event, which is the
control that rules out "this is real chamber physics the boundary merely
under-suppresses" — the plenum's own trace through the same window is smooth
and monotone, no spike, no sign flip.

### Field-level: the branch is dithering, with growing amplitude

`x_velocity`, mid-row, last 6 interior columns approaching the vent, `vent/`
run (fielddump on `pressure`, `x_velocity`, `Temp`, `rho_H2`, `density`):

| t [ms] | u [cm/s], 6 cols nearest vent | branch |
|---|---|---|
| 3.950 | +152 → +245 | forward (S_p active) |
| 4.075 | **−185 → −102** | reversal (S_p dropped) |
| 4.150 | −159 → −47 (weakening) | reversal |
| 4.200 | −136 → −28 (weakening further) | reversal |
| 4.250 | **+239 → +381** | forward again |
| 4.300 | **−962 → −861** (5× stronger) | reversal again |

Three crossings in 0.25 ms, each reversal excursion larger than the last.
`Y_H2` in these same cells is 1e-4-1e-2 throughout — active local combustion
right where the branch is flipping, i.e. exactly where `S_p` (dropped only on
the reversal side) would matter.

### Not a corruption event

Domain extrema at six times spanning before/during/after (min/max over the
full domain): pressure, density, temperature, `Y_H2` all stay finite and in a
physical range throughout (p peaks at 1.31e6 dyn/cm² lokally ≈ 1.3 atm,
T stays in [297, 1662] K, `Y_H2` decays monotonically to ~0 by t = 5.0 ms and
stays there). No NaN, no negative density/species. The event is a bounded,
self-terminating (once the last charge finishes burning) spurious resonance,
not a blow-up.

*(Caveat: the domain-averaged `burnt_frac` metric in `chamber_qoi.py` dips to
0.61 around t = 5.0-6.0 ms after tracking 0.96-0.99 through the event — this
looked alarming at first but is a metric artifact, not physics: the script's
`y_fresh` normalizer is recomputed per-snapshot as `max(Y_H2)` in the chamber
region, which becomes noise-dominated once the true `Y_H2` field is ~1e-4 and
falling. The underlying `Y_H2` field itself is small and monotonically
decaying at those times — see the domain-extrema table above. Fix the metric
before trusting `burnt_frac` past ~t = 4.7 ms in either script; not otherwise
consequential to this investigation.)*

## Code: the asymmetry

`Source/NSCBC.H`, forward outflow branch (`u_out_N ≥ 0`, active whenever an
outflow is behaving normally), lines ~763-799:

```cpp
amrex::Real L_in;
if (outflow_face) {
  const amrex::Real K = prm.sigma * one_m_M2 * c_N / prm.L_ref;
  L_in = K * (qN.p - tgt.p);
} else { ... }
L_in -= (1.0 - beta) * T_in;                    // transverse term
if (prm.beta_s < 1.0) {
  bool src_ok = true;
  const amrex::Real S_p = reaction_dpdt(qN, src_ok);
  if (src_ok) {
    L_in += (1.0 - prm.beta_s) * S_p;             // reaction-source venting
  }
}
Rm_g = Rm_N + layer * dRm + layer * dx_normal * L_in / (lambda_in * rho_c);
```
Ghost density/velocity then come from a full R₊/R₋ characteristic
reconstruction (further down in `apply()`), and `dRm` (material-slope
continuation, when `extrap_material` is on) is folded in too.

Reversal branch (`u_out_N < 0` at an outflow face — this run's `beta_s = 0`,
`extrap_material = 0`, so `dRm` would be 0 there regardless, but `S_p` and
`T_in` are architecturally unavailable), lines 604-638:

```cpp
const amrex::Real K_rev = prm.sigma * one_m_M2_rev(mach) * c_N / prm.L_ref;
amrex::Real p_g_rev =
  qN.p - layer * dx_normal * K_rev * (qN.p - tgt.p) / (2.0 * (c_N - u_out_N));
// ... floor, EOS T from (qN.rho, p_g_rev, qN.Y) ...
pack_ghost(qN.rho, qN.u, e_g, T_g, qN.Y, s_N, s_ghost);   // rho, u FROZEN from interior
```

Algebraically, `p_g_rev`'s correction is *exactly* the pressure change the
forward branch's `Rm_g` increment would produce **if `S_p`, `T_in`, and `dRm`
were all zero** — confirmed by hand (both reduce to
`Δp = -layer·dx_normal·K·(qN.p - tgt.p) / (2·(c_N - u_out_N))` in that limit).
So the commit's claim that the closure is "the same increment form as the
forward branch" is true for the bare pressure-relaxation term, but the
forward branch is *never* run with `S_p` and `T_in` forced to zero in this
case (`beta_s = 0` and `beta = 0.5` are the recipe) — meaning the two
branches compute measurably different ghost pressures for the same interior
state at `u_out_N = 0⁺` vs `0⁻` whenever there is local reaction or transverse
flow, which there is, right here.

Ghost density and velocity are also handled differently in kind, not just
degree: the forward branch reconstructs them from the characteristic system;
the reversal branch freezes them to the interior cell's values outright. That
compounds the discontinuity but was not isolated as a separate contributor in
this pass.

## What this is not (ruled out)

* **Not the σ-vs-reflection trade-off working as intended.** σ = 16 damps but
  does not remove the event; the *qualitative* signature (multiple growing
  zero-crossings) is not what a merely-weak-but-monotone relaxation produces.
* **Not real chamber gas dynamics the boundary fails to fully suppress.** The
  plenum control has the identical local combustion event at the identical
  location and time with no boundary present, and shows a smooth trace, no
  reversal, no spike.
* **Not a numerical blow-up / state corruption.** All fields stay finite and
  physical throughout; the run completes to `stop_time` normally in all three
  variants.
* **Not exercised by any existing gate.** The driver's C6 (`nscbc1d.cpp`,
  `check_fallbacks()`) calls the reversal branch exactly once, at a single
  static state, `layer = 1`, and asserts only finiteness + counter increment
  — it does not check restorative direction (the way C2 does for the forward
  branch) and has no dynamic/sustained-reversal scenario. This is a real gap,
  not an oversight in this write-up: queue item 2 ("Phase-1 coverage:
  backflow branch as local inflow with a sustained-recirculation test") is
  precisely the missing gate that would have caught this.

## Suggested next steps (not yet decided)

1. Carry `S_p` (and ideally `T_in`, `dRm`) into the reversal branch so the
   closure is symmetric across `u_out_N = 0` — the most direct fix, but needs
   care: the reversal branch currently has no characteristic reconstruction
   for density/velocity, so folding these terms in means either building
   that reconstruction for the reversal case too, or finding a reduced form
   that's provably continuous.
2. Alternatively, blend the two closures over a small Mach window straddling
   zero (e.g. a smoothstep in `u_out_N/c_N`) instead of a hard branch switch,
   which would remove the discontinuity without requiring the reversal branch
   to replicate every forward-branch term exactly.
3. Either direction needs a **dynamic** driver check before it can be trusted
   — a case that drives `u_out_N` through zero repeatedly under an active
   `S_p` (or synthetic surrogate) and asserts the ghost pressure/state does
   not develop growing oscillation. This is the natural home for queue item
   2's sustained-recirculation gate; it should probably be built to catch
   *this* mechanism specifically, not just generic backflow.
4. Regenerating the Lesson-9 comparison table (work queue item 3) should wait
   on a decision here — the vent variant's peak-aligned QoIs as currently
   computed are dominated by this artifact, not by the T7 physics the case is
   meant to measure (see the conversation this doc summarizes, or just rerun
   `chamber_qoi.py vent plenum --sigma16 vent-sigma16 --xvent 1.2` and note
   the peak lands at t ≈ 4.5 ms with the anomalous magnitude above).

## Artifacts

* `Exec/RegTests/NSCBC-Chamber/chamber_qoi.py` — new script (T7 QoI
  extraction: peak-aligned overpressure, V̇ balance, ring-down); not yet
  updated to exclude or separately report this window.
* `Exec/RegTests/NSCBC-Chamber/{vent,plenum,vent-sigma16}-traces.csv` — cached
  per-plotfile traces backing the tables above.
* Production plotfiles in `vent/`, `plenum/`, `vent-sigma16/` (321 each,
  git-ignored, reproducible from the inputs listed above).
