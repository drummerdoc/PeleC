# NSCBC on the `nscbc` branch — independent status evaluation, development plan, and a teaching curriculum

*Written 2026-08-21 against `nscbc` @ `84aa1ac` (23 commits over upstream `development` @ `ac17bd0`), the
`RESTART-PROMPT-nscbc.md`, `Docs/NSCBC-design-and-literature-review.md`, and the six documents in `PAPERS/`.
This document is a fresh reading, not a continuation of the earlier session's notes: where it disagrees with the
design document or the commit narrative it says so explicitly and gives the derivation or the line number.*

---

## 0. Verdict in one paragraph

The branch is a serious, well-instrumented piece of work: a GPU-clean, stateless ghost-cell characteristic
boundary kernel (`Source/NSCBC.H`, 1229 lines) wired into `BCfill` with a per-point problem hook, a 1-D
verification driver with twelve checks across four builds, five PeleC regression families (acoustic pulse in
2-D/3-D/AMR/EB/radial, vortex and circular pulse, sitting flame, exiting flame, reacting PMF), bit-identical
determinism under decomposition, MPI and restart, and a design document that is honest about what the
ghost-cell (GC) form cannot do. It is ready for acoustic waves leaving through a subsonic outflow and for steady
or slowly varying subsonic inflows. It is **not** yet ready for the application it targets (vented explosions,
fronts crossing the boundary): the shipped default σ = 0.25 crashes on a flame transit, the working recipe
(σ = 16 + `extrap_temperature`) buys stability with 28 % acoustic reflection, and the vortex test shows the
transverse terms behaving opposite to the published literature. Two defects found in this review — an inverted
sign on the reaction source term in the incoming wave (`NSCBC.H:1029`) and a tangential stencil that clamps
instead of wrapping in periodic directions (`BCfill.cpp:115–121, 179–182`) — are cheap to fix, sit directly
under the two worst measurements (flame exit, vortex β-response), and should be fixed and re-measured **before**
any of the roadmap items in the restart prompt (T5/T3, T7, boundary registers) are started, because those
items inherit their baselines from the current tables. After that, the literature in `PAPERS/` points at one
architectural addition — per-point boundary registers enabling NDNR/NRI — which is the right medium-term goal,
and the restart prompt's ordering is otherwise sound.

---

## 1. What the branch contains, and what it has actually established

### 1.1 Inventory

| Area | Files | State |
|---|---|---|
| Kernel | `Source/NSCBC.H` | Complete for subsonic in/outflow, supersonic copy/Dirichlet branches, flow-reversal fallbacks, transverse (Yoo–Im) terms, reaction source, two λ₀ closures (entropy / temperature), material continuation, far-field pin, diagnostics counters |
| Plumbing | `Source/BCfill.cpp`, `PeleC.cpp/.H`, `Params/*` | Per-direction `Params`, `bcnormal_nscbc` problem hook (`ProblemSpecificFunctions.H`), parameter validation, EB/AMR guards, counter reporting via `sum_interval` |
| Hydro/diffusion hooks | `Hydro.cpp`, `Diffusion.cpp` | Dead legacy Fortran NSCBC plumbing removed; hook-point comments only (the boundary heat flux is controlled through the `extrap_temperature` ghost closure, not here) |
| 1-D driver | `Verification/NSCBC1D/nscbc1d.cpp` (2442 lines) | Checks C1–C12; air/LiDryer/ADV=2/SRK builds; CI job `NSCBC-Driver` |
| Field tools | `Verification/NSCBCFields/` | `fielddump` + `metrics.py` |
| PeleC cases | `Exec/RegTests/NSCBC-{Acoustic,COVO,FlameOutflow,FlameOutflow-DRM,PMF}` | Inputs + READMEs with measured tables; registered in `Tests/CMakeLists.txt` as gold-file regressions |
| Docs | `Docs/sphinx/BoundaryConditions.rst`, `Docs/NSCBC-design-and-literature-review.md` | User guidance and design rationale |

### 1.2 What is proven (I re-read the evidence and accept these)

The acoustic contract holds: a planar pulse leaves a σ = 0.25 outflow with R ≈ 0.8 % of incident amplitude in
both the driver and PeleC 2-D/3-D; the circular pulse stays round to 0.07 % amplitude spread around the arc at
β = 0.8 versus 21 % with a hard boundary; the relaxation rate τ = 1/K is mesh-independent (C5); the inflow
reflection curve R = 2.3/4.8/19/57 % at `relax_u` = 0.5/2/10/50 is measured (C4). The fill is a pure function of
valid interior data and is bit-identical under box decomposition, 1-vs-2 ranks and checkpoint/restart. The
ghost-pressure bias Δp ≈ ρcL(∂u/∂n)/σ is derived, reproduced exactly in the driver (C9a) and shown to be a minor
(~5 %) part of the sitting-flame error; `extrap_temperature` is the correct ghost closure for the diffusion
operator (C8, C12). The EB body-state detectability problem is real and the `eb_zero_body_state=1` requirement is
justified.

### 1.3 What is asserted but, in my reading, not established

Three conclusions in the commit narrative rest on measurements that the defects in §2 contaminate, and should
be demoted from "settled" to "to be re-measured":

* *"β_s is measured-ineffective"* and *"an amplitude-side source/viscous term double-counts"* — both were
  measured with the source entering `L_in` with the wrong sign (§2.1). The refutation of the diffusive term
  (C12: 104 → −911) is very likely the same sign error seen through a different term.
* *"β = local Mach is on the wrong side of the minimum; β = 0 is near-unstable"* (COVO table) — measured on a
  y-periodic domain whose periodic corner ghosts and first/last-row transverse stencils are wrong (§2.2). The
  literature (Yoo–Im 2007, Lodato 2008, Motheau 2017) consistently shows β = M exiting a low-Mach vortex
  cleanly; the branch's β_opt(θ) formula is a plane-wave argument applied to a vortex and should not be
  published as guidance until the vortex case is re-run on a fixed stencil.
* *"The ghost-cell form is Motheau's NSCBC-GC"* (`BoundaryConditions.rst:49–50`) — it is not (§2.8). This
  matters for how the upstream PR is framed, not for correctness.

---

## 2. Independent audit — defects and risks, ranked

Line numbers refer to `84aa1ac`.

### 2.1 High — the reaction source enters the incoming wave with the wrong sign (`NSCBC.H:1029`)

Derivation, in the kernel's own variables (outward normal, `u_out`, frozen ρc). The LODI pressure and
normal-velocity equations with transverse (𝒯) and chemical source (S_p = dp/dt|_chem, > 0 for heat release,
which is what `reaction_dpdt` returns — check its `heat_term` sign at lines 671–674) read

```
∂p/∂t = −½(L5 + L1) − 𝒯_p + S_p
∂u/∂t = −(1/2ρc)(L5 − L1) − 𝒯_u
```

The incoming invariant is R₋ = u − p/ρc, so

```
∂R₋/∂t = (1/ρc) [ L1 + (𝒯_p − ρc 𝒯_u) − S_p ] = (1/ρc) [ L1 + T_in − S_p ]
```

where `T_in` is exactly the quantity assembled at lines 953–958. A perfectly non-reflecting boundary wants
∂R₋/∂t = 0, i.e. L1 = −T_in + S_p; with relaxation and the β weights,

```
L1 = K(p − p_t) − (1−β) T_in + (1−β_s) S_p .
```

The code has `L_in -= (1.0 - beta) * T_in` (correct) and `L_in -= (1.0 - prm.beta_s) * S_p` (sign inverted).
With β_s = 0 the boundary therefore *doubles* the heat-release push on R₋ instead of removing it. The
Sutherland–Kennedy (JCP 2003) result — that the correct source term cancels the steady-flame pressure offset
(1−α)S_p/K exactly — is the known cure for precisely the "ghost-pressure bias at a sitting flame" the branch
spent several commits characterising; with the present sign the offset is 2S_p/K, worse than with the term
off, which is what the tables show (PMF 65 → 80; FlameOutflow 2065 → 2074). Check C7 gates only the
thermodynamics of `reaction_dpdt`, never the sign of its use, so no gate could have caught this.

*Is this a convention disagreement (Poinsot vs Motheau) rather than an error?* No, for three reasons.
(i) There is no convention switch: `Params::validate()` (lines 365–383) requires β_s ≤ 1 and the sign is
hard-coded, so the control interface cannot select it. (ii) The kernel works in its own outward frame with
its own invariants, and checks the sign of its relaxation term itself (lines 908–914: p_N > p_t ⇒ ghost
pressure falls). Whatever the L-definitions differ by between authors (lo/hi faces, T1/T4), the outward frame
absorbs it; the source term only has to be consistent with the kernel's *own* relaxation term, and it is not.
(iii) Convention-free check from the exact steady state: heat release at constant pressure produces the
dilatation ∂u/∂n = S_p/(ρc²) (pressure equation with ∂p/∂t = 0), so the exact R₋ gradient at a steady flame
is +S_p/(ρc²); the kernel's mapping gives dR₋/dn = L_in/(ρc²) for u ≪ c, so the fill reproduces the exact
state with p = p_t only if L_in = K(p − p_t) **+** S_p. A 1-D linear-acoustics model with the kernel's fill
transcribed line for line (`_to_delete/source_sign_check.py`: rigid wall, source confined to the boundary
cell, minmod R₊ slope, Godunov fluxes) gives equilibrium offsets (p_N − p_t)K/S_p = 0.000 / 1.000 / 2.000 for
+S_p / term off / −S_p as written.

Fix: one character. Then re-run C12 (with the removed diffusive term restored under the corrected sign — it
was built with the same convention), the FlameOutflow σ/β_s table, and `nscbc-flameexit.inp` at σ = 0.25 and
σ = 1 with β_s = 0 before anything else.

### 2.2 High — tangential stencil clamps into the domain instead of wrapping in periodic directions

`nscbc_fill` clamps the stencil base and the transverse neighbours into `[domlo, domhi]` (`BCfill.cpp:115–121`,
`179–182`). For a wall-type tangential face that is right; for a **periodic** tangential direction the corner
ghost `(i < domlo, j < domlo)` must be the image of row `domhi`, and the transverse stencil at `j = domlo` must
wrap. The legacy path (`BCfill.cpp:261–265`) reads the already-periodic-filled tangential ghost and gets this
right; the NSCBC path does not. Every 2-D regression case except the circular pulse is periodic in the direction tangential to its
NSCBC faces (Acoustic, COVO, FlameOutflow in y; PMF in x), so every baseline contains it; it is invisible to uniform-in-y controls. The CTU transverse predictor and the
diffusion operator's transverse gradients read those corner ghosts in the first and last rows. Fix: use
`geom.isPeriodic(d)` and, for periodic `d`, read `dest` at the un-clamped index (it is valid after
`FillBoundary`, which runs before the physical-BC functor in `FillPatchSingleLevel`). Add a one-line gate:
after a fill on a periodic case, a ghost row must equal its wrapped image to round-off.

### 2.3 High (usability) — the shipped default does not survive a flame transit

`nscbc-flameexit.inp` at σ = 0.25: front pushed backwards, wrinkle +40 %, NaN. The working recipe (σ = 16 +
`extrap_temperature`) is documented, but a default that crashes the flagship use case is a release blocker
regardless of documentation. Re-evaluate after §2.1; if the transit still needs σ ≈ 16, the kernel should
either switch regimes automatically on the transit guard (`Diag::structure` already fires at the crossing
columns) or the default should move and the docs should say why.

### 2.4 Medium — COVO β-response contradicts the literature

See §1.3. Re-measure after §2.2 with Motheau's exact vortex setup (M = 0.2 … their Table/Fig. conditions) and
β ∈ {1, M, 0.5, 0}. If β = M is still worse than β = 1, the next suspect is double application of transverse
physics: the ghost column already carries `T_in(j)` row by row *and* the CTU predictor applies transverse flux
differences across that column — in the flux form the transverse term enters the boundary ODE once.

### 2.5 Medium — backflow closure is discontinuous and reflective (`NSCBC.H:853–880`)

An outflow cell with `u_out < 0` snaps to a pressure Dirichlet at `p_t` with the full interior velocity,
interior ρ and Y, and T recomputed from (ρ_N, Y_N, p_t); for `u_out → 0⁺` the ghost is `p_N − O(Δx·K)`. A cell hovering at zero normal velocity (the tail
of a rarefaction, a recirculation touching the face — both routine at a vent) chatters between closures on
successive fills. The 3-D radial case logged 11 456 reversals in 40 steps. A proper local-inflow branch
(targets = ambient T, Y, tangential velocity 0, same K family) is the standard treatment and is what a vent
plane needs.

### 2.6 Medium — TurbInflow fluctuations are silently dropped on NSCBC inlets

`add_turb` writes fluctuations into the first ghost layer; the NSCBC path then overwrites `dest(iv)` with no
`turb_fluc` argument to the hook (`BCfill.cpp:200–206`, `279–283`). The sphinx table (`BoundaryConditions.rst:454`)
says the combination works. It does not; either wire it (pass the fluctuation into `Target::u`) or say so.

### 2.7 Medium — box-decomposition dependence with β < 1 is weaker than claimed but not zero

The transverse neighbours are also clamped into the FAB (`BCfill.cpp:179–182`); because `dest` is the grown
FAB and interior ghosts are valid after `FillBoundary`, this clamp rarely binds, but `n_stencil` depends on
`depth_fab` for boxes thinner than three cells, and any future `blocking_factor < 3` case will see it. Note in
the docs; guard in `read_params`.

### 2.8 Medium (framing) — this is not Motheau's NSCBC-GC

Motheau, Almgren & Bell (AIAA J. 2017) compute the wave amplitudes L_i at the boundary with one-sided
high-order normal derivatives of the primitive field, replace the incoming L_i by the relaxation (with β = M
transverse terms), invert the characteristic relations into **target normal derivatives** at the boundary, and
fill ghosts by a finite-difference inversion consistent with the interior stencil. This branch extrapolates
linearised **invariants** R±, S with a minmod-limited per-cell slope, injects the relaxation as an R₋ gradient,
and closes with the EOS. The relaxation rate and the transverse expression are faithful to Poinsot–Lele/Yoo–Im;
the spatial realisation is new. Call it what it is ("invariant-extrapolation ghost-cell NSCBC") in the sphinx docs and
the upstream PR (the `Hydro.cpp:118` comment already attributes Motheau correctly to the deleted Fortran), and note that the deleted Fortran's corner treatment (both incoming families modelled) is not
reproduced — corners are owned by the lowest live `idir` and the other direction's stencil is clamped.

### 2.9 Lower-severity items

Wall/NSCBC corners: the NSCBC fill overwrites the wall's reflected ghost so the ghost column is y-constant
across a wall face (`BCfill.cpp:64–77`); untested. Frozen ρc and c² across a 3-cell stencil that spans a flame
(ρ ratio 4–6) makes δR₊, δS differences of inconsistent variables — acknowledged, but it is exactly the
configuration the flame cases measure. Supersonic inflow silently uses `qN.p` when `Target.p` is unset
(`NSCBC.H:836`). Boundary pressure is rebuilt from `UTEMP` rather than `e` (`NSCBC.H:572–578`), an O(Δt)
inconsistency after time interpolation. Counters only print when `sum_interval > 0` (default −1), so the
"silent fallback" problem the counters exist to solve is still silent by default. Per-thread local memory for
a 53-species mechanism is ≈ 2× the header's estimate. The AMR `static bool warned[][]` is shared across levels.

---

## 3. What the papers in `PAPERS/` change

First, the file map, since the names mislead: `1-s2.0-S0021999126003402` = Dupuy, Meziat Ramirez, Douasbin,
Poinsot, *JCP* 562 (2026) 114987 (NDNR); `1-s2.0-S0045793019301951` = Daviller, Oztarlik, Poinsot, *Comput.
Fluids* 190 (2019) 503 (NRI-NSCBC); `FINAL QUILLATRE.pdf` = Vermorel, Quillâtre, Poinsot, *Combust. Flame* 183
(2017) 207 (SydGex multiscale); `Multiscale_NewVersion_for_review.pdf` = Quillâtre, Vermorel, Poinsot, Ricoux,
*Ind. Eng. Chem. Res.* 52 (2013) 11414 (vented deflagration LES); `QUILLATRE.pdf` = Quillâtre et al., MCS7
2011; the `.docx` is an ALCF press article about the SydGex billion-cell run and contains nothing on boundary
conditions. The design document identifies them correctly.

The design document's Parts II–IV are accurate on the main points (NDNR band σ ≥ 10, τ∞ = 3t_a, f_c ≈ f₀/10;
NRI's u⁻ subtraction and I ≡ 1; factor 2 vs 1 injection; the plenum practice; the 650 Hz post-peak mode;
the CLR flame-exit artifact). Corrections and additions that affect the plan:

**σ conventions differ by a factor of two between Daviller and Dupuy/this code.** Daviller writes
L₅ = ρc[2K(u − uᵗ)] (his Eq. 18), relaxing u at rate K; Poinsot–Lele, Dupuy and `NSCBC.H:1017–1022` put K
(not 2K) in the wave, relaxing at rate K/2. Hence σ_Daviller = σ_repo/2: his σ = 5 (I > 20 at resonance),
17/170 (nozzle), 23 (turbulent channel), 15 (slot-flame DNS) correspond to repo σ ≈ 10, 34/340, 46, 30.
Dupuy's values transfer 1:1. The design document's line "literature σ values transfer directly" is true only
for the Dupuy/Poinsot–Lele convention; T3/T4/T6 targets must be doubled.

**NDNR is implementable as a small state machine, and the paper gives everything needed.** Per boundary point:
p⁺ = −½∫L⁺dt at an outlet (u⁻ = +(1/2ρc)∫L⁻dt at an inlet; mind the sign asymmetry in their Eq. 8); EMA
recursion Ṽ(t) = (1 − Δt/τ)Ṽ(t−Δt) + (Δt/τ)V(t) with Ṽ(t₀) = V(t₀); the **adaptive** start-up rule τ(t) =
t − t₀ until τ∞ (their Eq. 14 — without it the filter is biased for t < τ∞); outlet L⁻ = −K[pᵗ − (p − p⁺ + p̃⁺)]
(their A.1). Defaults from their Table 2: σ ∈ [10, L/(cΔt)], τ∞ = 3t_a with t_a = 2L/c, equivalently
f_{c,K} ≥ 3f₀ and f_{c,τ∞} ≈ f₀/10 with f₀ = c/4L; overestimating t_a is nearly harmless, underestimating it
is not. Their convergence criterion (|u − uᵗ| within 0.1 % of the perturbation for the rest of the run) and
viability threshold (t_conv < 30t_a) are the gates T1/T2 should use. Their DDE-validation tube (ρ = 1.17
kg/m³, c = 348 m/s, L = 0.435 m) and their wave-packet (Eq. 28: B = 2ρc, ε = 0.4, centre frequency f₀,
bandwidth f₀) are ready-made 1-D references.

**Daviller footnote 6** — "in practice the u⁻ integral must be high-pass filtered to remove any continuous
component" — is the bridge from NRI to NDNR and the main risk for any register implementation: an unfiltered
∫L dt drifts. Build NDNR, not NRI.

**The GC form has no amplitude slot for injection.** Daviller's injection terms (−2ρc ∂u_a/∂t acoustic,
−ρc ∂u_v/∂t vortical) enter L₅ directly; our inlet imposes a time-varying *value* through `relax_u`. The
transfer function of that value-inlet to an injected harmonic is unmeasured (T3/T5 in the design document).
Until it is, forced-response and thermoacoustic users should be warned that injected amplitudes are wrong by
a frequency-dependent factor of order 1/(1 + (ω/K)²)^{1/2}.

**The explosion papers validate on overpressure magnitude and trend only**; time-to-peak is explicitly
shifted out (C&F 2017 §2.2). The T7 QoI should therefore be the overpressure trace with peaks aligned, plus
the IECR V̇_comb − V̇_vent diagnostic (their Eqs. 3–5), whose sign change marks the peak and whose post-peak
oscillation is attributed to acoustic reflection modulating the vent flux — i.e. the boundary is *in* the
observable. The plenum-border BC in IECR/C&F is Granet et al. 2010 (3-D NSCBC with transverse terms), and the
chamber mesh coarsens in the last third of the *chamber*, not the plenum.

---

## 4. Development plan

The plan is in five phases. Each has an exit gate that is a number, not a feeling. Costs are for the 2-core
cloud sandbox documented in the restart prompt; your workstation will be faster.

### Phase 0 — Fix, re-baseline, re-decide (1–2 sessions)

1. Flip the sign at `NSCBC.H:1029`; restore the diffusive dp/dt term behind a `beta_d` with the corrected
   sign (it was built and verified exact; only its use was wrong). Update the `Params` comments.
2. Wrap periodic tangential directions in `nscbc_fill`; add the image-row gate.
3. Re-run, in this order: driver C1–C12 (expect C12's 104 → −911 to become 104 → ≲ 104); FlameOutflow
   σ/β_s/β_d table at σ ∈ {0.25, 1, 16}; `nscbc-flameexit.inp` at σ ∈ {0.25, 1, 16} × β_s ∈ {1, 0};
   COVO β ∈ {1, M, 0.5, 0} and the circular pulse. Refresh the README tables and the kernel's β guidance.
4. Decide the default σ and the transit behaviour from the new tables (§2.3). If the corrected source term
   lets σ ≈ 1 survive a transit, the NDNR motivation weakens at outlets and strengthens at inlets; if not,
   NDNR at outlets stays the top architectural item.
5. Rename the scheme in the docs (§2.8); document the Daviller σ factor; wire or un-document TurbInflow.

Gate: all regression gold files regenerated with a commit message that states which numbers moved and why.

### Phase 1 — Close the cheap coverage holes (1 session)

Backflow branch as local inflow (§2.5) with a sustained-recirculation test (a 2-D backward-facing step whose
recirculation bubble touches the outflow, or the 3-D radial case with a metric). Wall/NSCBC corner test
(2-D channel, no-slip walls, NSCBC outflow; gate: wall shear in the last column equals the interior column to
second order). Supersonic inflow `Target.p` handling. Counter reporting on by default when `bc_nscbc = 1`
(print at `plot_int`, not only at `sum_interval`).

### Phase 2 — Injection fidelity: T3 and T5 (1–2 sessions; driver first, PeleC 1-D/2-D after)

Harmonic inlet forcing against a hard-pressure outlet. T3 measures the deterioration index I(f) =
|L̂₅/L̂₅ᵗ| = 1/(1 + R₁e^{2ikL}), R₁ = K/(K − 2iω) in the repo convention, at three frequencies straddling the
first quarter-wave resonance for σ_repo ∈ {0, 4, 10}. T5 measures P_RMS(x) against Daviller's Eq. 33 with
k± = ω/(c ± u); analytic anti-node amplitude √2·ρc·u_a. These two gates are what the boundary-registers work
(Phase 3) must improve, so they come first.

### Phase 3 — Boundary registers → NDNR outlets, then NRI/NDNR inlets (3–4 sessions)

Design note first (half a session), then code. One `FabArray`-like register per characteristic face on each
level holding, per boundary point, `∫L_out dt`, its EMA, and `t₀`; updated exactly once per level advance at
the new-time fill (so SDC re-fills stay idempotent), regridded by copy-on-intersection with a fresh-start rule
(`Ṽ = V` where no parent data exists), written into and read from checkpoints, and passed to the fill as
additional inputs so `apply()` stays algebraic. The outgoing amplitude L_out is computed from the same minmod
slope the fill already forms (for an outlet L⁺ = (u + c)ρc·δR₊/Δx in the kernel's variables). Registers-off
is bit-identical CLR — that is the first gate. Second gate: T1/T2 reproduce Dupuy Figs. 8–9 (CLR optimum
σ ≈ 0.3 for a step; NDNR t_conv < 10t_a across σ ∈ [10, 100] for the wave packet). Third gate: `nscbc-flameexit`
at σ = 16 with τ∞ = 3t_a shows the post-exit residual unchanged and the peak transit error's acoustic component
(the 28 %) gone. Then the inlet side, and T3/T5 repeated: I ≡ 1 is the target.

### Phase 4 — Application capstone: T7 mini-SydGex, then T8 (2–3 sessions)

2-D Sydney-like chamber (25 × 5 cm scaled or as is), closed end, 1 cm burnt kernel, optional baffle, open
end; run (a) with a plenum and the boundary on the plenum's far side and (b) with the boundary at the vent,
Phase-0 recipe and Phase-3 NDNR. QoI: peak-aligned overpressure traces, V̇_comb − V̇_vent, and the
post-peak ringing frequency and decay. T8 (outwardly propagating laminar flame on a quarter circle, Markstein
slope independent of arc radius) is the boundary-independence gate for slowly arriving spherical dilatation.

### Phase 5 — Hardware and upstream

CUDA/HIP compile and CPU-vs-GPU bit comparison on NSCBC-Acoustic; register-spill profile at 53 species.
Curate the narrative into an AMReX-Combustion PR: kernel + BCfill + docs + Acoustic/COVO/PMF cases first;
the flame-exit case and registers as follow-ups.

### What I would *not* do

Reformulate to the flux form (the branch's own measurement that the GC bias is ~5 % stands, and the sign fix
removes the part of the residual that motivated it). Tune σ by formula in reacting cases (Dupuy's point: when
c evolves, t_a is unknowable — build the broad-band method instead). Publish the β_opt(θ) curve before §2.2 is
fixed.

---

## 5. Teaching curriculum — how the scheme handles subsonic inflows and outflows

Each lesson names the physics, gives a runnable recipe on cases that already exist in the branch (plus two
small new inputs to add), states what you should see and what the quantitative gate is, and names the knob
whose effect the lesson isolates. Do the lessons in order: every lesson's "expected" column assumes you have
seen the previous one. All runs are a few minutes on a workstation except where noted. Commands assume a
GNUmake build in each case directory (`make -j TPL` once for SUNDIALS, then `make -j`).

### Lesson 0 — The picture to hold in your head

At a subsonic boundary, of the five characteristic families (two acoustic u ± c, one entropy, two vorticity,
all travelling at u), the ones leaving the domain are *measured* from the interior and the ones entering
must be *modelled*. At an outflow one family enters (the u − c acoustic), so one quantity — pressure — must be
told what the outside world is. At an inflow four enter, so velocity (normal and tangential), temperature
and composition must be told. "Non-reflecting" means: model the entering family as if nothing were coming
back. "Relaxation" means: model it as a weak spring pulling the boundary toward a target, strength K = σc/L,
so that the mean does not drift. Every trade-off in the subject is the tension between those two sentences,
and σ is the dial. Rudy–Strikwerda's σ ≈ 0.27 minimises the reflection of a single pulse; Dupuy's DDE
analysis shows the same value minimises the step-response convergence time — and that both degrade quickly
on either side.

The ghost-cell twist in this branch: instead of integrating an ODE at the boundary node, the kernel paints,
in every ghost layer, the state the outgoing families extrapolate to plus the state the modelled incoming
family implies, and lets PeleC's Riemann solver and diffusion operator read those ghosts like any other
cells. The relaxation becomes a *spatial gradient* of the incoming invariant R₋ = u − p/ρc across the ghost
layers. Consequence you will see repeatedly: anything that makes the outgoing extrapolation wrong (a flame,
a strong dilatation, a vortex core) leaks into the pressure the boundary thinks it sees.

Read `Docs/NSCBC-design-and-literature-review.md` Part I alongside this lesson.

### Lesson 1 — A pressure pulse leaving through an outflow (1-D driver, then PeleC)

*Physics.* A Gaussian pressure pulse splits into two acoustic waves; the one reaching the outflow should leave
without trace. What comes back is the reflection coefficient R.

*Recipe (driver).* `cd Verification/NSCBC1D && make -j && ./nscbc1d` runs C1–C12; C4 is this lesson.
`./nscbc1d sweep` already prints the reflection coefficient over σ ∈ [0, 16] — that table is the lesson.
(The driver takes no other command-line parameters; Lessons 3 and 7 will need a ParmParse-style override
for σ, `relax_u`, L and the forcing — add it once, at the `Params` setup in `nscbc1d.cpp:~2380`.)

*Recipe (PeleC).* `Exec/RegTests/NSCBC-Acoustic`, `nscbc-acoustic.inp`; then the same with `pelec.bc_nscbc=0`
(hard boundary) and with `pelec.bc_nscbc_sigma = 0 / 1 / 4 / 16`. Measure R as the peak of |p − p_amb| in the
upstream half after the incident pulse has left, divided by the incident amplitude. The Acoustic README gives
the protocol; `Verification/NSCBCFields/metrics.py` has `residual`, `circularity` and `sphericity` modes but
no reflection mode yet — a ten-line addition (`reflection`: max |p − p_amb| over x < x₀ at the final plotfile
divided by `prob.amp`) that Lessons 1, 3 and 4 all use.

*Expect.* Hard: R ≈ 1. σ = 0.25: R ≈ 0.8 %. σ = 0 slightly lower but watch the mean pressure — it is now
unanchored and will drift in Lesson 3. σ = 16: R ≈ 28 %. Plot R(σ): it is monotone and it is the price list
for everything that follows.

*Knob isolated.* σ. Also run `bc_nscbc_order = 1` once: the reflection barely changes here (a pulse is smooth
and acoustic), which is the control for Lesson 5 where order matters a great deal.

### Lesson 2 — The same pulse at oblique incidence: transverse terms and corners (PeleC 2-D)

*Physics.* A circular pulse meets the faces at every angle. The 1-D model of the incoming wave is wrong by the
transverse terms 𝒯 = u_t∂_t p − ρc u_t∂_t u_n + γp∂_t u_t (Yoo–Im); β weights how much of them the incoming
model cancels. At the corners two faces must agree.

*Recipe.* `NSCBC-COVO`, `nscbc-pulse2d.inp`; β ∈ {1, 0.8, 0.5, 0}; and `pelec.bc_nscbc=0`. Plot pressure
contours at the time the front crosses the faces; compute the amplitude spread around the corner-directed arc
(`metrics.py`).

*Expect (current branch).* Hard boundary: the ring squares off and grows four reflection lobes, 21 % spread.
β = 1: 0.9 %. β = 0.8: 0.07 %. β = 0: 3 %. After the Phase-0 periodic fix nothing should change here (this case
is not periodic) — which makes it the control for Lesson 4.

*Knob isolated.* β. Note that β here is a *weight on a correction*, so "more" is not monotonically better.

### Lesson 3 — Holding the mean: drift, step response, and why σ has a floor (driver; new T1)

*Physics.* With σ → 0 the outflow is transparent and has no opinion about the mean pressure: any net mass or
heat input makes p drift without bound. With σ large the boundary answers every perturbation with a
reflection. Dupuy's step test shows the two failure modes on one axis.

*Recipe.* Add a driver mode (T1 of the design document): duct L = 0.435 m, air at 1 atm (c = 348 m/s),
hard-pressure (fully reflecting) far end, NSCBC inflow with a step in uᵗ of 1 m/s at t₀; σ ∈
{0.1, 0.3, 1, 3, 10, 30}. Record t_conv = last time |u − uᵗ| > 0.1 % of the step, in units of t_a = 2L/c.

*Expect.* A minimum near σ ≈ 0.3 and growth on both sides; at σ = 10 the trace is a slowly decaying
square-ish ringing at the quarter-wave period 4L/c, each bounce losing the fraction Lesson 1 told you. This
is the baseline NDNR (Phase 3) must beat: with an EMA at τ∞ = 3t_a Dupuy's Fig. 9 shows t_conv < 10t_a for
every σ in [10, 100].

*Knob isolated.* σ again, now on the *mean-holding* axis rather than the reflection axis — the same dial seen
from the other side.

### Lesson 4 — A vortex leaving (PeleC 2-D COVO)

*Physics.* An isentropic vortex carries vorticity and an entropy/pressure deficit but no acoustic content. An
ideal outflow lets it pass in silence; any pressure left in the domain afterwards was manufactured by the
boundary. The vortex also exercises the *inflow* (it is fed by one) and the transverse terms in earnest,
because at the core ∂_t u_n is large.

*Recipe.* `NSCBC-COVO`, `nscbc-covo.inp` (M = 0.2, strength 0.5, y-periodic); β ∈ {1, −1 (= local M), 0.5, 0};
the periodic no-boundary reference (`geometry.is_periodic = 1 1`, `lo_bc/hi_bc = Interior`); and the hard
reference. QoI: rms(p − p_amb)/δp_vortex over the domain after the vortex has left, versus the periodic floor.

*Expect (current branch).* Floor 1×; hard 14×; β = 1 10×; β = 0.5 3.8×; β = M 16×; β = 0 132×. **This table is
the one the Phase-0 periodic fix is expected to change.** Re-run it after the fix; if β = M moves to the good
side (as in Yoo–Im, Lodato, Motheau) the lesson becomes "β = M is the physical choice"; if it does not, it
stays "the transverse terms are being applied twice somewhere" and Phase 0 item 4 in §4 is the next step.
Either way, watch the movie of vorticity: with a hard outflow the core *stretches* along the boundary and a
pressure wave runs back upstream; with a good characteristic outflow the core crosses the plane undeformed.

*Knob isolated.* β, now on a vortical structure; and `relax_u` at the inlet — repeat the best β with
`bc_nscbc_relax_u = 0.5` and `10` and watch the upstream-running wave reflect off the inlet (R = 2.3 % vs 19 %).

### Lesson 5 — A flame *sitting* on the outflow: dilatation, the ghost-pressure bias and the closures

*Physics.* A flame is a sustained dilatation ∂u/∂n > 0 plus a 4–6× density drop plus heat release. The
outgoing invariant R₊ = u + p/ρc then has a slope dominated by dilatation, not acoustics; extrapolating it
manufactures a ghost pressure Δp ≈ ρcL(∂u/∂n)/σ. The λ₀ family (entropy/temperature/species) must also be
closed, and the *diffusion operator reads the ghosts*, so the temperature closure is the boundary heat flux.

*Recipe.* `NSCBC-FlameOutflow`, `nscbc-flameoutflow.inp` (H₂/air LiDryer, wrinkled sheet parked on the
outflow; shielded 480×64 reference for `measure.py`; absolute `prob.pmf_datafile` path). Run the matrix
σ ∈ {0.25, 1, 16} × `extrap_temperature` ∈ {0, 1} × β_s ∈ {1, 0}, first on the current code, then after the
Phase-0 sign fix. QoI: mean-pressure error vs the shielded reference at t = 2.4e-5 s, and the radical
profiles (`radicals.py`).

*Expect (current branch).* σ = 1: +2074 dyn/cm²; +1894 with `extrap_T`; σ = 16 + `extrap_T`: +67; hard: −527;
β_s = 0 moves 2065 → 2074 (the wrong-sign term, tiny and wrong way). **After the fix** β_s = 0 should reduce
the σ = 1 error by the Sutherland–Kennedy offset (1−α)S_p/K — expect the bulk of the +1200 "open physics
lead" in the restart prompt to be this term. `order = 1` here flips the sign of the error: the limiter
clamps the front's structure — this is the lesson that the hydro reconstruction is part of the boundary
condition in a GC scheme.

*Knobs isolated.* `extrap_temperature` (the diffusive closure), β_s (the reactive source), σ (anchoring vs
bias), `order`.

### Lesson 6 — A flame *leaving* (PeleC 2-D flame exit)

*Physics.* A front crossing the plane at U = 40 S_L is a travelling contact with large Δρ and modest Δu, at
M ≈ 0.01, carrying mass and enthalpy and essentially no acoustics. The boundary is asked to *anchor*, not to
be transparent: the error is set by how strongly the ghost pressure is held while the density jump crosses.

*Recipe.* `nscbc-flameexit.inp` as shipped (σ = 0.25, β = 0.5, β_s = 0, extrap_T and extrap_material on);
then σ = 0.25 with both extrap flags off (the bare default); then σ = 16 + `extrap_temperature` only (the
recipe); then hard `pelec.bc_nscbc=0`; then, after Phase 0, σ = 1 + β_s = 0.
`exit_metrics.py` tracks the front position (must advect at U − S_L), wrinkle amplitude (must stay frozen) and
the post-exit residual against the exact uniform fresh stream. The transit guard's `material structure`
counter should fire at the crossing columns — set `pelec.sum_interval = 50` to see it.

*Expect (current branch).* Bare σ = 0.25: front pushed *backwards*, wrinkle +40 %, NaN at the crossing.
Shipped file (σ = 0.25 with both extrap flags): −340 000 dyn/cm², flame quenched and expelled — the
material continuation turns the crossing into runaway venting. σ = 16 + extrap_T: peak transit error +1800,
post-exit +8.5 (0.18 % of ambient). Hard: −200 / +1.2 (the worst
acoustic boundary in the suite is near-perfect for a pure transit). `extrap_material` on during a transit:
runaway venting at low σ (it continues a structure that should leave). The lesson is that a ghost-cell
boundary has two regimes — acoustic (be transparent) and material (be an anchor) — and a front crossing
switches regime in a few cells; Dupuy's Appendix B shows the mirror image (CLR acoustics pulling a flame *out*
3× early). After Phase 0, re-run σ ∈ {0.25, 1} with β_s = 0: the question is whether the corrected source
term lets a small σ survive a transit; the answer decides how much of Phase 3 is needed at outlets.

### Lesson 7 — The inflow as a signal source: forcing, reflection and the factor of two (driver; new T3/T5)

*Physics.* An inlet that relaxes a *value* uᵗ(t) is a low-pass filter on injected signals with cut-off ~K and
a reflection coefficient R₁ = K/(K − 2iω) for the upstream-running wave; in a duct with a reflecting exit the
two combine into the deterioration index I(f) = 1/(1 + R₁e^{2ikL}), which peaks at the quarter-wave
resonances. Daviller shows I > 20 at σ_Dav = 5 (σ_repo = 10); the fix (subtract the outgoing wave u⁻ from the
deviation) needs the register state of Phase 3.

*Recipe.* Driver: duct, NSCBC inflow with uᵗ = ū + u_a sin 2πft, hard-pressure exit; f ∈ {0.7, 1.0, 1.3}×f₀
with f₀ = c/4L; `relax_u` ∈ {0.5, 2, 10}. Record the amplitude of the L₅-like quantity at the inlet (δR₊ in
the kernel's variables) relative to the injected 2ρc u_a, and the P_RMS(x) profile against Daviller Eq. 33.

*Expect.* I → 1 only as `relax_u` → 0 (and then Lesson 3's drift returns); at `relax_u` = 10 the anti-nodes are
over-predicted by ~50 % and the nodes filled in, as in Daviller Fig. 15. Then — and this is the demonstration
that motivates Phase 3 — the same plot with the NRI/NDNR inlet is flat at I = 1 for every `relax_u`.

*Knob isolated.* `relax_u`, now as a transfer function rather than a single R.

### Lesson 8 — Backflow and walls at the boundary (PeleC 2-D; Phase 1 cases)

*Physics.* Outflow faces that momentarily ingest fluid (rarefaction tails, recirculation, vent-plane
eddies) need a local-inflow model; corners with no-slip walls need the wall's reflection preserved. Both are
ordinary at a vent.

*Recipe.* The 3-D radial acoustic case (`nscbc-acoustic-3d-radial.inp`) already produces thousands of
reversal events — print the counters. After Phase 1, the backward-facing-step case with the recirculation
touching the outflow: QoI = pressure noise generated at the face versus an extended-domain reference.

*Expect (current branch).* Each reversal snaps the ghost to p_t; the residual improvement over hard drops from
120× to 1.5× in the radial case. After Phase 1 it should be continuous in u_n and silent.

### Lesson 9 (capstone) — The vented chamber two ways (PeleC 2-D; Phase 4)

Mini-SydGex with the boundary on a plenum versus at the vent. You will see, in one overpressure trace, every
earlier lesson: the laminar-phase dilatation (Lesson 5), the front crossing the vent plane (Lesson 6), the
post-peak ringing at the first longitudinal mode whose decay is set by the boundary's reflection (Lessons 1
and 3), and — if Phase 3 is in — the NDNR boundary holding the mean without ringing. The QoI is
V̇_comb − V̇_vent (IECR Eqs. 3–5): its zero crossing is the peak, and the boundary's impedance is visible in
what happens after it.

---

## 6. Parameter cheat-sheet (current code)

| Key | Default | What it does | Lesson |
|---|---|---|---|
| `pelec.bc_nscbc` | 0 | Master switch; acts only where `bcnormal_nscbc` returns a live type | all |
| `bc_nscbc_sigma` | 0.25 | Outflow pressure relaxation, K = σ(1−M²)c/L | 1, 3, 5, 6 |
| `bc_nscbc_relax_u` | 2.0 | Inflow normal-velocity relaxation (same K form ×ρc) | 4, 7 |
| `bc_nscbc_relax_t` | 0.2 | Inflow T and tangential-velocity relaxation (per-layer blend) | 4 |
| `bc_nscbc_beta` | 1.0 | Transverse weight: 1 off, 0 full, <0 local Mach | 2, 4 |
| `bc_nscbc_beta_s` | 1.0 | Reaction-source weight: 1 off, 0 full (**sign to be fixed**) | 5, 6 |
| `bc_nscbc_order` | 2 | 1 copy, 2 minmod-linear extrapolation | 1, 5 |
| `bc_nscbc_extrap_temperature` | 0 | Close λ₀ on T (correct diffusive flux) | 5, 6 |
| `bc_nscbc_extrap_material` | 0 | Entropy-bounded material slope on R₋ (fronts that sit, not leave) | 5, 6 |
| `bc_nscbc_pin_farfield` | 0 | Hard R₋ pin to (u = 0, p_t); ignores σ, β | 6 |
| `pelec.sum_interval` | −1 | Set > 0 to see the fallback/transit counters | 6, 8 |

---

## 7. Immediate next actions (for the next working session)

1. Apply the two Phase-0 fixes in the cloud clone, rebuild the driver (~2 min after the AMReX 1-D install),
   run C1–C12, and rebuild one PeleC exe (FlameOutflow; ~40 min on 2 cores).
2. Re-run the FlameOutflow and flame-exit matrices with β_s = 0; re-run COVO β ∈ {1, M, 0.5, 0}.
3. Update the three READMEs and the kernel's `Params` comments with the new tables; correct the Motheau
   attribution and add the Daviller σ-factor note to `BoundaryConditions.rst`.
4. Only then pick up the restart prompt's queue at T5/T3, with the driver σ-from-command-line and the
   T1/T3/T5 modes added together, since they share the duct setup.
