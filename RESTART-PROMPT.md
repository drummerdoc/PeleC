# Restart prompt — PeleC NSCBC work

Paste everything below into the first message of a new session.

---

We are continuing NSCBC development in my PeleC repo (`/Users/marcusd/src/PeleC`, branch
`nscbc`, based on upstream `development` @ `ac17bd0`). Read, in order: `CLAUDE.md` (the working
brief; its "Decisions settled by measurement" and "State" sections are authoritative),
`Docs/NSCBC-design-and-literature-review.md` Part I (the formulation), and the READMEs of
whatever you touch (`Exec/RegTests/NSCBC-*`, `Verification/NSCBC1D`). The git log carries the
full reasoning — the commit messages are the lab notebook.

## Environment (this machine)

* Build against `/Users/marcusd/src/PeleLMeX/Submodules/PelePhysics` by passing
  `PELE_PHYSICS_HOME=` on every make line (PeleC's own submodule is uninitialised; SUNDIALS is
  prebuilt under that tree's `ThirdParty/INSTALL`, no `make TPL` needed). `COMP=llvm`.
* The 1-D driver builds in seconds; a 2-D case exe in ~5 min; the flame matrices run in
  minutes to ~2 h. Run sims in the background and poll — ONE production job at a time,
  at 6 MPI ranks (Marc's standing rule: overlapped 8-rank runs browned the laptop out);
  throttle builds to `make -j4` while a run is live. `make realclean COMP=llvm` before
  rebuilding the driver after a `Chemistry_Model` switch OR a kernel edit — make will not
  notice either on its own (the exe target has no prerequisites), and realclean without
  `COMP=` cleans the wrong suffix and silently does nothing.
* `prob.pmf_datafile` is CWD-relative: always pass an absolute path.

## State (all measured; do not relitigate without new data)

* Reaction source enters the incoming wave with `+(1-beta_s)*S_p` (Sutherland–Kennedy); gate
  `Verification/NSCBC1D/source_sign_check.py`.
* Flame-on-boundary recipe: `bc_nscbc_extrap_temperature=1` + `bc_nscbc_beta_s=0` — holds
  hard-outflow accuracy at any sigma sitting, survives transits at any sigma. `sigma=16` +
  `beta_s=1` is best absolute anchoring at 28% reflection; do not stack sigma=16 with beta_s=0.
* Periodic tangential stencils CLAMP into the domain. Wrap (aperiodic, 2e-2) and strip-band
  (bitwise but ~100 lines of machinery for a 2e-4..1.6e-3 residual) were built and rejected;
  the seam gate reports the residual and aborts above 1e-2.
* beta: 1 default (safe), 0.5 optimum for broad structures leaving, beta=0 near-unstable,
  pointwise local-Mach worse than off. beta_opt(theta) = 1 − cos/(1+cos).
* `extrap_material` is for fronts that SIT, never during a transit.
* Determinism (decomposition, MPI, restart) is bit-identical and must stay so; the fill is a
  pure function of the interior state — no boundary history.
* Outflow reversal runs through the SAME closure as forward outflow (unified acoustics,
  restoring in p and u, `S_p`/`T_in` live) with the λ₀ material slopes upwinded off
  (`w_mat`) — no dedicated branch. Three dedicated branches died before this: the ambient
  pin (NaN via pressure step), the soft pressure-only relaxation (spurious 0.32 atm chamber
  spike via dropped terms + frozen ghost velocity), and unified-without-upwinding (cold
  runaway NaN via self-fed material extrapolation). Gate: driver C13 a/b/c; forensics in
  `Docs/NSCBC-reversal-branch-defect.md`. Forward-flow histories are bit-identical (recipe
  rows reproduce to the digit); failure-mode rows re-measured +14–22% honest. Treating
  SUSTAINED recirculation as local inflow with a target state is still queue item 2, but
  note the reversal counters are silent unless `pelec.sum_interval` is set.
* C11x (driver, reported not gated): profile-fit ghost closure (fitU + source-consistency
  bound) reaches ~oracle level on a sustained front, releases unsustainable structure, robust
  to 15% end-state and full shape distortion of the family. 2-3D candidate; open items in
  `Verification/NSCBC1D/README.md`.
* Driver: air 56/56, LiDryer 60/60, SRK 43/43 (statics; dynamics skipped by design);
  PELE_NUM_ADV=2 expected 56 (CI). CI job `NSCBC-Driver` runs all four; GNUmake takes
  `PELE_EOS=SRK`. C13 (reversal continuity/relaxation/material-freeze) is static and runs
  in all builds. The driver GNUmake needs `make realclean COMP=llvm` before rebuilds —
  it will NOT notice kernel edits or mechanism switches on its own.
* Static GPU audit of the branch code is clean (commit e6a8ad9): device paths, captures,
  counters, and every synchronize reviewed. The actual device build remains item 5.

## Local artifacts (this machine, not in git)

* `Exec/RegTests/NSCBC-{FlameOutflow,COVO}/runs/` hold the animation runs and MP4s
  (flame exit 3-way, vortex, pulse) — git-ignored, reproducible from the inputs.
* The teaching page "Waves at the Open Boundary" (σ/β tutorial with live figures) is a
  claude.ai artifact: https://claude.ai/code/artifact/4d4f11af-7932-4c80-95f7-cf9bd90b984a
* Untracked `PAPERS/` at the repo root is Marc's reference library — expected, leave it.

## Resolved 2026-08-25: the reversal closure is now unified and gated

The chamber production A/B found the `b764cfb` reversal branch dropping
`S_p`/`dR₊`/`T_in` (spurious 0.32 atm spike at the flame-vent crossing), and
the first fix attempt exposed a second defect (extrapolated material content
advected back in during backflow — the cold runaway, 385 → 89 K, NaN; the
original pin's own failure mode). Final closure in `Source/NSCBC.H`: ONE
outflow closure both sides of `u_out = 0` (unified acoustics — restoring in
both p and u, `S_p`/`T_in`/`dRm` always live), with the λ₀ material SLOPES
upwinded off during reversal (`w_mat`: 1 for `u_out ≥ 0` so every
forward-flow history is bit-identical — recipe rows reproduce to the digit,
measured; 0 under firm reversal; 10⁻³c continuity band). Gate: driver C13
a/b/c, static, red-measured against both defective closures. Full story:
`Docs/NSCBC-reversal-branch-defect.md`. Fallout re-measured: the
failure-mode rows in FlameOutflow/DRM READMEs shifted +14–22% (the old
branch was bleeding pressure through reversals it manufactured itself —
counters now show e.g. the σ=1/eT=0/β_s=1 sitting row holds its ENTIRE face
in reversal, which is why those rows were closure-sensitive); recipe rows
unchanged to the digit; chamber table regenerated (README).

## Work queue

1. **T5 + T3** — injection fidelity: DONE (2026-08-25), both halves. Driver:
   `./nscbc1d t1|t3|t5` duct modes, self-gated, tables in `Verification/NSCBC1D/README.md`
   (t1: Dupuy's CLR curve, optimum relax_u = 0.3 at 6.06 t_a; t3: relax_u = 0 injects
   nothing, off-resonance injection monotone but unfaithful — 0.14/0.95/2.71 of target
   at relax_u 0.5/2/5 — and ON-resonance injection collapses at any stiffness, a velocity
   relaxation cannot drive a velocity node; t5: standing-wave geometry exact, amplitude
   carries the coupling story). PeleC: `NSCBC-Acoustic` duct mode (`pulse_type = 2`,
   `nscbc-acoustic-duct.inp`, `duct_metrics.py`) reproduces the driver's full 12-entry
   deterioration table to 1–4% (several entries to the third digit) and the antinode
   amplitudes to 0.1–0.3% — the deterioration is the formulation's, not either code's.
   Measurement note: detrend before Fourier-projecting at small relax_u (mean wander
   leaks as spurious amplitude; measured up to I_u = 0.73 of pure leak). The t3 table
   is the standing motivation for queue item 4's NRI-style registers.
2. **Phase-1 coverage** — DONE (2026-08-25), all four sub-items. Sustained recirculation:
   `bc_nscbc_backflow_material` ramps the entering lambda_0 content to the Target's
   reservoir state over an outward-Mach band [1e-3, 1e-2] (breathing stays frozen);
   gate C14, the flush test — frozen closure holds a hot duct at 597 K on its own
   exhaust, the flag flushes it to 299.5 K against a 300 K reservoir, and the ramp is
   bit-inert at |M| ~ 1e-4. Supersonic inflow without `Target.p` is a COUNTED
   substitution (`target_incomplete`, C6-gated both ways). Counters can no longer be
   silent: with `sum_interval` unset, the report runs every 100 coarse steps and prints
   only when something counted. Wall/NSCBC corner: `nscbc-acoustic-corner.inp` (y-walls
   against the x-hi characteristic outflow) — R = 0.752%, identical to the periodic-y
   baseline to the third digit, y-structure 0.0000%, max|v| 0.000%: the corner is
   invisible. Driver: air 61/61, LiDryer 65/65, SRK 45/45.
3. **T7 mini-SydGex** — duct set DONE; EB follow-ups BUILT, production IN PROGRESS
   (2026-08-25, may have been interrupted by a power loss — CHECK
   `Exec/RegTests/NSCBC-Chamber/box/` and `baffle/` run.logs before trusting any
   plotfiles there). State: `chamber-box.inp` (chamber as interior EB in a 2.4 × 0.9
   box — the lateral plenum) and `chamber-baffle.inp` (+ Sydney baffle, 67% blockage);
   one registered geometry `nscbc-chamber-box`; `eb_zero_body_state = 1`. THREE
   production attempts failed and were each diagnosed and fixed in commits: (i) the
   Godunov path NaN'd at an outer cut-cell corner at 1.4 ms → both inps run
   `pelec.do_mol = 1` (the EB-supported hydro; box-vs-duct A/B therefore differs in
   scheme too); (ii) grid-aligned EB surfaces → placed off-grid at x.4·dx via new
   `prob.chamber_yc` (x0 = 0.205, yc = 0.455, kernel follows; dims still 1.2 × 0.3);
   (iii) THE REAL KILLER: PeleC's EB walls default to ISOTHERMAL AT 1 KELVIN
   (`eb_isothermal = true`, `eb_boundary_T = 1.0`) — the cut cells refrigerated
   (268 K by 25 µs → 77 K by 1.6 ms), the flame NaN'd at the coldest corner, under
   any scheme/alignment. Diagnosed by a 500-step A/B (identical cooling with
   bc_nscbc off and either body state; 298.0 K exactly with adiabatic walls). Both
   inps now set `pelec.eb_isothermal = 0`, matching the adiabatic domain walls.
   Production COMPLETE (2026-08-25): both variants ran the full 8 ms, zero NaN
   (box 168k steps/2.7 h, baffle 174k/3.4 h at 6 ranks; ignition IC clipped to the
   chamber interior after Marc caught the kernel poking through the closed-end
   slab). A/B tables in the case README, measured with the cavity-masked
   `chamber_qoi.py --xvent 1.405 --xlo 0.205 --ymask 0.305,0.605`. Headlines: the
   lateral plenum is a ~3× better vent than the straight duct (chamber dp a third,
   burn ~40% slower); the baffle produces the full Sydney phenomenology — 6× gap
   jet, 8.5× overpressure peak (31.7k dyn/cm² at 3.7 ms), sustained vent breathing
   (14.1M reversal fills vs the box's 1.5M, all through the unified closure).
   Neither un-baffled variant rings down inside 8 ms (the extension charge burns
   late); the baffle run settles by 8 ms. Building all this also found the
   outflow's equilibration micro-inflow (−0.28 cm/s — real, benignly counted) and
   put a 1e-9 c roundoff deadband on the reversal counter.
4. **Boundary registers** (design note before code): per-face EMA/integrated registers,
   checkpointed, updated once per advance outside the fill. NDNR motivation at outlets is
   WEAKENED by the flame-closure results; strongest remaining cases are inlets (NRI) and
   reflection removal where sigma=16 is still chosen. The C11x closure may want the trend
   register as its transit/decay gate.
5. Hardware-gated: CUDA/HIP compile, CPU-vs-GPU on NSCBC-Acoustic, register-spill profile.
6. Upstream: curate the branch into an AMReX-Combustion PR when ready.

Start by confirming `git status` is clean (untracked `PAPERS/` is expected) and the driver
runs green, then pick up the queue at item 1 — or at item 3's write-up if the Lesson-9
production plotfiles have arrived — unless I say otherwise.
