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
  minutes to ~2 h. Run sims in the background and poll. `make realclean COMP=llvm` before
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
2. **Phase-1 coverage**: backflow branch as local inflow with a sustained-recirculation test;
   wall/NSCBC corner test; supersonic-inflow `Target.p` handling; counters on by default.
3. **T7 mini-SydGex** — DONE (2026-08-25). Three-variant production set (vent, plenum,
   vent-σ16) completed with the final kernel; T7 table and phase-resolved analysis are in
   the case README's Status section; `chamber_qoi.py` (beside `chamber_metrics.py`) is the
   measurement script. The case's production runs killed three defective reversal closures
   along the way. Remaining follow-ups: laterally-expanding plenum and the baffle (both
   need EB).
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
