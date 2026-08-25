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
4. **Boundary registers** — DESIGN NOTE DONE with a measured prototype
   (2026-08-25): `Docs/NSCBC-boundary-registers-design.md`. Headlines: a register-only
   NRI transplant FAILS in the GC form (no amplitude slot to protect — I.3e, measured);
   the slot is `Target::dudt`, stateless and ALREADY IN THE KERNEL (usable by any
   problem hook today; zero = classical inlet identically); the register's true job is
   the learned reference mean (EMA of R₋, τ = 3 t_a). Register + dudt measured in the
   driver-held-register prototype: I_in = 0.89–1.08 over the full 3-stiffness ×
   3-frequency matrix including ON resonance, vs 0.07–4.7 classical — Daviller's I ≡ 1,
   transplanted. PHASES A–C IMPLEMENTED IN PELEC (2026-08-25): flat replicated-band
   store in `BCfill.cpp` (`pelec.bc_nscbc_nri` inlets, `pelec.bc_nscbc_ndnr` outlets;
   update once per level-0 advance from the PeleC::advance dispatcher, BCfill composes
   the effective Target, kernel stays pure; trend_dS advisory only). Gates green: PeleC
   duct NRI point I_in = 0.970 (driver 1.031, FF-only 2.42); planar σ=16 NDNR 30.10% →
   3.66% with the mean anchored to 37 ppm (driver 28.14% → 3.56%); restart bit-identity
   (`NSCBCRegisters` in the checkpoint); n1-vs-n6 bit-identity incl. tangentially-split
   register faces. Two paid-for lessons live in the design note's implementation
   record: the invariant needs a TIME-FROZEN impedance (`ema_rhoc`; local rho·c aliases
   the mean p/(ρc) into R_out at signal amplitude — degraded I_in 1.03→1.6 in BOTH
   codes), and `amrex::DefaultGeometry().Domain()` is not a reliable level test in the
   fill (the store now records its own domain; peek(dom) is the level test). Testing
   traps: NSCBC-Acoustic's GNUmakefile is `USE_MPI=FALSE` — six mpiexec clones race one
   plotfile; multi-rank plotfiles shard `Cell_D`, so compare fields, never files.
5. Hardware-gated: CUDA/HIP compile, CPU-vs-GPU on NSCBC-Acoustic, register-spill
   profile. STATIC GPU AUDIT DONE (2026-08-25, pre-hardware, NVIDIA-focused) — all
   branch device code surveyed and one sync gap fixed:
   * CLEAN: no printf/Abort/std::string/allocation on device paths; finiteness via
     `amrex::Math::isfinite`; every device lambda captures PODs by value, none touches
     a class member (no hidden `this`); `PCHypFillExtDir` is all-POD (pointers + POD
     Params in a GpuArray); reg-index helpers are `AMREX_GPU_HOST_DEVICE`, dimension-
     generic, and internal-linkage-safe (plain functions — nvcc's extended-lambda
     linkage restriction applies to lambdas, and every ParallelFor/ReduceOps lambda
     sits in an external-linkage member function); counters bump via
     `HostDevice::Atomic::Add(Long*)`; ProbParmDevice structs are POD with raw device
     pointers, containers live in ProbParmHost (chamber PMF follows the standard
     pattern); `std::abs` on device matches upstream precedent (Riemann.H/Godunov.H/
     EB.H) and is supported on CUDA/HIP; NO `AMREX_USE_GPU` ifdefs anywhere in branch
     code — CPU and GPU compile the identical path, so the A/B is meaningful.
   * SYNC AUDIT (also clean, one fix): update ends `streamSynchronize`; MFIter loops
     sync implicitly at destruction; all `Gpu::copy` calls are the blocking variant;
     wrap-gate uses `rd.value(op)` (syncs internally). FIXED: the counter reset after
     a report is a device fill on the current stream, unordered against the next
     advance's fills on MFIter's rotating streams — a zero could land on a fresh
     count. Now `Device::synchronize()` before the read+reset and `streamSynchronize`
     after (diagnostics-only race, but counters are evidence and must not lie).
   * TOMORROW'S PROTOCOL: (1) build — `make realclean; make -j TPL USE_CUDA=TRUE;
     make -j USE_CUDA=TRUE USE_MPI=TRUE CUDA_ARCH=<sm>` per case (Make.PeleC passes
     USE_CUDA through to the SUNDIALS TPL; TPL and app configs must match; remember
     the case GNUmakefiles default USE_MPI=FALSE). (2) correctness — NSCBC-Acoustic
     planar σ=0.25 row, σ=16 NDNR row, duct NRI point, then a short chamber-box run
     (EB + do_mol on GPU); compare CPU-vs-GPU field-by-field with fielddump (%.17e):
     expect agreement to roundoff, NOT bitwise (GPU FMA contraction differs — do not
     chase bit-identity across architectures; n1-vs-n6 bit-identity is the
     determinism gate and it is CPU-side). (3) profile — `CUDA_VERBOSE=TRUE` for
     ptxas register counts on the fill kernel (it carries Y[NUM_SPECIES] locals and
     the EOS; spill is the expected risk with LiDryer/drm19), then ncu occupancy if
     it spills. The 1-D driver is a host tool — no GPU work there.
6. Upstream: curate the branch into an AMReX-Combustion PR when ready.

Start by confirming `git status` is clean (untracked `PAPERS/` is expected) and the driver
runs green, then pick up the queue at item 5 (hardware-gated GPU work) or item 6
(upstream PR curation) — unless I say otherwise.
