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
  minutes to ~2 h. Run sims in the background and poll. `rm -rf tmp_build_dir *.ex` before
  switching `Chemistry_Model` — make will not notice on its own.
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
* Outflow reversal (transient backflow through an outflow) uses a soft σ-rated pressure
  relaxation in the ghosts, continuous through u = 0 — NOT a pin to ambient. The old hard pin
  NaN'd the first NSCBC-Chamber production run at 4.08 ms (vent breathing during ring-down)
  and was the cause of both old NaN rows in the flame-exit table. Every previously-completing
  table row in the suite reproduces to the digit under the new closure (healthy runs never
  reverse); only failure-mode rows shifted, and all READMEs carry post-fix numbers. A
  sustained-recirculation gate for this branch is still queue item 2.
* C11x (driver, reported not gated): profile-fit ghost closure (fitU + source-consistency
  bound) reaches ~oracle level on a sustained front, releases unsustainable structure, robust
  to 15% end-state and full shape distortion of the family. 2-3D candidate; open items in
  `Verification/NSCBC1D/README.md`.
* Driver: air 53/53, LiDryer 57/57, PELE_NUM_ADV=2 57/57, SRK 40/40 (statics; dynamics
  skipped by design). CI job `NSCBC-Driver` runs all four; GNUmake takes `PELE_EOS=SRK`.
* Static GPU audit of the branch code is clean (commit e6a8ad9): device paths, captures,
  counters, and every synchronize reviewed. The actual device build remains item 5.

## Local artifacts (this machine, not in git)

* `Exec/RegTests/NSCBC-{FlameOutflow,COVO}/runs/` hold the animation runs and MP4s
  (flame exit 3-way, vortex, pulse) — git-ignored, reproducible from the inputs.
* The teaching page "Waves at the Open Boundary" (σ/β tutorial with live figures) is a
  claude.ai artifact: https://claude.ai/code/artifact/4d4f11af-7932-4c80-95f7-cf9bd90b984a
* Untracked `PAPERS/` at the repo root is Marc's reference library — expected, leave it.

## Work queue

1. **T5 + T3** — injection-fidelity tests (harmonic inlet forcing vs analytic standing-wave
   P_RMS; Daviller's deterioration index): what the value-relaxation inlet does to injected
   signals. Driver first (add the sigma-from-CLI + T1/T3/T5 duct modes together), PeleC after.
2. **Phase-1 coverage**: backflow branch as local inflow with a sustained-recirculation test;
   wall/NSCBC corner test; supersonic-inflow `Target.p` handling; counters on by default.
3. **T7 mini-SydGex** — BUILT as `Exec/RegTests/NSCBC-Chamber` (Lesson 9 is runnable).
   Marc's first vent production run found the reversal-pin defect (fixed; see State) and
   its data is discarded; BOTH variants are to be regenerated with a post-fix executable
   (2026-08-24). When those plotfiles come back: process with `chamber_metrics.py`
   (`--xvent 1.2`), build the peak-aligned trace difference / V̇ balance / ring-down
   table into the case README's Status section, and consider the σ = 16 ring-down
   variant the README names. Follow-ups: laterally-expanding plenum and the baffle
   (both need EB).
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
