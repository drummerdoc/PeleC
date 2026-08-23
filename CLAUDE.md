# PeleC — working notes for Claude Code (branch `nscbc`)

PeleC is an AMReX-based compressible reacting-flow solver (C++17, CPU/GPU, CGS units).
Branch `nscbc` implements a ghost-cell Navier–Stokes Characteristic Boundary Condition.
Owner: Marc. Repo at `/Users/marcusd/src/PeleC`, fork `origin` = drummerdoc/PeleC,
`upstream` = AMReX-Combustion/PeleC. Base of the branch: upstream `development` @ `ac17bd0`.

## Read these before touching NSCBC code, in this order

1. `RESTART-PROMPT-nscbc.md` — the one-page current state and work queue. **Start here.**
   (`Docs/NSCBC-status-plan-and-curriculum.md` is the Phase-0-era audit and teaching curriculum;
   its defect list is executed and its tables pre-date the fixes — historical record plus curriculum.)
2. `Docs/NSCBC-design-and-literature-review.md` — formulation (Part I), papers (II–IV).
3. `Source/NSCBC.H` (kernel), `Source/BCfill.cpp` (plumbing, `nscbc_fill`), the READMEs in
   `Exec/RegTests/NSCBC-*` and `Verification/NSCBC1D/README.md` (measured tables and protocols).
4. `git log ac17bd0..nscbc` — the commit messages carry the reasoning.

## Decisions settled by measurement — do not relitigate without new data

- Ghost-cell architecture stays; flux-form reformulation rejected (GC bias is ~5 % of the error).
- `extrap_temperature` is the correct ghost closure for the diffusion operator (C8, C12).
- `extrap_material` is for fronts that SIT, never during a transit.
- Determinism (decomposition, MPI, restart) is bit-identical and must stay so.
- EB requires `pelec.eb_zero_body_state=1`; AMR face touching an NSCBC face only warns.

## State (Phase 0 complete; see git log from `1329779` for the record)

- The reaction source enters the incoming wave with `+(1-beta_s)*S_p` (Sutherland–Kennedy);
  gate: `Verification/NSCBC1D/source_sign_check.py`. All beta_s-sensitive tables were re-measured
  after the fix; READMEs carry only post-fix numbers.
- Flame-on-boundary recipe: `bc_nscbc_extrap_temperature=1` + `bc_nscbc_beta_s=0` holds
  hard-outflow accuracy at ANY sigma (sitting) and survives transits at any sigma. sigma=16+beta_s=1
  is the best absolute anchoring at 28% reflection. Do not stack sigma=16 with beta_s=0.
- Periodic tangential stencils CLAMP into the domain (the original scheme). Three designs measured:
  the clamp (stable everywhere; O(1e-4) seam residual under amrex's corner-strip protocol), a wrap
  through resident images (aperiodic at 2e-2), and a strip-consistent band (bitwise-periodic but
  destabilises NSCBC-FlameOutflow-DRM at beta<1: NaN in ~150 steps from a two-row stencil change).
  Stability won. Gate: `nscbc_check_periodic_wrap()` reports the residual, aborts only above 1e-3;
  exercised by `NSCBC-COVO/nscbc-wrapgate.inp`. Do not rebuild the wrap or the band without a DRM
  stability result.
- C11x (driver, reported not gated): a profile-fit ghost closure (fitU + source-consistency bound)
  holds a sustained front at ~oracle level, releases unsustainable structure, robust to family
  errors. Candidate for a 2-3D closure; open items listed in `Verification/NSCBC1D/README.md`.

## Build and run

- Per-case GNUmake in `Exec/RegTests/<case>`: `make -j TPL` once (SUNDIALS), then `make -j`.
  Regression inputs and measurement scripts live beside each case; READMEs give the protocol.
- 1-D driver: `Verification/NSCBC1D` (needs a 1-D AMReX; see its README). `./nscbc1d` runs C1–C12;
  `./nscbc1d sweep` prints the σ-reflection table. Four builds: air, LiDryer, PELE_NUM_ADV=2, PELE_EOS=SRK.
- `Verification/NSCBCFields/fielddump` + `metrics.py` (modes: residual, circularity, sphericity).
- Run long builds and simulations in the background and poll; do not block on them.
- `prob.pmf_datafile` is CWD-relative: always pass an absolute path.
- `measure.py`'s shielded reference is valid only to t = 2.4e-5 s and needs `prob.nscbc_inflow=0`.
- AMReX CLI args with spaces must be single argv tokens.
- Regression gates in `Tests/CMakeLists.txt` are gold-file comparisons; physics numbers live in READMEs.
  When a fix moves a number, update the README table and say in the commit message what moved and why.

## Style

- Commit messages in the existing narrative style: what was measured, the numbers, what it settles.
- Keep `Source/NSCBC.H` GPU-clean: no printf/Abort/std::string/allocation on device paths;
  `amrex::Math::isfinite`, not `std::isfinite`; POD captures only.
- `_to_delete/` at the repo root is scratch and may be purged.
