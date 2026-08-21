# PeleC — working notes for Claude Code (branch `nscbc`)

PeleC is an AMReX-based compressible reacting-flow solver (C++17, CPU/GPU, CGS units).
Branch `nscbc` implements a ghost-cell Navier–Stokes Characteristic Boundary Condition.
Owner: Marc. Repo at `/Users/marcusd/src/PeleC`, fork `origin` = drummerdoc/PeleC,
`upstream` = AMReX-Combustion/PeleC. Base of the branch: upstream `development` @ `ac17bd0`.

## Read these before touching NSCBC code, in this order

1. `Docs/NSCBC-status-plan-and-curriculum.md` — current status, ranked defects, the five-phase
   plan, the nine-lesson teaching curriculum. **Start here. It supersedes `RESTART-PROMPT-nscbc.md`
   where they disagree.**
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

## Known defects to fix first (Phase 0 of the plan) — see the status doc §2

- `Source/NSCBC.H:1029`: reaction source must be `L_in += (1-beta_s)*S_p`, not `-=`.
  Verified convention-free in `Verification/NSCBC1D/source_sign_check.py`. After the flip, re-measure
  C12, the FlameOutflow/flame-exit tables and COVO before trusting any old number.
- ~~`Source/BCfill.cpp`: tangential stencil in periodic directions~~ FIXED, but not by wrapping:
  amrex's corner protocol (`StateDataPhysBCFunct` recomputes corner ghosts on an image-band strip
  FAB) makes a centred difference across the periodic seam unachievable. The fill restricts its
  tangential stencil to the strip-consistent band (`tang_range`); the gate is
  `nscbc_check_periodic_wrap()` + `NSCBC-COVO/nscbc-wrapgate.inp` (both faces bitwise, both
  decompositions bitwise after 20 steps). Do not "fix" this back to a wrap.
- Scheme is NOT Motheau's NSCBC-GC (invariant extrapolation, not derivative targets) — fix the sphinx docs.
- Daviller 2019's σ is half of this code's σ (his 2K vs our K in the wave).

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
