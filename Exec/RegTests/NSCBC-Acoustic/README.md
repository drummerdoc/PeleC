# NSCBC-Acoustic

A right-running isentropic acoustic pulse in a quiescent gas, launched toward a
subsonic outflow at `x-hi`. `x-lo` is a wall and `y` is periodic, so the only
thing that can absorb the pulse is the outflow boundary. This is the AMReX-side
counterpart of check C4 in `Verification/NSCBC1D`: the same physical problem,
solved by a completely different code path, as a cross-check that the boundary
condition survives contact with the framework.

## Running it

```sh
./PeleC-NSCBC-Acoustic nscbc-acoustic.inp                      # characteristic
./PeleC-NSCBC-Acoustic nscbc-acoustic.inp pelec.bc_nscbc=0     # hard boundary
```

With `bc_nscbc = 0` the problem's `bcnormal` imposes the ambient pressure
directly in the ghost cells — deliberately crude, so the comparison is stark.
With `bc_nscbc = 1` the `bcnormal_nscbc` hook in `prob.H` returns an outflow
target and the characteristic treatment takes over. Note how little the hook
has to say: for a pure non-reflecting outflow the target pressure is the only
quantity that may be specified, so it is the only quantity it sets.

## Results

`L = 10 cm`, `c = 34719 cm/s`, so one acoustic transit is `2.88e-4 s`. Run to
`t = 4.6e-4 s` (1.6 transit times), by which point the pulse has crossed the
outflow. `R` is the peak residual pressure disturbance in the upstream half,
measured against the instantaneous domain mean so the σ-driven anchoring
adjustment is not miscounted as a reflected wave, divided by the incident
amplitude.

| boundary treatment | R [%] | mean p at end [dyn/cm²] |
|---|---|---|
| hard: ambient p imposed in the ghost cells | **97.2** | 1013178.3 |
| NSCBC, σ = 0.25 | **0.81** | 1013241.8 |

Target pressure is 1013250.0. The hard boundary reflects essentially the whole
pulse; the characteristic treatment reduces the reflection by a factor of 120
and holds the mean pressure to 8 parts per million of the target.

**Cross-check.** The standalone 1-D driver predicts `R = 0.758%` at σ = 0.25
(`Verification/NSCBC1D/README.md`). PeleC's 2-D Godunov solve gives 0.810%.
Two independent solvers, two independent implementations of the surrounding
machinery, agreeing to 7% on a quantity that spans two orders of magnitude
between the good and bad boundary conditions. That agreement is the point of
having both.

## Things worth trying

* `pelec.bc_nscbc_sigma = 0` — perfectly non-reflecting; R drops further, and
  the mean pressure is then unanchored and free to drift over a long run.
* `pelec.bc_nscbc_sigma = 2.0` — over-relaxed; the reflection grows roughly
  linearly in σ.
* `pelec.bc_nscbc_order = 1` — zeroth-order extrapolation of the outgoing
  invariant instead of minmod-limited linear.
* `pelec.bc_nscbc_pin_farfield = 1` — the value-pin formulation, which is both
  non-reflecting and anchored but does not converge under mesh refinement.
* Refine `amr.n_cell` and confirm R does not grow: the relaxation is a rate,
  not a per-cell value blend, so it is grid-converged.

## Three dimensions

The same sources build in 1-D, 2-D and 3-D — `AMREX_SPACEDIM`, `AMREX_D_DECL`,
`AMREX_D_TERM` and `AMREX_D_EXPR` throughout — so only the inputs file changes.
The 2-D planar result is unmoved: mean pressure at the end is 1013178.3 with a
hard boundary and 1013241.8 with the characteristic one, the same to the last
printed digit as before the conversion.

| inputs | what it is |
|---|---|
| `nscbc-acoustic.inp` | 2-D planar; the historical regression |
| `nscbc-acoustic-3d.inp` | 3-D planar; confirms the kernel builds and runs in 3-D |
| `nscbc-acoustic-3d-radial.inp` | 3-D radial, **all six faces characteristic** |

### Why the radial case exists

A planar pulse loads one face at normal incidence, which is what the 2-D run
already tests. It says nothing about a ghost cell that lies outside the domain
in two or three directions at once — and that is the part of `BCfill.cpp` with
no 1-D analogue, so the standalone driver cannot reach it either.

A radial pulse reaches the six faces at normal incidence, the twelve edges at
45°, and the eight corners along the body diagonal, in one run. Because the
initial condition is exactly isotropic about the box centre, **any** departure
from sphericity in the departing front is boundary-generated; no reference
solution is needed to say so.

`metrics.py sphericity` reports the spread in the front radius over ~4000
directions, binned by χ, the angle to the nearest face normal: χ = 0 is a face
centre, 45° an edge, 54.7° a corner. Rays whose front has already left are
dropped, so as time goes on the surviving rays are exactly the oblique ones.

### What it measures

`c = 34719 cm/s`, box 5 cm, 96³. The front reaches the faces at 7.2×10⁻⁵ s, the
edges at 1.02×10⁻⁴ and the corners at 1.25×10⁻⁴.

| t [s] | rays | χ range | radius spread % | amplitude spread % |
|---|---|---|---|---|
| | | | NSCBC / hard | NSCBC / hard |
| 3.0×10⁻⁵ | 4000 | all | 7.355 / 7.355 | 1.782 / 1.782 |
| 6.0×10⁻⁵ | 4000 | all | 2.891 / 2.891 | 5.158 / 5.157 |
| 7.5×10⁻⁵ | 2260 | edges + corners | 1.964 / 2.452 | **4.70 / 7.34** |
| 9.0×10⁻⁵ | 555 | corners | 1.291 / 1.612 | **7.69 / 40.29** |
| 1.05×10⁻⁴ | 25 | corners | 0.600 / 0.523 | **2.25 / 22.14** |

Read the first two rows first: before the front reaches any face the two runs
are **identical to six figures**, as they must be, since nothing has touched the
boundary yet. That is the metric's own sanity check, and it is why the later
rows can be believed.

After the face crossing they separate, and the discriminator is the **amplitude**
spread, not the radius: the front arrives at the right time either way, but a
hard boundary corrupts its strength. At 9×10⁻⁵ s — corner-bound rays only, the
face and edge parts of the wave already gone — the amplitude around the
surviving arc varies by 40% with a hard boundary and 7.7% with the
characteristic one.

So corner and edge ownership works. The `apply()` algebra is dimension-agnostic
by construction, and this says the plumbing around it is too.

## AMR and EB variants

`nscbc-acoustic-amr.inp` runs the planar pulse from inside a 2× refined patch
kept away from the outflow: mean pressure and upstream residual match the
single-level run (1013241.8 / 0.751% vs 0.752%). PeleC now warns — once per
face — when a refined level touches a Hard/UserBC face with `bc_nscbc = 1`,
because the fill's stencil is level-local and a fine patch on a
characteristic face makes the boundary condition level-dependent. Measured
here the on-face artefact is small (0.744%, +0.1 dyn/cm²); nothing
guarantees that at higher ratios or oblique incidence.

`nscbc-acoustic-eb.inp` seats an EB solid inside the fill's stencil at the
outflow. It found two things: PeleC's *default* body state is a sampled
fluid state the fill cannot detect (the counter read zero with a solid in
the stencil — the silent fallback the counters exist to prevent), so
`bc_nscbc` with EB geometry now *requires* `pelec.eb_zero_body_state = 1`
and aborts otherwise; and a body cutting the domain face itself NaNs under
the characteristic *and* the hard boundary — a pre-existing PeleC
EB-at-domain-boundary limitation. With the flag set, the run counts ~56000
body-state stencil degradations and completes cleanly.

### Two things this case turned up

**The flow-reversal fallback is not hypothetical.** Running with
`pelec.sum_interval > 0`, the counters read zero everywhere until the pulse
leaves and then report `flow reversal 11456` in a 40-step window. That is
correct physics, not a bug: the rarefaction behind an outgoing spherical wave
pulls the pressure below ambient and draws gas back in through faces configured
as outflows. The kernel detects it and switches to the ambient-pressure closure,
which keeps the interior composition and tangential velocity so a grazing flow
is not arrested — and the run is stable through it. Before this case, that path
had only ever been exercised by a synthetic state in check C6.

**It also caps the benefit.** For cells in reversal the boundary is effectively a
pressure Dirichlet, so it is not non-reflecting there. That is why the residual
after the wave has gone improves only ~1.5× over a hard boundary here
(L2 0.038 vs 0.056 of the incident amplitude) against 120× for the planar case.
If a problem spends a lot of its time in reversal at an outflow, the honest fix
is to give that face an inflow target rather than to expect the outflow model to
cope. The residual number is also not purely boundary error in this case: the
imploding half of the split pulse passes through the origin and is still in the
box at the final time.
