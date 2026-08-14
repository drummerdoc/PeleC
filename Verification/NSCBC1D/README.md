# nscbc1d — standalone verification of `Source/NSCBC.H`

This driver compiles the production NSCBC kernel **unmodified** against a
minimal 1-D finite-volume Euler solver. Its purpose is to separate two failure
modes that are otherwise indistinguishable in a full PeleC run:

> *"the boundary condition is wrong"* vs *"the AMReX plumbing around it is wrong"*

It is cheap enough to sweep parameters, so it also produces the reference curves
that the AMReX-side regression tests are checked against.

## Build and run

```sh
cmake -S . -B build -DAMReX_DIR=<amrex-install>/lib/cmake/AMReX
cmake --build build
./build/nscbc1d          # run all checks
./build/nscbc1d sweep    # also dump the sigma sweep
```

Any AMReX build with `AMReX_SPACEDIM=1` works; MPI, OpenMP, EB and particles are
all unnecessary. The default mechanism is `air` (2 species, Fuego); override
with `-DPELE_MECHANISM=`. Build it a second time against a reacting mechanism
to exercise check C7, which is skipped otherwise:

```sh
cmake -S . -B build_lidryer -DAMReX_DIR=... -DPELE_MECHANISM=LiDryer
cmake --build build_lidryer && ./build_lidryer/nscbc1d      # 26/26
```

Nothing in the driver may assume a particular mechanism. The first version of
`air_Y()` returned "0.233 for O2, else 0.767", which is right for the
two-species `air` mechanism and gives a composition summing to 6.6 for
LiDryer's nine — every check downstream then failed for reasons having nothing
to do with the boundary condition.

## Units

PeleC and PelePhysics work in **CGS**: cm, g, s, K, dyn/cm², erg/g.
`Constants::PATM = 1.01325e6`. The first version of this driver was written in
SI and every sound speed was 100× too large, which silently turned the
supersonic-outflow check into a subsonic one — that is exactly the class of
error this driver exists to catch, and it is why the checks assert on physical
relationships rather than on hard-coded numbers.

## What is checked

| | Check | What a failure means |
|---|---|---|
| **C1** | A uniform state is reproduced exactly in every ghost layer, at every σ, for both inflow and outflow | The kernel is manufacturing a gradient out of nothing; everything downstream is meaningless |
| **C2** | Every relaxation moves the boundary *toward* its target, on both the lo and hi face | A sign error. This is the assertion the legacy Fortran lacked, and its absence is why users had to be told "`relax_T` must be negative" |
| **C3** | Outflow *extrapolates* composition rather than imposing it; inflow imposes it exactly; `Σ Y = 1` and `UEDEN = UEINT + KE` hold to round-off | Species over-specification at outflow (the legacy defect), or a broken state identity |
| **C4** | Acoustic reflection of a pressure pulse: below 1% at σ=0.25, essentially zero at σ=0, and 2nd-order extrapolation no worse than 1st | The extrapolation or the invariant algebra is wrong |
| **C5** | The relaxation is a **rate**: grid-independent, and equal to `K = σ(1−M²)c/L` | The parameterisation has drifted to a value-blend, whose effective rate scales as `c/Δx` and therefore doubles when the mesh does |
| **C6** | Supersonic, reversed-flow and EB-body-state fallbacks each return a finite physical state and increment their counter | A silent fallback — i.e. a bug that will not be found in production |
| **C7** | The closed-form reaction source `dp/dt|_react` matches a directional finite difference along the reaction path, converging as τ→0; chemistry conserves mass; a cold state gives zero | The thermodynamics of `reaction_dpdt()` is wrong. Skipped automatically when the mechanism has no reactions |

## Reference results

Measured with `air`, `n = 400`, `L = 10 cm`, `order = 2`, a 0.1% Gaussian
pressure pulse, integrated for 1.6 acoustic transit times. `R` is the peak
residual wave amplitude in the upstream half after the pulse has left, measured
against the *instantaneous domain mean* so that the σ-driven anchoring
adjustment is not miscounted as a reflected wave.

| σ | R [%] | mean p drift [dyn/cm²] | τ_relax [s] |
|---|---|---|---|
| 0.00 | 0.00078 | −0.0076 | ∞ |
| 0.05 | 0.160 | −1.68 | 5.75e−3 |
| 0.10 | 0.315 | −3.30 | 2.87e−3 |
| 0.15 | 0.466 | −4.88 | 1.92e−3 |
| 0.20 | 0.614 | −6.44 | 1.44e−3 |
| **0.25** | **0.758** | **−7.96** | **1.15e−3** |
| 0.30 | 0.899 | −9.44 | 9.58e−4 |
| 0.50 | 1.43 | −15.0 | 5.75e−4 |
| 1.00 | 2.56 | −26.8 | 2.87e−4 |
| 2.00 | 4.18 | −43.5 | 1.44e−4 |
| 4.00 | 7.20 | −60.3 | 7.19e−5 |
| 8.00 | 14.97 | −69.3 | 3.59e−5 |
| 16.00 | 28.14 | −70.8 | 1.80e−5 |
| — (`pin_farfield`) | 0.015 | +0.010 | n/a (value pin) |

The sweep runs past σ = 2 because `Exec/RegTests/NSCBC-FlameOutflow` needs
σ ≈ 10 to anchor an outflow with a flame crossing it, and the price of that has
to be quotable: 20-30% reflection, and a mean drift that has saturated — beyond
σ ≈ 8 the anchoring stops improving while the reflection keeps growing, so
σ > 16 buys nothing at all.

Three things to read out of that table.

**The reflection/anchoring trade-off is real and monotone.** `R` grows very
nearly linearly in σ while the pressure anchoring strengthens in step. σ = 0 is
perfectly non-reflecting and completely unanchored. This is the curve that makes
"σ = 0.25 is a good default" a measured statement rather than a received one.

**The relaxation is a genuine rate.** Doubling and quadrupling the resolution
changes the measured `K` by 3% (1925 → 1861 s⁻¹ from n=200 to n=800), and the
measured value sits within 7% of `σc/L`. A boundary condition parameterised as a
*value blend* instead has an effective rate of
`c/Δx`, which doubles when the mesh does, so its σ is not transferable and its
behaviour is not grid-converged. Keeping this check green is what keeps
literature σ values meaningful in PeleC.

**`pin_farfield` is not simply "σ = ∞".** The hard far-field pin is *both*
nearly non-reflecting (R = 0.015%) *and* anchored (drift +0.01 vs −7.96 at
σ=0.25) — because it constrains the incoming invariant's *value* rather than its
gradient, and a value constraint on the incoming characteristic alone reflects
nothing. Its cost is that it anchors to `p_target + ρc·u_out` rather than to
`p_target`, and that it does not converge under refinement. It is the right
choice for an open boundary onto a large quiescent reservoir and the wrong one
for a duct exhausting into a plenum whose true mean pressure is not `p_target`.

## Adding a check

Checks are plain functions that call `check(bool, name, detail)`. Prefer
assertions on *physical relationships* (monotonicity, grid-independence,
conservation identities, direction of relaxation) over assertions on numbers:
the numbers move when the mechanism, the mesh or the interior scheme changes,
and the relationships do not.
