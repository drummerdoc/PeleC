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
