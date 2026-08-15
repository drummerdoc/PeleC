# NSCBC-FlameOutflow — a flame sitting *on* a characteristic outflow

A wrinkled, unanchored premixed sheet (LiDryer, φ = 0.4 H₂/air, S_L = 22.8 cm/s)
in a uniform stream at U = 2 S_L, with its mean position exactly on the outflow
plane. Half the boundary is in burnt gas, half in fresh, and the reaction zone
passes through the boundary cells at the two crossings — which is where the
wrinkle slope is largest, so the flame normal there is tilted 51° off the
boundary normal.

```sh
./PeleC2d.gnu.ex nscbc-flameoutflow.inp
./PeleC2d.gnu.ex nscbc-flameoutflow.inp pelec.bc_nscbc_sigma=16.0   # see below
```

The stream carries the sheet downstream at U − S_L = S_L, so nothing anchors it
and nothing has to be held in place. It is the configuration
`NSCBC-PMF/README.md` says is needed and cannot provide: a live heat-release
term, strong tangential gradients, and a multi-species diffusive structure, all
inside the boundary cells, from step 0.

What the case does **not** do is advect out within a measurement. The drift Mach
number is S_L/c ≈ 6×10⁻⁴, so clearing the boundary takes some 400 acoustic
transits. The sheet is initialised already straddling and is quasi-frozen over
the few relaxation times the boundary needs. That is not a compromise: the
quantity under test is the boundary's response to a reaction zone in its own
cells, and that is present throughout.

## The measurement protocol

Earlier attempts at this measurement were worthless, for two reasons worth
recording because both are easy to repeat.

**There is no useful self-referencing metric here.** The configuration has a
genuine mass imbalance — the front drifting downstream converts light burnt gas
into dense fresh gas, so the domain really does gain mass — and the resulting
pressure ramp is several times larger than any difference between boundary
settings. Comparing a run against its own domain mean measures the ramp, not the
boundary.

**A longer-domain reference is only a reference if it differs in one thing.**
The obvious reference — same problem, outflow moved downstream — differs in two
others: the physical ramp is spread over a larger volume, and, because
`L_ref = ProbHi − ProbLo`, *every* relaxation rate `K = σ(1−M²)c/L_ref` changes
with it, including the one at the inlet.

The protocol that works:

* `prob.nscbc_inflow = 0` in **every** run, test and reference. The inlet is a
  hard Dirichlet and the characteristic treatment is on the outflow alone, so
  the length-dependent inflow rate cannot contaminate anything.
* The reference puts its outflow at x = 3.0, i.e. 2.4 cm downstream of the test
  outflow. In burnt gas at c ≈ 8.6×10⁴ cm/s nothing from it can reach x < 0.6
  before t = 2.8×10⁻⁵ s. Runs stop at 2.4×10⁻⁵ s.
* Up to that time the reference restricted to x < 0.6 **is** the exact solution
  of the test problem: identical grid, identical scheme, identical inlet, and
  its outflow outside the domain of dependence. The difference field is the test
  run's outflow error.
* Because of that shielding the mean level is part of the error and must **not**
  be subtracted. Getting this wrong was the earlier mistake that made a hard
  Dirichlet outflow look better than it is.

`measure.py` implements it.

## What it measures

Errors against the shielded reference at t = 2.4×10⁻⁵ s (≈ 3.4 τ_relax at
σ = 1), in dyn/cm²; p_amb = 1.013×10⁶, so 1000 is 10⁻³ of ambient. `R` is the
acoustic reflection at the same σ from `Verification/NSCBC1D`.

| outflow | σ | β | β_s | order | mean Δp | L2(Δp) | R [%] |
|---|---|---|---|---|---|---|---|
| hard `p = p_amb` | — | — | — | — | −527 | 586 | — |
| characteristic | 0.25 | 0.5 | 0 | 2 | +2185 | 2209 | 0.76 |
| characteristic | 1 | 1 | 1 | 2 | +2301 | 2325 | 2.56 |
| characteristic | 1 | 0.5 | 1 | 2 | +2065 | 2088 | 2.56 |
| characteristic | 1 | 1 | 0 | 2 | +2260 | 2288 | 2.56 |
| characteristic | 1 | 0.5 | 0 | 2 | +2074 | 2098 | 2.56 |
| characteristic | 4 | 0.5 | 0 | 2 | +1450 | 1477 | 7.20 |
| **characteristic** | **16** | **0.5** | **0** | **2** | **+218** | **330** | **28.1** |
| characteristic | 1 | 0.5 | 0 | **1** | −2393 | 2408 | — |
| `pin_farfield` | — | 0.5 | 0 | 2 | +748 | 787 | 0.015 |

Five things come out of that table, and none of them is the thing this case was
built to look for.

**σ is the control that matters, and the default is far too small.** From
σ = 0.25 to σ = 16 the error falls by a factor of ten, and only at σ = 16 does
the characteristic outflow beat a hard pressure outflow. The default σ = 0.25,
which `Verification/NSCBC1D` shows is a good inert choice, is the *worst* entry
in the table.

**β = 0.5 helps by about 10%, consistently.** Every pairing that differs only in
β moves the same way and by about the same amount. That is the same β ≈ 0.5 the
COVO case measured, arrived at independently, which is worth something.

**β_s does essentially nothing.** 2301 → 2260 and 2065 → 2074: one improves by
2%, the other degrades by 0.4%. This is not a silent failure — the diagnostics
report `reaction source dropped 0` throughout — the correction is being applied
and its effect is below the noise. See below for why.

**`order = 2` is not a refinement, it is load-bearing.** Dropping to first order
flips the sign and gives an error seven times larger than σ = 16. The outgoing
invariant's extrapolation is what lets the front's structure leave; without it
the ghost clamps the front and the domain drains.

**`pin_farfield` is the interesting compromise.** It is nearly non-reflecting
(0.015%) and still anchors to within 748 — better than σ = 4 at a fraction of
the reflection. Its known cost is that it anchors to p_target + ρc·u_out rather
than p_target, which is exactly the residual seen here.

## Why β_s does nothing, and what actually goes wrong

The natural reading of the table is that the reaction source correction is
broken. It is not. Bisecting the configuration at σ = 1 (t = 10⁻⁵ s, domain-mean
pressure rise):

| configuration | mean Δp |
|---|---|
| uniform fresh gas, no flame in the domain | +25 |
| flame, `do_react = 0` (it just diffuses) | +417 |
| flame, reactions on, diffusion off | +793 |
| flame, reactions and diffusion | +1283 |
| flame, hard outflow | −23 |

**The error is already there with chemistry switched off.** Whatever is going
wrong is not the reaction source term, and β_s cannot fix it because it is not
what is being measured. Reactions roughly double an error that a purely
non-reacting density front had already created.

The mechanism is in the ghost pressure. At an outflow

```
p_g = ½ ρc (R₊,g − R₋,g)
R₊,g = R₊,N + ℓ δR₊                        (extrapolated, limited)
R₋,g = R₋,N + ℓ Δn L_in / ((c − u_out) ρc)  (relaxed)
```

so, per ghost layer ℓ,

```
p_g − p_N  =  ½ ρc ℓ δR₊  −  (ℓ σ / 2 n_x) (p_N − p_∞)
                 ^^^^                ^^^^
              extrapolation        anchoring
```

The anchoring increment carries a factor 1/n_x — with n_x = 96 and σ = 1 it is
0.5% of the pressure offset it is trying to remove. The extrapolation term
carries no such factor. Setting the two equal gives the equilibrium offset

```
Δp  ≈  ρ c L_ref (du_out/dn)|_boundary / σ
```

which is the whole story: a **normal velocity gradient at the boundary** biases
the ghost pressure, and the relaxation can only fight it in proportion to σ. A
flame crossing the outflow is a large `du_out/dn` — the gas accelerates from
45.6 to 123.8 cm/s across 0.07 cm — and it is dilatational, not acoustic, so
there is nothing wrong with the invariant algebra. LODI simply has no way to
tell the two apart.

Two checks that the formula is the right one rather than a plausible story:

* It predicts a σ⁻¹ trend. Measured: 2074 → 1450 → 218 for σ = 1 → 4 → 16.
* It predicts the offset is **grid-converged**, because δR₊ is a per-cell slope
  falling as 1/n_x while the anchoring carries 1/n_x explicitly. Measured
  (no chemistry, σ = 1, t = 10⁻⁵ s): n_x = 48 → +531, 96 → +409, 192 → +371.
  It converges to a finite value; it does not refine away.

## What to tell a user

This is the case that turns "place outflows away from flames" from received
wisdom into a number. If you must put one there:

* Raise σ. The inert default of 0.25 is the worst choice available. σ ≈ 10 is
  where the anchoring starts to win, and it costs you roughly 20% acoustic
  reflection, so you are trading one error for another and should decide which
  one your problem cares about.
* Or use `pin_farfield`, if a fixed offset of order ρc·u_out is acceptable and
  transparency is not.
* Set `bc_nscbc_beta = 0.5`. It is a free 10%.
* Leave `bc_nscbc_order = 2`.
* `bc_nscbc_beta_s` will not save you. The reaction source correction is real
  and correctly applied, but at a front-crossing outflow it is a second-order
  effect on top of a first-order problem.
* Run with `pelec.sum_interval > 0` at least once and read the fallback line.
  A `reaction source dropped` or `transverse dropped` count means a dial is
  silently doing nothing.

## Traps

**`pelec.allow_negative_energy = 0`** is harmless in an inert case and fatal
here — see `NSCBC-PMF/README.md`. The inputs file sets it explicitly to 1.

**`prob.probtype = 1`** retains a V-flame anchored on an inlet hot spot, from an
earlier attempt at this measurement. It compiles and initialises but has not
been used for a measurement; the wrinkled sheet superseded it because it needs
no anchor, no apex geometry and no steady state to wait for.

## Locating the reaction zone

Use a radical, not the temperature. The temperature midpoint sits in the preheat
zone, which is wide and diffusion-controlled; the radical peak sits in the
reaction zone, which is the thing that actually has to be inside the boundary
cells for `beta_s` to have anything to correct. `radicals.py` does both the
front-position tracking and the boundary-column profile.

Along the outflow column at t = 24 µs, H₂/air, σ = 1:

```
  y [cm]    T [K]         Y_H        Y_OH    q [erg/cm3/s]
  0.0031      299   1.494e-17   4.342e-14        9.387e+00
  0.1031      644   1.963e-07   4.007e-06        2.337e+08
  0.1281     1008   2.348e-05   2.532e-04        5.549e+09
  0.2531     1215   5.116e-05   9.520e-04        5.961e+09
  0.3781      301   6.946e-17   2.304e-13        8.336e+01
```

Twelve orders of magnitude in `Y_H` and nine in the heat release, along a single
column of boundary cells, with the two peaks at the two designed crossings. That
is the case doing what it claims.

The radical is also the better front tracker for the stability question below,
because the peak is 2–3 cells wide where the temperature midpoint is a ramp over
ten.

## Is the measurement contaminated by thermodiffusive instability?

A fair question of the H₂ case: at φ = 0.4, Le(H₂) ≈ 0.3, so an imposed wrinkle
is thermodiffusively unstable and would grow. It does not, because there is no
time for it to. Tracking the front by the H peak:

| t [s] | mean x_f [cm] | wrinkle half-amplitude [cm] |
|---|---|---|
| 0 | 0.65294 | 0.07990 |
| 9.78e−6 | 0.65262 | 0.07896 |
| 2.40e−5 | 0.65247 | 0.07809 |

The amplitude *decays* 2.3% over the run. The Darrieus–Landau e-folding time at
this wavenumber is ~9×10⁻³ s and the flame time δ/S_L is 3×10⁻³ s, against a
measurement window of 2.4×10⁻⁵ s — 130 to 390 times shorter than either. The
boundary equilibrates on the acoustic time, which is what makes the measurement
possible at all, and the flame is frozen on that time.

`NSCBC-FlameOutflow-DRM` runs the same problem with CH₄/air at φ = 0.75
(drm19, Le(CH₄) = 0.97, thermodiffusively neutral) and settles the question by
construction rather than by argument.

## A diffusive boundary condition (`bc_nscbc_extrap_temperature`)

The measurements above were made before the outflow had any diffusive treatment
at all. `Source/Diffusion.cpp` forms the conductive and species fluxes at a
physical boundary face from these ghost cells, and in the entropy closure the
ghost temperature is whatever the EOS returns from the extrapolated density and
the acoustically-set pressure — nothing chooses it with a heat flux in mind.

Check C8 in `Verification/NSCBC1D` puts a number on it. On a flame-like ramp
(2×10⁴ K/cm, 0.01 cm cells) the ghost overstates the face temperature gradient
by **40%**: T_ghost = 1680 K where the interior ramp continues to 1600 K.

`pelec.bc_nscbc_extrap_temperature = 1` closes the λ₀ family on temperature
instead — T is extrapolated on the same minmod slope as everything else and ρ
follows from the EOS — so the face gradient is the interior one, exactly. A
uniform state still comes back to 2×10⁻¹⁶, so nothing hyperbolic is given up.

Measured here, at t = 2.4×10⁻⁵ s, mean Δp against the shielded reference:

| outflow | entropy closure | temperature closure | change |
|---|---|---|---|
| σ = 1, β = 0.5, β_s = 0 | +2074 | +1894 | −9% |
| σ = 16, β = 0.5, β_s = 0 | +218 | **+67** | **−69%** |
| hard `p = p_amb` | −527 | — | — |

## The material continuation (`bc_nscbc_extrap_material`), and what it settles

The ghost-pressure bias of the previous section has a targeted fix:
continue the *material* part of the incoming invariant's slope into the
ghost, bounded through the entropy family so the incoming model never feeds
on the waves it launches (`Source/NSCBC.H`; gated by checks C9(a) and C11 in
`Verification/NSCBC1D`, where it removes the bias exactly and holds a
manufactured sustained front that the entropy closure walks away from).

Measured here, σ = 1, β = 0.5, β_s = 0, t = 2.4×10⁻⁵ s, mean Δp:

| closure | mean Δp | change |
|---|---|---|
| entropy | +2074 | — |
| + material | +1972 | −5% |
| + temperature | +1894 | −9% |
| + **material and temperature** | **+1200** | **−42%** |

Two conclusions, and the first one closes an open question.

**The extrapolation bias is real but is not what dominates this case.** The
fix removes the C9(a) mechanism *exactly* in the driver and buys 5% here, so
the σ-suppressed error of this configuration is mostly something else — by
the bisection above and the C9(b) result, the unmodelled diffusive and
reactive enthalpy deposition in the boundary cells. That also answers the
question the C9/C10 commits left open, **against** the flux-form
reformulation's premise: a flux-form NSCBC removes the same extrapolation
bias, so it too would buy ~5% here, at many times the cost.

**The two ghost closures compound.** Together they make the ghost a fully
consistent material continuation — p from the corrected acoustic pair, u
from the continued slope, T extrapolated, ρ from the EOS — and take 42%
where each alone takes single digits. Note two caveats: at U = 2 S_L the
quasi-steady bound in the material continuation underestimates the front's
true slope by 2× (the sheet drifts at S_L = U/2), so this configuration is
near its worst case; and the continuation is for *quasi-steady* structure —
see the exit test below before enabling it anywhere a front actually crosses.

## The exit test — `nscbc-flameexit.inp`

This case's original question was always "what happens when the flame
reaches the outflow", and the quasi-frozen configuration answers it only for
a front that *sits* there. `nscbc-flameexit.inp` makes the sheet actually
leave: U = 40 S_L, the sheet initialised fully inside, and the whole passage
— approach, crossing, exit, emptied domain — inside 6×10⁻⁴ s. The wrinkle
keeps the parent's slope (51° tilt at the crossings) and the grid keeps the
parent's dx. Physics fixes the truth during transit (the front must advect
at U − S_L ≈ 889 cm/s with its wrinkle frozen; the decay times are 100× the
window), and after the exit the exact solution is known outright: a uniform
fresh stream at (U, p_amb, T_in). `exit_metrics.py` tracks both.

| outflow | transit: peak mean Δp | front | wrinkle | post-exit residual Δp, rms(u−U) |
|---|---|---|---|---|
| hard `p = p_amb` | −200 | on schedule | held | **+1.2, 0.01** |
| **σ = 16 + extrap_T** | +1800 | on schedule | held | **+8.5, 0.20** |
| σ = 16 + extrap_T + extrap_material | −3520 | on schedule | held | +8.5, 0.20 |
| σ = 0.25 (the shipped default) | +37000 → **NaN, crash** | pushed *backwards* | grows 40% | — |
| σ = 0.25 + extrap_T + extrap_material | **−340000**, flame quenched | expelled | destroyed | recovering on τ = 1/K |
| `pin_farfield` + extrap_T | +12000, persistent | 40% slow | distorted | +17500, 460 |

What the table says, in order of importance.

**The shipped default destroys the flame and then the run.** At σ = 0.25
the anchoring time τ = L/(σc) is comparable to the transit itself, so the
boundary integrates the crossing's flux imbalance into a pressure ramp that
pushes the front back upstream, pumps the wrinkle amplitude 40%, and NaNs at
t = 2×10⁻⁴ s. This is not a tuning nuance; it is the difference between the
run completing and not.

**A fast dilatational transit wants anchoring, not transparency.** The
crossing is not an acoustic event: what leaves is mass and enthalpy at
M ~ 0.01, whose own pressure field is p_amb to one part in 10⁴. The hard
Dirichlet — the *worst* acoustic boundary in this suite, R = 97% — handles
it almost perfectly, and among characteristic treatments the quality ranking
is exactly the anchoring-strength ranking. σ = 16 + `extrap_temperature`
exits the flame with front kinematics and wrinkle amplitude
indistinguishable from the hard-outflow truth, a transient of 0.18% of
ambient, and a domain that returns to the exact uniform stream to 8 parts
per million. The ~28% acoustic reflection that σ = 16 costs elsewhere never
appears, because there is nothing acoustic to reflect.

**`extrap_material` is for fronts that stay, not fronts that leave.** At
σ = 0.25 the continuation turns the crossing into a runaway: the entropy
bound is wide open while the front is in the boundary cells, the continued
slopes over-vent the domain, the weak relaxation cannot answer, and the
domain drains to −0.34 atm — cold enough to quench the flame before it
finishes leaving. At σ = 16 the same mechanism is bounded (the transit dip
doubles relative to extrap_T alone, and the run is otherwise clean), so a
strong σ makes the flag safe but not helpful here. Its measured value is the
quasi-frozen table above.

**`pin_farfield` anchors a through-flow duct to the wrong state.** It is
stable, but it holds p + ρcu rather than p, so the whole operating point
shifts: the stream runs 460 cm/s slow, the front crosses 40% late, and the
domain settles +17500 dyn/cm² off ambient and stays there. Use it for the
quiescent-reservoir case it was built for, not for an exit with a mean flow.

The two rows say something worth reading carefully. At σ = 1 the ghost-pressure
bias of the previous section still dominates, so removing the diffusive error
buys only 9%. At σ = 16 that bias is largely suppressed, the diffusive error is
most of what remains, and removing it takes the error down by more than a factor
of three.

**σ = 16 with the temperature closure is the first setting in this case that
beats a hard pressure outflow on anchoring** — 67 against 527, a factor of eight
— while remaining a characteristic boundary. The reflection cost of σ = 16 is
unchanged at ~28%, so the trade-off of the previous section has not gone away;
what has changed is that the anchoring side of it is now genuinely good rather
than merely less bad.

The default is off, because it changes every outflow result. Turn it on when a
thermal or compositional structure is anywhere near the boundary.
