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
