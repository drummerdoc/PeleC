# NSCBC-COVO — two multi-dimensional tests of the characteristic boundary

One executable, two problems, selected with `prob.probtype`. Both are here to
answer questions the 1-D tests cannot: what happens when a wave meets the
boundary *obliquely*, and what happens when the thing crossing the boundary is
not a wave at all.

```sh
./PeleC-NSCBC-COVO nscbc-covo.inp        # probtype 0, convected vortex
./PeleC-NSCBC-COVO nscbc-pulse2d.inp     # probtype 1, circular pulse
# add pelec.bc_nscbc=0 to either for the hard-boundary reference
```

---

## probtype 1 — circular Gaussian pulse

After Motheau et al. (2017). A quiescent square, all four faces non-reflecting,
a Gaussian pressure bump at the centre. The outgoing circular wave reaches every
face and every corner at once. The exact wavefront is a circle of radius `ct`,
so **any departure from circularity is boundary error** — and unlike an
amplitude metric it cannot be hidden by a uniform pressure offset.

Two numbers, both measured over the polar rays whose front is still inside the
domain (i.e. the corner-directed ones, once the face-normal parts have exited):
the spread in front **radius**, and the spread in front **amplitude**, around
the arc.

| | radius spread | amplitude spread |
|---|---|---|
| hard boundary | 0.98% | **21.4%** |
| NSCBC, β = 1 (transverse off) | 0.47% | 0.90% |
| NSCBC, β = 0.8 | 0.40% | **0.066%** |
| NSCBC, β = 0.5 | 0.40% | 1.08% |
| NSCBC, β = 0 (full transverse) | 0.36% | 2.95% |

The amplitude column is the striking one. The hard boundary does not so much
bend the front as *chew* it: by the time the corner-directed arc is all that is
left, its strength varies by a factor of 1.2 around the arc, because the parts
that passed near the faces have already been corrupted by reflections. The
characteristic treatment holds it to 0.9% without transverse terms and to
**0.066% with them at β = 0.8** — a further factor of 14.

The radius column saturates: 0.4% is roughly one cell, so all the β < 1 rows
are "as circular as this mesh can represent". Amplitude uniformity is the
discriminating measure here, not radius.

The picture says it faster than the table. With the hard boundary the
rarefaction ring squares off, kinks at the corners, and grows four bright
reflection lobes; with the NSCBC it stays round and simply leaves.

`W` is the domain half-width. The 1.28% "before contact" figure is the
Cartesian-mesh discretisation error for a circular front and is the floor for
this metric, not a boundary effect — note that it is measured over all 720 rays
while the crossing-time rows use only the ~52 corner-directed ones, so the rows
are comparable to each other but not to the first row.

---

## probtype 0 — convected isentropic vortex (COVO)

A Yee/Shu vortex carried at M = 0.2 from a subsonic inflow at x-lo to a subsonic
outflow at x-hi, y periodic. A vortex is a pure vorticity/entropy structure: it
carries essentially no acoustic content, so **an ideal boundary would let it
leave in silence** and any pressure signal left in the domain afterwards was
manufactured by the boundary. It is also the only test in this suite that
exercises the **inflow** path end to end.

This is the test that shows the current implementation's limit honestly.

| configuration | rms/δp | rms vs floor |
|---|---|---|
| **no boundary at all** (periodic; discretisation only) | **0.0047** | 1× |
| hard inflow + hard outflow | 0.0652 | 14.0× |
| NSCBC, β = 1 (transverse off) | 0.0468 | 10.0× |
| NSCBC, β = 0.8 | 0.0341 | 7.3× |
| NSCBC, **β = 0.5** | **0.0176** | **3.8×** |
| NSCBC, β = 0.35 | 0.0187 | 4.0× |
| NSCBC, β = 0.2 (= local Mach) | 0.0754 | 16.1× |
| NSCBC, β = 0 (full transverse) | 0.6184 | 132× |

Dial sensitivity at β = 1, for reference: σ = 0 gives 15×, `relax_u` = 0.2 gives
14×, and hard-inflow + NSCBC-outflow gives 14× — i.e. none of them matters
compared to β.

`δp` is the vortex's own pressure deficit (1.90e4 dyn/cm², ~1.9% of ambient).
The periodic row is obtained by comparing the final field with the initial field
shifted by `u·t`; it includes the shift-interpolation error, so the true floor
is lower still and the ratios below are conservative.

Three things follow, and they are the point of the test.

**Without transverse terms the treatment helps by only 1.4×**, and the residual
stays an order of magnitude above the no-boundary floor — against 120× for a
*planar* pulse at normal incidence in NSCBC-Acoustic. That gap is the whole
content of "1-D-normal LODI", and none of σ, `relax_u` or the choice of which
boundary is characteristic closes it. When knobs span a decade and the answer
does not move, the error is not in those knobs.

**Turning the transverse terms on closes most of it.** β = 0.5 brings the
residual to 3.8× the floor: 2.7× better than β = 1 and 3.7× better than a hard
boundary. That is the term the dial sweep was pointing at.

**But the response to β is sharply asymmetric.** Too little correction costs a
factor of a few; too much is catastrophic — β = 0.2 is worse than no transverse
terms at all and β = 0 is close to unstable, at 132× the floor. That asymmetry
is why the shipped default is β = 1 rather than the optimum. A wrong β is far
more dangerous than an absent one.

**The optimum is predictable, not empirical.** For a plane wave meeting the
boundary at angle θ, the correction that exactly cancels the obliqueness error
in the 1-D characteristic decomposition is `(1−β) = cos θ / (1 + cos θ)`, so

    β_opt = 1 − cos θ / (1 + cos θ)

which is 0.5 at normal incidence and rises toward 1 at grazing incidence (0.67
at 60°, 0.79 at 75°). Both measured optima land on that curve: the vortex — a
broad, near-normal-incidence structure — optimises at **0.5**, and the circular
pulse, whose surviving arc crossed the faces near-tangentially, optimises at
**0.8**. Two independent problems, one formula, no fitting.

**The outflow was never the reflecting part.** Switching only the outflow to
the characteristic treatment cuts the upstream residual (mean |δp| over
x < 3 cm) from 1648 to 342 dyn/cm² — a factor of 4.8. It was not sending a
reflected wave back upstream; it was failing to absorb the vortex cleanly,
which is a different fault, and β is its fix.

Use `prob.nscbc_inflow = 0` to leave the inflow to the ordinary `bcnormal` and
isolate the outflow, which is how the third row above was obtained.

---

## Verification note

The control that makes all of this trustworthy: with `prob.vortex_strength = 0`
— uniform flow, nothing crossing anything — both boundary conditions produce
**exactly zero** pressure disturbance to machine precision, over 3000 steps.
That is the AMReX-side confirmation of check C1 in `Verification/NSCBC1D`, and
it means every number above is a response to the flow rather than noise the
boundary generates on its own.

## Analysis tooling

`Verification/NSCBC2D/fielddump.cpp` flattens one variable of a single-level 2-D
plotfile onto a regular array (`fextract` gives 1-D slices only, and the
quantity of interest here is a shape). `Verification/NSCBC2D/metrics.py`
computes the circularity and residual metrics from those dumps.
