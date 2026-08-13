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
| before the wave reaches the boundary (`ct` = 0.80 W) | 1.28% | 0.5% |
| hard boundary, crossing (`ct` = 1.21 W) | 0.87% | **21.4%** |
| NSCBC σ = 0.25, crossing (`ct` = 1.21 W) | **0.22%** | **0.89%** |

The amplitude column is the striking one. The hard boundary does not so much
bend the front as *chew* it: by the time the corner-directed arc is all that is
left, its strength varies by a factor of 1.2 around the arc, because the parts
that passed near the faces have already been corrupted by reflections. The
characteristic treatment holds it to 0.9% — within a factor of two of the
pre-contact discretisation baseline.

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

| configuration | max/δp | rms/δp | rms vs floor |
|---|---|---|---|
| **no boundary at all** (periodic; discretisation only) | 0.081 | **0.0047** | 1× |
| hard inflow + hard outflow | 0.151 | 0.0652 | 14× |
| hard inflow + NSCBC outflow | 0.170 | 0.0660 | 14× |
| NSCBC inflow + NSCBC outflow | 0.126 | **0.0468** | **10×** |
| … same, σ = 0 | 0.168 | 0.0702 | 15× |
| … same, `relax_u` = 0.2 | 0.170 | 0.0677 | 14× |

`δp` is the vortex's own pressure deficit (1.90e4 dyn/cm², ~1.9% of ambient).
The periodic row is obtained by comparing the final field with the initial field
shifted by `u·t`; it includes the shift-interpolation error, so the true floor
is lower still and the ratios below are conservative.

Three things follow, and they are the point of the test.

**The characteristic treatment helps, but only by about 1.4× in rms**, and the
residual stays an order of magnitude above the no-boundary floor. Compare the
NSCBC-Acoustic case, where a *planar* pulse at normal incidence reflects 120×
less. The gap between 120× and 1.4× is the whole content of "1-D-normal LODI".

**The dials are not what is limiting it.** Moving σ from 0.25 to 0, or `relax_u`
from 2 to 0.2 — a factor of ten each — changes the rms by less than the choice
of *which* boundary is characteristic. When a knob spans a decade and the answer
does not move, the error is not in that knob. It is in the terms that are not
there: the transverse contributions to the modelled incoming wave, Motheau's
β. A vortex crossing an outflow is the maximally transverse event, which is
exactly why this case was chosen.

**The outflow is nonetheless doing its job.** Switching only the outflow to the
characteristic treatment cuts the upstream residual (mean |δp| over x < 3 cm)
from 1648 to 342 dyn/cm² — a factor of 4.8. It is not sending a strong
reflected wave back upstream; it is failing to absorb the vortex cleanly, which
is a different fault with a different fix.

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
