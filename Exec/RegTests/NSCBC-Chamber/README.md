# NSCBC-Chamber — the T7 mini-SydGex capstone (Lesson 9)

A closed-end premixed chamber venting through its open end, run **two ways**
so that the difference of the two overpressure traces is the boundary
condition's contribution to an application-level answer:

* **`chamber-vent.inp`** — the characteristic boundary IS the vent plane
  (x = 1.2 cm), carrying the flame closures; the front crosses it live.
* **`chamber-plenum.inp`** — the same chamber, but the duct continues to
  x = 3.6 cm and the boundary sits two chamber lengths downstream. Identical
  boundary settings: *placement* is the only difference.

A half-disc of burnt products (r = 1 mm) on the closed-end wall ignites a
quiescent H₂/air φ = 0.4 charge (the shipped PMF profile; S_L = 22.8 cm/s).
The flame fingers down the duct, drives flow out of the open end, crosses the
vent plane, and the emptied chamber rings at its longitudinal modes — the
laminar-phase dilatation (Lesson 5), the crossing (Lesson 6), and the
boundary-set ring-down (Lessons 1 and 3) in one trace.

## Scaling, stated honestly

True SydGex is a 25 cm chamber and a multi-ms turbulent burn. This chamber is
**1.2 × 0.3 cm** (the Sydney L/H = 4) with a laminar burn, ~17 flame
thicknesses long and ~4 across, and the "plenum" is a straight duct
extension — the boundary-placement question T7 poses, without EB geometry.
Two follow-ups are deliberate omissions: a laterally-expanding plenum
(chamber inside a box) needs EB walls, and the baffle of the Sydney rig needs
either EB or an embedded thin wall. Scale up freely on a bigger machine: keep
dx ≤ 0.0125 cm (≥ 5.5 cells per thermal thickness; 0.00625 for the measured
11), grow `prob_hi`/`amr.n_cell` together, and grow `stop_time` with the
front-travel estimate below.

## Running

```sh
make -j TPL && make -j             # or the CMake build; LiDryer / Fuego / Simple
./PeleC2d.<comp>.ex chamber-vent.inp   prob.pmf_datafile=$PWD/LiDryer_H2_p1_phi0_4000tu0300.dat
./PeleC2d.<comp>.ex chamber-plenum.inp prob.pmf_datafile=$PWD/LiDryer_H2_p1_phi0_4000tu0300.dat
```

`prob.pmf_datafile` is CWD-relative — always pass an absolute path. For MPI
builds set `USE_MPI = TRUE` (GNUmake) and run under `mpiexec`; the case is
bit-identical across decompositions like every case in this suite.

Cost model (measured at the shipped scale, serial Apple M-class core):
dt ≈ 6.7×10⁻⁸ s, ~120k steps to `stop_time = 8 ms`, ≈ 2×10⁵ cell-steps/s —
about an hour serial for the vent variant, ~2.5× that for the plenum; minutes
on a many-core node. Front-travel estimate for sizing `stop_time`: the front
leaves the kernel at ~σ_exp·S_L ≈ 100 cm/s and accelerates with its area
ratio; allow (L_chamber − r_kernel)/(150 cm/s) plus ~1 ms of ring-down.

## Measuring

```sh
./chamber_metrics.py runs/vent runs/plenum --xvent 1.2
```

prints, per plotfile: closed-end overpressure, chamber-mean overpressure,
burnt volume fraction, and the volumetric flux through the vent plane. The T7
QoIs are built from the **difference of the two runs' traces**: peak-aligned
overpressure (peaks aligned, not clocks — ignition transients differ),
V̇_comb − V̇_vent whose zero crossing marks the peak, and the post-peak
ringing frequency and decay rate, which is where the boundary's reflection
coefficient is directly visible (Lesson 9's punchline: repeat the vent
variant at σ = 16 and watch the ring-down change while the peak barely
moves).

## Boundary settings

The vent carries the flame-crossing recipe measured in `NSCBC-FlameOutflow`:
`extrap_temperature = 1`, `beta_s = 0`, `beta = 0.5`, σ = 0.25 —
acoustically open, so the ring-down shows the chamber rather than the
boundary. `extrap_material` stays off (a front crosses). The plenum variant
uses identical settings so placement is the only difference; with the
boundary two chamber lengths away, its settings are nearly irrelevant, which
is the point.

## Status

Production-complete, and the case has now killed THREE defective reversal
closures on its way here — it is the only configuration in the suite whose
boundary genuinely reverses (the vent breathes), so it finds what nothing
else can. The history, in order: the original hard ambient pin NaN'd the
first production run at 4.08 ms (backflow −144 → −2628 cm/s, a 138 K cell);
the soft pressure-only replacement survived but manufactured a spurious
0.32 atm chamber spike at 4.0–4.6 ms by drawing inflow against a
+0.3 atm adverse gradient (dropped `S_p`/`dR₊`/`T_in`, frozen ghost
velocity); and the first unification of the closure NaN'd again at 4.05 ms
by extrapolating material content into the backflow (the cold runaway,
385 → 89 K in 42 µs). Full forensics:
`Docs/NSCBC-reversal-branch-defect.md`; static gate: driver C13 a/b/c. The
shipped closure — unified acoustics, `w_mat`-upwinded material slopes —
completes all three variants to 8 ms.

The T7 comparison, measured with `chamber_qoi.py` on the 2026-08-25
production set (identical boundary settings, placement the only difference):

| variant | V̇ plateau [cm²/s] | end-of-burn peak [dyn/cm²] | f_ring [kHz] | decay λ [1/s] |
|---|---|---|---|---|
| vent, σ = 0.25 | 72.2 | **+37,473** | 3.33 | 674 |
| vent, σ = 16 | 71.6 | **+6,875** | 14.16 | 504 |
| plenum (same window) | 72.1 | +1,951 | 0.51 | 184 |

Three phases, three answers. **Buildup**: placement costs nothing — both
variants drive identical venting (V̇ 72.2 vs 72.1 cm²/s); the plenum's
standing +~3×10³ dyn/cm² offset is the inertia of the 2.4 cm duct column it
pushes, physics any boundary would see. **End of burn** (t ≈ 4.0–4.6 ms):
the last ~3% of charge burns while the vent flow chokes — a real
constant-volume pressurisation — and how much of it the boundary lets
through is where placement bites: +37k at σ = 0.25 against the plenum's
+2k, halved and halved again by σ = 16. This is the honest remaining
boundary cost, down 8.5× from the defective closure's 321k, and with NO
counter-gradient inflow (V̇ dips to −9 cm²/s where the defective closure
drew −374 against the gradient). **Ring-down** is where σ is directly
legible: at σ = 0.25 the "mode" is the relaxation itself
(3.33 kHz ≈ K/2π = σc/2πL), while σ = 16 stiffens the end toward a real
acoustic termination and the chamber rings near its burnt-gas quarter-wave
(14.2 measured vs ~18.8 kHz ideal) — the Lesson-9 punchline, measured. The
σ = 16 run also takes brief genuine backflow breaths (V̇ to −110 cm²/s)
through the fixed closure without incident.

Not registered in `Tests/CMakeLists.txt`: at ~1 CPU-hour per variant it is
an application capstone, not a CI gate.
