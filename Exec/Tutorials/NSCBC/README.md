# Tutorial: characteristic (NSCBC) boundary conditions in PeleC

A ladder of small, preconfigured runs that shows what PeleC's characteristic
boundary condition does, how each of its settings changes the answer, and
where it stops being the right tool. Every rung builds an existing case,
runs a handful of variants, and plots one figure; the plot script prints
your numbers next to ours, so you know when you have reproduced it.

By the end you should be able to choose `σ`, `relax_u`, `relax_t` and `β`
for a new problem, recognise the signs of a boundary that is in the wrong
place, know which structures may cross an outflow and which may not, and
know when the right answer is a longer domain rather than a better setting.

The reference for every key is the *Boundary conditions* chapter of the
PeleC documentation; this tutorial is the hands-on companion to it.

## The ladder

| Rung | Question | Case | Cost | Status |
|---|---|---|---|---|
| [1](01-outflow-pulse/README.md) | Why a characteristic outflow? | NSCBC-Acoustic 2-D | 1 min | ready |
| [2](02-outflow-sigma/README.md) | σ: letting waves out against holding the mean | NSCBC-Acoustic 2-D | 5 min | ready |
| 3 | Subsonic inflow: `relax_u`, forcing (`dudt`) and the NRI register | NSCBC-Acoustic duct | — | planned |
| 4 | A turbulent inflow through the characteristic inlet | NSCBC-Acoustic 3-D duct + TurbInflow | — | planned |
| 5 | Oblique waves and vortices: the transverse weight β | NSCBC-COVO | — | planned |
| 6 | What may cross an outflow: vortices, entropy and composition fronts, flames; the flame closures | NSCBC-COVO, NSCBC-FlameOutflow | — | planned |
| 7 | Flow reversal at a vent and sustained backflow | NSCBC-Chamber | — | planned |
| 8 | When the answer is a longer domain | NSCBC-FlameOutflow, NSCBC-Chamber | — | planned |

## Running a rung

```sh
cd 01-outflow-pulse
./run.sh            # builds the case it needs, runs every variant into runs/
python plot.py      # writes the figure here and prints measured vs reference
```

`run.sh` skips any variant that already finished; delete its directory under
`runs/` to rerun it. Environment knobs (see `common.sh`): `USE_MPI`
(default `TRUE`), `NPROCS` (default 4), `MAKE_JOBS`, `MPIEXEC`, `PELE_HOME`.

Prerequisites: a working PeleC build environment (the cases build with
GNU make), MPI, and Python 3 with `numpy`, `matplotlib` and `yt`
(`pip install numpy matplotlib yt`).

## Conventions

* Lengths in cm, pressures in dyn/cm², the CGS units of PeleC. The target
  pressure is 1 atm = 1.01325 × 10⁶ dyn/cm².
* Times are given in acoustic transits, `t c / L`, with `L` the domain length
  normal to the boundary under study.
* The reflection coefficient `R` is the largest residual left upstream after a
  pulse has had time to leave, measured from the median of that region and
  divided by the incident amplitude (`nscbc_tut.reflection`).
* Reference figures live in `Docs/sphinx/images/nscbc_tutorial/`.
