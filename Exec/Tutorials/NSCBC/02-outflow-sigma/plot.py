"""Rung 2 figure: reflection and mean-pressure recovery against sigma."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import nscbc_tut as nt  # noqa: E402

AMP = 1.0e-3 * nt.P_AMB
CS, L, AREA = 34719.0, 10.0, 10.0 * 0.5
VARIANTS = [
    ("s0", r"$\sigma$ = 0"),
    ("s0.25", r"$\sigma$ = 0.25"),
    ("s1", r"$\sigma$ = 1"),
    ("s4", r"$\sigma$ = 4"),
    ("s16", r"$\sigma$ = 16"),
    ("s16_ndnr", r"$\sigma$ = 16 + NDNR"),
]
# Reference numbers: reflection R [%], and the largest |mean-pressure
# deviation| after 20 L/c as a fraction of the initial overpressure.
REFERENCE = {
    "s0": (0.00098, 1.0),
    "s0.25": (0.813, 0.00063),
    "s1": (2.89, 0.0016),
    "s4": (7.57, 0.21),
    "s16": (27.5, 0.74),
    "s16_ndnr": (0.915, 9.9e-05),
}

fig, (ax_r, ax_m) = nt.plt.subplots(
    1, 2, figsize=(13, 4.4), layout="constrained", width_ratios=[1, 1.6]
)
colors = nt.plt.cm.viridis(np.linspace(0, 0.9, len(VARIANTS)))
Rs, late = [], []
for (tag, label), col in zip(VARIANTS, colors):
    d = nt.load(nt.plotfiles(os.path.join("runs", "pulse_" + tag))[-1],
                ["pressure"])
    p = (nt.transverse_mean(d["pressure"]) - nt.P_AMB) / AMP
    Rs.append(100.0 * nt.reflection(p, 1.0))

    log = nt.read_datlog(os.path.join("runs", "ring_" + tag, "datlog"))
    t, pm = nt.mean_pressure(log, AREA)
    dev = (pm - nt.P_AMB) / AMP  # the initial overpressure is 1 in these units
    tl = t * CS / L
    ax_m.plot(tl, dev, color=col, label=label,
              lw=2.2 if tag == "s16_ndnr" else 1.2)
    late.append(np.abs(dev[tl > 20.0]).max())

for (tag, label), R, lt in zip(VARIANTS, Rs, late):
    label = label.replace(r"$\sigma$", "sigma")
    ref = REFERENCE.get(tag)
    if ref:
        nt.check(f"{label}: reflection R", R, ref[0], "{:.3g}", " %")
        nt.check(f"{label}: max |mean dev| after 20 L/c", lt, ref[1], "{:.2g}")
    else:
        print(f"  {label:<24s} R = {R:.3g} %   max |mean dev| after 20 L/c = {lt:.2g}")

ax_r.bar(range(len(Rs)), Rs, color=colors)
ax_r.set_yscale("log")
ax_r.set_xticks(range(len(Rs)), [v[1] for v in VARIANTS], rotation=30,
                ha="right")
ax_r.set_ylabel("pulse reflection R [%]")
ax_r.set_title("(A) letting a wave out")
for i, R in enumerate(Rs):
    ax_r.text(i, R * 1.15, f"{R:.2g}", ha="center", fontsize=8)

ax_m.axhline(0.0, color="k", lw=0.6)
ax_m.set_xlabel("t c / L")
ax_m.set_ylabel(r"$(\langle p \rangle - p_{target})$ / initial overpressure")
ax_m.set_title("(B) holding the mean: a vented duct returning to the target")
ax_m.set_ylim(-1.1, 1.1)
ax_m.legend(fontsize=8, ncol=2)
nt.save(fig, "outflow_sigma.png")
