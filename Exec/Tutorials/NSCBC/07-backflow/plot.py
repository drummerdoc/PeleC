"""Rung 7 figure: sustained backflow through an outflow, with and without
the reservoir closure."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import nscbc_tut as nt  # noqa: E402

AREA = 10.0 * 0.5
VARIANTS = [("frozen", "backflow_material = 0 (default)"),
            ("reservoir", "backflow_material = 1")]
# Reference: mean temperature [K] and inflow speed [cm/s] at t = 30 ms.
REFERENCE = {"frozen": (600.0, 2864.0), "reservoir": (300.0, 1769.0)}

fig, (ax_t, ax_x) = nt.plt.subplots(1, 2, figsize=(13, 4.4),
                                    layout="constrained")
for (tag, label), col in zip(VARIANTS, ["C3", "C0"]):
    log = nt.read_datlog(os.path.join("runs", tag, "datlog"))
    t = log["time"] * 1e3
    T = log["temp"] / AREA
    u = -log["xmom"] / log["mass"]  # speed of the stream drawn in at x = L
    ax_t.plot(t, T, color=col, label=label)
    ref = REFERENCE.get(tag)
    if ref:
        nt.check(f"{tag}: mean T at 30 ms [K]", T[-1], ref[0], "{:.0f}")
        nt.check(f"{tag}: inflow speed at 30 ms [cm/s]", u[-1], ref[1], "{:.0f}")
    else:
        print(f"  {tag}: mean T at 30 ms = {T[-1]:.0f} K, "
              f"inflow speed = {u[-1]:.0f} cm/s")
    files = nt.plotfiles(os.path.join("runs", tag))
    want = [5e-3, 10e-3, 15e-3] if tag == "reservoir" else [30e-3]
    for w in want:
        d = min((nt.load(f, ["Temp"]) for f in files),
                key=lambda dd: abs(dd["time"] - w))
        ax_x.plot(d["x"], nt.transverse_mean(d["Temp"]), color=col,
                  alpha=0.4 + 0.6 * w / 0.015 if tag == "reservoir" else 1.0,
                  label=f"{tag}, {1e3 * d['time']:.0f} ms")
ax_t.axhline(300.0, color="k", ls=":", lw=0.8)
ax_t.text(0.5, 310, "reservoir 300 K", fontsize=8)
ax_t.set_xlabel("t [ms]")
ax_t.set_ylabel("mean temperature in the duct [K]")
ax_t.set_title("a 600 K duct drawing 300 K air in through its outflow")
ax_t.legend(fontsize=8)
ax_x.set_xlabel("x [cm]   (gas enters at x = 10)")
ax_x.set_ylabel("temperature [K]")
ax_x.set_title("what is drawn in")
ax_x.legend(fontsize=8)
nt.save(fig, "backflow.png")
