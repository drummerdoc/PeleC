"""Rung 1 figure: x-t diagrams of the pulse for both outflows."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import nscbc_tut as nt  # noqa: E402

AMP = 1.0e-3 * nt.P_AMB  # the inputs' prob.amp
CS = 34719.0  # ambient sound speed, cm/s (GammaLaw air, 300 K)
L = 10.0  # domain length, cm
REFERENCE = {"hard": 97.29, "nscbc": 0.81}  # reflection R in %

fig, axes = nt.plt.subplots(1, 3, figsize=(14, 4.2), layout="constrained")
inset = axes[2].inset_axes([0.45, 0.12, 0.5, 0.45])
for ax, name, title in [
    (axes[0], "hard", "hard pressure ghost"),
    (axes[1], "nscbc", r"characteristic outflow ($\sigma$ = 0.25)"),
]:
    files = nt.plotfiles(os.path.join("runs", name))
    rows, times = [], []
    for f in files:
        d = nt.load(f, ["pressure"])
        rows.append((nt.transverse_mean(d["pressure"]) - nt.P_AMB) / AMP)
        times.append(d["time"] * CS / L)
    x = d["x"]
    im = ax.pcolormesh(
        x, times, np.array(rows), cmap="RdBu_r", vmin=-1, vmax=1, shading="auto"
    )
    ax.set_title(title)
    ax.set_xlabel("x [cm]")
    ax.set_ylabel("t c / L")

    # Reflection: what is left in the upstream half once the pulse has had
    # time to leave (nscbc_tut.reflection).
    last = (nt.transverse_mean(d["pressure"]) - nt.P_AMB) / AMP
    R = 100.0 * nt.reflection(last, 1.0)
    nt.check(f"reflection R, {name}", R, REFERENCE[name], "{:.2f}", " %")
    axes[2].plot(x, last, label=f"{name}: R = {R:.2f} %")
    inset.plot(x, last)
fig.colorbar(im, ax=axes[:2], label=r"$p'$ / pulse amplitude")
axes[2].set_title(r"$p'$ at t c / L = 1.6")
axes[2].set_xlabel("x [cm]")
axes[2].legend(loc="upper left")
inset.set_ylim(-0.02, 0.02)
inset.set_title("zoom: \u00b12%", fontsize=8)
inset.tick_params(labelsize=7)
nt.save(fig, "outflow_pulse.png")
