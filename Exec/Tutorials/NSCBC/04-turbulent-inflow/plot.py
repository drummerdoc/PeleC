"""Rung 4 figure: how much of a turbulent inflow gets in, against x."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import nscbc_tut as nt  # noqa: E402

SCALE = 100.0  # turbinflow.xlo.turb_scale_vel [cm/s]
T0 = 1.0e-3  # statistics from here on (the inlet region has filled by then)
VARIANTS = [
    ("defaults", "relax_u = 2, relax_t = 0.2 (defaults)"),
    ("u20", "relax_u = 20, relax_t = 0.2"),
    ("t20", "relax_u = 2, relax_t = 20"),
    ("ut20", "relax_u = 20, relax_t = 20"),
    ("ut20_nri", "relax_u = 20, relax_t = 20, NRI"),
]
# Reference: transmitted fraction at the first cell, (normal, tangential).
REFERENCE = {
    "defaults": (0.53, 0.29),
    "u20": (0.76, 0.32),
    "t20": (1.16, 0.90),
    "ut20": (0.93, 0.90),
    "ut20_nri": (1.40, 0.90),
}


def hit_rms(path):
    d = np.loadtxt(path, delimiter=",", skiprows=1)
    return np.sqrt(np.mean(d[:, 3:6] ** 2, axis=0))


def intensities(run):
    """u', and the tangential v', w' combined, against x: rms over the
    cross-section and over t >= T0."""
    U, V, W = [], [], []
    for f in nt.plotfiles(os.path.join("runs", run)):
        d = nt.load(f, ["x_velocity", "y_velocity", "z_velocity"])
        if d["time"] >= T0:
            U.append(d["x_velocity"])
            V.append(d["y_velocity"])
            W.append(d["z_velocity"])
    U, V, W = (np.array(a) for a in (U, V, W))  # (nt, nx, ny, nz)

    def rms(a):
        return np.sqrt(((a - a.mean(axis=(0, 2, 3), keepdims=True)) ** 2)
                       .mean(axis=(0, 2, 3)))

    return d["x"], rms(U), np.sqrt(0.5 * (rms(V) ** 2 + rms(W) ** 2))


r = SCALE * hit_rms(os.path.join("runs", "hit.dat"))
u_pre, t_pre = r[0], np.sqrt(0.5 * (r[1] ** 2 + r[2] ** 2))
print(f"prescribed: u' = {u_pre:.1f} cm/s, tangential = {t_pre:.1f} cm/s")

fig, (ax_n, ax_t) = nt.plt.subplots(1, 2, figsize=(13, 4.4),
                                    layout="constrained", sharey=True)
cols = nt.plt.cm.viridis(np.linspace(0, 0.9, len(VARIANTS)))
for (tag, label), col in zip(VARIANTS, cols):
    x, un, ut = intensities(tag)
    ls = "--" if tag.endswith("nri") else "-"
    ax_n.plot(x, un / u_pre, ls, color=col, label=label)
    ax_t.plot(x, ut / t_pre, ls, color=col, label=label)
    fn, ft = un[0] / u_pre, ut[0] / t_pre
    ref = REFERENCE.get(tag)
    if ref:
        nt.check(f"{tag}: normal fraction at inlet", fn, ref[0], "{:.2f}")
        nt.check(f"{tag}: tangential fraction at inlet", ft, ref[1], "{:.2f}")
    else:
        print(f"  {tag}: at the inlet, normal {fn:.2f}, tangential {ft:.2f}")
for ax, what in [(ax_n, "normal (streamwise) u'"),
                 (ax_t, "tangential v', w'")]:
    ax.axhline(1.0, color="k", lw=0.6)
    ax.set_xlabel("x [cm]  (inlet at 0)")
    ax.set_title(f"{what} / prescribed")
    ax.set_xlim(0, 2.0)
    ax.set_ylim(0, 1.5)
ax_n.set_ylabel("realised / prescribed intensity")
ax_t.legend(fontsize=8, loc="upper right")
nt.save(fig, "turbulent_inflow.png")
