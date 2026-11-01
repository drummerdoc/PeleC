"""Rung 5 figure: the transverse weight beta, on oblique waves and on a
vortex."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import nscbc_tut as nt  # noqa: E402

PULSE = [("hard", "hard"), ("b1.0", r"$\beta$ = 1"), ("b0.8", r"$\beta$ = 0.8"),
         ("b0.5", r"$\beta$ = 0.5"), ("b0.0", r"$\beta$ = 0")]
COVO = [("hard", "hard"), ("b1.0", r"$\beta$ = 1"), ("b0.5", r"$\beta$ = 0.5"),
        ("b0.2", r"$\beta$ = 0.2"), ("b0.0", r"$\beta$ = 0")]
# Reference: pulse front amplitude spread [%]; vortex residual / floor.
REFERENCE = {
    "pulse_hard": 21.4, "pulse_b1.0": 0.896, "pulse_b0.8": 0.0661,
    "pulse_b0.5": 1.08, "pulse_b0.0": 2.95,
    "covo_hard": 14.3, "covo_b1.0": 10.3, "covo_b0.5": 3.86,
    "covo_b0.2": 7.14, "covo_b0.0": 772.0,
}


def plain(label):
    return label.replace(r"$\beta$", "beta")


def crossing(run):
    """The plotfile at which only the corner-directed arc of the front is
    still inside the domain (front radius 5.6-6.4 cm), and its sound speed."""
    files = nt.plotfiles(os.path.join("runs", run))
    d0 = nt.load(files[0], ["pressure", "density"])
    c = np.sqrt(1.4 * d0["pressure"][0, 0] / d0["density"][0, 0])
    for f in files:
        d = nt.load(f, ["pressure"])
        if 5.6 < c * d["time"] < 6.4:
            return d, c
    raise SystemExit(f"{run}: no plotfile at the crossing time")


def vortex_residual(run):
    files = nt.plotfiles(os.path.join("runs", run))
    a0 = nt.load(files[0], ["pressure", "x_velocity"])
    a1 = nt.load(files[-1], ["pressure"])
    dp = np.abs(a0["pressure"] - nt.P_AMB).max()  # the vortex's own deficit
    if run != "covo_periodic":
        return np.sqrt(((a1["pressure"] - nt.P_AMB) ** 2).mean()) / dp
    # No boundary: compare with the initial field carried along by u t.
    x = a0["x"]
    dx, L = x[1] - x[0], x[-1] - x[0] + (x[1] - x[0])
    u = np.median(a0["x_velocity"])
    n = ((u * a1["time"]) % L) / dx
    i0, f = int(n), n - int(n)
    p0 = a0["pressure"]
    ref = (1 - f) * np.roll(p0, i0, axis=0) + f * np.roll(p0, i0 + 1, axis=0)
    return np.sqrt(((a1["pressure"] - ref) ** 2).mean()) / dp


fig, ax = nt.plt.subplot_mosaic(
    [["ph", "p8", "A"], ["v1", "v5", "B"]], figsize=(14, 8.4),
    layout="constrained", width_ratios=[1, 1, 1.15])

spreads = []
for tag, label in PULSE:
    d, c = crossing("pulse_" + tag)
    dp = d["pressure"] - nt.P_AMB
    s, nrays = nt.front_amplitude_spread(dp, d["x"], d["y"], c * d["time"])
    spreads.append(100.0 * s)
    if tag in ("hard", "b0.8"):
        a = ax[{"hard": "ph", "b0.8": "p8"}[tag]]
        lim = 1.0e-4 * nt.P_AMB
        im_a = a.pcolormesh(d["x"], d["y"], dp.T, cmap="RdBu_r", vmin=-lim,
                            vmax=lim, shading="auto")
        a.set_aspect("equal")
        a.set_title(f"(A) pulse front, {label}")
        a.set_xticks([])
        a.set_yticks([])

print("(A) circular pulse: spread of the front amplitude around the arc")
for (tag, label), s in zip(PULSE, spreads):
    ref = REFERENCE.get("pulse_" + tag)
    if ref is not None:
        nt.check(f"{plain(label)}: amplitude spread", s, ref, "{:.3g}", " %")
    else:
        print(f"  {plain(label)}: {s:.3g} %")
ax["A"].bar(range(len(spreads)), spreads,
            color=["grey"] + list(nt.plt.cm.viridis(np.linspace(0, 0.9, 4))))
ax["A"].set_yscale("log")
ax["A"].set_xticks(range(len(PULSE)), [p[1] for p in PULSE])
ax["A"].set_ylabel("front amplitude spread [%]")
ax["A"].set_title("front strength around the arc")
for i, s in enumerate(spreads):
    ax["A"].text(i, s * 1.15, f"{s:.2g}", ha="center", fontsize=8)

# What the vortex leaves behind: the pressure field once it has gone.
lim = None
for tag, key in [("b1.0", "v1"), ("b0.5", "v5")]:
    d = nt.load(nt.plotfiles(os.path.join("runs", "covo_" + tag))[-1],
                ["pressure"])
    dp = d["pressure"] - nt.P_AMB
    lim = lim or np.abs(dp).max()
    im_b = ax[key].pcolormesh(d["x"], d["y"], dp.T, cmap="RdBu_r", vmin=-lim,
                              vmax=lim, shading="auto")
    ax[key].set_aspect("equal")
    ax[key].set_xticks([])
    ax[key].set_yticks([])
    lab = dict(COVO)[tag]
    ax[key].set_title(f"(B) left behind by the vortex, {lab}")

floor = vortex_residual("covo_periodic")
res = [vortex_residual("covo_" + tag) / floor for tag, _ in COVO]
print(f"(B) vortex: residual after it has left, over the no-boundary floor "
      f"({floor:.4f} of the vortex's pressure deficit)")
for (tag, label), r in zip(COVO, res):
    ref = REFERENCE.get("covo_" + tag)
    if ref is not None:
        nt.check(f"{plain(label)}: residual / floor", r, ref, "{:.3g}", "x")
    else:
        print(f"  {plain(label)}: {r:.3g}x")
ax["B"].bar(range(len(res)), res,
            color=["grey"] + list(nt.plt.cm.viridis(np.linspace(0, 0.9, 4))))
ax["B"].set_yscale("log")
ax["B"].axhline(1.0, color="k", ls="--", lw=0.8)
ax["B"].text(len(res) - 0.5, 1.1, "no boundary", ha="right", fontsize=8)
ax["B"].set_xticks(range(len(COVO)), [p[1] for p in COVO])
ax["B"].set_ylabel("residual / no-boundary floor")
ax["B"].set_title("residual after the vortex has left")
for i, r in enumerate(res):
    ax["B"].text(i, r * 1.15, f"{r:.2g}x" if r < 100 else f"{r:.0f}x",
                 ha="center", fontsize=8)
fig.colorbar(im_a, ax=[ax["ph"], ax["p8"]], shrink=0.8,
             label="p − p_amb [dyn/cm²] (one scale)")
fig.colorbar(im_b, ax=[ax["v1"], ax["v5"]], shrink=0.8,
             label="p − p_amb [dyn/cm²] (one scale)")
nt.save(fig, "transverse.png")
