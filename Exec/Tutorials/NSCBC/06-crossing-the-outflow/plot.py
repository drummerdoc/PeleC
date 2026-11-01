"""Rung 6 figure: a hot front, inert and burning, leaving through an outflow."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import nscbc_tut as nt  # noqa: E402

U = 912.1643066  # stream velocity, 40 S_L (printed by the case at start-up)
S_L = 22.80410767  # laminar flame speed
X_OUT = 0.6  # the outflow plane [cm]

# (run, label, colour, line style)
INERT = [("inert_hard", "hard p = p_amb", "k", "-"),
         ("inert_defaults", "σ = 0.25", "C3", "-"),
         ("inert_closures", "σ = 0.25, eT = 1", "C0", "-"),
         ("inert_s16_defaults", "σ = 16", "C1", "-"),
         ("inert_s16_closures", "σ = 16, eT = 1", "C2", "-")]
FLAME = [("flame_hard", "hard p = p_amb", "k", "-"),
         ("flame_defaults", "σ = 0.25 (defaults)", "C3", "-"),
         ("flame_closures", "σ = 0.25, eT = 1, β_s = 0", "C0", "-"),
         ("flame_s16_defaults", "σ = 16", "C1", "-"),
         ("flame_s16_closures", "σ = 16, eT = 1, β_s = 0", "C2", "-"),
         ("flame_s16_closures_ndnr", "σ = 16, eT = 1, β_s = 0, NDNR", "C4",
          "-")]
# Reference: largest domain-mean pressure excursion [dyn/cm^2], the residual
# at 0.6 ms [dyn/cm^2], and the time the front has fully left [ms] (None:
# it never does).
REFERENCE = {
    "inert_hard": (-192, 1.2, 0.21),
    "inert_defaults": (20067, 1738.2, 0.29),
    "inert_closures": (-77562, -1895.3, 0.15),
    "inert_s16_defaults": (451, 8.6, 0.21),
    "inert_s16_closures": (-2997, 8.6, 0.21),
    "flame_hard": (-203, 1.2, 0.23),
    "flame_defaults": (58615, 55581.0, None),
    "flame_closures": (-8818, -250.3, 0.23),
    "flame_s16_defaults": (1561, 8.5, 0.23),
    "flame_s16_closures": (-432, 8.5, 0.22),
    "flame_s16_closures_ndnr": (-1062, -13.7, 0.22),
}


def crossings(T, x, T_mid):
    """Per transverse row, where T last rises through T_mid (NaN if not)."""
    xf = np.full(T.shape[1], np.nan)
    for j in range(T.shape[1]):
        r = T[:, j]
        i = np.where((r[:-1] < T_mid) & (r[1:] >= T_mid))[0]
        if len(i):
            i = i[-1]
            xf[j] = x[i] + (T_mid - r[i]) / (r[i + 1] - r[i]) * (x[i + 1] - x[i])
    return xf


def history(run, speed):
    """Mean pressure excursion, the front's lead over its schedule (mean over
    the rows that still hold it), and the exit time of one run."""
    files = nt.plotfiles(os.path.join("runs", run))
    if not files or not os.path.exists(os.path.join("runs", run, "done")):
        return None
    t, dp, lead = [], [], []
    xf0 = T_mid = None
    t_exit = None
    for f in files:
        d = nt.load(f, ["pressure", "Temp"])
        if xf0 is None:
            T_mid = 0.5 * (d["Temp"].min() + d["Temp"].max())
            xf0 = crossings(d["Temp"], d["x"], T_mid)
        xf = crossings(d["Temp"], d["x"], T_mid)
        t.append(d["time"])
        dp.append(d["pressure"].mean() - nt.P_AMB)
        held = np.isfinite(xf)
        lead.append(np.mean(xf[held] - xf0[held] - speed * d["time"])
                    if held.any() else np.nan)
        if not held.any() and t_exit is None:
            t_exit = d["time"]
    return np.array(t), np.array(dp), np.array(lead), t_exit


def report(run, h):
    t, dp, _, t_exit = h
    k = np.argmax(np.abs(dp))
    meas = (dp[k], dp[-1], None if t_exit is None else 1e3 * t_exit)
    ref = REFERENCE.get(run)
    exit_s = "never" if meas[2] is None else f"{meas[2]:.2f} ms"
    if ref is None:
        print(f"  {run}: extreme {meas[0]:.0f}, residual {meas[1]:.1f}, "
              f"front gone {exit_s}")
        return
    nt.check(f"{run}: extreme mean dp [dyn/cm2]", meas[0], ref[0], "{:.0f}")
    nt.check(f"{run}: residual at 0.6 ms [dyn/cm2]", meas[1], ref[1],
             "{:.1f}")
    ref_s = "never" if ref[2] is None else f"{ref[2]:.2f} ms"
    print(f"  {run}: front gone {exit_s}   (reference {ref_s})")


def symlog_axis(ax):
    ax.set_yscale("symlog", linthresh=100.0)
    ax.axhspan(-0.001 * nt.P_AMB, 0.001 * nt.P_AMB, color="0.9", zorder=0)
    ax.axhline(0.0, color="0.5", lw=0.6)
    ax.set_xlabel("t [ms]      (grey band: ±0.1% of p_amb)")
    ax.set_ylabel("domain-mean p − p_amb [dyn/cm²]")


fig, axes = nt.plt.subplots(2, 2, figsize=(13, 9), layout="constrained")
ax_i, ax_f, ax_s, ax_m = axes.ravel()

for group, speed, ax in [(INERT, U, ax_i), (FLAME, U - S_L, ax_f)]:
    for run, label, col, ls in group:
        h = history(run, speed)
        if h is None:
            print(f"  {run}: not run (or not finished)")
            continue
        report(run, h)
        t, dp, lead, _ = h
        ax.plot(1e3 * t, dp, color=col, ls=ls, label=label)
        if run in ("inert_defaults", "inert_closures"):
            ax_s.plot(1e3 * t, 10 * lead, color=col, ls="--",
                      label=f"inert, {label}")
        elif group is FLAME and "ndnr" not in run:
            ax_s.plot(1e3 * t, 10 * lead, color=col, ls=ls,
                      label=f"flame, {label}")
    symlog_axis(ax)
    ax.legend(fontsize=8, loc="lower left" if group is INERT else "best")
ax_i.set_title("(A) an inert hot front leaving (no chemistry)")
ax_f.set_title("(B) a flame leaving")

ax_s.axhline(0.0, color="0.5", lw=0.6)
ax_s.set_xlabel("t [ms]")
ax_s.set_ylabel("front position − its schedule [mm]")
ax_s.set_title("is the front on schedule? (mean over the rows still holding it)")
ax_s.legend(fontsize=8)

# What blocked looks like: the temperature at 0.3 ms, defaults against the
# closures.
maps = []
for run in ("flame_defaults", "flame_closures"):
    files = nt.plotfiles(os.path.join("runs", run))
    if not files:
        continue
    d = min((nt.load(f, ["Temp"]) for f in files),
            key=lambda dd: abs(dd["time"] - 3.0e-4))
    maps.append((run, d))
if maps:
    stack = np.concatenate([d["Temp"].T for _, d in maps][::-1], axis=0)
    im = ax_m.imshow(stack, origin="lower", aspect="equal", cmap="inferno",
                     extent=[0, X_OUT, 0, 0.2 * len(maps)], vmin=290,
                     vmax=2500)
    ax_m.axhline(0.2, color="w", lw=1.5)
    for k, (run, d) in enumerate(maps[::-1]):
        lab = [lbl for r, lbl, _, _ in FLAME if r == run][0]
        gone = "" if np.any(d["Temp"] > 600.0) else "  (flame gone)"
        ax_m.text(0.01, 0.2 * k + 0.17, lab + gone, color="w", fontsize=8)
    ax_m.set_title(f"temperature at t = {1e3 * maps[0][1]['time']:.2f} ms "
                   "(outflow at the right edge)")
    ax_m.set_xlabel("x [cm]")
    ax_m.set_yticks([])
    fig.colorbar(im, ax=ax_m, label="T [K]", shrink=0.7)

nt.save(fig, "crossing.png")
