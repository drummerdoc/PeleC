"""Rung 8 figure: a vented chamber with its boundary at the vent and with
the boundary moved two chamber lengths downstream."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import nscbc_tut as nt  # noqa: E402

X_VENT = 1.2  # the vent plane: the end of the chamber and of the charge [cm]
H = 0.3  # chamber height [cm]
BUILDUP = (2.0e-3, 3.5e-3)  # the flame is in the chamber, venting steadily
END_OF_BURN = (3.8e-3, 6.0e-3)  # the flame reaches the vent plane

VARIANTS = [("vent", "boundary at the vent, σ = 0.25", "C3"),
            ("vent_s16", "boundary at the vent, σ = 16", "C1"),
            ("plenum", "boundary 2.4 cm downstream, σ = 0.25", "C0"),
            ("plenum_s16", "boundary 2.4 cm downstream, σ = 16", "C2")]
# Reference, closed-end overpressure [dyn/cm^2]: mean during the buildup,
# largest value at the end of burn, and the dominant frequency after it [kHz]
# (None: no clear ring).
REFERENCE = {"vent": (595, 37178, 2.74), "vent_s16": (-15, 6874, 13.83),
             "plenum": (1817, 2886, None), "plenum_s16": (-72, 1133, 2.50)}
# Reference: rms difference [cm^2/s] between the sigma = 0.25 and sigma = 16
# vent-plane flux while the flame reaches the vent plane, per layout.
REFERENCE_PAIR = {"vent": 67.0, "plenum": 2.5}


def traces(run):
    """Closed-end overpressure and the volume flux through the vent plane,
    one row per plotfile.  Cached in runs/."""
    cache = os.path.join("runs", f"{run}-traces.npy")
    files = nt.plotfiles(os.path.join("runs", run))
    if os.path.exists(cache):
        rows = np.load(cache)
        if len(rows) == len(files):
            return rows
    rows = []
    for f in files:
        d = nt.load(f, ["pressure", "x_velocity"])
        x, p, u = d["x"], d["pressure"], d["x_velocity"]
        iv = np.nonzero(x < X_VENT)[0][-1]  # last column inside the chamber
        rows.append([d["time"], p[x < 0.1 * X_VENT].mean() - nt.P_AMB,
                     u[iv].sum() * H / p.shape[1]])
    rows = np.array(rows)
    np.save(cache, rows)
    return rows


def ring_frequency(t, dp):
    """Dominant frequency of a detrended trace (Hann-windowed FFT,
    parabolic peak)."""
    w = max(5, len(t) // 8) | 1
    s = dp - np.convolve(dp, np.ones(w) / w, mode="same")
    sp = np.abs(np.fft.rfft(s * np.hanning(len(s))))
    fr = np.fft.rfftfreq(len(s), np.median(np.diff(t)))
    i = 1 + int(np.argmax(sp[1:]))
    if i < len(sp) - 1:
        a, b, c = sp[i - 1], sp[i], sp[i + 1]
        i = i + 0.5 * (a - c) / (a - 2 * b + c)
    return i * fr[1]


def measures(tr):
    t, dp = tr[:, 0], tr[:, 1]
    b = (t >= BUILDUP[0]) & (t <= BUILDUP[1])
    e = (t >= END_OF_BURN[0]) & (t <= END_OF_BURN[1])
    k = np.nonzero(e)[0][np.argmax(dp[e])]
    late = t >= t[k] + 5e-4
    f = ring_frequency(t[late], dp[late]) if late.sum() > 16 else np.nan
    return dp[b].mean(), dp[k], 1e-3 * f


fig = nt.plt.figure(figsize=(14, 9), layout="constrained")
gs = fig.add_gridspec(3, 3, height_ratios=[1.0, 1.0, 0.75])
ax_p = fig.add_subplot(gs[0, :2])
ax_z = fig.add_subplot(gs[1, 0])
ax_r = fig.add_subplot(gs[1, 1])
ax_v = fig.add_subplot(gs[:2, 2])
ax_m = fig.add_subplot(gs[2, :])

flux = {}
for run, label, col in VARIANTS:
    files = nt.plotfiles(os.path.join("runs", run))
    if not files or not os.path.exists(os.path.join("runs", run, "done")):
        print(f"  {run}: not run (or not finished)")
        continue
    tr = traces(run)
    meas = measures(tr)
    ref = REFERENCE.get(run)
    if ref:
        nt.check(f"{run}: buildup mean [dyn/cm2]", meas[0], ref[0], "{:.0f}")
        nt.check(f"{run}: end-of-burn peak [dyn/cm2]", meas[1], ref[1],
                 "{:.0f}")
        if ref[2] is not None:
            nt.check(f"{run}: ring-down frequency [kHz]", meas[2], ref[2],
                     "{:.2f}")
    else:
        print(f"  {run}: buildup {meas[0]:.0f}, end-of-burn peak "
              f"{meas[1]:.0f} dyn/cm2, ring-down {meas[2]:.2f} kHz")
    flux[run] = tr[:, :3:2]
    t = 1e3 * tr[:, 0]
    ax_p.plot(t, tr[:, 1], color=col, lw=1.0, label=label)
    ax_z.plot(t, tr[:, 1], color=col, lw=1.0)
    ax_r.plot(t, tr[:, 1], color=col, lw=1.0)
    ax_v.plot(t, tr[:, 2], color=col, lw=1.0, label=label)

# The two-sigma test: does the flow out of the chamber depend on the
# boundary's setting while the flame reaches the vent plane?
for lay, ref in REFERENCE_PAIR.items():
    if lay in flux and f"{lay}_s16" in flux:
        a, b = flux[lay], flux[f"{lay}_s16"]
        m = (a[:, 0] >= 3.8e-3) & (a[:, 0] <= 4.8e-3)
        d = a[m, 1] - np.interp(a[m, 0], b[:, 0], b[:, 1])
        nt.check(f"{lay}: vent flux, sigma 0.25 vs 16, rms diff [cm2/s]",
                 np.sqrt(np.mean(d ** 2)), ref, "{:.1f}")

ax_p.set_ylabel("closed-end p − p_amb [dyn/cm²]")
ax_p.set_xlabel("t [ms]")
ax_p.set_title("overpressure at the closed end of the chamber")
ax_p.legend(fontsize=8, ncol=2)
for ax in (ax_p, ax_z):
    ax.axvspan(1e3 * END_OF_BURN[0], 1e3 * END_OF_BURN[1], color="0.93",
               zorder=0)
ax_z.set_xlim(3.6, 5.6)
ax_z.set_xlabel("t [ms]")
ax_z.set_ylabel("closed-end p − p_amb [dyn/cm²]")
ax_z.set_title("the flame reaches the vent plane")
ax_r.set_xlim(2.0, 3.0)
ax_r.set_ylim(-3000, 6000)
ax_r.set_xlabel("t [ms]")
ax_r.set_title("during the burn (zoom)")
ax_v.axhline(0.0, color="0.5", lw=0.6)
ax_v.set_xlabel("t [ms]")
ax_v.set_ylabel("volume flux through the vent plane [cm²/s]")
ax_v.set_title("flow out through the vent")

# The two layouts at 6 ms, on one length scale: the chamber alone, and the
# chamber with the duct, where the charge pushed out ahead of the flame is
# still burning.
maps = []
for run in ("vent", "plenum"):
    files = nt.plotfiles(os.path.join("runs", run))
    if files:
        d = min((nt.load(f, ["Temp"]) for f in files),
                key=lambda dd: abs(dd["time"] - 6.0e-3))
        maps.append((run, d))
for k, (run, d) in enumerate(maps):
    y0 = (len(maps) - 1 - k) * (H + 0.08)
    im = ax_m.imshow(d["Temp"].T, origin="lower", cmap="inferno", vmin=290,
                     vmax=1500, extent=[0, d["x"][-1] + 0.5 * (d["x"][1]
                                                              - d["x"][0]),
                                        y0, y0 + H])
    ax_m.text(0.02, y0 + H + 0.01, f"{run}: boundary at x = "
              f"{d['x'][-1] + 0.5 * (d['x'][1] - d['x'][0]):.1f} cm",
              fontsize=8)
if maps:
    ax_m.axvline(X_VENT, color="w", ls=":", lw=1.0)
    ax_m.set_xlim(0, 3.6)
    ax_m.set_ylim(0, len(maps) * (H + 0.08))
    ax_m.set_aspect("equal")
    ax_m.set_yticks([])
    ax_m.set_xlabel("x [cm]   (dotted: the vent plane)")
    ax_m.set_title(f"temperature at t = {1e3 * maps[0][1]['time']:.1f} ms")
    fig.colorbar(im, ax=ax_m, label="T [K]", shrink=0.8)

nt.save(fig, "longer_domain.png")
