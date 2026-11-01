"""Rung 3 figure: an inflow's step response and its injection fidelity."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import nscbc_tut as nt  # noqa: E402

CS, L = 34719.0, 10.0
T_A = 2.0 * L / CS  # acoustic round trip
U_INIT, U_TGT = 2000.0, 2100.0
FREQ, AMP = 692.1, 34.7  # forcing (0.8 f0) and its target amplitude [cm/s]
KS_STEP = ["0.1", "0.3", "1", "3", "10"]
KS_FORCE = ["0.5", "2", "5"]
MODES = [("plain", "relaxation only"), ("ff", "+ dudt feed-forward"),
         ("ffnri", "+ dudt + NRI register")]
# Reference: t_conv [t_a] per step run; I_in per forced run.
REFERENCE = {
    "step_k0.1": 20.0, "step_k0.3": 4.65, "step_k1": 7.64,
    "step_k3": 19.9, "step_k10": 20.0,
    "force_k0.5_plain": 0.223, "force_k2_plain": 1.48, "force_k5_plain": 4.56,
    "force_k0.5_ff": 1.13, "force_k2_ff": 2.34, "force_k5_ff": 5.08,
    "force_k0.5_ffnri": 0.960, "force_k2_ffnri": 0.963, "force_k5_ffnri": 1.01,
}


def step_response(k):
    log = nt.read_datlog(os.path.join("runs", f"step_k{k}", "datlog"))
    ubar = log["xmom"] / log["mass"]
    return log["time"] / T_A, (ubar - U_INIT) / (U_TGT - U_INIT)


def injection(run):
    """I_in (incoming wave over the ideal injector's) and I_u (inlet velocity
    amplitude over the target's), projected at the forcing frequency."""
    t, u, p, rho = [], [], [], []
    for f in nt.plotfiles(os.path.join("runs", run)):
        d = nt.load(f, ["x_velocity", "pressure", "density"])
        t.append(d["time"])
        u.append(d["x_velocity"][0].mean())  # the inlet column
        p.append(d["pressure"][0].mean())
        rho.append(d["density"][0].mean())
    t, u, p = map(np.array, (t, u, p))
    rc = np.mean(rho) * CS
    m = t >= 6.0 * T_A
    tm, w = t[m], np.gradient(t[m])
    om = 2.0 * np.pi * FREQ

    def amp(sig):
        s = sig[m] - np.polyval(np.polyfit(tm, sig[m], 1), tm)  # detrend
        a = (2.0 / w.sum()) * np.sum(s * np.sin(om * tm) * w)
        b = (2.0 / w.sum()) * np.sum(s * np.cos(om * tm) * w)
        return np.hypot(a, b)

    return amp(u + p / rc) / (2.0 * AMP), amp(u) / AMP


fig, (ax_s, ax_f) = nt.plt.subplots(1, 2, figsize=(13, 4.4), layout="constrained")
cols = nt.plt.cm.plasma(np.linspace(0, 0.85, len(KS_STEP)))
print("(A) step: t_conv = last time |u - target| > 1% of the step")
for k, col in zip(KS_STEP, cols):
    tt, r = step_response(k)
    ax_s.plot(tt, r, color=col, label=f"relax_u = {k}")
    off = np.where(np.abs(r - 1.0) > 0.01)[0]
    tconv = tt[off[-1]] if len(off) else 0.0
    ref = REFERENCE.get(f"step_k{k}")
    if ref is not None:
        nt.check(f"relax_u = {k}: t_conv [t_a]", tconv, ref, "{:.3g}")
    else:
        print(f"  relax_u = {k}: t_conv = {tconv:.3g} t_a")
ax_s.axhline(1.0, color="k", lw=0.6)
ax_s.set_xlabel(r"t / t$_a$  (t$_a$ = 2L/c, one round trip)")
ax_s.set_ylabel(r"($\bar u$ − 2000) / step")
ax_s.set_title("(A) the inlet target steps 2000 → 2100 cm/s")
ax_s.legend(fontsize=8)

print("(B) forced inlet at 0.8 f0: I_in (1 = the wave asked for), I_u")
width = 0.26
for j, (mode, label) in enumerate(MODES):
    Iin, Iu = [], []
    for k in KS_FORCE:
        a, b = injection(f"force_k{k}_{mode}")
        Iin.append(a)
        Iu.append(b)
        ref = REFERENCE.get(f"force_k{k}_{mode}")
        if ref is not None:
            nt.check(f"relax_u = {k}, {mode}: I_in", a, ref, "{:.3g}")
        else:
            print(f"  relax_u = {k}, {mode}: I_in = {a:.3g}, I_u = {b:.3g}")
    xs = np.arange(len(KS_FORCE)) + (j - 1) * width
    ax_f.bar(xs, Iin, width, label=label)
    for x, v in zip(xs, Iin):
        ax_f.text(x, v + 0.08, f"{v:.2f}", ha="center", fontsize=8)
ax_f.axhline(1.0, color="k", ls="--", lw=0.8)
ax_f.set_xticks(range(len(KS_FORCE)), [f"relax_u = {k}" for k in KS_FORCE])
ax_f.set_ylabel(r"injection index I$_{in}$  (1 = exact)")
ax_f.set_title("(B) forcing the inlet at 0.8 f$_0$")
ax_f.legend(fontsize=8)
nt.save(fig, "inflow.png")
