"""Shared helpers for the NSCBC tutorial plot scripts.

Needs numpy, matplotlib and yt (pip install numpy matplotlib yt).
"""

import glob
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

P_AMB = 1.01325e6  # dyn/cm^2, the ambient/target pressure of every rung
GAMMA = 1.4  # GammaLaw default, used to turn the data log's rho_e into <p>


def plotfiles(run_dir):
    """Sorted plotfile directories of one run."""
    return sorted(
        d for d in glob.glob(os.path.join(run_dir, "plt*")) if os.path.isdir(d)
    )


def load(pltfile, fields):
    """Level-0 arrays of the given fields, plus time and cell centres."""
    import yt

    yt.set_log_level(50)
    ds = yt.load(pltfile)
    cg = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    out = {f: np.array(cg["boxlib", f]) for f in fields}
    lo = np.array(ds.domain_left_edge)
    hi = np.array(ds.domain_right_edge)
    n = np.array(ds.domain_dimensions)
    dx = (hi - lo) / n
    out["x"] = lo[0] + (np.arange(n[0]) + 0.5) * dx[0]
    out["time"] = float(ds.current_time)
    return out


def transverse_mean(a):
    """Average over every direction but x."""
    return a.reshape(a.shape[0], -1).mean(axis=1)


def reflection(p_final, amp):
    """Reflection coefficient: the largest residual left in the upstream
    half of the domain, measured from that half's median (so a uniform
    anchoring offset is not counted as a wave, and a strong reflected pulse
    does not drag the baseline the way a mean would), over the incident
    amplitude."""
    up = p_final[: len(p_final) // 2]
    return np.abs(up - np.median(up)).max() / amp


def read_datlog(path):
    """The pelec data log (amr.data_log) as a dict of columns."""
    with open(path) as f:
        names = f.readline().split()
    data = np.loadtxt(path, skiprows=1, ndmin=2)
    return {n: data[:, i] for i, n in enumerate(names)}


def mean_pressure(datlog, volume):
    """Domain-mean pressure history from the data log (GammaLaw gas)."""
    return datlog["time"], (GAMMA - 1.0) * datlog["rho_e"] / volume


def check(label, value, reference, fmt="{:.3g}", unit=""):
    """Print a measured number next to the tutorial's reference value."""
    v = fmt.format(value)
    r = fmt.format(reference)
    print(f"  {label:<44s} {v:>10s}{unit}   (reference {r}{unit})")


def save(fig, path):
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path}")
