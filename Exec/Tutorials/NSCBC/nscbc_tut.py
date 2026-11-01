"""Shared helpers for the NSCBC tutorial plot scripts.

Needs numpy and matplotlib only; plotfiles are read directly.
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
    """Plotfile directories of one run, in step order (plt99999 comes before
    plt100000)."""
    dirs = [d for d in glob.glob(os.path.join(run_dir, "plt*"))
            if os.path.isdir(d) and os.path.basename(d)[3:].isdigit()]
    return sorted(dirs, key=lambda d: int(os.path.basename(d)[3:]))


def _parse_box(text):
    """'((lo...) (hi...) (type...))' -> (lo, hi) integer tuples."""
    import re

    nums = [list(map(int, g.split(","))) for g in re.findall(r"\(([-\d,]+)\)", text)]
    return nums[0], nums[1]


def load(pltfile, fields):
    """Level-0 arrays of the given fields (shape nx, ny[, nz]) from an AMReX
    plotfile, plus its time and the cell centres x, y[, z].  Reads the files
    directly with numpy: single-level, cell-centred data, double precision,
    which is what every rung writes."""
    with open(os.path.join(pltfile, "Header")) as f:
        lines = f.read().split("\n")
    ncomp = int(lines[1])
    names = [lines[2 + i].strip() for i in range(ncomp)]
    k = 2 + ncomp
    dim = int(lines[k])
    time = float(lines[k + 1])
    lo = np.array(lines[k + 3].split(), dtype=float)
    hi = np.array(lines[k + 4].split(), dtype=float)
    dlo, dhi = _parse_box(lines[k + 6])
    shape = tuple(dhi[d] - dlo[d] + 1 for d in range(dim))

    lev = os.path.join(pltfile, "Level_0")
    with open(os.path.join(lev, "Cell_H")) as f:
        h = f.read().split("\n")
    boxes, fabs = [], []
    for line in h:
        if line.startswith("((") and line.count("(") >= 4:  # a box
            boxes.append(_parse_box(line))
        elif line.startswith("FabOnDisk:"):
            _, fname, off = line.split()
            fabs.append((fname, int(off)))
    data = {n: np.empty(shape) for n in fields}
    for (blo, bhi), (fname, off) in zip(boxes, fabs):
        with open(os.path.join(lev, fname), "rb") as f:
            f.seek(off)
            header = f.readline().decode()
            nc = int(header.split()[-1])
            n = [bhi[d] - blo[d] + 1 for d in range(dim)]
            raw = np.fromfile(f, dtype="<f8", count=int(np.prod(n)) * nc)
        # Fortran order: x fastest, one component after another.
        raw = raw.reshape([nc] + n[::-1])
        sl = tuple(slice(blo[d] - dlo[d], bhi[d] - dlo[d] + 1) for d in range(dim))
        for fld in fields:
            data[fld][sl] = raw[names.index(fld)].T
    dx = (hi - lo) / np.array(shape)
    for d, name in enumerate("xyz"[:dim]):
        data[name] = lo[d] + (np.arange(shape[d]) + 0.5) * dx[d]
    data["time"] = time
    return data


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


def front_amplitude_spread(dp, x, y, radius, ntheta=720):
    """Spread of the peak |dp| along rays from the domain centre, over the
    rays whose front (at `radius`) is still well inside the domain: the
    circular front's strength should be the same in every direction, so any
    variation is boundary error.  Returns (spread, number of rays)."""
    x0, y0 = 0.5 * (x[0] + x[-1]), 0.5 * (y[0] + y[-1])
    th = np.linspace(0.0, 2.0 * np.pi, ntheta, endpoint=False)
    ct, st = np.cos(th), np.sin(th)
    with np.errstate(divide="ignore"):
        to_x = np.where(ct != 0, (x[-1] - x0) / np.abs(ct), np.inf)
        to_y = np.where(st != 0, (y[-1] - y0) / np.abs(st), np.inf)
    keep = radius < 0.9 * np.minimum(to_x, to_y)
    rr = np.linspace(0.55 * radius, 1.45 * radius, 2400)
    R, T = np.meshgrid(rr, th[keep], indexing="ij")
    fx = (x0 + R * np.cos(T) - x[0]) / (x[1] - x[0])
    fy = (y0 + R * np.sin(T) - y[0]) / (y[1] - y[0])
    i = np.clip(np.floor(fx).astype(int), 0, len(x) - 2)
    j = np.clip(np.floor(fy).astype(int), 0, len(y) - 2)
    tx, ty = np.clip(fx - i, 0, 1), np.clip(fy - j, 0, 1)
    v = ((1 - tx) * (1 - ty) * dp[i, j] + tx * (1 - ty) * dp[i + 1, j]
         + (1 - tx) * ty * dp[i, j + 1] + tx * ty * dp[i + 1, j + 1])
    amp = np.abs(v).max(axis=0)
    return (amp.max() - amp.min()) / amp.mean(), int(keep.sum())


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
