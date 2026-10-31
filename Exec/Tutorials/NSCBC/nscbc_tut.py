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
    """Sorted plotfile directories of one run."""
    return sorted(
        d for d in glob.glob(os.path.join(run_dir, "plt*")) if os.path.isdir(d)
    )


def _parse_box(text):
    """'((lo...) (hi...) (type...))' -> (lo, hi) integer tuples."""
    import re

    nums = [list(map(int, g.split(","))) for g in re.findall(r"\(([-\d,]+)\)", text)]
    return nums[0], nums[1]


def load(pltfile, fields):
    """Level-0 arrays of the given fields (shape nx, ny[, nz]) from an AMReX
    plotfile, plus its time and the x cell centres.  Reads the files
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
    data["x"] = lo[0] + (np.arange(shape[0]) + 0.5) * dx[0]
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
