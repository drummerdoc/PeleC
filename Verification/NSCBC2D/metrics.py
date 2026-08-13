#!/usr/bin/env python3
"""
Metrics for the 2-D NSCBC tests (Exec/RegTests/NSCBC-COVO).

Reads the regular arrays written by fielddump and reports:

  circularity  -- for the circular-pulse case, the spread in the radius of the
                  outgoing wavefront over polar angle, restricted to the rays
                  along which the front is still inside the domain.  The exact
                  solution is a circle, so any spread is boundary error.  This
                  is the diagnostic used by Motheau et al. (2017): a good
                  non-reflecting boundary lets the pressure contours stay
                  circular as the wave crosses; a poor one flattens and buckles
                  them, worst at the corners where the wave arrives obliquely.

  residual     -- what is left in the domain after the wave has gone, relative
                  to the incident amplitude.  For the vortex case this is the
                  whole story: a vortex carries no acoustic content, so any
                  pressure signal left behind was manufactured by the boundary.

Usage:
    metrics.py circularity <p_amb> <c> <file> [<file> ...]
    metrics.py residual    <p_amb> <ref_amp_file> <file> [<file> ...]
"""
import sys
import numpy as np


def load(fn):
    with open(fn) as f:
        f.readline()
        nx, ny, xlo, ylo, dx, dy, t = f.readline().split()
        nx, ny = int(nx), int(ny)
        a = np.loadtxt(f)
    a = a.reshape(ny, nx)
    x = float(xlo) + (np.arange(nx) + 0.5) * float(dx)
    y = float(ylo) + (np.arange(ny) + 0.5) * float(dy)
    return x, y, a, float(t)


def bilinear(x, y, a, xq, yq):
    """Sample a on a regular grid at scattered (xq, yq)."""
    dx, dy = x[1] - x[0], y[1] - y[0]
    fi = (xq - x[0]) / dx
    fj = (yq - y[0]) / dy
    i0 = np.clip(np.floor(fi).astype(int), 0, len(x) - 2)
    j0 = np.clip(np.floor(fj).astype(int), 0, len(y) - 2)
    tx = np.clip(fi - i0, 0.0, 1.0)
    ty = np.clip(fj - j0, 0.0, 1.0)
    return ((1 - tx) * (1 - ty) * a[j0, i0] + tx * (1 - ty) * a[j0, i0 + 1] +
            (1 - tx) * ty * a[j0 + 1, i0] + tx * ty * a[j0 + 1, i0 + 1])


def dist_to_boundary(theta, x, y):
    """Distance from the origin to the domain boundary along each ray."""
    xmax, ymax = x[-1], y[-1]
    xmin, ymin = x[0], y[0]
    ct, st = np.cos(theta), np.sin(theta)
    big = 1e30
    with np.errstate(divide="ignore", invalid="ignore"):
        dx1 = np.where(ct > 0, xmax / ct, np.where(ct < 0, xmin / ct, big))
        dy1 = np.where(st > 0, ymax / st, np.where(st < 0, ymin / st, big))
    return np.minimum(dx1, dy1)


def circularity(p_amb, c, files, ntheta=720):
    print()
    print("  Wavefront circularity -- radius of the outgoing front vs polar angle")
    print("  (exact solution is a circle; every deviation is boundary error)")
    print()
    print("  %-16s %8s %8s %10s %10s %10s %10s" %
          ("file", "t/tau", "rays", "r_mean", "spread%", "amp_sprd%", "peak|dp|"))
    print("  " + "-" * 82)
    for fn in files:
        x, y, a, t = load(fn)
        dp = a - p_amb
        r_th = c * t
        theta = np.linspace(0.0, 2 * np.pi, ntheta, endpoint=False)
        dbnd = dist_to_boundary(theta, x, y)
        # Only rays whose front is still comfortably inside the domain.
        keep = r_th < 0.90 * dbnd
        if keep.sum() < 8 or r_th <= 0:
            print("  %-16s %8.3f %8s" % (fn.split("/")[-1], t, "front gone"))
            continue
        th = theta[keep]
        # Scan each retained ray for the front.
        # The radial scan must be finer than the effect being measured: at 400
        # samples one increment is ~0.2 % of r_th, which quantises the radius
        # spread and makes two genuinely different boundaries report the same
        # number.  2400 puts the quantisation an order of magnitude below the
        # grid spacing.
        rr = np.linspace(0.55 * r_th, 1.45 * r_th, 2400)
        R, TH = np.meshgrid(rr, th, indexing="ij")
        vals = bilinear(x, y, dp, R * np.cos(TH), R * np.sin(TH))
        k = np.argmax(np.abs(vals), axis=0)
        r_peak = rr[k]
        amp = np.abs(vals[k, np.arange(len(th))])
        spread = (r_peak.max() - r_peak.min()) / r_peak.mean()
        amp_spread = (amp.max() - amp.min()) / amp.mean()
        print("  %-16s %8.3f %8d %10.4f %10.3f %10.3f %10.4e" %
              (fn.split("/")[-1], r_th / max(x[-1], 1e-30), keep.sum(),
               r_peak.mean(), 100 * spread, 100 * amp_spread, amp.mean()))


def residual(p_amb, ref_file, files):
    _, _, a0, _ = load(ref_file)
    inc = np.abs(a0 - p_amb).max()
    print()
    print("  Residual pressure disturbance, relative to the incident amplitude")
    print("  incident |dp|_max = %.5e dyn/cm^2" % inc)
    print()
    print("  %-16s %10s %12s %12s %14s" %
          ("file", "t [s]", "max|dp|/inc", "L2|dp|/inc", "mean p"))
    print("  " + "-" * 68)
    for fn in files:
        x, y, a, t = load(fn)
        dp = a - p_amb
        l2 = np.sqrt((dp ** 2).mean())
        print("  %-16s %10.4e %12.5f %12.6f %14.4f" %
              (fn.split("/")[-1], t, np.abs(dp).max() / inc, l2 / inc, a.mean()))


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "circularity":
        circularity(float(sys.argv[2]), float(sys.argv[3]), sys.argv[4:])
    elif mode == "residual":
        residual(float(sys.argv[2]), sys.argv[3], sys.argv[4:])
    else:
        print(__doc__)
        sys.exit(1)
