"""A periodic box of synthetic homogeneous isotropic turbulence, written in
the text format of PelePhysics' Support/TurbInflowGenerator/gen_hit_ic.py
(cell centres on [0, 2 pi]^3, Fortran order, u_rms = 1), from numpy alone.

    python make_hit.py <N> <k0> <outfile>

Divergence-free random Fourier modes with a Passot-Pouquet spectrum
E(k) ~ k^4 exp(-2 k^2 / k0^2), peaked at wavenumber k0.  The seed is fixed,
so every user gets the same field and can match the reference numbers.
"""

import sys

import numpy as np

N, k0, out = int(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
rng = np.random.default_rng(20261031)

k1 = np.fft.fftfreq(N, d=1.0 / N)  # integer wavenumbers on a 2 pi box
kx, ky, kz = np.meshgrid(k1, k1, k1, indexing="ij")
kk = np.sqrt(kx**2 + ky**2 + kz**2)
kk[0, 0, 0] = 1.0

a = rng.normal(size=(3, N, N, N)) + 1j * rng.normal(size=(3, N, N, N))
kdota = (kx * a[0] + ky * a[1] + kz * a[2]) / kk**2
a[0] -= kx * kdota  # project out the compressive part
a[1] -= ky * kdota
a[2] -= kz * kdota
spectrum = kk**4 * np.exp(-2.0 * kk**2 / k0**2)
a *= np.sqrt(spectrum / kk**2)
a[:, 0, 0, 0] = 0.0
a[:, kk > N / 3] = 0.0  # keep the field resolved on the box

u = np.real(np.fft.ifftn(a, axes=(1, 2, 3)))
u /= np.sqrt(np.mean(u**2, axis=(1, 2, 3)).mean())  # u_rms = 1 per component

xc = (np.arange(N) + 0.5) * 2.0 * np.pi / N
X, Y, Z = np.meshgrid(xc, xc, xc, indexing="ij")
cols = [c.reshape(-1, order="F") for c in (X, Y, Z, u[0], u[1], u[2])]
np.savetxt(out, np.vstack(cols).T, fmt="%.18e", delimiter=",",
           header="x, y, z, u, v, w")
rms = np.sqrt(np.mean(u**2, axis=(1, 2, 3)))
print(f"wrote {out}: N = {N}, k0 = {k0:g}, rms u v w = "
      + " ".join(f"{r:.3f}" for r in rms))
