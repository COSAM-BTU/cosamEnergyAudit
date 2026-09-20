#!/usr/bin/env python3
"""Checkerboard indicators of reconstructed OpenFOAM fields on the uniform N^3 blockMesh (review round 3, A6).
For a cell field q: (1) Nyquist fraction  E_Nyq = sum |q_hat|^2 over modes with any |k_i| = N/2  /  sum |q_hat|^2 (mean excluded);
(2) odd-even indicator  chi = rms( q_P - mean of the six face neighbours ) / rms( q - <q> ).  Both for p and for |U| components.
Usage: checkerboard.py <case> <N> <time> [<time> ...]"""
import sys, os, re, numpy as np
def read_field(path, ncomp):
    b = open(path, "rb").read(); i = b.index(b"internalField"); j = b.index(b"\n", i) + 1
    k = b.index(b"\n", j); n = int(b[j:k].strip()); l = b.index(b"(", k) + 1
    a = np.frombuffer(b[l:l + 8 * n * ncomp], dtype="<f8").reshape(n, ncomp) if ncomp > 1 else np.frombuffer(b[l:l + 8 * n], dtype="<f8")
    return a
def indicators(q, N):
    q = q.reshape(N, N, N, order="F") if q.size == N**3 else None   # blockMesh: i fastest
    q = q - q.mean(); qh = np.fft.fftn(q); E = np.abs(qh)**2; E[0, 0, 0] = 0.0
    kx = np.fft.fftfreq(N) * N; nyq = np.abs(kx) == N // 2
    mask = nyq[:, None, None] | nyq[None, :, None] | nyq[None, None, :]
    fN = E[mask].sum() / E.sum()
    hi = np.abs(kx) >= 3 * N // 8; band = hi[:, None, None] | hi[None, :, None] | hi[None, None, :]; fB = E[band].sum() / E.sum()
    nb = (np.roll(q, 1, 0) + np.roll(q, -1, 0) + np.roll(q, 1, 1) + np.roll(q, -1, 1) + np.roll(q, 1, 2) + np.roll(q, -1, 2)) / 6.0
    chi = np.sqrt(np.mean((q - nb)**2)) / np.sqrt(np.mean(q**2))
    return fN, fB, chi
case, N = sys.argv[1], int(sys.argv[2])
for t in sys.argv[3:]:
    d = os.path.join(case, t); p = read_field(os.path.join(d, "p"), 1); U = read_field(os.path.join(d, "U"), 3)
    fp, bp, cp = indicators(p, N); fu = []; bu = []; cu = []
    for c in range(3):
        f, b, ch = indicators(U[:, c].copy(), N); fu.append(f); bu.append(b); cu.append(ch)
    print("CB %-52s t=%-4s p: Nyq=%.2e band=%.3e chi=%.4f | U: Nyq=%.2e band=%.3e chi=%.4f" % (os.path.basename(case).replace("_nC2_dc1_pe10", ""), t, fp, bp, cp, np.mean(fu), np.mean(bu), np.mean(cu)), flush=True)
