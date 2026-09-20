#!/usr/bin/env python3
"""V0: compare a tgvSpectral.py output (or any 't Ek eps ...' file) with a reference file (UCL gdiag / Incompact3d .dat).
Reports peak dissipation & time, L2 differences of K(t) and eps(t) over the common interval, and max relative deviation.
Usage: compareRef.py <ours.dat> <ref.dat> [--epscol 2] [--enstcol 3]"""
import numpy as np, sys, argparse
ap = argparse.ArgumentParser(); ap.add_argument("ours"); ap.add_argument("ref"); ap.add_argument("--tmax", type=float, default=20.0)
a = ap.parse_args()
o = np.loadtxt(a.ours, comments="#"); r = np.loadtxt(a.ref, comments="#")
to, Ko, epso = o[:, 0], o[:, 1], o[:, 2]; tr, Kr, epsr = r[:, 0], r[:, 1], r[:, 2]
# ours: eps column may be -dEk/dt (finite difference) — use 2nu*enstrophy (col 4) if present as the smooth dissipation
if o.shape[1] >= 5: epso_s = o[:, 4]
else: epso_s = epso
tmax = min(a.tmax, to.max(), tr.max()); tt = np.linspace(0, tmax, 2001)
Ki, Kri = np.interp(tt, to, Ko), np.interp(tt, tr, Kr)
Ei, Eri = np.interp(tt, to, epso_s), np.interp(tt, tr, epsr)
print("ours : peak eps = %.6g at t = %.3f ; K(t_end) = %.8f" % (epso_s.max(), to[np.argmax(epso_s)], Ko[-1]))
print("ref  : peak eps = %.6g at t = %.3f ; K(t_end) = %.8f" % (epsr.max(), tr[np.argmax(epsr)], np.interp(tmax, tr, Kr)))
print("peak eps rel. diff = %.3e ; peak time diff = %.3f" % (epso_s.max()/epsr.max() - 1, to[np.argmax(epso_s)] - tr[np.argmax(epsr)]))
print("L2(K)   over [0,%g] = %.3e (rel. to K0: %.3e)" % (tmax, np.sqrt(np.mean((Ki - Kri)**2)), np.sqrt(np.mean((Ki - Kri)**2))/0.125))
print("L2(eps) over [0,%g] = %.3e (rel. to peak: %.3e) ; max|rel| eps = %.3e" % (tmax, np.sqrt(np.mean((Ei - Eri)**2)), np.sqrt(np.mean((Ei - Eri)**2))/epsr.max(), np.max(np.abs(Ei - Eri)/np.maximum(Eri, 1e-12))))
