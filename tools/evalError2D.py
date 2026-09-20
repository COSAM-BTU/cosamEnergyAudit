#!/usr/bin/env python3
"""L2/Linf error of U at the latest written time vs the analytic 2D TGV (Re=100 default): u = e^{-2 nu t}(sin x cos y, -cos x sin y).
Requires ASCII U and cell centres C (postProcess -func writeCellCentres -latestTime). Usage: evalError2D.py <case> [nu]"""
import re, sys, os, glob, numpy as np
case = sys.argv[1]; nu = float(sys.argv[2]) if len(sys.argv) > 2 else 0.01
times = sorted([float(d) for d in os.listdir(case) if re.match(r"^[0-9.]+$", d) and float(d) > 0]); t = times[-1]; td = os.path.join(case, "%g" % t)
def vec(f):
    s = open(f).read(); body = s[s.find("internalField"):]
    n = int(re.search(r"List<vector>\s*(\d+)", body).group(1))
    vals = re.findall(r"\(([-0-9.e+]+) ([-0-9.e+]+) ([-0-9.e+]+)\)", body)[:n]
    return np.array(vals, dtype=float)
U = vec(os.path.join(td, "U")); C = vec(os.path.join(td, "C"))
ue = np.exp(-2*nu*t)*np.column_stack([np.sin(C[:, 0])*np.cos(C[:, 1]), -np.cos(C[:, 0])*np.sin(C[:, 1]), np.zeros(len(C))])
err = U - ue; Umax = np.abs(ue).max()
print("t=%g N=%d L2=%.6e Linf=%.6e (relative to |u|max=%.4f)" % (t, len(U), np.sqrt(np.mean(np.sum(err**2, 1)))/Umax, np.abs(err).max()/Umax, Umax))
