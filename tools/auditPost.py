#!/usr/bin/env python3
"""Post-process a kineticEnergyAudit time series: CSV export, summary metrics, and PGFplots-ready column data.
Usage: auditPost.py <case_dir> [--ref referenceData/re1600/spectral_Re1600_512.gdiag] [--out post_<run_id>]"""
import numpy as np, sys, os, json, glob, argparse
ap = argparse.ArgumentParser(); ap.add_argument("case"); ap.add_argument("--ref", default=None); ap.add_argument("--out", default=None)
a = ap.parse_args()
cols = "t dt timeIndex CoMax K dKdt epsNu twoNuSS Omega twoNuOmega eTime eTimeDiss eTimeStore eConv eMesh eCont ePres eDev ePhi eIter sumE Rclosure cumClosure RcellMax chkConv chkDiff chkTime chkPres eConvOwner eMeshOwner epsNuOwner eTimeHalf eConvHalf eDiffHalf ePresHalf eIterHalf maskFrac sumFluxC sumFluxD sumFluxX eIterInstr chkInstr".split()
files = sorted(glob.glob(os.path.join(a.case, "postProcessing/kineticEnergyAudit/*/kineticEnergyAudit.dat")))
arrs = [x for x in (np.loadtxt(f, comments="#", ndmin=2) for f in files) if x.size]
if not arrs: sys.exit("no audit rows yet in %s" % a.case)
d = np.vstack(arrs); d = d[np.argsort(d[:, 0])]
n = d.shape[1]; c = {k: d[:, i] for i, k in enumerate(cols[:n])}
meta = json.load(open(os.path.join(a.case, "run_meta.json"))) if os.path.exists(os.path.join(a.case, "run_meta.json")) else {"run_id": os.path.basename(os.path.abspath(a.case))}
out = a.out or ("post_" + meta["run_id"]); os.makedirs(out, exist_ok=True)
np.savetxt(os.path.join(out, "audit_timeseries.csv"), d, delimiter=",", header=",".join(cols[:n]), comments="")
# metrics
t = c["t"]; dt = c["dt"]
eps_num = -c["dKdt"] - c["epsNu"]                       # total numerical dissipation rate (positive = loss)
ipk = int(np.argmax(-c["dKdt"]))
integ = lambda x: float(np.sum(x*dt))
m = dict(run_id=meta["run_id"], t_end=float(t[-1]), K_end=float(c["K"][-1]),
         peak_minus_dKdt=float(-c["dKdt"][ipk]), t_peak=float(t[ipk]), peak_epsNu=float(c["epsNu"].max()), t_peak_epsNu=float(t[np.argmax(c["epsNu"])]),
         int_epsNu=integ(c["epsNu"]), int_eNum=integ(eps_num), ratio_eNum_over_epsNu=integ(eps_num)/max(integ(c["epsNu"]), 1e-300),
         max_Rclosure=float(c["Rclosure"].max()), max_RcellMax=float(c["RcellMax"].max()), max_cumClosure=float(c["cumClosure"].max()),
         max_CoMax=float(c["CoMax"].max()),
         shares={k: integ(c[k])/max(integ(eps_num), 1e-300) for k in ("eTime", "eConv", "eMesh", "eCont", "ePres", "eDev")},
         share_eIter=-integ(c["eIter"])/max(integ(eps_num), 1e-300))
if a.ref and os.path.exists(a.ref):
    r = np.loadtxt(a.ref, comments="#"); rt, rE, rEps = r[:, 0], r[:, 1], r[:, 2]
    tt = t[t <= rt.max()]
    m["ref_peak_eps"] = float(rEps.max()); m["ref_t_peak"] = float(rt[np.argmax(rEps)])
    m["L2_K_vs_ref"] = float(np.sqrt(np.mean((np.interp(tt, t, c["K"]) - np.interp(tt, rt, rE))**2)))
    m["L2_eps_vs_ref"] = float(np.sqrt(np.mean((np.interp(tt, t, -c["dKdt"]) - np.interp(tt, rt, rEps))**2)))
json.dump(m, open(os.path.join(out, "metrics.json"), "w"), indent=1)
# PGFplots-ready stacked data (per unit volume rates) at a coarse stride
stride = max(1, len(t)//400)
hdr = "t epsNu eTime eConv eMesh eCont ePres eDev eIter minusdKdt eNum"
tab = np.column_stack([t, c["epsNu"], c["eTime"], c["eConv"], c["eMesh"], c["eCont"], c["ePres"], c["eDev"], c["eIter"], -c["dKdt"], eps_num])[::stride]
np.savetxt(os.path.join(out, "stacked.dat"), tab, header=hdr, comments="")
print(json.dumps(m, indent=1)[:1500])
