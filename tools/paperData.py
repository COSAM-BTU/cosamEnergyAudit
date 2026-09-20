#!/usr/bin/env python3
"""Paper data extraction (runs ON super00). Writes ~/tgvAudit/paperData/: metrics_all.csv (one row per run), ts_<run_id>.csv
(strided time series), ref_*.csv (references), v1p_diff.csv. Columns follow the FO header (INSTRUCTIONS §6)."""
import numpy as np, glob, os, json, csv, re
R="/scratch/soscfd00/tgvAudit"; OUT=os.path.expanduser("~/tgvAudit/paperData"); os.makedirs(OUT, exist_ok=True)
REF=os.path.expanduser("~/tgvAudit/stage1/spectral/ref/spectral_Re1600_512.gdiag"); SP=R+"/spectral"
cols="t dt timeIndex CoMax K dKdt epsNu twoNuSS Omega twoNuOmega eTime eTimeDiss eTimeStore eConv eMesh eCont ePres eDev ePhi eIter sumE Rclosure cumClosure RcellMax chkConv chkDiff chkTime chkPres eConvOwner eMeshOwner epsNuOwner eTimeHalf eConvHalf eDiffHalf ePresHalf eIterHalf maskFrac sumFluxC sumFluxD sumFluxX eIterInstr chkInstr".split()
def load(case):
    fs=sorted(glob.glob(case+"/postProcessing/kineticEnergyAudit/*/kineticEnergyAudit.dat")); arrs=[np.loadtxt(f,comments="#",ndmin=2) for f in fs]
    d=np.vstack([a for a in arrs if a.size]); d=d[np.argsort(d[:,0])]; return {k:d[:,i] for i,k in enumerate(cols[:d.shape[1]])}
def ref_for(Re):
    if Re==1600: u=np.loadtxt(REF,comments="#"); return u[:,0],u[:,1],u[:,2],2*u[:,3]/1600.,"UCL512"
    n={100:128,280:192}[Re]; u=np.loadtxt("%s/spec_Re%d_N%d.dat"%(SP,Re,n),comments="#"); return u[:,0],u[:,1],u[:,2],u[:,2],"specN%d"%n
# references
for Re,(tr,Kr,er,e2,nm) in [(Re,ref_for(Re)) for Re in (100,280,1600)]:
    np.savetxt(OUT+"/ref_Re%d_%s.csv"%(Re,nm), np.column_stack([tr,Kr,er])[::max(1,len(tr)//400)], delimiter=",", header="t,K,eps", comments="")
u=np.loadtxt(SP+"/spec_Re1600_N256.dat",comments="#"); np.savetxt(OUT+"/ref_Re1600_ownspec256.csv", np.column_stack([u[:,0],u[:,1],u[:,2]])[::10], delimiter=",", header="t,K,eps", comments="")
rows=[]
for blk in ("V2","S","Sint","F","T","V1p","Vadapt","Lchk","V3","Vpar"):
    for case in sorted(glob.glob(R+"/"+blk+"/*_h*")):
        if not glob.glob(case+"/postProcessing/kineticEnergyAudit/*/kineticEnergyAudit.dat"): continue
        meta=json.load(open(case+"/run_meta.json")) if os.path.exists(case+"/run_meta.json") else {}
        rid=meta.get("run_id",os.path.basename(case)); c=load(case); t=c["t"]; dt=c["dt"]
        Re=int(meta.get("Re",1600)); N=int(meta.get("N",0)); nproc=max(1,len(glob.glob(case+"/processor*")))
        I=lambda x: float(np.sum(x*dt))
        ipk=int(np.argmax(-c["dKdt"])); epsnum=-c["dKdt"]-c["epsNu"]; c["eNum"]=epsnum   # numerical dissipation rate = sum_i e_i - e_iter
        tr,Kr,er,e2,rn=ref_for(Re); m=tr<=t[-1]+1e-9; ir=int(np.argmax(er[m])); ei=np.interp(t,tr,er); Ki=np.interp(t,tr,Kr)
        log=case+"/log.pimpleFoam"; wall=None; steps=0
        if os.path.exists(log):
            for line in open(log):
                if line.startswith("Time = "): steps+=1
                if line.startswith("ExecutionTime"): wall=float(line.split()[2])
        cm=case+"/log.checkMesh"; nonorth=None; cells=None
        if os.path.exists(cm):
            s=open(cm).read(); mo=re.search(r"non-orthogonality Max: ([0-9.e+-]+)",s); nonorth=float(mo.group(1)) if mo else None
            mo=re.search(r"cells:\s+(\d+)",s); cells=int(mo.group(1)) if mo else None
        row=dict(run_id=rid, block=blk, Re=Re, N=N, scheme=meta.get("scheme"), ddt=meta.get("ddt"), dt=meta.get("dt"), nCorr=meta.get("nCorr"), nOuter=meta.get("nOuter"),
                 ddtCorr=meta.get("ddtCorr"), ptol=meta.get("ptol"), mesh=meta.get("mesh"), momPred=meta.get("momPred"), tag=meta.get("tag"), fo_version=meta.get("fo_version"),
                 ranks=nproc, cells=cells, nonorth_max=nonorth, steps=steps, wall_s=wall, corehours=(wall*nproc/3600. if wall else None), t_end=float(t[-1]),
                 K0=float(c["K"][0]), K_end=float(c["K"][-1]), peak=float(-c["dKdt"][ipk]), t_peak=float(t[ipk]), epsNu_peak=float(c["epsNu"].max()), t_epsNu_peak=float(t[np.argmax(c["epsNu"])]),
                 twoNuOmega_peak=float(c["twoNuOmega"].max()), ref=rn, ref_peak=float(er[m][ir]), ref_t_peak=float(tr[m][ir]), ref_K_end=float(np.interp(t[-1],tr,Kr)),
                 peak_rel_err=float(-c["dKdt"][ipk]/er[m][ir]-1), L2K_rel=float(np.sqrt(np.mean((c["K"]-Ki)**2))/c["K"][0]), L2eps_rel=float(np.sqrt(np.mean((-c["dKdt"]-ei)**2))/er[m][ir]),
                 int_epsNu=I(c["epsNu"]), int_sumE=I(c["sumE"]), int_eNum=I(epsnum), rho=I(epsnum)/I(c["epsNu"]), rho_abs=I(np.abs(epsnum))/I(c["epsNu"]),
                 eNum_over_epsNu_at_peak=float(epsnum[ipk]/c["epsNu"][ipk]),
                 maxRcl=float(c["Rclosure"].max()), maxCum=float(c["cumClosure"].max()), maxRcell=float(c["RcellMax"].max()), steps_Rcl_gt_1e12=int((c["Rclosure"]>1e-12).sum()),
                 CoMax=float(c["CoMax"].max()), dt_min=float(dt.min()), dt_max=float(dt.max()), maskFrac_mean=float(c["maskFrac"].mean()) if "maskFrac" in c else None,
                 max_sumFluxC=float(np.abs(c["sumFluxC"]).max()) if "sumFluxC" in c else None)
        for k in ("eTime","eTimeDiss","eTimeStore","eConv","eMesh","eCont","ePres","eDev","eIter","ePhi"):
            row["int_"+k]=I(c[k]); sgn=-1.0 if k=="eIter" else 1.0
            row["share_"+k]=sgn*I(c[k])/I(epsnum) if I(epsnum)!=0 else None   # shares of e_num (sum to 1; s_iter = -int e_iter / int e_num)
        for k in ("eTimeHalf","eConvHalf","eDiffHalf","ePresHalf","eIterHalf","eConvOwner","eMeshOwner","epsNuOwner"):
            if k in c: row["int_"+k]=I(c[k])
        rows.append(row)
        sel=["t","dt","CoMax","K","dKdt","epsNu","twoNuOmega","eTime","eTimeDiss","eTimeStore","eConv","eMesh","eCont","ePres","eDev","eIter","sumE","eNum","Rclosure","cumClosure","RcellMax","chkPres","sumFluxC"]
        sel=[k for k in sel if k in c]; stride=max(1,len(t)//400)
        arr=np.column_stack([c[k] for k in sel])[::stride]
        if stride>1 and (len(t)-1)%stride: arr=np.vstack([arr, np.array([c[k][-1] for k in sel])])
        np.savetxt(OUT+"/ts_%s.csv"%rid, arr, delimiter=",", header=",".join(sel), comments="", fmt="%.10e")
        # full-resolution closure series for closure figures (Rclosure needs every step)
        np.savetxt(OUT+"/rcl_%s.csv"%rid, np.column_stack([t,c["Rclosure"],c["cumClosure"],c["RcellMax"]]), delimiter=",", header="t,Rclosure,cumClosure,RcellMax", comments="", fmt="%.6e")
keys=sorted(set(k for r in rows for k in r), key=lambda k: (k not in ("run_id","block","Re","N"), k))
w=csv.DictWriter(open(OUT+"/metrics_all.csv","w",newline=""), fieldnames=keys); w.writeheader(); [w.writerow(r) for r in rows]
# V1p differences vs r16
try:
    ref=load(glob.glob(R+"/V1p/*_r16")[0]); out=[ref["t"]]; hdr=["t"]
    for r in (1,32,64,128):
        a=load(glob.glob(R+"/V1p/*_r%d"%r)[0]); n=min(len(a["t"]),len(ref["t"]))
        for k in ("K","dKdt","ePres"): out.append(np.abs(a[k][:n]-ref[k][:n])); hdr.append("d%s_r%d"%(k,r))
    n=min(len(x) for x in out); np.savetxt(OUT+"/v1p_diff.csv", np.column_stack([x[:n] for x in out]), delimiter=",", header=",".join(hdr), comments="", fmt="%.6e")
except Exception as e: print("V1p diff skipped:", e)
print("runs:", len(rows), "→", OUT); print(sorted(set(r["block"] for r in rows)))
