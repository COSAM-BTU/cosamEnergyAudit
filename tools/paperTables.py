#!/usr/bin/env python3
"""Generate the LaTeX tables of the paper (tabular* with booktabs) from the extracted metrics_all.csv.
Usage: paperTables.py [--data <csv dir>] [--out <table dir>]   (defaults: ../paperData and ../tables)"""
import csv, os, math
import argparse
_ap=argparse.ArgumentParser(description="Generate the LaTeX tables of the paper.")
_ap.add_argument("--data", default=None, help="directory holding metrics_all.csv (default: ../paperData)")
_ap.add_argument("--out", default=None, help="output directory for the tables (default: ../tables)")
_a=_ap.parse_args()
HERE=os.path.dirname(os.path.abspath(__file__))
D=os.path.abspath(_a.data) if _a.data else os.path.join(HERE,"..","paperData")
OUT=os.path.abspath(_a.out) if _a.out else os.path.join(HERE,"..","tables"); os.makedirs(OUT,exist_ok=True)
rows=list(csv.DictReader(open(os.path.join(D,"metrics_all.csv"))))
f=lambda x,d=3: ("%%.%df"%d)%float(x)
pct=lambda x: "%+.1f"%(100*float(x))
sci=lambda x: ("$%s\\times10^{%d}$"%(("%.0f"%(float(x)/10**math.floor(math.log10(abs(float(x)))))), math.floor(math.log10(abs(float(x)))))) if float(x)!=0 else "0"
def sciE(x):
    v=float(x)
    if v==0: return "0"
    e=math.floor(math.log10(abs(v))); m=v/10**e
    if round(m,1)>=10: m/=10; e+=1
    return "$%.1f\\times10^{%d}$"%(m,e)
def tab(name, caption, label, colspec, header, body, wide=False, size="\\footnotesize", notes=None):
    env="table*" if wide else "table"; width="\\textwidth" if wide else "\\columnwidth"
    s="\\begin{%s}[htbp]\n\\centering\n%s\n\\caption{%s}\n\\label{%s}\n\\begin{tabular*}{%s}{@{\\extracolsep{\\fill}}%s@{}}\n\\toprule\n%s \\\\\n\\midrule\n"%(env,size,caption,label,width,colspec,header)
    s+="".join(r+" \\\\\n" for r in body)+"\\bottomrule\n\\end{tabular*}\n"
    if notes: s+="\\smallskip\n\n{%s %s}\n"%(size,notes)
    s+="\\end{%s}\n"%env; open(os.path.join(OUT,name+".tex"),"w").write(s); print("wrote",name, len(body),"rows")
def get(**kw):
    out=[r for r in rows if all(str(r.get(k))==str(v) for k,v in kw.items())]; return out
def one(**kw):
    o=get(**kw); assert len(o)==1, (kw,len(o)); return o[0]
def label(r):
    s=[]
    if r["scheme"]!="linear": s.append(r["scheme"])
    if r["ddt"]!="backward": s.append({"Euler":"Euler","CN":"CN $\\psi{=}0.9$"}.get(r["ddt"],r["ddt"]))
    dt0={"64":0.05,"128":0.025,"256":0.0125,"32":0.1}[r["N"]]
    if abs(float(r["dt"])-dt0)>1e-12: s.append("$\\Delta t{=}%g$"%float(r["dt"]))
    if r["nCorr"]!="2": s.append("nCorr %s"%r["nCorr"])
    if r["nOuter"]!="1": s.append("nOuter %s"%r["nOuter"])
    if r["ddtCorr"]!="1": s.append("ddtCorr off")
    if r["momPred"]!="1": s.append("No predictor")
    if r["ptol"]!="1e-10": s.append("$p$ tolerance $10^{-6}$")
    if r["mesh"]!="hex": s.append(r["mesh"])
    return ", ".join(s) if s else "B0"
# --- validation
val=[("100","V2",32),("100","V2",64),("100","V2",128),("280","S",64),("280","S",128),("1600","S",64),("1600","S",128),("1600","F",256)]
body=[]
for Re,blk,N in val:
    r=[x for x in get(Re=Re,block=blk,N=str(N),scheme="linear",ddt="backward",ddtCorr="1",mesh="hex",nCorr="2",nOuter="1",ptol="1e-10",momPred="1") if abs(float(x["dt"])-{32:0.1,64:0.05,128:0.025,256:0.0125}[N])<1e-12][0]
    body.append("%s & $%d^3$ & %s & %s & %s [%s] & %s & %s & %s & %s [%s]"%(Re if N in (32,) or (Re=="280" and N==64) or (Re=="1600" and N==64) else "", N, f(r["peak"],5), pct(r["peak_rel_err"]), f(r["t_peak"],2), f(r["ref_t_peak"],2), pct(float(r["epsNu_peak"])/float(r["ref_peak"])-1), pct(float(r["twoNuOmega_peak"])/float(r["ref_peak"])-1), sciE(r["L2K_rel"]), f(r["K_end"],6), f(r["ref_K_end"],6)))
tab("tab_validation","Validation of the baseline setting B0 against the spectral references (Re$\\,=\\,$100: own $128^3$; Re$\\,=\\,$280: own $192^3$; Re$\\,=\\,$1600: workshop $512^3$). Peak of the total dissipation rate $-\\mathrm{d}K/\\mathrm{d}t$ and its relative error $\\Delta_\\mathrm{peak}$, peak time (reference in brackets), relative errors of the peaks of the resolved viscous dissipation $\\varepsilon_\\nu$ and of the enstrophy-based dissipation $2\\nu\\Omega$ with respect to the reference peak, root-mean-square difference over $0 \\le t \\le 20$ of the kinetic energy normalized by $K(0)$, and the final kinetic energy (reference in brackets). All quantities are non-dimensional.",
    "tab:validation","llrrrrrrr","Re & Mesh & Peak & $\\Delta_\\mathrm{peak}$ (\\%) & $t_\\mathrm{peak}$ & $\\varepsilon_\\nu$ (\\%) & $2\\nu\\Omega$ (\\%) & $L_2(K)/K_0$ & $K(20)$",body,wide=True)
# --- screening tables
def screening(Re,N,blocks,name,cap,lab):
    rs=[r for r in rows if r["Re"]==Re and r["N"]==str(N) and r["block"] in blocks and r["mesh"]=="hex"]
    rs.sort(key=lambda r: float(r["rho"]))
    body=[]
    for r in rs:
        ill=abs(float(r["rho"]))<0.01   # shares undefined when the net numerical dissipation vanishes
        sh=lambda k: "--" if ill else f(r[k],2)
        body.append("%s & %s & %s & %s & %s & %s & %s & %s & %s"%(label(r), f(r["rho"],3), sh("share_eTime"), sh("share_eConv"), sh("share_ePres"), sh("share_eIter"), pct(r["peak_rel_err"]), f(r["t_peak"],2), sciE(r["maxRcl"])))
    tab(name,cap,lab,"lrrrrrrrr","Setting & $\\rho$ & $s_\\mathrm{time}$ & $s_\\mathrm{conv}$ & $s_\\mathrm{pres}$ & $s_\\mathrm{iter}$ & $\\Delta\\varepsilon_\\mathrm{peak}$ (\\%) & $t_\\mathrm{peak}$ & $\\max R_\\mathrm{closure}$",body,wide=True,
        notes="$\\rho = \\int_0^{20} e_\\mathrm{num}\\,\\mathrm{d}t / \\int_0^{20}\\varepsilon_\\nu\\,\\mathrm{d}t$ with $e_\\mathrm{num} = -\\mathrm{d}K/\\mathrm{d}t - \\varepsilon_\\nu$; shares $s_i = \\int e_i\\,\\mathrm{d}t/\\int e_\\mathrm{num}\\,\\mathrm{d}t$ and $s_\\mathrm{iter} = -\\int e_\\mathrm{iter}\\,\\mathrm{d}t/\\int e_\\mathrm{num}\\,\\mathrm{d}t$, so that the shares of all terms sum to one and a negative share adds energy (shares are not reported when $|\\rho| < 0.01$, where the net numerical dissipation vanishes and individual terms cancel; the mesh, continuity and explicit-stress shares are below 0.01 in every row and are not listed, so the tabulated shares need not sum exactly to one); $\\Delta\\varepsilon_\\mathrm{peak}$ is the relative error of the peak of $-\\mathrm{d}K/\\mathrm{d}t$ with respect to the reference. B0: linear, backward, $\\Delta t = %s$, nCorrectors 2, ddtCorr on, $p$ tolerance $10^{-10}$. nCorr and nOuter are the numbers of pressure and outer correctors per step; the $p$ tolerance is the absolute tolerance of the pressure solver ($10^{-10}$ in B0)."%{128:"0.025",64:"0.05"}[N])
screening("1600",128,("S","Sint"),"tab_screening1600","Screening at Re$\\,=\\,$1600 on the $128^3$ hexahedral mesh: numerical-to-viscous dissipation ratio $\\rho$, shares of the numerical dissipation carried by the time scheme, the convective scheme, the pressure--velocity coupling and the algebraic residual, and the peak error against the workshop $512^3$ reference. Rows sorted by $\\rho$.","tab:screening1600")
screening("280",128,("S",),"tab_screening280","Screening at Re$\\,=\\,$280 on the $128^3$ hexahedral mesh (reference: own pseudo-spectral $192^3$). Same quantities as Table~\\ref{tab:screening1600}.","tab:screening280")
rs=[r for r in rows if r["N"]=="64" and r["block"]=="S" and r["mesh"]=="hex"]; rs.sort(key=lambda r:(r["Re"],float(r["rho"])))
body=["%s & %s & %s & %s & %s & %s & %s & %s & %s & %s"%(r["Re"], label(r), f(r["rho"],3), f(r["share_eTime"],2), f(r["share_eConv"],2), f(r["share_ePres"],2), f(r["share_eIter"],2), pct(r["peak_rel_err"]), f(r["t_peak"],2), sciE(r["maxRcl"])) for r in rs]
tab("tab_screening64","Convective schemes on the $64^3$ mesh at Re$\\,=\\,$280 and 1600 (B0 otherwise; upwind is the known-answer test).","tab:screening64","llrrrrrrrr","Re & Setting & $\\rho$ & $s_\\mathrm{time}$ & $s_\\mathrm{conv}$ & $s_\\mathrm{pres}$ & $s_\\mathrm{iter}$ & $\\Delta\\varepsilon_\\mathrm{peak}$ (\\%) & $t_\\mathrm{peak}$ & $\\max R_\\mathrm{closure}$",body,wide=True,notes="Shares as in Table~\\ref{tab:screening1600}; for the cubic scheme the algebraic term adds energy ($s_\\mathrm{iter} < 0$), which is why $s_\\mathrm{pres}$ exceeds one.")
# --- topology
body=[]
for r in sorted([x for x in rows if x["Re"]=="1600" and x["N"]=="128" and (x["block"]=="T" or (x["block"]=="S" and x["scheme"] in ("linear","LUST") and label(x) in ("B0","LUST")))], key=lambda x:(x["mesh"]!="hex",{"hex":0,"sin":1,"rand":2,"poly":3}[x["mesh"]],x["scheme"])):
    body.append("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s"%({"hex":"Hexahedral","sin":"Sinusoidal","rand":"Random","poly":"Polyhedral"}[r["mesh"]], f(int(r["cells"])/1e6,2), f(r["nonorth_max"],1), r["scheme"], f(r["rho"],3), f(r["share_ePres"],2), f(r["share_eMesh"],3), f(r["share_eConv"],2), pct(r["peak_rel_err"]), sciE(r["maxRcl"])))
for r in sorted([x for x in rows if x["Re"]=="1600" and x["N"]=="64" and x["block"] in ("T",) or (x["Re"]=="1600" and x["N"]=="64" and x["block"]=="S" and label(x)=="B0")], key=lambda x:{"hex":0,"sin":1,"rand":2,"poly":3}[x["mesh"]]):
    body.append("%s ($64^3$-eq.) & %s & %s & %s & %s & %s & %s & %s & %s & %s"%({"hex":"Hexahedral","sin":"Sinusoidal","rand":"Random","poly":"Polyhedral"}[r["mesh"]], f(int(r["cells"])/1e6,3), f(r["nonorth_max"],1), r["scheme"], f(r["rho"],3), f(r["share_ePres"],2), f(r["share_eMesh"],3), f(r["share_eConv"],2), pct(r["peak_rel_err"]), sciE(r["maxRcl"])))
tab("tab_topology","Mesh topology at Re$\\,=\\,$1600 (T block; $\\Delta t = 0.025$ on the $128^3$-equivalent meshes, $0.05$ on the $64^3$-equivalent ones): number of cells, maximum non-orthogonality from \\texttt{checkMesh}, dissipation ratio $\\rho$, shares of $e_\\mathrm{pres}$, $e_\\mathrm{mesh}$ and $e_\\mathrm{conv}$, peak error and maximum closure residual (the $10^{-8}$--$10^{-7}$ values on the random and polyhedral meshes are the transient coupled-patch flux effect discussed in Sections~\\ref{sec:transient} and \\ref{sec:parallel}).","tab:topology","lrrlrrrrrr","Mesh & Cells ($10^6$) & Non-orth. (deg) & Scheme & $\\rho$ & $s_\\mathrm{pres}$ & $s_\\mathrm{mesh}$ & $s_\\mathrm{conv}$ & $\\Delta\\varepsilon_\\mathrm{peak}$ (\\%) & $\\max R_\\mathrm{closure}$",body,wide=True,notes="The hexahedral rows are the corresponding runs of the screening block, repeated for comparison; the T block proper comprises the nine runs on the perturbed and polyhedral meshes.")
# --- convention sensitivity
body=[]
for r in [one(run_id="S1600_h128_linear_bac_dt0025_nC2_dc1_pe10"), one(run_id="S1600_h128_LUST_bac_dt0025_nC2_dc1_pe10"), one(run_id="T1600_h128_poly_linear_bac_dt0025_nC2_dc1_pe10"), one(run_id="S1600_h128_linear_eul_dt0025_nC2_dc1_pe10")]:
    E=float(r["int_epsNu"]); g=lambda k: f(float(r[k])/E,4) if r.get(k) not in (None,"") else "--"
    body.append("%s, %s, %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s"%({"hex":"Hex","poly":"Poly"}[r["mesh"]], "linear" if r["scheme"]=="linear" else r["scheme"], {"backward":"backward","Euler":"Euler"}.get(r["ddt"],r["ddt"]), g("int_eTime"), g("int_eTimeHalf"), g("int_eConv"), g("int_eConvHalf"), g("int_ePres"), g("int_ePresHalf"), g("int_eIter"), g("int_eIterHalf"), g("int_eMesh"), g("int_eMeshOwner")))
tab("tab_convention","Sensitivity of the attribution to the convention (Re$\\,=\\,$1600, $128^3$, $\\Delta t = 0.025$): time integrals of the terms normalized by $\\int\\varepsilon_\\nu\\,\\mathrm{d}t$ for the primary convention (multiplier $\\mathbf{U}^{n+1}$, half--half face-to-cell distribution) and for the alternatives computed in the same runs (multiplier $\\mathbf{U}^{n+1/2}$, superscript $\\frac12$; owner-side distribution, superscript own). For the Euler scheme (last row) the $\\mathbf{U}^{n+1/2}$ multiplier removes the time term exactly.","tab:convention","lrrrrrrrrrr","Case & $e_\\mathrm{time}$ & $e^{\\frac12}_\\mathrm{time}$ & $e_\\mathrm{conv}$ & $e^{\\frac12}_\\mathrm{conv}$ & $e_\\mathrm{pres}$ & $e^{\\frac12}_\\mathrm{pres}$ & $e_\\mathrm{iter}$ & $e^{\\frac12}_\\mathrm{iter}$ & $e_\\mathrm{mesh}$ & $e^\\mathrm{own}_\\mathrm{mesh}$",body,wide=True)
# --- V1p
body=[]; ref=one(run_id="V1p1600_h64_linear_bac_dt005_nC2_dc1_pe10_r16")
import numpy as np
dd=np.loadtxt(os.path.join(D,"v1p_diff.csv"),delimiter=",",skiprows=1); hdr=open(os.path.join(D,"v1p_diff.csv")).readline().strip().split(",")
for rk in (1,16,32,64,128):
    r=one(run_id="V1p1600_h64_linear_bac_dt005_nC2_dc1_pe10_r%d"%rk); dK="--" if rk==16 else sciE(dd[:,hdr.index("dK_r%d"%rk)].max()/0.125)
    body.append("%d & %.10f & %s & %s & %s"%(rk, float(r["K_end"]), dK, sciE(r["maxRcl"]), f(float(r["wall_s"])/60,1)))
tab("tab_v1p","Parallel reproducibility (V1p): $64^3$ B0 run at Re$\\,=\\,$1600 on 1 to 128 ranks with the fields initialized in parallel. Final kinetic energy, maximum deviation $\\Delta K_{\\max}$ of $K(t)$ from the 16-rank run normalized by $K(0)$, maximum closure residual and wall time.","tab:v1p","rrrrr","Ranks & $K(20)$ & $\\Delta K_{\\max}/K_0$ & $\\max R_\\mathrm{closure}$ & Wall (min)",body)
# --- campaign summary
blocks=[("V2","Resolved limit, Re$\\,=\\,$100 ($32^3$--$128^3$, two schemes)"),("S","Screening, Re$\\,=\\,$280 and 1600 ($64^3$, $128^3$)"),("Sint","Interactions, Re$\\,=\\,$1600 ($128^3$)"),("F","Focused, Re$\\,=\\,$1600 ($256^3$)"),("T","Mesh topology, Re$\\,=\\,$1600"),("V1p","Parallel reproducibility ($64^3$)"),("Vadapt","Adaptive time step ($128^3$)"),("Lchk","Box-length check ($64^3$)"),("Vpar","Decomposition and corrector controls ($64^3$ and $128^3$ random and polyhedral meshes)")]
body=[]; tot=0; nt=0
for b,desc in blocks:
    rs=get(block=b); ch=round(sum(float(r["corehours"]) for r in rs if r["corehours"])); tot+=ch; nt+=len(rs)
    body.append("%s & %s & %d & %s & %.0f"%(b,desc,len(rs),"/".join(sorted(set(str(max(1,int(r["ranks"]))) for r in rs),key=int)),ch))
body.append("\\midrule\nTotal & & %d & & %.0f"%(nt,tot))
tab("tab_campaign","Simulation campaign: blocks, number of runs, MPI ranks (one rank per physical core) and physical core-hours (solver wall time $\\times$ ranks; the pseudo-spectral references, the timing runs, the two-dimensional verification runs and the ten-step instrumented-solver (V4) and built-in-function-object (V3) checks are not included; the total is the sum of the rounded block values).","tab:campaign","llrrr","Block & Content & Runs & Ranks & Core-hours",body,wide=True)
print("campaign core-hours (solver only):", round(tot))
