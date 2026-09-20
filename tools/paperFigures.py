#!/usr/bin/env python3
"""Generate the paper's PGFplots figures (PGFPLOTS_RULES.md compliant) from stage1/paperData into paper/figures/<name>/.
Each figure: standalone .tex + its CSV data copied into the subfolder; compiled with pdflatex; symlink figures/<name>.pdf."""
import os, shutil, subprocess, csv, math, re, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); D=os.path.join(HERE,"..","paperData"); FIG=os.path.join(HERE,"..","..","paper","figures")
COL=["blue","red","black","orange","brown","green!20!black","gray"]
# ---- global colour maps: ONE mapping used by every figure, assigned explicitly per series
#      (never through a cycle list, so a panel missing one member keeps the others' colours)
SCHEMECOL={"linear":"blue","LUST":"red","linearUpwind":"black","limitedLinear":"orange","cubic":"brown","upwind":"green!20!black"}
MESHCOL={32:"orange",64:"blue",128:"red",256:"brown"}
REFCOL="black"   # reference (pseudo-spectral) solution: ALWAYS black, dashed, line width=1pt
OFF="dashdot"    # a "ddtCorr off" variant keeps the colour of its mesh and is densely dashdotted
def rid(block,Re,N,scheme="linear",ddt="bac",dt=None,nC=2,dc=1,pe=10,mesh=None,tag=None):
    dts={0.1:"dt01",0.05:"dt005",0.025:"dt0025",0.0125:"dt00125"}; dt=dt or {32:0.1,64:0.05,128:0.025,256:0.0125}[N]
    s="%s%d_h%d"%(block,Re,N)+("_%s"%mesh if mesh else "")+"_%s_%s_%s_nC%d_dc%d_pe%d"%(scheme,ddt,dts[dt],nC,dc,pe)+("_%s"%tag if tag else ""); return s
def ticks(lo,hi,step,dec):
    n=int(round((hi-lo)/step)); vals=[lo+i*step for i in range(n+1)]; f="%%.%df"%dec
    return "xtick={%s}, xticklabels={%s}"%(",".join(f%v for v in vals),",".join("$"+f%v+"$" for v in vals)), vals
def yticks(lo,hi,step,dec):
    n=int(round((hi-lo)/step)); vals=[lo+i*step for i in range(n+1)]; f="%%.%df"%dec
    return "ytick={%s}, yticklabels={%s}"%(",".join(f%v for v in vals),",".join("$"+f%v+"$" for v in vals))
def logticks(axis,lo,hi,every=2):
    e=list(range(int(math.log10(lo)),int(math.log10(hi))+1,every)); return "%stick={%s}"%(axis,",".join("1e%d"%k for k in e))
def figure(name, series, xlabel, ylabel, xlim, ylim, xt, yt, legend_cols=-1, ymode=None, xmode=None, extra="", legend_offset=-0.35, panel=None):
    """series: list of dict(file, x, y (col or expr), color, style ('data'|'ref'|'marks'|'guide'), legend, [csvsrc])"""
    d=os.path.join(FIG,name); os.makedirs(d,exist_ok=True)
    for s in series:
        src=s.get("csvsrc") or os.path.join(D,s["file"])
        if os.path.exists(src): shutil.copy(src, os.path.join(d,os.path.basename(s["file"])))
    grid="major grid style={line width=0.3pt, gray},\n    minor grid style={line width=0.12pt, black!15},\n    minor tick num=1,"
    if ymode=="log" and xmode!="log": grid="major grid style={line width=0.3pt, gray},\n    minor grid style={line width=0.12pt, black!15},\n    minor x tick num=1, minor y tick num=9,"
    if ymode=="log" and xmode=="log": grid="major grid style={line width=0.3pt, gray},\n    minor grid style={line width=0.08pt, black!10},\n    minor tick num=1,"
    opts=["width=7.5cm","height=6cm","xlabel={%s}"%xlabel,"ylabel={%s}"%ylabel,"xmin=%s, xmax=%s"%xlim,"ymin=%s, ymax=%s"%ylim,"grid=both",grid,
          "legend style={\n        at={(0.5,%s)},\n        anchor=north,\n        legend columns=%d,\n        draw=none,\n    }"%(legend_offset,legend_cols),xt,yt,"tick align=inside","tick pos=left","scaled x ticks=false","scaled y ticks=false"]
    if ymode=="log": opts.append("ymode=log"); 
    if xmode=="log": opts.append("xmode=log")
    if extra: opts.append(extra)
    body=""
    for s in series:
        st={"data":"%s, line width=0.6pt"%s["color"],"ref":"%s, line width=1pt, dashed"%s["color"],"marks":"%s, only marks, mark=*, mark size=2.5pt"%s["color"],"omarks":"%s, only marks, mark=o, mark size=2.5pt"%s["color"],"guide":"%s, line width=1pt, dashed"%s["color"],"dotted":"%s, line width=0.6pt, densely dotted"%s["color"],"dash":"%s, line width=0.6pt, densely dashed"%s["color"],"dashdot":"%s, line width=0.6pt, densely dashdotted"%s["color"]}[s["style"]]
        if "domain" in s: body+="\\addplot[%s, domain=%s, samples=2] {%s};\n"%(st,s["domain"],s["expr"])
        else:
            y=s["y"]; ycol=("y expr={%s}"%y) if ("\\thisrow" in y) else ("y=%s"%y)
            body+="\\addplot[%s]\n    table[x=%s, %s, col sep=comma] {%s};\n"%(st,s["x"],ycol,os.path.basename(s["file"]))
        body+="\\addlegendentry{%s}\n\n"%s["legend"] if s.get("legend") else "\\addlegendimage{empty legend}\\addlegendentry{}\n\n" if False else ""
    if panel: body+="\\node[anchor=north west, inner sep=2pt] at (rel axis cs:0.02,0.98) {(%s)};\n"%panel
    tex="\\documentclass[tikz]{standalone}\n\\usepackage{pgfplots}\n\\pgfplotsset{compat=newest}\n\n\\begin{document}\n\\begin{tikzpicture}\n\n\\begin{axis}[\n    "+",\n    ".join(opts)+",\n]\n\n"+body+"\\end{axis}\n\\end{tikzpicture}\n\\end{document}\n"
    open(os.path.join(d,name+".tex"),"w").write(tex)
    r=subprocess.run(["pdflatex","-interaction=nonstopmode",name+".tex"],cwd=d,capture_output=True,text=True)
    ok=os.path.exists(os.path.join(d,name+".pdf")) and r.returncode==0
    link=os.path.join(FIG,name+".pdf")
    if os.path.lexists(link): os.remove(link)
    os.symlink(os.path.join(name,name+".pdf"),link)
    for ext in (".aux",".log"):
        p=os.path.join(d,name+ext); os.path.exists(p) and os.remove(p)
    print("%-22s %s"%(name,"OK" if ok else "FAILED: "+r.stdout[-400:]))
    return ok
T="$t$ $[-]$"; TX="xtick={0,5,10,15,20}, xticklabels={$0$, $5$, $10$, $15$, $20$}"
ts=lambda r:"ts_%s.csv"%r; rc=lambda r:"rcl_%s.csv"%r
def S(file,y,color,legend,style="data",x="t"): return dict(file=file,x=x,y=y,color=color,legend=legend,style=style)
MD="-\\thisrow{dKdt}"
# ---- validation: K(t) and -dK/dt(t) vs reference, per Re
def valfig(Re, runs, refcsv, reflab):
    """runs: list of (runid, legend, colour, style) -- colours from MESHCOL, reference always REFCOL dashed."""
    ser=[S(ts(r),"K",c,l,st) for r,l,c,st in runs]+[S(refcsv,"K",REFCOL,reflab,"ref")]
    figure("valRe%d_K"%Re, ser, T, "$K$ $[-]$", (0,20), (0,0.125), TX, yticks(0,0.125,0.025,3), legend_cols=2 if len(ser)>3 else -1, panel="a")
    ser=[S(ts(r),MD,c,l,st) for r,l,c,st in runs]+[S(refcsv,"eps",REFCOL,reflab,"ref")]
    figure("valRe%d_eps"%Re, ser, T, "$-\\mathrm{d}K/\\mathrm{d}t$ $[-]$", (0,20), (0,0.015), TX, yticks(0,0.015,0.005,3), legend_cols=2 if len(ser)>3 else -1, panel="b")
valfig(100,[(rid("V2",100,32),"$32^3$",MESHCOL[32],"data"),(rid("V2",100,64),"$64^3$",MESHCOL[64],"data"),(rid("V2",100,128),"$128^3$",MESHCOL[128],"data")],"ref_Re100_specN128.csv","Spectral $128^3$")
valfig(280,[(rid("S",280,64),"$64^3$",MESHCOL[64],"data"),(rid("S",280,128),"$128^3$",MESHCOL[128],"data"),(rid("S",280,128,dc=0),"$128^3$, ddtCorr off",MESHCOL[128],OFF)],"ref_Re280_specN192.csv","Spectral $192^3$")
valfig(1600,[(rid("S",1600,64),"$64^3$",MESHCOL[64],"data"),(rid("S",1600,128),"$128^3$",MESHCOL[128],"data"),(rid("F",1600,256),"$256^3$",MESHCOL[256],"data")],"ref_Re1600_UCL512.csv","Spectral $512^3$")
# ---- dissipation split (error cancellation)
r128=rid("S",1600,128); r256=rid("F",1600,256); r256dc0=rid("F",1600,256,dc=0,tag="F2"); r128dc0=rid("S",1600,128,dc=0)
lab3=["$128^3$","$256^3$","$256^3$, ddtCorr off"]
figure("split_minusdKdt",[S(ts(r128),MD,MESHCOL[128],lab3[0]),S(ts(r256),MD,MESHCOL[256],lab3[1]),S(ts(r256dc0),MD,MESHCOL[256],lab3[2],OFF),S("ref_Re1600_UCL512.csv","eps",REFCOL,"Spectral $512^3$","ref")],T,"$-\\mathrm{d}K/\\mathrm{d}t$ $[-]$",(0,20),(0,0.015),TX,yticks(0,0.015,0.005,3),legend_cols=2, panel="a")
figure("split_epsNu",[S(ts(r128),"epsNu",MESHCOL[128],lab3[0]),S(ts(r256),"epsNu",MESHCOL[256],lab3[1]),S(ts(r256dc0),"epsNu",MESHCOL[256],lab3[2],OFF),S("ref_Re1600_UCL512.csv","eps",REFCOL,"Spectral $512^3$","ref")],T,"$\\varepsilon_\\nu$ $[-]$",(0,20),(0,0.015),TX,yticks(0,0.015,0.005,3),legend_cols=2, panel="b")
figure("split_eNum",[S(ts(r128),"eNum",MESHCOL[128],lab3[0]),S(ts(r256),"eNum",MESHCOL[256],lab3[1]),S(ts(r256dc0),"eNum",MESHCOL[256],lab3[2],OFF),S(ts(r128dc0),"eNum",MESHCOL[128],"$128^3$, ddtCorr off",OFF)],T,"$e_\\mathrm{num}$ $[-]$",(0,20),(0,0.006),TX,yticks(0,0.006,0.002,3),legend_cols=2, panel="c")
# ---- signature figures: components vs t
def sig(name, r, withConv, panel=None):
    ser=[S(ts(r),MD,COL[0],"$-\\mathrm{d}K/\\mathrm{d}t$"),S(ts(r),"epsNu",COL[1],"$\\varepsilon_\\nu$"),S(ts(r),"eNum",COL[2],"$e_\\mathrm{num}$"),S(ts(r),"ePres",COL[3],"$e_\\mathrm{pres}$")]
    ser+=[S(ts(r),"eConv",COL[4],"$e_\\mathrm{conv}$")] if withConv else [S(ts(r),"eTime",COL[4],"$e_\\mathrm{time}$")]
    ser+=[S(ts(r),"eIter",COL[5],"$e_\\mathrm{iter}$","dotted"),S(ts(r),"\\thisrow{eMesh}+\\thisrow{eCont}+\\thisrow{eDev}"+("+\\thisrow{eTime}" if withConv else ""),COL[6],"Other")]
    figure(name,ser,T,"$e_i$ $[-]$",(0,20),(-0.005,0.015),TX,yticks(-0.005,0.015,0.005,3),legend_cols=3,legend_offset=-0.32,panel=panel)
sig("sig_128B0",r128,False,"a"); sig("sig_128LUST",rid("S",1600,128,scheme="LUST"),True,"b"); sig("sig_256B0",r256,False,"d")
# ---- ddtCorr: e_pres(t) and rho vs dt
figure("ddtCorr_ePres",[S(ts(r128),"ePres",MESHCOL[128],"$128^3$"),S(ts(r128dc0),"ePres",MESHCOL[128],"$128^3$, ddtCorr off",OFF),S(ts(r256),"ePres",MESHCOL[256],"$256^3$"),S(ts(r256dc0),"ePres",MESHCOL[256],"$256^3$, ddtCorr off",OFF)],T,"$e_\\mathrm{pres}$ $[-]$",(0,20),(-0.001,0.006),TX,yticks(-0.001,0.006,0.001,3),legend_cols=2, panel="a")
rows=list(csv.DictReader(open(os.path.join(D,"metrics_all.csv"))))
def m(cond): return [r for r in rows if cond(r)]
def wcsv(name,hdr,data):
    p=os.path.join(D,name); w=csv.writer(open(p,"w",newline="")); w.writerow(hdr); [w.writerow(x) for x in data]; return name
sel=lambda dc: sorted([(float(r["dt"]),float(r["rho"])) for r in rows if r["Re"]=="1600" and r["N"]=="128" and r["scheme"]=="linear" and r["ddt"]=="backward" and r["nCorr"]=="2" and r["nOuter"]=="1" and r["ptol"]=="1e-10" and r["mesh"]=="hex" and r["momPred"]=="1" and r["ddtCorr"]==str(dc) and r["block"] in ("S","Sint")])
f1=wcsv("rho_dc1.csv",["dt","rho"],sel(1)); f0=wcsv("rho_dc0.csv",["dt","rho"],sel(0)); print("rho points:",sel(1),sel(0))
figure("ddtCorr_rho",[S(f1,"rho",MESHCOL[128],"$128^3$, ddtCorr on","marks",x="dt"),S(f0,"rho",MESHCOL[128],"$128^3$, ddtCorr off","omarks",x="dt")],"$\\Delta t$ $[-]$","$\\rho$ $[-]$",(0,0.06),(0,0.5),"xtick={0,0.01,0.02,0.03,0.04,0.05,0.06}, xticklabels={$0$, $0.01$, $0.02$, $0.03$, $0.04$, $0.05$, $0.06$}",yticks(0,0.5,0.1,1), panel="b")
# ---- schemes
sch=lambda s: rid("S",1600,128,scheme=s)
figure("schemes_eConv",[S(ts(sch("LUST")),"eConv",SCHEMECOL["LUST"],"LUST"),S(ts(sch("linearUpwind")),"eConv",SCHEMECOL["linearUpwind"],"linearUpwind"),S(ts(sch("limitedLinear")),"eConv",SCHEMECOL["limitedLinear"],"limitedLinear"),S(ts(sch("cubic")),"eConv",SCHEMECOL["cubic"],"cubic")],T,"$e_\\mathrm{conv}$ $[-]$",(0,20),(-0.002,0.008),TX,yticks(-0.002,0.008,0.002,3),legend_cols=2, panel="a")
figure("schemes_ePres",[S(ts(r128),"ePres",SCHEMECOL["linear"],"linear"),S(ts(sch("LUST")),"ePres",SCHEMECOL["LUST"],"LUST"),S(ts(sch("linearUpwind")),"ePres",SCHEMECOL["linearUpwind"],"linearUpwind"),S(ts(sch("limitedLinear")),"ePres",SCHEMECOL["limitedLinear"],"limitedLinear"),S(ts(sch("cubic")),"ePres",SCHEMECOL["cubic"],"cubic")],T,"$e_\\mathrm{pres}$ $[-]$",(0,20),(-0.002,0.008),TX,yticks(-0.002,0.008,0.002,3),legend_cols=2,legend_offset=-0.32, panel="b")
# ---- closure (semilog)
RCL="max(abs(\\thisrow{Rclosure}),1e-16)"
figure("closure_Rcl",[S(rc(r128),RCL,COL[0],"Hex $128^3$"),S(rc(rid("T",1600,128,mesh="sin")),RCL,COL[1],"Sinusoidal"),S(rc(rid("T",1600,128,mesh="rand")),RCL,COL[2],"Random"),S(rc(rid("T",1600,128,mesh="poly")),RCL,COL[3],"Polyhedral"),S(rc(rid("Vadapt",1600,128,tag="adaptive")),RCL,COL[4],"Adaptive $\\Delta t$","dash"),S(rc(r256),RCL,COL[5],"Hex $256^3$","dotted")],T,"$R_\\mathrm{closure}$ $[-]$",(0,20),(1e-16,1e-6),TX,logticks("y",1e-16,1e-6,2),legend_cols=2,ymode="log",legend_offset=-0.32)
rr=rid("T",1600,128,mesh="rand")
_tsr=list(csv.DictReader(open(os.path.join(D,ts(rr)))))
fRat=wcsv("rand_fluxRatio.csv",["t","fluxRatio"],[[x["t"],abs(float(x["sumFluxC"]))/max(abs(float(x["dKdt"])),1e-300)] for x in _tsr])
figure("rand_transient",[S(rc(rr),RCL,COL[0],"$R_\\mathrm{closure}$"),S(fRat,"max(\\thisrow{fluxRatio},1e-18)",COL[1],"$|\\sum_f F^C_f|/|\\mathrm{d}K/\\mathrm{d}t|$"),S(rc(rr),"max(abs(\\thisrow{RcellMax}),1e-16)",COL[2],"$R_\\mathrm{cell,max}$")],T,"residual $[-]$",(0,20),(1e-18,1e-6),TX,logticks("y",1e-18,1e-6,2),legend_cols=1,ymode="log",legend_offset=-0.32)
# ---- V1p
V1P="max(\\thisrow{%s},1e-15)"
figure("v1p_diff",[S("v1p_diff.csv",V1P%"dK_r1",COL[0],"1 rank"),S("v1p_diff.csv",V1P%"dK_r32",COL[1],"32 ranks"),S("v1p_diff.csv",V1P%"dK_r64",COL[2],"64 ranks"),S("v1p_diff.csv",V1P%"dK_r128",COL[3],"128 ranks")],T,"$|K-K_{16}|$ $[-]$",(0,20),(1e-15,1e-9),TX,logticks("y",1e-15,1e-9,2),legend_cols=2,ymode="log",legend_offset=-0.32)
# ---- V2 resolved limit: sumE/epsNu
figure("v2_ratio",[S(ts(rid("V2",100,32)),"\\thisrow{eNum}/\\thisrow{epsNu}",MESHCOL[32],"$32^3$"),S(ts(rid("V2",100,64)),"\\thisrow{eNum}/\\thisrow{epsNu}",MESHCOL[64],"$64^3$"),S(ts(rid("V2",100,128)),"\\thisrow{eNum}/\\thisrow{epsNu}",MESHCOL[128],"$128^3$")],T,"$e_\\mathrm{num}/\\varepsilon_\\nu$ $[-]$",(0,20),(0,0.15),TX,yticks(0,0.15,0.05,2))
# ---- V1 (2D) order plots from V1_summary.txt
rows1={}
for line in open(os.path.join(HERE,"..","results","V1_summary.txt")):
    mo=re.match(r"(\S+)\s+N=\s*(\d+) dt=\s*(\S+) ddt=(.+?)\s+mP=.*L2=(\S+).*eTime=(\S+) .*ePres=(\S+) .*eIter=(\S+)",line)
    if mo: rows1[mo.group(1)]=dict(N=int(mo.group(2)),dt=float(mo.group(3)),ddt=mo.group(4).strip(),L2=float(mo.group(5)),eTime=abs(float(mo.group(6))),ePres=abs(float(mo.group(7))),eIter=abs(float(mo.group(8))))
mesh=[rows1["mesh_N%d"%N] for N in (32,64,128,256,512)]
wcsv("v1_mesh.csv",["h","L2","eTime","ePres","eIter"],[[2*math.pi/r["N"],r["L2"],r["eTime"],r["ePres"],r["eIter"]] for r in mesh])
h0=2*math.pi/32
figure("v1_mesh",[S("v1_mesh.csv","L2",COL[0],"$\\|u-u_\\mathrm{exact}\\|_2$","marks",x="h"),S("v1_mesh.csv","ePres",COL[1],"$|\\int e_\\mathrm{pres}\\,\\mathrm{d}t|$","marks",x="h"),S("v1_mesh.csv","eTime",COL[2],"$|\\int e_\\mathrm{time}\\,\\mathrm{d}t|$","marks",x="h"),
    dict(file="v1_mesh.csv",style="guide",color=COL[3],legend="$h^2$",domain="0.01:0.2",expr="%g*(x/%g)^2"%(mesh[0]["L2"],h0)),dict(file="v1_mesh.csv",style="guide",color=COL[4],legend="$h^4$",domain="0.01:0.2",expr="%g*(x/%g)^4"%(mesh[0]["ePres"],h0))],
    "$h$ $[-]$","Error $[-]$",(0.01,0.2),(1e-7,1e-1),"xtick={0.01,0.02,0.05,0.1,0.2}, xticklabels={$0.01$, $0.02$, $0.05$, $0.1$, $0.2$}",logticks("y",1e-7,1e-1,1),legend_cols=2,ymode="log",xmode="log",legend_offset=-0.32, panel="a")
dts=[0.05,0.025,0.0125,0.00625]
for tag,lab in (("Euler","Euler"),("backward","backward"),("CrankNicolson09","CN")):
    wcsv("v1_dt_%s.csv"%tag,["dt","eTime","L2"],[[rows1["dt%g_%s"%(d,tag)]["dt"],rows1["dt%g_%s"%(d,tag)]["eTime"],rows1["dt%g_%s"%(d,tag)]["L2"]] for d in dts if "dt%g_%s"%(d,tag) in rows1])
e0=rows1["dt0.05_Euler"]["eTime"]; b0=rows1["dt0.05_backward"]["eTime"]
figure("v1_dt",[S("v1_dt_Euler.csv","eTime",COL[0],"Euler","marks",x="dt"),S("v1_dt_backward.csv","eTime",COL[1],"backward","marks",x="dt"),S("v1_dt_CrankNicolson09.csv","eTime",COL[2],"CN $\\psi=0.9$","marks",x="dt"),
    dict(file="v1_dt_Euler.csv",style="guide",color=COL[3],legend="$\\Delta t^1$",domain="0.005:0.06",expr="%g*(x/0.05)"%e0),dict(file="v1_dt_Euler.csv",style="guide",color=COL[4],legend="$\\Delta t^2$",domain="0.005:0.06",expr="%g*(x/0.05)^2"%b0)],
    "$\\Delta t$ $[-]$","$|\\int e_\\mathrm{time}\\,\\mathrm{d}t|$ $[-]$",(0.005,0.06),(1e-7,1e-3),"xtick={0.005,0.01,0.02,0.05}, xticklabels={$0.005$, $0.01$, $0.02$, $0.05$}",logticks("y",1e-7,1e-3,1),legend_cols=3,ymode="log",xmode="log",legend_offset=-0.32, panel="b")
# ---- adaptive dt
rv=rid("Vadapt",1600,128,tag="adaptive")
figure("vadapt_dt",[S(ts(rv),"dt",COL[0],"$\\Delta t$")],T,"$\\Delta t$ $[-]$",(0,20),(0.010,0.035),TX,yticks(0.010,0.035,0.005,3), panel="a")
figure("vadapt_Co",[S(ts(rv),"CoMax",COL[0],"$\\mathrm{Co}_\\mathrm{max}$"),S(ts(r128),"CoMax",COL[1],"$\\mathrm{Co}_\\mathrm{max}$, fixed $\\Delta t$")],T,"$\\mathrm{Co}_\\mathrm{max}$ $[-]$",(0,20),(0,1.2),TX,yticks(0,1.2,0.2,1), panel="b")
# ---- Pareto: |peak error| vs core-hours (Re=1600, all runs with finite corehours)
pts={64:[],128:[],256:[]}
for r in rows:
    if r["Re"]=="1600" and r["corehours"] and r["N"] in ("64","128","256") and r["block"] in ("S","Sint","F","T"): pts[int(r["N"])].append((float(r["corehours"]),100*abs(float(r["peak_rel_err"]))))
for N in pts: wcsv("pareto_%d.csv"%N,["coreh","err"],sorted(pts[N]))
figure("pareto",[S("pareto_64.csv","err",MESHCOL[64],"$64^3$","marks",x="coreh"),S("pareto_128.csv","err",MESHCOL[128],"$128^3$","marks",x="coreh"),S("pareto_256.csv","err",MESHCOL[256],"$256^3$","marks",x="coreh")],"$t_\\mathrm{CPU}$ $[\\mathrm{core\\,h}]$","$|\\Delta\\varepsilon_\\mathrm{peak}|$ $[\\%]$",(0.1,10000),(0,15),"xtick={1e-1,1e0,1e1,1e2,1e3,1e4}",yticks(0,15,5,0),xmode="log",legend_cols=3)
print("done")
# ---- V3: audit vs built-in function objects and the classical enstrophy-based estimate (64^3 Re=1600 B0)
figure("v3_compare",[S("v3_compare.csv","eNum_FO",COL[0],"$e_\\mathrm{num}$ (audit)"),S("v3_compare.csv","eNum_classical",COL[1],"$-\\mathrm{d}K/\\mathrm{d}t-2\\nu\\Omega$"),S("v3_compare.csv","epsNu",COL[2],"$\\varepsilon_\\nu$"),S("v3_compare.csv","twoNuOmega",COL[3],"$2\\nu\\Omega$")],T,"rate $[-]$",(0,20),(0,0.012),TX,yticks(0,0.012,0.002,3),legend_cols=2)

# ---- alignment pass (2026-09-20): panels placed side by side in the manuscript must have identical PDF heights, otherwise the
# baseline alignment of \includegraphics shifts the axis frames when the legends have different numbers of rows. For every row
# of panels the shorter PDFs get their legend moved down (larger gap between axis and legend) until the heights match.
import re as _re
ROWS = [("v1_mesh", "v1_dt"), ("vadapt_dt", "vadapt_Co"), ("valRe100_K", "valRe100_eps"), ("valRe280_K", "valRe280_eps"),
        ("valRe1600_K", "valRe1600_eps"), ("sig_128B0", "sig_128LUST"), ("schemes_eConv", "schemes_ePres"), ("ddtCorr_ePres", "ddtCorr_rho"),
        ("split_minusdKdt", "split_epsNu"), ("split_eNum", "sig_256B0")]
def _pdfh(n):
    out = subprocess.run(["pdfinfo", os.path.join(FIG, n, n + ".pdf")], capture_output=True, text=True).stdout
    return float(_re.search(r"Page size:\s+[\d.]+ x ([\d.]+)", out).group(1))
def _compile(n):
    d = os.path.join(FIG, n); r = subprocess.run(["pdflatex", "-interaction=nonstopmode", n + ".tex"], cwd=d, capture_output=True, text=True)
    for ext in (".aux", ".log"):
        p = os.path.join(d, n + ext); os.path.exists(p) and os.remove(p)
    return r.returncode == 0
def _shift_legend(n, dpt, factor):
    texf = os.path.join(FIG, n, n + ".tex"); s = open(texf).read(); m = _re.search(r"at=\{\(0\.5,(-?[\d.]+)\)\}", s)
    off = float(m.group(1)); new = off - dpt * factor
    open(texf, "w").write(s.replace(m.group(0), "at={(0.5,%.5f)}" % new, 1)); _compile(n); return new
for row in ROWS:
    factor = {n: 1.0 / (4.5 * 72.27 / 2.54) for n in row}   # initial guess: 1 rel unit ~ 4.5 cm axis box; refined from the observed response
    for it in range(4):
        hs = {n: _pdfh(n) for n in row}; hmax = max(hs.values()); moved = False
        for n in row:
            dh = hmax - hs[n]
            if dh > 0.3:
                moved = True; h0 = hs[n]; _shift_legend(n, dh, factor[n]); h1 = _pdfh(n)
                if abs(h1 - h0) > 0.05: factor[n] = factor[n] * dh / (h1 - h0)   # calibrate the rel-unit-per-pt response
        if not moved: break
    print("aligned %-32s heights: %s" % ("/".join(row), ", ".join("%.1f" % _pdfh(n) for n in row)))
