#!/usr/bin/env python3
"""Campaign case generator (INSTRUCTIONS §8): builds a 3D TGV case directory from parameters, with the run_id naming
convention  <block><Re>_h<N>[_<mesh>]_<scheme>_<ddt>_dt<dt>_nC<n>_dc<0|1>_p<tol>[_<tag>].
Usage: makeCase.py --block S --Re 1600 --N 128 --scheme linear --ddt backward --dt 0.025 [--nCorr 2] [--nOuter 1]
       [--ddtCorr 1] [--ptol 1e-10] [--mesh hex|sin|rand|poly] [--momPred 1] [--endTime 20] [--writeInterval 0.5] [--tag x] [--out DIR]"""
import re, argparse, os, shutil, math, json, subprocess, sys
ap = argparse.ArgumentParser()
ap.add_argument("--block", default="S"); ap.add_argument("--Re", type=float, required=True); ap.add_argument("--N", type=int, required=True)
ap.add_argument("--scheme", default="linear"); ap.add_argument("--ddt", default="backward"); ap.add_argument("--dt", type=float, required=True)
ap.add_argument("--nCorr", type=int, default=2); ap.add_argument("--nOuter", type=int, default=1); ap.add_argument("--ddtCorr", type=int, default=1)
ap.add_argument("--ptol", default="1e-10"); ap.add_argument("--mesh", default="hex"); ap.add_argument("--momPred", type=int, default=1)
ap.add_argument("--endTime", type=float, default=20.0); ap.add_argument("--writeInterval", type=float, default=0.5)
ap.add_argument("--tag", default=""); ap.add_argument("--out", default="."); ap.add_argument("--template", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "tgv3D_16"))
a = ap.parse_args()
schemeMap = {"linear": "Gauss linear", "LUST": "Gauss LUST grad(U)", "linearUpwind": "Gauss linearUpwind grad(U)",
             "limitedLinear": "Gauss limitedLinear 1", "cubic": "Gauss cubic", "upwind": "Gauss upwind", "midPoint": "Gauss midPoint"}
ddtMap = {"Euler": "Euler", "backward": "backward", "CN": "CrankNicolson 0.9"}
ddtEntry = ddtMap[a.ddt] + (" 0" if (a.ddtCorr == 0 and a.ddt == "backward") else "")
if a.ddtCorr == 0 and a.ddt != "backward": sys.exit("ddtCorr off is only possible with backward (ddtPhiCoeff)")
dtTag = ("%g" % a.dt).replace("0.", "0").replace(".", "")
rid = "%s%d_h%d%s_%s_%s_dt%s_nC%d_dc%d_p%s%s" % (a.block, int(a.Re), a.N, ("" if a.mesh == "hex" else "_" + a.mesh), a.scheme,
        a.ddt.lower()[:3] if a.ddt != "CN" else "cn", dtTag, a.nCorr, a.ddtCorr, a.ptol.replace("1e-", "e"), ("_" + a.tag if a.tag else ""))
if a.nOuter != 1: rid += "_nO%d" % a.nOuter
if a.momPred == 0: rid += "_noMP"
dst = os.path.join(a.out, rid)
if os.path.exists(dst): sys.exit("exists: " + dst)
shutil.copytree(a.template, dst)
def fd(entry, value, f):
    subprocess.run(["foamDictionary", "-entry", entry, "-set", value, os.path.join(dst, f)], check=True, stdout=subprocess.DEVNULL)
# blockMeshDict: edit the N macro textually (NOT via foamDictionary, which rewrites the file with 6-digit precision -> L = 6.28319 != 2*pi;
# found 2026-09-17: perturbed-mesh cyclics then mismatch (closure 1e-9) and the rand perturbation moves far-plane points)
bm = os.path.join(dst, "system/blockMeshDict"); txt = open(bm).read(); txt2, nsub = re.subn(r"^N\s+\d+;", "N %d;" % a.N, txt, flags=re.M)
assert nsub == 1 and "6.283185307179586" in txt2, "blockMeshDict template: expected 'N <n>;' line and full-precision L"
open(bm, "w").write(txt2)
fd("nu", "%.10g" % (1.0 / a.Re), "constant/transportProperties")
fd("divSchemes/div(phi,U)", schemeMap[a.scheme], "system/fvSchemes")
fd("ddtSchemes/default", ddtEntry, "system/fvSchemes")
fd("deltaT", "%g" % a.dt, "system/controlDict"); fd("endTime", "%g" % a.endTime, "system/controlDict")
fd("writeControl", "timeStep", "system/controlDict"); fd("writeInterval", "%d" % round(a.writeInterval / a.dt), "system/controlDict")
fd("PIMPLE/nCorrectors", str(a.nCorr), "system/fvSolution"); fd("PIMPLE/nOuterCorrectors", str(a.nOuter), "system/fvSolution")
fd("PIMPLE/momentumPredictor", "yes" if a.momPred else "no", "system/fvSolution")
fd("PIMPLE/nNonOrthogonalCorrectors", "0" if a.mesh == "hex" else "1", "system/fvSolution")
fd("solvers/p/tolerance", a.ptol, "system/fvSolution"); fd("solvers/pFinal/tolerance", a.ptol, "system/fvSolution")
fd("solvers/p/solver", "GAMG", "system/fvSolution"); fd("solvers/p/smoother", "GaussSeidel", "system/fvSolution"); fd("solvers/pFinal/solver", "GAMG", "system/fvSolution"); fd("solvers/pFinal/smoother", "GaussSeidel", "system/fvSolution")
vf = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "VERSION")
meta = dict(run_id=rid, fo_version=(open(vf).read().strip() if os.path.exists(vf) else "unversioned"), **vars(a)); json.dump(meta, open(os.path.join(dst, "run_meta.json"), "w"), indent=1)
open(os.path.join(dst, "Allrun.pre"), "w").write("""#!/bin/bash
# pre-processing (serial mesh; fields set IN PARALLEL after decomposition — see INSTRUCTIONS §7)
cd "${0%%/*}" || exit 1
source /usr/lib/openfoam/openfoam2406/etc/bashrc
MESH=%s; N=%d; DX=$(python3 -c "import math; print(2*math.pi/%d)")
foamDictionary -entry writeFormat -set ascii system/controlDict > /dev/null
if [ "$MESH" = "poly" ]; then
  PM=/scratch/soscfd00/tgvAudit/meshes/poly$N/constant/polyMesh
  [ -d $PM ] || { echo "POLY MESH MISSING: $PM"; exit 1; }
  rm -rf constant/polyMesh && cp -r $PM constant/polyMesh && echo "poly mesh copied from $PM" > log.blockMesh
else
  blockMesh > log.blockMesh 2>&1 || exit 1
fi
case $MESH in
  sin)  python3 ../perturbPoints3D.py sin $(python3 -c "print(0.2*$DX)") > log.perturb 2>&1;;
  rand) python3 ../perturbPoints3D.py rand 0.2 $DX > log.perturb 2>&1;;
esac
foamDictionary -entry writeFormat -set binary system/controlDict > /dev/null
checkMesh > log.checkMesh 2>&1
R=run_record_$(date +%%Y%%m%%d); mkdir -p $R && cp -r system constant/transportProperties constant/turbulenceProperties $R/ && foamVersion > $R/foamVersion.txt 2>&1; cp log.checkMesh $R/ 2>/dev/null
echo PRE_OK
""" % (a.mesh, a.N, a.N))
os.chmod(os.path.join(dst, "Allrun.pre"), 0o755)
print(rid)
