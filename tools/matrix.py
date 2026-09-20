#!/usr/bin/env python3
"""Campaign matrix: the list of run specifications passed to makeCase.py.
Blocks: V2 (resolved limit at Re=100), S (screening at Re=280 and 1600), Sint (selected interactions), F (focused 256^3),
T (mesh topology), V1p (parallel reproducibility), Vadapt (adaptive time step), Lchk (box-length check), Vpar (controls).
The rank counts below are those used for the published campaign; adjust RANKS for a different machine.
Usage: matrix.py [--block V2|S|Sint|T|all] [--print]"""
import argparse, json
DT = {32: 0.1, 64: 0.05, 128: 0.025, 256: 0.0125}
RANKS = {32: 8, 64: 16, 128: 64, 256: 128}   # MPI ranks per mesh size
def spec(block, Re, N, scheme="linear", ddt="backward", dt=None, nCorr=2, nOuter=1, ddtCorr=1, ptol="1e-10", mesh="hex", momPred=1, tag="", ranks=None):
    return dict(block=block, Re=Re, N=N, scheme=scheme, ddt=ddt, dt=dt if dt else DT[N], nCorr=nCorr, nOuter=nOuter,
                ddtCorr=ddtCorr, ptol=ptol, mesh=mesh, momPred=momPred, tag=tag, ranks=ranks or RANKS[N])
def build(block):
    runs = []
    if block in ("V2", "all"):
        for N in (32, 64, 128): runs.append(spec("V2", 100, N))
        for sch in ("LUST", "linearUpwind"): runs.append(spec("V2", 100, 64, scheme=sch))
    if block in ("S", "all"):
        for Re in (1600, 280):
            for N in (64, 128): runs.append(spec("S", Re, N))                                    # baseline B0
            for sch in ("LUST", "linearUpwind", "limitedLinear", "cubic"):
                runs.append(spec("S", Re, 128, scheme=sch)); runs.append(spec("S", Re, 64, scheme=sch))
            for ddt in ("Euler", "CN"): runs.append(spec("S", Re, 128, ddt=ddt))
            for dt in (0.0125, 0.05): runs.append(spec("S", Re, 128, dt=dt))
            runs.append(spec("S", Re, 128, nCorr=3)); runs.append(spec("S", Re, 128, nOuter=3))
            runs.append(spec("S", Re, 128, ddtCorr=0)); runs.append(spec("S", Re, 128, momPred=0))
            runs.append(spec("S", Re, 128, ptol="1e-6"))
        runs.append(spec("S", 1600, 64, scheme="upwind", tag="known"))                            # known-answer
    if block in ("Sint", "all"):   # selected interactions, Re=1600, 128^3
        for sch in ("LUST", "linearUpwind"):
            for dt in (0.0125, 0.05): runs.append(spec("Sint", 1600, 128, scheme=sch, dt=dt))
            runs.append(spec("Sint", 1600, 128, scheme=sch, ddt="Euler"))
        for dt in (0.0125, 0.05): runs.append(spec("Sint", 1600, 128, ddtCorr=0, dt=dt))
    if block in ("T", "all"):
        for mesh in ("sin", "rand", "poly"):
            runs.append(spec("T", 1600, 128, mesh=mesh)); runs.append(spec("T", 1600, 128, mesh=mesh, scheme="LUST"))
        runs.append(spec("T", 1600, 128, mesh="sin", scheme="midPoint", tag="emesh0"))
        runs.append(spec("T", 1600, 64, mesh="sin")); runs.append(spec("T", 1600, 64, mesh="poly"))
    if block in ("Vadapt", "all"): runs.append(spec("Vadapt", 1600, 128, tag="adaptive"))
    if block in ("V1p", "all"):   # parallel reproducibility: 64^3 Re=1600 B0 at 1/32/64/128 ranks (16-rank run exists in S)
        for r in (1, 16, 32, 64, 128): runs.append(spec("V1p", 1600, 64, tag="r%d" % r, ranks=r))
    if block in ("Lchk", "all"): runs.append(spec("Lchk", 1600, 64, tag="exactL"))   # box length exactly 2*pi
    if block in ("F", "all"):
        runs.append(spec("F", 1600, 256))                                                          # F1: baseline at 256^3
        runs.append(spec("F", 1600, 256, ddtCorr=0, tag="F2"))                                     # F2: the least dissipative setting of the screening, run to t = 12
    return runs
if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--block", default="all"); ap.add_argument("--print", action="store_true"); a = ap.parse_args()
    runs = build(a.block)
    ch = {N: 0.0 for N in RANKS}
    for r in runs: ch[r["N"]] += 1
    print(json.dumps(runs, indent=0) if a.print else "")
    print("runs:", len(runs), "by N:", {k: int(v) for k, v in ch.items() if v})
