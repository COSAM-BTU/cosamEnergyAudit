#!/usr/bin/env python3
"""Campaign matrix (INSTRUCTIONS §8) → list of run specifications for makeCase.py.
Blocks: V2 (Re=100 resolved limit), S (screening Re=280/1600), F (focused 256^3, filled after screening), T (topology).
Ranks per mesh (timing test 2026-09-15: 128^3 at 64 ranks = 3.90 s/step, 55 core-h/run; at 128 ranks 2.90 s/step,
82 core-h/run (67% efficiency) → 128^3 runs use 64 ranks; 256^3 → 128 ranks (measured separately).
Usage: matrix.py [--block V2|S|Sint|T|all] [--print]"""
import argparse, json
DT = {32: 0.1, 64: 0.05, 128: 0.025, 256: 0.0125}
RANKS = {32: 8, 64: 16, 128: 64, 256: 128}   # 128^3: 64 ranks is the efficient choice (see timing)
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
        for r in (1, 16, 32, 64, 128): runs.append(spec("V1p", 1600, 64, tag="r%d" % r, ranks=r))   # r16 added 2026-09-17 (exact-L cases; S run has L=6.28319)
    if block in ("Lchk", "all"): runs.append(spec("Lchk", 1600, 64, tag="exactL"))   # 2026-09-17: quantify the L=6.28319 (S/V2/Sint/F) vs exact 2*pi effect (compare with S1600_h64 B0)
    if block in ("F", "all"):
        runs.append(spec("F", 1600, 256))                                                          # F1: 256^3 B0, t=20 (tag "")
        runs.append(spec("F", 1600, 256, ddtCorr=0, tag="F2"))                                     # F2 (chosen 2026-09-17 from interim report 2: least dissipative
        # Re=1600 128^3 setting = ddtCorr off, |∫Σe|/∫ε_ν = 0.10 vs B0 0.38; one-factor contrast to F1), run to t=12 (--endTime 12)
    return runs
if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--block", default="all"); ap.add_argument("--print", action="store_true"); a = ap.parse_args()
    runs = build(a.block)
    ch = {N: 0.0 for N in RANKS}
    for r in runs: ch[r["N"]] += 1
    print(json.dumps(runs, indent=0) if a.print else "")
    print("runs:", len(runs), "by N:", {k: int(v) for k, v in ch.items() if v})
