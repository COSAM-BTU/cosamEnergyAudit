#!/usr/bin/env python3
"""Statistics of decomposed binary volScalarFields (for fixed ParaView colour ranges). Usage: fieldStats.py <case> <time> <field> [<field>...]"""
import sys, glob, re, struct, numpy as np, os
case, time = sys.argv[1], sys.argv[2]
def read_scalar(path):
    b=open(path,"rb").read(); i=b.find(b"internalField")
    if i<0: return None
    j=b.find(b"nonuniform List<scalar>", i)
    if j<0:   # uniform
        m=re.search(rb"internalField\s+uniform\s+([0-9.eE+-]+)", b); return np.array([float(m.group(1))])
    k=b.find(b"(", j); nstr=b[j+len(b"nonuniform List<scalar>"):k].strip(); n=int(nstr)
    return np.frombuffer(b[k+1:k+1+8*n], dtype="<f8")
V=None
for f in sys.argv[3:]:
    vals=[]; 
    for p in sorted(glob.glob(case+"/processor*")):
        fp=os.path.join(p,time,f)
        if os.path.exists(fp): a=read_scalar(fp); vals.append(a)
    if not vals: fp=os.path.join(case,time,f); vals=[read_scalar(fp)]
    x=np.concatenate(vals); ax=np.abs(x)
    q=np.percentile(x,[0.5,1,5,25,50,75,95,99,99.5])
    print("%-8s t=%s n=%d mean=%.3e std=%.3e min=%.3e max=%.3e | p0.5=%.3e p1=%.3e p5=%.3e p50=%.3e p95=%.3e p99=%.3e p99.5=%.3e | frac>0=%.3f | |x| p99=%.3e" % (f,time,len(x),x.mean(),x.std(),x.min(),x.max(),q[0],q[1],q[2],q[4],q[6],q[7],q[8],(x>0).mean(),np.percentile(ax,99)))
