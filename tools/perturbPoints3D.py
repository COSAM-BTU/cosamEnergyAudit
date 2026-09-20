#!/usr/bin/env python3
"""Periodic perturbation of an ASCII polyMesh/points file for the [0,2pi]^3 TGV box.
Recipe (a) sinusoidal (INSTRUCTIONS §8.4): x'_i = x_i + A sin(x) sin(y) sin(z)  -> vanishes on all box planes,
  so cyclic faces are unchanged; smooth, keeps 2nd order.
Recipe (b) random: interior points displaced by U(-alpha, alpha)*dx with a fixed seed; boundary-plane points fixed.
Usage: perturbPoints3D.py sin A | perturbPoints3D.py rand alpha dx [seed]"""
import re, sys, math, random
mode = sys.argv[1]
L = 2*math.pi; tol = 1e-9
p = 'constant/polyMesh/points'
s = open(p).read()
head, rest = s.split('(', 1)
body, tail = rest.rsplit(')', 1)
pts = re.findall(r'\(([^()]*)\)', body)
out = []
if mode == 'sin':
    A = float(sys.argv[2])
    for t in pts:
        x, y, z = map(float, t.split())
        d = A*math.sin(x)*math.sin(y)*math.sin(z)
        out.append('(%.15g %.15g %.15g)' % (x + d, y + d, z + d))
elif mode == 'rand':
    alpha = float(sys.argv[2]); dx = float(sys.argv[3]); seed = int(sys.argv[4]) if len(sys.argv) > 4 else 12345
    rng = random.Random(seed)
    for t in pts:
        x, y, z = map(float, t.split())
        onB = any(abs(c) < tol or abs(c - L) < tol for c in (x, y, z))
        if onB:
            out.append('(%.15g %.15g %.15g)' % (x, y, z))
        else:
            out.append('(%.15g %.15g %.15g)' % (x + rng.uniform(-alpha, alpha)*dx, y + rng.uniform(-alpha, alpha)*dx, z + rng.uniform(-alpha, alpha)*dx))
else:
    raise SystemExit('mode must be sin or rand')
open(p, 'w').write(head + '(\n' + '\n'.join(out) + '\n)' + tail)
print('perturbed', len(out), 'points, mode', mode)
