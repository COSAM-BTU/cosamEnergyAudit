#!/usr/bin/env python3
"""Single-node pseudo-spectral (Fourier–Galerkin) incompressible Navier–Stokes solver for the 3D Taylor–Green vortex
on [0,2pi]^3, used to generate the reference solutions. Rotational form u x omega, 2/3-rule dealiasing,
explicit RK4 (viscous term explicit), real-to-complex FFTs (scipy.fft with threads; pyfftw if available).
Output columns match the UCL reference file: t  Ek  -dEk/dt  enstrophy  (+ 2*nu*enstrophy, |u|max).
Usage: tgvSpectral.py --N 128 --Re 1600 --dt 0.005 --T 20 --out spec_Re1600_N128.dat [--every 1] [--threads 32]"""
import numpy as np, argparse, time, os, sys
ap = argparse.ArgumentParser()
ap.add_argument("--N", type=int, default=64); ap.add_argument("--Re", type=float, default=1600.0)
ap.add_argument("--dt", type=float, default=0.005); ap.add_argument("--T", type=float, default=20.0)
ap.add_argument("--out", default=None); ap.add_argument("--every", type=int, default=1); ap.add_argument("--threads", type=int, default=os.cpu_count())
ap.add_argument("--dealias", default="2/3")
a = ap.parse_args()
N, nu, dt = a.N, 1.0/a.Re, a.dt
backend = "numpy"
try:
    import pyfftw
    pyfftw.interfaces.cache.enable(); pyfftw.interfaces.cache.set_keepalive_time(3600)
    from pyfftw.interfaces import scipy_fft as sf
    backend = "pyfftw"
except Exception:
    try:
        import scipy.fft as sf; backend = "scipy"
    except Exception:
        sf = None
if sf is not None:
    rfftn = lambda x: sf.rfftn(x, workers=a.threads); irfftn = lambda x: sf.irfftn(x, s=(N, N, N), workers=a.threads)
else:
    rfftn = np.fft.rfftn; irfftn = lambda x: np.fft.irfftn(x, s=(N, N, N))
print("FFT backend:", backend, "threads:", a.threads, file=sys.stderr)
x = 2*np.pi*np.arange(N)/N
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
u0 = np.stack([np.sin(X)*np.cos(Y)*np.cos(Z), -np.cos(X)*np.sin(Y)*np.cos(Z), np.zeros_like(X)])
k1 = np.fft.fftfreq(N, 1.0/N); k3 = np.fft.rfftfreq(N, 1.0/N)
KX, KY, KZ = np.meshgrid(k1, k1, k3, indexing="ij"); K = np.stack([KX, KY, KZ]); K2 = KX**2 + KY**2 + KZ**2
K2inv = np.where(K2 > 0, 1.0/np.where(K2 > 0, K2, 1), 0.0)
kmax = N/3.0 if a.dealias == "2/3" else N/2.0
mask = (np.abs(KX) < kmax) & (np.abs(KY) < kmax) & (np.abs(KZ) < kmax)
Ninv = 1.0/N**3
# Parseval weights for r2c layout (kz=0 and kz=N/2 planes counted once)
w = np.full(K2.shape, 2.0); w[..., 0] = 1.0
if N % 2 == 0: w[..., -1] = 1.0
def fft3(u): return np.stack([rfftn(u[i])*Ninv for i in range(3)])
def ifft3(uh): return np.stack([irfftn(uh[i]/Ninv) for i in range(3)])
def curl_hat(uh):
    return 1j*np.stack([K[1]*uh[2] - K[2]*uh[1], K[2]*uh[0] - K[0]*uh[2], K[0]*uh[1] - K[1]*uh[0]])
def rhs(uh):
    u = ifft3(uh); om = ifft3(curl_hat(uh))
    cr = np.stack([u[1]*om[2] - u[2]*om[1], u[2]*om[0] - u[0]*om[2], u[0]*om[1] - u[1]*om[0]])
    ch = fft3(cr)*mask
    kdotc = (K[0]*ch[0] + K[1]*ch[1] + K[2]*ch[2])*K2inv
    ch -= K*kdotc                       # Leray projection (pressure)
    return ch - nu*K2*uh
def energy(uh):   # 0.5<|u|^2>
    return 0.5*np.sum(w*(np.abs(uh[0])**2 + np.abs(uh[1])**2 + np.abs(uh[2])**2))
def enstrophy(uh):
    oh = curl_hat(uh); return 0.5*np.sum(w*(np.abs(oh[0])**2 + np.abs(oh[1])**2 + np.abs(oh[2])**2))
uh = fft3(u0)*mask
kdot = (K[0]*uh[0] + K[1]*uh[1] + K[2]*uh[2])*K2inv; uh -= K*kdot
nsteps = int(round(a.T/dt)); t = 0.0
out = open(a.out, "w") if a.out else sys.stdout
out.write("# tgvSpectral.py N=%d Re=%g dt=%g dealias=%s RK4 rotational-form; columns: t Ek -dEk/dt enstrophy 2nu*enstrophy\n" % (N, a.Re, dt, a.dealias))
Eprev = energy(uh); Om = enstrophy(uh)
out.write("%.6f %.12g %.12g %.12g %.12g\n" % (0.0, Eprev, 2*nu*Om, Om, 2*nu*Om)); out.flush()
t0 = time.time()
for n in range(1, nsteps + 1):
    k1_ = rhs(uh); k2_ = rhs(uh + 0.5*dt*k1_); k3_ = rhs(uh + 0.5*dt*k2_); k4_ = rhs(uh + dt*k3_)
    uh = uh + dt*(k1_ + 2*k2_ + 2*k3_ + k4_)/6.0
    t = n*dt
    if n % a.every == 0 or n == nsteps:
        E = energy(uh); Om = enstrophy(uh)
        out.write("%.6f %.12g %.12g %.12g %.12g\n" % (t, E, (Eprev - E)/(dt*a.every), Om, 2*nu*Om)); out.flush()
        Eprev = E
    if n % max(1, nsteps//20) == 0:
        print("step %d/%d t=%.3f Ek=%.10f Omega=%.8f  %.1f s/step" % (n, nsteps, t, energy(uh), Om, (time.time()-t0)/n), file=sys.stderr, flush=True)
