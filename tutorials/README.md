# Tutorials

Two small Taylor–Green cases that exercise the function object and print its closure residual.

| Case | Description | Runtime |
|---|---|---|
| `tgv2D` | Two-dimensional Taylor–Green vortex at Re = 100 on a 32² mesh; the analytic solution `u = e^{-2 nu t}(sin x cos y, -cos x sin y)` is known, so `tools/evalError2D.py` can be used to measure the velocity error | seconds |
| `tgv3D` | Three-dimensional Taylor–Green vortex at Re = 1600 on a 16³ mesh | seconds |

Run either case with the OpenFOAM environment loaded:

```bash
cd tgv2D
./Allrun
```

The audit writes `postProcessing/kineticEnergyAudit/0/kineticEnergyAudit.dat`; the closure residual `Rclosure` in the last column
group should stay at the level of 1e-14. Set `writeFields true` in `system/controlDict` to obtain the per-cell fields at the write
times.
