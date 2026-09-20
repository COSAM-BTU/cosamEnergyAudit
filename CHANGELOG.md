# Changelog

## 1.1.1 (2026-09-20)

- Documentation: the derivation document is rewritten in English and is self-contained.
- Tools: campaign directory, output directory and reference solution are command-line arguments; the polyhedral mesh directory is
  read from `POLYMESH_DIR`; the figure and table scripts take `--data` and `--out`.
- Tutorials: the run scripts use the standard OpenFOAM `RunFunctions` instead of a fixed installation path.
- No change to the function object source; results are unaffected.

## 1.1.0 (2026-09-18)

- Variable time step: the BDF2 time term is evaluated with the operator route `a . T(a) V - dK/dt` whenever the step changes,
  instead of the constant-step closed form. Bit-identical to 1.0.0 for a constant time step. Found by the closure check itself
  during an adaptive-time-step run, which reported residuals of up to 3e-2 at exactly the steps at which the step changed.
- The first old-time step `deltaT0` is read from `Time::deltaT0Value()`.
- Documentation: derivation addendum for the variable-step scheme; note on matching coupled patches.

## 1.0.0 (2026-09-15)

- First release: exact discrete kinetic-energy identity for the collocated finite-volume discretization of OpenFOAM, evaluated at
  run time with the unmodified `pimpleFoam` executable; global time series, per-cell fields, alternative conventions, diagnostics
  and offline recomputation.
