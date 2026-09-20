# cosamEnergyAudit — runtime discrete kinetic-energy audit for OpenFOAM (v2406)

`kineticEnergyAudit` is a function object that reconstructs, from the registered final fields (U, p, phi and their
old-time levels) and the case's own `fvSchemes` operators, the discrete momentum equation actually assembled by
`pimpleFoam` (laminar, PISO/PIMPLE), dots it with U^{n+1} and closes the discrete kinetic-energy budget **exactly**:

    dK/dt = -eps_nu - e_time - e_conv - e_mesh - e_cont - e_pres - e_dev + e_iter

| term | meaning | closed form / definition |
|---|---|---|
| eps_nu | face-based viscous dissipation consistent with the implicit Laplacian (incl. non-orthogonal correction) | ½Σ_f [g_f − (a_O−a_N)·ν|S_f|k_f·(∇a)_f], g_f = ν|S_f|δ_f|a_N−a_O|² |
| e_time | time-scheme error; backward split into dissipative (≥0) and storage parts | Euler ½|a−b|²V/Δt; BDF2 [X(a,b)−X(b,d)]V/Δt + |a−2b+d|²V/(4Δt); CN(ψ) ψ[ΔK/Δt − a·ddt0 V] + ((1+ψ)/2Δt)|a−b|²V |
| e_conv | convection-scheme error relative to `linear` | ½Σ_f (d_f^{sch} − m_f), d_f = s φ_f (a_P−a_N)·(U_f^{sch}−U_f^{mid}) |
| e_mesh | `linear` vs `midPoint` (mesh non-uniformity) | ½Σ_f m_f (zero on uniform hex) |
| e_cont | continuity error | ½|a_P|² (div φ)_P V |
| e_pres | pressure–velocity coupling (Rhie–Chow + ddtCorr) | p_P (G^T a)_P V, (G^T a)_P V = Σ_f ω_P S_f·(a_O − a_N) |
| e_dev | explicit `dev2(T(grad U))` term of `linearViscousStress` | −a·D2 V |
| e_iter | algebraic/iterative residual of the momentum system as solved (φ^n) | a·r V |
| e_phi | flux-lag term between φ^n and φ^{n+1} (reported, not in the closure) | a·[C(φ^{n+1}) − C(φ^n)]a V |

Primary convention: multiplier U^{n+1}, midPoint reference, φ^n (as solved in PISO), half–half face-to-cell
distribution. Alternatives (U^{n+1/2} multiplier, owner distribution) are written as global time series in the same run.
Full derivation: `DERIVATION.md` (study repository); numerically verified to round-off with `stage1/toyClosure.py`.

## Usage
```
functions
{
    kineticEnergyAudit
    {
        type            kineticEnergyAudit;
        libs            (cosamEnergyAuditFunctionObjects);
        writeFields     true;      // per-cell eNum, nuNum (masked), eTime, eConv, eMesh, eCont, ePres, eDev, eIter, epsNu, SS
        maskCoeff       1e-3;      // nu_num mask: cells with S:S <= maskCoeff*<S:S> are masked
        // instrumentedResidual rTrueV;   // optional: compare with a solver-written as-solved residual field (V4)
        // offline true;                  // recompute from written raw inputs (U_prev, phi_prev, U_prevprev, ddt0(U))
        writeControl    writeTime;
    }
}
```
Time series: `postProcessing/kineticEnergyAudit/<time>/kineticEnergyAudit.dat` (one row per step; columns documented in
the header; includes `Rclosure`, `cumClosure`, `RcellMax`, consistency checks `chkConv/chkDiff/chkTime/chkPres`, the
telescoping-flux sums, the alternative-convention totals and the mask fraction).

Offline recomputation of any written time: `pimpleFoam -postProcess -dict system/auditOffline -time <t1>,<t2>`
(bit-identical to the runtime rows).

## Requirements and restrictions (v0.1)
- OpenFOAM v2406, `pimpleFoam` with `simulationType laminar`; all non-empty patches must be coupled (cyclic/processor).
- `ddtSchemes`: Euler, backward (fixed or variable Δt — variable Δt exact since v1.1), CrankNicolson ψ. First-step rules of the schemes are mirrored.
- Convection schemes: any `Gauss <scheme>` (the scheme's own face interpolation is used; e_conv is measured against
  `linear`). `laplacianSchemes Gauss linear corrected`, `gradSchemes Gauss linear`.
- Use `relaxationFactors { equations { ".*" 1; } }`, `pRefCell/pRefValue`, `momentumPredictor yes`.
- **Coupled patches must match exactly.** On perturbed/periodic meshes the closure holds to round-off only if the cyclic face pairs are
  geometrically identical (weights w and 1−w complementary). A box length rounded to 6 digits (e.g. `foamDictionary` rewriting
  `blockMeshDict`, L = 6.28319 ≠ 2π) plus a periodic point perturbation left mismatched cyclics and R_closure ≈ 1e-9 (2026-09-17);
  the pressure-transpose check `chkPres` flags this.
- **Parallel runs: initialise the fields in parallel** (`decomposePar` → `mpirun setExprFields -parallel`). Fields set
  in serial and decomposed carry unevaluated processor-patch values into the first momentum assembly; the audit flags
  this as R_closure ≈ 6e-3 at step 1 (a reproducibility note in the paper).

## Build
```
source /usr/lib/openfoam/openfoam2406/etc/bashrc
cd cosamEnergyAudit && wmake      # -> $FOAM_USER_LIBBIN/libcosamEnergyAuditFunctionObjects.so
```

## Changes
- **v1.1 (2026-09-18):** `backward` with a *variable* time step: v1.0 evaluated e_time with the constant-step BDF2 closed form
  (X(a,b) − X(b,d) + |a−2b+d|²/4Δt), which is not an identity when Δt ≠ Δt₀ (adaptive Δt run: R_closure up to 3e-2 at every
  step where Δt changed; caught by `chkTime`). v1.1 uses the exact operator route e_time = a·T(a) V − ΔK/Δt with OpenFOAM's
  variable-step coefficients whenever Δt ≠ Δt₀; the storage part keeps the constant-step G-norm change (telescopes exactly),
  the dissipative part is the remainder (sign-definite only for Δt = Δt₀). Constant-step results are bit-identical to v1.0.
  Δt₀ is taken from `Time::deltaT0` (also available offline from `<time>/uniform/time`).

## Verification status (2026-09-15, re-run 2026-09-18 with v1.1: identical)
2D TGV 32² (Euler/backward/backward-0/CN0.9; linear/upwind/LUST/linearUpwind/limitedLinear/cubic; uniform and
perturbed mesh): R_closure ≤ 5e-14 incl. step 1; 3D 16³ uniform/perturbed (non-orthogonality up to 23°): ≤ 3e-15;
polyhedral (Gmsh tet → polyDualMesh) cyclic mesh: 5e-15; 4 ranks vs serial 1.6e-12; offline and restart bit-identical;
instrumented-solver residual (A U − H(U) + ∇p) agrees with e_iter for linear/upwind, and isolates the explicit
limiter/gradient lag for LUST and limitedLinear.

License: GPL-3.0-or-later. COSAM, Bursa Technical University.
