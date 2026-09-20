# Derivation of the discrete kinetic-energy identity

This document derives, term by term, the identity that the `kineticEnergyAudit` function object evaluates, and it fixes the
attribution convention. The reference implementation is OpenFOAM v2406; the operator definitions below are those of its source.
All identities were verified numerically to round-off on two-dimensional periodic meshes (uniform and randomly perturbed) with
random fields, with relative errors below 2e-14.

## 0. Notation

- Cell `P` (volume `V_P`, centre `x_P`); face `f` (area vector `S_f` pointing from the owner `O(f)` to the neighbour `N(f)`,
  magnitude `|S_f|`, centre `x_f`); `s_Pf = +1` if `P = O(f)` and `-1` otherwise. Cyclic and processor faces are treated as
  interior faces, the neighbour value being obtained from the coupled patch.
- Face flux `phi_f` (volumetric, positive out of the owner): the Rhie-Chow flux produced by the pressure equation.
  `(div phi)_P V_P = sum_f s_Pf phi_f`.
- Time levels: `a = U^{n+1}`, `b = U^n`, `d = U^{n-1}`; cell energy `K_P = 1/2 |a_P|^2 V_P`; domain energy `K = sum_P K_P / V`.
- Linear interpolation weight `w_f` (owner side): `U_f^lin = w_f U_O + (1-w_f) U_N`, with
  `w_f = |S_f.(x_N - x_f)| / (|S_f.(x_f - x_O)| + |S_f.(x_N - x_f)|)`; mid-point interpolation `U_f^mid = 1/2 (U_O + U_N)`.
  On a uniform hexahedral mesh `w_f = 1/2`.
- Selected convection scheme: `U_f^sch = U_f^mid + delta_f`, where `delta_f` is the scheme-specific deviation (for limited or
  gradient-based schemes it depends on `U`).

## 1. The discrete momentum equation assembled by the solver

For `pimpleFoam` (laminar, PISO mode, one outer corrector) the momentum predictor is
`fvm::ddt(U) + fvm::div(phi, U) + divDevReff(U) = 0`, solved as `solve(UEqn == -fvc::grad(p))`, with
`divDevReff = -fvc::div(nu dev2(T(grad U))) - fvm::laplacian(nu, U)`. In cell form (multiplied by `V_P`):

    T(a)_P V_P + sum_f s phi^n_f U_f^sch - sum_f s J_f(a) - D2_P V_P + sum_f s S_f p_f = r_P V_P        (1)

- The convective flux is `phi^n`: the momentum matrix is assembled before the pressure correctors update the flux. With more than
  one outer corrector the flux of the last outer iteration approaches `phi^{n+1}`.
- Diffusive face flux `J_f = nu_f |S_f| [delta_f (U_N - U_O) + k_f.(grad U)_f]` (`Gauss linear corrected`), where `delta_f` is the
  non-orthogonal delta coefficient, `k_f` the non-orthogonal correction vector (zero on an orthogonal mesh) and `(grad U)_f` the
  linearly interpolated Gauss gradient.
- `D2_P = [fvc::div(nu dev2(T(grad U)))]_P` is explicit; it vanishes in the continuous limit for a divergence-free field but not
  discretely.
- `r_P` is the momentum residual: what remains when (1) is reassembled with the final `a`, `p^{n+1}` and `phi^n`. It contains the
  PISO splitting error, the linear-solver residuals, and the lag of every nonlinear ingredient (limiters, explicit gradients,
  non-orthogonal corrections, `dev2`) that the solver evaluates with a provisional velocity and the audit with `a`.

Primary convention: (1) is dotted with `a_P` and summed over the cells.

## 2. Time term: `e_time := a.T(a) V - (K^{n+1} - K^n)/dt`

### 2.1 Euler, `T(a) = (a - b)/dt`

`a.(a-b) = 1/2 |a|^2 - 1/2 |b|^2 + 1/2 |a-b|^2`, hence `e_time = 1/2 |a-b|^2 V/dt >= 0`: exactly dissipative.

### 2.2 Backward (BDF2), constant step, `T(a) = (3a - 4b + d)/(2 dt)`

The G-stability identity is
`2 a.(3a - 4b + d) = |a|^2 + |2a - b|^2 - |b|^2 - |2b - d|^2 + |a - 2b + d|^2`.
With `E_G(a,b) := 1/4 (|a|^2 + |2a-b|^2)` and `X(a,b) := E_G(a,b) - 1/2 |a|^2 = 1/4 (a-b).(3a-b)`,

    a.T(a) = [E_G(a,b) - E_G(b,d)]/dt + |a - 2b + d|^2/(4 dt)
    e_time = [X(a,b) - X(b,d)] V/dt  +  |a - 2b + d|^2 V/(4 dt)  =: e_time^store + e_time^diss

`e_time^diss >= 0` is strictly dissipative; `e_time^store` telescopes in time (its integral reduces to boundary terms) and has no
definite sign. Check: for `a - b = b - d = v` the term reduces to `1/2 |v|^2 V/dt`, storage only.

For a variable step OpenFOAM uses `T(a) = (c a - c0 b + c00 d)/dt` with `r = dt/dt0`, `c = 1 + r/(1+r)`, `c00 = r^2/(1+r)`,
`c0 = c + c00`, which reduces to the constant-step coefficients `3/2, 2, 1/2` for `r = 1`. The time term is then still evaluated
exactly from its definition; see the addendum below. In the first step of a run the previous level is undefined and the scheme
reduces to Euler.

### 2.3 Crank-Nicolson with off-centering `psi`, `T(a) = ((1+psi)/dt)(a - b) - psi ddt0^n`

`ddt0^n` is the stored time derivative of the previous step, which OpenFOAM updates at the first `ddt` call of a step. With the
Euler identity,

    e_time = psi [(K^{n+1} - K^n)/dt - a.ddt0^n V] + ((1+psi)/(2 dt)) |a-b|^2 V

For `psi = 0` this is the Euler term; the first step of a run uses Euler.

### 2.4 Alternative multiplier `U^{n+1/2} = 1/2 (a + b)`

`1/2 (a+b).(a-b)/dt = (K^{n+1} - K^n)/dt` exactly, so the Euler time term vanishes identically with this multiplier and the
temporal error migrates into the other terms. It is reported as a global time series only; the per-cell fields use `a`.

## 3. Convective term: `W_P := a_P . sum_f s phi_f U_f^sch`

Writing `U_f^sch = U_f^mid + delta_f` and `a_P.delta_f = 1/2 (a_P - a_N').delta_f + U_f^mid.delta_f`, the contribution of one face
to the owner cell (for the neighbour, `s` changes sign) is

    s phi_f a_P.U_f^sch = s phi_f 1/2 |a_P|^2  +  1/2 d_f  +  s F^C_f
    d_f   := s_Pf phi_f (a_P - a_P').delta_f          (the same value in both cells)
    F^C_f := phi_f [1/2 a_O.a_N + U_f^mid.delta_f]    (the same in both cells, hence telescoping)

where `P'` is the cell on the other side of `f`. This defines, per cell,

    e_cont(P) := 1/2 |a_P|^2 (div phi)_P V_P
    m_f       := d_f for linear interpolation = s phi_f (w_f - 1/2)(a_P - a_P').(a_O - a_N)   (zero on a uniform hex mesh)
    e_mesh(P) := 1/2 sum_f m_f
    e_conv(P) := 1/2 sum_f (d_f^sch - m_f)            (deviation of the scheme from linear; identically zero for linear)

so that `W_P = e_cont + e_conv + e_mesh + sum_f s F^C_f` and, globally, `sum_P W_P = sum_P (e_cont + e_conv + e_mesh)`.

Particular schemes: for `upwind`, `delta_f = +-1/2 (a_O - a_N)` depending on the flow direction, so that
`d_f = 1/2 |phi_f| |a_O - a_N|^2 >= 0` face by face, which provides a known-answer test. For `linear`,
`delta_f = (w_f - 1/2)(a_O - a_N)`. For `cubic`, `LUST`, `linearUpwind` and `limitedLinear`, `delta_f` involves gradients or
limiters and `d_f` has no definite sign.

Face-to-cell distribution: half to each cell (primary); the owner-side alternative changes the local field but not the global sum.
Flux convention: `phi^n` (primary); the effect of `phi^{n+1}` is reported as the flux-lag term `e_phi`.

## 4. Pressure term: `a_P.(G p)_P V_P = sum_f s S_f.a_P p_f`, `p_f = w_f p_O + (1-w_f) p_N`

Summation by parts gives `sum_P a_P.(Gp)_P V_P = sum_P p_P (G^T a)_P V_P` with

    (G^T a)_P V_P := sum_{f in P} omega_Pf S_f.(a_O(f) - a_N(f)),    omega_Pf = w_f if P = O(f), else 1 - w_f

The same expression `omega_Pf S_f.(a_O - a_N)` is used on both sides of a face; writing it with the opposite sign on the neighbour
side breaks the identity. Equivalently `(G^T a)_P V_P = -sum_f s S_f.U~_f` with the weight-transposed interpolation
`U~_f = (1 - omega_Pf) a_P + omega_Pf a_P'`, which on a uniform hexahedral mesh is minus the mid-point divergence of the
cell-centred velocity. The cell identity is exact:

    a_P.(Gp)_P V_P = p_P (G^T a)_P V_P + sum_f s X_f,    X_f := w_f p_O (S_f.a_N) + (1 - w_f) p_N (S_f.a_O)

with `X_f` the same in both cells, hence telescoping. This defines `e_pres(P) := p_P (G^T a)_P V_P`.

Interpretation: the pressure equation drives the divergence of the face flux to zero, but not the weight-transposed divergence of
the cell-centred velocity. The difference is the Rhie-Chow filter together with the transient flux correction, since
`phi = S.interp(HbyA) + ddtCorr - rAU_f grad_f p . S` while `a = HbyA - rAU grad p`. The share of the transient correction cannot
be separated within one run; it is obtained by comparing paired runs with the correction on and off (`ddtCorr` set to zero in the
`ddtSchemes` entry). An alternative local form, `e_pres'(P) := sum_f s p_f S_f.(a_P - U_f^mid)`, has the same global sum and a
different local distribution.

## 5. Diffusive term: `a_P . sum_f s J_f`

Splitting `a_P.J_f` into a symmetric and an antisymmetric part gives, for the owner,
`a_O.J = 1/2 (a_O + a_N).J + 1/2 (a_O - a_N).J`, and for the neighbour
`-a_N.J = -1/2 (a_O + a_N).J + 1/2 (a_O - a_N).J`. The symmetric part `F^D_f := 1/2 (a_O + a_N).J_f` telescopes, and because the
diffusive term enters (1) with a minus sign every cell receives `-1/2 (a_O - a_N).J_f` from it. With the orthogonal part of `J_f`
equal to `nu |S_f| delta_f (a_N - a_O)`,

    eps_nu(P) := 1/2 sum_f [ g_f - (a_O - a_N) . nu_f |S_f| k_f.(grad a)_f ],   g_f := nu_f |S_f| delta_f |a_N - a_O|^2 >= 0

so that `a_P . sum_f s J_f = -eps_nu(P) + sum_f s F^D_f`. This face-based dissipation is consistent with the implicit Laplacian
that the solver inverts; on a polyhedral mesh its local value is not guaranteed to be non-negative. Two further quantities are
defined: the diagnostic `e_diff(P) := eps_nu(P) - 2 nu (S:S)_P V_P`, which does not enter the closure, and the explicit stress
term `e_dev(P) := -a_P.D2_P V_P`, which does.

## 6. Algebraic residual: `e_iter(P) := a_P.r_P V_P`

The residual is that of the equation as solved, with `phi^n`. It contains the PISO splitting error, the linear-solver residuals,
and the lag of the nonlinear ingredients. For the linear scheme on an orthogonal hexahedral mesh the last contribution reduces to
the `dev2` term, so that `e_iter` is essentially the splitting residual; with LUST, limited schemes or polyhedral meshes it grows,
which is what the comparison with an instrumented solver measures. The term is called an algebraic residual, not a dissipation.

## 7. Closure

Dotting (1) with `a_P V_P` and collecting Sections 2 to 6 gives, per cell,

    (K_P^{n+1} - K_P^n)/dt = -e_time - e_conv - e_mesh - e_cont - eps_nu - e_dev - e_pres + e_iter
                             - sum_f s (F^C_f - F^D_f + X_f)

and globally, since the face fluxes telescope,

    dK/dt = -eps_nu - e_time - e_conv - e_mesh - e_cont - e_pres - e_dev + e_iter                     (2)

The closure residual is `R_closure := |LHS - RHS| / max(|LHS|, |RHS|, eps0)` with a floor `eps0 = 1e-15`, and the cumulative
closure is `|K(t) - K(0) + int_0^t (eps_nu + sum_i e_i - e_iter) dt'| / K(0)`. Equation (2) is an algebraic identity; its residual
tests the implementation (face telescoping, coupled-patch handling, floating-point reduction), not the solver. Whether the audit
reconstructs the equation that the solver actually solved is a separate question, answered with an instrumented copy of the solver
that writes the residual of its own momentum system.

## 8. Convention pairs

Each alternative is written as a global time series in the same run; the per-cell fields use the primary convention.

| Choice | Primary | Alternative |
|---|---|---|
| Multiplier | `U^{n+1}` | `U^{n+1/2}` (Section 2.4) |
| Convective reference | mid-point (`e_conv` = scheme - linear; `e_mesh` = linear - mid-point) | linear (`e_mesh` identically zero) |
| Flux | `phi^n` | `phi^{n+1}` (reported as `e_phi`) |
| Face-to-cell distribution | half-half | owner side |
| Local form of `e_pres` | `p_P (G^T a)_P V_P` | `sum_f s p_f S_f.(a_P - U_f^mid)` |

## 9. Summary of definitions and signs

| Term | Definition (per cell, times `V_P`) | Sign |
|---|---|---|
| `e_time` | `a.T(a) V - dK/dt` (closed forms in Section 2) | Euler >= 0; BDF2 dissipative part >= 0, storage part either |
| `e_conv` | `1/2 sum_f (d_f^sch - m_f)` | either (upwind: non-negative in total) |
| `e_mesh` | `1/2 sum_f m_f` | either (zero on a uniform hexahedral mesh) |
| `e_cont` | `1/2 \|a\|^2 (div phi) V` | either (at the level of the pressure tolerance) |
| `e_pres` | `p (G^T a) V` | either |
| `eps_nu` | `1/2 sum_f [g_f - (a_O - a_N).nu \|S_f\| k_f.(grad a)_f]` | non-negative on an orthogonal mesh; either on a polyhedral one |
| `e_dev` | `-a.D2 V` | either, small |
| `e_iter` | `a.r V` | either |
| `e_phi` | flux-lag term | either (reported, not part of (2)) |
| `e_diff` | `eps_nu - 2 nu (S:S) V` | diagnostic |

## Addendum: BDF2 with a variable time step

With `r = dt/dt0` the `backward` scheme of OpenFOAM uses `T(a) = (c a - c0 b + c00 d)/dt`, `c = 1 + r/(1+r)`, `c00 = r^2/(1+r)`,
`c0 = c + c00`. The time term is evaluated from its definition, `e_time = a.T(a) V - dK/dt` (the operator route, checked against
`fvc::ddt(U)` at every step). The constant-step split into a storage and a non-negative dissipative part is an identity for
`r = 1` only; for `r != 1` the function object keeps the change of the G-norm as the storage part and reports the remainder as the
dissipative part, which is then not sign-definite, because a sign-definite split for variable steps requires a step-dependent
G-norm that does not telescope. In a 16^3 test with the time step changing at every step by up to 17 per cent the closure residual
stays below 3e-15, whereas the constant-step closed form gives residuals of up to 5e-2 at the steps at which the step changes.
