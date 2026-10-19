# Plan: polar vorticity-streamfunction solver, validated on the cylinder at Re = 40

> **Status (2026-10-18): implemented** and merged (PR #57): `src/polar.c`,
> `src/cylinder.c`, `tests/test_polar.c`, `tools/cylinder_validate.py`; documented
> in the methodology manual (`docs/methodology.tex`, section "Polar Grid: Flow
> Past a Circular Cylinder"). How it differs from the plan below:
>
> - **The polar grid is opt-in**, with `--geometry polar`; Cartesian stays the default. The only new
>   `solver_config` field is `const polar_grid *polar`, NULL for Cartesian, and every structure that
>   is filled field by field is now zero-initialised. All 196 existing checks, the CLI tests and the
>   stored-result regression (bit-identical) pass unchanged.
> - **Far-field omega.** The outer omega is zero-gradient on outflow, the same as `freestream`; the
>   analytic Schlichting omega is not imposed. The far-field model enters through psi_out only.
> - **Phase 6 items.**
>   - The compact fourth-order xi Poisson solver (`--poisson-order 4`) is implemented.
>   - Local pseudo-time is not implemented; the runs were fast enough without it, 44 s for
>     128 x 90 to steady state.
>   - The paper's Test 1 with the reference field is not possible: the field isn't public.
> - **Manufactured solution.** Orders 4 and 6 converge at about 2.7, not 3. The tests assert
>   at least 2.5.

Reference: R. Gautier, D. Biau, E. Lamballais, *A reference solution of the flow
over a circular cylinder at Re = 40*, arXiv:1310.6641 (`docs/1310.6641v1.pdf` on
`main`).

## 0. What the paper gives us to validate against

The paper's full reference field (`refsolcyl.tgz` with `interpol.m`) is **not** in
the PDF. It is "available on the journal website, also on demand to the
corresponding author". Without it we cannot reproduce their Test 1 (Dirichlet
outer boundary interpolated from the reference field at r = 10.5D, with an L2
convergence study). What we can compare against:

| Quantity | Reference (fine grid 200x1024, r_inf = 40) | Source |
|---|---|---|
| C_D | 1.4931 | Table 2 |
| L_w / D (from the rear stagnation point) | 2.2360 | Table 2 |
| theta_s (from the front stagnation point) | 126.3945 deg | Table 2 |
| recirculation centre a/D, b/D | 0.71, 0.59 | Table 1 |
| u_x, u_y, omega at (x, y) = (3, 3) | 1.0670111347, 1.5134054e-2, -6.5846e-6 | Appendix A |
| r_inf sensitivity (r_inf = 30/40/50) | C_D 1.4906/1.4931/1.4943, L_w 2.2346/2.2360/2.2369, theta_s 126.406/126.395/126.389 | Table 2 |
| centreline velocity u_x(r, theta = 0) | rises monotonically to 1, with no inflection near the outflow | Fig. 1 |
| literature scatter | 1.48 < C_D < 1.62, 2.13 < L_w < 2.35 | Table 1 |

Conventions (D = 1, U_inf = 1, cylinder radius R = 1/2, Re = U D / nu = 40):
theta = 0 points downstream (the wake axis), as in the paper's u_1 formula.
The paper's theta_s is measured from the **front** stagnation point:
theta_s = 180 deg - (our angle of separation).

**Action (optional, recommended):** email the corresponding author
(eric.lamballais@univ-poitiers.fr) for `refsolcyl.tgz`. With it, Phase 6 adds
the paper's full convergence test (Section 4 / Fig. 4). The rest of the plan
does not depend on it.

Store the numbers above in `tests/reference/cylinder_re40.csv` so that tests and
tools read them from one place.

## 1. Formulation

### Coordinates

Log-polar (conformal) mapping, uniform in (xi, theta):

```
r = R e^xi,   xi in [0, xi_inf],  xi_inf = ln(r_inf / R)
theta in [0, 2 pi), periodic
```

Square cells (d xi = d theta = 2 pi / N_theta) give N_r = xi_inf / d xi + 1. For
r_inf = 40: N_theta = 128, 256, 512 gives N_r = 90, 180, 358. The wall cell is
R d xi = 0.012 at N_theta = 256, so the Re = 40 boundary layer (about D/sqrt(Re) = 0.16) gets
about 13 points. The cells grow linearly with r. That is what lets r_inf reach 30 to 50
cheaply.

### Equations

With h = r the scale factor of the conformal map:

```
psi_xixi + psi_thth = -r^2 omega                                  (Poisson)
u_r  =  (1/r) psi_th,    u_th = -(1/r) psi_xi                     (velocity)
d omega/dt = -(u_r omega_xi + u_th omega_th) / r
             + (omega_xixi + omega_thth) / (Re r^2)               (transport)
```

Free stream: psi = r sin(theta). Continuity is
`(1/r^2) [D_xi (D_th psi) - D_th (D_xi psi)]`. This is zero to round-off for the
Kronecker operators, because D_xi and D_th commute, the same property the Cartesian solver
relies on.

### Mapping onto the existing data layout

The code stores fields as `ny` rows of `nx` values, with `x` contiguous. Use:

| code | polar |
|---|---|
| column j, "x", `nx` | theta_j = 2 pi j / N_theta, **periodic** (`SDiff1_periodic`, `SDiff2_periodic`) |
| row i, "y", `ny` | xi_i = i d xi; row 0 = cylinder, row ny-1 = outer boundary (`SDiff1`, `SDiff2`, one-sided end rows) |
| field `u` | u_theta |
| field `v` | u_r |

The transport right-hand side is then exactly the existing loop of `dwdt()`, with
the advective part multiplied by `1/r_i` and the viscous part by `1/r_i^2` (one
factor per row):

```
out_ij = (1/r_i) (-u dwdx - v dwdy) + (1/(Re r_i^2)) (d2wdx2 + d2wdy2)
```

**Sign trap:** (x, y) = (theta, xi) is left-handed relative to (r, theta). The
polar velocity recovery is therefore u = -(1/r) DY psi and v = +(1/r) DX psi, the
opposite of the Cartesian u = DY psi, v = -DX psi. The MMS test in Phase 2 must be in
place before anything else depends on this.

### Boundary conditions

**Cylinder (row 0):** psi = 0 and psi_xi = 0 (no slip). The wall vorticity is
omega_w = -(1/R^2) psi_xixi. That is the existing Briley formula in xi:
`omega_w = briley(psi_0..psi_3, -R u_th_wall, d xi) / R^2`, with u_th_wall = 0 for
a fixed cylinder. A non-zero value leaves room for a rotating cylinder later. The
analogue of `wall_closure 0` is omega_w = (1/R) D_xi u_th on the wall row, using the
one-sided DY row. Implement both. Briley is the default; compare them in Phase 5.

**Outer boundary (row ny-1).** The FFT-in-theta Poisson solver needs psi to have
**one boundary type for all theta**: Dirichlet, with values that can vary in theta. omega is
only ever *set* on the boundary row, so its type can change with theta.

- `freestream` (Phase 2): psi = r_inf sin(theta). omega = 0 on inflow
  (cos theta < 0). Zero normal gradient on outflow (omega_N = 2 omega_{N-1} - omega_{N-2}, or a
  plain copy). This is the paper's "Dirichlet" case. It should show the
  spurious inflection in the centreline velocity in Fig. 1.
- `farfield` (Phase 4): the paper's triple decomposition u = u_inf + u_1 + u_2,
  converted to psi and omega:
  - u_1 is the Schlichting wake deficit (paper's formula; the paper folds U_inf into
    its u_1 . e_x = 1 - ...). Define the deficit d(r, theta) =
    C_D sqrt(Re) / (4 sqrt(pi)) * exp(-Re r (1 - cos theta)/2) / sqrt(r), with
    u_1 . e_x = -d and u_1 . e_y = -(sin theta / 2) d.
  - u_2 is the mass correction: a uniform radial outflow (Hasan et al.) chosen so that the net flux through r_inf is zero, which is
    what makes psi single-valued.
  - psi_out(theta) = r_inf sin(theta) + integral_0^theta r_inf (u_1 + u_2).e_r
    d theta'. Compute it on the theta grid by FFT: remove the mean of the integrand
    (that *is* u_2), then divide mode m by i m. psi_out comes out odd in theta,
    and psi_out(0) = 0 matches psi = 0 on the symmetry line.
  - omega_out = -d(u_1.e_x)/dy written in polar form, which is non-zero only in the wake.
    Zero-gradient outflow is an alternative, if the analytic form proves fragile.
  - C_D appears in the boundary condition. Use the running surface C_D (Section 3, item 1),
    updated every K steps with under-relaxation, C_D <- C_D + 0.5 (C_D_new - C_D),
    starting from 1.5, and held fixed within a step so that RK4 stages see one
    boundary. The run stays deterministic. The paper notes the same coupling between resolution and the boundary condition.

### Poisson solver (new): real FFT in theta, banded solve in xi

Per theta mode m (`nx/2 + 1` modes, `fftw_plan_many_dft_r2c` over rows):

```
(delta_xixi - lambda_m) psi_m = -(r^2 omega)_m,   psi_m(0) = 0,  psi_m(xi_inf) = (psi_out)_m
```

- Take lambda_m = -(symbol of `SDiff2_periodic(N_theta, order)` at m). The theta
  operator then matches DX2 exactly, as `periodic_eigenvalues()` already does for the
  periodic solver.
- `poisson_order 2`: 3-point delta_xixi, one Thomas solve per mode. The matrix is constant
  in xi, so factor it once per mode at setup.
- `poisson_order 4` (Phase 6): compact (Numerov) operator
  delta^2 - lambda (1 + h^2/12 delta^2) applied to psi, with right-hand side (1 + h^2/12 delta^2) f. It is still
  tridiagonal and fourth order in xi. This is the analogue of the existing compact
  9-point option.
- The xi problem is Dirichlet at both ends, so m = 0 is non-singular.
- OpenMP: parallelise over modes.
- Why not DST in xi: the inhomogeneous outer Dirichlet values and any
  future non-conformal stretching r(xi), with variable coefficients, fit a banded solve per mode directly.

## 2. Code structure

Extend the existing solver; do not fork it. RK4, Euler, OpenMP loops, the backend and
the time loop all carry over.

**`include/polar.h`, `src/polar.c` (new):**

- `polar_grid`: R, r_inf, d xi, `r[i]`, `1/r[i]`, `1/r[i]^2`, `cos/sin theta_j`,
  outer BC type, far-field state (lagged C_D, `psi_out[j]`). It is owned by the caller
  and must outlive the solver, like the operators.
- `polar_grid polar_grid_make(int ntheta, int nr, double R, double r_inf)`.
- `polar_poisson_setup/solve/cleanup`, holding the FFTW plans with `fftw_plans_hold()`
  and `fftw_plans_release()`.
- `polar_set_wall_vorticity(w, psi, u, v, grid, closure)`, for the cylinder row and the
  outer row.
- `polar_apply_bc(u, v, grid)`: u = v = 0 on row 0, far-field velocity on row ny-1.
- `polar_farfield_update(grid, C_D)`: recompute `psi_out` and `omega_out`.
- Diagnostics (Section 3).
- `polar_time_step_limits()`.
- `printvtk_polar()`: VTK `STRUCTURED_GRID` with Cartesian points
  x = r cos theta, y = r sin theta. Repeat column 0 at theta = 2 pi so that ParaView closes the
  annulus. It also writes the velocity as Cartesian (u_x, u_y).

**`solver_config` (`include/fluiddyn.h`):** add `int geometry` (0 Cartesian,
1 polar) and `const struct polar_grid *polar`. Polar runs set `periodic = 0`.

**Hooks in `src/fluiddyn.c`.** Each is one `if (cfg.geometry == 1)` branch calling into
`polar.c`:

| where | change |
|---|---|
| `rk4_alloc()` | polar Poisson plans instead of `fft_setup()` |
| `solve_psi()` | right-hand side is -r_i^2 w, and the polar solve uses outer Dirichlet `psi_out` |
| `velocity_from_vorticity()` | u = -(1/r_i) DY psi, v = (1/r_i) DX psi (sign trap above) |
| `psi_is_current()` | compare **all** columns. Today `edge = 1` skips columns 0 and nx-1, which are interior in theta. As written, a stale psi would be reused. |
| `dwdt()`, `euler()` | per-row metric factors in the right-hand side loop. `euler()` needs the metric passed in, or the polar branch builds `out` and does w += dt*out. |
| `set_wall_vorticity()`, `apply_wall_bc()`, `step()` | polar versions: two boundary rows, no left or right walls. Includes the final refresh at the end of `step()`. |
| `max_stable_dt()` / `time_step_limits()` | `polar_time_step_limits()`. Viscous rate = (lambda_max(D_xixi) + lambda_max(D_thth)) / (Re R^2), Courant length R min(d xi, d theta), u_max = 2 (the potential-flow peak). The current versions assume uniform h and would accept an unstable dt. |

**Backend (`src/backend.c`):** a polar config always runs on the CPU. If the GPU is
asked for, print the same fallback message as today. This is mandatory: otherwise a
CUDA build silently runs Cartesian kernels on polar data. `backend_continuity()` and
`backend_integrals()` get polar branches with the area weight r^2 d xi d theta and the
continuity form above.

**`src/main.c`:**
- `--case cylinder`, with Re defaulting to 40 for this case. Use the "negative: not given"
  pattern already in `main()`.
- `--ntheta N`, `--nr N` (default from square cells), `--r-inf R` (default 40),
  `--outer-bc freestream|farfield`, `--steady-tol TOL` (stop when
  max|omega^{n+1} - omega^n| / dt < TOL), `--diag-interval N`.
- Reject the options that do not apply (forcing, `--operators fourier`, `--velocity-order 4`) the way `main()`
  already refuses mismatched options.
- Initial condition: impulsive start, omega = 0. The first Poisson solve gives potential
  flow and the wall closure creates the vortex sheet.
- Outputs:
  - `output/cylinder_history.csv`: t, C_D, C_Dp, C_Df, L_w, theta_s, residual
  - `output/cylinder_surface.csv`: theta, omega_w, p_w, C_p
  - `output/cylinder_centreline.csv`: r, u_x on theta = 0 and pi
  - a summary block printed at the end
  - VTK frames

## 3. Diagnostics: the validation quantities

1. **C_D from the cylinder surface.** The paper uses the momentum flux, eq. (1), but
   that needs interior pressure, which psi-omega does not have.
   - On the no-slip wall, theta-momentum gives dp/dtheta = (1/Re) omega_xi at r = R.
     So p_w(theta) is the integral of (1/Re) omega_xi in theta. Do it spectrally by FFT,
     dividing by i m. The m = 0 part of dp/dtheta should be zero; report it as a consistency check.
   - C_Dp = -2R times the integral over 0..2 pi of p_w cos(theta),
     C_Df = -(2R/Re) times the integral of omega_w sin(theta), and C_D = C_Dp + C_Df.
     Trapezoidal sums are spectrally accurate for periodic integrands.
   - Optional cross-check: the momentum deficit across the wake at a few radii.
2. **Separation angle:** the zero of omega_w(theta) on (0, pi), located by
   Fourier interpolation of omega_w and a root bracket. Report 180 deg - theta_zero.
3. **Wake length:** the zero crossing of u_r on the theta = 0 column (u_r < 0 in
   the bubble), by cubic interpolation in xi. L_w = r_zero - R.
4. **Recirculation centre:** the extremum of psi in 0 < theta < pi inside the bubble, by
   a local quadratic fit in (xi, theta), converted to Cartesian. Use Coutanceau &
   Bouard's convention: a is the streamwise distance from the rear stagnation point, b
   the spacing between the two centres. Confirm this against the sketch in Table 1
   before comparing.
5. **Centreline velocity** u_x(r, theta = 0) out to r_inf. Check that it is monotonic
   past the bubble (the paper's Fig. 1 criterion for the outer boundary condition).
6. **Point value at (3, 3):** bicubic interpolation in (xi, theta) of u_x, u_y, omega.
7. **Symmetry:** max |psi(theta) + psi(-theta)|. The scheme preserves the mirror
   symmetry, so this should stay at round-off. Re = 40 is below the Hopf bifurcation
   (Re about 47), so the steady symmetric state is the attractor.
8. **Steady residual:** max|omega^{n+1} - omega^n| / dt, plus the C_D history.

## 4. Verification (before any cylinder run)

These go in `tests/test_solver.c` and `tests/mms.c`, following their existing style:

1. **Operators:** `skronecker` of a periodic theta operator with a wall xi operator. Check it is exact on
   polynomials in xi times Fourier modes in theta, and that DX DY = DY DX to round-off.
2. **Polar Poisson, exact discrete test:** omega = 0 with psi(R) = 0 and psi(r_inf) =
   r_inf sin theta has the exact solution psi = A (r - R^2/r) sin theta, with
   A = r_inf / (r_inf - R^2/r_inf). Check the error order in xi: 2 for `poisson_order 2`, 4 for order 4.
3. **Polar Poisson, manufactured:** psi = (r - R)^2 (r_inf - r) (sin theta +
   0.5 cos 3 theta). Check the order under refinement in xi and theta, and that the solve
   agrees with an iterative or dense solve on a small grid.
4. **Full-solver MMS in polar.** Extend `mms_case` with `geometry`:
   - psi = g(t) f(r) T(theta), with f = (r - R)^2 (r_inf - r)^2 (no slip at both radii)
     and T a sum of 2 to 3 modes, so that the nonlinear term is non-zero.
   - Get the source in closed form, as `mms.c` does now. Generate the expressions once
     offline, for example with sympy, and paste them in.
   - Check the observed spatial order for orders 2, 4 and 6 with each closure. This catches the sign trap, the metric factors
     and the wall closure.
   - Reuse `test_temporal_order` for RK4 time order on this case.
5. **Continuity** at round-off, and the **Poisson reuse** count (the `psi_is_current` fix).
6. **Diagnostics on synthetic fields:**
   - C_Df on omega_w = c sin theta has a known integral.
   - The pressure integration on a pure mode.
   - The L_w and theta_s root finders on analytic profiles.
7. **Backend guard:** a CUDA build reports CPU for a polar config.
8. Optional: **circular Couette flow** (outer cylinder rotating), an exact steady
   solution. It exercises m = 0 and a moving-wall Briley term.

## 5. Validation campaign at Re = 40

Run to steady state with the residual below 1e-8, or 1e-6 for exploratory runs.

| Run set | Grids N_theta x N_r | r_inf | Outer BC | Purpose |
|---|---|---|---|---|
| A | 128x90, 256x180, 512x358 | 40 | freestream | grid convergence, observed order and Richardson-extrapolated C_D, L_w, theta_s |
| B | 256x(scaled), 512x(scaled) | 30, 40, 50 | farfield | match Table 2 layout, r_inf sensitivity |
| C | 512 | 40 | freestream vs farfield | reproduce Fig. 1: inflection with Dirichlet, monotone with far field |
| D | 512 | 40 | farfield | closure (Briley vs velocity), FD order 2/4/6, poisson_order 2/4 ablation |

**Acceptance targets** (finest grid, `farfield`, r_inf = 40):
- C_D within 1% of 1.4931
- L_w/D within 1% of 2.2360
- theta_s within 0.3 deg of 126.39
- a/D and b/D within 0.01 of 0.71 and 0.59
- u_x, u_y at (3, 3) within 1e-3
- centreline monotone
- symmetry error at round-off

The `freestream` runs should fall inside the literature scatter of Table 1, biased
by blockage. A miss of more than 2% against the far-field target is a bug to find, not a
tolerance to loosen. The paper's own coarse-grid spread at r_inf = 40 is 1.5%.

**Cost (estimate, check on the first run).** The step is limited by viscosity at the wall:
RK4 dt ~ 14 (R d xi)^2. That is about 2e-3 at N_theta = 256 and 5e-4 at 512. Reaching steady state
takes t ~ 100 to 200 D/U, so roughly 1e5 to 4e5 steps on the fine grid. The paper reports the same
azimuthal time-step penalty in Section 5. If that is too slow for the campaign, add the
optional **local pseudo-time** of Phase 6. The steady state is unchanged and the
convergence is much faster.

`tools/cylinder_validate.py` (new) reads the output CSVs and
`tests/reference/cylinder_re40.csv`. It prints the comparison tables and the observed orders, and
plots the centreline (Fig. 1 analogue), omega_w(theta) and C_p(theta). Use matplotlib and the
style in `tools/`.

**CI:** add to `tests/regression.sh`:
- 100 steps of a coarse cylinder (128x90, r_inf = 20) against stored output, like the
  cavity short run
- one steady coarse run with loose bounds: C_D in [1.45, 1.60], L_w in [2.1, 2.4]. It should
  take about 30 s; measure it first.

## 6. Phases (one PR each, branched from `main`, not `paper/draft`)

1. **Polar grid and Poisson.** `polar.c/h`, the FFT-theta / Thomas-xi solver, tests 4.2 and 4.3,
   `tests/reference/cylinder_re40.csv`.
2. **Solver integration.**
   - `geometry` in `solver_config`, the fluiddyn hooks including the
     `psi_is_current` column fix, the Briley and velocity closures on the cylinder, the
     `freestream` outer boundary, time-step limits, the backend guard
   - tests 4.1, 4.4, 4.5, 4.7
3. **Cylinder case.** CLI, initial condition, steady stop, the diagnostics of Section 3 (test 4.6),
   CSV and VTK output, the coarse regression in CI.
4. **Far-field boundary.** Schlichting plus mass correction for psi_out and omega_out, the
   lagged C_D update, run set C.
5. **Validation campaign.** Run sets A to D, `tools/cylinder_validate.py`, results in
   `docs/methodology.tex` (new section: polar formulation, boundary conditions,
   validation tables) and a README section with a vorticity figure.
6. **Optional:**
   - fourth-order compact xi Poisson
   - local pseudo-time for steady runs: scale each row's right-hand side by
     s_i = min(r_i^2/R^2, beta r_i/R). The first term is the local viscous limit, the second the local advective
     one. Steady runs only, and refused with unsteady output.
   - the paper's Test 1 with the reference field, if the authors send it: Dirichlet at
     r = 10.5D from `interpol.m` output, and the L2 error norm of eq. (3) under refinement
   - next cases: a rotating cylinder, and Re = 100 vortex shedding (Strouhal number about 0.165)

## 7. Risks and traps

- **Handedness sign** of the (theta, xi) layout. Catch it with the MMS test before anything else.
- **Four-wall assumptions** in `psi_is_current`, `apply_wall_bc`,
  `set_wall_vorticity`, `compute_integrals` (trapezoid weights), `printvtk`
  (uniform spacing) and the time-step limits. Each needs a polar branch. The
  `psi_is_current` one silently reuses a stale psi.
- **GPU path:** must refuse polar configs explicitly.
- **psi single-valuedness:** the outer psi must carry zero net flux. Otherwise there is a jump at
  theta = 2 pi and Gibbs ringing in the FFT solve.
- **C_D feedback** in the far-field boundary: keep it lagged and relaxed. Report C_D both
  with the boundary's C_D frozen and with it free, to show the coupling has converged.
- **Slow approach to steady state** near the Hopf threshold. Monitor the residual and the C_D
  history, not only the final time.
