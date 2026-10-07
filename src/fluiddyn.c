#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fluiddyn.h"
#include "linearalg.h"

void euler(mtrx w, mtrx dwdx, mtrx dwdy, mtrx d2wdx2, mtrx d2wdy2, mtrx u, mtrx v, const double *f,
           double Re, double dt)
{
    int i;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (w.m * w.n >= OMP_MIN_WORK)
#endif
    for (i = 0; i < w.m; i++)
    {
        int j;
        for (j = 0; j < w.n; j++)
        {
            double rhs = -MAt(u, i, j) * MAt(dwdx, i, j) - MAt(v, i, j) * MAt(dwdy, i, j) + (1. / Re) * (MAt(d2wdx2, i, j) + MAt(d2wdy2, i, j));
            if (f) rhs += f[i * w.n + j];
            MAt(w, i, j) = rhs * dt + MAt(w, i, j);
        }
    }
}

// ---------------------------------------------------------------------------
// RK4 time integration
// ---------------------------------------------------------------------------

#include <stdlib.h>
#include "poisson.h"

rk4_ctx rk4_alloc(const solver_config *cfg)
{
    rk4_ctx ctx;
    int nx = cfg->nx, ny = cfg->ny;

    ctx.cfg = *cfg;

    ctx.dwdx = initm(ny, nx);
    ctx.dwdy = initm(ny, nx);
    ctx.d2wdx2 = initm(ny, nx);
    ctx.d2wdy2 = initm(ny, nx);
    ctx.psi = initm(ny, nx);
    ctx.psi_scratch = initm(ny, nx);
    ctx.k1 = initm(ny, nx);
    ctx.k2 = initm(ny, nx);
    ctx.k3 = initm(ny, nx);
    ctx.k4 = initm(ny, nx);
    ctx.w_tmp = initm(ny, nx);
    ctx.rhs = initm(ny, nx);
    if ((cfg->poisson_order != 2 && cfg->poisson_order != 4) || (cfg->wall_closure != 0 && cfg->wall_closure != 1))
    {
        printf("** Error: poisson_order must be 2 or 4 and wall_closure 0 or 1 **\n");
        exit(1);
    }
    if (cfg->poisson_order == 4 && cfg->poisson_type != 3 && !cfg->periodic)
    {
        printf("** Error: the fourth-order Poisson operator needs the FFT solver (poisson_type 3) **\n");
        exit(1);
    }
    if (cfg->periodic && cfg->poisson_type != 3)
    {
        printf("** Error: periodic boundaries need the FFT Poisson solver (poisson_type 3) **\n");
        exit(1);
    }
    ctx.fft = cfg->poisson_type == 3 && !cfg->periodic ? fft_setup(nx, ny) : NULL;
    ctx.periodic = cfg->periodic ? periodic_setup(nx, ny, cfg->DX2, cfg->DY2) : NULL;
    ctx.source = cfg->vorticity_source ? initm(ny, nx) : (mtrx){0};
    ctx.steps = 0;
    return ctx;
}

void rk4_free(rk4_ctx *ctx)
{
    freem(&ctx->dwdx);
    freem(&ctx->dwdy);
    freem(&ctx->d2wdx2);
    freem(&ctx->d2wdy2);
    freem(&ctx->psi);
    freem(&ctx->psi_scratch);
    freem(&ctx->k1);
    freem(&ctx->k2);
    freem(&ctx->k3);
    freem(&ctx->k4);
    freem(&ctx->w_tmp);
    freem(&ctx->rhs);
    fft_cleanup(ctx->fft);
    ctx->fft = NULL;
    periodic_cleanup(ctx->periodic);
    ctx->periodic = NULL;
    if (ctx->source.M) freem(&ctx->source);
}

// Solve nabla^2 psi = -w and recover u = dpsi/dy, v = -dpsi/dx.
static void velocity_from_vorticity(mtrx w, mtrx u, mtrx v, rk4_ctx *ctx)
{
    // Poisson solve: nabla^2 psi = -w
    negcpy(ctx->rhs, w);
    if (ctx->cfg.periodic)
        poisson_periodic(ctx->periodic, ctx->rhs, ctx->psi);
    else if (ctx->cfg.poisson_type == 1)
        poisson(ctx->rhs, ctx->psi, ctx->psi_scratch, ctx->cfg.dx, ctx->cfg.dy,
                ctx->cfg.poisson_max_it, ctx->cfg.poisson_tol);
    else if (ctx->cfg.poisson_type == 2)
        poisson_SOR(ctx->rhs, ctx->psi, ctx->psi_scratch, ctx->cfg.dx, ctx->cfg.dy,
                    ctx->cfg.poisson_max_it, ctx->cfg.poisson_tol, ctx->cfg.beta);
    else if (ctx->cfg.poisson_type == 3 && ctx->fft)
        poisson_FFT_order(ctx->fft, ctx->rhs, ctx->psi, ctx->cfg.dx, ctx->cfg.dy, ctx->cfg.poisson_order);
    else if (ctx->cfg.poisson_type == 3)
    {
        printf("** Error: the workspace has no FFT plans; it was not allocated by rk4_alloc() "
               "for poisson_type 3 **\n");
        exit(1);
    }
    else
    {
        printf("** Error: valid Poisson solver types are 1, 2 or 3 **\n");
        exit(1);
    }

    // Recover u = dpsi/dy, v = -dpsi/dx
    spmv(*ctx->cfg.DY, ctx->psi.M, u.M);
    spmv(*ctx->cfg.DX, ctx->psi.M, v.M);
    negcpy(v, v);
}

// First and second derivatives of w into the workspace
static void derivatives(mtrx w, rk4_ctx *ctx)
{
    spmv(*ctx->cfg.DX, w.M, ctx->dwdx.M);
    spmv(*ctx->cfg.DY, w.M, ctx->dwdy.M);
    spmv(*ctx->cfg.DX2, w.M, ctx->d2wdx2.M);
    spmv(*ctx->cfg.DY2, w.M, ctx->d2wdy2.M);
}

// Row r of A x
static double csr_row(const smtrx *A, const double *x, int r)
{
    int k;
    double sum = 0.0;
    for (k = A->row_ptr[r]; k < A->row_ptr[r + 1]; k++)
        sum += A->values[k] * x[A->col_idx[k]];
    return sum;
}

// Third-order wall vorticity from the stream function (Briley 1971). With
// psi = 0 along a wall, w = -d2psi/ds2 there, s the distance from the wall;
// from psi_0..psi_3 at s = 0, h, 2h, 3h and the wall-tangential velocity
// U = dpsi/ds(0):
//   w = (85 psi_0 - 108 psi_1 + 27 psi_2 - 4 psi_3) / (18 h^2) + 11 U / (3 h)
static double briley(double p0, double p1, double p2, double p3, double U, double h)
{
    return (85.0 * p0 - 108.0 * p1 + 27.0 * p2 - 4.0 * p3) / (18.0 * h * h) + 11.0 * U / (3.0 * h);
}

// Wall vorticity from psi (wall_closure 1). dpsi/ds is u at the bottom wall,
// -u at the top, -v at the left and v at the right. Corners take the bottom
// or top value, as in the velocity-based formula.
static void set_wall_vorticity_psi(mtrx w, const rk4_ctx *ctx)
{
    int i, j, nx = ctx->cfg.nx, ny = ctx->cfg.ny;
    mtrx p = ctx->psi;
    const wall_bc *bc = &ctx->cfg.bc;
    double dx = ctx->cfg.dx, dy = ctx->cfg.dy;

    for (j = 0; j < nx; j++)
    {
        MAt(w, 0, j) = briley(MAt(p, 0, j), MAt(p, 1, j), MAt(p, 2, j), MAt(p, 3, j), bc->u[2], dy);
        MAt(w, ny - 1, j) = briley(MAt(p, ny - 1, j), MAt(p, ny - 2, j), MAt(p, ny - 3, j), MAt(p, ny - 4, j),
                                   -bc->u[3], dy);
    }
    for (i = 1; i < ny - 1; i++)
    {
        MAt(w, i, 0) = briley(MAt(p, i, 0), MAt(p, i, 1), MAt(p, i, 2), MAt(p, i, 3), -bc->v[0], dx);
        MAt(w, i, nx - 1) = briley(MAt(p, i, nx - 1), MAt(p, i, nx - 2), MAt(p, i, nx - 3), MAt(p, i, nx - 4),
                                   bc->v[1], dx);
    }
}

// Wall vorticity: w = dv/dx - du/dy on the boundary nodes (wall_closure 0),
// or the third-order formula from psi (wall_closure 1), which reads ctx->psi
// and ignores u and v
static void set_wall_vorticity(mtrx w, mtrx u, mtrx v, const rk4_ctx *ctx)
{
    int i, j, nx = ctx->cfg.nx, ny = ctx->cfg.ny;
    const smtrx *DX = ctx->cfg.DX, *DY = ctx->cfg.DY;

    if (ctx->cfg.wall_closure == 1)
    {
        set_wall_vorticity_psi(w, ctx);
        return;
    }

    for (j = 0; j < nx; j++)
    {
        int b = j, t = (ny - 1) * nx + j; // bottom and top rows
        w.M[b] = csr_row(DX, v.M, b) - csr_row(DY, u.M, b);
        w.M[t] = csr_row(DX, v.M, t) - csr_row(DY, u.M, t);
    }
    for (i = 1; i < ny - 1; i++)
    {
        int l = i * nx, r = i * nx + nx - 1; // left and right columns
        w.M[l] = csr_row(DX, v.M, l) - csr_row(DY, u.M, l);
        w.M[r] = csr_row(DX, v.M, r) - csr_row(DY, u.M, r);
    }
}

// The vorticity source at time t into ctx->source; NULL without cfg.vorticity_source
static const double *vorticity_source_at(double t, rk4_ctx *ctx)
{
    if (!ctx->cfg.vorticity_source) return NULL;
    ctx->cfg.vorticity_source(t, ctx->source, ctx->cfg.source_data);
    return ctx->source.M;
}

// out = -u*(dw/dx) - v*(dw/dy) + (1/Re)*(d2w/dx2 + d2w/dy2) + f(t). Also sets
// u, v from w's interior and the wall velocities, and w's boundary from u, v.
void dwdt(mtrx w, mtrx u, mtrx v, mtrx out, double t, rk4_ctx *ctx)
{
    int i;
    int nx = ctx->cfg.nx, ny = ctx->cfg.ny;

    // Velocity of this stage, then the wall vorticity that goes with it. The
    // wall vorticity follows the velocity, so it has to be set again at every
    // RK4 stage: set once per step, it lags the interior and limits RK4 to
    // first order in time. A periodic grid has no walls.
    velocity_from_vorticity(w, u, v, ctx);
    if (!ctx->cfg.periodic)
    {
        apply_wall_bc(u, v, &ctx->cfg.bc);
        set_wall_vorticity(w, u, v, ctx);
    }

    derivatives(w, ctx);

    // RHS: dw/dt = -u*dw/dx - v*dw/dy + (1/Re)*(d2w/dx2 + d2w/dy2)
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (ny * nx >= OMP_MIN_WORK)
#endif
    for (i = 0; i < ny; i++)
    {
        int j;
        for (j = 0; j < nx; j++)
            MAt(out, i, j) = -MAt(u, i, j) * MAt(ctx->dwdx, i, j) - MAt(v, i, j) * MAt(ctx->dwdy, i, j) + (1.0 / ctx->cfg.Re) * (MAt(ctx->d2wdx2, i, j) + MAt(ctx->d2wdy2, i, j));
    }

    const double *f = vorticity_source_at(t, ctx);
    if (f)
        for (i = 0; i < nx * ny; i++)
            out.M[i] += f[i];
}

// Classical RK4: w_{n+1} = w_n + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
// On return u and v are the velocity of the interior of w_{n+1}. The wall
// entries of w_{n+1} are not meaningful (the combination advances them with
// the transport RHS at the walls); step() replaces them.
void rk4(mtrx w, mtrx u, mtrx v, rk4_ctx *ctx)
{
    double dt = ctx->cfg.dt, t = ctx->cfg.t0 + (double)ctx->steps * dt;
    int i;
    int ny = ctx->cfg.ny, nx = ctx->cfg.nx;

    // k1 = f(w_n)
    dwdt(w, u, v, ctx->k1, t, ctx);

    // k2 = f(w_n + dt/2 * k1)
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (ny * nx >= OMP_MIN_WORK)
#endif
    for (i = 0; i < ny; i++)
    {
        int j;
        for (j = 0; j < nx; j++)
            MAt(ctx->w_tmp, i, j) = MAt(w, i, j) + 0.5 * dt * MAt(ctx->k1, i, j);
    }
    dwdt(ctx->w_tmp, u, v, ctx->k2, t + 0.5 * dt, ctx);

    // k3 = f(w_n + dt/2 * k2)
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (ny * nx >= OMP_MIN_WORK)
#endif
    for (i = 0; i < ny; i++)
    {
        int j;
        for (j = 0; j < nx; j++)
            MAt(ctx->w_tmp, i, j) = MAt(w, i, j) + 0.5 * dt * MAt(ctx->k2, i, j);
    }
    dwdt(ctx->w_tmp, u, v, ctx->k3, t + 0.5 * dt, ctx);

    // k4 = f(w_n + dt * k3)
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (ny * nx >= OMP_MIN_WORK)
#endif
    for (i = 0; i < ny; i++)
    {
        int j;
        for (j = 0; j < nx; j++)
            MAt(ctx->w_tmp, i, j) = MAt(w, i, j) + dt * MAt(ctx->k3, i, j);
    }
    dwdt(ctx->w_tmp, u, v, ctx->k4, t + dt, ctx);

    // Combine: w_{n+1} = w_n + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (ny * nx >= OMP_MIN_WORK)
#endif
    for (i = 0; i < ny; i++)
    {
        int j;
        for (j = 0; j < nx; j++)
            MAt(w, i, j) += (dt / 6.0) * (MAt(ctx->k1, i, j) + 2.0 * MAt(ctx->k2, i, j) + 2.0 * MAt(ctx->k3, i, j) + MAt(ctx->k4, i, j));
    }

    // Final Poisson solve so u, v are consistent with w_{n+1}
    velocity_from_vorticity(w, u, v, ctx);
}

// ---------------------------------------------------------------------------
// Time-step limit
// ---------------------------------------------------------------------------

// Sum of |coefficients| in row r of A
static double row_abs_sum(const smtrx *A, int r)
{
    int k;
    double s = 0.0;
    for (k = A->row_ptr[r]; k < A->row_ptr[r + 1]; k++)
        s += fabs(A->values[k]);
    return s;
}

double max_stable_dt(const smtrx *dxx, const smtrx *dyy, double Re, int time_scheme)
{
    // Largest |eigenvalue| of each second-derivative operator, taken from the
    // middle row. The interior stencils are centered with coefficients of
    // alternating sign, so the sum of their magnitudes is the value of the
    // stencil's symbol at the highest grid frequency.
    double lambda = (row_abs_sum(dxx, dxx->m / 2) + row_abs_sum(dyy, dyy->m / 2)) / Re;

    // Forward Euler is stable for real eigenvalues in [-2, 0], classical RK4
    // down to -2.785293...
    return (time_scheme == 1 ? 2.0 : 2.785293563405282) / lambda;
}

double euler_advection_dt(double Re, double u_max)
{
    return u_max > 0.0 ? 2.0 / (Re * u_max * u_max) : HUGE_VAL;
}

dt_limits time_step_limits(const smtrx *dxx, const smtrx *dyy, double h, double Re, double u_max,
                           double max_co, int time_scheme)
{
    dt_limits l;
    l.courant = u_max > 0. ? max_co * h / u_max : HUGE_VAL;
    l.viscous = max_stable_dt(dxx, dyy, Re, time_scheme);
    l.advection = time_scheme == 1 ? euler_advection_dt(Re, u_max) : HUGE_VAL;
    l.accept = fmin(l.courant, l.viscous);
    l.suggest = fmin(l.accept, l.advection);
    return l;
}

double round_down_3(double x)
{
    if (!(x > 0.0) || !isfinite(x)) return x;
    double scale = pow(10.0, 2.0 - floor(log10(x)));
    double r = floor(x * scale) / scale;
    // floor(log10) can be off by one next to a power of ten
    return r > x ? (floor(x * scale) - 1.0) / scale : r;
}

// ---------------------------------------------------------------------------
// Full timestep
// ---------------------------------------------------------------------------

void apply_wall_bc(mtrx u, mtrx v, const wall_bc *bc)
{
    int i, j;
    int ny = u.m, nx = u.n; // ny rows (y), nx columns (x)

    for (j = 0; j < nx; j++)
    {
        MAt(u, 0, j) = bc->u[2];
        MAt(v, 0, j) = bc->v[2];
        MAt(u, ny - 1, j) = bc->u[3];
        MAt(v, ny - 1, j) = bc->v[3];
    }
    for (i = 0; i < ny; i++)
    {
        MAt(u, i, 0) = bc->u[0];
        MAt(v, i, 0) = bc->v[0];
        MAt(u, i, nx - 1) = bc->u[1];
        MAt(v, i, nx - 1) = bc->v[1];
    }
}

void step(mtrx w, mtrx u, mtrx v, rk4_ctx *ctx)
{
    double dt = ctx->cfg.dt;
    const wall_bc *bc = &ctx->cfg.bc;

    // Vorticity BCs: w = dv/dx - du/dy evaluated at boundaries (walls only).
    // The formula from psi needs psi, which a first step does not have yet.
    if (!ctx->cfg.periodic && ctx->cfg.wall_closure == 1 && ctx->steps == 0)
        velocity_from_vorticity(w, u, v, ctx);
    if (!ctx->cfg.periodic)
    {
        apply_wall_bc(u, v, bc);
        set_wall_vorticity(w, u, v, ctx);
    }

    if (ctx->cfg.time_scheme == 1)
    {
        // Euler: single RHS evaluation, then one Poisson solve
        derivatives(w, ctx);

        euler(w, ctx->dwdx, ctx->dwdy, ctx->d2wdx2, ctx->d2wdy2, u, v,
              vorticity_source_at(ctx->cfg.t0 + (double)ctx->steps * dt, ctx), ctx->cfg.Re, dt);
        velocity_from_vorticity(w, u, v, ctx);
    }
    else
    {
        // RK4: four RHS evaluations, each with a Poisson solve and its own
        // wall vorticity
        rk4(w, u, v, ctx);
    }

    // The update advanced the wall entries of w with the transport equation,
    // which does not hold there. Replace them with the wall vorticity of the
    // new velocity, so that the w the caller sees (and writes out) is
    // consistent with u and v. u and v themselves stay the velocity of the
    // Poisson solution, whose divergence the continuity check measures; the
    // wall velocities are imposed on copies in the stage buffers.
    if (!ctx->cfg.periodic)
    {
        mtrxcpy(ctx->k1, u);
        mtrxcpy(ctx->k2, v);
        apply_wall_bc(ctx->k1, ctx->k2, bc);
        set_wall_vorticity(w, ctx->k1, ctx->k2, ctx);
    }
    ctx->steps++;
}
