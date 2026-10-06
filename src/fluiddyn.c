#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fluiddyn.h"
#include "linearalg.h"

void euler(mtrx w, mtrx dwdx, mtrx dwdy, mtrx d2wdx2, mtrx d2wdy2, mtrx u, mtrx v, double Re, double dt)
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
            MAt(w, i, j) = (-MAt(u, i, j) * MAt(dwdx, i, j) - MAt(v, i, j) * MAt(dwdy, i, j) + (1. / Re) * (MAt(d2wdx2, i, j) + MAt(d2wdy2, i, j))) * dt + MAt(w, i, j);
        }
    }
}

// ---------------------------------------------------------------------------
// RK4 time integration
// ---------------------------------------------------------------------------

#include <stdlib.h>
#include "poisson.h"

rk4_ctx rk4_alloc(int nx, int ny)
{
    rk4_ctx ctx;

    ctx.nx = nx; ctx.ny = ny;
    ctx.dwdx    = initm(ny, nx); ctx.dwdy    = initm(ny, nx);
    ctx.d2wdx2  = initm(ny, nx); ctx.d2wdy2  = initm(ny, nx);
    ctx.dpsidx  = initm(ny, nx); ctx.dpsidy  = initm(ny, nx);
    ctx.psi     = initm(ny, nx); ctx.psi_scratch = initm(ny, nx);
    ctx.k1      = initm(ny, nx); ctx.k2      = initm(ny, nx);
    ctx.k3      = initm(ny, nx); ctx.k4      = initm(ny, nx);
    ctx.w_tmp   = initm(ny, nx); ctx.rhs     = initm(ny, nx);
    return ctx;
}

void rk4_free(rk4_ctx *ctx)
{
    freem(&ctx->dwdx);   freem(&ctx->dwdy);
    freem(&ctx->d2wdx2); freem(&ctx->d2wdy2);
    freem(&ctx->dpsidx); freem(&ctx->dpsidy);
    freem(&ctx->psi);    freem(&ctx->psi_scratch);
    freem(&ctx->k1);     freem(&ctx->k2);
    freem(&ctx->k3);     freem(&ctx->k4);
    freem(&ctx->w_tmp);  freem(&ctx->rhs);
}

// Solve nabla^2 psi = -w and recover u = dpsi/dy, v = -dpsi/dx.
static void velocity_from_vorticity(mtrx w, mtrx u, mtrx v, rk4_ctx *ctx)
{
    // Poisson solve: nabla^2 psi = -w
    negcpy(ctx->rhs, w);
    if (ctx->poisson_type == 1)
        poisson(ctx->rhs, ctx->psi, ctx->psi_scratch, ctx->dx, ctx->dy,
                ctx->poisson_max_it, ctx->poisson_tol);
    else if (ctx->poisson_type == 2)
        poisson_SOR(ctx->rhs, ctx->psi, ctx->psi_scratch, ctx->dx, ctx->dy,
                    ctx->poisson_max_it, ctx->poisson_tol, ctx->beta);
    else if (ctx->poisson_type == 3)
        poisson_FFT(ctx->rhs, ctx->psi, ctx->dx, ctx->dy);
    else
    {
        printf("** Error: valid Poisson solver types are 1, 2 or 3 **\n");
        exit(1);
    }

    // Recover u = dpsi/dy, v = -dpsi/dx
    spmv(*ctx->DY, ctx->psi.M, u.M);
    spmv(*ctx->DX, ctx->psi.M, v.M);
    negcpy(v, v);
}

// First and second derivatives of w into the workspace
static void derivatives(mtrx w, rk4_ctx *ctx)
{
    spmv(*ctx->DX,  w.M, ctx->dwdx.M);
    spmv(*ctx->DY,  w.M, ctx->dwdy.M);
    spmv(*ctx->DX2, w.M, ctx->d2wdx2.M);
    spmv(*ctx->DY2, w.M, ctx->d2wdy2.M);
}

// Evaluate dw/dt and update u, v consistent with w via Poisson solve.
// out = -u*(dw/dx) - v*(dw/dy) + (1/Re)*(d2w/dx2 + d2w/dy2)
void dwdt(mtrx w, mtrx u, mtrx v, mtrx out, rk4_ctx *ctx)
{
    int i;
    int nx = ctx->nx, ny = ctx->ny;

    derivatives(w, ctx);

    velocity_from_vorticity(w, u, v, ctx);

    // RHS: dw/dt = -u*dw/dx - v*dw/dy + (1/Re)*(d2w/dx2 + d2w/dy2)
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (ny * nx >= OMP_MIN_WORK)
#endif
    for (i = 0; i < ny; i++)
    {
        int j;
        for (j = 0; j < nx; j++)
            MAt(out, i, j) = - MAt(u, i, j) * MAt(ctx->dwdx,   i, j)
                             - MAt(v, i, j) * MAt(ctx->dwdy,   i, j)
                             + (1.0 / ctx->Re) * (MAt(ctx->d2wdx2, i, j)
                                                + MAt(ctx->d2wdy2, i, j));
    }
}

// Classical RK4: w_{n+1} = w_n + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
// u and v are updated to be consistent with w_{n+1} on return.
void rk4(mtrx w, mtrx u, mtrx v, double dt, rk4_ctx *ctx)
{
    int i;
    int ny = ctx->ny, nx = ctx->nx;

    // k1 = f(w_n)
    dwdt(w, u, v, ctx->k1, ctx);

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
    dwdt(ctx->w_tmp, u, v, ctx->k2, ctx);

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
    dwdt(ctx->w_tmp, u, v, ctx->k3, ctx);

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
    dwdt(ctx->w_tmp, u, v, ctx->k4, ctx);

    // Combine: w_{n+1} = w_n + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (ny * nx >= OMP_MIN_WORK)
#endif
    for (i = 0; i < ny; i++)
    {
        int j;
        for (j = 0; j < nx; j++)
            MAt(w, i, j) += (dt / 6.0) * (MAt(ctx->k1, i, j)
                                       + 2.0 * MAt(ctx->k2, i, j)
                                       + 2.0 * MAt(ctx->k3, i, j)
                                           + MAt(ctx->k4, i, j));
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
    l.courant   = u_max > 0. ? max_co * h / u_max : HUGE_VAL;
    l.viscous   = max_stable_dt(dxx, dyy, Re, time_scheme);
    l.advection = time_scheme == 1 ? euler_advection_dt(Re, u_max) : HUGE_VAL;
    l.accept    = fmin(l.courant, l.viscous);
    l.suggest   = fmin(l.accept, l.advection);
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
        MAt(u, 0, j)    = bc->u[2];  MAt(v, 0, j)    = bc->v[2];
        MAt(u, ny-1, j) = bc->u[3];  MAt(v, ny-1, j) = bc->v[3];
    }
    for (i = 0; i < ny; i++)
    {
        MAt(u, i, 0)    = bc->u[0];  MAt(v, i, 0)    = bc->v[0];
        MAt(u, i, nx-1) = bc->u[1];  MAt(v, i, nx-1) = bc->v[1];
    }
}

void step(mtrx w, mtrx u, mtrx v, double dt, int time_scheme, const wall_bc *bc, rk4_ctx *ctx)
{
    int i, j;
    int nx = ctx->nx, ny = ctx->ny;
    // dpsidx/dpsidy are not needed until velocity recovery — use them as scratch
    mtrx dvdx = ctx->dpsidx, dudy = ctx->dpsidy;

    apply_wall_bc(u, v, bc);

    // Vorticity BCs: w = dv/dx - du/dy evaluated at boundaries
    spmv(*ctx->DY, u.M, dudy.M);
    spmv(*ctx->DX, v.M, dvdx.M);

    for (j = 0; j < nx; j++)
    {
        MAt(w, 0, j)    = MAt(dvdx, 0, j)    - MAt(dudy, 0, j);
        MAt(w, ny-1, j) = MAt(dvdx, ny-1, j) - MAt(dudy, ny-1, j);
    }
    for (i = 0; i < ny; i++)
    {
        MAt(w, i, 0)    = MAt(dvdx, i, 0)    - MAt(dudy, i, 0);
        MAt(w, i, nx-1) = MAt(dvdx, i, nx-1) - MAt(dudy, i, nx-1);
    }

    if (time_scheme == 1)
    {
        // Euler: single RHS evaluation, then one Poisson solve
        derivatives(w, ctx);

        euler(w, ctx->dwdx, ctx->dwdy, ctx->d2wdx2, ctx->d2wdy2, u, v, ctx->Re, dt);
        velocity_from_vorticity(w, u, v, ctx);
    }
    else
    {
        // RK4: four RHS evaluations, each with a Poisson solve
        // u and v are updated to be consistent with w on return
        rk4(w, u, v, dt, ctx);
    }
}
