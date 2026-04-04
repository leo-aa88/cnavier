#include <stdio.h>
#include <stdlib.h>
#include "fluiddyn.h"
#include "linearalg.h"

void euler(mtrx w, mtrx dwdx, mtrx dwdy, mtrx d2wdx2, mtrx d2wdy2, mtrx u, mtrx v, double Re, double dt)
{
    int i, j;

    for (i = 0; i < w.m; i++)
    {
        for (j = 0; j < w.n; j++)
        {
            MAt(w, i, j) = (-MAt(u, i, j) * MAt(dwdx, i, j) - MAt(v, i, j) * MAt(dwdy, i, j) + (1. / Re) * (MAt(d2wdx2, i, j) + MAt(d2wdy2, i, j))) * dt + MAt(w, i, j);
        }
    }
}

mtrx continuity(mtrx dudx, mtrx dvdy)
{
    int i, j;
    mtrx temp;
    temp = initm(dudx.m, dudx.n);

    for (i = 0; i < temp.m; i++)
    {
        for (j = 0; j < temp.n; j++)
        {
            MAt(temp, i, j) = MAt(dudx, i, j) + MAt(dvdy, i, j);
        }
    }
    return temp;
}

mtrx vorticity(mtrx dudy, mtrx dvdx)
{
    int i, j;
    mtrx temp;
    temp = initm(dudy.m, dudy.n);

    for (i = 0; i < temp.m; i++)
    {
        for (j = 0; j < temp.n; j++)
        {
            MAt(temp, i, j) = MAt(dvdx, i, j) - MAt(dudy, i, j);
        }
    }
    return temp;
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
    ctx.dwdx    = initm(nx, ny); ctx.dwdy    = initm(nx, ny);
    ctx.d2wdx2  = initm(nx, ny); ctx.d2wdy2  = initm(nx, ny);
    ctx.dpsidx  = initm(nx, ny); ctx.dpsidy  = initm(nx, ny);
    ctx.psi     = initm(nx, ny); ctx.psi_scratch = initm(nx, ny);
    ctx.k1      = initm(nx, ny); ctx.k2      = initm(nx, ny);
    ctx.k3      = initm(nx, ny); ctx.k4      = initm(nx, ny);
    ctx.w_tmp   = initm(nx, ny);
    ctx.flat_w   = (double *)malloc((size_t)nx * ny * sizeof(double));
    ctx.flat_psi = (double *)malloc((size_t)nx * ny * sizeof(double));
    ctx.flat_tmp = (double *)malloc((size_t)nx * ny * sizeof(double));
    if (!ctx.flat_w || !ctx.flat_psi || !ctx.flat_tmp)
    {
        printf("** Error: insufficient memory for RK4 workspace **\n");
        exit(1);
    }
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
    freem(&ctx->w_tmp);
    free(ctx->flat_w);
    free(ctx->flat_psi);
    free(ctx->flat_tmp);
}

// Evaluate dw/dt and update u, v consistent with w via Poisson solve.
// out = -u*(dw/dx) - v*(dw/dy) + (1/Re)*(d2w/dx2 + d2w/dy2)
void dwdt(mtrx w, mtrx u, mtrx v, mtrx out, rk4_ctx *ctx)
{
    int i, j;
    int nx = ctx->nx, ny = ctx->ny;

    // Derivatives of w
    flatten(w, ctx->flat_w, nx, ny);
    spmv(*ctx->DX,  ctx->flat_w, ctx->flat_tmp); unflatten(ctx->flat_tmp, ctx->dwdx,   nx, ny);
    spmv(*ctx->DY,  ctx->flat_w, ctx->flat_tmp); unflatten(ctx->flat_tmp, ctx->dwdy,   nx, ny);
    spmv(*ctx->DX2, ctx->flat_w, ctx->flat_tmp); unflatten(ctx->flat_tmp, ctx->d2wdx2, nx, ny);
    spmv(*ctx->DY2, ctx->flat_w, ctx->flat_tmp); unflatten(ctx->flat_tmp, ctx->d2wdy2, nx, ny);

    // Poisson solve: nabla^2 psi = -w  =>  invert sign, solve, restore
    invsig(w);
    if (ctx->poisson_type == 1)
        poisson(w, ctx->psi, ctx->psi_scratch, ctx->dx, ctx->dy,
                ctx->poisson_max_it, ctx->poisson_tol);
    else if (ctx->poisson_type == 2)
        poisson_SOR(w, ctx->psi, ctx->psi_scratch, ctx->dx, ctx->dy,
                    ctx->poisson_max_it, ctx->poisson_tol, ctx->beta);
    else
        poisson_FFT(w, ctx->psi, ctx->dx, ctx->dy);
    invsig(w);

    // Recover u = dpsi/dy, v = -dpsi/dx
    flatten(ctx->psi, ctx->flat_psi, nx, ny);
    spmv(*ctx->DY, ctx->flat_psi, ctx->flat_tmp); unflatten(ctx->flat_tmp, ctx->dpsidy, nx, ny);
    spmv(*ctx->DX, ctx->flat_psi, ctx->flat_tmp); unflatten(ctx->flat_tmp, ctx->dpsidx, nx, ny);
    mtrxcpy(u, ctx->dpsidy);
    invsig(ctx->dpsidx);
    mtrxcpy(v, ctx->dpsidx);

    // RHS: dw/dt = -u*dw/dx - v*dw/dy + (1/Re)*(d2w/dx2 + d2w/dy2)
    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
            MAt(out, i, j) = - MAt(u, i, j) * MAt(ctx->dwdx,   i, j)
                             - MAt(v, i, j) * MAt(ctx->dwdy,   i, j)
                             + (1.0 / ctx->Re) * (MAt(ctx->d2wdx2, i, j)
                                                + MAt(ctx->d2wdy2, i, j));
}

// Classical RK4: w_{n+1} = w_n + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
// u and v are updated to be consistent with w_{n+1} on return.
void rk4(mtrx w, mtrx u, mtrx v, double dt, rk4_ctx *ctx)
{
    int i, j;
    int nx = ctx->nx, ny = ctx->ny;

    // k1 = f(w_n)
    dwdt(w, u, v, ctx->k1, ctx);

    // k2 = f(w_n + dt/2 * k1)
    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
            MAt(ctx->w_tmp, i, j) = MAt(w, i, j) + 0.5 * dt * MAt(ctx->k1, i, j);
    dwdt(ctx->w_tmp, u, v, ctx->k2, ctx);

    // k3 = f(w_n + dt/2 * k2)
    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
            MAt(ctx->w_tmp, i, j) = MAt(w, i, j) + 0.5 * dt * MAt(ctx->k2, i, j);
    dwdt(ctx->w_tmp, u, v, ctx->k3, ctx);

    // k4 = f(w_n + dt * k3)
    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
            MAt(ctx->w_tmp, i, j) = MAt(w, i, j) + dt * MAt(ctx->k3, i, j);
    dwdt(ctx->w_tmp, u, v, ctx->k4, ctx);

    // Combine: w_{n+1} = w_n + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
            MAt(w, i, j) += (dt / 6.0) * (MAt(ctx->k1, i, j)
                                       + 2.0 * MAt(ctx->k2, i, j)
                                       + 2.0 * MAt(ctx->k3, i, j)
                                           + MAt(ctx->k4, i, j));

    // Final Poisson solve so u, v are consistent with w_{n+1}
    dwdt(w, u, v, ctx->k1, ctx); // k1 reused as scratch — u, v updated as side effect
}
