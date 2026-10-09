#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fluiddyn.h"
#include "linearalg.h"
#include "fourier.h"

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
    ctx.w_solved = initm(ny, nx);
    ctx.psi_valid = 0;
    ctx.solves = 0;
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
    if ((cfg->fourier && !cfg->periodic) || (cfg->dealias && !cfg->fourier))
    {
        printf("** Error: the Fourier operators need a periodic grid, and dealiasing the Fourier operators **\n");
        exit(1);
    }
    // The operators in Fourier space, or the periodic Poisson solver for the
    // sparse ones
    ctx.fourier = cfg->fourier ? fourier_setup(cfg) : NULL;
    ctx.periodic = cfg->periodic && !cfg->fourier ? periodic_setup(nx, ny, cfg->DX2, cfg->DY2) : NULL;
    ctx.source = cfg->vorticity_source ? initm(ny, nx) : (mtrx){0};
    ctx.kolmogorov = kolmogorov_rows(&cfg->forcing, ny, cfg->dy, cfg->dy * (cfg->periodic ? ny : ny - 1));
    if (cfg->forcing.random_rate > 0.0 && !cfg->periodic)
    {
        printf("** Error: random forcing needs a periodic grid **\n");
        exit(1);
    }
    ctx.hyp1 = cfg->forcing.hyperviscosity > 0.0 ? initm(ny, nx) : (mtrx){0};
    ctx.hyp2 = cfg->forcing.hyperviscosity > 0.0 ? initm(ny, nx) : (mtrx){0};
    ctx.uw = cfg->advection == 1 ? initm(ny, nx) : (mtrx){0};
    ctx.vw = cfg->advection == 1 ? initm(ny, nx) : (mtrx){0};
    if (cfg->advection == 1 && !cfg->periodic)
    {
        printf("** Error: the skew-symmetric nonlinear term needs a periodic grid **\n");
        exit(1);
    }
    if (cfg->forcing.hyperviscosity > 0.0 && cfg->forcing.hyper_order < 2)
    {
        printf("** Error: the hyperviscosity order must be at least 2 **\n");
        exit(1);
    }
    ctx.kicks = NULL;
    if (cfg->forcing.random_rate > 0.0)
    {
        periodic_symbols sym;
        periodic_symbols_of(cfg, &sym);
        ctx.kicks = random_forcing_setup(&cfg->forcing, nx, ny, cfg->dx, cfg->dy, cfg->dt, &sym);
        periodic_symbols_free(&sym);
    }
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
    freem(&ctx->w_solved);
    fft_cleanup(ctx->fft);
    ctx->fft = NULL;
    periodic_cleanup(ctx->periodic);
    ctx->periodic = NULL;
    fourier_free(ctx->fourier);
    ctx->fourier = NULL;
    if (ctx->source.M) freem(&ctx->source);
    free(ctx->kolmogorov);
    ctx->kolmogorov = NULL;
    random_forcing_free(ctx->kicks);
    ctx->kicks = NULL;
    if (ctx->hyp1.M) freem(&ctx->hyp1);
    if (ctx->hyp2.M) freem(&ctx->hyp2);
    if (ctx->uw.M) freem(&ctx->uw);
    if (ctx->vw.M) freem(&ctx->vw);
}

// Does w hold, wherever the Poisson solve reads it, the values ctx->psi was
// solved for? The FFT solvers read all of w on a periodic grid and its
// interior with walls (the compact right-hand side extrapolates its wall
// values from the interior). The iterative solvers start from the previous
// psi, so a second solve may converge further: never skipped for them.
static int psi_is_current(mtrx w, const rk4_ctx *ctx)
{
    int i, nx = ctx->cfg.nx, ny = ctx->cfg.ny, edge = ctx->cfg.periodic ? 0 : 1;

    if (!ctx->psi_valid || (!ctx->cfg.periodic && ctx->cfg.poisson_type != 3)) return 0;
    for (i = edge; i < ny - edge; i++)
        if (memcmp(&MAt(w, i, edge), &MAt(ctx->w_solved, i, edge), (size_t)(nx - 2 * edge) * sizeof(double)) != 0)
            return 0;
    return 1;
}

// Solve nabla^2 psi = -w and recover u = dpsi/dy, v = -dpsi/dx. The solve is
// skipped when psi is already that of w (psi_is_current()): RK4's first stage
// would otherwise repeat the last solve of the step before. The velocity is
// recomputed from psi either way. remember: record w for that comparison; only
// the solves of the step's own w that the next RK4 first stage may repeat do,
// so the stage solves and runs with random kicks (which change w every step)
// do not pay for a copy of w.
static void solve_psi(mtrx w, int remember, rk4_ctx *ctx)
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
    // Only a solve that a later call can reuse is worth a copy of w (never
    // one of the iterative solvers, psi_is_current())
    remember = remember && (ctx->cfg.periodic || ctx->cfg.poisson_type == 3);
    if (remember) mtrxcpy(ctx->w_solved, w);
    ctx->psi_valid = remember;
    ctx->solves++;
}

static void velocity_from_vorticity(mtrx w, mtrx u, mtrx v, int remember, rk4_ctx *ctx)
{
    if (ctx->fourier)
    {
        // The solve and the velocity in one pass through Fourier space
        if (psi_is_current(w, ctx))
            fourier_velocity(ctx->fourier, ctx->psi.M, u.M, v.M);
        else
        {
            fourier_poisson(ctx->fourier, w.M, ctx->psi.M, u.M, v.M);
            if (remember) mtrxcpy(ctx->w_solved, w);
            ctx->psi_valid = remember;
            ctx->solves++;
        }
        return;
    }
    if (!psi_is_current(w, ctx)) solve_psi(w, remember, ctx);

    // Recover u = dpsi/dy, v = -dpsi/dx
    spmv(ctx->cfg.DYv ? *ctx->cfg.DYv : *ctx->cfg.DY, ctx->psi.M, u.M);
    spmv(ctx->cfg.DXv ? *ctx->cfg.DXv : *ctx->cfg.DX, ctx->psi.M, v.M);
    negcpy(v, v);
}

// First and second derivatives of w into the workspace
static void derivatives(mtrx w, rk4_ctx *ctx)
{
    if (ctx->fourier)
    {
        // The Laplacian whole into d2wdx2; d2wdy2 stays zero
        fourier_derivatives(ctx->fourier, w.M, ctx->dwdx.M, ctx->dwdy.M, ctx->d2wdx2.M);
        return;
    }
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

// out += -drag w + the Kolmogorov source (forcing.h)
// out += -drag w + the Kolmogorov source - nu_h (-L)^p w - alpha_h psi
// (forcing.h). psi must be that of w: dwdt() and step() solve for it first.
static void add_forcing_terms(mtrx out, mtrx w, rk4_ctx *ctx)
{
    int i, k, n = ctx->cfg.nx * ctx->cfg.ny, nx = ctx->cfg.nx, ny = ctx->cfg.ny;
    const forcing_config *fc = &ctx->cfg.forcing;
    double drag = fc->drag;

    if (fc->hyperviscosity > 0.0 && ctx->fourier)
    {
        fourier_power(ctx->fourier, w.M, fc->hyper_order, ctx->hyp1.M);
        for (k = 0; k < n; k++)
            out.M[k] -= fc->hyperviscosity * ctx->hyp1.M[k];
    }
    else if (fc->hyperviscosity > 0.0)
    {
        // (-L)^p w by p applications of -L = -(DX2 + DY2), alternating between
        // hyp1 and hyp2; rhs is free scratch outside the Poisson solve
        double *a = w.M, *bufs[2] = {ctx->hyp1.M, ctx->hyp2.M};
        for (i = 0; i < fc->hyper_order; i++)
        {
            double *b = bufs[i % 2];
            spmv(*ctx->cfg.DX2, a, b);
            spmv(*ctx->cfg.DY2, a, ctx->rhs.M);
            for (k = 0; k < n; k++)
                b[k] = -(b[k] + ctx->rhs.M[k]);
            a = b;
        }
        for (k = 0; k < n; k++)
            out.M[k] -= fc->hyperviscosity * a[k];
    }
    if (fc->hypodrag > 0.0)
        for (k = 0; k < n; k++)
            out.M[k] -= fc->hypodrag * ctx->psi.M[k];
    if (drag == 0.0 && !ctx->kolmogorov) return;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (nx * ny >= OMP_MIN_WORK)
#endif
    for (i = 0; i < ny; i++)
        for (int j = 0; j < nx; j++)
            MAt(out, i, j) += -drag * MAt(w, i, j) + (ctx->kolmogorov ? ctx->kolmogorov[i] : 0.0);
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
    velocity_from_vorticity(w, u, v, 0, ctx);
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
    if (ctx->cfg.advection == 1)
        skew_correction(&ctx->cfg, ctx->fourier, u.M, v.M, w.M, ctx->dwdx.M, ctx->dwdy.M, out.M, ctx->uw.M,
                        ctx->vw.M);

    const double *f = vorticity_source_at(t, ctx);
    if (f)
        for (i = 0; i < nx * ny; i++)
            out.M[i] += f[i];
    add_forcing_terms(out, w, ctx);
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

    // With dealiasing, w_{n+1} without the modes the 2/3 rule cuts
    if (ctx->fourier && ctx->cfg.dealias) fourier_filter(ctx->fourier, w.M);

    // Final Poisson solve so u, v are consistent with w_{n+1}; the next
    // step's first stage can reuse it unless a kick changes w first
    velocity_from_vorticity(w, u, v, !ctx->kicks, ctx);
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

void skew_correction(const solver_config *cfg, struct fourier_ops *four, const double *u, const double *v,
                     const double *w, const double *wx, const double *wy, double *out, double *s1, double *s2)
{
    int k, n = cfg->nx * cfg->ny;

    // Every read of wx and wy comes before the first write to s1 and s2
    for (k = 0; k < n; k++)
        out[k] += 0.5 * (u[k] * wx[k] + v[k] * wy[k]);
    if (four)
    {
        // DX(u w) + DY(v w) in one pass through Fourier space
        for (k = 0; k < n; k++)
        {
            s1[k] = u[k] * w[k];
            s2[k] = v[k] * w[k];
        }
        fourier_divergence(four, s1, s2, s1);
        for (k = 0; k < n; k++)
            out[k] -= 0.5 * s1[k];
        return;
    }
    for (k = 0; k < n; k++)
        s1[k] = u[k] * w[k];
    spmv(*cfg->DX, s1, s2);
    for (k = 0; k < n; k++)
        out[k] -= 0.5 * s2[k];
    for (k = 0; k < n; k++)
        s1[k] = v[k] * w[k];
    spmv(*cfg->DY, s1, s2);
    for (k = 0; k < n; k++)
        out[k] -= 0.5 * s2[k];
}

// Smallest non-zero eigenvalue of -A, A a periodic second-derivative
// operator of n points: its symbol at the lowest wavenumber, from the middle
// row with the column offsets taken modulo n
static double lowest_eigenvalue(const smtrx *A)
{
    int k, r = A->m / 2, n = A->m;
    double s = 0.0;
    for (k = A->row_ptr[r]; k < A->row_ptr[r + 1]; k++)
    {
        int off = A->col_idx[k] - r;
        if (2 * off > n) off -= n;
        if (2 * off < -n) off += n;
        s -= A->values[k] * cos(2.0 * 3.14159265358979323846 * off / n);
    }
    return s;
}

double damping_rate_max(const smtrx *dxx, const smtrx *dyy, double Re, const forcing_config *fc)
{
    // Largest eigenvalue of -(DX2 + DY2), from the middle rows. The interior
    // stencils are centered with coefficients of alternating sign, so the sum
    // of their magnitudes is the value of the stencil's symbol at the highest
    // grid frequency.
    double qmax = row_abs_sum(dxx, dxx->m / 2) + row_abs_sum(dyy, dyy->m / 2);
    double rate = qmax / Re;

    if (!fc) return rate;
    // sigma(Q) = Q/Re + nu_h Q^p + alpha + alpha_h/Q is convex for Q > 0, so
    // its maximum over the eigenvalues lies at the largest or the smallest
    rate += fc->drag;
    if (fc->hyperviscosity > 0.0) rate += fc->hyperviscosity * pow(qmax, fc->hyper_order);
    if (fc->hypodrag > 0.0)
    {
        double qmin = fmin(lowest_eigenvalue(dxx), lowest_eigenvalue(dyy));
        double low = qmin / Re + fc->drag + fc->hypodrag / qmin;
        if (fc->hyperviscosity > 0.0) low += fc->hyperviscosity * pow(qmin, fc->hyper_order);
        rate = fmax(rate + fc->hypodrag / qmax, low);
    }
    return rate;
}

double max_stable_dt_forced(const smtrx *dxx, const smtrx *dyy, double Re, const forcing_config *fc,
                            int time_scheme)
{
    // Forward Euler is stable for real eigenvalues in [-2, 0], classical RK4
    // down to -2.785293...
    return (time_scheme == 1 ? 2.0 : 2.785293563405282) / damping_rate_max(dxx, dyy, Re, fc);
}

double max_stable_dt(const smtrx *dxx, const smtrx *dyy, double Re, int time_scheme)
{
    return max_stable_dt_forced(dxx, dyy, Re, NULL, time_scheme);
}

double euler_advection_dt(double Re, double u_max)
{
    return u_max > 0.0 ? 2.0 / (Re * u_max * u_max) : HUGE_VAL;
}

dt_limits time_step_limits(const smtrx *dxx, const smtrx *dyy, double h, double Re, double u_max,
                           double max_co, int time_scheme)
{
    return time_step_limits_forced(dxx, dyy, h, Re, u_max, max_co, time_scheme, NULL);
}

dt_limits time_step_limits_forced(const smtrx *dxx, const smtrx *dyy, double h, double Re, double u_max,
                                  double max_co, int time_scheme, const forcing_config *fc)
{
    dt_limits l;
    l.courant = u_max > 0. ? max_co * h / u_max : HUGE_VAL;
    l.viscous = max_stable_dt_forced(dxx, dyy, Re, fc, time_scheme);
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
        velocity_from_vorticity(w, u, v, 1, ctx);
    if (!ctx->cfg.periodic)
    {
        apply_wall_bc(u, v, bc);
        set_wall_vorticity(w, u, v, ctx);
    }

    // Random forcing: a kick with new phases at the start of every step.
    // RK4's first stage derives the velocity from the kicked w; Euler uses
    // the velocity it is given, so it is updated here.
    if (ctx->kicks)
    {
        random_forcing_draw(ctx->kicks);
        random_forcing_add(ctx->kicks, w, ctx->cfg.dy);
        if (ctx->cfg.time_scheme == 1) velocity_from_vorticity(w, u, v, 0, ctx);
    }

    if (ctx->cfg.time_scheme == 1)
    {
        // Euler: single RHS evaluation, then one Poisson solve. The source,
        // the forcing and damping terms and the skew-symmetric correction go
        // in as one field (in w_tmp, which Euler does not otherwise use).
        const double *f = vorticity_source_at(ctx->cfg.t0 + (double)ctx->steps * dt, ctx);
        const forcing_config *fc = &ctx->cfg.forcing;
        // hypodrag needs the psi of w, which a first step does not have yet
        if (fc->hypodrag > 0.0 && ctx->steps == 0 && !ctx->kicks) velocity_from_vorticity(w, u, v, 0, ctx);
        derivatives(w, ctx);
        if (fc->drag != 0.0 || ctx->kolmogorov || fc->hyperviscosity > 0.0 || fc->hypodrag > 0.0 ||
            ctx->cfg.advection == 1)
        {
            int k, n = ctx->cfg.nx * ctx->cfg.ny;
            for (k = 0; k < n; k++)
                ctx->w_tmp.M[k] = f ? f[k] : 0.0;
            add_forcing_terms(ctx->w_tmp, w, ctx);
            if (ctx->cfg.advection == 1)
                skew_correction(&ctx->cfg, ctx->fourier, u.M, v.M, w.M, ctx->dwdx.M, ctx->dwdy.M, ctx->w_tmp.M,
                                ctx->uw.M, ctx->vw.M);
            f = ctx->w_tmp.M;
        }

        euler(w, ctx->dwdx, ctx->dwdy, ctx->d2wdx2, ctx->d2wdy2, u, v, f, ctx->cfg.Re, dt);
        if (ctx->fourier && ctx->cfg.dealias) fourier_filter(ctx->fourier, w.M);
        velocity_from_vorticity(w, u, v, 0, ctx);
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
