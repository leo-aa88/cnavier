#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linearalg.h"
#include "poisson.h"

double error(mtrx u1, mtrx u2)
{
    double e = 0;
    int i, j;

    for (i = 0; i < u1.m; i++)
        for (j = 0; j < u1.n; j++)
            e += sqrt(pow(MAt(u2, i, j) - MAt(u1, i, j), 2));
    return e;
}

// Gauss-Seidel Poisson solver.
// u and u0 are pre-allocated by the caller; u holds the result on return.
int poisson(mtrx f, mtrx u, mtrx u0, double dx, double dy, int itmax, double tol)
{
    int i, k;
    int ny = f.m, nx = f.n; // ny rows (y), nx columns (x)
    double e;
    double dx2 = dx * dx, dy2 = dy * dy;
    double denom = 2.0 * (dx2 + dy2);

    zerosm(u);

    for (k = 0; k < itmax; k++)
    {
        mtrxcpy(u0, u);
#ifdef _OPENMP
        /* Red–black ordering: two parallel phases (sequential GS is not safe to omp parallel for). */
#pragma omp parallel for schedule(static) if (ny * nx >= OMP_MIN_WORK)
        for (i = 1; i < ny - 1; i++)
        {
            int j0 = (i & 1) ? 1 : 2;
            int j;
            for (j = j0; j < nx - 1; j += 2)
                MAt(u, i, j) = (dx2 * (MAt(u, i + 1, j) + MAt(u, i - 1, j)) + dy2 * (MAt(u, i, j + 1) + MAt(u, i, j - 1)) - dx2 * dy2 * MAt(f, i, j)) / denom;
        }
#pragma omp parallel for schedule(static) if (ny * nx >= OMP_MIN_WORK)
        for (i = 1; i < ny - 1; i++)
        {
            int j0 = (i & 1) ? 2 : 1;
            int j;
            for (j = j0; j < nx - 1; j += 2)
                MAt(u, i, j) = (dx2 * (MAt(u, i + 1, j) + MAt(u, i - 1, j)) + dy2 * (MAt(u, i, j + 1) + MAt(u, i, j - 1)) - dx2 * dy2 * MAt(f, i, j)) / denom;
        }
#else
        for (i = 1; i < ny - 1; i++)
            for (int j = 1; j < nx - 1; j++)
                MAt(u, i, j) = (dx2 * (MAt(u, i + 1, j) + MAt(u, i - 1, j)) + dy2 * (MAt(u, i, j + 1) + MAt(u, i, j - 1)) - dx2 * dy2 * MAt(f, i, j)) / denom;
#endif
        e = error(u, u0);
        if (e < tol)
        {
            printf("Poisson solved in %d iterations - RSS error: %E\n", k, e);
            return k;
        }
    }
    printf("Error: max iterations reached for Poisson solver.\n");
    exit(1);
}

double sor_beta(int nx, int ny, double dx, double dy)
{
    double cx = cos(PI / (nx - 1)), cy = cos(PI / (ny - 1));
    // With dx == dy the weights are equal; average directly so that square
    // grids keep the exact value they always had
    double rho = (dx == dy) ? 0.5 * (cx + cy)
                            : (dy * dy * cx + dx * dx * cy) / (dx * dx + dy * dy);
    return 2.0 / (1.0 + sqrt(1.0 - rho * rho));
}

// SOR Poisson solver.
// u and u0 are pre-allocated by the caller; u holds the result on return.
int poisson_SOR(mtrx f, mtrx u, mtrx u0, double dx, double dy, int itmax, double tol, double beta)
{
    int i, k;
    int ny = f.m, nx = f.n; // ny rows (y), nx columns (x)
    double e;
    double dx2 = dx * dx, dy2 = dy * dy;
    double denom = 2.0 * (dx2 + dy2);

    zerosm(u);

    for (k = 0; k < itmax; k++)
    {
        mtrxcpy(u0, u);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (ny * nx >= OMP_MIN_WORK)
        for (i = 1; i < ny - 1; i++)
        {
            int j0 = (i & 1) ? 1 : 2;
            int j;
            for (j = j0; j < nx - 1; j += 2)
                MAt(u, i, j) = beta * (dx2 * (MAt(u, i + 1, j) + MAt(u, i - 1, j)) + dy2 * (MAt(u, i, j + 1) + MAt(u, i, j - 1)) - dx2 * dy2 * MAt(f, i, j)) / denom + (1.0 - beta) * MAt(u0, i, j);
        }
#pragma omp parallel for schedule(static) if (ny * nx >= OMP_MIN_WORK)
        for (i = 1; i < ny - 1; i++)
        {
            int j0 = (i & 1) ? 2 : 1;
            int j;
            for (j = j0; j < nx - 1; j += 2)
                MAt(u, i, j) = beta * (dx2 * (MAt(u, i + 1, j) + MAt(u, i - 1, j)) + dy2 * (MAt(u, i, j + 1) + MAt(u, i, j - 1)) - dx2 * dy2 * MAt(f, i, j)) / denom + (1.0 - beta) * MAt(u0, i, j);
        }
#else
        for (i = 1; i < ny - 1; i++)
            for (int j = 1; j < nx - 1; j++)
                MAt(u, i, j) = beta * (dx2 * (MAt(u, i + 1, j) + MAt(u, i - 1, j)) + dy2 * (MAt(u, i, j + 1) + MAt(u, i, j - 1)) - dx2 * dy2 * MAt(f, i, j)) / denom + (1.0 - beta) * MAt(u0, i, j);
#endif
        e = error(u, u0);
        if (e < tol)
        {
            printf("Poisson SOR solved in %d iterations - RSS error: %E\n", k, e);
            return k;
        }
    }
    printf("Error: max iterations reached for Poisson SOR solver.\n");
    exit(1);
}

// ---------------------------------------------------------------------------
// FFT-based direct Poisson solver
// ---------------------------------------------------------------------------
// Solves the 5-point discretisation of nabla^2 u = f with u = 0 on the wall
// nodes (first and last row and column), the same problem the Gauss-Seidel
// and SOR solvers solve. The unknowns are the (nx-2) x (ny-2) interior nodes.
//
// The 2D Discrete Sine Transform (DST-I) of the interior diagonalises that
// operator: a DST-I of length m implies zeros at indices -1 and m, which are
// exactly the wall nodes. The solution is exact (to floating-point precision)
// in a single pass:
//   1. Forward DST-I of the interior right-hand side
//   2. Divide each mode by its eigenvalue
//   3. Inverse DST-I (= forward DST-I / (4*(mx+1)*(my+1)))
//
// FFTW's RODFT00 plan is the DST-I.
// ---------------------------------------------------------------------------

#include <fftw3.h>

// The 2D DST-I is done as two passes of 1D DST-I: along the rows, then along
// the columns. Each pass works in batches of DST_BATCH rows (or columns), one
// FFTW plan per batch size, and OpenMP builds spread the batches over the
// threads. Every row is therefore transformed by the same plan however many
// threads there are, so serial and OpenMP builds give bit-identical results.
// FFTW_ESTIMATE makes the plans themselves deterministic; FFTW_MEASURE times
// candidate plans and can pick a different one, with different round-off, on
// every run.
//
// DST_BATCH is a multiple of 8, so every batch starts 64 bytes after the
// previous one and has the alignment the plans were made for.
#define DST_BATCH 8

typedef struct
{
    fftw_plan full; // exactly DST_BATCH transforms, or NULL if count < DST_BATCH
    fftw_plan rest; // the last count % DST_BATCH transforms, or NULL if none
    int count;
    ptrdiff_t step; // offset between batches, in doubles
} dst_pass;

struct fft_solver
{
    int nx, ny;
    double *buf; // interior work buffer, (ny-2) rows of (nx-2), row-major
    dst_pass pass_rows, pass_cols;
};

// fftw_cleanup() may only run once no plan is left, i.e. after the last
// solver is freed
static int live_solvers = 0;

// count transforms of length len, elements stride apart, transforms dist apart
static dst_pass make_pass(double *buf, int len, int count, int stride, int dist)
{
    dst_pass p;
    fftw_r2r_kind kind = FFTW_RODFT00;

    p.count = count;
    p.step = (ptrdiff_t)DST_BATCH * dist;
    p.full = NULL;
    p.rest = NULL;
    if (count >= DST_BATCH)
        p.full = fftw_plan_many_r2r(1, &len, DST_BATCH, buf, NULL, stride, dist,
                                    buf, NULL, stride, dist, &kind, FFTW_ESTIMATE);
    if (count % DST_BATCH)
        p.rest = fftw_plan_many_r2r(1, &len, count % DST_BATCH, buf, NULL, stride, dist,
                                    buf, NULL, stride, dist, &kind, FFTW_ESTIMATE);
    if ((count >= DST_BATCH && !p.full) || (count % DST_BATCH && !p.rest))
    {
        printf("** Error: FFTW could not plan the sine transform **\n");
        exit(1);
    }
    return p;
}

static void free_pass(dst_pass *p)
{
    if (p->full) fftw_destroy_plan(p->full);
    if (p->rest) fftw_destroy_plan(p->rest);
}

static void run_pass(const dst_pass *p, double *buf, int parallel)
{
    int b, batches = (p->count + DST_BATCH - 1) / DST_BATCH;
    (void)parallel;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (parallel)
#endif
    for (b = 0; b < batches; b++)
    {
        // Every batch but a short last one holds exactly DST_BATCH transforms
        double *x = buf + b * p->step;
        fftw_execute_r2r(p->count - b * DST_BATCH >= DST_BATCH ? p->full : p->rest, x, x);
    }
}

// 2D DST-I of the interior buffer, in place
static void dst2d(fft_solver *s)
{
    int parallel = s->nx * s->ny >= OMP_MIN_WORK;
    run_pass(&s->pass_rows, s->buf, parallel);
    run_pass(&s->pass_cols, s->buf, parallel);
}

fft_solver *fft_setup(int nx, int ny)
{
    int rows = ny - 2, cols = nx - 2; // interior nodes; fields are ny rows of nx
    fft_solver *s = (fft_solver *)malloc(sizeof(fft_solver));

    if (!s)
    {
        printf("** Error: insufficient memory **\n");
        exit(1);
    }
    s->nx = nx;
    s->ny = ny;
    s->buf = (double *)fftw_malloc((size_t)rows * cols * sizeof(double));
    if (!s->buf)
    {
        printf("** Error: fftw_malloc failed **\n");
        exit(1);
    }

    s->pass_rows = make_pass(s->buf, cols, rows, 1, cols); // rows of cols contiguous values (x)
    s->pass_cols = make_pass(s->buf, rows, cols, cols, 1); // columns, values cols apart (y)
    live_solvers++;
    return s;
}

void fft_cleanup(fft_solver *s)
{
    if (!s) return;
    free_pass(&s->pass_rows);
    free_pass(&s->pass_cols);
    fftw_free(s->buf);
    free(s);
    // Release FFTW's planner state once nothing uses it, so leak checkers see
    // nothing left
    if (--live_solvers == 0)
        fftw_cleanup();
}

void poisson_FFT(fft_solver *s, mtrx f, mtrx u, double dx, double dy)
{
    poisson_FFT_order(s, f, u, dx, dy, 2);
}

// f at interior node (i, j), or, for a wall node, the cubic extrapolation of
// the interior along the wall normal (f is not used on the walls, so the
// solve still depends on the interior of f only)
static double f_or_extrapolated(mtrx f, int i, int j)
{
    int ny = f.m, nx = f.n;
    if (j == 0) return 3.0 * MAt(f, i, 1) - 3.0 * MAt(f, i, 2) + MAt(f, i, 3);
    if (j == nx - 1) return 3.0 * MAt(f, i, nx - 2) - 3.0 * MAt(f, i, nx - 3) + MAt(f, i, nx - 4);
    if (i == 0) return 3.0 * MAt(f, 1, j) - 3.0 * MAt(f, 2, j) + MAt(f, 3, j);
    if (i == ny - 1) return 3.0 * MAt(f, ny - 2, j) - 3.0 * MAt(f, ny - 3, j) + MAt(f, ny - 4, j);
    return MAt(f, i, j);
}

void poisson_FFT_order(fft_solver *s, mtrx f, mtrx u, double dx, double dy, int order)
{
    int i;
    int nx = s->nx, ny = s->ny;
    double *fft_buf = s->buf;
    int rows = ny - 2, cols = nx - 2; // interior: row index i is y, column j is x
    // Order 4: the compact 9-point operator
    //   [dxx + dyy + (dx^2 + dy^2)/12 dxx dyy] u = [1 + (dx^2 dxx + dy^2 dyy)/12] f
    // with dxx, dyy the second differences. Its right-hand side is
    // f + (f_E + f_W + f_N + f_S - 4 f)/12, for any dx, dy.
    double cross = order == 4 ? (dx * dx + dy * dy) / 12.0 : 0.0;

    // Copy the interior right-hand side into the work buffer
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (nx * ny >= OMP_MIN_WORK)
#endif
    for (i = 0; i < rows; i++)
    {
        int j;
        for (j = 0; j < cols; j++)
        {
            double c = MAt(f, i + 1, j + 1);
            if (order == 4)
                c += (f_or_extrapolated(f, i + 1, j) + f_or_extrapolated(f, i + 1, j + 2) +
                      f_or_extrapolated(f, i, j + 1) + f_or_extrapolated(f, i + 2, j + 1) - 4.0 * c) /
                     12.0;
            fft_buf[i * cols + j] = c;
        }
    }

    // Forward DST-I
    dst2d(s);

    // Divide by eigenvalues of the 2D Laplacian under DST-I:
    //   λ_ij = (2*cos(π*(i+1)/(rows+1)) - 2) / dy²
    //         + (2*cos(π*(j+1)/(cols+1)) - 2) / dx²
    double inv_norm = 1.0 / (4.0 * (double)(rows + 1) * (double)(cols + 1));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (nx * ny >= OMP_MIN_WORK)
#endif
    for (i = 0; i < rows; i++)
    {
        double lambda_i = (2.0 * cos(PI * (i + 1) / (double)(rows + 1)) - 2.0) / (dy * dy);
        int j;
        for (j = 0; j < cols; j++)
        {
            double lambda_j = (2.0 * cos(PI * (j + 1) / (double)(cols + 1)) - 2.0) / (dx * dx);
            fft_buf[i * cols + j] /= (lambda_i + lambda_j + cross * lambda_i * lambda_j);
        }
    }

    // Inverse DST-I (same transform; normalise by 1/(2(rows+1)) * 1/(2(cols+1)))
    dst2d(s);

    // Write the normalised interior into u; u = 0 on the wall nodes
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (nx * ny >= OMP_MIN_WORK)
#endif
    for (i = 0; i < ny; i++)
    {
        int j;
        for (j = 0; j < nx; j++)
            MAt(u, i, j) = (i == 0 || i == ny - 1 || j == 0 || j == nx - 1)
                               ? 0.0
                               : fft_buf[(i - 1) * cols + (j - 1)] * inv_norm;
    }
}

// ---------------------------------------------------------------------------
// Periodic Poisson solver
// ---------------------------------------------------------------------------
// Solves (DX2 + DY2) u = f on a doubly periodic grid, with DX2, DY2 the
// solver's own second-derivative operators, so the Poisson operator has the
// same order as the transport derivatives. Both are circulant, so a 2D real
// FFT diagonalises them: mode (ky, kx) has eigenvalue lx[kx] + ly[ky], the
// symbol of each operator's stencil. The mean (mode 0, eigenvalue 0) is not
// determined by the equation; u is returned with zero mean, and the mean of f
// is ignored.
// ---------------------------------------------------------------------------

struct periodic_solver
{
    int nx, ny, kx;     // kx = nx/2 + 1 complex values per row of the spectrum
    double *in;         // ny x nx real field
    fftw_complex *spec; // ny x kx spectrum
    double *lx, *ly;    // eigenvalues of DX2 by kx and of DY2 by ky
    fftw_plan forward, backward;
};

// Symbol of a circulant operator at wavenumber k: sum over the entries of its
// first row of value * cos(2 pi k offset / n), offset counted in units of
// `step` columns. The stencils are symmetric, so the sine part vanishes.
static double symbol(const smtrx *A, int k, int n, int step)
{
    int e;
    double s = 0.0;
    for (e = A->row_ptr[0]; e < A->row_ptr[1]; e++)
    {
        int offset = A->col_idx[e] / step; // whole grid points along the axis
        s += A->values[e] * cos(2.0 * PI * (double)k * (double)offset / (double)n);
    }
    return s;
}

void periodic_eigenvalues(int nx, int ny, const smtrx *DX2, const smtrx *DY2, double *lx, double *ly)
{
    int k;
    // DX2 = I (x) dxx: its first row is the x stencil, in columns 0..nx-1.
    // DY2 = dyy (x) I: its first row is the y stencil, in columns m*nx.
    for (k = 0; k < nx / 2 + 1; k++)
        lx[k] = symbol(DX2, k, nx, 1);
    for (k = 0; k < ny; k++)
        ly[k] = symbol(DY2, k, ny, nx);
}

periodic_solver *periodic_setup(int nx, int ny, const smtrx *DX2, const smtrx *DY2)
{
    periodic_solver *s = (periodic_solver *)malloc(sizeof(periodic_solver));

    if (!s)
    {
        printf("** Error: insufficient memory **\n");
        exit(1);
    }
    s->nx = nx;
    s->ny = ny;
    s->kx = nx / 2 + 1;
    s->in = (double *)fftw_malloc((size_t)nx * ny * sizeof(double));
    s->spec = (fftw_complex *)fftw_malloc((size_t)s->kx * ny * sizeof(fftw_complex));
    s->lx = (double *)malloc((size_t)s->kx * sizeof(double));
    s->ly = (double *)malloc((size_t)ny * sizeof(double));
    if (!s->in || !s->spec || !s->lx || !s->ly)
    {
        printf("** Error: insufficient memory **\n");
        exit(1);
    }
    periodic_eigenvalues(nx, ny, DX2, DY2, s->lx, s->ly);

    s->forward = fftw_plan_dft_r2c_2d(ny, nx, s->in, s->spec, FFTW_ESTIMATE);
    s->backward = fftw_plan_dft_c2r_2d(ny, nx, s->spec, s->in, FFTW_ESTIMATE);
    if (!s->forward || !s->backward)
    {
        printf("** Error: FFTW could not plan the periodic transform **\n");
        exit(1);
    }
    live_solvers++;
    return s;
}

void periodic_cleanup(periodic_solver *s)
{
    if (!s) return;
    fftw_destroy_plan(s->forward);
    fftw_destroy_plan(s->backward);
    fftw_free(s->in);
    fftw_free(s->spec);
    free(s->lx);
    free(s->ly);
    free(s);
    if (--live_solvers == 0)
        fftw_cleanup();
}

void poisson_periodic(periodic_solver *s, mtrx f, mtrx u)
{
    int i, n = s->nx * s->ny;
    double inv_n = 1.0 / (double)n;

    for (i = 0; i < n; i++)
        s->in[i] = f.M[i];
    fftw_execute(s->forward);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n >= OMP_MIN_WORK)
#endif
    for (i = 0; i < s->ny; i++)
    {
        int j;
        for (j = 0; j < s->kx; j++)
        {
            double lambda = s->lx[j] + s->ly[i];
            double scale = (i == 0 && j == 0) ? 0.0 : inv_n / lambda;
            s->spec[i * s->kx + j][0] *= scale;
            s->spec[i * s->kx + j][1] *= scale;
        }
    }
    fftw_execute(s->backward);
    for (i = 0; i < n; i++)
        u.M[i] = s->in[i];
}
