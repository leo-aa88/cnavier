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
void poisson(mtrx f, mtrx u, mtrx u0, double dx, double dy, int itmax, double tol)
{
    int i, k;
    int nx = f.m, ny = f.n;
    double e;
    double dx2 = dx * dx, dy2 = dy * dy;
    double denom = 2.0 * (dx2 + dy2);

    zerosm(u);

    for (k = 0; k < itmax; k++)
    {
        mtrxcpy(u0, u);
#ifdef _OPENMP
        /* Red–black ordering: two parallel phases (sequential GS is not safe to omp parallel for). */
#pragma omp parallel for schedule(static) if (nx * ny >= OMP_MIN_WORK)
        for (i = 1; i < nx - 1; i++)
        {
            int j0 = (i & 1) ? 1 : 2;
            int j;
            for (j = j0; j < ny - 1; j += 2)
                MAt(u, i, j) = (dy2 * (MAt(u, i+1, j) + MAt(u, i-1, j))
                            + dx2 * (MAt(u, i, j+1) + MAt(u, i, j-1))
                            - dx2 * dy2 * MAt(f, i, j)) / denom;
        }
#pragma omp parallel for schedule(static) if (nx * ny >= OMP_MIN_WORK)
        for (i = 1; i < nx - 1; i++)
        {
            int j0 = (i & 1) ? 2 : 1;
            int j;
            for (j = j0; j < ny - 1; j += 2)
                MAt(u, i, j) = (dy2 * (MAt(u, i+1, j) + MAt(u, i-1, j))
                            + dx2 * (MAt(u, i, j+1) + MAt(u, i, j-1))
                            - dx2 * dy2 * MAt(f, i, j)) / denom;
        }
#else
        for (i = 1; i < nx - 1; i++)
            for (int j = 1; j < ny - 1; j++)
                MAt(u, i, j) = (dy2 * (MAt(u, i+1, j) + MAt(u, i-1, j))
                            + dx2 * (MAt(u, i, j+1) + MAt(u, i, j-1))
                            - dx2 * dy2 * MAt(f, i, j)) / denom;
#endif
        e = error(u, u0);
        if (e < tol)
        {
            printf("Poisson solved in %d iterations - RSS error: %E\n", k, e);
            return;
        }
    }
    printf("Error: max iterations reached for Poisson solver.\n");
    exit(1);
}

// SOR Poisson solver.
// u and u0 are pre-allocated by the caller; u holds the result on return.
void poisson_SOR(mtrx f, mtrx u, mtrx u0, double dx, double dy, int itmax, double tol, double beta)
{
    int i, k;
    int nx = f.m, ny = f.n;
    double e;
    double dx2 = dx * dx, dy2 = dy * dy;
    double denom = 2.0 * (dx2 + dy2);

    zerosm(u);

    for (k = 0; k < itmax; k++)
    {
        mtrxcpy(u0, u);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (nx * ny >= OMP_MIN_WORK)
        for (i = 1; i < nx - 1; i++)
        {
            int j0 = (i & 1) ? 1 : 2;
            int j;
            for (j = j0; j < ny - 1; j += 2)
                MAt(u, i, j) = beta  * (dy2 * (MAt(u, i+1, j) + MAt(u, i-1, j))
                                    + dx2 * (MAt(u, i, j+1) + MAt(u, i, j-1))
                                    - dx2 * dy2 * MAt(f, i, j)) / denom
                           + (1.0 - beta) * MAt(u0, i, j);
        }
#pragma omp parallel for schedule(static) if (nx * ny >= OMP_MIN_WORK)
        for (i = 1; i < nx - 1; i++)
        {
            int j0 = (i & 1) ? 2 : 1;
            int j;
            for (j = j0; j < ny - 1; j += 2)
                MAt(u, i, j) = beta  * (dy2 * (MAt(u, i+1, j) + MAt(u, i-1, j))
                                    + dx2 * (MAt(u, i, j+1) + MAt(u, i, j-1))
                                    - dx2 * dy2 * MAt(f, i, j)) / denom
                           + (1.0 - beta) * MAt(u0, i, j);
        }
#else
        for (i = 1; i < nx - 1; i++)
            for (int j = 1; j < ny - 1; j++)
                MAt(u, i, j) = beta  * (dy2 * (MAt(u, i+1, j) + MAt(u, i-1, j))
                                    + dx2 * (MAt(u, i, j+1) + MAt(u, i, j-1))
                                    - dx2 * dy2 * MAt(f, i, j)) / denom
                           + (1.0 - beta) * MAt(u0, i, j);
#endif
        e = error(u, u0);
        if (e < tol)
        {
            printf("Poisson SOR solved in %d iterations - RSS error: %E\n", k, e);
            return;
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
    int       count;
    ptrdiff_t step; // offset between batches, in doubles
} dst_pass;

static dst_pass pass_rows, pass_cols;
static double  *fft_buf;  // interior work buffer, (nx-2) x (ny-2), row-major
static int      fft_nx;
static int      fft_ny;

// count transforms of length len, elements stride apart, transforms dist apart
static dst_pass make_pass(int len, int count, int stride, int dist)
{
    dst_pass p;
    fftw_r2r_kind kind = FFTW_RODFT00;

    p.count = count;
    p.step  = (ptrdiff_t)DST_BATCH * dist;
    p.full  = NULL;
    p.rest  = NULL;
    if (count >= DST_BATCH)
        p.full = fftw_plan_many_r2r(1, &len, DST_BATCH, fft_buf, NULL, stride, dist,
                                    fft_buf, NULL, stride, dist, &kind, FFTW_ESTIMATE);
    if (count % DST_BATCH)
        p.rest = fftw_plan_many_r2r(1, &len, count % DST_BATCH, fft_buf, NULL, stride, dist,
                                    fft_buf, NULL, stride, dist, &kind, FFTW_ESTIMATE);
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

static void run_pass(const dst_pass *p, int parallel)
{
    int b, batches = (p->count + DST_BATCH - 1) / DST_BATCH;
    (void)parallel;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (parallel)
#endif
    for (b = 0; b < batches; b++)
    {
        // Every batch but a short last one holds exactly DST_BATCH transforms
        double *x = fft_buf + b * p->step;
        fftw_execute_r2r(p->count - b * DST_BATCH >= DST_BATCH ? p->full : p->rest, x, x);
    }
}

// 2D DST-I of the interior buffer, in place
static void dst2d(void)
{
    int parallel = fft_nx * fft_ny >= OMP_MIN_WORK;
    run_pass(&pass_rows, parallel);
    run_pass(&pass_cols, parallel);
}

void fft_setup(int nx, int ny)
{
    int mx = nx - 2, my = ny - 2; // interior nodes

    fft_nx  = nx;
    fft_ny  = ny;
    fft_buf = (double *)fftw_malloc((size_t)mx * my * sizeof(double));
    if (!fft_buf) { printf("** Error: fftw_malloc failed **\n"); exit(1); }

    pass_rows = make_pass(my, mx, 1, my); // mx rows of my contiguous values
    pass_cols = make_pass(mx, my, my, 1); // my columns, values my apart
}

void fft_cleanup(void)
{
    free_pass(&pass_rows);
    free_pass(&pass_cols);
    fftw_free(fft_buf);
    fftw_cleanup(); // release FFTW's planner state, so leak checkers see nothing left
}

void poisson_FFT(mtrx f, mtrx u, double dx, double dy)
{
    int i;
    int nx = fft_nx, ny = fft_ny;
    int mx = nx - 2, my = ny - 2;

    // Copy the interior right-hand side into the work buffer
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (nx * ny >= OMP_MIN_WORK)
#endif
    for (i = 0; i < mx; i++)
    {
        int j;
        for (j = 0; j < my; j++)
            fft_buf[i * my + j] = MAt(f, i + 1, j + 1);
    }

    // Forward DST-I
    dst2d();

    // Divide by eigenvalues of the 2D Laplacian under DST-I:
    //   λ_ij = (2*cos(π*(i+1)/(mx+1)) - 2) / dx²
    //         + (2*cos(π*(j+1)/(my+1)) - 2) / dy²
    double inv_norm = 1.0 / (4.0 * (double)(mx + 1) * (double)(my + 1));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (nx * ny >= OMP_MIN_WORK)
#endif
    for (i = 0; i < mx; i++)
    {
        double lambda_i = (2.0 * cos(PI * (i + 1) / (double)(mx + 1)) - 2.0)
                          / (dx * dx);
        int j;
        for (j = 0; j < my; j++)
        {
            double lambda_j = (2.0 * cos(PI * (j + 1) / (double)(my + 1)) - 2.0)
                              / (dy * dy);
            fft_buf[i * my + j] /= (lambda_i + lambda_j);
        }
    }

    // Inverse DST-I (same transform; normalise by 1/(2(mx+1)) * 1/(2(my+1)))
    dst2d();

    // Write the normalised interior into u; u = 0 on the wall nodes
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (nx * ny >= OMP_MIN_WORK)
#endif
    for (i = 0; i < nx; i++)
    {
        int j;
        for (j = 0; j < ny; j++)
            MAt(u, i, j) = (i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
                         ? 0.0 : fft_buf[(i - 1) * my + (j - 1)] * inv_norm;
    }
}
