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
#pragma omp parallel for schedule(static)
        for (i = 1; i < nx - 1; i++)
        {
            int j0 = (i & 1) ? 1 : 2;
            int j;
            for (j = j0; j < ny - 1; j += 2)
                MAt(u, i, j) = (dy2 * (MAt(u, i+1, j) + MAt(u, i-1, j))
                            + dx2 * (MAt(u, i, j+1) + MAt(u, i, j-1))
                            - dx2 * dy2 * MAt(f, i, j)) / denom;
        }
#pragma omp parallel for schedule(static)
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
#pragma omp parallel for schedule(static)
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
#pragma omp parallel for schedule(static)
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
// Solves ∇²u = f on [0,Lx] x [0,Ly] with homogeneous Dirichlet BCs (u=0
// on all boundaries) using the 2D Discrete Sine Transform (DST-I).
//
// The DST-I diagonalises the second-derivative finite-difference operator,
// so the solution is exact (to floating-point precision) in a single pass:
//   1. Forward DST-I of the interior RHS
//   2. Divide each mode by its eigenvalue
//   3. Inverse DST-I (= forward DST-I / (2*(nx+1)*(ny+1)))
//
// FFTW's RODFT00 plan is the DST-I.
// ---------------------------------------------------------------------------

#include <fftw3.h>

static fftw_plan plan_fwd;
static fftw_plan plan_inv;
static double   *fft_buf;  // shared work buffer, size (nx)*(ny)
static int       fft_nx;
static int       fft_ny;

void fft_setup(int nx, int ny)
{
    fft_nx  = nx;
    fft_ny  = ny;
    fft_buf = (double *)fftw_malloc((size_t)nx * ny * sizeof(double));
    if (!fft_buf) { printf("** Error: fftw_malloc failed **\n"); exit(1); }

    // FFTW_RODFT00 = DST-I in both dimensions
    // The transform operates on an nx x ny array stored row-major.
    plan_fwd = fftw_plan_r2r_2d(nx, ny,
                                 fft_buf, fft_buf,
                                 FFTW_RODFT00, FFTW_RODFT00,
                                 FFTW_MEASURE);
    plan_inv = fftw_plan_r2r_2d(nx, ny,
                                 fft_buf, fft_buf,
                                 FFTW_RODFT00, FFTW_RODFT00,
                                 FFTW_MEASURE);
}

void fft_cleanup(void)
{
    fftw_destroy_plan(plan_fwd);
    fftw_destroy_plan(plan_inv);
    fftw_free(fft_buf);
}

void poisson_FFT(mtrx f, mtrx u, double dx, double dy)
{
    int i;
    int nx = fft_nx;
    int ny = fft_ny;

    // Copy interior RHS into the work buffer (boundaries stay 0 by Dirichlet)
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (i = 0; i < nx; i++)
    {
        int j;
        for (j = 0; j < ny; j++)
            fft_buf[i * ny + j] = MAt(f, i, j);
    }

    // Forward DST-I
    fftw_execute(plan_fwd);

    // Divide by eigenvalues of the 2D Laplacian under DST-I:
    //   λ_ij = (2*cos(π*(i+1)/(nx+1)) - 2) / dx²
    //         + (2*cos(π*(j+1)/(ny+1)) - 2) / dy²
    double inv_norm = 1.0 / (4.0 * (double)(nx + 1) * (double)(ny + 1));
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (i = 0; i < nx; i++)
    {
        double lambda_i = (2.0 * cos(PI * (i + 1) / (double)(nx + 1)) - 2.0)
                          / (dx * dx);
        int j;
        for (j = 0; j < ny; j++)
        {
            double lambda_j = (2.0 * cos(PI * (j + 1) / (double)(ny + 1)) - 2.0)
                              / (dy * dy);
            fft_buf[i * ny + j] /= (lambda_i + lambda_j);
        }
    }

    // Inverse DST-I (same transform; normalise by 1/(2(nx+1)) * 1/(2(ny+1)))
    fftw_execute(plan_inv);

    // Write normalised result into u
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (i = 0; i < nx; i++)
    {
        int j;
        for (j = 0; j < ny; j++)
            MAt(u, i, j) = fft_buf[i * ny + j] * inv_norm;
    }
}