#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "forcing.h"

#define FORCING_PI 3.14159265358979323846

// Symbol of the circulant operator A at wavenumber index k along an axis of n
// points: the entries of its first row, offset counted in units of `step`
// columns. Returns the real and imaginary parts.
static void symbol(const smtrx *A, int k, int n, int step, double *re, double *im)
{
    int e;
    *re = 0.0;
    *im = 0.0;
    for (e = A->row_ptr[0]; e < A->row_ptr[1]; e++)
    {
        int offset = A->col_idx[e] / step;
        double th = 2.0 * FORCING_PI * (double)k * (double)offset / (double)n;
        *re += A->values[e] * cos(th);
        *im += A->values[e] * sin(th);
    }
}

// splitmix64: a small, well-mixed generator; uniform in [0, 1)
static double uniform(unsigned long long *state)
{
    unsigned long long z = (*state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    z ^= z >> 31;
    return (double)(z >> 11) * (1.0 / 9007199254740992.0);
}

double *kolmogorov_rows(const forcing_config *f, int ny, double dy, double Ly)
{
    int i;
    double k = 2.0 * FORCING_PI * f->kolmogorov_n / Ly, *rows;

    if (f->kolmogorov_amp == 0.0) return NULL;
    if (!(rows = (double *)malloc(ny * sizeof(double))))
    {
        printf("** Error: insufficient memory **\n");
        exit(1);
    }
    for (i = 0; i < ny; i++)
        rows[i] = -f->kolmogorov_amp * k * cos(k * (i * dy));
    return rows;
}

random_forcing *random_forcing_setup(const forcing_config *f, int nx, int ny, double dx, double dy, double dt,
                                     const smtrx *DX, const smtrx *DY, const smtrx *DX2, const smtrx *DY2)
{
    int pass, m, n, count = 0;
    double Lx = nx * dx, Ly = ny * dy, dk0 = 2.0 * FORCING_PI / (Lx > Ly ? Lx : Ly);
    random_forcing *rf;

    if (f->random_rate <= 0.0) return NULL;
    rf = (random_forcing *)calloc(1, sizeof(random_forcing));
    if (!rf)
    {
        printf("** Error: insufficient memory **\n");
        exit(1);
    }

    // Two passes: count the modes in the shell, then fill them in. Half plane
    // (n > 0, or n = 0 and m > 0) so that no mode is the conjugate of another,
    // and |m| < nx/2, |n| < ny/2 so that none is its own conjugate.
    for (pass = 0; pass < 2; pass++)
    {
        count = 0;
        for (n = 0; 2 * n < ny; n++)
            for (m = -((nx - 1) / 2); 2 * m < nx; m++)
            {
                double kx = 2.0 * FORCING_PI * m / Lx, ky = 2.0 * FORCING_PI * n / Ly;
                double k = sqrt(kx * kx + ky * ky) / dk0;
                if ((n == 0 && m <= 0) || fabs(k - f->random_kf) > f->random_dk) continue;
                if (pass == 1)
                {
                    rf->kx[count] = kx;
                    rf->ky[count] = ky;
                }
                count++;
            }
        if (pass == 0)
        {
            if (count == 0)
            {
                printf("** Error: no wavevector of the grid lies in the random-forcing shell **\n");
                exit(1);
            }
            rf->modes = count;
            rf->kx = (double *)calloc(count, sizeof(double));
            rf->ky = (double *)calloc(count, sizeof(double));
            rf->amp = (double *)malloc(count * sizeof(double));
            rf->phase = (double *)malloc(count * sizeof(double));
            if (!rf->kx || !rf->ky || !rf->amp || !rf->phase)
            {
                printf("** Error: insufficient memory **\n");
                exit(1);
            }
        }
    }

    // A mode w = a cos(k.x + phase) has psi^ = w^ / Q and velocity amplitude
    // a sqrt(A) / Q, with A = |symbol of DX|^2 + |symbol of DY|^2 and
    // Q = -(symbol of DX2 + DY2), so its discrete energy is A a^2 / (4 Q^2).
    // Distinct modes are orthogonal on the grid, so with each mode carrying
    // eps dt / modes, a kick carries eps dt.
    for (m = 0; m < rf->modes; m++)
    {
        int jx = (int)lround(rf->kx[m] * Lx / (2.0 * FORCING_PI)), jy = (int)lround(rf->ky[m] * Ly / (2.0 * FORCING_PI));
        double xr, xi, yr, yi, x2r, x2i, y2r, y2i, A, Q;
        symbol(DX, (jx + nx) % nx, nx, 1, &xr, &xi);
        symbol(DY, (jy + ny) % ny, ny, nx, &yr, &yi);
        symbol(DX2, (jx + nx) % nx, nx, 1, &x2r, &x2i);
        symbol(DY2, (jy + ny) % ny, ny, nx, &y2r, &y2i);
        A = xr * xr + xi * xi + yr * yr + yi * yi;
        Q = -(x2r + y2r);
        rf->amp[m] = 2.0 * Q * sqrt(f->random_rate * dt / (rf->modes * A));
    }
    rf->state = f->random_seed;
    return rf;
}

void random_forcing_free(random_forcing *rf)
{
    if (!rf) return;
    free(rf->kx);
    free(rf->ky);
    free(rf->amp);
    free(rf->phase);
    free(rf);
}

void random_forcing_draw(random_forcing *rf)
{
    int m;
    for (m = 0; m < rf->modes; m++)
        rf->phase[m] = 2.0 * FORCING_PI * uniform(&rf->state);
}

void random_forcing_add(const random_forcing *rf, mtrx w, double dx, double dy)
{
    int i;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (w.m * w.n >= OMP_MIN_WORK)
#endif
    for (i = 0; i < w.m; i++)
        for (int j = 0; j < w.n; j++)
        {
            double s = 0.0;
            for (int m = 0; m < rf->modes; m++)
                s += rf->amp[m] * cos(rf->kx[m] * (j * dx) + rf->ky[m] * (i * dy) + rf->phase[m]);
            MAt(w, i, j) += s;
        }
}
