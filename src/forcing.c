#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "forcing.h"
#include "poisson.h"

#define FORCING_PI 3.14159265358979323846

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
                                     const periodic_symbols *sym)
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
        int ix = (jx + nx) % nx, iy = (jy + ny) % ny;
        double A = sym->d1x_re[ix] * sym->d1x_re[ix] + sym->d1x_im[ix] * sym->d1x_im[ix] +
                   sym->d1y_re[iy] * sym->d1y_re[iy] + sym->d1y_im[iy] * sym->d1y_im[iy];
        double Q = -(sym->d2x[ix] + sym->d2y[iy]);
        rf->amp[m] = 2.0 * Q * sqrt(f->random_rate * dt / (rf->modes * A));
    }
    rf->nx = nx;
    rf->ny = ny;
    rf->cx = (double *)malloc((size_t)rf->modes * nx * sizeof(double));
    rf->sx = (double *)malloc((size_t)rf->modes * nx * sizeof(double));
    rf->cy = (double *)malloc((size_t)rf->modes * ny * sizeof(double));
    rf->sy = (double *)malloc((size_t)rf->modes * ny * sizeof(double));
    if (!rf->cx || !rf->sx || !rf->cy || !rf->sy)
    {
        printf("** Error: insufficient memory **\n");
        exit(1);
    }
    for (m = 0; m < rf->modes; m++)
        for (n = 0; n < nx; n++)
        {
            rf->cx[(size_t)m * nx + n] = rf->amp[m] * cos(rf->kx[m] * (n * dx));
            rf->sx[(size_t)m * nx + n] = rf->amp[m] * sin(rf->kx[m] * (n * dx));
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
    free(rf->cx);
    free(rf->sx);
    free(rf->cy);
    free(rf->sy);
    free(rf);
}

void random_forcing_draw(random_forcing *rf)
{
    int m;
    for (m = 0; m < rf->modes; m++)
        rf->phase[m] = 2.0 * FORCING_PI * uniform(&rf->state);
}

void random_forcing_add(random_forcing *rf, mtrx w, double dy)
{
    int i, m, nx = rf->nx, ny = rf->ny;

    for (m = 0; m < rf->modes; m++)
        for (i = 0; i < ny; i++)
        {
            rf->cy[(size_t)m * ny + i] = cos(rf->ky[m] * (i * dy) + rf->phase[m]);
            rf->sy[(size_t)m * ny + i] = sin(rf->ky[m] * (i * dy) + rf->phase[m]);
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (w.m * w.n >= OMP_MIN_WORK)
#endif
    for (i = 0; i < ny; i++)
        for (int j = 0; j < nx; j++)
        {
            double s = 0.0;
            for (int k = 0; k < rf->modes; k++)
                s += rf->cx[(size_t)k * nx + j] * rf->cy[(size_t)k * ny + i] -
                     rf->sx[(size_t)k * nx + j] * rf->sy[(size_t)k * ny + i];
            MAt(w, i, j) += s;
        }
}

// A phase in [0, 2 pi) for wavevector (m, n) and the seed, the same on every
// grid that holds the mode
static double mode_phase(unsigned long long seed, int m, int n)
{
    unsigned long long state = seed ^ ((unsigned long long)(unsigned)m * 0x9E3779B97F4A7C15ULL) ^
                               ((unsigned long long)(unsigned)n * 0xC2B2AE3D27D4EB4FULL);
    uniform(&state);
    return 2.0 * FORCING_PI * uniform(&state);
}

void random_initial_field(mtrx w, mtrx u, mtrx v, double dx, double dy, double k0, double energy,
                          unsigned long long seed)
{
    int i, j, nx = w.n, ny = w.m, h = nx / 2 + 1, n = nx * ny;
    double Lx = nx * dx, Ly = ny * dy, dk0 = 2.0 * FORCING_PI / (Lx > Ly ? Lx : Ly), e = 0.0, scale;
    double *real = fftw_alloc_real((size_t)n);
    fftw_complex *psi = fftw_alloc_complex((size_t)ny * h), *tmp = fftw_alloc_complex((size_t)ny * h);
    fftw_plan inv;
    mtrx *out[3] = {&w, &u, &v};

    if (!real || !psi || !tmp)
    {
        printf("** Error: insufficient memory **\n");
        exit(1);
    }
    fftw_plans_hold();
    inv = fftw_plan_dft_c2r_2d(ny, nx, tmp, real, FFTW_ESTIMATE);

    // psi^ of mode (m, n) has the phase mode_phase(m, n), the same on every
    // grid, and the amplitude of the spectrum: a shell of radius k holds ~k
    // modes, each of energy k^2 |psi^|^2 / 2, so |psi^| ~ k^(1/2) exp(-(k/k0)^2)
    // up to the constant fixed below. The column m = 0 of the half spectrum
    // holds (0, n) and its conjugate (0, -n). No Nyquist modes.
    for (i = 0; i < ny; i++)
        for (j = 0; j < h; j++)
        {
            int jy = i <= ny / 2 ? i : i - ny;
            double kx = 2.0 * FORCING_PI * j / Lx, ky = 2.0 * FORCING_PI * jy / Ly;
            double k = sqrt(kx * kx + ky * ky) / dk0, *c = psi[i * h + j];
            int nyquist = (nx % 2 == 0 && j == nx / 2) || (ny % 2 == 0 && i == ny / 2);
            double a = (k == 0.0 || nyquist) ? 0.0 : sqrt(k) * exp(-(k / k0) * (k / k0));
            double ph = (j == 0 && jy < 0) ? -mode_phase(seed, 0, -jy) : mode_phase(seed, j, jy);
            c[0] = a * cos(ph);
            c[1] = a * sin(ph);
        }

    // w = -lap psi, u = dpsi/dy, v = -dpsi/dx: multiply by k^2, i ky, -i kx
    for (int f = 0; f < 3; f++)
    {
        for (i = 0; i < ny; i++)
            for (j = 0; j < h; j++)
            {
                int jy = i <= ny / 2 ? i : i - ny;
                double kx = 2.0 * FORCING_PI * j / Lx, ky = 2.0 * FORCING_PI * jy / Ly;
                double *c = psi[i * h + j], *d = tmp[i * h + j];
                if (f == 0)
                {
                    d[0] = (kx * kx + ky * ky) * c[0];
                    d[1] = (kx * kx + ky * ky) * c[1];
                }
                else
                {
                    double s = f == 1 ? ky : -kx; // multiply by i s
                    d[0] = -s * c[1];
                    d[1] = s * c[0];
                }
            }
        fftw_execute(inv);
        for (i = 0; i < n; i++)
            out[f]->M[i] = real[i];
    }
    for (i = 0; i < n; i++)
        e += 0.5 * (u.M[i] * u.M[i] + v.M[i] * v.M[i]);
    scale = e > 0.0 ? sqrt(energy / (e / n)) : 0.0;
    for (int f = 0; f < 3; f++)
        for (i = 0; i < n; i++)
            out[f]->M[i] *= scale;

    fftw_destroy_plan(inv);
    fftw_plans_release();
    fftw_free(real);
    fftw_free(psi);
    fftw_free(tmp);
}
