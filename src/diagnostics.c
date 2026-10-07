#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "diagnostics.h"
#include "poisson.h"
#include "utils.h"

// Velocity at node (i, j) for the integrals: the wall velocity on a wall node
// (left and right walls win at the corners, as in apply_wall_bc), else the
// field's value
static void velocity_at(const solver_config *cfg, mtrx u, mtrx v, int i, int j, double *uu, double *vv)
{
    int nx = cfg->nx, ny = cfg->ny, wall = -1;

    if (!cfg->periodic)
        wall = j == 0 ? 0 : j == nx - 1 ? 1
                        : i == 0        ? 2
                        : i == ny - 1   ? 3
                                        : -1;
    *uu = wall < 0 ? MAt(u, i, j) : cfg->bc.u[wall];
    *vv = wall < 0 ? MAt(v, i, j) : cfg->bc.v[wall];
}

flow_integrals compute_integrals(const solver_config *cfg, mtrx u, mtrx v, mtrx w, double *wx, double *wy)
{
    int i, j, nx = cfg->nx, ny = cfg->ny;
    double sE = 0.0, sZ = 0.0, sP = 0.0, norm;
    flow_integrals r;

    spmv(*cfg->DX, w.M, wx);
    spmv(*cfg->DY, w.M, wy);
    for (i = 0; i < ny; i++)
        for (j = 0; j < nx; j++)
        {
            int k = i * nx + j;
            double uu, vv, weight = 1.0;
            if (!cfg->periodic)
                weight = ((i == 0 || i == ny - 1) ? 0.5 : 1.0) * ((j == 0 || j == nx - 1) ? 0.5 : 1.0);
            velocity_at(cfg, u, v, i, j, &uu, &vv);
            sE += weight * (uu * uu + vv * vv);
            sZ += weight * w.M[k] * w.M[k];
            sP += weight * (wx[k] * wx[k] + wy[k] * wy[k]);
        }
    norm = cfg->periodic ? (double)nx * ny : (double)(nx - 1) * (ny - 1);
    r.E = 0.5 * sE / norm;
    r.Z = 0.5 * sZ / norm;
    r.P = 0.5 * sP / norm;
    return r;
}

// ---------------------------------------------------------------------------
// Spectra
// ---------------------------------------------------------------------------

struct spectra
{
    solver_config cfg;
    int nx, ny, kx, bins;
    double dk;
    int *bin;        // shell of each of the ny*kx modes
    double *weight;  // 1 or 2: modes kx and -kx in the half spectrum
    double *lap;     // eigenvalue of DX2 + DY2 of each mode
    double *in, *nl; // real field, nonlinear term
    double *wx, *wy; // derivatives of w
    fftw_complex *uh, *vh, *wh, *nh;
    fftw_plan plan; // r2c of `in` into a spectrum
};

spectra *spectra_setup(const solver_config *cfg)
{
    int i, j, nx = cfg->nx, ny = cfg->ny;
    double Lx = nx * cfg->dx, Ly = ny * cfg->dy, kmax = 0.0;
    spectra *s = (spectra *)calloc(1, sizeof(spectra));
    double *lx, *ly;

    if (!s || !cfg->periodic)
    {
        printf("** Error: spectra need a periodic grid **\n");
        exit(1);
    }
    s->cfg = *cfg;
    s->nx = nx;
    s->ny = ny;
    s->kx = nx / 2 + 1;
    s->dk = 2.0 * PI / (Lx > Ly ? Lx : Ly);
    s->bin = (int *)malloc((size_t)s->kx * ny * sizeof(int));
    s->weight = (double *)malloc((size_t)s->kx * ny * sizeof(double));
    s->lap = (double *)malloc((size_t)s->kx * ny * sizeof(double));
    s->in = (double *)fftw_malloc((size_t)nx * ny * sizeof(double));
    s->nl = (double *)malloc((size_t)nx * ny * sizeof(double));
    s->wx = (double *)malloc((size_t)nx * ny * sizeof(double));
    s->wy = (double *)malloc((size_t)nx * ny * sizeof(double));
    s->uh = (fftw_complex *)fftw_malloc((size_t)s->kx * ny * sizeof(fftw_complex));
    s->vh = (fftw_complex *)fftw_malloc((size_t)s->kx * ny * sizeof(fftw_complex));
    s->wh = (fftw_complex *)fftw_malloc((size_t)s->kx * ny * sizeof(fftw_complex));
    s->nh = (fftw_complex *)fftw_malloc((size_t)s->kx * ny * sizeof(fftw_complex));
    lx = (double *)malloc((size_t)s->kx * sizeof(double));
    ly = (double *)malloc((size_t)ny * sizeof(double));
    if (!s->bin || !s->weight || !s->lap || !s->in || !s->nl || !s->wx || !s->wy || !s->uh || !s->vh || !s->wh ||
        !s->nh || !lx || !ly)
    {
        printf("** Error: insufficient memory **\n");
        exit(1);
    }
    periodic_eigenvalues(nx, ny, cfg->DX2, cfg->DY2, lx, ly);
    for (i = 0; i < ny; i++)
        for (j = 0; j < s->kx; j++)
        {
            int m = i <= ny / 2 ? i : i - ny; // signed wavenumber index in y
            double kxv = 2.0 * PI * j / Lx, kyv = 2.0 * PI * m / Ly, k = sqrt(kxv * kxv + kyv * kyv);
            int idx = i * s->kx + j;
            s->bin[idx] = (int)floor(k / s->dk + 0.5);
            // Modes j and nx - j are conjugate; the half spectrum holds one of
            // each pair, except j = 0 and, for even nx, j = nx/2
            s->weight[idx] = (j == 0 || (nx % 2 == 0 && j == nx / 2)) ? 1.0 : 2.0;
            s->lap[idx] = lx[j] + ly[i];
            if (k > kmax) kmax = k;
        }
    s->bins = (int)floor(kmax / s->dk + 0.5) + 1;
    free(lx);
    free(ly);
    s->plan = fftw_plan_dft_r2c_2d(ny, nx, s->in, s->uh, FFTW_ESTIMATE);
    if (!s->plan)
    {
        printf("** Error: FFTW could not plan the spectra **\n");
        exit(1);
    }
    fftw_plans_hold();
    return s;
}

int spectra_bins(const spectra *s)
{
    return s->bins;
}

double spectra_dk(const spectra *s)
{
    return s->dk;
}

// Spectrum of the nx*ny field x
static void transform(spectra *s, const double *x, fftw_complex *out)
{
    int k, n = s->nx * s->ny;
    for (k = 0; k < n; k++)
        s->in[k] = x[k];
    fftw_execute_dft_r2c(s->plan, s->in, out);
}

void spectra_compute(spectra *s, mtrx u, mtrx v, mtrx w, double *E, double *Z, double *PE, double *PZ)
{
    int b, k, n = s->nx * s->ny, modes = s->kx * s->ny;
    double inv = 1.0 / ((double)n * n); // Parseval: mean(f^2) = sum |f^|^2 / n^2
    double *TE = (double *)calloc(s->bins, sizeof(double)), *TZ = (double *)calloc(s->bins, sizeof(double));

    if (!TE || !TZ)
    {
        printf("** Error: insufficient memory **\n");
        exit(1);
    }
    if (E)
        for (b = 0; b < s->bins; b++)
            E[b] = 0.0;
    if (Z)
        for (b = 0; b < s->bins; b++)
            Z[b] = 0.0;

    // The solver's nonlinear term N = -(u DX w + v DY w)
    spmv(*s->cfg.DX, w.M, s->wx);
    spmv(*s->cfg.DY, w.M, s->wy);
    for (k = 0; k < n; k++)
        s->nl[k] = -(u.M[k] * s->wx[k] + v.M[k] * s->wy[k]);

    transform(s, u.M, s->uh);
    transform(s, v.M, s->vh);
    transform(s, w.M, s->wh);
    transform(s, s->nl, s->nh);
    for (k = 0; k < modes; k++)
    {
        double c = s->weight[k] * inv;
        double wr = s->wh[k][0], wi = s->wh[k][1];
        // psi^ = -w^ / lap; the mean mode carries no energy
        double pr = k == 0 ? 0.0 : -wr / s->lap[k], pi = k == 0 ? 0.0 : -wi / s->lap[k];
        int bb = s->bin[k];
        if (E)
            E[bb] += 0.5 * c * (s->uh[k][0] * s->uh[k][0] + s->uh[k][1] * s->uh[k][1] + s->vh[k][0] * s->vh[k][0] + s->vh[k][1] * s->vh[k][1]);
        if (Z) Z[bb] += 0.5 * c * (wr * wr + wi * wi);
        TE[bb] += c * (pr * s->nh[k][0] + pi * s->nh[k][1]);
        TZ[bb] += c * (wr * s->nh[k][0] + wi * s->nh[k][1]);
    }
    // Fluxes: what the nonlinear term removes from all bins up to k
    for (b = 0; b < s->bins; b++)
    {
        double prevE = b ? (PE ? PE[b - 1] : 0.0) : 0.0, prevZ = b ? (PZ ? PZ[b - 1] : 0.0) : 0.0;
        if (PE) PE[b] = prevE - TE[b];
        if (PZ) PZ[b] = prevZ - TZ[b];
    }
    free(TE);
    free(TZ);
}

void spectra_write(spectra *s, mtrx u, mtrx v, mtrx w, double t)
{
    int b, frame = output_frame("spectrum", ".csv");
    char name[96];
    FILE *f;
    double *E = (double *)malloc((size_t)s->bins * sizeof(double));
    double *Z = (double *)malloc((size_t)s->bins * sizeof(double));
    double *PE = (double *)malloc((size_t)s->bins * sizeof(double));
    double *PZ = (double *)malloc((size_t)s->bins * sizeof(double));

    if (!E || !Z || !PE || !PZ)
    {
        printf("** Error: insufficient memory **\n");
        exit(1);
    }
    spectra_compute(s, u, v, w, E, Z, PE, PZ);
    snprintf(name, sizeof(name), "./output/spectrum-1-%d.csv", frame);
    if (!(f = fopen(name, "w")))
    {
        printf("\nError while opening file\n");
        exit(1);
    }
    fprintf(f, "# t = %.17g\nk,E,Z,Pi_E,Pi_Z\n", t);
    for (b = 0; b < s->bins; b++)
        fprintf(f, "%.17g,%.17g,%.17g,%.17g,%.17g\n", b * s->dk, E[b], Z[b], PE[b], PZ[b]);
    fclose(f);
    free(E);
    free(Z);
    free(PE);
    free(PZ);
}

void spectra_free(spectra *s)
{
    if (!s) return;
    fftw_destroy_plan(s->plan);
    fftw_plans_release();
    free(s->bin);
    free(s->weight);
    free(s->lap);
    fftw_free(s->in);
    free(s->nl);
    free(s->wx);
    free(s->wy);
    fftw_free(s->uh);
    fftw_free(s->vh);
    fftw_free(s->wh);
    fftw_free(s->nh);
    free(s);
}
