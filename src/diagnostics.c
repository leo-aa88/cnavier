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

double kolmogorov_factor(const solver_config *cfg)
{
    int e, n = cfg->forcing.kolmogorov_n, ny = cfg->ny, nx = cfg->nx;
    double k, k1 = 0.0, Q = 0.0;

    if (!cfg->periodic || cfg->forcing.kolmogorov_amp == 0.0) return 1.0;
    k = 2.0 * PI * n / (ny * cfg->dy);
    // DY = dyy (x) I and DY2 likewise: their first rows hold the y stencils
    // in columns m*nx. DY e^{iky} = i k1 e^{iky}, DY2 e^{iky} = -Q e^{iky}.
    for (e = cfg->DY->row_ptr[0]; e < cfg->DY->row_ptr[1]; e++)
    {
        int m = cfg->DY->col_idx[e] / nx; // whole rows: the offset along y
        k1 += cfg->DY->values[e] * sin(2.0 * PI * n * (double)m / ny);
    }
    for (e = cfg->DY2->row_ptr[0]; e < cfg->DY2->row_ptr[1]; e++)
    {
        int m = cfg->DY2->col_idx[e] / nx;
        Q -= cfg->DY2->values[e] * cos(2.0 * PI * n * (double)m / ny);
    }
    return k * k1 / Q;
}

flow_integrals compute_integrals(const solver_config *cfg, mtrx u, mtrx v, mtrx w, double *wx, double *wy)
{
    int i, j, nx = cfg->nx, ny = cfg->ny;
    double sE = 0.0, sZ = 0.0, sP = 0.0, sI = 0.0, norm;
    double Ly = cfg->dy * (cfg->periodic ? ny : ny - 1), kk = 2.0 * PI * cfg->forcing.kolmogorov_n / Ly;
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
            sI += weight * uu * sin(kk * (i * cfg->dy));
        }
    norm = cfg->periodic ? (double)nx * ny : (double)(nx - 1) * (ny - 1);
    r.E = 0.5 * sE / norm;
    r.Z = 0.5 * sZ / norm;
    r.P = 0.5 * sP / norm;
    // Energy input: the work of the Kolmogorov force A sin(k y) on u, and the
    // rate the random kicks inject by construction
    r.I = cfg->forcing.kolmogorov_amp * sI / norm + cfg->forcing.random_rate;
    r.I_disc = kolmogorov_factor(cfg) * cfg->forcing.kolmogorov_amp * sI / norm + cfg->forcing.random_rate;
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
    double *ratio;   // A/Q of each mode: |symbol of DX|^2 + |symbol of DY|^2 over -lap
    double *in, *nl; // real field, nonlinear term
    double *wx, *wy; // derivatives of w
    fftw_complex *uh, *vh, *wh, *nh;
    fftw_plan plan; // r2c of `in` into a spectrum
};

// |symbol|^2 of a circulant first-derivative operator at wavenumber k: the
// entries of its first row, offset counted in units of `step` columns
static double symbol_sq(const smtrx *A, int k, int n, int step)
{
    int e;
    double re = 0.0, im = 0.0;
    for (e = A->row_ptr[0]; e < A->row_ptr[1]; e++)
    {
        int offset = A->col_idx[e] / step;
        double th = 2.0 * PI * (double)k * (double)offset / (double)n;
        re += A->values[e] * cos(th);
        im += A->values[e] * sin(th);
    }
    return re * re + im * im;
}

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
    s->ratio = (double *)malloc((size_t)s->kx * ny * sizeof(double));
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
    if (!s->bin || !s->weight || !s->lap || !s->ratio || !s->in || !s->nl || !s->wx || !s->wy || !s->uh || !s->vh || !s->wh ||
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
            s->ratio[idx] = idx == 0 ? 0.0
                                     : (symbol_sq(cfg->DXv ? cfg->DXv : cfg->DX, j, nx, 1) +
                                        symbol_sq(cfg->DYv ? cfg->DYv : cfg->DY, i, ny, nx)) /
                                           -s->lap[idx];
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

    // The solver's nonlinear term N = -(u DX w + v DY w), or its
    // skew-symmetric form (skew_correction())
    spmv(*s->cfg.DX, w.M, s->wx);
    spmv(*s->cfg.DY, w.M, s->wy);
    for (k = 0; k < n; k++)
        s->nl[k] = -(u.M[k] * s->wx[k] + v.M[k] * s->wy[k]);
    if (s->cfg.advection == 1) skew_correction(&s->cfg, u.M, v.M, w.M, s->wx, s->wy, s->nl, s->wx, s->wy);

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
        // The energy of a mode is 1/2 A |psi^|^2 with A from the first-
        // derivative operators that give u and v, and w^ = Q psi^ with Q from
        // the Laplacian, so the nonlinear term changes it at (A/Q) Re(psi^* N^)
        TE[bb] += c * s->ratio[k] * (pr * s->nh[k][0] + pi * s->nh[k][1]);
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

void spectra_dissipation(spectra *s, mtrx w, double *DE, double *DZ, double *FE, double *FZ)
{
    int b, k, n = s->nx * s->ny, modes = s->kx * s->ny;
    const forcing_config *fc = &s->cfg.forcing;
    double inv = 1.0 / ((double)n * n), nu = 1.0 / s->cfg.Re;

    for (b = 0; b < s->bins; b++)
    {
        if (DE) DE[b] = 0.0;
        if (DZ) DZ[b] = 0.0;
        if (FE) FE[b] = 0.0;
        if (FZ) FZ[b] = 0.0;
    }
    transform(s, w.M, s->wh);
    // The mean vorticity (zero in a physical periodic flow) carries no energy
    // and feels only the drag
    if (FZ) FZ[0] = fc->drag * (s->wh[0][0] * s->wh[0][0] + s->wh[0][1] * s->wh[0][1]) * inv;
    for (k = 1; k < modes; k++)
    {
        // A damping term of symbol -sigma changes the mode's energy
        // 1/2 A |psi^|^2 = 1/2 (A/Q^2) |w^|^2 at -sigma (A/Q^2) |w^|^2 and its
        // enstrophy at -sigma |w^|^2. Viscosity has sigma = nu Q, the
        // hyperviscosity nu_h Q^p, the drag alpha and the hypodrag alpha_h / Q.
        double c = s->weight[k] * inv, w2 = s->wh[k][0] * s->wh[k][0] + s->wh[k][1] * s->wh[k][1];
        double Q = -s->lap[k], small = nu * Q, large = fc->drag + fc->hypodrag / Q;
        if (fc->hyperviscosity > 0.0) small += fc->hyperviscosity * pow(Q, fc->hyper_order);
        if (DE) DE[s->bin[k]] += small * c * s->ratio[k] / Q * w2;
        if (DZ) DZ[s->bin[k]] += small * c * w2;
        if (FE) FE[s->bin[k]] += large * c * s->ratio[k] / Q * w2;
        if (FZ) FZ[s->bin[k]] += large * c * w2;
    }
}

void spectra_all(spectra *s, mtrx u, mtrx v, mtrx w, double *out)
{
    size_t B = (size_t)s->bins;
    spectra_compute(s, u, v, w, out, out + B, out + 2 * B, out + 3 * B);
    spectra_dissipation(s, w, out + 4 * B, out + 5 * B, out + 6 * B, out + 7 * B);
}

void spectra_tables(const spectra *s, const int **bin, const double **weight, const double **lap,
                    const double **ratio)
{
    *bin = s->bin;
    *weight = s->weight;
    *lap = s->lap;
    *ratio = s->ratio;
}

void spectra_write_frame(const spectra *s, const double *out, double t)
{
    int b, q, B = s->bins, frame = output_frame("spectrum", ".csv");
    char name[96];
    FILE *f;

    snprintf(name, sizeof(name), "./output/spectrum-1-%d.csv", frame);
    if (!(f = fopen(name, "w")))
    {
        printf("\nError while opening file\n");
        exit(1);
    }
    fprintf(f, "# t = %.17g\nk,E,Z,Pi_E,Pi_Z,D_E,D_Z,F_E,F_Z\n", t);
    for (b = 0; b < B; b++)
    {
        fprintf(f, "%.17g", b * s->dk);
        for (q = 0; q < SPECTRA_COLUMNS; q++)
            fprintf(f, ",%.17g", out[q * B + b]);
        fprintf(f, "\n");
    }
    fclose(f);
}

void spectra_write(spectra *s, mtrx u, mtrx v, mtrx w, double t)
{
    double *out = (double *)malloc((size_t)SPECTRA_COLUMNS * s->bins * sizeof(double));
    if (!out)
    {
        printf("** Error: insufficient memory **\n");
        exit(1);
    }
    spectra_all(s, u, v, w, out);
    spectra_write_frame(s, out, t);
    free(out);
}

void spectra_free(spectra *s)
{
    if (!s) return;
    fftw_destroy_plan(s->plan);
    fftw_plans_release();
    free(s->bin);
    free(s->weight);
    free(s->lap);
    free(s->ratio);
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
