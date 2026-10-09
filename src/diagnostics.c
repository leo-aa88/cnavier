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
    int n = cfg->forcing.kolmogorov_n, ny = cfg->ny;
    double k, k1, Q;
    periodic_symbols sym;

    if (!cfg->periodic || cfg->forcing.kolmogorov_amp == 0.0) return 1.0;
    k = 2.0 * PI * n / (ny * cfg->dy);
    // DY e^{iky} = i k1 e^{iky}, DY2 e^{iky} = -Q e^{iky}
    periodic_symbols_of(cfg, &sym);
    k1 = sym.d1y_im[n % ny];
    Q = -sym.d2y[n % ny];
    periodic_symbols_free(&sym);
    return k * k1 / Q;
}

flow_integrals compute_integrals(const solver_config *cfg, mtrx u, mtrx v, mtrx w, double *wx, double *wy)
{
    if (cfg->fourier)
    {
        fourier_ops *f = fourier_setup(cfg);
        flow_integrals r = compute_integrals_with(cfg, f, u, v, w, wx, wy);
        fourier_free(f);
        return r;
    }
    return compute_integrals_with(cfg, NULL, u, v, w, wx, wy);
}

flow_integrals compute_integrals_with(const solver_config *cfg, fourier_ops *four, mtrx u, mtrx v, mtrx w, double *wx,
                                      double *wy)
{
    int i, j, nx = cfg->nx, ny = cfg->ny;
    double sE = 0.0, sZ = 0.0, sP = 0.0, sI = 0.0, norm;
    double Ly = cfg->dy * (cfg->periodic ? ny : ny - 1), kk = 2.0 * PI * cfg->forcing.kolmogorov_n / Ly;
    flow_integrals r;

    if (four)
        fourier_derivatives(four, w.M, wx, wy, NULL);
    else
    {
        spmv(*cfg->DX, w.M, wx);
        spmv(*cfg->DY, w.M, wy);
    }
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
    fftw_plan plan;       // r2c of `in` into a spectrum
    unsigned long id;     // distinct for every spectra_setup() of the run
    fourier_ops *fourier; // the operators in Fourier space (cfg.fourier only)
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
    static unsigned long next_id = 0;
    s->cfg = *cfg;
    s->id = ++next_id;
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
    periodic_symbols sym;
    periodic_symbols_of(cfg, &sym);
    for (j = 0; j < s->kx; j++)
        lx[j] = sym.d2x[j];
    for (i = 0; i < ny; i++)
        ly[i] = sym.d2y[i];
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
                                     : (sym.d1x_re[j] * sym.d1x_re[j] + sym.d1x_im[j] * sym.d1x_im[j] +
                                        sym.d1y_re[i] * sym.d1y_re[i] + sym.d1y_im[i] * sym.d1y_im[i]) /
                                           -s->lap[idx];
            if (k > kmax) kmax = k;
        }
    s->bins = (int)floor(kmax / s->dk + 0.5) + 1;
    free(lx);
    free(ly);
    periodic_symbols_free(&sym);
    s->fourier = cfg->fourier ? fourier_setup(cfg) : NULL;
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
    // skew-symmetric or 3/2-padded form (nonlinear_correction())
    if (s->fourier)
        fourier_derivatives(s->fourier, w.M, s->wx, s->wy, NULL);
    else
    {
        spmv(*s->cfg.DX, w.M, s->wx);
        spmv(*s->cfg.DY, w.M, s->wy);
    }
    for (k = 0; k < n; k++)
        s->nl[k] = -(u.M[k] * s->wx[k] + v.M[k] * s->wy[k]);
    if (nonlinear_corrected(&s->cfg))
        nonlinear_correction(&s->cfg, s->fourier, u.M, v.M, w.M, s->wx, s->wy, s->nl, s->wx, s->wy);

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

// spectra_dissipation() from the spectrum of w already in s->wh
static void dissipation_from_spectrum(spectra *s, double *DE, double *DZ, double *FE, double *FZ)
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

void spectra_dissipation(spectra *s, mtrx w, double *DE, double *DZ, double *FE, double *FZ)
{
    transform(s, w.M, s->wh);
    dissipation_from_spectrum(s, DE, DZ, FE, FZ);
}

void spectra_all(spectra *s, mtrx u, mtrx v, mtrx w, double *out)
{
    size_t B = (size_t)s->bins;
    // spectra_compute() leaves the spectrum of w in s->wh: four transforms, not five
    spectra_compute(s, u, v, w, out, out + B, out + 2 * B, out + 3 * B);
    dissipation_from_spectrum(s, out + 4 * B, out + 5 * B, out + 6 * B, out + 7 * B);
}

unsigned long spectra_id(const spectra *s)
{
    return s->id;
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
    fourier_free(s->fourier);
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
