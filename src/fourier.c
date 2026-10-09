#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "fourier.h"
#include "fluiddyn.h"
#include "poisson.h"

// Symbol of a circulant operator from its first row: sum over the entries of
// value * exp(i 2 pi k offset / n), offset = column / step (whole grid points
// along the axis: step 1 for x, nx for y in a 2-D operator)
static void row_symbol(const smtrx *A, int n, int step, double *re, double *im)
{
    for (int k = 0; k < n; k++)
    {
        double sr = 0.0, si = 0.0;
        for (int e = A->row_ptr[0]; e < A->row_ptr[1]; e++)
        {
            int offset = A->col_idx[e] / step; // whole grid points along the axis
            double th = 2.0 * PI * (double)k * (double)offset / (double)n;
            sr += A->values[e] * cos(th);
            si += A->values[e] * sin(th);
        }
        re[k] = sr;
        if (im) im[k] = si;
    }
}

static double *alloc_doubles(int n)
{
    double *p = (double *)malloc((size_t)n * sizeof(double));
    if (!p)
    {
        printf("** Error: insufficient memory **\n");
        exit(1);
    }
    return p;
}

void periodic_symbols_of(const solver_config *cfg, periodic_symbols *s)
{
    int nx = cfg->nx, ny = cfg->ny, k;
    int one_d = cfg->D1x && cfg->D1y && cfg->D2x && cfg->D2y;

    if (!one_d && !(cfg->DX && cfg->DY && cfg->DX2 && cfg->DY2))
    {
        printf("** Error: the periodic symbols need the 1-D or the 2-D operators **\n");
        exit(1);
    }
    s->nx = nx;
    s->ny = ny;
    s->d1x_re = alloc_doubles(nx);
    s->d1x_im = alloc_doubles(nx);
    s->d2x = alloc_doubles(nx);
    s->d1y_re = alloc_doubles(ny);
    s->d1y_im = alloc_doubles(ny);
    s->d2y = alloc_doubles(ny);
    s->maskx = alloc_doubles(nx);
    s->masky = alloc_doubles(ny);
    row_symbol(one_d ? cfg->D1x : cfg->DX, nx, 1, s->d1x_re, s->d1x_im);
    row_symbol(one_d ? cfg->D2x : cfg->DX2, nx, 1, s->d2x, NULL);
    row_symbol(one_d ? cfg->D1y : cfg->DY, ny, one_d ? 1 : nx, s->d1y_re, s->d1y_im);
    row_symbol(one_d ? cfg->D2y : cfg->DY2, ny, one_d ? 1 : nx, s->d2y, NULL);
    // 2/3 rule: keep |k| < n/3 along each axis, so that the products of two
    // kept fields alias only onto modes that are cut
    for (k = 0; k < nx; k++)
        s->maskx[k] = (cfg->dealias != 1 || 3 * (k <= nx / 2 ? k : nx - k) < nx) ? 1.0 : 0.0;
    for (k = 0; k < ny; k++)
        s->masky[k] = (cfg->dealias != 1 || 3 * (k <= ny / 2 ? k : ny - k) < ny) ? 1.0 : 0.0;
}

void periodic_symbols_free(periodic_symbols *s)
{
    free(s->d1x_re);
    free(s->d1x_im);
    free(s->d2x);
    free(s->d1y_re);
    free(s->d1y_im);
    free(s->d2y);
    free(s->maskx);
    free(s->masky);
}

struct fourier_ops
{
    periodic_symbols sym;
    int nx, ny, kx;
    double *in;               // nx * ny real field
    fftw_complex *hat, *hat2; // spectra of the inputs, ny * kx
    fftw_complex *work;       // the spectrum an inverse transform consumes
    fftw_plan r2c, c2r;

    // 3/2 padding (cfg.dealias 2): the padded grid mx * my, its fields and spectra
    int mx, my, mkx;
    double *pf[4];      // u, v, DX w, DY w on the padded grid, then the product
    fftw_complex *phat; // a padded half spectrum, my * mkx
    fftw_plan pr2c, pc2r;
};

fourier_ops *fourier_setup(const solver_config *cfg)
{
    fourier_ops *f = (fourier_ops *)calloc(1, sizeof(fourier_ops));
    int nx = cfg->nx, ny = cfg->ny, kx = nx / 2 + 1;

    if (!f || !cfg->periodic)
    {
        printf("** Error: the Fourier operators need a periodic grid **\n");
        exit(1);
    }
    periodic_symbols_of(cfg, &f->sym);
    f->nx = nx;
    f->ny = ny;
    f->kx = kx;
    f->in = (double *)fftw_malloc((size_t)nx * ny * sizeof(double));
    f->hat = (fftw_complex *)fftw_malloc((size_t)kx * ny * sizeof(fftw_complex));
    f->hat2 = (fftw_complex *)fftw_malloc((size_t)kx * ny * sizeof(fftw_complex));
    f->work = (fftw_complex *)fftw_malloc((size_t)kx * ny * sizeof(fftw_complex));
    if (!f->in || !f->hat || !f->hat2 || !f->work)
    {
        printf("** Error: insufficient memory **\n");
        exit(1);
    }
    fftw_plans_hold();
    f->r2c = fftw_plan_dft_r2c_2d(ny, nx, f->in, f->hat, FFTW_ESTIMATE);
    f->c2r = fftw_plan_dft_c2r_2d(ny, nx, f->work, f->in, FFTW_ESTIMATE);
    if (!f->r2c || !f->c2r)
    {
        printf("** Error: FFTW could not plan the Fourier operators **\n");
        exit(1);
    }
    if (cfg->dealias == 2)
    {
        // Products of two fields with |k| < n/2 alias onto |k| >= m - n, so
        // m >= 3n/2 keeps every mode of the n grid free of aliasing
        f->mx = (3 * nx + 1) / 2;
        f->my = (3 * ny + 1) / 2;
        f->mkx = f->mx / 2 + 1;
        for (int q = 0; q < 4; q++)
            f->pf[q] = (double *)fftw_malloc((size_t)f->mx * f->my * sizeof(double));
        f->phat = (fftw_complex *)fftw_malloc((size_t)f->mkx * f->my * sizeof(fftw_complex));
        if (!f->pf[0] || !f->pf[1] || !f->pf[2] || !f->pf[3] || !f->phat)
        {
            printf("** Error: insufficient memory **\n");
            exit(1);
        }
        f->pr2c = fftw_plan_dft_r2c_2d(f->my, f->mx, f->pf[0], f->phat, FFTW_ESTIMATE);
        f->pc2r = fftw_plan_dft_c2r_2d(f->my, f->mx, f->phat, f->pf[0], FFTW_ESTIMATE);
        if (!f->pr2c || !f->pc2r)
        {
            printf("** Error: FFTW could not plan the padded transforms **\n");
            exit(1);
        }
    }
    return f;
}

void fourier_free(fourier_ops *f)
{
    if (!f) return;
    fftw_destroy_plan(f->r2c);
    fftw_destroy_plan(f->c2r);
    if (f->mx)
    {
        fftw_destroy_plan(f->pr2c);
        fftw_destroy_plan(f->pc2r);
        for (int q = 0; q < 4; q++)
            fftw_free(f->pf[q]);
        fftw_free(f->phat);
    }
    fftw_free(f->in);
    fftw_free(f->hat);
    fftw_free(f->hat2);
    fftw_free(f->work);
    periodic_symbols_free(&f->sym);
    free(f);
    fftw_plans_release();
}

const periodic_symbols *fourier_symbols(const fourier_ops *f)
{
    return &f->sym;
}

// out = spectrum of x
static void forward(fourier_ops *f, const double *x, fftw_complex *out)
{
    int n = f->nx * f->ny;
    for (int k = 0; k < n; k++)
        f->in[k] = x[k];
    fftw_execute_dft_r2c(f->r2c, f->in, out);
}

// out = the field of f->work, normalised (FFTW does not); consumes f->work
static void inverse(fourier_ops *f, double *out)
{
    int n = f->nx * f->ny;
    double inv = 1.0 / (double)n;
    fftw_execute_dft_c2r(f->c2r, f->work, f->in);
    for (int k = 0; k < n; k++)
        out[k] = f->in[k] * inv;
}

// What an operation multiplies mode (row i, column j) of the half spectrum by
enum op
{
    OP_DX,
    OP_DY,
    OP_LAP,
    OP_PSI, // -(DX2 + DY2)^-1, zero for the mean
    OP_FILTER
};

static void multiplier(const fourier_ops *f, enum op op, int i, int j, double *re, double *im)
{
    const periodic_symbols *s = &f->sym;
    double mask = s->maskx[j] * s->masky[i], lap = s->d2x[j] + s->d2y[i];
    *im = 0.0;
    switch (op)
    {
    case OP_DX:
        *re = s->d1x_re[j] * mask;
        *im = s->d1x_im[j] * mask;
        break;
    case OP_DY:
        *re = s->d1y_re[i] * mask;
        *im = s->d1y_im[i] * mask;
        break;
    case OP_LAP:
        *re = lap * mask;
        break;
    case OP_PSI:
        *re = (i == 0 && j == 0) ? 0.0 : -mask / lap;
        break;
    default:
        *re = mask;
        break;
    }
}

// f->work = multiplier(op) * src
static void apply(fourier_ops *f, const fftw_complex *src, enum op op)
{
    for (int i = 0; i < f->ny; i++)
        for (int j = 0; j < f->kx; j++)
        {
            int k = i * f->kx + j;
            double mr, mi;
            multiplier(f, op, i, j, &mr, &mi);
            f->work[k][0] = mr * src[k][0] - mi * src[k][1];
            f->work[k][1] = mr * src[k][1] + mi * src[k][0];
        }
}

void fourier_derivatives(fourier_ops *f, const double *w, double *wx, double *wy, double *lap)
{
    forward(f, w, f->hat);
    if (wx)
    {
        apply(f, f->hat, OP_DX);
        inverse(f, wx);
    }
    if (wy)
    {
        apply(f, f->hat, OP_DY);
        inverse(f, wy);
    }
    if (lap)
    {
        apply(f, f->hat, OP_LAP);
        inverse(f, lap);
    }
}

// u = DY psi, v = -DX psi from the spectrum of psi in f->hat2
static void velocity_from_hat(fourier_ops *f, double *u, double *v)
{
    int k, n = f->nx * f->ny;
    if (u)
    {
        apply(f, f->hat2, OP_DY);
        inverse(f, u);
    }
    if (v)
    {
        apply(f, f->hat2, OP_DX);
        inverse(f, v);
        for (k = 0; k < n; k++)
            v[k] = -v[k];
    }
}

void fourier_poisson(fourier_ops *f, const double *w, double *psi, double *u, double *v)
{
    forward(f, w, f->hat);
    apply(f, f->hat, OP_PSI);
    for (int k = 0; k < f->kx * f->ny; k++)
    {
        f->hat2[k][0] = f->work[k][0];
        f->hat2[k][1] = f->work[k][1];
    }
    if (psi) inverse(f, psi);
    velocity_from_hat(f, u, v);
}

void fourier_velocity(fourier_ops *f, const double *psi, double *u, double *v)
{
    forward(f, psi, f->hat2);
    velocity_from_hat(f, u, v);
}

void fourier_divergence(fourier_ops *f, const double *a, const double *b, double *out)
{
    int k;
    forward(f, a, f->hat);
    forward(f, b, f->hat2);
    apply(f, f->hat2, OP_DY);
    // work = DY b; add DX a mode by mode
    for (int i = 0; i < f->ny; i++)
        for (int j = 0; j < f->kx; j++)
        {
            double mr, mi;
            k = i * f->kx + j;
            multiplier(f, OP_DX, i, j, &mr, &mi);
            f->work[k][0] += mr * f->hat[k][0] - mi * f->hat[k][1];
            f->work[k][1] += mr * f->hat[k][1] + mi * f->hat[k][0];
        }
    inverse(f, out);
}

void fourier_power(fourier_ops *f, const double *w, int p, double *out)
{
    forward(f, w, f->hat);
    apply(f, f->hat, OP_LAP);
    // work = lap w^; multiply by (-lap)^(p-1) and the sign
    for (int i = 0; i < f->ny; i++)
        for (int j = 0; j < f->kx; j++)
        {
            int k = i * f->kx + j;
            double q = -(f->sym.d2x[j] + f->sym.d2y[i]), fac = -1.0;
            for (int r = 1; r < p; r++)
                fac *= q;
            f->work[k][0] *= fac;
            f->work[k][1] *= fac;
        }
    inverse(f, out);
}

void fourier_filter(fourier_ops *f, double *w)
{
    forward(f, w, f->hat);
    apply(f, f->hat, OP_FILTER);
    inverse(f, w);
}

// Row of the padded spectrum that holds row i (signed wavenumber m) of the n
// grid's, or -1 for the Nyquist row, which the padded product leaves out
static int padded_row(int i, int n, int m)
{
    int s = 2 * i <= n ? i : i - n;
    if (2 * (s < 0 ? -s : s) >= n) return -1;
    return s >= 0 ? s : s + m;
}

void fourier_nonlinear_padded(fourier_ops *f, const double *w, double *out)
{
    const periodic_symbols *s = &f->sym;
    int i, j, q, k, nx = f->nx, ny = f->ny, mx = f->mx, my = f->my, mkx = f->mkx, pn = mx * my;
    double inv = 1.0 / ((double)nx * ny);

    forward(f, w, f->hat);
    // u = DY psi, v = -DX psi, DX w, DY w on the padded grid: the n grid's
    // modes (the Nyquist ones left out) times the symbols, zero elsewhere
    for (q = 0; q < 4; q++)
    {
        for (k = 0; k < mkx * my; k++)
            f->phat[k][0] = f->phat[k][1] = 0.0;
        for (i = 0; i < ny; i++)
        {
            int r = padded_row(i, ny, my);
            if (r < 0) continue;
            for (j = 0; j < f->kx; j++)
            {
                double mr, mi, psi;
                if (2 * j >= nx) continue;
                psi = (i == 0 && j == 0) ? 0.0 : -1.0 / (s->d2x[j] + s->d2y[i]);
                if (q == 0)
                {
                    mr = s->d1y_re[i] * psi;
                    mi = s->d1y_im[i] * psi;
                }
                else if (q == 1)
                {
                    mr = -s->d1x_re[j] * psi;
                    mi = -s->d1x_im[j] * psi;
                }
                else if (q == 2)
                {
                    mr = s->d1x_re[j];
                    mi = s->d1x_im[j];
                }
                else
                {
                    mr = s->d1y_re[i];
                    mi = s->d1y_im[i];
                }
                const double *a = f->hat[i * f->kx + j];
                double *b = f->phat[r * mkx + j];
                b[0] = (mr * a[0] - mi * a[1]) * inv;
                b[1] = (mr * a[1] + mi * a[0]) * inv;
            }
        }
        fftw_execute_dft_c2r(f->pc2r, f->phat, f->pf[q]);
    }
    // The product on the padded grid, which represents it exactly up to the
    // n grid's modes; then back to the n grid, without the Nyquist modes
    for (k = 0; k < pn; k++)
        f->pf[0][k] = -(f->pf[0][k] * f->pf[2][k] + f->pf[1][k] * f->pf[3][k]);
    fftw_execute_dft_r2c(f->pr2c, f->pf[0], f->phat);
    for (i = 0; i < ny; i++)
    {
        int r = padded_row(i, ny, my);
        for (j = 0; j < f->kx; j++)
        {
            double *d = f->work[i * f->kx + j];
            if (r < 0 || 2 * j >= nx)
                d[0] = d[1] = 0.0;
            else
            {
                // padded r2c sums over mx*my points; inverse() divides by nx*ny
                d[0] = f->phat[r * mkx + j][0] * ((double)nx * ny / pn);
                d[1] = f->phat[r * mkx + j][1] * ((double)nx * ny / pn);
            }
        }
    }
    inverse(f, out);
}
