// Polar grid around a circular cylinder (polar.h)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include "polar.h"
#include "poisson.h"
#include "utils.h"

static void *xmalloc(size_t bytes)
{
    void *p = malloc(bytes);
    if (!p)
    {
        printf("** Error: insufficient memory **\n");
        exit(1);
    }
    return p;
}

int polar_square_nr(int ntheta, double R, double r_inf)
{
    return (int)lround(log(r_inf / R) / (2.0 * PI / ntheta)) + 1;
}

polar_grid polar_grid_make(int ntheta, int nr, double R, double r_inf)
{
    polar_grid g;
    int i, j;

    if (nr <= 0) nr = polar_square_nr(ntheta, R, r_inf);
    if (ntheta < 8 || nr < POLAR_MIN_NR || !(R > 0.0) || !(r_inf > R))
    {
        printf("** Error: a polar grid needs ntheta >= 8, nr >= %d (here %d) and 0 < R < r_inf **\n", POLAR_MIN_NR,
               nr);
        exit(1);
    }
    g.ntheta = ntheta;
    g.nr = nr;
    g.R = R;
    g.r_inf = r_inf;
    g.dtheta = 2.0 * PI / ntheta;
    g.dxi = log(r_inf / R) / (nr - 1);
    g.r = (double *)xmalloc((size_t)nr * sizeof(double));
    g.inv_r = (double *)xmalloc((size_t)nr * sizeof(double));
    g.inv_r2 = (double *)xmalloc((size_t)nr * sizeof(double));
    for (i = 0; i < nr; i++)
    {
        g.r[i] = i == nr - 1 ? r_inf : R * exp(i * g.dxi);
        g.inv_r[i] = 1.0 / g.r[i];
        g.inv_r2[i] = g.inv_r[i] * g.inv_r[i];
    }
    g.cos_t = (double *)xmalloc((size_t)ntheta * sizeof(double));
    g.sin_t = (double *)xmalloc((size_t)ntheta * sizeof(double));
    // Mirror-symmetric to the last bit, cos(theta_j) = cos(theta_{-j}) and
    // sin(theta_j) = -sin(theta_{-j}), with cos = 0 exactly at theta = pi/2
    // and 3 pi/2: rounded values there would treat one as inflow and the
    // other as outflow on the outer boundary, and break the symmetry of the flow
    for (j = 0; j <= ntheta / 2; j++)
    {
        g.cos_t[j] = 4 * j == ntheta ? 0.0 : cos(j * g.dtheta);
        g.sin_t[j] = 2 * j == ntheta ? 0.0 : sin(j * g.dtheta);
    }
    for (j = ntheta / 2 + 1; j < ntheta; j++)
    {
        g.cos_t[j] = g.cos_t[ntheta - j];
        g.sin_t[j] = -g.sin_t[ntheta - j];
    }
    g.psi_out = (double *)xmalloc((size_t)ntheta * sizeof(double));
    g.w_out = (double *)xmalloc((size_t)ntheta * sizeof(double));
    g.wall_ut = 0.0;
    g.cd = 0.0;
    g.Re = 0.0;
    g.version = 0;
    polar_freestream_set(&g);
    return g;
}

void polar_grid_free(polar_grid *g)
{
    free(g->r);
    free(g->inv_r);
    free(g->inv_r2);
    free(g->cos_t);
    free(g->sin_t);
    free(g->psi_out);
    free(g->w_out);
    memset(g, 0, sizeof(*g));
}

void polar_freestream_set(polar_grid *g)
{
    int j;
    g->outer_bc = POLAR_FREESTREAM;
    for (j = 0; j < g->ntheta; j++)
    {
        g->psi_out[j] = g->r_inf * g->sin_t[j];
        g->w_out[j] = 0.0;
    }
    g->version++;
}

void polar_dirichlet_set(polar_grid *g, const double *psi_out, const double *w_out)
{
    memcpy(g->psi_out, psi_out, (size_t)g->ntheta * sizeof(double));
    memcpy(g->w_out, w_out, (size_t)g->ntheta * sizeof(double));
    g->outer_bc = POLAR_DIRICHLET;
    g->version++;
}

void polar_farfield_set(polar_grid *g, double cd, double Re)
{
    int j, m, n = g->ntheta, kt = n / 2 + 1;
    double r = g->r_inf, amp = cd * sqrt(Re) / (4.0 * sqrt(PI)) / sqrt(r), mean = 0.0;
    double *f = (double *)fftw_malloc((size_t)n * sizeof(double));
    fftw_complex *fh = (fftw_complex *)fftw_malloc((size_t)kt * sizeof(fftw_complex));
    fftw_plan fwd, bwd;

    if (!f || !fh)
    {
        printf("** Error: insufficient memory **\n");
        exit(1);
    }
    fftw_plans_hold();
    fwd = fftw_plan_dft_r2c_1d(n, f, fh, FFTW_ESTIMATE);
    bwd = fftw_plan_dft_c2r_1d(n, fh, f, FFTW_ESTIMATE);

    // r_inf times the radial velocity of the deficit, u1.e_r = u1x cos + u1y sin
    for (j = 0; j < n; j++)
    {
        double c = g->cos_t[j], s = g->sin_t[j], d = amp * exp(-0.5 * Re * r * (1.0 - c));
        f[j] = r * (-d * c - 0.5 * s * d * s);
        mean += f[j] / n;
    }
    // The mean is the flux the deficit removes; the uniform outflow returns
    // it. The rest is integrated in theta mode by mode (divided by i m), with
    // zero mean, so that psi = 0 at theta = 0 (the integrand is even).
    for (j = 0; j < n; j++)
        f[j] -= mean;
    fftw_execute(fwd);
    fh[0][0] = fh[0][1] = 0.0;
    for (m = 1; m < kt; m++)
    {
        double re = fh[m][0], im = fh[m][1];
        if (2 * m == n)
            fh[m][0] = fh[m][1] = 0.0; // the Nyquist mode has no derivative
        else
        {
            fh[m][0] = im / m / n;
            fh[m][1] = -re / m / n;
        }
    }
    fftw_execute(bwd);
    for (j = 0; j < n; j++)
    {
        g->psi_out[j] = r * g->sin_t[j] + f[j];
        g->w_out[j] = 0.0;
    }
    fftw_destroy_plan(fwd);
    fftw_destroy_plan(bwd);
    fftw_free(f);
    fftw_free(fh);
    fftw_plans_release();
    g->outer_bc = POLAR_FARFIELD;
    g->cd = cd;
    g->Re = Re;
    g->version++;
}

double polar_potential_psi(const polar_grid *g, double r, double theta)
{
    double R2 = g->R * g->R, A = g->r_inf / (g->r_inf - R2 / g->r_inf);
    return A * (r - R2 / r) * sin(theta);
}

// ---------------------------------------------------------------------------
// Poisson solver
// ---------------------------------------------------------------------------

struct polar_poisson
{
    const polar_grid *g;
    int nt, nr, kt, order;
    double *in;         // nr + 1 rows of nt: f (rows 0..nr-1), psi_out (row nr)
    fftw_complex *spec; // their transforms, nr + 1 rows of kt
    fftw_plan forward;  // the nr + 1 rows
    fftw_plan backward; // nr rows back
    double *off;        // off-diagonal of each mode's tridiagonal system, [kt]
    double *cp, *dinv;  // Thomas factors of each mode: modified super-diagonal and inverse pivots, [kt * nr]
    double *work;       // elimination of each mode, [kt * nr]
};

polar_poisson *polar_poisson_setup(const polar_grid *g, const smtrx *DX2, int order)
{
    polar_poisson *s = (polar_poisson *)xmalloc(sizeof(polar_poisson));
    int nt = g->ntheta, nr = g->nr, kt = nt / 2 + 1, m, i, e;
    double h2 = g->dxi * g->dxi;

    if ((order != 2 && order != 4) || g->nr < POLAR_MIN_NR)
    {
        printf("** Error: the polar Poisson operator is of order 2 or 4, on nr >= %d **\n", POLAR_MIN_NR);
        exit(1);
    }
    s->g = g;
    s->nt = nt;
    s->nr = nr;
    s->kt = kt;
    s->order = order;
    s->in = (double *)fftw_malloc((size_t)nt * (nr + 1) * sizeof(double));
    s->spec = (fftw_complex *)fftw_malloc((size_t)kt * (nr + 1) * sizeof(fftw_complex));
    if (!s->in || !s->spec)
    {
        printf("** Error: insufficient memory **\n");
        exit(1);
    }
    s->off = (double *)xmalloc((size_t)kt * sizeof(double));
    s->cp = (double *)xmalloc((size_t)kt * nr * sizeof(double));
    s->dinv = (double *)xmalloc((size_t)kt * nr * sizeof(double));
    s->work = (double *)xmalloc((size_t)kt * nr * sizeof(double));

    fftw_plans_hold();
    s->forward = fftw_plan_many_dft_r2c(1, &nt, nr + 1, s->in, NULL, 1, nt, s->spec, NULL, 1, kt, FFTW_ESTIMATE);
    s->backward = fftw_plan_many_dft_c2r(1, &nt, nr, s->spec, NULL, 1, kt, s->in, NULL, 1, nt, FFTW_ESTIMATE);
    if (!s->forward || !s->backward)
    {
        printf("** Error: FFTW could not plan the polar transform **\n");
        exit(1);
    }

    for (m = 0; m < kt; m++)
    {
        // The symbol of DX2 on exp(i m theta), from its first row: theta
        // stencil in columns 0..nt-1
        double sym = 0.0, off, diag, denom;
        for (e = DX2->row_ptr[0]; e < DX2->row_ptr[1]; e++)
            sym += DX2->values[e] * cos(2.0 * PI * m * DX2->col_idx[e] / nt);
        if (order == 2)
        {
            off = 1.0 / h2;
            diag = -2.0 / h2 + sym;
        }
        else
        {
            // delta^2 psi + sym (1 + h^2 delta^2 / 12) psi = (1 + h^2 delta^2 / 12) f
            off = 1.0 / h2 + sym / 12.0;
            diag = -2.0 / h2 + 10.0 * sym / 12.0;
        }
        s->off[m] = off;
        // Unknowns on rows 1..nr-2
        for (i = 1; i < nr - 1; i++)
        {
            denom = diag - (i > 1 ? off * s->cp[m * nr + i - 1] : 0.0);
            s->dinv[m * nr + i] = 1.0 / denom;
            s->cp[m * nr + i] = off / denom;
        }
    }
    return s;
}

void polar_poisson_free(polar_poisson *s)
{
    if (!s) return;
    fftw_destroy_plan(s->forward);
    fftw_destroy_plan(s->backward);
    fftw_free(s->in);
    fftw_free(s->spec);
    free(s->off);
    free(s->cp);
    free(s->dinv);
    free(s->work);
    free(s);
    fftw_plans_release();
}

// The solve with s->in rows 1..nr-2 holding f
static void solve_in(polar_poisson *s, const double *psi_out, double *psi)
{
    int nt = s->nt, nr = s->nr, kt = s->kt, i, j, m;
    double scale = 1.0 / nt;

    if (s->order == 4)
        for (j = 0; j < nt; j++)
        {
            // f on the boundary rows, extrapolated (cubic) from the interior
            double *c = s->in + j;
            size_t r = (size_t)nt; // one row
            c[0] = 4.0 * c[r] - 6.0 * c[2 * r] + 4.0 * c[3 * r] - c[4 * r];
            c[(nr - 1) * r] = 4.0 * c[(nr - 2) * r] - 6.0 * c[(nr - 3) * r] + 4.0 * c[(nr - 4) * r] - c[(nr - 5) * r];
        }
    else
        for (j = 0; j < nt; j++)
            s->in[j] = s->in[(size_t)(nr - 1) * nt + j] = 0.0;
    memcpy(s->in + (size_t)nr * nt, psi_out, (size_t)nt * sizeof(double));
    fftw_execute(s->forward);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (nt * nr >= OMP_MIN_WORK)
#endif
    for (m = 0; m < kt; m++)
    {
        int c, k;
        const double *cp = s->cp + (size_t)m * nr, *dinv = s->dinv + (size_t)m * nr;
        double off = s->off[m], *d = s->work + (size_t)m * nr;
        fftw_complex *col = s->spec + m; // row k at col[k * kt]

        for (c = 0; c < 2; c++)
        {
            double g = col[(size_t)nr * kt][c];
            // Right-hand side of rows 1..nr-2 into d, then forward elimination
            for (k = 1; k < nr - 1; k++)
                d[k] = s->order == 2 ? col[(size_t)k * kt][c]
                                     : (col[(size_t)(k - 1) * kt][c] + 10.0 * col[(size_t)k * kt][c] +
                                        col[(size_t)(k + 1) * kt][c]) /
                                           12.0;
            d[nr - 2] -= off * g;
            d[1] *= dinv[1];
            for (k = 2; k < nr - 1; k++)
                d[k] = (d[k] - off * d[k - 1]) * dinv[k];
            // Back substitution
            for (k = nr - 3; k >= 1; k--)
                d[k] -= cp[k] * d[k + 1];
            col[0][c] = 0.0;
            for (k = 1; k < nr - 1; k++)
                col[(size_t)k * kt][c] = d[k];
            col[(size_t)(nr - 1) * kt][c] = g;
        }
    }

    fftw_execute(s->backward);
    for (j = 0; j < nt; j++)
    {
        psi[j] = 0.0;
        psi[(size_t)(nr - 1) * nt + j] = psi_out[j];
    }
    for (i = 1; i < nr - 1; i++)
        for (j = 0; j < nt; j++)
            psi[(size_t)i * nt + j] = s->in[(size_t)i * nt + j] * scale;
}

void polar_poisson_solve_f(polar_poisson *s, const double *f, const double *psi_out, double *psi)
{
    memcpy(s->in, f, (size_t)s->nt * s->nr * sizeof(double));
    solve_in(s, psi_out, psi);
}

void polar_poisson_solve(polar_poisson *s, const double *w, const double *psi_out, double *psi)
{
    int i, j, nt = s->nt;
    for (i = 0; i < s->nr; i++)
    {
        double r2 = s->g->r[i] * s->g->r[i];
        for (j = 0; j < nt; j++)
            s->in[(size_t)i * nt + j] = -r2 * w[(size_t)i * nt + j];
    }
    solve_in(s, psi_out, psi);
}

// ---------------------------------------------------------------------------
// Velocity and boundary conditions
// ---------------------------------------------------------------------------

void polar_velocity(const polar_grid *g, const smtrx *DX, const smtrx *DY, const double *psi, double *u, double *v)
{
    int i, j, nt = g->ntheta;

    spmv(*DY, (double *)psi, u);
    spmv(*DX, (double *)psi, v);
    for (i = 0; i < g->nr; i++)
        for (j = 0; j < nt; j++)
        {
            u[i * nt + j] *= -g->inv_r[i];
            v[i * nt + j] *= g->inv_r[i];
        }
}

void polar_apply_bc(const polar_grid *g, mtrx u, mtrx v)
{
    int j;
    for (j = 0; j < g->ntheta; j++)
    {
        MAt(u, 0, j) = g->wall_ut;
        MAt(v, 0, j) = 0.0;
    }
}

// Row r of A applied to f x, with x on every row of the grid scaled by s[row]
static double csr_row_scaled(const smtrx *A, const double *x, const double *s, int nt, int r)
{
    int k;
    double sum = 0.0;
    for (k = A->row_ptr[r]; k < A->row_ptr[r + 1]; k++)
        sum += A->values[k] * s[A->col_idx[k] / nt] * x[A->col_idx[k]];
    return sum;
}

void polar_boundary_vorticity(const polar_grid *g, int closure, const smtrx *DX, const smtrx *DY, mtrx w,
                              const mtrx *psi, mtrx u, mtrx v)
{
    int j, nt = g->ntheta, N = g->nr - 1;
    double R = g->R, h = g->dxi;

    for (j = 0; j < nt; j++)
    {
        if (closure == 1)
        {
            // Briley's third-order formula in xi (as in fluiddyn.c): -psi_xixi
            // from psi on rows 0..3 and psi_xi = -R u_theta on the wall
            const double *p = psi->M + j;
            double U = -R * g->wall_ut;
            size_t r = (size_t)nt; // one row
            MAt(w, 0, j) = ((85.0 * p[0] - 108.0 * p[r] + 27.0 * p[2 * r] - 4.0 * p[3 * r]) / (18.0 * h * h) +
                            11.0 * U / (3.0 * h)) *
                           g->inv_r2[0];
        }
        else
            MAt(w, 0, j) = g->inv_r2[0] *
                           (csr_row_scaled(DY, u.M, g->r, nt, j) - csr_row_scaled(DX, v.M, g->r, nt, j));

        // Outer boundary
        if (g->outer_bc == POLAR_DIRICHLET || g->cos_t[j] < 0.0)
            MAt(w, N, j) = g->w_out[j];
        else
            MAt(w, N, j) = MAt(w, N - 1, j); // outflow: zero normal gradient
    }
}

void polar_divergence(const polar_grid *g, const smtrx *DX, const smtrx *DY, const double *u, const double *v,
                      double *div, double *s1, double *s2)
{
    int i, j, nt = g->ntheta, n = g->ntheta * g->nr;

    for (i = 0; i < g->nr; i++)
        for (j = 0; j < nt; j++)
        {
            s1[i * nt + j] = g->r[i] * v[i * nt + j];
            s2[i * nt + j] = g->r[i] * u[i * nt + j];
        }
    spmv(*DY, s1, div);
    spmv(*DX, s2, s1);
    for (i = 0; i < n; i++)
        div[i] = (div[i] + s1[i]) * g->inv_r2[i / nt];
}

// ---------------------------------------------------------------------------
// Diagnostics
// ---------------------------------------------------------------------------

polar_forces polar_wall_forces(const polar_grid *g, const smtrx *DY, const double *w, double Re, double *p_wall)
{
    int j, nt = g->ntheta;
    double nu = 1.0 / Re, dt = g->dtheta;
    double *dwdxi = (double *)xmalloc((size_t)nt * sizeof(double));
    polar_forces f = {0};

    // On the wall: friction tau = nu w; theta-momentum gives dp/dtheta =
    // R nu dw/dr = nu dw/dxi. The force per unit length over (1/2) U^2 D,
    // D = 2R, is the integral over theta of (-p cos - nu w sin, -p sin +
    // nu w cos), and by parts -p cos -> p' sin, -p sin -> -p' cos. The
    // trapezoidal rule is spectrally accurate for these periodic integrands.
    for (j = 0; j < nt; j++)
    {
        double ww = w[j], dp;
        int k;
        dwdxi[j] = 0.0;
        for (k = DY->row_ptr[j]; k < DY->row_ptr[j + 1]; k++)
            dwdxi[j] += DY->values[k] * w[DY->col_idx[k]];
        dp = nu * dwdxi[j];
        f.cdp += dp * g->sin_t[j] * dt;
        f.cdf += -nu * ww * g->sin_t[j] * dt;
        f.cl += (-dp * g->cos_t[j] + nu * ww * g->cos_t[j]) * dt;
        f.dpdt_mean += dp / nt;
    }
    f.cd = f.cdp + f.cdf;

    if (p_wall)
    {
        int m, kt = nt / 2 + 1;
        double *in = (double *)fftw_malloc((size_t)nt * sizeof(double));
        fftw_complex *sp = (fftw_complex *)fftw_malloc((size_t)kt * sizeof(fftw_complex));
        fftw_plan fwd, bwd;
        if (!in || !sp)
        {
            printf("** Error: insufficient memory **\n");
            exit(1);
        }
        fftw_plans_hold();
        fwd = fftw_plan_dft_r2c_1d(nt, in, sp, FFTW_ESTIMATE);
        bwd = fftw_plan_dft_c2r_1d(nt, sp, in, FFTW_ESTIMATE);
        for (j = 0; j < nt; j++)
            in[j] = nu * dwdxi[j];
        fftw_execute(fwd);
        // p = integral of dp/dtheta: mode m divided by i m, zero mean
        sp[0][0] = sp[0][1] = 0.0;
        for (m = 1; m < kt; m++)
        {
            double re = sp[m][0], im = sp[m][1];
            if (2 * m == nt)
                sp[m][0] = sp[m][1] = 0.0;
            else
            {
                sp[m][0] = im / m / nt;
                sp[m][1] = -re / m / nt;
            }
        }
        fftw_execute(bwd);
        memcpy(p_wall, in, (size_t)nt * sizeof(double));
        fftw_destroy_plan(fwd);
        fftw_destroy_plan(bwd);
        fftw_free(in);
        fftw_free(sp);
        fftw_plans_release();
    }
    free(dwdxi);
    return f;
}

// Lagrange interpolation of degree 3 through (x[k], y[k]), k = 0..3, at t
static double lagrange4(const double *x, const double *y, double t)
{
    int a, b;
    double sum = 0.0;
    for (a = 0; a < 4; a++)
    {
        double l = 1.0;
        for (b = 0; b < 4; b++)
            if (b != a) l *= (t - x[b]) / (x[a] - x[b]);
        sum += l * y[a];
    }
    return sum;
}

// The root in [x[k], x[k+1]] of the cubic through the four points (y[k] and
// y[k+1] of opposite signs), by bisection
static double cubic_root(const double *x, const double *y, int k)
{
    double lo = x[k], hi = x[k + 1], flo = y[k];
    int it;
    for (it = 0; it < 100; it++)
    {
        double mid = 0.5 * (lo + hi), fm = lagrange4(x, y, mid);
        if ((fm > 0.0) == (flo > 0.0))
        {
            lo = mid;
            flo = fm;
        }
        else
            hi = mid;
    }
    return 0.5 * (lo + hi);
}

double polar_separation_angle(const polar_grid *g, const double *w)
{
    int j, k, nt = g->ntheta;
    double x[4], y[4];

    // In the bubble the flow along the rear of the upper surface runs
    // upstream (w > 0 on the wall); ahead of separation w < 0
    if (!(w[1] > 0.0)) return -1.0;
    for (j = 2; j <= nt / 2; j++)
        if (!(w[j] > 0.0)) break;
    if (j > nt / 2) return -1.0;
    for (k = 0; k < 4; k++)
    {
        int jj = j - 2 + k;
        x[k] = jj * g->dtheta;
        y[k] = w[(jj % nt + nt) % nt];
    }
    return 180.0 - cubic_root(x, y, 1) * 180.0 / PI;
}

double polar_wake_length(const polar_grid *g, const double *v)
{
    int i, k, nt = g->ntheta, nr = g->nr;
    double x[4], y[4];

    // u_r on theta = 0 is column 0
    if (!(v[nt] < 0.0)) return 0.0;
    for (i = 2; i < nr; i++)
        if (!(v[(size_t)i * nt] < 0.0)) break;
    if (i >= nr) return (g->r_inf - g->R) / (2.0 * g->R);
    // The sign changes between rows i - 1 and i: the cubic through rows
    // i - 2 .. i + 1, or through the last four rows when i is the last
    int first = i + 1 < nr ? i - 2 : nr - 4;
    for (k = 0; k < 4; k++)
    {
        x[k] = (first + k) * g->dxi;
        y[k] = v[(size_t)(first + k) * nt];
    }
    return (g->R * exp(cubic_root(x, y, i - 1 - first)) - g->R) / (2.0 * g->R);
}

int polar_vortex_centre(const polar_grid *g, const double *psi, double *a, double *b)
{
    int i, j, nt = g->ntheta, bi = -1, bj = -1;
    double best = 0.0, xi, th, r, fm, f0, fp, den;

    for (i = 1; i < g->nr - 1; i++)
        for (j = 1; j * g->dtheta < 0.5 * PI; j++)
            if (psi[i * nt + j] < best)
            {
                best = psi[i * nt + j];
                bi = i;
                bj = j;
            }
    if (bi < 0) return 0;
    // Vertex of the parabola through the three points in each direction
    xi = bi * g->dxi;
    fm = psi[(bi - 1) * nt + bj], f0 = psi[bi * nt + bj], fp = psi[(bi + 1) * nt + bj];
    den = fm - 2.0 * f0 + fp;
    if (den != 0.0) xi += 0.5 * g->dxi * (fm - fp) / den;
    th = bj * g->dtheta;
    fm = psi[bi * nt + bj - 1], fp = psi[bi * nt + bj + 1];
    den = fm - 2.0 * f0 + fp;
    if (den != 0.0) th += 0.5 * g->dtheta * (fm - fp) / den;
    r = g->R * exp(xi);
    *a = (r * cos(th) - g->R) / (2.0 * g->R);
    *b = 2.0 * r * sin(th) / (2.0 * g->R);
    return 1;
}

int polar_interpolate(const polar_grid *g, const double *u, const double *v, const double *w, double x, double y,
                      double *ux, double *uy, double *wv)
{
    int nt = g->ntheta, nr = g->nr, i0, j0, a, b;
    double r = sqrt(x * x + y * y), th = atan2(y, x), xi, wx[4], wt[4], sut = 0.0, sur = 0.0, sw = 0.0;

    if (th < 0.0) th += 2.0 * PI;
    if (!(r >= g->R) || !(r <= g->r_inf)) return 0;
    xi = log(r / g->R);
    i0 = (int)floor(xi / g->dxi) - 1;
    if (i0 < 0) i0 = 0;
    if (i0 > nr - 4) i0 = nr - 4;
    j0 = (int)floor(th / g->dtheta) - 1;
    // Lagrange weights in each direction
    for (a = 0; a < 4; a++)
    {
        double lx = 1.0, lt = 1.0;
        for (b = 0; b < 4; b++)
            if (b != a)
            {
                lx *= (xi - (i0 + b) * g->dxi) / ((a - b) * g->dxi);
                lt *= (th - (j0 + b) * g->dtheta) / ((a - b) * g->dtheta);
            }
        wx[a] = lx;
        wt[a] = lt;
    }
    for (a = 0; a < 4; a++)
        for (b = 0; b < 4; b++)
        {
            int k = (i0 + a) * nt + ((j0 + b) % nt + nt) % nt;
            double c = wx[a] * wt[b];
            sut += c * u[k];
            sur += c * v[k];
            sw += c * w[k];
        }
    if (ux) *ux = sur * cos(th) - sut * sin(th);
    if (uy) *uy = sur * sin(th) + sut * cos(th);
    if (wv) *wv = sw;
    return 1;
}

double polar_asymmetry(const polar_grid *g, const double *w)
{
    int i, j, nt = g->ntheta;
    double e = 0.0;
    for (i = 0; i < g->nr; i++)
        for (j = 0; j < nt; j++)
            e = fmax(e, fabs(w[i * nt + j] + w[i * nt + (nt - j) % nt]));
    return e;
}

void printvtk_polar(const polar_grid *g, const mtrx *f, const mtrx *u, const mtrx *v, const char *title)
{
    int i, j, nt = g->ntheta, nr = g->nr, count;
    char name[96];
    FILE *pf;

    count = output_frame(title, ".vtk");
    snprintf(name, sizeof(name), "./output/%s-1-%d.vtk", title, count);
    if ((pf = fopen(name, "w")) == NULL)
    {
        printf("\nError while opening file\n");
        exit(1);
    }
    printf("%s\n", name);
    fprintf(pf, "# vtk DataFile Version 2.0\n%s\nASCII\nDATASET STRUCTURED_GRID\n", title);
    fprintf(pf, "DIMENSIONS %d %d 1\nPOINTS %d double\n", nt + 1, nr, (nt + 1) * nr);
    for (i = 0; i < nr; i++)
        for (j = 0; j <= nt; j++)
            fprintf(pf, "%.9g %.9g 0\n", g->r[i] * g->cos_t[j % nt], g->r[i] * g->sin_t[j % nt]);
    fprintf(pf, "POINT_DATA %d\nSCALARS %s double\nLOOKUP_TABLE default\n", (nt + 1) * nr, title);
    for (i = 0; i < nr; i++)
        for (j = 0; j <= nt; j++)
            fprintf(pf, "%.9g\n", MAt(*f, i, j % nt));
    if (u && v)
    {
        fprintf(pf, "VECTORS velocity double\n");
        for (i = 0; i < nr; i++)
            for (j = 0; j <= nt; j++)
            {
                double ut = MAt(*u, i, j % nt), ur = MAt(*v, i, j % nt), c = g->cos_t[j % nt], s = g->sin_t[j % nt];
                fprintf(pf, "%.9g %.9g 0\n", ur * c - ut * s, ur * s + ut * c);
            }
    }
    fclose(pf);
}
