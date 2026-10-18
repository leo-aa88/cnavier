// Tests of the polar grid (polar.h): the Poisson solver, the velocity and
// its divergence, a manufactured solution of the full solver, the
// diagnostics, and the backend's CPU fallback

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "linearalg.h"
#include "finitediff.h"
#include "poisson.h"
#include "fluiddyn.h"
#include "backend.h"
#include "diagnostics.h"
#include "polar.h"
#include "test_polar.h"

static check_fn check;

// Operators of a polar grid: DX, DX2 in theta (periodic), DY, DY2 in xi
typedef struct
{
    polar_grid g;
    smtrx DX, DY, DX2, DY2;
} polar_ops;

static polar_ops ops_make(int ntheta, int nr, double R, double r_inf, int order)
{
    polar_ops p;
    p.g = polar_grid_make(ntheta, nr, R, r_inf);
    smtrx dt = SDiff1_periodic(ntheta, order, p.g.dtheta), dt2 = SDiff2_periodic(ntheta, order, p.g.dtheta);
    smtrx dx = SDiff1(p.g.nr, order, p.g.dxi), dx2 = SDiff2(p.g.nr, order, p.g.dxi);
    smtrx It = seye(ntheta), Ix = seye(p.g.nr);
    p.DX = skronecker(Ix, dt);
    p.DY = skronecker(dx, It);
    p.DX2 = skronecker(Ix, dt2);
    p.DY2 = skronecker(dx2, It);
    freesm(dt);
    freesm(dt2);
    freesm(dx);
    freesm(dx2);
    freesm(It);
    freesm(Ix);
    return p;
}

static void ops_free(polar_ops *p)
{
    freesm(p->DX);
    freesm(p->DY);
    freesm(p->DX2);
    freesm(p->DY2);
    polar_grid_free(&p->g);
}

// ---------------------------------------------------------------------------
// Poisson solver
// ---------------------------------------------------------------------------

// Error of the solve for w = 0 (potential flow) against the exact solution
static double potential_error(int ntheta, int order)
{
    polar_ops p = ops_make(ntheta, 0, 0.5, 4.0, 6);
    int n = p.g.ntheta * p.g.nr, i, j;
    double *w = (double *)calloc((size_t)n, sizeof(double)), *psi = (double *)malloc((size_t)n * sizeof(double));
    double e = 0.0;
    polar_poisson *s = polar_poisson_setup(&p.g, &p.DX2, order);

    polar_poisson_solve(s, w, p.g.psi_out, psi);
    for (i = 0; i < p.g.nr; i++)
        for (j = 0; j < p.g.ntheta; j++)
            e = fmax(e, fabs(psi[i * p.g.ntheta + j] - polar_potential_psi(&p.g, p.g.r[i], j * p.g.dtheta)));
    polar_poisson_free(s);
    free(w);
    free(psi);
    ops_free(&p);
    return e;
}

static void test_polar_poisson(void)
{
    printf("\nPolar Poisson solver\n");
    double e2a = potential_error(64, 2), e2b = potential_error(128, 2);
    double e4a = potential_error(64, 4), e4b = potential_error(128, 4);
    check("potential flow, order 2: error on 64 x 20", e2a, 1E-3);
    check("  observed order >= 1.9 (2 - log2 ratio)", 2.0 - log2(e2a / e2b), 0.1);
    check("potential flow, order 4: error on 64 x 20", e4a, 1E-5);
    check("  observed order >= 3.8 (4 - log2 ratio)", 4.0 - log2(e4a / e4b), 0.2);

    // Order 2 solves its discrete system exactly: the 3-point operator in xi
    // plus DX2 applied to the solution gives back f on the interior
    polar_ops p = ops_make(48, 17, 0.5, 6.0, 4);
    int nt = p.g.ntheta, nr = p.g.nr, n = nt * nr, i, j;
    double *f = (double *)calloc((size_t)n, sizeof(double)), *psi = (double *)calloc((size_t)n, sizeof(double));
    double *lap = (double *)calloc((size_t)n, sizeof(double)), res = 0.0, scale = 0.0, h2 = p.g.dxi * p.g.dxi;
    for (i = 0; i < nr; i++)
        for (j = 0; j < nt; j++)
            f[i * nt + j] = cos(3.0 * j * p.g.dtheta + 0.3) * exp(-p.g.r[i]) + sin(j * p.g.dtheta) * i;
    for (j = 0; j < nt; j++)
        p.g.psi_out[j] = 2.0 * sin(2.0 * j * p.g.dtheta) + 0.5;
    polar_poisson *s = polar_poisson_setup(&p.g, &p.DX2, 2);
    polar_poisson_solve_f(s, f, p.g.psi_out, psi);
    spmv(p.DX2, psi, lap);
    for (i = 1; i < nr - 1; i++)
        for (j = 0; j < nt; j++)
        {
            int k = i * nt + j;
            double a = lap[k] + (psi[k - nt] - 2.0 * psi[k] + psi[k + nt]) / h2;
            res = fmax(res, fabs(a - f[k]));
            scale = fmax(scale, fabs(f[k]));
        }
    check("order 2 solves its discrete system (residual / max |f|)", res / scale, 1E-12);
    double bmax = 0.0;
    for (j = 0; j < nt; j++)
        bmax = fmax(bmax, fabs(psi[j]) + fabs(psi[(nr - 1) * nt + j] - p.g.psi_out[j]));
    check("boundary values: psi = 0 on the wall, psi_out outside", bmax, 0.0);
    polar_poisson_free(s);

    // The velocity of any psi is divergence-free to round-off (DX and DY commute)
    double *u = (double *)malloc((size_t)n * sizeof(double)), *v = (double *)malloc((size_t)n * sizeof(double));
    double *s1 = (double *)malloc((size_t)n * sizeof(double)), *s2 = (double *)malloc((size_t)n * sizeof(double));
    double dmax = 0.0, umax = 0.0;
    polar_velocity(&p.g, &p.DX, &p.DY, psi, u, v);
    polar_divergence(&p.g, &p.DX, &p.DY, u, v, lap, s1, s2);
    for (i = 0; i < n; i++)
    {
        dmax = fmax(dmax, fabs(lap[i]));
        umax = fmax(umax, fabs(u[i]) + fabs(v[i]));
    }
    check("divergence of the velocity of psi / max |u|", dmax / umax, 1E-12);

    // Potential flow: u_r = A (1 - R^2/r^2) cos(theta) and the free stream far out
    {
        polar_ops q = ops_make(128, 0, 0.5, 8.0, 6);
        int qn = q.g.ntheta * q.g.nr, qt = q.g.ntheta;
        double *qw = (double *)calloc((size_t)qn, sizeof(double)), *qp = (double *)malloc((size_t)qn * sizeof(double));
        double *qu = (double *)malloc((size_t)qn * sizeof(double)), *qv = (double *)malloc((size_t)qn * sizeof(double));
        double A = q.g.r_inf / (q.g.r_inf - 0.25 / q.g.r_inf), e = 0.0, e_near = 0.0;
        polar_poisson *qs = polar_poisson_setup(&q.g, &q.DX2, 4);
        polar_poisson_solve(qs, qw, q.g.psi_out, qp);
        polar_velocity(&q.g, &q.DX, &q.DY, qp, qu, qv);
        // Interior rows: on the boundary rows DY is one-sided (the solver
        // imposes the wall velocity and does not use the outer row's)
        for (i = 1; i < q.g.nr - 1; i++)
            for (j = 0; j < qt; j++)
            {
                double r = q.g.r[i], ur = A * (1.0 - 0.25 / (r * r)) * q.g.cos_t[j],
                       ut = -A * (1.0 + 0.25 / (r * r)) * q.g.sin_t[j];
                double d = fabs(qu[i * qt + j] - ut) + fabs(qv[i * qt + j] - ur);
                // Rows 1, 2 and nr-3, nr-2 have the biased rows of DY
                if (i >= 3 && i < q.g.nr - 3)
                    e = fmax(e, d);
                else
                    e_near = fmax(e_near, d);
            }
        check("potential-flow velocity (u_r, u_theta signs), 128 x 57, interior", e, 1E-6);
        check("  and on the rows next to the boundaries (biased stencils)", e_near, 2E-3);
        polar_poisson_free(qs);
        free(qw);
        free(qp);
        free(qu);
        free(qv);
        ops_free(&q);
    }

    free(f);
    free(psi);
    free(lap);
    free(u);
    free(v);
    free(s1);
    free(s2);
    ops_free(&p);
}

// ---------------------------------------------------------------------------
// Manufactured solution of the full solver:
//   psi = g(t) F(r) T(theta),  F = (r - R)^2 (r_inf - r)^3,
//   T = sin(theta) + 0.4 cos(2 theta + 0.3) + 0.25 sin(3 theta)
// F(R) = F'(R) = 0 (no slip) and F, F', F'' vanish at r_inf, so psi = w = 0
// there (POLAR_DIRICHLET with zeros). The source makes it exact.
// ---------------------------------------------------------------------------

typedef struct
{
    double R, ri, Re;
    double c[6]; // F = sum c[k] r^k
    const polar_grid *g;
} mms_polar;

static double gt(double t)
{
    return 1.0 + 0.5 * sin(2.0 * t);
}
static double gdot(double t)
{
    return cos(2.0 * t);
}

// d-th derivative of T
static double Tth(double th, int d)
{
    // derivative d of sin(m th + p) = m^d sin(m th + p + d pi/2)
    return pow(1.0, d) * sin(th + d * 0.5 * PI) + 0.4 * pow(2.0, d) * sin(2.0 * th + 0.3 + 0.5 * PI + d * 0.5 * PI) +
           0.25 * pow(3.0, d) * sin(3.0 * th + d * 0.5 * PI);
}

// d-th derivative of F at r
static double Fr(const mms_polar *m, double r, int d)
{
    int k, q;
    double s = 0.0;
    for (k = d; k < 6; k++)
    {
        double c = m->c[k];
        for (q = 0; q < d; q++)
            c *= (k - q);
        s += c * pow(r, k - d);
    }
    return s;
}

static void mms_polar_init(mms_polar *m, double R, double ri, double Re)
{
    // (r - R)^2 (ri - r)^3 expanded: multiply the factors one at a time
    double p[6] = {1, 0, 0, 0, 0, 0}, q[6];
    int k, f;
    double roots[5] = {R, R, ri, ri, ri}, sign = -1.0; // (ri - r)^3 = -(r - ri)^3
    for (f = 0; f < 5; f++)
    {
        for (k = 0; k < 6; k++)
            q[k] = (k > 0 ? p[k - 1] : 0.0) - roots[f] * p[k];
        memcpy(p, q, sizeof(p));
    }
    for (k = 0; k < 6; k++)
        m->c[k] = sign * p[k];
    m->R = R;
    m->ri = ri;
    m->Re = Re;
}

// Exact w, u = u_theta, v = u_r at (r, th, t)
static void mms_polar_exact(const mms_polar *m, double r, double th, double t, double *w, double *u, double *v)
{
    double F = Fr(m, r, 0), F1 = Fr(m, r, 1), F2 = Fr(m, r, 2), g = gt(t);
    double L1 = F2 + F1 / r, G = F / (r * r);
    if (w) *w = -g * (L1 * Tth(th, 0) + G * Tth(th, 2));
    if (u) *u = -g * F1 * Tth(th, 0);
    if (v) *v = g * F * Tth(th, 1) / r;
}

static void mms_polar_source(double t, mtrx f, void *data)
{
    const mms_polar *m = (const mms_polar *)data;
    const polar_grid *pg = m->g;
    int i, j, d, nt = pg->ntheta;
    double g = gt(t), nu = 1.0 / m->Re;
    // T and its derivatives 0..4 on every column (the solution is separable)
    double *T = (double *)malloc((size_t)5 * nt * sizeof(double));
    for (j = 0; j < nt; j++)
        for (d = 0; d < 5; d++)
            T[d * nt + j] = Tth(j * pg->dtheta, d);
    for (i = 0; i < pg->nr; i++)
    {
        double r = pg->r[i], F = Fr(m, r, 0), F1 = Fr(m, r, 1), F2 = Fr(m, r, 2), F3 = Fr(m, r, 3), F4 = Fr(m, r, 4);
        double L1 = F2 + F1 / r, L1p = F3 + F2 / r - F1 / (r * r),
               L1pp = F4 + F3 / r - 2.0 * F2 / (r * r) + 2.0 * F1 / (r * r * r);
        double G = F / (r * r), Gp = F1 / (r * r) - 2.0 * F / (r * r * r),
               Gpp = F2 / (r * r) - 4.0 * F1 / (r * r * r) + 6.0 * F / (r * r * r * r);
        for (j = 0; j < nt; j++)
        {
            double T0 = T[j], T1 = T[nt + j], T2 = T[2 * nt + j], T3 = T[3 * nt + j], T4 = T[4 * nt + j];
            // Omega = -(L1 T + G T''), w = g Omega
            double W = -(L1 * T0 + G * T2), Wr = -(L1p * T0 + Gp * T2), Wt = -(L1 * T1 + G * T3);
            // lap(A(r) S(theta)) = (A'' + A'/r) S + A S''/r^2
            double lapW = -((L1pp + L1p / r) * T0 + L1 * T2 / (r * r) + (Gpp + Gp / r) * T2 + G * T4 / (r * r));
            double ur = g * F * T1 / r, ut = -g * F1 * T0;
            MAt(f, i, j) = gdot(t) * W + ur * g * Wr + ut / r * g * Wt - nu * g * lapW;
        }
    }
    free(T);
}

// Max error of w on the interior after time T from the exact state, and of
// the wall vorticity, on an ntheta grid with square cells
static void polar_mms_run_scheme(int ntheta, int order, int poisson_order, int closure, int time_scheme, double dt_frac,
                                 double *e_in, double *e_wall)
{
    const double R = 0.5, ri = 2.0, Re = 20.0, T = 0.1;
    polar_ops p = ops_make(ntheta, 0, R, ri, order);
    int nt = p.g.ntheta, nr = p.g.nr, i, j, steps;
    mms_polar m;
    mms_polar_init(&m, R, ri, Re);
    m.g = &p.g;
    {
        double *zero = (double *)calloc((size_t)nt, sizeof(double));
        polar_dirichlet_set(&p.g, zero, zero);
        free(zero);
    }

    smtrx dt2 = SDiff2_periodic(nt, order, p.g.dtheta), dx2 = SDiff2(nr, order, p.g.dxi);
    double dt = dt_frac * max_stable_dt(&dt2, &dx2, Re * R * R, time_scheme);
    freesm(dt2);
    freesm(dx2);
    steps = (int)ceil(T / dt);
    dt = T / steps;

    solver_config cfg = {0};
    cfg.nx = nt;
    cfg.ny = nr;
    cfg.dx = p.g.dtheta;
    cfg.dy = p.g.dxi;
    cfg.Re = Re;
    cfg.dt = dt;
    cfg.time_scheme = time_scheme;
    cfg.poisson_type = 3;
    cfg.poisson_order = poisson_order;
    cfg.wall_closure = closure;
    cfg.DX = &p.DX;
    cfg.DY = &p.DY;
    cfg.DX2 = &p.DX2;
    cfg.DY2 = &p.DY2;
    cfg.polar = &p.g;
    cfg.vorticity_source = mms_polar_source;
    cfg.source_data = &m;

    backend *b = backend_create(&cfg, 0);
    mtrx *u, *v, *w;
    backend_fields(b, &u, &v, &w);
    for (i = 0; i < nr; i++)
        for (j = 0; j < nt; j++)
            mms_polar_exact(&m, p.g.r[i], j * p.g.dtheta, 0.0, &MAt(*w, i, j), &MAt(*u, i, j), &MAt(*v, i, j));
    backend_set_fields(b, u, v, w);
    for (i = 0; i < steps; i++)
        backend_step(b);
    backend_fields(b, &u, &v, &w);
    *e_in = *e_wall = 0.0;
    for (i = 0; i < nr; i++)
        for (j = 0; j < nt; j++)
        {
            double we;
            mms_polar_exact(&m, p.g.r[i], j * p.g.dtheta, T, &we, NULL, NULL);
            if (i == 0)
                *e_wall = fmax(*e_wall, fabs(MAt(*w, i, j) - we));
            else if (i < nr - 1)
                *e_in = fmax(*e_in, fabs(MAt(*w, i, j) - we));
        }
    backend_free(b);
    ops_free(&p);
}

static void polar_mms_run(int ntheta, int order, int poisson_order, int closure, double *e_in, double *e_wall)
{
    polar_mms_run_scheme(ntheta, order, poisson_order, closure, 2, 0.5, e_in, e_wall);
}

static void test_polar_mms(void)
{
    double a, aw, b, bw, c, cw;
    char name[128];
    printf("\nPolar grid: manufactured solution of the full solver (RK4, r from 0.5 to 2, Re = 20, t = 0.1)\n");

    polar_mms_run(48, 2, 2, 0, &a, &aw);
    polar_mms_run(96, 2, 2, 0, &b, &bw);
    polar_mms_run(192, 2, 2, 0, &c, &cw);
    printf("  order 2, velocity closure: max |w - exact| interior %.3e %.3e %.3e (orders %.2f %.2f); wall %.3e %.3e "
           "%.3e\n",
           a, b, c, log2(a / b), log2(b / c), aw, bw, cw);
    check("  order 2 converges at second order in the interior (2 - observed)", 2.0 - log2(b / c), 0.25);
    check("  error on 192 x 46 / max |w| (about 30)", c / 30.0, 1E-3);

    polar_mms_run(48, 4, 4, 1, &a, &aw);
    polar_mms_run(96, 4, 4, 1, &b, &bw);
    polar_mms_run(192, 4, 4, 1, &c, &cw);
    double b4 = b;
    printf("  order 4, compact Poisson, Briley: interior %.3e %.3e %.3e (orders %.2f %.2f); wall %.3e %.3e %.3e\n", a,
           b, c, log2(a / b), log2(b / c), aw, bw, cw);
    snprintf(name, sizeof(name), "  order 4 + Briley: observed order %.2f (96 -> 192), 2.5 - order", log2(b / c));
    check(name, 2.5 - log2(b / c), 0.0);
    check("  and below the order-2 error by 5x on 192 x 46", c / 4.008e-3, 0.2);

    polar_mms_run(96, 6, 4, 1, &b, &bw);
    polar_mms_run(192, 6, 4, 1, &c, &cw);
    printf("  order 6, compact Poisson, Briley: interior %.3e %.3e (order %.2f)\n", b, c, log2(b / c));
    snprintf(name, sizeof(name), "  order 6 + Briley: observed order %.2f (96 -> 192), 2.5 - order", log2(b / c));
    check(name, 2.5 - log2(b / c), 0.0);

    // Euler at a step 20 times below its stability limit: the time error is
    // negligible, so it must reproduce the spatial error of RK4 on the same
    // grid (order 4, compact Poisson, Briley on 96 x 22; b above)
    double e1, e1w;
    polar_mms_run_scheme(96, 4, 4, 1, 1, 0.05, &e1, &e1w);
    printf("  Euler, 96 x 22, order 4: interior %.4e (RK4 %.4e)\n", e1, b4);
    check("  Euler on the polar grid matches RK4's spatial error (relative)", fabs(e1 - b4) / b4, 0.01);
}

// ---------------------------------------------------------------------------
// Diagnostics on fields whose answers are known
// ---------------------------------------------------------------------------

static void test_polar_diagnostics(void)
{
    polar_ops p = ops_make(128, 0, 0.5, 10.0, 6);
    int nt = p.g.ntheta, nr = p.g.nr, n = nt * nr, i, j;
    double *w = (double *)malloc((size_t)n * sizeof(double)), *v = (double *)malloc((size_t)n * sizeof(double));
    double *u = (double *)malloc((size_t)n * sizeof(double)), *pw = (double *)malloc((size_t)nt * sizeof(double));
    double Re = 40.0, nu = 1.0 / Re, alpha = -3.0, beta = 5.0, th0 = 0.9;

    printf("\nPolar grid: diagnostics\n");
    // w = (alpha + beta xi) sin(theta): dw/dxi = beta sin on the wall, so
    // C_Dp = nu beta pi, C_Df = -nu alpha pi, C_L = 0, p = -nu beta cos
    for (i = 0; i < nr; i++)
        for (j = 0; j < nt; j++)
            w[i * nt + j] = (alpha + beta * i * p.g.dxi) * p.g.sin_t[j];
    polar_forces f = polar_wall_forces(&p.g, &p.DY, w, Re, pw);
    check("C_D pressure part of a known wall field", fabs(f.cdp - nu * beta * PI), 1E-12);
    check("C_D friction part of a known wall field", fabs(f.cdf + nu * alpha * PI), 1E-12);
    check("C_L and mean dp/dtheta vanish", fabs(f.cl) + fabs(f.dpdt_mean), 1E-13);
    double pe = 0.0;
    for (j = 0; j < nt; j++)
        pe = fmax(pe, fabs(pw[j] + nu * beta * p.g.cos_t[j]));
    check("wall pressure integrated from dp/dtheta", pe, 1E-13);
    check("odd field: asymmetry", polar_asymmetry(&p.g, w), 1E-14);

    // Wall vorticity sin(theta) (cos theta - cos th0): separation at th0
    for (j = 0; j < nt; j++)
        w[j] = p.g.sin_t[j] * (p.g.cos_t[j] - cos(th0));
    check("separation angle (deg)", fabs(polar_separation_angle(&p.g, w) - (180.0 - th0 * 180.0 / PI)), 1E-4);
    for (j = 0; j < nt; j++)
        w[j] = -p.g.sin_t[j];
    check("no separation gives -1", fabs(polar_separation_angle(&p.g, w) + 1.0), 0.0);

    // u_r = r - 1.6 on the wake axis: the bubble ends at r = 1.6, L_w/D = 1.1
    for (i = 0; i < nr; i++)
        for (j = 0; j < nt; j++)
            v[i * nt + j] = p.g.r[i] - 1.6;
    check("wake length of a known profile", fabs(polar_wake_length(&p.g, v) - 1.1), 1E-4);
    // The sign change between the last two rows: the root is bracketed by
    // them (not extrapolated from the rows before)
    double rc = 0.5 * (p.g.r[nr - 2] + p.g.r[nr - 1]);
    for (i = 0; i < nr; i++)
        v[(size_t)i * nt] = p.g.r[i] - rc;
    check("wake length with the sign change on the last row",
          fabs(polar_wake_length(&p.g, v) - (rc - 0.5) / 1.0), 1E-3);
    for (i = 0; i < nr; i++)
        v[(size_t)i * nt] = -1.0;
    check("reverse flow out to the boundary: (r_inf - R) / D", fabs(polar_wake_length(&p.g, v) - 9.5), 0.0);

    // Interpolation of a smooth field: u_r = cos, u_theta = -sin is the
    // uniform stream u_x = 1, u_y = 0; w = x y / 10
    for (i = 0; i < nr; i++)
        for (j = 0; j < nt; j++)
        {
            double x = p.g.r[i] * p.g.cos_t[j], y = p.g.r[i] * p.g.sin_t[j];
            v[i * nt + j] = p.g.cos_t[j];
            u[i * nt + j] = -p.g.sin_t[j];
            w[i * nt + j] = x * y / 10.0;
        }
    double ux, uy, wv;
    int in = polar_interpolate(&p.g, u, v, w, 3.0, 3.0, &ux, &uy, &wv);
    check("interpolation at (3, 3): inside, uniform stream and x y / 10",
          in ? fabs(ux - 1.0) + fabs(uy) + fabs(wv - 0.9) : 1.0, 1E-5);
    check("interpolation outside the grid is refused", polar_interpolate(&p.g, u, v, w, 0.1, 0.1, &ux, &uy, &wv), 0.0);

    // psi with a minimum at a known place inside the upper half
    double xc = 1.2, yc = 0.3, a, b;
    for (i = 0; i < nr; i++)
        for (j = 0; j < nt; j++)
        {
            double x = p.g.r[i] * p.g.cos_t[j], y = p.g.r[i] * p.g.sin_t[j];
            w[i * nt + j] = -exp(-((x - xc) * (x - xc) + (y - yc) * (y - yc)) / 0.1);
        }
    int found = polar_vortex_centre(&p.g, w, &a, &b);
    check("vortex centre of a Gaussian at (1.2, 0.3): a/D, b/D", found ? fabs(a - 0.7) + fabs(b - 0.6) : 1.0, 5E-3);

    free(w);
    free(v);
    free(u);
    free(pw);
    ops_free(&p);
}

// The far-field condition: single-valued psi with the free stream's sign
// pattern, and freestream with cd = 0
static void test_polar_farfield(void)
{
    polar_grid g = polar_grid_make(256, 0, 0.5, 40.0);
    int j, nt = g.ntheta;
    double e = 0.0, odd = 0.0;

    printf("\nPolar grid: outer boundary conditions\n");
    polar_farfield_set(&g, 0.0, 40.0);
    for (j = 0; j < nt; j++)
        e = fmax(e, fabs(g.psi_out[j] - g.r_inf * g.sin_t[j]));
    check("far field with C_D = 0 is the free stream", e, 1E-12);
    polar_farfield_set(&g, 1.5, 40.0);
    for (j = 0; j < nt; j++)
        odd = fmax(odd, fabs(g.psi_out[j] + g.psi_out[(nt - j) % nt]));
    check("far-field psi is odd in theta (symmetric flow)", odd, 1E-12);
    // Outside the wake psi falls short of the free stream by C_D / 4 (half
    // the deficit flux C_D / 2), less the uniform outflow that returns it,
    // (C_D / 2) theta / (2 pi): at theta = pi/4, 0.375 - 0.09375
    double jump = g.r_inf * g.sin_t[nt / 8] - g.psi_out[nt / 8];
    check("deficit outside the wake, C_D/4 - (C_D/2) theta/(2 pi) (relative)", fabs(jump - 0.28125) / 0.28125, 0.02);
    polar_grid_free(&g);
}

// A polar configuration runs on the CPU even when the GPU is asked for
static void test_polar_backend(void)
{
    polar_ops p = ops_make(32, 0, 0.5, 4.0, 4);
    solver_config cfg = {0};
    printf("\nPolar grid: backend\n");
    cfg.nx = p.g.ntheta;
    cfg.ny = p.g.nr;
    cfg.dx = p.g.dtheta;
    cfg.dy = p.g.dxi;
    cfg.Re = 40.0;
    cfg.dt = 1E-4;
    cfg.time_scheme = 2;
    cfg.poisson_type = 3;
    cfg.poisson_order = 2;
    cfg.wall_closure = 1;
    cfg.DX = &p.DX;
    cfg.DY = &p.DY;
    cfg.DX2 = &p.DX2;
    cfg.DY2 = &p.DY2;
    cfg.polar = &p.g;
    backend *b = backend_create(&cfg, 1);
    check("polar configuration on the CPU backend", strcmp(backend_name(b), "CPU") != 0, 0.0);
    // A few steps of the impulsive start stay finite and symmetric
    int k, finite = 1;
    mtrx *w;
    for (k = 0; k < 20; k++)
        backend_step(b);
    backend_fields(b, NULL, NULL, &w);
    for (k = 0; k < w->m * w->n; k++)
        finite = finite && isfinite(w->M[k]);
    check("impulsive start: finite", !finite, 0.0);
    double wmax = 0.0;
    for (k = 0; k < w->m * w->n; k++)
        wmax = fmax(wmax, fabs(w->M[k]));
    check("impulsive start: symmetric (asymmetry / max |w|)", polar_asymmetry(&p.g, w->M) / wmax, 1E-12);
    backend_free(b);
    ops_free(&p);
}

// A polar configuration of the solver: RK4, compact Poisson, Briley
static solver_config polar_cfg(polar_ops *p)
{
    solver_config cfg = {0};
    cfg.nx = p->g.ntheta;
    cfg.ny = p->g.nr;
    cfg.dx = p->g.dtheta;
    cfg.dy = p->g.dxi;
    cfg.Re = 40.0;
    cfg.dt = 1E-3;
    cfg.time_scheme = 2;
    cfg.poisson_type = 3;
    cfg.poisson_order = 4;
    cfg.wall_closure = 1;
    cfg.DX = &p->DX;
    cfg.DY = &p->DY;
    cfg.DX2 = &p->DX2;
    cfg.DY2 = &p->DY2;
    cfg.polar = &p->g;
    return cfg;
}

// Every option the polar right-hand side does not implement is refused, not
// ignored; vorticity_source is honoured
static void test_polar_config(void)
{
    polar_ops p = ops_make(32, 0, 0.5, 4.0, 4);
    solver_config ok = polar_cfg(&p), c;
    int cases = 0, refused = 0;
    char name[96];

    printf("\nPolar grid: configuration\n");
    check("a valid polar configuration is accepted", polar_config_error(&ok) != NULL, 0.0);
    c = ok;
    c.vorticity_source = mms_polar_source;
    check("  with a vorticity source", polar_config_error(&c) != NULL, 0.0);
    c = ok;
    c.polar = NULL;
    check("a Cartesian configuration has no polar error", polar_config_error(&c) != NULL, 0.0);
#define REFUSED(stmt)                              \
    do                                             \
    {                                              \
        c = ok;                                    \
        stmt;                                      \
        cases++;                                   \
        refused += polar_config_error(&c) != NULL; \
    } while (0)
    REFUSED(c.forcing.drag = 0.1);
    REFUSED(c.forcing.kolmogorov_amp = 1.0);
    REFUSED(c.forcing.random_rate = 0.1);
    REFUSED(c.forcing.hyperviscosity = 1E-6);
    REFUSED(c.forcing.hypodrag = 1.0);
    REFUSED(c.fourier = 1);
    REFUSED(c.dealias = 1);
    REFUSED(c.advection = 1);
    REFUSED(c.periodic = 1);
    REFUSED(c.poisson_type = 2);
    REFUSED(c.DXv = &p.DX);
    REFUSED(c.DYv = &p.DY);
    REFUSED(c.nx = c.nx - 1);
    REFUSED(c.ny = c.ny + 1);
#undef REFUSED
    snprintf(name, sizeof(name), "each of %d unsupported options refused (count not refused)", cases);
    check(name, cases - refused, 0.0);
    ops_free(&p);
}

// The area means of the polar branch of compute_integrals() against fields
// with known integrals between R = 1/2 and r_inf = 4, all three second order
// (the trapezoidal rule in xi, for the integrals and the area alike). E and
// Z: u_theta = r - R, u_r = 0 (no slip on the cylinder), w = r. P: w = xi,
// which every row of DY (the one-sided ones too) differentiates exactly, so
// |grad w|^2 = 1/r^2 and P = (1/2) 2 pi ln(r_inf/R) / area.
static double polar_integrals_error(int ntheta, double *eE, double *eZ, double *eP)
{
    const double R = 0.5, ri = 4.0;
    polar_ops p = ops_make(ntheta, 0, R, ri, 6);
    solver_config cfg = polar_cfg(&p);
    int i, j, nt = p.g.ntheta, nr = p.g.nr;
    mtrx u = initm(nr, nt), v = initm(nr, nt), w = initm(nr, nt);
    double *wx = (double *)calloc((size_t)nt * nr, sizeof(double)), *wy = (double *)calloc((size_t)nt * nr, sizeof(double));
    for (i = 0; i < nr; i++)
        for (j = 0; j < nt; j++)
        {
            MAt(u, i, j) = p.g.r[i] - R;
            MAt(w, i, j) = p.g.r[i];
        }
    flow_integrals fi = compute_integrals(&cfg, u, v, w, wx, wy);
    for (i = 0; i < nr; i++)
        for (j = 0; j < nt; j++)
            MAt(w, i, j) = i * p.g.dxi;
    flow_integrals fp = compute_integrals(&cfg, u, v, w, wx, wy);
    // Area pi (ri^2 - R^2); E = (1/2) 2 pi int (r - R)^2 r dr / area,
    // Z = (1/2) 2 pi int r^3 dr / area, P = 1/2
    double area = PI * (ri * ri - R * R), P = PI * log(ri / R) / area;
#define FE(r) ((r) * (r) * (r) * (r) / 4.0 - 2.0 * R * (r) * (r) * (r) / 3.0 + R * R * (r) * (r) / 2.0)
    double E = PI * (FE(ri) - FE(R)) / area, Z = PI * (ri * ri * ri * ri - R * R * R * R) / 4.0 / area;
#undef FE
    *eE = fabs(fi.E - E) / E;
    *eZ = fabs(fi.Z - Z) / Z;
    *eP = fabs(fp.P - P) / P;
    freem(&u);
    freem(&v);
    freem(&w);
    free(wx);
    free(wy);
    ops_free(&p);
    return fmax(*eE, fmax(*eZ, *eP));
}

static void test_polar_integrals(void)
{
    double a[3], b[3];
    char name[128];
    printf("\nPolar grid: energy, enstrophy, palinstrophy (area means)\n");

    // The mean of a constant is exact: the area is summed by the same rule.
    // u_theta = 1 with the cylinder moving at 1 (its row takes wall_ut), w = 1
    {
        polar_ops p = ops_make(48, 0, 0.5, 4.0, 6);
        solver_config cfg = polar_cfg(&p);
        int k, n = p.g.ntheta * p.g.nr;
        mtrx u = initm(p.g.nr, p.g.ntheta), v = initm(p.g.nr, p.g.ntheta), w = initm(p.g.nr, p.g.ntheta);
        double *wx = (double *)calloc((size_t)n, sizeof(double)), *wy = (double *)calloc((size_t)n, sizeof(double));
        p.g.wall_ut = 1.0;
        for (k = 0; k < n; k++)
            u.M[k] = w.M[k] = 1.0;
        for (k = 0; k < p.g.ntheta; k++)
            u.M[k] = 0.0; // replaced by wall_ut
        flow_integrals fi = compute_integrals(&cfg, u, v, w, wx, wy);
        check("constants: E = Z = 1/2, P = 0 to round-off", fabs(fi.E - 0.5) + fabs(fi.Z - 0.5) + fabs(fi.P), 1E-14);
        freem(&u);
        freem(&v);
        freem(&w);
        free(wx);
        free(wy);
        ops_free(&p);
    }

    polar_integrals_error(64, &a[0], &a[1], &a[2]);
    polar_integrals_error(128, &b[0], &b[1], &b[2]);
    printf("  relative errors on 64 / 128 points: E %.2e / %.2e, Z %.2e / %.2e, P %.2e / %.2e\n", a[0], b[0], a[1], b[1],
           a[2], b[2]);
    check("E, Z, P within 5e-3 of the exact means on 128 points", fmax(b[0], fmax(b[1], b[2])), 5E-3);
    snprintf(name, sizeof(name), "E, Z, P converge at second order: %.2f, %.2f, %.2f; 1.9 - min", log2(a[0] / b[0]),
             log2(a[1] / b[1]), log2(a[2] / b[2]));
    check(name, 1.9 - fmin(log2(a[0] / b[0]), fmin(log2(a[1] / b[1]), log2(a[2] / b[2]))), 0.0);
}

// A change of the outer boundary values between steps reaches the next
// step: a solver whose boundary changes after 20 steps takes the same next
// step, bitwise, as a fresh solver started from its state with the new
// boundary (without the version check in psi_is_current() the first RK4
// stage would reuse the stream function of the old boundary). And the
// solver's stream function (backend_psi()) is the one that gave its velocity.
static void test_polar_boundary_change(void)
{
    polar_ops p = ops_make(32, 0, 0.5, 4.0, 4), q = ops_make(32, 0, 0.5, 4.0, 4);
    solver_config cp = polar_cfg(&p), cq = polar_cfg(&q);
    int i, k, nt = p.g.ntheta, nr = p.g.nr, n = nt * nr;
    double *psi_b = (double *)malloc((size_t)nt * sizeof(double)), *w_b = (double *)calloc((size_t)nt, sizeof(double));
    mtrx u0 = initm(nr, nt), v0 = initm(nr, nt), w0 = initm(nr, nt), *u, *v, *w, *wq;
    double *us = (double *)malloc((size_t)n * sizeof(double)), *vs = (double *)malloc((size_t)n * sizeof(double));

    printf("\nPolar grid: outer boundary changed between steps\n");
    for (k = 0; k < nt; k++)
        psi_b[k] = 1.2 * p.g.r_inf * p.g.sin_t[k] + 0.3 * p.g.cos_t[(2 * k) % nt];
    backend *a = backend_create(&cp, 0);
    for (k = 0; k < 20; k++)
        backend_step(a);
    backend_get_fields(a, &u0, &v0, &w0);
    polar_dirichlet_set(&p.g, psi_b, w_b);
    polar_dirichlet_set(&q.g, psi_b, w_b);
    backend_step(a);
    backend *b = backend_create(&cq, 0);
    backend_set_fields(b, &u0, &v0, &w0);
    backend_step(b);
    backend_fields(a, &u, &v, &w);
    backend_fields(b, NULL, NULL, &wq);
    double d = 0.0;
    for (k = 0; k < n; k++)
        d = fmax(d, fabs(w->M[k] - wq->M[k]));
    check("the step after the change equals a fresh solver's (w, bitwise)", d, 0.0);

    // u, v from the solver's own psi, on the interior rows (the cylinder row
    // of u, v holds the wall velocity only after polar_apply_bc())
    const mtrx *psi = backend_psi(a);
    polar_velocity(&p.g, &p.DX, &p.DY, psi->M, us, vs);
    d = 0.0;
    for (i = 1; i < nr; i++)
        for (k = 0; k < nt; k++)
            d = fmax(d, fabs(us[i * nt + k] - MAt(*u, i, k)) + fabs(vs[i * nt + k] - MAt(*v, i, k)));
    for (k = 0; k < nt; k++)
        d = fmax(d, fabs(MAt(*psi, nr - 1, k) - psi_b[k]));
    check("backend_psi() gives the solver's u, v and the new boundary (bitwise)", d, 0.0);

    backend_free(a);
    backend_free(b);
    freem(&u0);
    freem(&v0);
    freem(&w0);
    free(psi_b);
    free(w_b);
    free(us);
    free(vs);
    ops_free(&p);
    ops_free(&q);
}

void test_polar(check_fn fn)
{
    check = fn;
    test_polar_poisson();
    test_polar_farfield();
    test_polar_diagnostics();
    test_polar_backend();
    test_polar_mms();
    test_polar_config();
    test_polar_integrals();
    test_polar_boundary_change();
}
