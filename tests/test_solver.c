// Test suite: make test            (CPU solver)
//             make CUDA=1 test     (CPU solver + GPU backend checked against it)
//
// Exit status: 0 = every check ran and passed, 1 = a check failed,
// 77 = nothing failed but the GPU tests could not run (no usable device).

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include "linearalg.h"
#include "finitediff.h"
#include "poisson.h"
#include "fluiddyn.h"
#include "threads.h"
#include "mms.h"
#include "diagnostics.h"
#include "backend.h"
#include "utils.h"
#ifdef USE_CUDA
#include "cudasolver.h"
#endif

#define EXIT_SKIPPED 77 // conventional "skipped" status: neither a pass nor a failure

static int n_checks = 0;
static int n_failed = 0;
static int n_skipped = 0;       // GPU tests that could not run
static int n_skipped_other = 0; // other checks this machine cannot run

// Pass when value <= limit (a NaN value fails)
static void check(const char *name, double value, double limit)
{
    int ok = (value <= limit);
    n_checks++;
    if (!ok) n_failed++;
    printf("  [%s] %-52s %.3E (limit %.1E)\n", ok ? " ok " : "FAIL", name, value, limit);
}

// ---------------------------------------------------------------------------
// Problem setup — same construction as main.c
// ---------------------------------------------------------------------------

typedef struct
{
    int nx, ny;
    smtrx DX, DY, DX2, DY2;
    smtrx DXv, DYv;    // velocity operators (velocity_order 4), else empty
    solver_config cfg; // shared by the CPU and GPU solvers
    rk4_ctx ctx;       // CPU workspace
    mtrx u, v, w;
} problem;

// A lid-driven cavity on a unit square with nx x ny nodes (bc == NULL), or
// the given wall velocities. Fields are ny rows (y) of nx values (x), as in
// main.c. The configuration is complete, operators included, before the
// workspace takes its copy, and the configuration points into *p, so p must
// stay where it is until problem_free().
// problem_init_ext() also sets the derivative order, Re, the domain size,
// periodic boundaries (then dx = Lx/nx, and bc is unused), the start time, an
// optional vorticity source, the Poisson operator order (FFT solver), the
// wall-vorticity closure and the order of the velocity rows next to the walls.
static void problem_init_ext(problem *p, int nx, int ny, double Lx, double Ly, int periodic, int order,
                             double Re, int time_scheme, int poisson_type, double dt, double poisson_tol,
                             const wall_bc *bc, double t0, void (*vorticity_source)(double, mtrx, void *),
                             void *source_data, int poisson_order, int wall_closure, int velocity_order)
{
    double dx = Lx / (periodic ? nx : nx - 1), dy = Ly / (periodic ? ny : ny - 1);
    wall_bc lid = {{0., 0., 0., 1.}, {0., 0., 0., 0.}};

    p->nx = nx;
    p->ny = ny;

    smtrx sd_x = periodic ? SDiff1_periodic(nx, order, dx) : SDiff1(nx, order, dx);
    smtrx sd_y = periodic ? SDiff1_periodic(ny, order, dy) : SDiff1(ny, order, dy);
    smtrx sd_x2 = periodic ? SDiff2_periodic(nx, order, dx) : SDiff2(nx, order, dx);
    smtrx sd_y2 = periodic ? SDiff2_periodic(ny, order, dy) : SDiff2(ny, order, dy);
    smtrx sIx = seye(nx);
    smtrx sIy = seye(ny);

    p->DX = skronecker(sIy, sd_x);
    p->DY = skronecker(sd_y, sIx);
    p->DX2 = skronecker(sIy, sd_x2);
    p->DY2 = skronecker(sd_y2, sIx);

    freesm(sd_x);
    freesm(sd_y);
    freesm(sd_x2);
    freesm(sd_y2);
    freesm(sIx);
    freesm(sIy);

    p->cfg.nx = nx;
    p->cfg.ny = ny;
    p->cfg.dx = dx;
    p->cfg.dy = dy;
    p->cfg.Re = Re;
    p->cfg.dt = dt;
    p->cfg.time_scheme = time_scheme;
    p->cfg.poisson_type = poisson_type;
    p->cfg.poisson_max_it = 200000;
    p->cfg.poisson_order = poisson_order;
    p->cfg.wall_closure = wall_closure;
    p->cfg.poisson_tol = poisson_tol;
    p->cfg.beta = sor_beta(nx, ny, dx, dy);
    p->cfg.periodic = periodic;
    p->cfg.advection = 0;
    p->cfg.bc = bc ? *bc : lid;
    p->cfg.DX = &p->DX;
    p->cfg.DY = &p->DY;
    p->cfg.DX2 = &p->DX2;
    p->cfg.DY2 = &p->DY2;
    p->DXv = (smtrx){0};
    p->DYv = (smtrx){0};
    if (velocity_order == 4)
    {
        smtrx vx = SDiff1_wall4(nx, order, dx), vy = SDiff1_wall4(ny, order, dy), Ix = seye(nx), Iy = seye(ny);
        p->DXv = skronecker(Iy, vx);
        p->DYv = skronecker(vy, Ix);
        freesm(vx);
        freesm(vy);
        freesm(Ix);
        freesm(Iy);
    }
    p->cfg.DXv = velocity_order == 4 ? &p->DXv : NULL;
    p->cfg.DYv = velocity_order == 4 ? &p->DYv : NULL;
    p->cfg.t0 = t0;
    p->cfg.vorticity_source = vorticity_source;
    p->cfg.source_data = source_data;
    p->cfg.forcing = (forcing_config){0};
    p->ctx = rk4_alloc(&p->cfg);

    p->u = initm(ny, nx);
    p->v = initm(ny, nx);
    p->w = initm(ny, nx);
}

static void problem_init(problem *p, int nx, int ny, int time_scheme, int poisson_type, double dt,
                         double poisson_tol, const wall_bc *bc)
{
    problem_init_ext(p, nx, ny, 1.0, 1.0, 0, 6, 100., time_scheme, poisson_type, dt, poisson_tol, bc, 0.0, NULL, NULL, 2, 0, 2);
}

static void problem_free(problem *p)
{
    freem(&p->u);
    freem(&p->v);
    freem(&p->w);
    freesm(p->DX);
    freesm(p->DY);
    freesm(p->DX2);
    freesm(p->DY2);
    if (p->cfg.DXv)
    {
        freesm(p->DXv);
        freesm(p->DYv);
    }
    rk4_free(&p->ctx);
}

// Deterministic values in [-1, 1]
static void fill_pseudo_random(double *x, int n, unsigned int seed)
{
    int i;
    for (i = 0; i < n; i++)
    {
        seed = seed * 1664525u + 1013904223u;
        x[i] = (double)(seed >> 8) / (double)(1u << 23) - 1.0;
    }
}

static double max_abs(const double *x, int n)
{
    int i;
    double m = 0.0;
    for (i = 0; i < n; i++)
    {
        if (isnan(x[i])) return NAN;
        if (fabs(x[i]) > m) m = fabs(x[i]);
    }
    return m;
}

// max |a - b| / max |b|
static double rel_diff(const double *a, const double *b, int n)
{
    int i;
    double m = 0.0;
    for (i = 0; i < n; i++)
    {
        if (isnan(a[i]) || isnan(b[i])) return NAN;
        if (fabs(a[i] - b[i]) > m) m = fabs(a[i] - b[i]);
    }
    return m / max_abs(b, n);
}

// Max and min of du/dx + dv/dy on the CPU
static void continuity_range(problem *p, double *cmax, double *cmin)
{
    int k, n = p->nx * p->ny;
    double *dudx = (double *)malloc(n * sizeof(double));
    double *dvdy = (double *)malloc(n * sizeof(double));

    spmv(p->DX, p->u.M, dudx);
    spmv(p->DY, p->v.M, dvdy);
    *cmax = -__DBL_MAX__;
    *cmin = __DBL_MAX__;
    for (k = 0; k < n; k++)
    {
        double c = dudx[k] + dvdy[k];
        if (c > *cmax) *cmax = c;
        if (c < *cmin) *cmin = c;
    }
    free(dudx);
    free(dvdy);
}

// ---------------------------------------------------------------------------
// Unit tests of the building blocks
// ---------------------------------------------------------------------------

static mtrx csr_to_dense(smtrx S)
{
    int i, k;
    mtrx D = initm(S.m, S.n);
    for (i = 0; i < S.m; i++)
        for (k = S.row_ptr[i]; k < S.row_ptr[i + 1]; k++)
            MAt(D, i, S.col_idx[k]) = S.values[k];
    return D;
}

// Sparse operations against their dense counterparts
static void test_linearalg(void)
{
    int i, n = 12;
    double dx = 1.0 / (n - 1);
    mtrx x = initm(n, 1), y_dense, d, k_dense, k_sparse, e_dense, e_sparse, a, b;
    double y_sparse[12];
    smtrx s, s2, s1, se;

    printf("Unit: linear algebra\n");

    // y = A x, dense and CSR
    s = SDiff1(n, 6, dx);
    d = csr_to_dense(s);
    fill_pseudo_random(x.M, n, 3u);
    y_dense = mtrxmul(d, x);
    spmv(s, x.M, y_sparse);
    check("spmv vs dense matrix product", rel_diff(y_sparse, y_dense.M, n), 1E-13);

    // Kronecker products
    s2 = SDiff2(5, 2, 0.25);
    s1 = SDiff1(4, 2, 1.0 / 3.0);
    a = csr_to_dense(s2);
    b = csr_to_dense(s1);
    k_dense = kronecker(a, b);
    smtrx ks = skronecker(s2, s1);
    k_sparse = csr_to_dense(ks);
    check("skronecker vs kronecker", rel_diff(k_sparse.M, k_dense.M, 20 * 20), 0.0);

    // Identity
    se = seye(7);
    e_sparse = csr_to_dense(se);
    e_dense = eye(7);
    check("seye vs eye", rel_diff(e_sparse.M, e_dense.M, 49), 0.0);

    // Element-wise helpers
    mtrx m = initm(3, 4), m2 = initm(3, 4);
    double flat[12];
    for (i = 0; i < 12; i++)
        m.M[i] = (i % 5) - 2.5 * (i == 7) + 3.0 * (i == 2);
    flatten(m, flat, 3, 4);
    unflatten(flat, m2, 3, 4);
    check("flatten / unflatten round trip", rel_diff(m2.M, m.M, 12), 0.0);
    check("maxel", fabs(maxel(m) - 5.0), 0.0);
    check("minel", fabs(minel(m) - (-0.5)), 0.0);
    negcpy(m2, m);
    for (i = 0; i < 12; i++)
        m2.M[i] += m.M[i];
    check("negcpy", max_abs(m2.M, 12), 0.0);

    freem(&x);
    freem(&y_dense);
    freem(&d);
    freem(&a);
    freem(&b);
    freem(&k_dense);
    freem(&k_sparse);
    freem(&e_dense);
    freem(&e_sparse);
    freem(&m);
    freem(&m2);
    freesm(s);
    freesm(s1);
    freesm(s2);
    freesm(ks);
    freesm(se);
}

// ---------------------------------------------------------------------------
// Finite-difference operators
// ---------------------------------------------------------------------------

// Every row of SDiff1/SDiff2, boundary rows included, must differentiate
// polynomials exactly up to the degree its stencil is built for. A row b nodes
// from the nearer end uses, for the first derivative, a one-sided first-order
// stencil (b = 0) or a centered stencil of order min(order, 2b); for the second
// derivative a one-sided stencil exact for cubics (b = 0) or a centered one of
// order min(order, 2b), which is exact one degree higher.
static void test_finitediff_exactness(void)
{
    int o, op, r, p, n = 16;
    int orders[3] = {2, 4, 6};
    double dx = 0.1;
    char name[96];

    printf("CPU: finite-difference operators are exact on polynomials\n");
    for (o = 0; o < 3; o++)
        for (op = 1; op <= 2; op++)
        {
            int order = orders[o];
            smtrx S = op == 1 ? SDiff1(n, order, dx) : SDiff2(n, order, dx);
            double worst = 0.0;
            for (r = 0; r < n; r++)
            {
                int b = r < n - 1 - r ? r : n - 1 - r;
                int deg = (b == 0) ? (op == 1 ? 1 : 3)
                                   : (2 * b < order ? 2 * b : order) + (op == 2);
                for (p = 0; p <= deg; p++)
                {
                    // f(x) = (x - x_r)^p, whose derivative at x_r is known
                    double sum = 0.0, expect = (op == 1) ? (p == 1) : 2.0 * (p == 2);
                    int k;
                    for (k = S.row_ptr[r]; k < S.row_ptr[r + 1]; k++)
                        sum += S.values[k] * pow((S.col_idx[k] - r) * dx, p);
                    if (isnan(sum))
                        worst = NAN;
                    else if (!isnan(worst) && fabs(sum - expect) > worst)
                        worst = fabs(sum - expect);
                }
            }
            snprintf(name, sizeof(name), "order %d, %s derivative: every row exact", order,
                     op == 1 ? "first" : "second");
            check(name, isnan(worst) ? INFINITY : worst, 1E-8);
            freesm(S);
        }
}

// Largest |difference| between row r of S and the expected (columns, values);
// infinite if the row has a different set of non-zeros
static double row_diff(smtrx S, int r, int count, const int *cols, const double *vals)
{
    int k;
    double d = 0.0;
    if (S.row_ptr[r + 1] - S.row_ptr[r] != count) return INFINITY;
    for (k = 0; k < count; k++)
    {
        if (S.col_idx[S.row_ptr[r] + k] != cols[k]) return INFINITY;
        d = fmax(d, fabs(S.values[S.row_ptr[r] + k] - vals[k]));
    }
    return d;
}

// A few rows written out, from each block of the 6th-order operators: the
// first rows, an interior row, and the last rows that are copied from the
// first ones. Checks the column order, that zero coefficients are not stored,
// and the mirrored rows (built by reading back earlier entries).
static void test_finitediff_rows(void)
{
    int n = 10;
    double dx = 0.5, h = 1.0 / dx, h2 = 1.0 / (dx * dx);
    smtrx d1 = SDiff1(n, 6, dx), d2 = SDiff2(n, 6, dx);
    double worst = 0.0;

    printf("CPU: finite-difference operators, rows written out\n");
    {
        int c0[] = {0, 1};
        double v0[] = {-h, h};
        int c1[] = {0, 2};
        double v1[] = {-0.5 * h, 0.5 * h};
        int c2[] = {0, 1, 3, 4};
        double v2[] = {1.0 / 12.0 * h, -2.0 / 3.0 * h, 2.0 / 3.0 * h, -1.0 / 12.0 * h};
        int c5[] = {2, 3, 4, 6, 7, 8};
        double v5[] = {-1.0 / 60.0 * h, 3.0 / 20.0 * h, -3.0 / 4.0 * h, 3.0 / 4.0 * h, -3.0 / 20.0 * h, 1.0 / 60.0 * h};
        int c8[] = {7, 9};
        double v8[] = {-0.5 * h, 0.5 * h};
        int c9[] = {8, 9};
        double v9[] = {-h, h};
        worst = fmax(worst, row_diff(d1, 0, 2, c0, v0));
        worst = fmax(worst, row_diff(d1, 1, 2, c1, v1));
        worst = fmax(worst, row_diff(d1, 2, 4, c2, v2));
        worst = fmax(worst, row_diff(d1, 5, 6, c5, v5));
        worst = fmax(worst, row_diff(d1, 8, 2, c8, v8));
        worst = fmax(worst, row_diff(d1, 9, 2, c9, v9));
    }
    check("first derivative, order 6: rows 0, 1, 2, 5, 8, 9", worst, 1E-14);

    worst = 0.0;
    {
        int c0[] = {0, 1, 2, 3};
        double v0[] = {2 * h2, -5 * h2, 4 * h2, -1 * h2};
        int c1[] = {0, 1, 2};
        double v1[] = {h2, -2 * h2, h2};
        int c5[] = {2, 3, 4, 5, 6, 7, 8};
        double v5[] = {1.0 / 90.0 * h2, -3.0 / 20.0 * h2, 3.0 / 2.0 * h2, -49.0 / 18.0 * h2,
                       3.0 / 2.0 * h2, -3.0 / 20.0 * h2, 1.0 / 90.0 * h2};
        int c9[] = {6, 7, 8, 9};
        double v9[] = {-1 * h2, 4 * h2, -5 * h2, 2 * h2};
        worst = fmax(worst, row_diff(d2, 0, 4, c0, v0));
        worst = fmax(worst, row_diff(d2, 1, 3, c1, v1));
        worst = fmax(worst, row_diff(d2, 5, 7, c5, v5));
        worst = fmax(worst, row_diff(d2, 9, 4, c9, v9));
    }
    check("second derivative, order 6: rows 0, 1, 5, 9", worst, 1E-12);

    freesm(d1);
    freesm(d2);
}

// Run the cavity from rest to time T and return w
static mtrx run_to(int n, int scheme, double dt, double T)
{
    int t, steps = (int)floor(T / dt + 0.5);
    problem p;
    mtrx w = initm(n, n);
    problem_init(&p, n, n, scheme, 3, dt, 1E-3, NULL);
    for (t = 0; t < steps; t++)
        step(p.w, p.u, p.v, &p.ctx);
    mtrxcpy(w, p.w);
    problem_free(&p);
    return w;
}

// Largest |a - b| over all nodes, walls included
static double field_diff(mtrx a, mtrx b)
{
    int i;
    double m = 0.0;
    for (i = 0; i < a.m * a.n; i++)
    {
        double d = fabs(a.M[i] - b.M[i]);
        if (isnan(d)) return NAN;
        if (d > m) m = d;
    }
    return m;
}

// Observed order of the time integration: halve dt and compare the errors in
// w against a run of the same scheme with a much smaller step. Euler is first
// order, RK4 fourth order. The whole field is compared: step() returns the
// wall vorticity of the new velocity, so the walls converge with the
// interior. RK4 is fourth order only because the wall vorticity is set again
// at every stage; set once per step it was first order with about Euler's
// error.
static void test_temporal_order(void)
{
    int n = 17, s, k;
    double T = 0.2, order, err_at[3];
    char name[96];

    printf("Unit: observed order of the time integration, %dx%d grid, t = %g\n", n, n, T);
    for (s = 1; s <= 2; s++)
    {
        double err[2];
        mtrx ref = run_to(n, s, T / 4096, T);
        for (k = 0; k < 2; k++)
        {
            mtrx w = run_to(n, s, T / (64 << k), T);
            err[k] = field_diff(w, ref);
            freem(&w);
        }
        freem(&ref);
        err_at[s] = err[0];
        order = log2(err[0] / err[1]);
        if (s == 1)
        {
            snprintf(name, sizeof(name), "Euler: observed order %.2f, |order - 1|", order);
            check(name, isnan(order) ? INFINITY : fabs(order - 1.0), 0.15);
        }
        else
        {
            snprintf(name, sizeof(name), "RK4: observed order %.2f, |order - 4|", order);
            check(name, isnan(order) ? INFINITY : fabs(order - 4.0), 0.3);
        }
    }
    snprintf(name, sizeof(name), "RK4 / Euler error at dt = T/64 (%.1e / %.1e)", err_at[2], err_at[1]);
    check(name, isnan(err_at[2] / err_at[1]) ? INFINITY : err_at[2] / err_at[1], 1E-3);
}

// The manufactured solution of mms.h on an n x n unit square, run with its
// vorticity source from the exact state at t = 0 to T; returns w
static mtrx run_with_source(int n, int scheme, double dt, double T, mms_case *c)
{
    int t, steps = (int)floor(T / dt + 0.5);
    wall_bc walls = {{0., 0., 0., 0.}, {0., 0., 0., 0.}};
    mtrx w = initm(n, n);
    problem p;

    problem_init_ext(&p, n, n, 1.0, 1.0, 0, 6, c->Re, scheme, 3, dt, 1E-3, &walls, 0.0, mms_source, c, 2, 0, 2);
    mms_exact(c, 0.0, &p.w, &p.u, &p.v, NULL);
    for (t = 0; t < steps; t++)
        step(p.w, p.u, p.v, &p.ctx);
    mtrxcpy(w, p.w);
    problem_free(&p);
    return w;
}

// With a time-dependent vorticity source the schemes must keep their order, which
// they do only if each RK4 stage evaluates the source at its own time
static void test_source_temporal_order(void)
{
    int n = 17, s, k;
    double T = 0.25;
    mms_case c = {1.0, 1.0, 100.0, 1.0 / (n - 1), 1.0 / (n - 1), 0};
    char name[96];

    printf("Unit: observed order in time with a vorticity source, %dx%d grid, t = %g\n", n, n, T);
    for (s = 1; s <= 2; s++)
    {
        int coarse = s == 1 ? 64 : 16;
        double err[2], order;
        mtrx ref = run_with_source(n, s, T / (64 * coarse), T, &c);
        for (k = 0; k < 2; k++)
        {
            mtrx w = run_with_source(n, s, T / (coarse << k), T, &c);
            err[k] = field_diff(w, ref);
            freem(&w);
        }
        freem(&ref);
        order = log2(err[0] / err[1]);
        snprintf(name, sizeof(name), "%s with a source: observed order %.2f, |order - %d|",
                 s == 1 ? "Euler" : "RK4", order, s == 1 ? 1 : 4);
        check(name, isnan(order) ? INFINITY : fabs(order - (s == 1 ? 1.0 : 4.0)), s == 1 ? 0.15 : 0.3);
    }
}

// Spatial convergence against the manufactured solution (make convergence
// runs the full study). psi, u and v converge at second order, set by the
// five-point Poisson operator; w converges at least as fast on these grids.
// The bounds on the errors themselves catch a wrong scale that would leave
// the orders intact.
static void test_mms_spatial_order(void)
{
    int k, ns[3] = {17, 33, 65};
    mms_errors e[3];
    char name[96];

    printf("Unit: spatial order against a manufactured solution, RK4 + FFT, order 6\n");
    for (k = 0; k < 3; k++)
        e[k] = mms_run(ns[k], ns[k], 1.0, 1.0, 100.0, 6, 2, 3, 2.5E-3, 0.0, 0.25);

#define ORDER(field) log2(e[1].field.max / e[2].field.max)
    snprintf(name, sizeof(name), "psi: observed order %.2f (33 -> 65), |order - 2|", ORDER(psi));
    check(name, isnan(ORDER(psi)) ? INFINITY : fabs(ORDER(psi) - 2.0), 0.15);
    snprintf(name, sizeof(name), "u: observed order %.2f, |order - 2|", ORDER(u));
    check(name, isnan(ORDER(u)) ? INFINITY : fabs(ORDER(u) - 2.0), 0.15);
    snprintf(name, sizeof(name), "v: observed order %.2f, |order - 2|", ORDER(v));
    check(name, isnan(ORDER(v)) ? INFINITY : fabs(ORDER(v) - 2.0), 0.15);
    snprintf(name, sizeof(name), "w, walls included: observed order %.2f, 2 - order", ORDER(w));
    check(name, isnan(ORDER(w)) ? INFINITY : 2.0 - ORDER(w), 0.1);
#undef ORDER
    check("psi max error at 65x65", e[2].psi.max, 3E-4);
    check("w max error at 65x65", e[2].w.max, 5E-3);

    // The ablation study (make ablation) runs a copy of the RK4 step. With
    // nothing replaced it must compute exactly what step() computes.
    mms_errors copy = mms_run_ablated(33, 33, 1.0, 1.0, 100.0, 6, 0, 2.5E-3, 0.0, 0.25);
    check("ablation copy of step(), nothing replaced: same errors, bitwise",
          (copy.psi.max != e[1].psi.max) + (copy.u.max != e[1].u.max) + (copy.v.rms != e[1].v.rms) +
              (copy.w.max != e[1].w.max) + (copy.w.rms != e[1].w.rms) + (copy.w_interior.rms != e[1].w_interior.rms),
          0.0);

    // A run that starts at t0 = 0.4 must evaluate the source at the times of
    // that run, so its error is no larger than from t = 0
    mms_errors late = mms_run(33, 33, 1.0, 1.0, 100.0, 6, 2, 3, 2.5E-3, 0.4, 0.25);
    snprintf(name, sizeof(name), "start at t0 = 0.4: psi error %.2e vs %.2e from t = 0", late.psi.max, e[1].psi.max);
    check(name, late.psi.max / e[1].psi.max, 1.5);
}

// The ablation (make ablation, issue #25): with derivative order 6, the
// Poisson operator alone and the wall-vorticity closure alone each hold the
// order of w to two; with both replaced by the exact solution it rises to
// about four
static void test_order_ablation(void)
{
    static const struct
    {
        const char *name;
        int flags;
        double target, tol; // order between 33 and 65: |order - target| <= tol
    } cases[] = {
        {"exact wall w (Poisson limits)", MMS_EXACT_WALL_W, 2.0, 0.15},
        {"exact psi (wall closure limits)", MMS_EXACT_PSI, 2.0, 0.15},
        {"exact psi and wall w", MMS_EXACT_PSI | MMS_EXACT_WALL_W, 4.0, 0.3},
    };
    int k;
    char name[96];

    printf("Unit: what limits the spatial order (ablation), order 6, 33x33 -> 65x65\n");
    for (k = 0; k < (int)(sizeof(cases) / sizeof(cases[0])); k++)
    {
        mms_errors a = mms_run_ablated(33, 33, 1.0, 1.0, 100.0, 6, cases[k].flags, 2.5E-3, 0.0, 0.25);
        mms_errors b = mms_run_ablated(65, 65, 1.0, 1.0, 100.0, 6, cases[k].flags, 2.5E-3, 0.0, 0.25);
        double order = log2(a.w.max / b.w.max);
        snprintf(name, sizeof(name), "%s: w order %.2f, |order - %.0f|", cases[k].name, order, cases[k].target);
        check(name, isnan(order) ? INFINITY : fabs(order - cases[k].target), cases[k].tol);
    }
}

// Periodic operators: the wrapped centered stencils differentiate sin(2 pi x)
// at their nominal order, and annihilate constants
static void test_periodic_operators(void)
{
    int o, n, i;
    char name[96];

    printf("Unit: periodic derivative operators\n");
    for (o = 2; o <= 6; o += 2)
    {
        double err[2][2];
        for (n = 32; n <= 64; n *= 2)
        {
            double dx = 1.0 / n, e1 = 0.0, e2 = 0.0, ones = 0.0;
            smtrx d1 = SDiff1_periodic(n, o, dx), d2 = SDiff2_periodic(n, o, dx);
            double *f = (double *)malloc(n * sizeof(double)), *d = (double *)malloc(n * sizeof(double));
            for (i = 0; i < n; i++)
                f[i] = sin(2.0 * PI * i * dx);
            spmv(d1, f, d);
            for (i = 0; i < n; i++)
                e1 = fmax(e1, fabs(d[i] - 2.0 * PI * cos(2.0 * PI * i * dx)));
            spmv(d2, f, d);
            for (i = 0; i < n; i++)
                e2 = fmax(e2, fabs(d[i] + 4.0 * PI * PI * f[i]));
            for (i = 0; i < n; i++)
                f[i] = 1.0;
            spmv(d1, f, d);
            for (i = 0; i < n; i++)
                ones = fmax(ones, fabs(d[i]));
            spmv(d2, f, d);
            for (i = 0; i < n; i++)
                ones = fmax(ones, fabs(d[i]) * dx);
            err[n == 64][0] = e1;
            err[n == 64][1] = e2;
            if (n == 64)
            {
                snprintf(name, sizeof(name), "order %d: constants differentiate to zero", o);
                check(name, ones, 1E-10);
            }
            free(f);
            free(d);
            freesm(d1);
            freesm(d2);
        }
        double p1 = log2(err[0][0] / err[1][0]), p2 = log2(err[0][1] / err[1][1]);
        snprintf(name, sizeof(name), "order %d: first derivative order %.2f", o, p1);
        check(name, fabs(p1 - o), 0.1);
        snprintf(name, sizeof(name), "order %d: second derivative order %.2f", o, p2);
        check(name, fabs(p2 - o), 0.1);
    }
}

// The periodic Poisson solver solves (DX2 + DY2) psi = f exactly, for the
// operators of every order, on grids with odd and even sizes, and returns
// psi with zero mean
static void test_periodic_poisson(int nx, int ny, int order)
{
    int k, N = nx * ny;
    double dx = 1.0 / nx, dy = 2.0 / ny, mean = 0.0, res = 0.0, scale = 0.0;
    smtrx d2x = SDiff2_periodic(nx, order, dx), d2y = SDiff2_periodic(ny, order, dy);
    smtrx Ix = seye(nx), Iy = seye(ny);
    smtrx DX2 = skronecker(Iy, d2x), DY2 = skronecker(d2y, Ix);
    mtrx f = initm(ny, nx), psi = initm(ny, nx);
    double *a = (double *)malloc(N * sizeof(double)), *b = (double *)malloc(N * sizeof(double));
    periodic_solver *s = periodic_setup(nx, ny, &DX2, &DY2);

    fill_pseudo_random(f.M, N, 11);
    for (k = 0; k < N; k++)
        mean += f.M[k] / N;
    for (k = 0; k < N; k++)
        f.M[k] -= mean; // a periodic Poisson problem needs a zero-mean right-hand side
    poisson_periodic(s, f, psi);
    spmv(DX2, psi.M, a);
    spmv(DY2, psi.M, b);
    mean = 0.0;
    for (k = 0; k < N; k++)
    {
        res = fmax(res, fabs(a[k] + b[k] - f.M[k]));
        scale = fmax(scale, fabs(f.M[k]));
        mean += psi.M[k] / N;
    }
    printf("Unit: periodic Poisson solver, %dx%d grid, order %d\n", nx, ny, order);
    check("residual of (DX2 + DY2) psi = f, relative", res / scale, 1E-10);
    check("psi has zero mean", fabs(mean), 1E-12);

    periodic_cleanup(s);
    free(a);
    free(b);
    freem(&f);
    freem(&psi);
    freesm(d2x);
    freesm(d2y);
    freesm(Ix);
    freesm(Iy);
    freesm(DX2);
    freesm(DY2);
}

// Without walls the method reaches the nominal order of its stencils: the
// periodic manufactured solution converges at 2, 4 and 6
static void test_periodic_mms_order(void)
{
    int o;
    char name[96];

    printf("Unit: spatial order on a periodic grid (manufactured solution), RK4\n");
    for (o = 2; o <= 6; o += 2)
    {
        mms_errors a = mms_run_periodic(32, 32, 1.0, 1.0, 100.0, o, 2, 2.5E-3, 0.0, 0.25);
        mms_errors b = mms_run_periodic(64, 64, 1.0, 1.0, 100.0, o, 2, 2.5E-3, 0.0, 0.25);
        double pw = log2(a.w.max / b.w.max), pp = log2(a.psi.max / b.psi.max);
        snprintf(name, sizeof(name), "order %d: w order %.2f (32 -> 64)", o, pw);
        check(name, isnan(pw) ? INFINITY : fabs(pw - o), 0.15);
        snprintf(name, sizeof(name), "order %d: psi order %.2f", o, pp);
        check(name, isnan(pp) ? INFINITY : fabs(pp - o), 0.15);
    }
    printf("Unit: ... with the skew-symmetric nonlinear term\n");
    for (o = 2; o <= 6; o += 2)
    {
        // The products u w, v w carry twice the wavenumbers of the fields, so
        // the asymptotic range starts on finer grids than for the advective form
        mms_errors a = mms_run_periodic_advection(64, 64, 1.0, 1.0, 100.0, o, 1, 2.5E-3, 0.0, 0.25);
        mms_errors b = mms_run_periodic_advection(128, 128, 1.0, 1.0, 100.0, o, 1, 2.5E-3, 0.0, 0.25);
        double pw = log2(a.w.max / b.w.max);
        snprintf(name, sizeof(name), "skew, order %d: w order %.2f (64 -> 128)", o, pw);
        check(name, isnan(pw) ? INFINITY : fabs(pw - o), 0.15);
    }
    // A 2 x 1 domain with dx != dy
    mms_errors a = mms_run_periodic(48, 32, 2.0, 1.0, 100.0, 6, 2, 2.5E-3, 0.0, 0.25);
    mms_errors b = mms_run_periodic(96, 64, 2.0, 1.0, 100.0, 6, 2, 2.5E-3, 0.0, 0.25);
    double pw = log2(a.w.max / b.w.max);
    snprintf(name, sizeof(name), "2x1 domain, 48x32 -> 96x64, order 6: w order %.2f", pw);
    check(name, isnan(pw) ? INFINITY : fabs(pw - 6.0), 0.3);
}

// The compact Poisson operator (poisson_order 4) is at least fourth order, the
// 5-point one second order, on psi = sin(pi x) sin(pi y) e^(x + 2y). Its Laplacian is
// not zero on the walls, so the extrapolated wall values of the compact
// right-hand side are exercised.
static double compact_poisson_error(int nx, int ny, int order)
{
    int i, j;
    double dx = 1.0 / (nx - 1), dy = 1.0 / (ny - 1), err = 0.0;
    mtrx f = initm(ny, nx), u = initm(ny, nx);
    fft_solver *fs = fft_setup(nx, ny);

    for (i = 0; i < ny; i++)
        for (j = 0; j < nx; j++)
        {
            double x = j * dx, y = i * dy;
            double S = sin(PI * x) * exp(x), Sxx = ((1.0 - PI * PI) * sin(PI * x) + 2.0 * PI * cos(PI * x)) * exp(x);
            double T = sin(PI * y) * exp(2.0 * y), Tyy = ((4.0 - PI * PI) * sin(PI * y) + 4.0 * PI * cos(PI * y)) * exp(2.0 * y);
            MAt(f, i, j) = Sxx * T + S * Tyy;
        }
    poisson_FFT_order(fs, f, u, dx, dy, order);
    for (i = 0; i < ny; i++)
        for (j = 0; j < nx; j++)
            err = fmax(err, fabs(MAt(u, i, j) - sin(PI * j * dx) * exp(j * dx) * sin(PI * i * dy) * exp(2.0 * i * dy)));
    fft_cleanup(fs);
    freem(&f);
    freem(&u);
    return err;
}

static void test_compact_poisson(void)
{
    int order;
    char name[96];

    printf("Unit: FFT Poisson solver, 5-point and compact operators\n");
    for (order = 2; order <= 4; order += 2)
    {
        double sq = log2(compact_poisson_error(33, 33, order) / compact_poisson_error(65, 65, order));
        double rect = log2(compact_poisson_error(33, 17, order) / compact_poisson_error(65, 33, order));
        // The compact operator converges at about five here (fourth order is
        // its guaranteed rate), so it is checked from below
        snprintf(name, sizeof(name), "poisson_order %d: order %.2f (33 -> 65 squared)", order, sq);
        check(name, order == 2 ? fabs(sq - 2.0) : 3.8 - sq, order == 2 ? 0.2 : 0.0);
        snprintf(name, sizeof(name), "poisson_order %d: order %.2f (33x17 -> 65x33)", order, rect);
        check(name, order == 2 ? fabs(rect - 2.0) : 3.8 - rect, order == 2 ? 0.2 : 0.0);
    }
}

// SDiff1_wall4: its rows at and next to both ends differentiate quartics
// exactly; the other rows are SDiff1's
static void test_velocity_rows(void)
{
    int n = 21, o, i, p;
    double dx = 0.05, worst = 0.0, inner = 0.0;
    double f[21], d[21], d0[21];

    printf("Unit: fourth-order velocity rows next to the walls\n");
    for (o = 4; o <= 6; o += 2)
    {
        smtrx A = SDiff1_wall4(n, o, dx), B = SDiff1(n, o, dx);
        for (p = 0; p <= 4; p++)
        {
            for (i = 0; i < n; i++)
                f[i] = pow(i * dx - 0.3, p);
            spmv(A, f, d);
            spmv(B, f, d0);
            for (i = 0; i < n; i++)
            {
                double exact = p == 0 ? 0.0 : p * pow(i * dx - 0.3, p - 1);
                if (i <= 1 || i >= n - 2)
                    worst = fmax(worst, fabs(d[i] - exact));
                else
                    inner = fmax(inner, fabs(d[i] - d0[i]));
            }
        }
        freesm(A);
        freesm(B);
    }
    check("rows 0, 1, n-2, n-1 exact on polynomials up to degree 4", worst, 1E-10);
    check("other rows identical to SDiff1", inner, 0.0);
}

// Raising the wall-bounded order needs both the compact Poisson operator and
// the third-order wall closure (issues #26, #29): together the manufactured
// solution converges at about four; the wall closure alone stays at two
static void test_wall_order_pairing(void)
{
    char name[96];
    mms_errors a, b;
    double order;

    printf("Unit: compact Poisson operator and third-order wall closure, order 6, 33 -> 65\n");
    a = mms_run_closures(33, 33, 1.0, 1.0, 100.0, 6, 4, 1, 2, 2.5E-3, 0.0, 0.25);
    b = mms_run_closures(65, 65, 1.0, 1.0, 100.0, 6, 4, 1, 2, 2.5E-3, 0.0, 0.25);
    order = log2(a.w.max / b.w.max);
    snprintf(name, sizeof(name), "both: w order %.2f, 3.7 - order", order);
    check(name, isnan(order) ? INFINITY : 3.7 - order, 0.0);
    order = log2(a.psi.max / b.psi.max);
    snprintf(name, sizeof(name), "both: psi order %.2f, 3.7 - order", order);
    check(name, isnan(order) ? INFINITY : 3.7 - order, 0.0);
    check("both: w max error at 65x65", b.w.max, 2E-4);
    a = mms_run_closures(33, 33, 1.0, 1.0, 100.0, 6, 2, 1, 2, 2.5E-3, 0.0, 0.25);
    b = mms_run_closures(65, 65, 1.0, 1.0, 100.0, 6, 2, 1, 2, 2.5E-3, 0.0, 0.25);
    order = log2(a.w.max / b.w.max);
    snprintf(name, sizeof(name), "wall closure alone: w order %.2f, |order - 2|", order);
    check(name, isnan(order) ? INFINITY : fabs(order - 2.0), 0.2);

    // With fourth-order velocity rows as well, u and v converge at four too
    a = mms_run_closures(33, 33, 1.0, 1.0, 100.0, 6, 4, 1, 4, 2.5E-3, 0.0, 0.25);
    b = mms_run_closures(65, 65, 1.0, 1.0, 100.0, 6, 4, 1, 4, 2.5E-3, 0.0, 0.25);
    order = log2(a.u.max / b.u.max);
    snprintf(name, sizeof(name), "all three: u order %.2f, 3.7 - order", order);
    check(name, isnan(order) ? INFINITY : 3.7 - order, 0.0);
    order = log2(a.w.max / b.w.max);
    snprintf(name, sizeof(name), "all three: w order %.2f, 3.7 - order", order);
    check(name, isnan(order) ? INFINITY : 3.7 - order, 0.0);
    check("all three: u max error at 65x65", b.u.max, 3E-6);
}

// Integrals of the Taylor-Green vortex psi = sin(kx) sin(ky) e^(-2 nu k^2 t)/k,
// k = 2 pi: E = 1/4 e^(-4 nu k^2 t), Z = k^2/2 and P = k^4 times that factor
// over 1/4; and its whole spectrum in the shell |k| = sqrt(2) 2 pi
static void test_integrals_taylor_green(void)
{
    int i, j, t, n = 32, steps = 100, B;
    double k = 2.0 * PI, nu = 0.01, dt = 1E-3, *E, total = 0.0;
    problem p;
    spectra *sp;
    char name[96];

    printf("Diagnostics: Taylor-Green vortex, %dx%d periodic grid\n", n, n);
    problem_init_ext(&p, n, n, 1.0, 1.0, 1, 6, 1.0 / nu, 2, 3, dt, 1E-10, NULL, 0.0, NULL, NULL, 2, 0, 2);
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
        {
            double x = j * p.cfg.dx, y = i * p.cfg.dy;
            MAt(p.w, i, j) = 2.0 * k * sin(k * x) * sin(k * y);
            MAt(p.u, i, j) = sin(k * x) * cos(k * y);
            MAt(p.v, i, j) = -cos(k * x) * sin(k * y);
        }
    flow_integrals f = compute_integrals(&p.cfg, p.u, p.v, p.w, p.ctx.k1.M, p.ctx.k2.M);
    check("E at t = 0 is 1/4", fabs(f.E - 0.25), 1E-14);
    check("Z at t = 0 is k^2/2", fabs(f.Z - k * k / 2.0) / (k * k / 2.0), 1E-14);
    check("P at t = 0 is k^4 (finite differences, order 6)", fabs(f.P - k * k * k * k) / (k * k * k * k), 1E-5);
    for (t = 0; t < steps; t++)
        step(p.w, p.u, p.v, &p.ctx);
    f = compute_integrals(&p.cfg, p.u, p.v, p.w, p.ctx.k1.M, p.ctx.k2.M);
    snprintf(name, sizeof(name), "E at t = %g against 1/4 e^(-4 nu k^2 t)", steps * dt);
    check(name, fabs(f.E - 0.25 * exp(-4.0 * nu * k * k * steps * dt)) / f.E, 1E-6);

    sp = spectra_setup(&p.cfg);
    B = spectra_bins(sp);
    E = (double *)calloc(B, sizeof(double));
    spectra_compute(sp, p.u, p.v, p.w, E, NULL, NULL, NULL);
    for (i = 0; i < B; i++)
        if (i != 1) total += E[i];
    check("all energy in the shell of |k| = sqrt(2) 2 pi (bin 1)", fabs(E[1] - f.E) / f.E + total / f.E, 1E-12);
    free(E);
    spectra_free(sp);
    problem_free(&p);
}

// On the unforced three-mode periodic flow of mms.h: the spectra sum to the
// integrals (Parseval), the nonlinear term conserves energy and enstrophy up
// to the order of the scheme, and dE/dt = -2 nu Z up to that order too
static void budget_run(int n, double out[4])
{
    int t, B, k;
    double nu = 0.01, dt = 1E-4, *E, *Z, *PE, *PZ, sE = 0.0, sZ = 0.0, mE = 0.0, mZ = 0.0;
    mms_case c = {1.0, 1.0, 1.0 / nu, 1.0 / n, 1.0 / n, 1};
    problem p;
    spectra *sp;

    problem_init_ext(&p, n, n, 1.0, 1.0, 1, 6, 1.0 / nu, 2, 3, dt, 1E-10, NULL, 0.0, NULL, NULL, 2, 0, 2);
    mms_exact(&c, 0.0, &p.w, &p.u, &p.v, NULL);
    for (t = 0; t < 10; t++)
        step(p.w, p.u, p.v, &p.ctx);
    sp = spectra_setup(&p.cfg);
    B = spectra_bins(sp);
    E = (double *)calloc(B, sizeof(double));
    Z = (double *)calloc(B, sizeof(double));
    PE = (double *)calloc(B, sizeof(double));
    PZ = (double *)calloc(B, sizeof(double));
    flow_integrals f0 = compute_integrals(&p.cfg, p.u, p.v, p.w, p.ctx.k1.M, p.ctx.k2.M);
    spectra_compute(sp, p.u, p.v, p.w, E, Z, PE, PZ);
    for (k = 0; k < B; k++)
    {
        sE += E[k];
        sZ += Z[k];
        mE = fmax(mE, fabs(PE[k]));
        mZ = fmax(mZ, fabs(PZ[k]));
    }
    step(p.w, p.u, p.v, &p.ctx);
    flow_integrals f1 = compute_integrals(&p.cfg, p.u, p.v, p.w, p.ctx.k1.M, p.ctx.k2.M);
    out[0] = fabs(sE - f0.E) / f0.E + fabs(sZ - f0.Z) / f0.Z;                      // Parseval
    out[1] = fabs(PE[B - 1]) / mE;                                                 // net nonlinear energy transfer
    out[2] = fabs(PZ[B - 1]) / mZ;                                                 // ... and enstrophy transfer
    out[3] = fabs((f1.E - f0.E) / dt + nu * (f0.Z + f1.Z)) / (nu * (f0.Z + f1.Z)); // budget
    free(E);
    free(Z);
    free(PE);
    free(PZ);
    spectra_free(sp);
    problem_free(&p);
}

// The spectral fluxes are those of the discrete equations, exactly: the net
// nonlinear transfer summed over all shells equals the rate at which the
// solver's nonlinear term changes E = 1/2 <u^2 + v^2> and Z = 1/2 <w^2>,
// computed in physical space with the solver's own operators, to round-off.
// On coarse grids, where the first-derivative symbols differ from the
// Laplacian's, this needs the A/Q weight in the energy transfer.
static void flux_consistency(int nx, int ny, int order, double out[2])
{
    int k, B, N = nx * ny;
    double rE = 0.0, rZ = 0.0, mE = 0.0, mZ = 0.0;
    problem p;
    problem_init_ext(&p, nx, ny, 1.0, 1.0, 1, order, 100., 2, 3, 1E-3, 1E-10, NULL, 0.0, NULL, NULL, 2, 0, 2);
    periodic_solver *ps = periodic_setup(nx, ny, &p.DX2, &p.DY2);
    mtrx f = initm(ny, nx), psi = initm(ny, nx), nl = initm(ny, nx), un = initm(ny, nx), vn = initm(ny, nx);
    double *a = (double *)malloc(N * sizeof(double)), *b = (double *)malloc(N * sizeof(double));
    spectra *sp = spectra_setup(&p.cfg);
    double *PE, *PZ;

    // A pseudo-random w, with energy at every wavenumber (a few smooth modes
    // would not exchange energy at all: their products fall outside the
    // field's modes); psi, u and v as the solver makes them
    fill_pseudo_random(p.w.M, N, 7u);
    negcpy(f, p.w);
    poisson_periodic(ps, f, psi);
    spmv(p.DY, psi.M, p.u.M);
    spmv(p.DX, psi.M, p.v.M);
    negcpy(p.v, p.v);
    // The nonlinear term, and the velocity change it causes
    spmv(p.DX, p.w.M, a);
    spmv(p.DY, p.w.M, b);
    for (k = 0; k < N; k++)
        nl.M[k] = -(p.u.M[k] * a[k] + p.v.M[k] * b[k]);
    negcpy(f, nl);
    poisson_periodic(ps, f, psi);
    spmv(p.DY, psi.M, un.M);
    spmv(p.DX, psi.M, vn.M);
    negcpy(vn, vn);
    for (k = 0; k < N; k++)
    {
        rE += (p.u.M[k] * un.M[k] + p.v.M[k] * vn.M[k]) / N;
        rZ += p.w.M[k] * nl.M[k] / N;
    }
    B = spectra_bins(sp);
    PE = (double *)calloc(B, sizeof(double));
    PZ = (double *)calloc(B, sizeof(double));
    spectra_compute(sp, p.u, p.v, p.w, NULL, NULL, PE, PZ);
    for (k = 0; k < B; k++)
    {
        mE = fmax(mE, fabs(PE[k]));
        mZ = fmax(mZ, fabs(PZ[k]));
    }
    // The flux through the last shell is minus the net transfer
    out[0] = fabs(rE + PE[B - 1]) / mE;
    out[1] = fabs(rZ + PZ[B - 1]) / mZ;

    free(a);
    free(b);
    free(PE);
    free(PZ);
    freem(&f);
    freem(&psi);
    freem(&nl);
    freem(&un);
    freem(&vn);
    spectra_free(sp);
    periodic_cleanup(ps);
    problem_free(&p);
}

// The energy budget of the discrete equations closes: on a coarse grid with a
// random field, where A/Q is far from 1, the centred difference of E over two
// RK4 steps equals the nonlinear transfer minus the discrete viscous
// dissipation sum D_E to the time-stepping error, while the continuum budget
// dE/dt = -2 nu Z is off by far more
static void discrete_budget(int n, int order, double out[2])
{
    int k, B, N = n * n;
    double dt = 5E-5, nu = 0.01, rate, Em, Ep, sumD = 0.0, *PE, *DE;
    problem p;
    spectra *sp;

    problem_init_ext(&p, n, n, 1.0, 1.0, 1, order, 1.0 / nu, 2, 3, dt, 1E-10, NULL, 0.0, NULL, NULL, 2, 0, 2);
    fill_pseudo_random(p.w.M, N, 3u);
    step(p.w, p.u, p.v, &p.ctx); // makes u, v consistent with w
    Em = compute_integrals(&p.cfg, p.u, p.v, p.w, p.ctx.k1.M, p.ctx.k2.M).E;
    step(p.w, p.u, p.v, &p.ctx);
    flow_integrals f0 = compute_integrals(&p.cfg, p.u, p.v, p.w, p.ctx.k1.M, p.ctx.k2.M);
    sp = spectra_setup(&p.cfg);
    B = spectra_bins(sp);
    PE = (double *)calloc(B, sizeof(double));
    DE = (double *)calloc(B, sizeof(double));
    spectra_compute(sp, p.u, p.v, p.w, NULL, NULL, PE, NULL);
    spectra_dissipation(sp, p.w, DE, NULL, NULL, NULL);
    for (k = 0; k < B; k++)
        sumD += DE[k];
    rate = -PE[B - 1] - sumD;
    step(p.w, p.u, p.v, &p.ctx);
    Ep = compute_integrals(&p.cfg, p.u, p.v, p.w, p.ctx.k1.M, p.ctx.k2.M).E;
    out[0] = fabs((Ep - Em) / (2.0 * dt) - rate) / fabs(rate);
    out[1] = fabs((Ep - Em) / (2.0 * dt) + 2.0 * nu * f0.Z) / fabs(rate);
    free(PE);
    free(DE);
    spectra_free(sp);
    problem_free(&p);
}

// The same with every damping term: drag, hypodrag and hyperviscosity of order
// p. dE/dt = transfer - sum (D_E + F_E) and dZ/dt = transfer - sum (D_Z + F_Z),
// each to the time-stepping error. out: the energy and the enstrophy residual.
static void damped_budget(int n, int hyper_order, double out[2])
{
    int k, B, N = n * n;
    double dt = 2E-5, nu = 0.01, rate[2] = {0.0, 0.0}, Em, Ep, Zm, Zp, *PE, *PZ, *D[4];
    problem p;
    spectra *sp;

    problem_init_ext(&p, n, n, 1.0, 1.0, 1, 4, 1.0 / nu, 2, 3, dt, 1E-10, NULL, 0.0, NULL, NULL, 2, 0, 2);
    rk4_free(&p.ctx);
    p.cfg.forcing.drag = 0.3;
    p.cfg.forcing.hypodrag = 20.0;
    p.cfg.forcing.hyperviscosity = hyper_order == 2 ? 1E-6 : 1E-10;
    p.cfg.forcing.hyper_order = hyper_order;
    p.ctx = rk4_alloc(&p.cfg);
    fill_pseudo_random(p.w.M, N, 5u);
    step(p.w, p.u, p.v, &p.ctx);
    flow_integrals fm = compute_integrals(&p.cfg, p.u, p.v, p.w, p.ctx.k1.M, p.ctx.k2.M);
    Em = fm.E;
    Zm = fm.Z;
    step(p.w, p.u, p.v, &p.ctx);
    sp = spectra_setup(&p.cfg);
    B = spectra_bins(sp);
    PE = (double *)calloc(B, sizeof(double));
    PZ = (double *)calloc(B, sizeof(double));
    for (k = 0; k < 4; k++)
        D[k] = (double *)calloc(B, sizeof(double));
    spectra_compute(sp, p.u, p.v, p.w, NULL, NULL, PE, PZ);
    spectra_dissipation(sp, p.w, D[0], D[1], D[2], D[3]);
    rate[0] = -PE[B - 1];
    rate[1] = -PZ[B - 1];
    for (k = 0; k < B; k++)
    {
        rate[0] -= D[0][k] + D[2][k];
        rate[1] -= D[1][k] + D[3][k];
    }
    step(p.w, p.u, p.v, &p.ctx);
    flow_integrals fp = compute_integrals(&p.cfg, p.u, p.v, p.w, p.ctx.k1.M, p.ctx.k2.M);
    Ep = fp.E;
    Zp = fp.Z;
    out[0] = fabs((Ep - Em) / (2.0 * dt) - rate[0]) / fabs(rate[0]);
    out[1] = fabs((Zp - Zm) / (2.0 * dt) - rate[1]) / fabs(rate[1]);
    free(PE);
    free(PZ);
    for (k = 0; k < 4; k++)
        free(D[k]);
    spectra_free(sp);
    problem_free(&p);
}

// The net nonlinear enstrophy transfer, -Pi_Z through the last shell, on a
// random field, relative to the largest transfer in a shell: of order one for
// the advective form, round-off for the skew-symmetric one
static double enstrophy_defect(int nx, int ny, int order, int advection)
{
    int k, B, N = nx * ny;
    double *PZ, peak = 0.0, net;
    problem p;
    spectra *sp;

    problem_init_ext(&p, nx, ny, 1.0, 1.0, 1, order, 100.0, 2, 3, 1E-4, 1E-10, NULL, 0.0, NULL, NULL, 2, 0, 2);
    rk4_free(&p.ctx);
    p.cfg.advection = advection;
    p.ctx = rk4_alloc(&p.cfg);
    fill_pseudo_random(p.w.M, N, 9u);
    step(p.w, p.u, p.v, &p.ctx);
    sp = spectra_setup(&p.cfg);
    B = spectra_bins(sp);
    PZ = (double *)calloc(B, sizeof(double));
    spectra_compute(sp, p.u, p.v, p.w, NULL, NULL, NULL, PZ);
    for (k = 0; k < B; k++)
        peak = fmax(peak, fabs(PZ[k] - (k ? PZ[k - 1] : 0.0)));
    net = fabs(PZ[B - 1]) / peak;
    free(PZ);
    spectra_free(sp);
    problem_free(&p);
    return net;
}

static void test_budgets(void)
{
    double a[4], b[4], o;
    char name[96];

    printf("Diagnostics: net spectral transfer = the nonlinear rate of change of E and Z\n");
    for (int o = 2; o <= 6; o += 4)
    {
        double r[2];
        flux_consistency(16, 16, o, r);
        snprintf(name, sizeof(name), "order %d, 16x16: energy %.1e, enstrophy %.1e", o, r[0], r[1]);
        check(name, r[0] + r[1], 1E-12);
        flux_consistency(24, 20, o, r);
        snprintf(name, sizeof(name), "order %d, 24x20: energy %.1e, enstrophy %.1e", o, r[0], r[1]);
        check(name, r[0] + r[1], 1E-12);
    }

    printf("Diagnostics: the discrete energy budget closes, random field\n");
    for (int o = 2; o <= 6; o += 4)
    {
        double r[2];
        discrete_budget(24, o, r);
        snprintf(name, sizeof(name), "order %d, 24x24: discrete budget %.1e (continuum %.1e)", o, r[0], r[1]);
        check(name, r[0], 1E-5);
    }

    printf("Diagnostics: the skew-symmetric nonlinear term conserves enstrophy, random field\n");
    for (int o = 2; o <= 6; o += 2)
    {
        double adv = enstrophy_defect(24, 20, o, 0), skew = enstrophy_defect(24, 20, o, 1);
        snprintf(name, sizeof(name), "order %d, 24x20: net transfer %.1e (advective form %.1e)", o, skew, adv);
        check(name, skew, 1E-13);
    }

    printf("Diagnostics: ... with drag, hypodrag and hyperviscosity, random field, order 4\n");
    for (int hp = 2; hp <= 3; hp++)
    {
        double r[2];
        damped_budget(24, hp, r);
        snprintf(name, sizeof(name), "hyperviscosity order %d: energy %.1e, enstrophy %.1e", hp, r[0], r[1]);
        check(name, fmax(r[0], r[1]), 2E-6);
    }

    printf("Diagnostics: spectra and budgets, unforced periodic flow, order 6, 32x32 -> 64x64\n");
    budget_run(32, a);
    budget_run(64, b);
    check("Parseval: spectra sum to E and Z", fmax(a[0], b[0]), 1E-13);
    o = log2(a[1] / b[1]);
    snprintf(name, sizeof(name), "net nonlinear energy transfer %.1e, order %.2f", b[1], o);
    check(name, isnan(o) ? INFINITY : 5.0 - o, 0.0);
    o = log2(a[2] / b[2]);
    snprintf(name, sizeof(name), "net nonlinear enstrophy transfer %.1e, order %.2f", b[2], o);
    check(name, isnan(o) ? INFINITY : 5.0 - o, 0.0);
    o = log2(a[3] / b[3]);
    snprintf(name, sizeof(name), "dE/dt + 2 nu Z, relative %.1e, order %.2f", b[3], o);
    check(name, isnan(o) ? INFINITY : 5.0 - o, 0.0);
}

// Kolmogorov forcing with drag reaches the laminar state, a single mode in y on
// which the nonlinear term vanishes: w = f_K / (nu Q + drag), Q = -(symbol of
// DY2). The energy input balances the dissipation, I = 2 nu Z + 2 drag E.
static void test_kolmogorov_laminar(void)
{
    int i, j, t, n = 32, kn = 2;
    double A = 1.0, nu = 0.05, drag = 0.1, k = 2.0 * PI * kn, err_d = 0.0, err_c = 0.0, peak = 0.0, Q;
    problem p;
    char name[96];

    printf("Forcing: laminar Kolmogorov flow, %dx%d periodic grid, order 6\n", n, n);
    problem_init_ext(&p, n, n, 1.0, 1.0, 1, 6, 1.0 / nu, 2, 3, 2E-3, 1E-10, NULL, 0.0, NULL, NULL, 2, 0, 2);
    // The forcing is part of the copied configuration, so set it before the
    // workspace is made: reallocate the workspace with it
    rk4_free(&p.ctx);
    p.cfg.forcing.kolmogorov_amp = A;
    p.cfg.forcing.kolmogorov_n = kn;
    p.cfg.forcing.drag = drag;
    p.ctx = rk4_alloc(&p.cfg);
    for (t = 0; t < 1500; t++)
        step(p.w, p.u, p.v, &p.ctx);

    // Q from the operator itself: DY2 cos(k y) = -Q cos(k y)
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            p.ctx.k3.M[i * n + j] = cos(k * i * p.cfg.dy);
    spmv(p.DY2, p.ctx.k3.M, p.ctx.k4.M);
    Q = -p.ctx.k4.M[0];
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
        {
            double f = -A * k * cos(k * i * p.cfg.dy);
            err_d = fmax(err_d, fabs(MAt(p.w, i, j) - f / (nu * Q + drag)));
            err_c = fmax(err_c, fabs(MAt(p.w, i, j) - f / (nu * k * k + drag)));
            peak = fmax(peak, fabs(f / (nu * k * k + drag)));
        }
    check("steady w = f / (nu Q + drag), Q of the discrete operator", err_d / peak, 1E-10);
    check("... and within O(h^6) of the continuum f / (nu k^2 + drag)", err_c / peak, 1E-5);
    flow_integrals fi = compute_integrals(&p.cfg, p.u, p.v, p.w, p.ctx.k1.M, p.ctx.k2.M);
    // I is the continuum work <u A sin(k y)> and 2 nu Z the continuum
    // dissipation, so this balance holds to the order of the scheme
    snprintf(name, sizeof(name), "continuum balance I = 2 nu Z + 2 drag E (I = %.4f)", fi.I);
    check(name, fabs(fi.I - 2.0 * nu * fi.Z - 2.0 * drag * fi.E) / fi.I, 1E-4);
    // With the discrete injection and dissipation it holds to the steady
    // state's convergence
    spectra *sp = spectra_setup(&p.cfg);
    int b, B = spectra_bins(sp);
    double *DE = (double *)calloc(B, sizeof(double)), sumD = 0.0;
    spectra_dissipation(sp, p.w, DE, NULL, NULL, NULL);
    for (b = 0; b < B; b++)
        sumD += DE[b];
    snprintf(name, sizeof(name), "discrete balance I_disc = sum D_E + 2 drag E (I_disc / I = %.6f)", fi.I_disc / fi.I);
    check(name, fabs(fi.I_disc - sumD - 2.0 * drag * fi.E) / fi.I_disc, 1E-9);
    free(DE);
    spectra_free(sp);
    problem_free(&p);
}

// Growth of max |w| over 300 steps of a single Fourier mode (on which the
// nonlinear term vanishes) at `factor` times the stability limit of the
// damping terms in f: (mx, my) = (1, 0) is the mode of the smallest
// eigenvalue, (n/2, n/2) that of the largest
static double damped_mode_growth(forcing_config f, int checkerboard, int time_scheme, double factor)
{
    int i, j, t, n = 16;
    double h = 1.0 / n, Re = 1E8, w0 = 0.0, w1 = 0.0;
    smtrx d2 = SDiff2_periodic(n, 6, h);
    double dt = factor * max_stable_dt_forced(&d2, &d2, Re, &f, time_scheme);
    problem p;

    problem_init_ext(&p, n, n, 1.0, 1.0, 1, 6, Re, time_scheme, 3, dt, 1E-10, NULL, 0.0, NULL, NULL, 2, 0, 2);
    rk4_free(&p.ctx);
    p.cfg.forcing = f;
    p.ctx = rk4_alloc(&p.cfg);
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
        {
            MAt(p.w, i, j) = checkerboard ? ((i + j) % 2 ? -1.0 : 1.0) : sin(2.0 * PI * j * h);
            MAt(p.u, i, j) = MAt(p.v, i, j) = 0.0;
            w0 = fmax(w0, fabs(MAt(p.w, i, j)));
        }
    for (t = 0; t < 300; t++)
        step(p.w, p.u, p.v, &p.ctx);
    // Past the limit every mode grows, and the round-off in the others soon
    // makes the nonlinear term blow the field up to inf or NaN
    for (i = 0; i < n * n; i++)
        w1 = isfinite(p.w.M[i]) ? fmax(w1, fabs(p.w.M[i])) : INFINITY;
    freesm(d2);
    problem_free(&p);
    return w1 / w0;
}

// The stability limit includes drag, hypodrag and hyperviscosity, each sharp:
// the mode that sets it decays at 0.95 times the limit and grows at 1.05
static void test_damping_stability(void)
{
    static const char *names[] = {"drag", "hypodrag", "hyperviscosity"};
    char name[96];

    printf("Stability: the time-step limit of the drag, hypodrag and hyperviscosity\n");
    for (int term = 0; term < 3; term++)
        for (int scheme = 1; scheme <= 2; scheme++)
        {
            forcing_config f = {0};
            if (term == 0) f.drag = 1000.0;
            if (term == 1) f.hypodrag = 1E4;
            if (term == 2)
            {
                f.hyperviscosity = 1E-6;
                f.hyper_order = 3;
            }
            double below = damped_mode_growth(f, term == 2, scheme, 0.95);
            double above = damped_mode_growth(f, term == 2, scheme, 1.05);
            snprintf(name, sizeof(name), "%s, %s: growth %.1e at 0.95x the limit, %.1e at 1.05x", names[term],
                     scheme == 1 ? "Euler" : "RK4", below, above);
            check(name, below <= 1.0 && above > 1E3 ? 0.0 : 1.0, 0.0);
        }
}

// The decaying-turbulence initial field: energy 1/2 <u^2 + v^2> as asked,
// u, v the velocity of w (u_x + v_y = 0 and v_x - u_y = w, spectrally exact,
// so to the order of the stencils with DX, DY), and the same field on a grid
// twice as fine
static void test_random_initial_field(void)
{
    int i, j, n = 64, N = n * n;
    double err = 0.0, peak = 0.0, e = 0.0, div = 0.0, curl = 0.0, umax = 0.0;
    mtrx w = initm(n, n), u = initm(n, n), v = initm(n, n);
    mtrx w2 = initm(2 * n, 2 * n), u2 = initm(2 * n, 2 * n), v2 = initm(2 * n, 2 * n);
    problem p;
    char name[96];

    printf("Initial field: decaying turbulence, %dx%d and %dx%d\n", n, n, 2 * n, 2 * n);
    random_initial_field(w, u, v, 1.0 / n, 1.0 / n, 4.0, 0.5, 3);
    random_initial_field(w2, u2, v2, 0.5 / n, 0.5 / n, 4.0, 0.5, 3);
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
        {
            err = fmax(err, fabs(MAt(w, i, j) - MAt(w2, 2 * i, 2 * j)));
            peak = fmax(peak, fabs(MAt(w, i, j)));
        }
    for (i = 0; i < N; i++)
        e += 0.5 * (u.M[i] * u.M[i] + v.M[i] * v.M[i]) / N;
    check("energy 1/2 <u^2 + v^2> = 0.5", fabs(e - 0.5), 1E-13);
    snprintf(name, sizeof(name), "same field on the grid twice as fine (max |w| %.0f)", peak);
    check(name, err / peak, 1E-12);

    problem_init_ext(&p, n, n, 1.0, 1.0, 1, 6, 100.0, 2, 3, 1E-3, 1E-10, NULL, 0.0, NULL, NULL, 2, 0, 2);
    spmv(p.DX, u.M, p.ctx.k1.M);
    spmv(p.DY, v.M, p.ctx.k2.M);
    for (i = 0; i < N; i++)
    {
        div = fmax(div, fabs(p.ctx.k1.M[i] + p.ctx.k2.M[i]));
        umax = fmax(umax, fabs(u.M[i]));
    }
    spmv(p.DX, v.M, p.ctx.k1.M);
    spmv(p.DY, u.M, p.ctx.k2.M);
    for (i = 0; i < N; i++)
        curl = fmax(curl, fabs(p.ctx.k1.M[i] - p.ctx.k2.M[i] - w.M[i]));
    snprintf(name, sizeof(name), "DX u + DY v = 0, DX v - DY u = w (order 6): %.1e, %.1e", div / umax, curl / peak);
    // Truncation error of order 6 on the spectrum (2 % at |k|/2 pi = 12):
    // this checks signs and axes, which would be wrong by O(1)
    check(name, fmax(div / umax, curl / peak), 1E-2);
    problem_free(&p);
    freem(&w);
    freem(&u);
    freem(&v);
    freem(&w2);
    freem(&u2);
    freem(&v2);
}

// Drag on the Taylor-Green vortex: E = 1/4 exp(-(4 nu k^2 + 2 drag) t)
static void test_drag_decay(void)
{
    int i, j, t, n = 32, steps = 200;
    double k = 2.0 * PI, nu = 0.01, drag = 0.5, dt = 1E-3;
    problem p;

    printf("Forcing: linear drag on the Taylor-Green vortex\n");
    problem_init_ext(&p, n, n, 1.0, 1.0, 1, 6, 1.0 / nu, 2, 3, dt, 1E-10, NULL, 0.0, NULL, NULL, 2, 0, 2);
    rk4_free(&p.ctx);
    p.cfg.forcing.drag = drag;
    p.ctx = rk4_alloc(&p.cfg);
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
        {
            double x = j * p.cfg.dx, y = i * p.cfg.dy;
            MAt(p.w, i, j) = 2.0 * k * sin(k * x) * sin(k * y);
            MAt(p.u, i, j) = sin(k * x) * cos(k * y);
            MAt(p.v, i, j) = -cos(k * x) * sin(k * y);
        }
    for (t = 0; t < steps; t++)
        step(p.w, p.u, p.v, &p.ctx);
    flow_integrals fi = compute_integrals(&p.cfg, p.u, p.v, p.w, p.ctx.k1.M, p.ctx.k2.M);
    check("E against 1/4 exp(-(4 nu k^2 + 2 drag) t)",
          fabs(fi.E - 0.25 * exp(-(4.0 * nu * k * k + 2.0 * drag) * steps * dt)) / fi.E, 1E-6);
    problem_free(&p);
}

// A random kick carries exactly eps dt of the solver's discrete energy, its
// modes lie in the shell, and two generators with the same seed agree
static void test_random_kick(void)
{
    int k, n = 48, N = n * n;
    double eps = 0.7, dt = 1E-3, kmin = INFINITY, kmax = 0.0;
    forcing_config fc = {0};
    problem p;
    random_forcing *a, *b;
    mtrx f = initm(n, n), psi = initm(n, n), u = initm(n, n), v = initm(n, n), w = initm(n, n), w2 = initm(n, n);
    double E = 0.0;

    printf("Forcing: random kicks, %dx%d periodic grid\n", n, n);
    problem_init_ext(&p, n, n, 1.0, 1.0, 1, 6, 100., 2, 3, dt, 1E-10, NULL, 0.0, NULL, NULL, 2, 0, 2);
    fc.random_rate = eps;
    fc.random_kf = 6.0;
    fc.random_dk = 1.0;
    fc.random_seed = 42;
    a = random_forcing_setup(&fc, n, n, p.cfg.dx, p.cfg.dy, dt, &p.DX, &p.DY, &p.DX2, &p.DY2);
    b = random_forcing_setup(&fc, n, n, p.cfg.dx, p.cfg.dy, dt, &p.DX, &p.DY, &p.DX2, &p.DY2);
    random_forcing_draw(a);
    random_forcing_add(a, w, p.cfg.dy);
    random_forcing_draw(b);
    random_forcing_add(b, w2, p.cfg.dy);
    periodic_solver *ps = periodic_setup(n, n, &p.DX2, &p.DY2);
    negcpy(f, w);
    poisson_periodic(ps, f, psi);
    spmv(p.DY, psi.M, u.M);
    spmv(p.DX, psi.M, v.M);
    for (k = 0; k < N; k++)
        E += 0.5 * (u.M[k] * u.M[k] + v.M[k] * v.M[k]) / N;
    for (k = 0; k < a->modes; k++)
    {
        double kk = sqrt(a->kx[k] * a->kx[k] + a->ky[k] * a->ky[k]) / (2.0 * PI);
        kmin = fmin(kmin, kk);
        kmax = fmax(kmax, kk);
    }
    check("one kick carries eps dt of discrete energy", fabs(E - eps * dt) / (eps * dt), 1E-12);
    check("its modes lie in the shell 5 <= |k| / 2 pi <= 7", (kmin < 5.0 - 1E-12) + (kmax > 7.0 + 1E-12), 0.0);
    check("same seed, same kick", rel_diff(w.M, w2.M, N), 0.0);
    random_forcing_draw(a);
    check("the next kick has new phases", a->phase[0] == b->phase[0], 0.0);

    random_forcing_free(a);
    random_forcing_free(b);
    periodic_cleanup(ps);
    freem(&f);
    freem(&psi);
    freem(&u);
    freem(&v);
    freem(&w);
    freem(&w2);
    problem_free(&p);
}

// Wall velocities go to the right nodes; at the corners the walls x = 0 and
// x = Lx win, on the CPU as on the GPU
static void test_wall_bc(void)
{
    int i, j, nx = 5, ny = 4;
    wall_bc bc = {{1., 2., 3., 4.}, {5., 6., 7., 8.}};
    mtrx u = initm(ny, nx), v = initm(ny, nx);
    double bad = 0.0;

    printf("Unit: wall velocities\n");
    for (i = 0; i < nx * ny; i++)
        u.M[i] = v.M[i] = -1.0;
    apply_wall_bc(u, v, &bc);
    for (i = 0; i < ny; i++)
        for (j = 0; j < nx; j++)
        {
            int wall = (j == 0) ? 0 : (j == nx - 1) ? 1
                                  : (i == 0)        ? 2
                                  : (i == ny - 1)   ? 3
                                                    : -1;
            double eu = wall < 0 ? -1.0 : bc.u[wall], ev = wall < 0 ? -1.0 : bc.v[wall];
            bad += fabs(MAt(u, i, j) - eu) + fabs(MAt(v, i, j) - ev);
        }
    check("every node has its wall's velocity, interior untouched", bad, 0.0);
    freem(&u);
    freem(&v);
}

// Count the numbers in a file after the line starting with `after`
static int count_values_after(const char *path, const char *after, int *found)
{
    char line[4096];
    int count = 0, on = 0;
    FILE *f = fopen(path, "r");
    *found = 0;
    if (!f) return -1;
    while (fgets(line, sizeof(line), f))
    {
        if (on)
        {
            char *p = line, *end;
            for (;;)
            {
                strtod(p, &end);
                if (end == p) break;
                count++;
                p = end;
            }
        }
        else if (strncmp(line, after, strlen(after)) == 0)
            on = *found = 1;
    }
    fclose(f);
    return count;
}

// The VTK and centerline writers, in a scratch directory
static void test_output(void)
{
    int i, j, nx = 4, ny = 3, found;
    char cwd[4096], tmpl[] = "/tmp/cnavier-test-XXXXXX", line[256];
    mtrx a = initm(ny, nx), u = initm(ny, nx), v = initm(ny, nx);
    FILE *f;
    double worst = 0.0;

    printf("Unit: output files\n");
    if (!getcwd(cwd, sizeof(cwd)) || !mkdtemp(tmpl) || chdir(tmpl) != 0 || mkdir("output", 0700) != 0)
    {
        check("scratch directory for the output tests", 1.0, 0.0);
        return;
    }

    // Files of an earlier run: one this run overwrites, one it would not reach
    // (a later frame of a longer run), and one of another series
    f = fopen("output/unittest-1-0.vtk", "w");
    if (f)
    {
        fprintf(f, "stale\n1 2 3\n");
        fclose(f);
    }
    f = fopen("output/unittest-1-5.vtk", "w");
    if (f)
    {
        fprintf(f, "stale\n");
        fclose(f);
    }
    f = fopen("output/other-1-5.vtk", "w");
    if (f)
    {
        fprintf(f, "keep\n");
        fclose(f);
    }
    for (i = 0; i < nx * ny; i++)
        a.M[i] = 0.5 * i;
    printvtk(a, "unittest", 1.0 / (nx - 1), 1.0 / (ny - 1));
    f = fopen("output/unittest-1-0.vtk", "r");
    int header = f && fgets(line, sizeof(line), f) && strncmp(line, "# vtk DataFile", 14) == 0;
    int dims = 0, spacing = 0;
    if (f)
    {
        while (fgets(line, sizeof(line), f))
        {
            if (strcmp(line, "DIMENSIONS 4 3 1\n") == 0) dims = 1;
            if (strncmp(line, "SPACING 0.33333333333333331 0.5 1", 33) == 0) spacing = 1;
        }
        fclose(f);
    }
    check("VTK file starts with its header (stale file replaced)", !header, 0.0);
    check("VTK DIMENSIONS lists x then y", !dims, 0.0);
    check("VTK SPACING is dx dy", !spacing, 0.0);
    check("earlier run's later frame removed, other series kept",
          (access("output/unittest-1-5.vtk", F_OK) == 0) + (access("output/other-1-5.vtk", F_OK) != 0), 0.0);
    check("VTK file has one value per node",
          (double)abs(count_values_after("output/unittest-1-0.vtk", "LOOKUP_TABLE", &found) - nx * ny) + !found, 0.0);

    // u = column index, v = row index: on the vertical centerline (halfway
    // between columns 1 and 2) u is 1.5, on the horizontal one (row 1) v is 1
    for (i = 0; i < ny; i++)
        for (j = 0; j < nx; j++)
        {
            MAt(u, i, j) = j;
            MAt(v, i, j) = i;
        }
    print_centerline(u, v, nx, ny, 1.0 / (nx - 1), 1.0 / (ny - 1));
    for (int c = 0; c < 2; c++)
    {
        double coord, val;
        int rows = 0;
        f = fopen(c == 0 ? "output/centerline_u_sim.csv" : "output/centerline_v_sim.csv", "r");
        if (!f || !fgets(line, sizeof(line), f))
        {
            worst = INFINITY;
            break;
        }
        while (fgets(line, sizeof(line), f))
        {
            char *end, *comma;
            coord = strtod(line, &comma);
            if (comma == line || *comma != ',') break;
            val = strtod(comma + 1, &end);
            if (end == comma + 1) break;
            double expect_coord = rows * (c == 0 ? 1.0 / (ny - 1) : 1.0 / (nx - 1));
            double d = fabs(val - (c == 0 ? 1.5 : 1.0)) + fabs(coord - expect_coord);
            // fmax would drop a NaN; count it as a failure
            worst = isnan(d) ? INFINITY : fmax(worst, d);
            rows++;
        }
        fclose(f);
        if (rows != (c == 0 ? ny : nx)) worst = INFINITY;
    }
    check("centerline CSVs: node coordinates, values on the centerline", worst, 1E-6);

    // Clean up the scratch directory
    const char *names[] = {"output/unittest-1-0.vtk", "output/other-1-5.vtk", "output/centerline_u_sim.csv",
                           "output/centerline_v_sim.csv", "output/centerline_u_ghia.csv",
                           "output/centerline_v_ghia.csv"};
    for (i = 0; i < 6; i++)
        remove(names[i]);
    rmdir("output");
    if (chdir(cwd) == 0) rmdir(tmpl);
    freem(&a);
    freem(&u);
    freem(&v);
}

// ---------------------------------------------------------------------------
// CPU tests
// ---------------------------------------------------------------------------

// Each derivative operator must differentiate along its own axis. The test
// functions are polynomials every row of the operator (boundary rows
// included) differentiates exactly, and their derivatives along the other axis
// are different, so swapping x and y would show up as an O(1) error.
static void test_cpu_operator_axes(int nx, int ny)
{
    int i, j, k, N = nx * ny;
    char name[96];
    problem p;
    problem_init(&p, nx, ny, 2, 3, 0.002, 1E-3, NULL);
    double *f = (double *)malloc(N * sizeof(double));
    double *d = (double *)malloc(N * sizeof(double));
    double *e = (double *)malloc(N * sizeof(double));
    double err;

    printf("CPU: derivative operators act along x and y, %dx%d grid\n", nx, ny);
    // DX of x(1 + y^2) is 1 + y^2; DY of y(1 + x^2) is 1 + x^2;
    // DX2 of x^2(1 + y) is 2(1 + y); DY2 of y^2(1 + x) is 2(1 + x)
    for (int op = 0; op < 4; op++)
    {
        for (i = 0; i < ny; i++)
            for (j = 0; j < nx; j++)
            {
                double x = j * p.cfg.dx, y = i * p.cfg.dy;
                k = i * nx + j;
                if (op == 0)
                {
                    f[k] = x * (1 + y * y);
                    e[k] = 1 + y * y;
                }
                if (op == 1)
                {
                    f[k] = y * (1 + x * x);
                    e[k] = 1 + x * x;
                }
                if (op == 2)
                {
                    f[k] = x * x * (1 + y);
                    e[k] = 2 * (1 + y);
                }
                if (op == 3)
                {
                    f[k] = y * y * (1 + x);
                    e[k] = 2 * (1 + x);
                }
            }
        spmv(op == 0 ? p.DX : op == 1 ? p.DY
                          : op == 2   ? p.DX2
                                      : p.DY2,
             f, d);
        err = rel_diff(d, e, N);
        snprintf(name, sizeof(name), "%s exact on its test polynomial", op == 0 ? "DX" : op == 1 ? "DY"
                                                                                     : op == 2   ? "DX2"
                                                                                                 : "DY2");
        check(name, err, 1E-9);
    }

    free(f);
    free(d);
    free(e);
    problem_free(&p);
}

// A sine mode that vanishes on the wall nodes is an eigenvector of the discrete
// Laplacian, so the FFT solver must return it divided by its eigenvalue.
// Fields are ny rows (y) of nx values (x). Sizes 8, 9, 10 and 17 cover a
// transform pass with fewer than one batch of 8 rows, exactly one batch, and
// a full batch plus a short one.
static void test_cpu_poisson_fft(int nx, int ny)
{
    int i, j, p = 3 < ny - 2 ? 3 : 1, q = 5 < nx - 2 ? 5 : 2;
    char name[96];
    double dx = 1.0 / (nx - 1), dy = 1.0 / (ny - 1);
    mtrx f = initm(ny, nx), psi = initm(ny, nx), expected = initm(ny, nx);
    double lambda = (2.0 * cos(PI * p / (double)(ny - 1)) - 2.0) / (dy * dy) + (2.0 * cos(PI * q / (double)(nx - 1)) - 2.0) / (dx * dx);

    printf("CPU: FFT Poisson solver against an exact eigenmode, %dx%d grid\n", nx, ny);
    for (i = 1; i < ny - 1; i++)
        for (j = 1; j < nx - 1; j++)
        {
            MAt(f, i, j) = sin(PI * i * p / (double)(ny - 1)) * sin(PI * j * q / (double)(nx - 1));
            MAt(expected, i, j) = MAt(f, i, j) / lambda;
        }

    fft_solver *fs = fft_setup(nx, ny);
    poisson_FFT(fs, f, psi, dx, dy);
    fft_cleanup(fs);
    snprintf(name, sizeof(name), "psi vs eigenmode / eigenvalue, %dx%d", nx, ny);
    check(name, rel_diff(psi.M, expected.M, nx * ny), 1E-10);

    freem(&f);
    freem(&psi);
    freem(&expected);
}

// Largest |value| on the wall nodes
static double wall_max(mtrx a)
{
    int i, j;
    double m = 0.0;
    for (i = 0; i < a.m; i++)
        for (j = 0; j < a.n; j++)
            if ((i == 0 || j == 0 || i == a.m - 1 || j == a.n - 1) && fabs(MAt(a, i, j)) > m)
                m = fabs(MAt(a, i, j));
    return m;
}

// Largest residual of the 5-point Laplacian at the interior points. Row
// neighbours (i +- 1) are y-neighbours, column neighbours (j +- 1) x-neighbours.
static double poisson_residual(mtrx f, mtrx psi, double dx, double dy)
{
    int i, j;
    double r, rmax = 0.0;

    for (i = 1; i < f.m - 1; i++)
        for (j = 1; j < f.n - 1; j++)
        {
            r = (MAt(psi, i + 1, j) - 2.0 * MAt(psi, i, j) + MAt(psi, i - 1, j)) / (dy * dy) + (MAt(psi, i, j + 1) - 2.0 * MAt(psi, i, j) + MAt(psi, i, j - 1)) / (dx * dx) - MAt(f, i, j);
            if (isnan(r)) return NAN;
            if (fabs(r) > rmax) rmax = fabs(r);
        }
    return rmax;
}

// All three solvers must solve the same discrete problem: once the iterative
// ones have converged they agree with the direct one and satisfy the 5-point
// equations, and all three put psi = 0 on the wall nodes.
static void test_cpu_poisson_agree(int nx, int ny)
{
    double dx = 1.0 / (nx - 1), dy = 1.0 / (ny - 1);
    double beta = sor_beta(nx, ny, dx, dy);
    mtrx f = initm(ny, nx), fft = initm(ny, nx), sor = initm(ny, nx), gs = initm(ny, nx);
    mtrx scratch = initm(ny, nx);

    printf("CPU: FFT, SOR and Gauss-Seidel solve the same problem, %dx%d grid\n", nx, ny);
    fill_pseudo_random(f.M, nx * ny, 5u);

    fft_solver *fs = fft_setup(nx, ny);
    poisson_FFT(fs, f, fft, dx, dy);
    fft_cleanup(fs);
    poisson_SOR(f, sor, scratch, dx, dy, 400000, 1E-13, beta);
    poisson(f, gs, scratch, dx, dy, 400000, 1E-13);

    check("FFT residual", poisson_residual(f, fft, dx, dy), 1E-8);
    check("SOR vs FFT", rel_diff(sor.M, fft.M, nx * ny), 1E-9);
    check("Gauss-Seidel vs FFT", rel_diff(gs.M, fft.M, nx * ny), 1E-9);
    check("psi on the wall nodes, all three solvers",
          wall_max(fft) + wall_max(sor) + wall_max(gs), 0.0);

    freem(&f);
    freem(&fft);
    freem(&sor);
    freem(&gs);
    freem(&scratch);
}

// On a strongly anisotropic grid (dx != dy) the SOR parameter must account
// for the spacings: the isotropic formula took about 5x the iterations.
static void test_cpu_sor_anisotropic(int nx, int ny, int max_iterations)
{
    double dx = 1.0 / (nx - 1), dy = 1.0 / (ny - 1);
    mtrx f = initm(ny, nx), psi = initm(ny, nx), scratch = initm(ny, nx);
    char name[96];

    printf("CPU: SOR iterations on an anisotropic %dx%d grid\n", nx, ny);
    fill_pseudo_random(f.M, nx * ny, 9u);
    int k = poisson_SOR(f, psi, scratch, dx, dy, 100000, 1E-3, sor_beta(nx, ny, dx, dy));
    snprintf(name, sizeof(name), "SOR iterations at the shipped tolerance (%d)", k);
    check(name, k, max_iterations);
    freem(&f);
    freem(&psi);
    freem(&scratch);
}

static void test_cpu_poisson_iterative(void)
{
    int n = 24;
    double dx = 1.0 / (n - 1);
    double beta = sor_beta(n, n, dx, dx);
    mtrx f = initm(n, n), psi = initm(n, n), scratch = initm(n, n);

    printf("CPU: Gauss-Seidel and SOR residuals\n");
    fill_pseudo_random(f.M, n * n, 7u);

    poisson(f, psi, scratch, dx, dx, 200000, 1E-9);
    check("Gauss-Seidel residual", poisson_residual(f, psi, dx, dx), 1E-4);
    poisson_SOR(f, psi, scratch, dx, dx, 200000, 1E-9, beta);
    check("SOR residual", poisson_residual(f, psi, dx, dx), 1E-4);

    freem(&f);
    freem(&psi);
    freem(&scratch);
}

// Short lid-driven cavity run: the fields must stay finite, the lid must
// have spun up the flow, and the velocity field must be divergence-free.
static void test_cpu_step(int nx, int ny, int time_scheme, const char *label)
{
    int t, N = nx * ny;
    double cmax, cmin;
    char name[96];
    problem p;
    problem_init(&p, nx, ny, time_scheme, 3, 0.002, 1E-3, NULL);

    printf("CPU: 50 steps of the lid-driven cavity, %s + FFT, %dx%d grid\n", label, nx, ny);
    for (t = 0; t < 50; t++)
        step(p.w, p.u, p.v, &p.ctx);
    continuity_range(&p, &cmax, &cmin);

    snprintf(name, sizeof(name), "%s: fields finite (max |w|)", label);
    check(name, max_abs(p.w.M, N), 1E6);
    snprintf(name, sizeof(name), "%s: lid drives the flow (-max |u|)", label);
    check(name, -max_abs(p.u.M, N), -1E-3);
    snprintf(name, sizeof(name), "%s: max |du/dx + dv/dy|", label);
    check(name, fmax(fabs(cmax), fabs(cmin)), 1E-9);

    problem_free(&p);
}

// The time-step limit must be sharp: just below it a run stays bounded, just
// above it the run diverges.
static void test_cpu_stability_limit(int time_scheme, const char *label)
{
    int t, n = 32;
    double frac[2] = {0.95, 1.05};
    char name[96];
    int f;

    printf("CPU: stability limit of %s, %dx%d grid\n", label, n, n);
    for (f = 0; f < 2; f++)
    {
        problem p;
        smtrx d2 = SDiff2(n, 6, 1.0 / (n - 1));
        double limit = max_stable_dt(&d2, &d2, 100.0, time_scheme); // problem_init()'s Re
        freesm(d2);
        double wmax = 0.0;

        problem_init(&p, n, n, time_scheme, 3, frac[f] * limit, 1E-3, NULL);
        for (t = 0; t < 3000 && !(wmax > 1E6); t++)
        {
            step(p.w, p.u, p.v, &p.ctx);
            wmax = max_abs(p.w.M, n * n);
        }
        if (f == 0)
        {
            snprintf(name, sizeof(name), "%s at 0.95 x limit: max |w| after 3000 steps", label);
            check(name, wmax, 1E3);
        }
        else
        {
            // A NaN counts as diverged
            snprintf(name, sizeof(name), "%s at 1.05 x limit: diverges (-max |w|)", label);
            check(name, isnan(wmax) ? -INFINITY : -wmax, -1E6);
        }
        problem_free(&p);
    }
}

// The other limits the start-up check uses
static void test_cpu_time_step_limits(void)
{
    int k, n = 64;
    double dx = 1.0 / n;
    smtrx d2 = SDiff2(n, 6, dx);
    char name[96];

    printf("CPU: time-step limits and the suggested dt\n");

    // Forward Euler's centered-advection limit is 2 nu / u^2 = 2 / (Re u^2).
    // At Re = 1000 it is the binding one, far below the viscous limit.
    check("Euler advection limit at Re = 1000, u = 1 is 0.002",
          fabs(euler_advection_dt(1000.0, 1.0) - 0.002), 1E-15);
    check("Euler advection limit at Re = 100, u = 2 is 0.005",
          fabs(euler_advection_dt(100.0, 2.0) - 0.005), 1E-15);
    check("Re = 1000: advection limit binds (below the viscous one)",
          euler_advection_dt(1000.0, 1.0) >= max_stable_dt(&d2, &d2, 1000.0, 1), 0.0);
    check("no advection limit when the walls are at rest", isinf(euler_advection_dt(100.0, 0.0)) ? 0.0 : 1.0, 0.0);

    // What main() refuses and what it suggests. Euler at Re = 1000 on 64^2:
    // Courant 1/64 < viscous, so the run is refused above 1/64, but the
    // suggestion is the advection limit 0.002, which the warning also uses
    dt_limits l = time_step_limits(&d2, &d2, dx, 1000.0, 1.0, 1.0, 1);
    check("Euler, Re = 1000: refused above min(Courant, viscous)",
          fabs(l.accept - fmin(dx, max_stable_dt(&d2, &d2, 1000.0, 1))), 1E-15);
    check("Euler, Re = 1000: Courant binds the refusal", fabs(l.accept - dx), 1E-15);
    check("Euler, Re = 1000: suggestion is the advection limit", fabs(l.suggest - 0.002), 1E-15);
    l = time_step_limits(&d2, &d2, dx, 1000.0, 1.0, 1.0, 2);
    check("RK4: no advection limit, suggestion = refusal threshold",
          (isinf(l.advection) ? 0.0 : 1.0) + fabs(l.suggest - l.accept), 0.0);
    l = time_step_limits(&d2, &d2, dx, 100.0, 1.0, 1.0, 1);
    check("Euler, Re = 100: viscous binds, suggestion = refusal threshold",
          fabs(l.accept - max_stable_dt(&d2, &d2, 100.0, 1)) + fabs(l.suggest - l.accept), 1E-15);

    // The dt printed as a suggestion must itself be accepted: rounded down,
    // and still within 1% of the limit
    double limits[] = {0.0056250, 0.0040391, 0.0158730, 0.0058050, 1.0, 0.001, 9.9999, 0.00999999,
                       max_stable_dt(&d2, &d2, 100.0, 1), max_stable_dt(&d2, &d2, 100.0, 2)};
    double worst = 0.0;
    for (k = 0; k < (int)(sizeof(limits) / sizeof(limits[0])); k++)
    {
        char text[32];
        snprintf(text, sizeof(text), "%.3g", round_down_3(limits[k]));
        double back = strtod(text, NULL);
        if (back > limits[k] || back < 0.99 * limits[k]) worst = fmax(worst, fabs(back - limits[k]) / limits[k]);
        if (back > limits[k]) worst = INFINITY;
    }
    snprintf(name, sizeof(name), "suggested dt <= limit and within 1%%, %d cases", k);
    check(name, worst, 0.0);
    freesm(d2);
}

// ---------------------------------------------------------------------------
// OpenMP
// ---------------------------------------------------------------------------

#ifdef _OPENMP
#include <omp.h>

// Run `steps` steps with `threads` threads and return w, u, v in out[0..2]
static void run_threads(int n, int scheme, int poisson_type, int steps, int threads, mtrx out[3])
{
    int t, saved = omp_get_max_threads();
    problem p;
    problem_init(&p, n, n, scheme, poisson_type, 0.002, 1E-3, NULL);
    omp_set_num_threads(threads);
    for (t = 0; t < steps; t++)
        step(p.w, p.u, p.v, &p.ctx);
    omp_set_num_threads(saved);
    out[0] = initm(n, n);
    out[1] = initm(n, n);
    out[2] = initm(n, n);
    mtrxcpy(out[0], p.w);
    mtrxcpy(out[1], p.u);
    mtrxcpy(out[2], p.v);
    problem_free(&p);
}

// Above OMP_MIN_WORK every parallel loop and the batched transforms run on
// several threads. The results must not depend on how many: with the FFT
// solver they are bitwise those of one thread, and so are the red-black
// Gauss-Seidel/SOR sweeps.
static void test_openmp_thread_count(int scheme, int poisson_type, const char *label)
{
    int k, n = 65, threads = omp_get_num_procs() > 4 ? 4 : (omp_get_num_procs() > 1 ? omp_get_num_procs() : 2);
    mtrx one[3], many[3];
    char name[96];

    printf("OpenMP: 1 vs %d threads, %s, %dx%d grid (above OMP_MIN_WORK = %d points)\n",
           threads, label, n, n, OMP_MIN_WORK);
    run_threads(n, scheme, poisson_type, 10, 1, one);
    run_threads(n, scheme, poisson_type, 10, threads, many);
    double diff = 0.0;
    for (k = 0; k < 3; k++)
    {
        diff += memcmp(one[k].M, many[k].M, (size_t)n * n * sizeof(double)) != 0;
        freem(&one[k]);
        freem(&many[k]);
    }
    snprintf(name, sizeof(name), "%s: w, u, v bitwise identical", label);
    check(name, diff, 0.0);
}

// Path of this executable, or "" where /proc/self/exe does not exist
static const char *self_path(void)
{
    static char self[4096];
    ssize_t len = readlink("/proc/self/exe", self, sizeof(self) - 1);
    self[len > 0 ? len : 0] = '\0';
    return self;
}

// The thread count default_threads() picks, in a fresh process started with
// the given environment (OpenMP reads its variables at start-up). -1 if the
// child fails.
static int child_threads(const char *env)
{
    char cmd[4400], out[32];
    int threads = -1;
    FILE *f;

    snprintf(cmd, sizeof(cmd),
             "env -u OMP_NUM_THREADS -u OMP_PROC_BIND -u OMP_PLACES -u GOMP_CPU_AFFINITY %s '%s' --default-threads",
             env, self_path());
    if (!(f = popen(cmd, "r"))) return -1;
    if (fgets(out, sizeof(out), f))
    {
        char *end;
        long t = strtol(out, &end, 10);
        if (end != out && t > 0 && t < 100000) threads = (int)t;
    }
    pclose(f);
    return threads;
}

// One thread per physical core by default, but a thread placement or count
// the user sets is respected. With a placement set, the runtime binds the
// initial thread to one place before main(), so counting cores from its
// affinity mask would give 1.
static void test_openmp_default_threads(void)
{
    int procs = omp_get_num_procs(), cores = physical_cores();
    int expect = cores > 0 && cores < procs ? cores : procs;
    char name[96];

    printf("OpenMP: default thread count (%d CPUs, %d physical cores)\n", procs, cores);
    // The children inherit this thread's CPU mask, which the runtime has
    // already narrowed to one place if a placement is set here. The core count
    // and the child processes need Linux.
    // Only these reasons skip the checks; a child that fails is a failure.
    if (omp_get_proc_bind() != omp_proc_bind_false || cores == 0 || self_path()[0] == '\0')
    {
        printf("  SKIPPED: %s\n", omp_get_proc_bind() != omp_proc_bind_false
                                      ? "this process has a thread placement (OMP_PROC_BIND, OMP_PLACES "
                                        "or GOMP_CPU_AFFINITY)"
                                      : "needs Linux sysfs and /proc/self/exe");
        n_skipped_other += 4;
        return;
    }
    snprintf(name, sizeof(name), "no OpenMP settings: %d threads", expect);
    check(name, child_threads("") != expect, 0.0);
    snprintf(name, sizeof(name), "OMP_PROC_BIND=close: all %d CPUs", procs);
    check(name, child_threads("OMP_PROC_BIND=close") != procs, 0.0);
    snprintf(name, sizeof(name), "OMP_PLACES=cores: all %d CPUs", procs);
    check(name, child_threads("OMP_PLACES=cores") != procs, 0.0);
    check("OMP_NUM_THREADS=3: 3 threads", child_threads("OMP_NUM_THREADS=3") != 3, 0.0);
}
#endif

// Every solver copies the configuration when it is created. Changing the
// caller's struct afterwards (dt and the lid speed here) must have no effect,
// on any of them: the CPU workspace, the CPU backend and the GPU solver.
static void test_config_copied(void)
{
    int t, k, nx = 32, ny = 32, N = nx * ny;
    problem p;
    problem_init(&p, nx, ny, 2, 3, 0.002, 1E-3, NULL);
    solver_config caller = p.cfg;
    backend *cpu = backend_create(&caller, 0);
    mtrx *bw;
    double diff = 0.0;
#ifdef USE_CUDA
    gpu_solver *gpu = gpu_init(&caller);
    mtrx gw = initm(ny, nx);
#endif

    printf("Configuration is copied when a solver is created\n");
    for (k = 0; k < 2; k++)
    {
        for (t = 0; t < 10; t++)
        {
            step(p.w, p.u, p.v, &p.ctx);
            backend_step(cpu);
#ifdef USE_CUDA
            if (gpu) gpu_step(gpu);
#endif
        }
        // After the first 10 steps, change everything the solvers read each step
        caller.dt = 0.001;
        caller.bc.u[3] = -1.0;
        p.cfg = caller;
    }
    backend_fields(cpu, NULL, NULL, &bw);
    diff = rel_diff(bw->M, p.w.M, N);
    check("CPU backend and step() both ignore the change (w, bitwise)", diff, 0.0);
#ifdef USE_CUDA
    if (gpu)
    {
        gpu_get_fields(gpu, NULL, NULL, &gw);
        check("GPU solver ignores it too (w vs CPU)", rel_diff(gw.M, p.w.M, N), 1E-11);
        gpu_free(gpu);
    }
    else
        n_skipped++;
    freem(&gw);
#endif
    // Twenty steps at the original dt and lid speed, not ten and ten
    problem q;
    problem_init(&q, nx, ny, 2, 3, 0.002, 1E-3, NULL);
    for (t = 0; t < 20; t++)
        step(q.w, q.u, q.v, &q.ctx);
    check("same result as twenty unchanged steps", rel_diff(p.w.M, q.w.M, N), 0.0);
    problem_free(&q);
    backend_free(cpu);
    problem_free(&p);
}

// Driving the solver through the backend interface must give exactly what
// calling step() directly gives. In CUDA builds the GPU backend must also
// agree with the CPU one.
static void test_backend(void)
{
    int t, nx = 40, ny = 24, N = nx * ny, steps = 20;
    double cmax, cmin, cmax_ref, cmin_ref;
    problem p;
    problem_init(&p, nx, ny, 2, 3, 0.002, 1E-3, NULL);
    mtrx u = initm(ny, nx), v = initm(ny, nx), w = initm(ny, nx);
    backend *b = backend_create(&p.cfg, 0);

    printf("Backend interface, %dx%d grid\n", nx, ny);
    check("CPU backend chosen when the GPU is not wanted", strcmp(backend_name(b), "CPU") != 0, 0.0);
    backend_set_fields(b, &p.u, &p.v, &p.w);
    for (t = 0; t < steps; t++)
    {
        backend_step(b);
        step(p.w, p.u, p.v, &p.ctx);
    }
    backend_get_fields(b, &u, &v, &w);
    backend_continuity(b, &cmax, &cmin);
    continuity_range(&p, &cmax_ref, &cmin_ref);
    check("CPU backend vs step(): u, v, w",
          rel_diff(u.M, p.u.M, N) + rel_diff(v.M, p.v.M, N) + rel_diff(w.M, p.w.M, N), 0.0);
    check("CPU backend vs step(): continuity", fabs(cmax - cmax_ref) + fabs(cmin - cmin_ref), 0.0);
    backend_free(b);

#ifdef USE_CUDA
    b = backend_create(&p.cfg, 1);
    if (strcmp(backend_name(b), "CUDA") == 0)
    {
        mtrx u0 = initm(ny, nx), v0 = initm(ny, nx), w0 = initm(ny, nx);
        backend_set_fields(b, &u0, &v0, &w0);
        for (t = 0; t < steps; t++)
            backend_step(b);
        backend_get_fields(b, &u, &v, &w);
        check("GPU backend vs CPU: w", rel_diff(w.M, p.w.M, N), 1E-11);
        check("GPU backend vs CPU: u", rel_diff(u.M, p.u.M, N), 1E-11);
        freem(&u0);
        freem(&v0);
        freem(&w0);
    }
    else
        n_skipped++;
    backend_free(b);
#endif

    freem(&u);
    freem(&v);
    freem(&w);
    problem_free(&p);
}

// ---------------------------------------------------------------------------
// GPU tests — every check compares the device result with the CPU one
// ---------------------------------------------------------------------------

#ifdef USE_CUDA

static gpu_solver *gpu_for(problem *p)
{
    gpu_solver *g = gpu_init(&p->cfg);
    if (!g)
    {
        printf("** Error: CUDA device disappeared during the tests **\n");
        exit(1);
    }
    return g;
}

static void test_gpu_spmv(int nx, int ny)
{
    int op, N = nx * ny;
    char name[96];
    const char *op_name[] = {"DX", "DY", "DX2", "DY2"};
    problem p;
    problem_init(&p, nx, ny, 2, 3, 0.002, 1E-3, NULL);
    gpu_solver *g = gpu_for(&p);
    smtrx *ops[] = {&p.DX, &p.DY, &p.DX2, &p.DY2};
    double *x = (double *)malloc(N * sizeof(double));
    double *y_cpu = (double *)malloc(N * sizeof(double));
    double *y_gpu = (double *)malloc(N * sizeof(double));

    printf("GPU: CSR SpMV, %dx%d grid\n", nx, ny);
    fill_pseudo_random(x, N, 11u);
    for (op = 0; op < 4; op++)
    {
        spmv(*ops[op], x, y_cpu);
        gpu_spmv(g, op, x, y_gpu);
        snprintf(name, sizeof(name), "%s*x vs CPU", op_name[op]);
        check(name, rel_diff(y_gpu, y_cpu, N), 1E-12);
    }

    free(x);
    free(y_cpu);
    free(y_gpu);
    gpu_free(g);
    problem_free(&p);
}

static void test_gpu_poisson(int nx, int ny, int poisson_type, const char *label, double limit)
{
    int N = nx * ny;
    char name[96];
    problem p;
    problem_init(&p, nx, ny, 2, poisson_type, 0.002, 1E-10, NULL);
    gpu_solver *g = gpu_for(&p);
    mtrx w = initm(ny, nx), f = initm(ny, nx), psi_cpu = initm(ny, nx), psi_gpu = initm(ny, nx);

    printf("GPU: %s Poisson solver, %dx%d grid\n", label, nx, ny);
    fill_pseudo_random(w.M, N, 23u);

    // CPU solvers take the right-hand side f = -w
    mtrxcpy(f, w);
    invsig(f);
    if (poisson_type == 1)
        poisson(f, psi_cpu, p.ctx.psi_scratch, p.cfg.dx, p.cfg.dy,
                p.cfg.poisson_max_it, p.cfg.poisson_tol);
    else if (poisson_type == 2)
        poisson_SOR(f, psi_cpu, p.ctx.psi_scratch, p.cfg.dx, p.cfg.dy,
                    p.cfg.poisson_max_it, p.cfg.poisson_tol, p.cfg.beta);
    else
        poisson_FFT(p.ctx.fft, f, psi_cpu, p.cfg.dx, p.cfg.dy);

    gpu_poisson(g, w.M, psi_gpu.M);
    snprintf(name, sizeof(name), "%s: psi vs CPU", label);
    check(name, rel_diff(psi_gpu.M, psi_cpu.M, N), limit);
    if (poisson_type != 3)
    {
        snprintf(name, sizeof(name), "%s: residual", label);
        check(name, poisson_residual(f, psi_gpu, p.cfg.dx, p.cfg.dy), 1E-4);
    }

    freem(&w);
    freem(&f);
    freem(&psi_cpu);
    freem(&psi_gpu);
    gpu_free(g);
    problem_free(&p);
}

// Run the same case on both backends and compare the fields afterwards.
// bc == NULL keeps the lid-driven cavity.
static void test_gpu_step(int nx, int ny, int steps, double dt, int time_scheme, int poisson_type,
                          double poisson_tol, const wall_bc *bc, const char *label, double limit)
{
    int t, N = nx * ny;
    char name[96];
    problem p;
    problem_init(&p, nx, ny, time_scheme, poisson_type, dt, poisson_tol, bc);
    gpu_solver *g = gpu_for(&p);
    mtrx u = initm(ny, nx), v = initm(ny, nx), w = initm(ny, nx);

    printf("GPU: %d steps, %s, %dx%d grid\n", steps, label, nx, ny);
    gpu_set_fields(g, &p.u, &p.v, &p.w);
    for (t = 0; t < steps; t++)
    {
        step(p.w, p.u, p.v, &p.ctx);
        gpu_step(g);
    }
    gpu_get_fields(g, &u, &v, &w);

    snprintf(name, sizeof(name), "%s: w vs CPU", label);
    check(name, rel_diff(w.M, p.w.M, N), limit);
    snprintf(name, sizeof(name), "%s: u vs CPU", label);
    check(name, rel_diff(u.M, p.u.M, N), limit);
    snprintf(name, sizeof(name), "%s: v vs CPU", label);
    check(name, rel_diff(v.M, p.v.M, N), limit);

    freem(&u);
    freem(&v);
    freem(&w);
    gpu_free(g);
    problem_free(&p);
}

// A run with a vorticity source (the manufactured solution of mms.h), started
// at t0 = 0.3, must match the CPU: the GPU fills the source on the host at
// each stage's time and adds it
static void test_gpu_source(int nx, int ny, int steps, int time_scheme, const char *label)
{
    int t, N = nx * ny;
    char name[96];
    mms_case c = {1.0, 1.0, 100.0, 1.0 / (nx - 1), 1.0 / (ny - 1), 0};
    wall_bc walls = {{0., 0., 0., 0.}, {0., 0., 0., 0.}};
    problem p;
    problem_init_ext(&p, nx, ny, 1.0, 1.0, 0, 6, c.Re, time_scheme, 3, 0.002, 1E-10, &walls, 0.3, mms_source, &c, 2, 0, 2);
    gpu_solver *g = gpu_for(&p);
    mtrx u = initm(ny, nx), v = initm(ny, nx), w = initm(ny, nx);

    printf("GPU: %d steps with a vorticity source, %s, %dx%d grid\n", steps, label, nx, ny);
    mms_exact(&c, 0.3, &p.w, &p.u, &p.v, NULL);
    gpu_set_fields(g, &p.u, &p.v, &p.w);
    for (t = 0; t < steps; t++)
    {
        step(p.w, p.u, p.v, &p.ctx);
        gpu_step(g);
    }
    gpu_get_fields(g, &u, &v, &w);

    snprintf(name, sizeof(name), "%s with a source: w vs CPU", label);
    check(name, rel_diff(w.M, p.w.M, N), 1E-11);
    snprintf(name, sizeof(name), "%s with a source: u, v vs CPU", label);
    check(name, rel_diff(u.M, p.u.M, N) + rel_diff(v.M, p.v.M, N), 1E-11);

    freem(&u);
    freem(&v);
    freem(&w);
    gpu_free(g);
    problem_free(&p);
}

// A periodic run (the periodic manufactured solution, with its source) must
// match the CPU: periodic Poisson solver by cuFFT, no wall updates
static void test_gpu_periodic(int nx, int ny, int steps, int time_scheme, const char *label)
{
    int t, N = nx * ny;
    char name[96];
    mms_case c = {2.0, 1.0, 100.0, 2.0 / nx, 1.0 / ny, 1};
    problem p;
    problem_init_ext(&p, nx, ny, 2.0, 1.0, 1, 6, c.Re, time_scheme, 3, 0.002, 1E-10, NULL, 0.0, mms_source, &c, 2, 0, 2);
    gpu_solver *g = gpu_for(&p);
    mtrx u = initm(ny, nx), v = initm(ny, nx), w = initm(ny, nx);

    printf("GPU: %d steps on a periodic grid, %s, %dx%d grid\n", steps, label, nx, ny);
    mms_exact(&c, 0.0, &p.w, &p.u, &p.v, NULL);
    gpu_set_fields(g, &p.u, &p.v, &p.w);
    for (t = 0; t < steps; t++)
    {
        step(p.w, p.u, p.v, &p.ctx);
        gpu_step(g);
    }
    gpu_get_fields(g, &u, &v, &w);

    snprintf(name, sizeof(name), "%s, periodic: w vs CPU", label);
    check(name, rel_diff(w.M, p.w.M, N), 1E-11);
    snprintf(name, sizeof(name), "%s, periodic: u, v vs CPU", label);
    check(name, rel_diff(u.M, p.u.M, N) + rel_diff(v.M, p.v.M, N), 1E-11);

    freem(&u);
    freem(&v);
    freem(&w);
    gpu_free(g);
    problem_free(&p);
}

// The compact Poisson operator and the wall closure from psi on the GPU:
// a lid-driven cavity run must match the CPU, and the compact Poisson solve
// on its own too
static void test_gpu_closures(int nx, int ny, int steps, int time_scheme, const char *label)
{
    int t, N = nx * ny;
    char name[96];
    problem p;
    problem_init_ext(&p, nx, ny, 1.0, 1.0, 0, 6, 100., time_scheme, 3, 0.002, 1E-10, NULL, 0.0, NULL, NULL, 4, 1, 4);
    gpu_solver *g = gpu_for(&p);
    mtrx u = initm(ny, nx), v = initm(ny, nx), w = initm(ny, nx), f = initm(ny, nx), psi = initm(ny, nx);

    printf("GPU: %d steps, compact Poisson + wall closure from psi + fourth-order velocity rows, %s, %dx%d grid\n",
           steps, label, nx, ny);
    gpu_set_fields(g, &p.u, &p.v, &p.w);
    for (t = 0; t < steps; t++)
    {
        step(p.w, p.u, p.v, &p.ctx);
        gpu_step(g);
    }
    gpu_get_fields(g, &u, &v, &w);
    snprintf(name, sizeof(name), "%s, compact + psi closure: w vs CPU", label);
    check(name, rel_diff(w.M, p.w.M, N), 1E-11);
    snprintf(name, sizeof(name), "%s, compact + psi closure: u, v vs CPU", label);
    check(name, rel_diff(u.M, p.u.M, N) + rel_diff(v.M, p.v.M, N), 1E-11);

    // The Poisson solve alone, on the final w: the GPU takes w, the CPU -w
    mtrxcpy(f, p.w);
    invsig(f);
    poisson_FFT_order(p.ctx.fft, f, psi, p.cfg.dx, p.cfg.dy, 4);
    gpu_poisson(g, p.w.M, u.M);
    snprintf(name, sizeof(name), "%s, compact Poisson: psi vs CPU", label);
    check(name, rel_diff(u.M, psi.M, N), 1E-12);

    freem(&u);
    freem(&v);
    freem(&w);
    freem(&f);
    freem(&psi);
    gpu_free(g);
    problem_free(&p);
}

// The GPU integrals match compute_integrals() on the same fields, with walls
// (lid-driven cavity, after some steps, so the wall nodes carry the wall
// velocity) and on a periodic grid
static void test_gpu_integrals(void)
{
    int t, periodic;
    char name[96];

    printf("GPU: energy, enstrophy and palinstrophy\n");
    for (periodic = 0; periodic <= 1; periodic++)
    {
        int nx = periodic ? 48 : 40, ny = periodic ? 33 : 25;
        mms_case c = {1.0, 1.0, 100.0, 1.0 / nx, 1.0 / ny, 1};
        problem p;
        problem_init_ext(&p, nx, ny, 1.0, 1.0, periodic, 6, 100., 2, 3, 0.002, 1E-10, NULL, 0.0, NULL, NULL, 2, 0, 2);
        gpu_solver *g = gpu_for(&p);
        double E, Z, P;
        if (periodic) mms_exact(&c, 0.0, &p.w, &p.u, &p.v, NULL);
        gpu_set_fields(g, &p.u, &p.v, &p.w);
        for (t = 0; t < 20; t++)
        {
            step(p.w, p.u, p.v, &p.ctx);
            gpu_step(g);
        }
        flow_integrals f = compute_integrals(&p.cfg, p.u, p.v, p.w, p.ctx.k1.M, p.ctx.k2.M);
        double I;
        gpu_integrals(g, &E, &Z, &P, &I);
        snprintf(name, sizeof(name), "%s: E, Z, P vs CPU", periodic ? "periodic" : "walls");
        check(name, fabs(E - f.E) / f.E + fabs(Z - f.Z) / f.Z + fabs(P - f.P) / f.P, 1E-11);
        gpu_free(g);
        problem_free(&p);
    }
}

// Drag, hypodrag, hyperviscosity, Kolmogorov forcing and random kicks, with the
// skew-symmetric nonlinear term, on the GPU match the CPU (same seed, so the
// same kicks), with RK4 and Euler, and so does the energy input
static void test_gpu_forcing(int time_scheme, const char *label)
{
    int t, nx = 48, ny = 40, N = nx * ny;
    char name[96];
    problem p;
    problem_init_ext(&p, nx, ny, 1.0, 1.0, 1, 6, 1000., time_scheme, 3, 5E-4, 1E-10, NULL, 0.0, NULL, NULL, 2, 0, 2);
    rk4_free(&p.ctx);
    p.cfg.forcing.drag = 0.1;
    p.cfg.forcing.kolmogorov_amp = 0.5;
    p.cfg.forcing.kolmogorov_n = 3;
    p.cfg.forcing.random_rate = 0.2;
    p.cfg.forcing.random_kf = 5.0;
    p.cfg.forcing.random_dk = 1.0;
    p.cfg.forcing.random_seed = 7;
    p.cfg.forcing.hypodrag = 0.5;
    p.cfg.forcing.hyperviscosity = 1E-10;
    p.cfg.forcing.hyper_order = 2;
    p.cfg.advection = 1;
    p.ctx = rk4_alloc(&p.cfg);
    gpu_solver *g = gpu_for(&p);
    mtrx u = initm(ny, nx), v = initm(ny, nx), w = initm(ny, nx);
    double E, Z, P, I;

    printf("GPU: %s, skew-symmetric, with drag, hypodrag, hyperviscosity, Kolmogorov and random forcing, %dx%d "
           "periodic grid\n",
           label, nx, ny);
    gpu_set_fields(g, &p.u, &p.v, &p.w);
    for (t = 0; t < 40; t++)
    {
        step(p.w, p.u, p.v, &p.ctx);
        gpu_step(g);
    }
    gpu_get_fields(g, &u, &v, &w);
    snprintf(name, sizeof(name), "%s, forced: w vs CPU", label);
    check(name, rel_diff(w.M, p.w.M, N), 1E-11);
    snprintf(name, sizeof(name), "%s, forced: u, v vs CPU", label);
    check(name, rel_diff(u.M, p.u.M, N) + rel_diff(v.M, p.v.M, N), 1E-11);
    flow_integrals f = compute_integrals(&p.cfg, p.u, p.v, p.w, p.ctx.k1.M, p.ctx.k2.M);
    gpu_integrals(g, &E, &Z, &P, &I);
    snprintf(name, sizeof(name), "%s, forced: energy input I vs CPU", label);
    check(name, fabs(I - f.I) / fabs(f.I), 1E-11);

    freem(&u);
    freem(&v);
    freem(&w);
    gpu_free(g);
    problem_free(&p);
}

// Hyperviscosity strong enough to dominate near the cutoff: the GPU applies it
// in spectral space, the CPU by p sparse products, and the two agree
static void test_gpu_hyperviscosity(int p_order)
{
    int t, nx = 32, ny = 24, N = nx * ny;
    char name[96];
    problem p;
    problem_init_ext(&p, nx, ny, 1.0, 1.0, 1, 6, 1000., 2, 3, 2E-4, 1E-10, NULL, 0.0, NULL, NULL, 2, 0, 2);
    rk4_free(&p.ctx);
    p.cfg.forcing.hyper_order = p_order;
    // Damping rate 2/dt at the largest eigenvalue: a tenth of the stability limit
    smtrx dx2 = SDiff2_periodic(nx, 6, 1.0 / nx), dy2 = SDiff2_periodic(ny, 6, 1.0 / ny);
    p.cfg.forcing.hyperviscosity = 1.0;
    double rate = damping_rate_max(&dx2, &dy2, 1E300, &p.cfg.forcing);
    p.cfg.forcing.hyperviscosity = 2.0 / (p.cfg.dt * rate);
    freesm(dx2);
    freesm(dy2);
    p.ctx = rk4_alloc(&p.cfg);
    fill_pseudo_random(p.w.M, N, 11u);
    gpu_solver *g = gpu_for(&p);
    mtrx w = initm(ny, nx);
    double w0 = 0.0, w1 = 0.0;

    printf("GPU: hyperviscosity of order %d in spectral space, %dx%d periodic grid\n", p_order, nx, ny);
    for (t = 0; t < N; t++)
        w0 += p.w.M[t] * p.w.M[t];
    gpu_set_fields(g, &p.u, &p.v, &p.w);
    for (t = 0; t < 20; t++)
    {
        step(p.w, p.u, p.v, &p.ctx);
        gpu_step(g);
    }
    gpu_get_fields(g, NULL, NULL, &w);
    for (t = 0; t < N; t++)
        w1 += p.w.M[t] * p.w.M[t];
    snprintf(name, sizeof(name), "w vs CPU (<w^2> fell to %.2f of its start)", w1 / w0);
    check(name, rel_diff(w.M, p.w.M, N), 1E-12);
    freem(&w);
    gpu_free(g);
    problem_free(&p);
}

// The spectra, fluxes and dissipation computed on the device match
// spectra_all() on the host, column by column, with either nonlinear term and
// every damping term
static void test_gpu_spectra(int advection)
{
    static const char *cols[SPECTRA_COLUMNS] = {"E", "Z", "Pi_E", "Pi_Z", "D_E", "D_Z", "F_E", "F_Z"};
    int t, q, b, nx = 40, ny = 32, N = nx * ny;
    char name[96];
    problem p;
    problem_init_ext(&p, nx, ny, 1.0, 1.0, 1, 6, 500., 2, 3, 5E-4, 1E-10, NULL, 0.0, NULL, NULL, 2, 0, 2);
    rk4_free(&p.ctx);
    p.cfg.advection = advection;
    p.cfg.forcing.drag = 0.2;
    p.cfg.forcing.hypodrag = 3.0;
    p.cfg.forcing.hyperviscosity = 1E-12;
    p.cfg.forcing.hyper_order = 3;
    p.ctx = rk4_alloc(&p.cfg);
    fill_pseudo_random(p.w.M, N, 13u);
    for (t = 0; t < N; t++)
        p.w.M[t] += 0.3; // a mean vorticity, which only the drag damps
    gpu_solver *g = gpu_for(&p);
    spectra *sp = spectra_setup(&p.cfg);
    int B = spectra_bins(sp);
    double *dev = (double *)malloc((size_t)SPECTRA_COLUMNS * B * sizeof(double));
    double *host = (double *)malloc((size_t)SPECTRA_COLUMNS * B * sizeof(double));
    mtrx u = initm(ny, nx), v = initm(ny, nx), w = initm(ny, nx);

    printf("GPU: spectra and fluxes on the device, %s form, %dx%d periodic grid\n",
           advection ? "skew-symmetric" : "advective", nx, ny);
    gpu_set_fields(g, &p.u, &p.v, &p.w);
    for (t = 0; t < 5; t++)
        gpu_step(g);
    gpu_spectra(g, sp, dev);
    gpu_get_fields(g, &u, &v, &w);
    spectra_all(sp, u, v, w, host);
    for (q = 0; q < SPECTRA_COLUMNS; q++)
    {
        double diff = 0.0, peak = 0.0;
        for (b = 0; b < B; b++)
        {
            diff = fmax(diff, fabs(dev[q * B + b] - host[q * B + b]));
            peak = fmax(peak, fabs(host[q * B + b]));
        }
        snprintf(name, sizeof(name), "%s vs host (largest %.2e)", cols[q], peak);
        check(name, peak > 0.0 ? diff / peak : INFINITY, 1E-12);
    }
    freem(&u);
    freem(&v);
    freem(&w);
    free(dev);
    free(host);
    spectra_free(sp);
    gpu_free(g);
    problem_free(&p);
}

// Fields written to the device must come back unchanged, and the continuity
// diagnostic (a max and a min reduction) must match the CPU on a field that
// is far from divergence-free.
static void test_gpu_fields(int n)
{
    int k, N = n * n;
    double cmax_cpu, cmin_cpu, cmax_gpu, cmin_gpu;
    problem p;
    problem_init(&p, n, n, 2, 3, 0.002, 1E-3, NULL);
    gpu_solver *g = gpu_for(&p);
    mtrx u = initm(n, n), v = initm(n, n), w = initm(n, n);

    printf("GPU: field transfer and continuity diagnostic, %dx%d grid\n", n, n);
    fill_pseudo_random(p.u.M, N, 1u);
    fill_pseudo_random(p.v.M, N, 2u);
    fill_pseudo_random(p.w.M, N, 3u);
    // Put the extremes in the last rows, so a reduction that stops short of
    // the end of the array cannot get them right.
    for (k = N - 4 * n; k < N; k++)
    {
        p.u.M[k] *= 100.0;
        p.v.M[k] *= 100.0;
    }
    gpu_set_fields(g, &p.u, &p.v, &p.w);
    gpu_get_fields(g, &u, &v, &w);
    check("u, v, w identical after round trip",
          rel_diff(u.M, p.u.M, N) + rel_diff(v.M, p.v.M, N) + rel_diff(w.M, p.w.M, N), 0.0);

    continuity_range(&p, &cmax_cpu, &cmin_cpu);
    gpu_continuity(g, &cmax_gpu, &cmin_gpu);
    check("continuity max vs CPU", fabs(cmax_gpu - cmax_cpu) / fabs(cmax_cpu), 1E-12);
    check("continuity min vs CPU", fabs(cmin_gpu - cmin_cpu) / fabs(cmin_cpu), 1E-12);

    freem(&u);
    freem(&v);
    freem(&w);
    gpu_free(g);
    problem_free(&p);
}

// Run a GPU test, or count it as skipped when there is no device to run it on
#define GPU_TEST(call)   \
    do                   \
    {                    \
        if (have_gpu)    \
            (call);      \
        else             \
            n_skipped++; \
    } while (0)

static void run_gpu_tests(void)
{
    // A different velocity on every wall, so a wall mix-up cannot go unnoticed
    wall_bc four_walls = {{0.3, -0.2, 0.5, 1.0}, {0.1, -0.4, 0.2, -0.3}};

    // Probe for a device first
    problem probe;
    problem_init(&probe, 16, 16, 2, 3, 0.002, 1E-3, NULL);
    gpu_solver *g = gpu_init(&probe.cfg);
    int have_gpu = (g != NULL);
    problem_free(&probe);
    if (have_gpu)
        gpu_free(g);
    else
        printf("GPU: no usable CUDA device - GPU tests SKIPPED\n");

    // 300x300 is large enough for the reductions to span several passes
    GPU_TEST(test_gpu_fields(16));
    GPU_TEST(test_gpu_fields(300));
    GPU_TEST(test_gpu_spmv(32, 32));
    GPU_TEST(test_gpu_spmv(45, 45));
    GPU_TEST(test_gpu_spmv(300, 300));
    GPU_TEST(test_gpu_spmv(37, 21));
    GPU_TEST(test_gpu_poisson(32, 32, 3, "FFT", 1E-12));
    GPU_TEST(test_gpu_poisson(45, 45, 3, "FFT", 1E-12));
    GPU_TEST(test_gpu_poisson(33, 20, 3, "FFT", 1E-12));
    GPU_TEST(test_gpu_poisson(20, 33, 3, "FFT", 1E-12));
    GPU_TEST(test_gpu_poisson(26, 15, 2, "SOR", 1E-7));
    GPU_TEST(test_gpu_poisson(24, 24, 2, "SOR", 1E-7));
    GPU_TEST(test_gpu_poisson(24, 24, 1, "Gauss-Seidel", 1E-7));
    GPU_TEST(test_gpu_step(32, 32, 50, 0.002, 2, 3, 1E-10, NULL, "RK4 + FFT", 1E-11));
    GPU_TEST(test_gpu_step(32, 32, 50, 0.002, 1, 3, 1E-10, NULL, "Euler + FFT", 1E-11));
    GPU_TEST(test_gpu_step(45, 45, 50, 0.002, 2, 3, 1E-10, NULL, "RK4 + FFT", 1E-11));
    GPU_TEST(test_gpu_step(128, 128, 10, 0.0005, 2, 3, 1E-10, NULL, "RK4 + FFT", 1E-11));
    GPU_TEST(test_gpu_step(40, 24, 50, 0.002, 2, 3, 1E-10, NULL, "RK4 + FFT", 1E-11));
    GPU_TEST(test_gpu_step(24, 40, 50, 0.002, 1, 3, 1E-10, NULL, "Euler + FFT", 1E-11));
    GPU_TEST(test_gpu_source(33, 33, 50, 2, "RK4 + FFT"));
    GPU_TEST(test_gpu_periodic(48, 33, 50, 2, "RK4"));
    GPU_TEST(test_gpu_integrals());
    GPU_TEST(test_gpu_spectra(0));
    GPU_TEST(test_gpu_spectra(1));
    GPU_TEST(test_gpu_hyperviscosity(2));
    GPU_TEST(test_gpu_hyperviscosity(4));
    GPU_TEST(test_gpu_forcing(2, "RK4"));
    GPU_TEST(test_gpu_forcing(1, "Euler"));
    GPU_TEST(test_gpu_closures(40, 25, 30, 2, "RK4"));
    GPU_TEST(test_gpu_closures(32, 32, 30, 1, "Euler"));
    GPU_TEST(test_gpu_periodic(32, 32, 50, 1, "Euler"));
    GPU_TEST(test_gpu_source(40, 24, 50, 1, "Euler + FFT"));
    GPU_TEST(test_gpu_step(40, 24, 20, 0.002, 2, 3, 1E-10, &four_walls,
                           "RK4 + FFT, four moving walls", 1E-11));
    GPU_TEST(test_gpu_step(26, 15, 3, 0.002, 2, 2, 1E-10, NULL, "RK4 + SOR", 1E-7));
    GPU_TEST(test_gpu_step(32, 32, 20, 0.002, 2, 3, 1E-10, &four_walls,
                           "RK4 + FFT, four moving walls", 1E-11));

    // The iterative solvers are run to a tight tolerance here. The default
    // CPU build sweeps lexicographically and the GPU red-black, so the two
    // only agree once the iteration has converged.
    GPU_TEST(test_gpu_step(24, 24, 3, 0.002, 2, 2, 1E-10, NULL, "RK4 + SOR", 1E-7));
    GPU_TEST(test_gpu_step(24, 24, 3, 0.002, 1, 1, 1E-10, NULL, "Euler + Gauss-Seidel", 1E-7));
#ifdef _OPENMP
    // With OPENMP=1 the CPU sweeps red-black too, so the backends must also
    // agree at the shipped tolerance, where the iteration stops far from
    // convergence.
    GPU_TEST(test_gpu_step(64, 64, 10, 0.005, 2, 2, 1E-3, NULL,
                           "RK4 + SOR, shipped tolerance", 1E-9));
    GPU_TEST(test_gpu_step(24, 24, 10, 0.002, 1, 1, 1E-3, NULL,
                           "Euler + Gauss-Seidel, shipped tolerance", 1E-9));
#endif
}

#endif // USE_CUDA

int main(int argc, char **argv)
{
#ifdef _OPENMP
    // Used by test_openmp_default_threads()
    if (argc > 1 && strcmp(argv[1], "--default-threads") == 0)
    {
        default_threads();
        printf("%d\n", omp_get_max_threads());
        return 0;
    }
#else
    (void)argc;
    (void)argv;
#endif
    test_finitediff_exactness();
    test_finitediff_rows();
    test_linearalg();
    test_wall_bc();
    test_output();
    test_temporal_order();
    test_source_temporal_order();
    test_mms_spatial_order();
    test_order_ablation();
    test_compact_poisson();
    test_velocity_rows();
    test_wall_order_pairing();
    test_periodic_operators();
    test_periodic_poisson(32, 24, 6);
    test_periodic_poisson(15, 21, 4);
    test_periodic_poisson(16, 16, 2);
    test_periodic_mms_order();
    test_integrals_taylor_green();
    test_kolmogorov_laminar();
    test_drag_decay();
    test_random_initial_field();
    test_damping_stability();
    test_random_kick();
    test_budgets();
    test_cpu_operator_axes(13, 9);
    test_cpu_operator_axes(9, 13);
    test_cpu_poisson_fft(8, 8);
    test_cpu_poisson_fft(9, 9);
    test_cpu_poisson_fft(10, 10);
    test_cpu_poisson_fft(17, 17);
    test_cpu_poisson_fft(32, 32);
    test_cpu_poisson_fft(33, 20);
    test_cpu_poisson_fft(9, 17);
    test_cpu_poisson_agree(24, 24);
    test_cpu_poisson_agree(26, 15);
    test_cpu_poisson_iterative();
    test_cpu_sor_anisotropic(10, 200, 200);
    test_cpu_sor_anisotropic(200, 10, 200);
    test_cpu_step(32, 32, 2, "RK4");
    test_cpu_step(32, 32, 1, "Euler");
    test_cpu_step(40, 24, 2, "RK4");
    test_cpu_step(24, 40, 2, "RK4");
    test_cpu_stability_limit(1, "Euler");
    test_cpu_stability_limit(2, "RK4");
    test_cpu_time_step_limits();
#ifdef _OPENMP
    test_openmp_thread_count(2, 3, "RK4 + FFT");
    test_openmp_thread_count(2, 2, "RK4 + SOR");
    test_openmp_thread_count(1, 1, "Euler + Gauss-Seidel");
    test_openmp_default_threads();
#endif
    test_backend();
    test_config_copied();
#ifdef USE_CUDA
    run_gpu_tests();
#endif

    printf("\n%d checks, %d failed", n_checks, n_failed);
    if (n_skipped_other) printf(", %d OpenMP thread-count checks SKIPPED (see above)", n_skipped_other);
#ifdef USE_CUDA
    if (n_skipped)
        printf(", %d GPU tests SKIPPED (no usable CUDA device)", n_skipped);
#else
    printf(" (CPU only: the GPU tests need CUDA=1)");
#endif
    printf("\n");

    if (n_failed) return 1;
    return n_skipped ? EXIT_SKIPPED : 0;
}
