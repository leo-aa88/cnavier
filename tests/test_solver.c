// Test suite: make test            (CPU solver)
//             make CUDA=1 test     (CPU solver + GPU backend checked against it)
//
// Exit status: 0 = every check ran and passed, 1 = a check failed,
// 77 = nothing failed but the GPU tests could not run (no usable device).

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "linearalg.h"
#include "finitediff.h"
#include "poisson.h"
#include "fluiddyn.h"
#ifdef USE_CUDA
#include "cudasolver.h"
#endif

#define EXIT_SKIPPED 77 // conventional "skipped" status: neither a pass nor a failure

static int n_checks = 0;
static int n_failed = 0;
static int n_skipped = 0; // GPU tests that could not run

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
    int     nx, ny;
    double  dt;
    int     time_scheme;
    smtrx   DX, DY, DX2, DY2;
    rk4_ctx ctx;
    wall_bc bc;
    mtrx    u, v, w;
} problem;

static problem problem_alloc(int n, int time_scheme, int poisson_type, double dt, double poisson_tol)
{
    problem p;
    int order = 6;
    double dx = 1.0 / n, dy = 1.0 / n;
    double rho = 0.5 * (cos(PI / n) + cos(PI / n));
    wall_bc lid = {{0., 0., 0., 1.}, {0., 0., 0., 0.}};

    p.nx = n; p.ny = n; p.dt = dt; p.time_scheme = time_scheme;
    p.bc = lid;

    smtrx sd_x  = SDiff1(n, order, dx);
    smtrx sd_y  = SDiff1(n, order, dy);
    smtrx sd_x2 = SDiff2(n, order, dx);
    smtrx sd_y2 = SDiff2(n, order, dy);
    smtrx sIx   = seye(n);
    smtrx sIy   = seye(n);

    p.DX  = skronecker(sIy,   sd_x);
    p.DY  = skronecker(sd_y,  sIx);
    p.DX2 = skronecker(sIy,   sd_x2);
    p.DY2 = skronecker(sd_y2, sIx);

    freesm(sd_x); freesm(sd_y); freesm(sd_x2); freesm(sd_y2);
    freesm(sIx);  freesm(sIy);

    p.ctx = rk4_alloc(n, n);
    p.ctx.Re = 100.; p.ctx.dx = dx; p.ctx.dy = dy;
    p.ctx.poisson_type = poisson_type;
    p.ctx.poisson_max_it = 200000; p.ctx.poisson_tol = poisson_tol;
    p.ctx.beta = 2.0 / (1.0 + sqrt(1.0 - rho * rho));

    p.u = initm(n, n);
    p.v = initm(n, n);
    p.w = initm(n, n);

    if (poisson_type == 3) fft_setup(n, n);
    return p;
}

// The operators are referenced through pointers, so bind them once the
// problem sits at its final address.
static void problem_bind(problem *p)
{
    p->ctx.DX = &p->DX; p->ctx.DY = &p->DY; p->ctx.DX2 = &p->DX2; p->ctx.DY2 = &p->DY2;
}

static void problem_free(problem *p)
{
    if (p->ctx.poisson_type == 3) fft_cleanup();
    freem(&p->u); freem(&p->v); freem(&p->w);
    freesm(p->DX); freesm(p->DY); freesm(p->DX2); freesm(p->DY2);
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
    *cmax = -__DBL_MAX__; *cmin = __DBL_MAX__;
    for (k = 0; k < n; k++)
    {
        double c = dudx[k] + dvdy[k];
        if (c > *cmax) *cmax = c;
        if (c < *cmin) *cmin = c;
    }
    free(dudx); free(dvdy);
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
                    if (isnan(sum)) worst = NAN;
                    else if (!isnan(worst) && fabs(sum - expect) > worst) worst = fabs(sum - expect);
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
        int c0[] = {0, 1};         double v0[] = {-h, h};
        int c1[] = {0, 2};         double v1[] = {-0.5 * h, 0.5 * h};
        int c2[] = {0, 1, 3, 4};   double v2[] = {1.0 / 12.0 * h, -2.0 / 3.0 * h, 2.0 / 3.0 * h, -1.0 / 12.0 * h};
        int c5[] = {2, 3, 4, 6, 7, 8};
        double v5[] = {-1.0 / 60.0 * h, 3.0 / 20.0 * h, -3.0 / 4.0 * h, 3.0 / 4.0 * h, -3.0 / 20.0 * h, 1.0 / 60.0 * h};
        int c8[] = {7, 9};         double v8[] = {-0.5 * h, 0.5 * h};
        int c9[] = {8, 9};         double v9[] = {-h, h};
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
        int c0[] = {0, 1, 2, 3};   double v0[] = {2 * h2, -5 * h2, 4 * h2, -1 * h2};
        int c1[] = {0, 1, 2};      double v1[] = {h2, -2 * h2, h2};
        int c5[] = {2, 3, 4, 5, 6, 7, 8};
        double v5[] = {1.0 / 90.0 * h2, -3.0 / 20.0 * h2, 3.0 / 2.0 * h2, -49.0 / 18.0 * h2,
                       3.0 / 2.0 * h2, -3.0 / 20.0 * h2, 1.0 / 90.0 * h2};
        int c9[] = {6, 7, 8, 9};   double v9[] = {-1 * h2, 4 * h2, -5 * h2, 2 * h2};
        worst = fmax(worst, row_diff(d2, 0, 4, c0, v0));
        worst = fmax(worst, row_diff(d2, 1, 3, c1, v1));
        worst = fmax(worst, row_diff(d2, 5, 7, c5, v5));
        worst = fmax(worst, row_diff(d2, 9, 4, c9, v9));
    }
    check("second derivative, order 6: rows 0, 1, 5, 9", worst, 1E-12);

    freesm(d1); freesm(d2);
}

// ---------------------------------------------------------------------------
// CPU tests
// ---------------------------------------------------------------------------

// A sine mode is an eigenvector of the discrete Laplacian, so the FFT solver
// must return it divided by its eigenvalue.
static void test_cpu_poisson_fft(void)
{
    int i, j, n = 32, p = 3, q = 5;
    double dx = 1.0 / n;
    mtrx f = initm(n, n), psi = initm(n, n), expected = initm(n, n);
    double lambda = (2.0 * cos(PI * p / (double)(n + 1)) - 2.0) / (dx * dx)
                  + (2.0 * cos(PI * q / (double)(n + 1)) - 2.0) / (dx * dx);

    printf("CPU: FFT Poisson solver against an exact eigenmode\n");
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
        {
            MAt(f, i, j) = sin(PI * (i + 1) * p / (double)(n + 1))
                         * sin(PI * (j + 1) * q / (double)(n + 1));
            MAt(expected, i, j) = MAt(f, i, j) / lambda;
        }

    fft_setup(n, n);
    poisson_FFT(f, psi, dx, dx);
    fft_cleanup();
    check("psi vs eigenmode / eigenvalue", rel_diff(psi.M, expected.M, n * n), 1E-10);

    freem(&f); freem(&psi); freem(&expected);
}

// Largest residual of the 5-point Laplacian at the interior points
static double poisson_residual(mtrx f, mtrx psi, double dx, double dy)
{
    int i, j;
    double r, rmax = 0.0;

    for (i = 1; i < f.m - 1; i++)
        for (j = 1; j < f.n - 1; j++)
        {
            r = (MAt(psi, i+1, j) - 2.0 * MAt(psi, i, j) + MAt(psi, i-1, j)) / (dx * dx)
              + (MAt(psi, i, j+1) - 2.0 * MAt(psi, i, j) + MAt(psi, i, j-1)) / (dy * dy)
              - MAt(f, i, j);
            if (isnan(r)) return NAN;
            if (fabs(r) > rmax) rmax = fabs(r);
        }
    return rmax;
}

static void test_cpu_poisson_iterative(void)
{
    int n = 24;
    double dx = 1.0 / n;
    double rho = cos(PI / n);
    double beta = 2.0 / (1.0 + sqrt(1.0 - rho * rho));
    mtrx f = initm(n, n), psi = initm(n, n), scratch = initm(n, n);

    printf("CPU: Gauss-Seidel and SOR residuals\n");
    fill_pseudo_random(f.M, n * n, 7u);

    poisson(f, psi, scratch, dx, dx, 200000, 1E-9);
    check("Gauss-Seidel residual", poisson_residual(f, psi, dx, dx), 1E-4);
    poisson_SOR(f, psi, scratch, dx, dx, 200000, 1E-9, beta);
    check("SOR residual", poisson_residual(f, psi, dx, dx), 1E-4);

    freem(&f); freem(&psi); freem(&scratch);
}

// Short lid-driven cavity run: the fields must stay finite, the lid must
// have spun up the flow, and the velocity field must be divergence-free.
static void test_cpu_step(int time_scheme, const char *label)
{
    int t, n = 32;
    double cmax, cmin;
    char name[96];
    problem p = problem_alloc(n, time_scheme, 3, 0.002, 1E-3);
    problem_bind(&p);

    printf("CPU: 50 steps of the lid-driven cavity, %s + FFT\n", label);
    for (t = 0; t < 50; t++)
        step(p.w, p.u, p.v, p.dt, p.time_scheme, &p.bc, &p.ctx);
    continuity_range(&p, &cmax, &cmin);

    snprintf(name, sizeof(name), "%s: fields finite (max |w|)", label);
    check(name, max_abs(p.w.M, n * n), 1E6);
    snprintf(name, sizeof(name), "%s: lid drives the flow (-max |u|)", label);
    check(name, -max_abs(p.u.M, n * n), -1E-3);
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
        problem p = problem_alloc(n, time_scheme, 3, 0.0, 1E-3);
        problem_bind(&p);
        smtrx d2 = SDiff2(n, 6, p.ctx.dx);
        double limit = max_stable_dt(&d2, &d2, p.ctx.Re, time_scheme);
        freesm(d2);
        double wmax = 0.0;

        p.dt = frac[f] * limit;
        for (t = 0; t < 3000 && !(wmax > 1E6); t++)
        {
            step(p.w, p.u, p.v, p.dt, p.time_scheme, &p.bc, &p.ctx);
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
// GPU tests — every check compares the device result with the CPU one
// ---------------------------------------------------------------------------

#ifdef USE_CUDA

static gpu_solver *gpu_for(problem *p)
{
    gpu_solver *g = gpu_init(&p->ctx, p->dt, p->time_scheme, &p->bc);
    if (!g)
    {
        printf("** Error: CUDA device disappeared during the tests **\n");
        exit(1);
    }
    return g;
}

static void test_gpu_spmv(int n)
{
    int op, N = n * n;
    char name[96];
    const char *op_name[] = {"DX", "DY", "DX2", "DY2"};
    problem p = problem_alloc(n, 2, 3, 0.002, 1E-3);
    problem_bind(&p);
    gpu_solver *g = gpu_for(&p);
    smtrx *ops[] = {&p.DX, &p.DY, &p.DX2, &p.DY2};
    double *x = (double *)malloc(N * sizeof(double));
    double *y_cpu = (double *)malloc(N * sizeof(double));
    double *y_gpu = (double *)malloc(N * sizeof(double));

    printf("GPU: CSR SpMV, %dx%d grid\n", n, n);
    fill_pseudo_random(x, N, 11u);
    for (op = 0; op < 4; op++)
    {
        spmv(*ops[op], x, y_cpu);
        gpu_spmv(g, op, x, y_gpu);
        snprintf(name, sizeof(name), "%s*x vs CPU", op_name[op]);
        check(name, rel_diff(y_gpu, y_cpu, N), 1E-12);
    }

    free(x); free(y_cpu); free(y_gpu);
    gpu_free(g);
    problem_free(&p);
}

static void test_gpu_poisson(int n, int poisson_type, const char *label, double limit)
{
    int N = n * n;
    char name[96];
    problem p = problem_alloc(n, 2, poisson_type, 0.002, 1E-10);
    problem_bind(&p);
    gpu_solver *g = gpu_for(&p);
    mtrx w = initm(n, n), f = initm(n, n), psi_cpu = initm(n, n), psi_gpu = initm(n, n);

    printf("GPU: %s Poisson solver, %dx%d grid\n", label, n, n);
    fill_pseudo_random(w.M, N, 23u);

    // CPU solvers take the right-hand side f = -w
    mtrxcpy(f, w);
    invsig(f);
    if (poisson_type == 1)
        poisson(f, psi_cpu, p.ctx.psi_scratch, p.ctx.dx, p.ctx.dy,
                p.ctx.poisson_max_it, p.ctx.poisson_tol);
    else if (poisson_type == 2)
        poisson_SOR(f, psi_cpu, p.ctx.psi_scratch, p.ctx.dx, p.ctx.dy,
                    p.ctx.poisson_max_it, p.ctx.poisson_tol, p.ctx.beta);
    else
        poisson_FFT(f, psi_cpu, p.ctx.dx, p.ctx.dy);

    gpu_poisson(g, w.M, psi_gpu.M);
    snprintf(name, sizeof(name), "%s: psi vs CPU", label);
    check(name, rel_diff(psi_gpu.M, psi_cpu.M, N), limit);
    if (poisson_type != 3)
    {
        snprintf(name, sizeof(name), "%s: residual", label);
        check(name, poisson_residual(f, psi_gpu, p.ctx.dx, p.ctx.dy), 1E-4);
    }

    freem(&w); freem(&f); freem(&psi_cpu); freem(&psi_gpu);
    gpu_free(g);
    problem_free(&p);
}

// Run the same case on both backends and compare the fields afterwards.
// bc == NULL keeps the lid-driven cavity.
static void test_gpu_step(int n, int steps, double dt, int time_scheme, int poisson_type,
                          double poisson_tol, const wall_bc *bc, const char *label, double limit)
{
    int t, N = n * n;
    char name[96];
    problem p = problem_alloc(n, time_scheme, poisson_type, dt, poisson_tol);
    problem_bind(&p);
    if (bc) p.bc = *bc;
    gpu_solver *g = gpu_for(&p);
    mtrx u = initm(n, n), v = initm(n, n), w = initm(n, n);

    printf("GPU: %d steps, %s, %dx%d grid\n", steps, label, n, n);
    gpu_set_fields(g, &p.u, &p.v, &p.w);
    for (t = 0; t < steps; t++)
    {
        step(p.w, p.u, p.v, p.dt, p.time_scheme, &p.bc, &p.ctx);
        gpu_step(g);
    }
    gpu_get_fields(g, &u, &v, &w);

    snprintf(name, sizeof(name), "%s: w vs CPU", label);
    check(name, rel_diff(w.M, p.w.M, N), limit);
    snprintf(name, sizeof(name), "%s: u vs CPU", label);
    check(name, rel_diff(u.M, p.u.M, N), limit);
    snprintf(name, sizeof(name), "%s: v vs CPU", label);
    check(name, rel_diff(v.M, p.v.M, N), limit);

    freem(&u); freem(&v); freem(&w);
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
    problem p = problem_alloc(n, 2, 3, 0.002, 1E-3);
    problem_bind(&p);
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

    freem(&u); freem(&v); freem(&w);
    gpu_free(g);
    problem_free(&p);
}

// Run a GPU test, or count it as skipped when there is no device to run it on
#define GPU_TEST(call) do { if (have_gpu) call; else n_skipped++; } while (0)

static void run_gpu_tests(void)
{
    // A different velocity on every wall, so a wall mix-up cannot go unnoticed
    wall_bc four_walls = {{0.3, -0.2, 0.5, 1.0}, {0.1, -0.4, 0.2, -0.3}};

    // Probe for a device first
    problem probe = problem_alloc(16, 2, 3, 0.002, 1E-3);
    problem_bind(&probe);
    gpu_solver *g = gpu_init(&probe.ctx, probe.dt, probe.time_scheme, &probe.bc);
    int have_gpu = (g != NULL);
    problem_free(&probe);
    if (have_gpu)
        gpu_free(g);
    else
        printf("GPU: no usable CUDA device - GPU tests SKIPPED\n");

    // 300x300 is large enough for the reductions to span several passes
    GPU_TEST(test_gpu_fields(16));
    GPU_TEST(test_gpu_fields(300));
    GPU_TEST(test_gpu_spmv(32));
    GPU_TEST(test_gpu_spmv(45));
    GPU_TEST(test_gpu_spmv(300));
    GPU_TEST(test_gpu_poisson(32, 3, "FFT", 1E-12));
    GPU_TEST(test_gpu_poisson(45, 3, "FFT", 1E-12));
    GPU_TEST(test_gpu_poisson(24, 2, "SOR", 1E-7));
    GPU_TEST(test_gpu_poisson(24, 1, "Gauss-Seidel", 1E-7));
    GPU_TEST(test_gpu_step(32, 50, 0.002, 2, 3, 1E-10, NULL, "RK4 + FFT", 1E-11));
    GPU_TEST(test_gpu_step(32, 50, 0.002, 1, 3, 1E-10, NULL, "Euler + FFT", 1E-11));
    GPU_TEST(test_gpu_step(45, 50, 0.002, 2, 3, 1E-10, NULL, "RK4 + FFT", 1E-11));
    GPU_TEST(test_gpu_step(128, 10, 0.0005, 2, 3, 1E-10, NULL, "RK4 + FFT", 1E-11));
    GPU_TEST(test_gpu_step(32, 20, 0.002, 2, 3, 1E-10, &four_walls,
                           "RK4 + FFT, four moving walls", 1E-11));

    // The iterative solvers are run to a tight tolerance here. The default
    // CPU build sweeps lexicographically and the GPU red-black, so the two
    // only agree once the iteration has converged.
    GPU_TEST(test_gpu_step(24, 3, 0.002, 2, 2, 1E-10, NULL, "RK4 + SOR", 1E-7));
    GPU_TEST(test_gpu_step(24, 3, 0.002, 1, 1, 1E-10, NULL, "Euler + Gauss-Seidel", 1E-7));
#ifdef _OPENMP
    // With OPENMP=1 the CPU sweeps red-black too, so the backends must also
    // agree at the shipped tolerance, where the iteration stops far from
    // convergence.
    GPU_TEST(test_gpu_step(64, 10, 0.005, 2, 2, 1E-3, NULL,
                           "RK4 + SOR, shipped tolerance", 1E-9));
    GPU_TEST(test_gpu_step(24, 10, 0.002, 1, 1, 1E-3, NULL,
                           "Euler + Gauss-Seidel, shipped tolerance", 1E-9));
#endif
}

#endif // USE_CUDA

int main(void)
{
    test_finitediff_exactness();
    test_finitediff_rows();
    test_cpu_poisson_fft();
    test_cpu_poisson_iterative();
    test_cpu_step(2, "RK4");
    test_cpu_step(1, "Euler");
    test_cpu_stability_limit(1, "Euler");
    test_cpu_stability_limit(2, "RK4");
    test_cpu_time_step_limits();
#ifdef USE_CUDA
    run_gpu_tests();
#endif

    printf("\n%d checks, %d failed", n_checks, n_failed);
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
