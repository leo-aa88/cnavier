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
    int           nx, ny;
    smtrx         DX, DY, DX2, DY2;
    solver_config cfg; // shared by the CPU and GPU solvers
    rk4_ctx       ctx; // CPU workspace
    mtrx          u, v, w;
} problem;

// A lid-driven cavity on a unit square with nx x ny nodes (bc == NULL), or
// the given wall velocities. Fields are ny rows (y) of nx values (x), as in
// main.c. The configuration is complete, operators included, before the
// workspace takes its copy, and the configuration points into *p, so p must
// stay where it is until problem_free().
static void problem_init(problem *p, int nx, int ny, int time_scheme, int poisson_type, double dt,
                         double poisson_tol, const wall_bc *bc)
{
    int order = 6;
    double dx = 1.0 / (nx - 1), dy = 1.0 / (ny - 1);
    wall_bc lid = {{0., 0., 0., 1.}, {0., 0., 0., 0.}};

    p->nx = nx; p->ny = ny;

    smtrx sd_x  = SDiff1(nx, order, dx);
    smtrx sd_y  = SDiff1(ny, order, dy);
    smtrx sd_x2 = SDiff2(nx, order, dx);
    smtrx sd_y2 = SDiff2(ny, order, dy);
    smtrx sIx   = seye(nx);
    smtrx sIy   = seye(ny);

    p->DX  = skronecker(sIy,   sd_x);
    p->DY  = skronecker(sd_y,  sIx);
    p->DX2 = skronecker(sIy,   sd_x2);
    p->DY2 = skronecker(sd_y2, sIx);

    freesm(sd_x); freesm(sd_y); freesm(sd_x2); freesm(sd_y2);
    freesm(sIx);  freesm(sIy);

    p->cfg.nx = nx; p->cfg.ny = ny; p->cfg.dx = dx; p->cfg.dy = dy;
    p->cfg.Re = 100.; p->cfg.dt = dt; p->cfg.time_scheme = time_scheme;
    p->cfg.poisson_type = poisson_type;
    p->cfg.poisson_max_it = 200000; p->cfg.poisson_tol = poisson_tol;
    p->cfg.beta = sor_beta(nx, ny, dx, dy);
    p->cfg.bc = bc ? *bc : lid;
    p->cfg.DX = &p->DX; p->cfg.DY = &p->DY; p->cfg.DX2 = &p->DX2; p->cfg.DY2 = &p->DY2;
    p->ctx = rk4_alloc(&p->cfg);

    p->u = initm(ny, nx);
    p->v = initm(ny, nx);
    p->w = initm(ny, nx);
}

static void problem_free(problem *p)
{
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
    for (i = 0; i < 12; i++) m.M[i] = (i % 5) - 2.5 * (i == 7) + 3.0 * (i == 2);
    flatten(m, flat, 3, 4);
    unflatten(flat, m2, 3, 4);
    check("flatten / unflatten round trip", rel_diff(m2.M, m.M, 12), 0.0);
    check("maxel", fabs(maxel(m) - 5.0), 0.0);
    check("minel", fabs(minel(m) - (-0.5)), 0.0);
    negcpy(m2, m);
    for (i = 0; i < 12; i++) m2.M[i] += m.M[i];
    check("negcpy", max_abs(m2.M, 12), 0.0);

    freem(&x); freem(&y_dense); freem(&d); freem(&a); freem(&b);
    freem(&k_dense); freem(&k_sparse); freem(&e_dense); freem(&e_sparse);
    freem(&m); freem(&m2);
    freesm(s); freesm(s1); freesm(s2); freesm(ks); freesm(se);
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

// Largest |a - b| over the interior nodes. The wall rows of w are left
// stale by step() (it sets the wall vorticity at the start of the next step),
// so they say nothing about the time integration.
static double interior_diff(mtrx a, mtrx b)
{
    int i, j;
    double m = 0.0;
    for (i = 1; i < a.m - 1; i++)
        for (j = 1; j < a.n - 1; j++)
        {
            double d = fabs(MAt(a, i, j) - MAt(b, i, j));
            if (isnan(d)) return NAN;
            if (d > m) m = d;
        }
    return m;
}

// Observed order of the time integration: halve dt and compare the interior
// errors against a run of the same scheme with a much smaller step. Euler is
// first order. RK4 is fourth order in its stages, but the wall vorticity is
// set once per step, from the velocities at its start, and not between the
// stages, so the coupled scheme is first order too, and at the same step its
// interior error is about Euler's (0.024 against 0.025 here). The RK4 check
// therefore only guards against it getting worse than first order; it cannot
// tell RK4 from Euler. Once the wall vorticity is updated at every stage, the
// interior error converges at fourth order and this test should require it.
static void test_temporal_order(void)
{
    int n = 17, s, k;
    double T = 0.2, order;
    char name[96];

    printf("Unit: observed order of the time integration, %dx%d grid, t = %g\n", n, n, T);
    for (s = 1; s <= 2; s++)
    {
        double err[2];
        mtrx ref = run_to(n, s, T / 4096, T);
        for (k = 0; k < 2; k++)
        {
            mtrx w = run_to(n, s, T / (64 << k), T);
            err[k] = interior_diff(w, ref);
            freem(&w);
        }
        freem(&ref);
        order = log2(err[0] / err[1]);
        if (s == 1)
        {
            snprintf(name, sizeof(name), "Euler: observed order %.2f, |order - 1|", order);
            check(name, isnan(order) ? INFINITY : fabs(order - 1.0), 0.15);
        }
        else
        {
            snprintf(name, sizeof(name), "RK4: observed order %.2f, at least first", order);
            check(name, isnan(order) ? INFINITY : 1.0 - order, 0.1);
        }
    }
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
    for (i = 0; i < nx * ny; i++) u.M[i] = v.M[i] = -1.0;
    apply_wall_bc(u, v, &bc);
    for (i = 0; i < ny; i++)
        for (j = 0; j < nx; j++)
        {
            int wall = (j == 0) ? 0 : (j == nx - 1) ? 1 : (i == 0) ? 2 : (i == ny - 1) ? 3 : -1;
            double eu = wall < 0 ? -1.0 : bc.u[wall], ev = wall < 0 ? -1.0 : bc.v[wall];
            bad += fabs(MAt(u, i, j) - eu) + fabs(MAt(v, i, j) - ev);
        }
    check("every node has its wall's velocity, interior untouched", bad, 0.0);
    freem(&u); freem(&v);
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
    if (f) { fprintf(f, "stale\n1 2 3\n"); fclose(f); }
    f = fopen("output/unittest-1-5.vtk", "w");
    if (f) { fprintf(f, "stale\n"); fclose(f); }
    f = fopen("output/other-1-5.vtk", "w");
    if (f) { fprintf(f, "keep\n"); fclose(f); }
    for (i = 0; i < nx * ny; i++) a.M[i] = 0.5 * i;
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
        if (!f || !fgets(line, sizeof(line), f)) { worst = INFINITY; break; }
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
    for (i = 0; i < 6; i++) remove(names[i]);
    rmdir("output");
    if (chdir(cwd) == 0) rmdir(tmpl);
    freem(&a); freem(&u); freem(&v);
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
                if (op == 0) { f[k] = x * (1 + y * y); e[k] = 1 + y * y; }
                if (op == 1) { f[k] = y * (1 + x * x); e[k] = 1 + x * x; }
                if (op == 2) { f[k] = x * x * (1 + y); e[k] = 2 * (1 + y); }
                if (op == 3) { f[k] = y * y * (1 + x); e[k] = 2 * (1 + x); }
            }
        spmv(op == 0 ? p.DX : op == 1 ? p.DY : op == 2 ? p.DX2 : p.DY2, f, d);
        err = rel_diff(d, e, N);
        snprintf(name, sizeof(name), "%s exact on its test polynomial", op == 0 ? "DX" : op == 1 ? "DY" : op == 2 ? "DX2" : "DY2");
        check(name, err, 1E-9);
    }

    free(f); free(d); free(e);
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
    double lambda = (2.0 * cos(PI * p / (double)(ny - 1)) - 2.0) / (dy * dy)
                  + (2.0 * cos(PI * q / (double)(nx - 1)) - 2.0) / (dx * dx);

    printf("CPU: FFT Poisson solver against an exact eigenmode, %dx%d grid\n", nx, ny);
    for (i = 1; i < ny - 1; i++)
        for (j = 1; j < nx - 1; j++)
        {
            MAt(f, i, j) = sin(PI * i * p / (double)(ny - 1))
                         * sin(PI * j * q / (double)(nx - 1));
            MAt(expected, i, j) = MAt(f, i, j) / lambda;
        }

    fft_solver *fs = fft_setup(nx, ny);
    poisson_FFT(fs, f, psi, dx, dy);
    fft_cleanup(fs);
    snprintf(name, sizeof(name), "psi vs eigenmode / eigenvalue, %dx%d", nx, ny);
    check(name, rel_diff(psi.M, expected.M, nx * ny), 1E-10);

    freem(&f); freem(&psi); freem(&expected);
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
            r = (MAt(psi, i+1, j) - 2.0 * MAt(psi, i, j) + MAt(psi, i-1, j)) / (dy * dy)
              + (MAt(psi, i, j+1) - 2.0 * MAt(psi, i, j) + MAt(psi, i, j-1)) / (dx * dx)
              - MAt(f, i, j);
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

    freem(&f); freem(&fft); freem(&sor); freem(&gs); freem(&scratch);
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
    freem(&f); freem(&psi); freem(&scratch);
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

    freem(&f); freem(&psi); freem(&scratch);
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
    out[0] = initm(n, n); out[1] = initm(n, n); out[2] = initm(n, n);
    mtrxcpy(out[0], p.w); mtrxcpy(out[1], p.u); mtrxcpy(out[2], p.v);
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
        freem(&one[k]); freem(&many[k]);
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
        freem(&u0); freem(&v0); freem(&w0);
    }
    else
        n_skipped++;
    backend_free(b);
#endif

    freem(&u); freem(&v); freem(&w);
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

    free(x); free(y_cpu); free(y_gpu);
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

    freem(&w); freem(&f); freem(&psi_cpu); freem(&psi_gpu);
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

    freem(&u); freem(&v); freem(&w);
    gpu_free(g);
    problem_free(&p);
}

// Run a GPU test, or count it as skipped when there is no device to run it on
#define GPU_TEST(call) do { if (have_gpu) (call); else n_skipped++; } while (0)

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
