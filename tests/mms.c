#include <math.h>
#include "mms.h"
#include "finitediff.h"
#include "fluiddyn.h"
#include "poisson.h"

#define MMS_PI 3.14159265358979323846

// Time factor of the solution, and its derivative
static double g(double t)
{
    return 0.2 + 0.1 * sin(2.0 * MMS_PI * t);
}
static double dg(double t)
{
    return 0.2 * MMS_PI * cos(2.0 * MMS_PI * t);
}

// S(s) = sin^2(a s) and its derivatives up to the fourth
static void shape(double a, double s, double d[5])
{
    double sn = sin(a * s), s2 = sin(2.0 * a * s), c2 = cos(2.0 * a * s);
    d[0] = sn * sn;
    d[1] = a * s2;
    d[2] = 2.0 * a * a * c2;
    d[3] = -4.0 * a * a * a * s2;
    d[4] = -8.0 * a * a * a * a * c2;
}

// Steady parts at one point: psi_s = X Y, u_s = X Y', v_s = -X' Y,
// w_s = -(X'' Y + X Y''), its gradient and its Laplacian
typedef struct
{
    double psi, u, v, w, wx, wy, lapw;
} point;

// The periodic solution: modes (a, b, amplitude, phase)
static const double modes[3][4] = {
    {1.0, 1.0, 1.0, 0.3},
    {2.0, -1.0, 0.5, 1.1},
    {1.0, 3.0, 0.25, 2.0},
};

static point at_periodic(const mms_case *c, double x, double y)
{
    int m;
    point p = {0, 0, 0, 0, 0, 0, 0};

    for (m = 0; m < 3; m++)
    {
        double kx = 2.0 * MMS_PI * modes[m][0] / c->Lx, ky = 2.0 * MMS_PI * modes[m][1] / c->Ly;
        double A = modes[m][2] / (2.0 * MMS_PI), K2 = kx * kx + ky * ky;
        double th = kx * x + ky * y + modes[m][3], cs = cos(th), sn = sin(th);
        p.psi += A * cs;
        p.u += -A * ky * sn;
        p.v += A * kx * sn;
        p.w += A * K2 * cs;
        p.wx += -A * K2 * kx * sn;
        p.wy += -A * K2 * ky * sn;
        p.lapw += -A * K2 * K2 * cs;
    }
    return p;
}

static point at(const mms_case *c, double x, double y)
{
    double X[5], Y[5];
    point p;

    if (c->periodic) return at_periodic(c, x, y);

    shape(MMS_PI / c->Lx, x, X);
    shape(MMS_PI / c->Ly, y, Y);
    p.psi = X[0] * Y[0];
    p.u = X[0] * Y[1];
    p.v = -X[1] * Y[0];
    p.w = -(X[2] * Y[0] + X[0] * Y[2]);
    p.wx = -(X[3] * Y[0] + X[1] * Y[2]);
    p.wy = -(X[2] * Y[1] + X[0] * Y[3]);
    p.lapw = -(X[4] * Y[0] + 2.0 * X[2] * Y[2] + X[0] * Y[4]);
    return p;
}

void mms_exact(const mms_case *c, double t, mtrx *w, mtrx *u, mtrx *v, mtrx *psi)
{
    int i, j;
    double gt = g(t);
    mtrx *f = w ? w : u ? u
                  : v   ? v
                        : psi;

    if (!f) return;
    for (i = 0; i < f->m; i++)
        for (j = 0; j < f->n; j++)
        {
            point p = at(c, j * c->dx, i * c->dy);
            if (w) MAt(*w, i, j) = gt * p.w;
            if (u) MAt(*u, i, j) = gt * p.u;
            if (v) MAt(*v, i, j) = gt * p.v;
            if (psi) MAt(*psi, i, j) = gt * p.psi;
        }
}

void mms_source(double t, mtrx f, void *data)
{
    const mms_case *c = (const mms_case *)data;
    int i, j;
    double gt = g(t), dgt = dg(t);

    // w = g w_s, u = g u_s, v = g v_s:
    // f = g' w_s + g^2 (u_s w_s,x + v_s w_s,y) - (g / Re) lap(w_s)
    for (i = 0; i < f.m; i++)
        for (j = 0; j < f.n; j++)
        {
            point p = at(c, j * c->dx, i * c->dy);
            MAt(f, i, j) = dgt * p.w + gt * gt * (p.u * p.wx + p.v * p.wy) - gt * p.lapw / c->Re;
        }
}

// Largest and rms of |a - b| over all nodes (part 0), the wall nodes (1) or
// the interior nodes (2)
static mms_norm norm_of(mtrx a, mtrx b, int part)
{
    int i, j, count = 0;
    mms_norm r = {0.0, 0.0};

    for (i = 0; i < a.m; i++)
        for (j = 0; j < a.n; j++)
        {
            int wall = i == 0 || j == 0 || i == a.m - 1 || j == a.n - 1;
            double d = fabs(MAt(a, i, j) - MAt(b, i, j));
            if ((part == 1 && !wall) || (part == 2 && wall)) continue;
            if (isnan(d)) d = INFINITY;
            if (d > r.max) r.max = d;
            r.rms += d * d;
            count++;
        }
    r.rms = sqrt(r.rms / count);
    return r;
}

static mms_errors run(int nx, int ny, double Lx, double Ly, double Re, int order, int time_scheme,
                      int poisson_type, double dt, double t0, double T, int periodic)
{
    mms_case c = {Lx, Ly, Re, Lx / (periodic ? nx : nx - 1), Ly / (periodic ? ny : ny - 1), periodic};
    wall_bc walls = {{0., 0., 0., 0.}, {0., 0., 0., 0.}};
    int t, steps;
    mms_errors e;

    smtrx d1x = periodic ? SDiff1_periodic(nx, order, c.dx) : SDiff1(nx, order, c.dx);
    smtrx d1y = periodic ? SDiff1_periodic(ny, order, c.dy) : SDiff1(ny, order, c.dy);
    smtrx d2x = periodic ? SDiff2_periodic(nx, order, c.dx) : SDiff2(nx, order, c.dx);
    smtrx d2y = periodic ? SDiff2_periodic(ny, order, c.dy) : SDiff2(ny, order, c.dy);
    dt = fmin(dt, 0.5 * max_stable_dt(&d2x, &d2y, Re, time_scheme));
    steps = (int)ceil(T / dt - 1E-9);
    smtrx Ix = seye(nx), Iy = seye(ny);
    smtrx DX = skronecker(Iy, d1x), DY = skronecker(d1y, Ix);
    smtrx DX2 = skronecker(Iy, d2x), DY2 = skronecker(d2y, Ix);
    freesm(d1x);
    freesm(d1y);
    freesm(d2x);
    freesm(d2y);
    freesm(Ix);
    freesm(Iy);

    solver_config cfg;
    cfg.nx = nx;
    cfg.ny = ny;
    cfg.dx = c.dx;
    cfg.dy = c.dy;
    cfg.Re = Re;
    cfg.dt = T / steps;
    cfg.t0 = t0;
    cfg.time_scheme = time_scheme;
    cfg.poisson_type = poisson_type;
    cfg.poisson_max_it = 1000000;
    cfg.poisson_tol = 1E-13;
    cfg.beta = sor_beta(nx, ny, c.dx, c.dy);
    cfg.periodic = periodic;
    cfg.bc = walls;
    cfg.DX = &DX;
    cfg.DY = &DY;
    cfg.DX2 = &DX2;
    cfg.DY2 = &DY2;
    cfg.vorticity_source = mms_source;
    cfg.source_data = &c;
    rk4_ctx ctx = rk4_alloc(&cfg);

    mtrx w = initm(ny, nx), u = initm(ny, nx), v = initm(ny, nx);
    mtrx we = initm(ny, nx), ue = initm(ny, nx), ve = initm(ny, nx), psie = initm(ny, nx);

    mms_exact(&c, t0, &w, &u, &v, NULL);
    for (t = 0; t < steps; t++)
        step(w, u, v, &ctx);

    mms_exact(&c, t0 + T, &we, &ue, &ve, &psie);
    if (!periodic) apply_wall_bc(u, v, &walls);
    e.psi = norm_of(ctx.psi, psie, 0);
    e.u = norm_of(u, ue, 0);
    e.v = norm_of(v, ve, 0);
    e.w = norm_of(w, we, 0);
    e.w_wall = norm_of(w, we, 1);
    e.w_interior = norm_of(w, we, 2);
    e.steps = steps;

    freem(&w);
    freem(&u);
    freem(&v);
    freem(&we);
    freem(&ue);
    freem(&ve);
    freem(&psie);
    rk4_free(&ctx);
    freesm(DX);
    freesm(DY);
    freesm(DX2);
    freesm(DY2);
    return e;
}

mms_errors mms_run(int nx, int ny, double Lx, double Ly, double Re, int order, int time_scheme,
                   int poisson_type, double dt, double t0, double T)
{
    return run(nx, ny, Lx, Ly, Re, order, time_scheme, poisson_type, dt, t0, T, 0);
}

mms_errors mms_run_periodic(int nx, int ny, double Lx, double Ly, double Re, int order, int time_scheme,
                            double dt, double t0, double T)
{
    return run(nx, ny, Lx, Ly, Re, order, time_scheme, 3, dt, t0, T, 1);
}

// ---------------------------------------------------------------------------
// Ablation: a copy of step() (RK4, FFT solver) with parts replaced by the
// exact solution
// ---------------------------------------------------------------------------

typedef struct
{
    mms_case c;
    int flags;
    smtrx DX, DY, DX2, DY2;
    fft_solver *fft;
    mtrx psi, rhs, u, v, wx, wy, wxx, wyy, dvx, duy, src, exact, k1, k2, k3, k4, w_tmp;
} ablation;

// u, v of the stage: from the exact solution, or from psi (exact, or solved
// from the interior of w) with the wall velocities (zero) imposed
static void ab_velocity(ablation *a, mtrx w, double t)
{
    static const wall_bc walls = {{0., 0., 0., 0.}, {0., 0., 0., 0.}};

    if (a->flags & MMS_EXACT_VELOCITY)
    {
        mms_exact(&a->c, t, NULL, &a->u, &a->v, NULL);
        return;
    }
    if (a->flags & MMS_EXACT_PSI)
        mms_exact(&a->c, t, NULL, NULL, NULL, &a->psi);
    else
    {
        negcpy(a->rhs, w);
        poisson_FFT(a->fft, a->rhs, a->psi, a->c.dx, a->c.dy);
    }
    spmv(a->DY, a->psi.M, a->u.M);
    spmv(a->DX, a->psi.M, a->v.M);
    negcpy(a->v, a->v);
    apply_wall_bc(a->u, a->v, &walls);
}

// The boundary entries of w: exact, or D_x v - D_y u as set_wall_vorticity()
// computes them
static void ab_wall_vorticity(ablation *a, mtrx w, double t)
{
    int i, j, ny = w.m, nx = w.n;

    if (a->flags & MMS_EXACT_WALL_W)
        mms_exact(&a->c, t, &a->exact, NULL, NULL, NULL);
    else
    {
        spmv(a->DX, a->v.M, a->dvx.M);
        spmv(a->DY, a->u.M, a->duy.M);
    }
    for (i = 0; i < ny; i++)
        for (j = 0; j < nx; j++)
            if (i == 0 || j == 0 || i == ny - 1 || j == nx - 1)
                MAt(w, i, j) = (a->flags & MMS_EXACT_WALL_W) ? MAt(a->exact, i, j)
                                                             : MAt(a->dvx, i, j) - MAt(a->duy, i, j);
}

// dwdt(): velocity, wall vorticity, derivatives, right-hand side
static void ab_rhs(ablation *a, mtrx w, double t, mtrx out)
{
    int k, n = w.m * w.n;

    ab_velocity(a, w, t);
    ab_wall_vorticity(a, w, t);
    spmv(a->DX, w.M, a->wx.M);
    spmv(a->DY, w.M, a->wy.M);
    spmv(a->DX2, w.M, a->wxx.M);
    spmv(a->DY2, w.M, a->wyy.M);
    for (k = 0; k < n; k++)
        out.M[k] = -a->u.M[k] * a->wx.M[k] - a->v.M[k] * a->wy.M[k] + (1.0 / a->c.Re) * (a->wxx.M[k] + a->wyy.M[k]);
    mms_source(t, a->src, &a->c);
    for (k = 0; k < n; k++)
        out.M[k] += a->src.M[k];
}

// step() with RK4: four stages, the final velocity and wall vorticity
static void ab_step(ablation *a, mtrx w, double t, double dt)
{
    int k, n = w.m * w.n;

    ab_rhs(a, w, t, a->k1);
    for (k = 0; k < n; k++)
        a->w_tmp.M[k] = w.M[k] + 0.5 * dt * a->k1.M[k];
    ab_rhs(a, a->w_tmp, t + 0.5 * dt, a->k2);
    for (k = 0; k < n; k++)
        a->w_tmp.M[k] = w.M[k] + 0.5 * dt * a->k2.M[k];
    ab_rhs(a, a->w_tmp, t + 0.5 * dt, a->k3);
    for (k = 0; k < n; k++)
        a->w_tmp.M[k] = w.M[k] + dt * a->k3.M[k];
    ab_rhs(a, a->w_tmp, t + dt, a->k4);
    for (k = 0; k < n; k++)
        w.M[k] += (dt / 6.0) * (a->k1.M[k] + 2.0 * a->k2.M[k] + 2.0 * a->k3.M[k] + a->k4.M[k]);
    ab_velocity(a, w, t + dt);
    ab_wall_vorticity(a, w, t + dt);
}

mms_errors mms_run_ablated(int nx, int ny, double Lx, double Ly, double Re, int order, int flags,
                           double dt, double t0, double T)
{
    ablation a;
    mtrx *fields[] = {&a.psi, &a.rhs, &a.u, &a.v, &a.wx, &a.wy, &a.wxx, &a.wyy, &a.dvx, &a.duy,
                      &a.src, &a.exact, &a.k1, &a.k2, &a.k3, &a.k4, &a.w_tmp};
    int f, t, steps, nf = (int)(sizeof(fields) / sizeof(fields[0]));
    mms_errors e;

    a.c = (mms_case){Lx, Ly, Re, Lx / (nx - 1), Ly / (ny - 1), 0};
    a.flags = flags;
    smtrx d1x = SDiff1(nx, order, a.c.dx), d1y = SDiff1(ny, order, a.c.dy);
    smtrx d2x = SDiff2(nx, order, a.c.dx), d2y = SDiff2(ny, order, a.c.dy);
    dt = fmin(dt, 0.5 * max_stable_dt(&d2x, &d2y, Re, 2));
    steps = (int)ceil(T / dt - 1E-9);
    dt = T / steps;
    smtrx Ix = seye(nx), Iy = seye(ny);
    a.DX = skronecker(Iy, d1x);
    a.DY = skronecker(d1y, Ix);
    a.DX2 = skronecker(Iy, d2x);
    a.DY2 = skronecker(d2y, Ix);
    freesm(d1x);
    freesm(d1y);
    freesm(d2x);
    freesm(d2y);
    freesm(Ix);
    freesm(Iy);
    a.fft = fft_setup(nx, ny);
    for (f = 0; f < nf; f++)
        *fields[f] = initm(ny, nx);

    mtrx w = initm(ny, nx), we = initm(ny, nx), ue = initm(ny, nx), ve = initm(ny, nx), psie = initm(ny, nx);
    mms_exact(&a.c, t0, &w, NULL, NULL, NULL);
    for (t = 0; t < steps; t++)
        ab_step(&a, w, t0 + (double)t * dt, dt);

    mms_exact(&a.c, t0 + T, &we, &ue, &ve, &psie);
    e.psi = norm_of(a.psi, psie, 0);
    e.u = norm_of(a.u, ue, 0);
    e.v = norm_of(a.v, ve, 0);
    e.w = norm_of(w, we, 0);
    e.w_wall = norm_of(w, we, 1);
    e.w_interior = norm_of(w, we, 2);
    e.steps = steps;

    freem(&w);
    freem(&we);
    freem(&ue);
    freem(&ve);
    freem(&psie);
    for (f = 0; f < nf; f++)
        freem(fields[f]);
    fft_cleanup(a.fft);
    freesm(a.DX);
    freesm(a.DY);
    freesm(a.DX2);
    freesm(a.DY2);
    return e;
}
