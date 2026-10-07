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

static point at(const mms_case *c, double x, double y)
{
    double X[5], Y[5];
    point p;

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

mms_errors mms_run(int nx, int ny, double Lx, double Ly, double Re, int order, int time_scheme,
                   int poisson_type, double dt, double t0, double T)
{
    mms_case c = {Lx, Ly, Re, Lx / (nx - 1), Ly / (ny - 1)};
    wall_bc walls = {{0., 0., 0., 0.}, {0., 0., 0., 0.}};
    int t, steps;
    mms_errors e;

    smtrx d1x = SDiff1(nx, order, c.dx), d1y = SDiff1(ny, order, c.dy);
    smtrx d2x = SDiff2(nx, order, c.dx), d2y = SDiff2(ny, order, c.dy);
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
    apply_wall_bc(u, v, &walls);
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
