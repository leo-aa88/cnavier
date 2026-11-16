// Time integrators in decaying 2-D turbulence (#48): the largest stable time
// step of RK4, AB3 and RK46, and their error against the cost, as the issue
// asks.
//
//   make integrators                       256^2, Re = 2e4, on the GPU if built with CUDA=1
//   ./integrators_study --n 128 --cpu      another grid, or the CPU
//
// The case is --case decaying (random phases, spectrum peak at |k|/2pi = 10,
// E = 1/2, order 6, skew-symmetric nonlinear term), as in the resolution
// study. Stability: the largest dt at which max |w| stays below twice its
// initial value up to t = 0.2, by bisection. Accuracy: each scheme at
// fractions of its own limit up to t = 0.1, before the runs decorrelate,
// against RK4 with a step 64 times below its limit: the relative L2 error of
// w and the largest relative error of the enstrophy spectrum Z(k) over the
// shells up to two thirds of the Nyquist wavenumber, with the wall-clock time.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "linearalg.h"
#include "finitediff.h"
#include "fluiddyn.h"
#include "forcing.h"
#include "backend.h"
#include "diagnostics.h"

typedef struct
{
    int n, gpu;
    double Re;
    smtrx d1, d2, DX, DY, DX2, DY2;
} study;

static void study_init(study *s, int n, double Re, int gpu)
{
    double h = 1.0 / n;
    smtrx I = seye(n);
    s->n = n;
    s->Re = Re;
    s->gpu = gpu;
    s->d1 = SDiff1_periodic(n, 6, h);
    s->d2 = SDiff2_periodic(n, 6, h);
    s->DX = skronecker(I, s->d1);
    s->DY = skronecker(s->d1, I);
    s->DX2 = skronecker(I, s->d2);
    s->DY2 = skronecker(s->d2, I);
    freesm(I);
}

static void study_free(study *s)
{
    freesm(s->d1);
    freesm(s->d2);
    freesm(s->DX);
    freesm(s->DY);
    freesm(s->DX2);
    freesm(s->DY2);
}

static solver_config config(const study *s, int scheme, double dt)
{
    solver_config c = {0};
    c.nx = c.ny = s->n;
    c.dx = c.dy = 1.0 / s->n;
    c.Re = s->Re;
    c.dt = dt;
    c.time_scheme = scheme;
    c.poisson_type = 3;
    c.poisson_order = 2;
    c.periodic = 1;
    c.advection = 1;
    c.DX = &s->DX;
    c.DY = &s->DY;
    c.DX2 = &s->DX2;
    c.DY2 = &s->DY2;
    c.D1x = c.D1y = &s->d1;
    c.D2x = c.D2y = &s->d2;
    return c;
}

static double max_abs(const double *x, int n)
{
    // A NaN anywhere makes the result NaN (fmax() alone would ignore it, and
    // a field that blew up would look bounded)
    double m = 0.0;
    for (int k = 0; k < n; k++)
    {
        if (isnan(x[k])) return NAN;
        m = fmax(m, fabs(x[k]));
    }
    return m;
}

// Run to T with a step at most dt (shortened so that a whole number of steps
// reaches T); w (n^2 values, or NULL) receives the final vorticity. Returns the
// wall-clock seconds of the steps; *growth is max |w| over the run divided by
// its initial value (infinite for a non-finite field), checked every check
// steps (0: at the end only).
static double run(const study *s, int scheme, double dt, double T, int check, double *w_out, double *growth,
                  int *steps_out)
{
    int n = s->n, N = n * n, steps = (int)ceil(T / dt - 1E-9), k;
    solver_config c = config(s, scheme, T / steps);
    backend *b = backend_create(&c, s->gpu);
    mtrx *u, *v, *w;
    double w0, peak, seconds;
    struct timespec t0, t1;

    backend_fields(b, &u, &v, &w);
    random_initial_field(*w, *u, *v, c.dx, c.dy, 10.0, 0.5, 1);
    backend_set_fields(b, u, v, w);
    w0 = peak = max_abs(w->M, N);
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (k = 1; k <= steps; k++)
    {
        backend_step(b);
        if ((check > 0 && k % check == 0) || k == steps)
        {
            backend_fields(b, NULL, NULL, &w);
            double m = max_abs(w->M, N);
            peak = isfinite(m) ? fmax(peak, m) : INFINITY;
            if (!(peak <= 1E3 * w0)) break; // diverged: no need to go on
        }
    }
    backend_fields(b, NULL, NULL, &w);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    seconds = (double)(t1.tv_sec - t0.tv_sec) + 1E-9 * (double)(t1.tv_nsec - t0.tv_nsec);
    if (w_out) memcpy(w_out, w->M, (size_t)N * sizeof(double));
    *growth = peak / w0;
    if (steps_out) *steps_out = steps;
    backend_free(b);
    return seconds;
}

// Largest stable dt by bisection between a stable and an unstable step
static double stable_limit(const study *s, int scheme, double lo, double hi)
{
    double g;
    for (int it = 0; it < 14; it++)
    {
        double mid = sqrt(lo * hi);
        run(s, scheme, mid, 0.2, 20, NULL, &g, NULL);
        if (g < 2.0)
            lo = mid;
        else
            hi = mid;
    }
    return lo;
}

// Relative L2 error of w and the largest relative error of Z(k) over the
// shells 1 .. two thirds of Nyquist
static void errors(const study *s, spectra *sp, const double *w, const double *ref, const double *E_ref, double *e2,
                   double *espec)
{
    int n = s->n, N = n * n, bins = spectra_bins(sp), k;
    double num = 0.0, den = 0.0, *E = (double *)malloc((size_t)bins * sizeof(double));
    mtrx mw = initm(n, n), mu = initm(n, n), mv = initm(n, n);
    for (k = 0; k < N; k++)
    {
        num += (w[k] - ref[k]) * (w[k] - ref[k]);
        den += ref[k] * ref[k];
    }
    *e2 = sqrt(num / den);
    memcpy(mw.M, w, (size_t)N * sizeof(double));
    spectra_compute(sp, mu, mv, mw, NULL, E, NULL, NULL); // Z(k) of w (u, v unused for it)
    *espec = 0.0;
    for (k = 1; k < bins && k <= n / 3; k++)
        if (E_ref[k] > 0.0) *espec = fmax(*espec, fabs(E[k] / E_ref[k] - 1.0));
    free(E);
    freem(&mw);
    freem(&mu);
    freem(&mv);
}

int main(int argc, char **argv)
{
    int n = 256, gpu = 1, s1, k;
    double Re = 2E4, T = 0.1;
    for (k = 1; k < argc; k++)
    {
        char *end = NULL;
        int ok = 1;
        if (strcmp(argv[k], "--n") == 0 && k + 1 < argc)
        {
            long val = strtol(argv[++k], &end, 10);
            ok = *end == '\0' && val >= 16 && val <= 4096;
            n = (int)val;
        }
        else if (strcmp(argv[k], "--re") == 0 && k + 1 < argc)
        {
            Re = strtod(argv[++k], &end);
            ok = *end == '\0' && Re > 0.0 && isfinite(Re);
        }
        else if (strcmp(argv[k], "--cpu") == 0)
            gpu = 0;
        else
            ok = 0;
        if (!ok)
        {
            printf("Usage: %s [--n N] [--re RE] [--cpu]\n", argv[0]);
            return 1;
        }
    }
    study s;
    study_init(&s, n, Re, gpu);
    solver_config c = config(&s, 2, 1E-3);
    spectra *sp = spectra_setup(&c);
    int N = n * n, bins = spectra_bins(sp);
    double *ref = (double *)malloc((size_t)N * sizeof(double)), *w = (double *)malloc((size_t)N * sizeof(double));
    double *E_ref = (double *)malloc((size_t)bins * sizeof(double)), limit[5] = {0}, g; // E_ref: Z(k) of the reference
    static const int schemes[3] = {2, 3, 4};

    printf("Decaying turbulence, %d^2, Re = %g, order 6, skew-symmetric; %s\n", n, Re, gpu ? "GPU if available" : "CPU");
    // The bracket: the viscous limits of the schemes are far above the
    // advective ones here, so start from them and the Courant estimate
    for (s1 = 0; s1 < 3; s1++)
    {
        int sc = schemes[s1];
        double hi = max_stable_dt(&s.d2, &s.d2, Re, sc), lo = hi / 1E3;
        limit[sc] = stable_limit(&s, sc, lo, hi);
        printf("  %s: largest stable dt %.4e (viscous limit %.4e)\n", time_scheme_name(sc), limit[sc], hi);
    }
    printf("  AB3 / RK4 limit %.4f (imaginary-axis intervals: 0.7236 / 2.8284 = 0.2558)\n", limit[3] / limit[2]);
    printf("  RK46 / RK4 limit %.4f (imaginary-axis intervals: 3.8160 / 2.8284 = 1.3492)\n", limit[4] / limit[2]);

    // Reference: RK4 at 1/64 of its limit
    run(&s, 2, limit[2] / 64.0, T, 0, ref, &g, NULL);
    {
        mtrx mw = initm(n, n), mu = initm(n, n), mv = initm(n, n);
        memcpy(mw.M, ref, (size_t)N * sizeof(double));
        spectra_compute(sp, mu, mv, mw, NULL, E_ref, NULL, NULL);
        freem(&mw);
        freem(&mu);
        freem(&mv);
    }
    printf("\nError at t = %g against RK4 at 1/64 of its limit (columns: scheme, dt / own limit, dt, steps, "
           "seconds, relative L2 error of w, max relative error of Z(k) over K <= n/3):\n",
           T);
    for (s1 = 0; s1 < 3; s1++)
    {
        int sc = schemes[s1];
        for (int f = 1; f <= 16; f *= 2)
        {
            int steps;
            double e2, es, dt = 0.9 * limit[sc] / f, sec = run(&s, sc, dt, T, 0, w, &g, &steps);
            errors(&s, sp, w, ref, E_ref, &e2, &es);
            printf("  %-4s %6.4f %10.4e %6d %8.3f %10.3e %10.3e\n", time_scheme_name(sc), 0.9 / f, T / steps, steps,
                   sec, e2, es);
        }
    }
    spectra_free(sp);
    free(ref);
    free(w);
    free(E_ref);
    study_free(&s);
    return 0;
}
