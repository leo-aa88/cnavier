#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <limits.h>
#include <getopt.h>
#include "linearalg.h"
#include "finitediff.h"
#include "utils.h"
#include "poisson.h"
#include "fluiddyn.h"
#include "threads.h"
#include "backend.h"
#include "diagnostics.h"
#ifdef _OPENMP
#include <omp.h>
#endif

// Largest supported grid: keeps nx*ny, the CSR non-zero count (up to 7 per
// row) and the FFT extension 4*(nx-1)*(ny-1) within int range
#define MAX_GRID 16384

// Parse a whole string as an int. Returns 0 if it is not one.
static int parse_int(const char *s, int *out)
{
    char *end;
    long val;

    errno = 0;
    val = strtol(s, &end, 10);
    if (errno != 0 || end == s || *end != '\0' || val < INT_MIN || val > INT_MAX)
        return 0;
    *out = (int)val;
    return 1;
}

// Parse a whole string as a finite double. Returns 0 if it is not one.
static int parse_double(const char *s, double *out)
{
    char *end;
    double val;

    errno = 0;
    val = strtod(s, &end);
    if (errno != 0 || end == s || *end != '\0' || !isfinite(val))
        return 0;
    *out = val;
    return 1;
}

// Flow cases: the lid-driven cavity has four walls, the others are doubly
// periodic on the unit square
enum flow_case
{
    CASE_CAVITY,
    CASE_TAYLOR_GREEN,
    CASE_SHEAR_LAYER,
    CASE_KOLMOGOROV,
    CASE_FORCED,
    CASE_DECAYING
};

#define CASE_PI 3.14159265358979323846

// Taylor-Green vortex: psi = sin(kx) sin(ky) e^(-2 k^2 t / Re) / k, k = 2 pi,
// an exact solution (the nonlinear term vanishes), with |u| <= 1
static void taylor_green(double x, double y, double t, double Re, double *w, double *u, double *v)
{
    double k = 2.0 * CASE_PI, decay = exp(-2.0 * k * k * t / Re);
    *w = 2.0 * k * sin(k * x) * sin(k * y) * decay;
    *u = sin(k * x) * cos(k * y) * decay;
    *v = -cos(k * x) * sin(k * y) * decay;
}

// Double shear layer (Bell, Colella and Glaz 1989): two tanh layers of
// thickness 1/30 at y = 1/4 and y = 3/4 with a small sinusoidal v that makes
// them roll up into vortices
static void shear_layer(double x, double y, double *w, double *u, double *v)
{
    const double delta = 1.0 / 30.0, eps = 0.05;
    double s = y <= 0.5 ? (y - 0.25) / delta : (0.75 - y) / delta, sech = 1.0 / cosh(s);
    *u = tanh(s);
    *v = eps * sin(2.0 * CASE_PI * x);
    // w = dv/dx - du/dy
    *w = 2.0 * CASE_PI * eps * cos(2.0 * CASE_PI * x) - (y <= 0.5 ? 1.0 : -1.0) * sech * sech / delta;
}

static void usage(const char *prog)
{
    printf("Usage: %s [options]\n", prog);
    printf("  --nx N               grid points in x\n");
    printf("  --ny N               grid points in y\n");
    printf("  --n N                grid points in both x and y\n");
    printf("  --dt DT              time step\n");
    printf("  --tf TF              final time\n");
    printf("  --output-interval N  write VTK every N iterations (0 = never)\n");
    printf("  --re RE              Reynolds number\n");
    printf("  --order N            order of the finite differences: 2, 4 or 6 (default 6),\n"
           "                       or compact6 (periodic cases: Lele's sixth-order compact)\n");
    printf("  --integrals-interval N  write E, Z, P to output/integrals.csv every N steps\n"
           "                       (default 1; 0 = never)\n");
    printf("  --poisson-order N    2 (5-point, default) or 4 (compact 9-point) Poisson\n"
           "                       operator of the FFT solver with walls\n");
    printf("  --wall-closure NAME  wall vorticity: velocity (dv/dx - du/dy, default) or\n"
           "                       briley (third order, from the stream function)\n");
    printf("  --velocity-order N   2 (default) or 4: order of the derivative rows next to\n"
           "                       the walls that give u, v from the stream function\n");
    printf("  --case NAME          cavity (default; lid-driven, four walls), or on a\n"
           "                       doubly periodic unit square: taylor-green (decaying\n"
           "                       vortex, compared with the exact solution at the end)\n"
           "                       shear-layer (double shear layer that rolls up),\n"
           "                       kolmogorov (Kolmogorov forcing from a small random\n"
           "                       perturbation), forced (random forcing from rest) or\n"
           "                       decaying (decaying turbulence from random phases)\n");
    printf("  --drag ALPHA         linear drag -ALPHA w (default 0; 0.1 for forced)\n");
    printf("  --kolmogorov-amp A   Kolmogorov body force A sin(k y) in x (default 1 for\n"
           "  --kolmogorov-n N     kolmogorov, else 0), k = 2 pi N (default 4)\n");
    printf("  --forcing-rate EPS   random forcing injecting energy at rate EPS (default\n"
           "                       0.1 for forced, else 0) on the wavenumbers\n"
           "  --forcing-k KF       |k| / 2 pi within KF +- DK (defaults 8 and 1)\n"
           "  --forcing-width DK\n");
    printf("  --advection NAME     periodic cases: nonlinear term in advective or skew\n"
           "                       (skew-symmetric, conserves enstrophy) form (default:\n"
           "                       skew for forced and decaying, else advective)\n");
    printf("  --hyperviscosity NU  hyperviscosity -NU (-lap)^P w (default 0)\n"
           "  --hyper-order P      its order P >= 2 (default 4)\n");
    printf("  --hypodrag ALPHA     large-scale drag -ALPHA psi (default 0)\n");
    printf("  --peak-k K0          decaying: the initial spectrum envelope\n"
           "                       (k/k0)^4 exp(-2 (k/k0)^2) peaks at |k| / 2 pi = K0\n"
           "                       (default 10); energy 1/2\n");
    printf("  --seed S             seed of the random forcing and initial field (default 1)\n");
    printf("  --spectrum-interval N  periodic cases: write output/spectrum-1-<n>.csv every N\n"
           "                       steps (default: with the VTK frames; 0 = never)\n");
#ifdef USE_CUDA
    printf("  --cpu                run on the CPU instead of the GPU\n");
#endif
    printf("  --help               show this message\n");
    printf("Unset options keep the defaults at the top of src/main.c.\n");
}

int main(int argc, char *argv[])
{
    int i, j, t;

    // Physical parameters
    double Re = 100.;
    int Lx = 1;
    int Ly = 1;

    // Numerical parameters
    int nx = 64;
    int ny = 64;
    double dt = 0.005;
    double tf = 30;
    double max_co = 1.;
    int order = 6;
    int poisson_max_it = 10000;
    double poisson_tol = 1E-3;
    int output_interval = 20;
    int poisson_type = 3; // 1=Gauss-Seidel  2=SOR  3=FFT (direct, exact)
    int time_scheme = 2;  // 1=Euler  2=RK4
    int use_gpu = 1;      // only meaningful when built with CUDA=1
    enum flow_case flow = CASE_CAVITY;
    int poisson_order = 2;      // 2 = 5-point, 4 = compact 9-point (FFT solver, walls)
    int wall_closure = 0;       // 0 = dv/dx - du/dy, 1 = Briley's formula from psi
    int velocity_order = 2;     // 4 = fourth-order rows next to the walls for u, v
    int integrals_interval = 1; // steps between lines of output/integrals.csv, 0 = none
    // Forcing and drag (forcing.h); negative: not given, take the case's default
    double drag = -1., kolmogorov_amp = -1., forcing_rate = -1., forcing_k = 8., forcing_width = 1.;
    int kolmogorov_n = 4;
    double hyperviscosity = 0., hypodrag = 0., peak_k = 10.;
    int hyper_order = 4;
    int spectrum_interval = -1; // negative: with the VTK frames
    int advection = -1;         // 0 = advective, 1 = skew-symmetric nonlinear term; -1: the case's default
    unsigned long long seed = 1;
    int forcing_given = 0;

    // Command-line overrides
    static struct option long_opts[] = {
        {"n", required_argument, 0, 'n'},
        {"nx", required_argument, 0, 'x'},
        {"ny", required_argument, 0, 'y'},
        {"dt", required_argument, 0, 'd'},
        {"tf", required_argument, 0, 'f'},
        {"output-interval", required_argument, 0, 'o'},
        {"re", required_argument, 0, 'r'},
        {"case", required_argument, 0, 'k'},
        {"poisson-order", required_argument, 0, 'p'},
        {"wall-closure", required_argument, 0, 'w'},
        {"velocity-order", required_argument, 0, 'u'},
        {"integrals-interval", required_argument, 0, 'i'},
        {"drag", required_argument, 0, 'a'},
        {"kolmogorov-amp", required_argument, 0, 'K'},
        {"kolmogorov-n", required_argument, 0, 'N'},
        {"forcing-rate", required_argument, 0, 'e'},
        {"forcing-k", required_argument, 0, 'q'},
        {"forcing-width", required_argument, 0, 'W'},
        {"seed", required_argument, 0, 's'},
        {"hyperviscosity", required_argument, 0, 'H'},
        {"hyper-order", required_argument, 0, 'P'},
        {"hypodrag", required_argument, 0, 'D'},
        {"peak-k", required_argument, 0, 'g'},
        {"spectrum-interval", required_argument, 0, 'S'},
        {"advection", required_argument, 0, 'A'},
        {"order", required_argument, 0, 'O'},
        {"cpu", no_argument, 0, 'c'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}};
    int opt, opt_index;
    while ((opt = getopt_long(argc, argv, "", long_opts, &opt_index)) != -1)
    {
        int ok = 1;
        switch (opt)
        {
        case 'n':
            ok = parse_int(optarg, &nx);
            ny = nx;
            break;
        case 'x':
            ok = parse_int(optarg, &nx);
            break;
        case 'y':
            ok = parse_int(optarg, &ny);
            break;
        case 'd':
            ok = parse_double(optarg, &dt);
            break;
        case 'f':
            ok = parse_double(optarg, &tf);
            break;
        case 'o':
            ok = parse_int(optarg, &output_interval);
            break;
        case 'r':
            ok = parse_double(optarg, &Re) && Re > 0.;
            break;
        case 'k':
            if (strcmp(optarg, "cavity") == 0)
                flow = CASE_CAVITY;
            else if (strcmp(optarg, "taylor-green") == 0)
                flow = CASE_TAYLOR_GREEN;
            else if (strcmp(optarg, "shear-layer") == 0)
                flow = CASE_SHEAR_LAYER;
            else if (strcmp(optarg, "kolmogorov") == 0)
                flow = CASE_KOLMOGOROV;
            else if (strcmp(optarg, "forced") == 0)
                flow = CASE_FORCED;
            else if (strcmp(optarg, "decaying") == 0)
                flow = CASE_DECAYING;
            else
                ok = 0;
            break;
        case 'p':
            ok = parse_int(optarg, &poisson_order) && (poisson_order == 2 || poisson_order == 4);
            break;
        case 'w':
            if (strcmp(optarg, "velocity") == 0)
                wall_closure = 0;
            else if (strcmp(optarg, "briley") == 0)
                wall_closure = 1;
            else
                ok = 0;
            break;
        case 'u':
            ok = parse_int(optarg, &velocity_order) && (velocity_order == 2 || velocity_order == 4);
            break;
        case 'i':
            ok = parse_int(optarg, &integrals_interval) && integrals_interval >= 0;
            break;
        case 'a':
            ok = parse_double(optarg, &drag) && drag >= 0.;
            forcing_given = 1;
            break;
        case 'K':
            ok = parse_double(optarg, &kolmogorov_amp);
            forcing_given = 1;
            break;
        case 'N':
            ok = parse_int(optarg, &kolmogorov_n) && kolmogorov_n >= 1;
            forcing_given = 1;
            break;
        case 'e':
            ok = parse_double(optarg, &forcing_rate) && forcing_rate >= 0.;
            forcing_given = 1;
            break;
        case 'q':
            ok = parse_double(optarg, &forcing_k) && forcing_k > 0.;
            forcing_given = 1;
            break;
        case 'W':
            ok = parse_double(optarg, &forcing_width) && forcing_width > 0.;
            forcing_given = 1;
            break;
        case 'H':
            ok = parse_double(optarg, &hyperviscosity) && hyperviscosity >= 0.;
            forcing_given = 1;
            break;
        case 'P':
            ok = parse_int(optarg, &hyper_order) && hyper_order >= 2 && hyper_order <= 8;
            forcing_given = 1;
            break;
        case 'D':
            ok = parse_double(optarg, &hypodrag) && hypodrag >= 0.;
            forcing_given = 1;
            break;
        case 'g':
            ok = parse_double(optarg, &peak_k) && peak_k > 0.;
            break;
        case 'O':
            if (strcmp(optarg, "compact6") == 0)
                order = FD_COMPACT6;
            else
                ok = parse_int(optarg, &order) && (order == 2 || order == 4 || order == 6);
            break;
        case 'A':
            if (strcmp(optarg, "advective") == 0)
                advection = 0;
            else if (strcmp(optarg, "skew") == 0)
                advection = 1;
            else
                ok = 0;
            break;
        case 'S':
            ok = parse_int(optarg, &spectrum_interval) && spectrum_interval >= 0;
            break;
        case 's':
        {
            int sv = 0;
            ok = parse_int(optarg, &sv) && sv >= 0;
            seed = (unsigned long long)sv;
            break;
        }
        case 'c':
            use_gpu = 0;
            break;
        case 'h':
            usage(argv[0]);
            return 0;
        default:
            usage(argv[0]);
            return 1;
        }
        if (!ok)
        {
            printf("** Error: invalid value '%s' for --%s **\n", optarg, long_opts[opt_index].name);
            return 1;
        }
    }
    if (optind < argc)
    {
        printf("** Error: unexpected argument '%s' **\n", argv[optind]);
        usage(argv[0]);
        return 1;
    }

    // Reject what the solver cannot handle. The tests on dt and tf are
    // written so that a NaN fails them.
    if (nx < 8 || nx > MAX_GRID || ny < 8 || ny > MAX_GRID)
    {
        printf("** Error: grid size must be between 8 and %d **\n", MAX_GRID);
        return 1;
    }
    if (!(dt > 0.) || !isfinite(dt) || !(tf > 0.) || !isfinite(tf))
    {
        printf("** Error: dt and tf must be positive and finite **\n");
        return 1;
    }
    if (output_interval < 0)
    {
        printf("** Error: output interval must not be negative **\n");
        return 1;
    }
    // The number of timesteps is checked as a double, before it becomes an int
    if (!(tf / dt >= 1.) || tf / dt > (double)INT_MAX)
    {
        printf("** Error: tf/dt gives %g timesteps; it must be between 1 and %d **\n",
               tf / dt, INT_MAX);
        return 1;
    }

    // Host memory estimate: four CSR operators with at most 7 non-zeros per
    // row, plus the arrays of the backend asked for (the CPU workspace, or on
    // the GPU only host copies of the fields), compared with the memory
    // available now. A run that would not fit stops here with a message
    // instead of being killed by the kernel once the pages are touched. Other
    // processes can still take memory after this check, so it is a guard
    // against the clear cases only; backend_create() checks again for the
    // backend it actually uses.
    // The compact operators are banded circulants about 65 entries wide on
    // average (81 for the first derivatives, 49 for the second)
    double row_entries = order == FD_COMPACT6 ? fmin(65.0, (double)(nx > ny ? nx : ny)) : 7.0;
    double mem_needed = (double)nx * ny * (velocity_order == 4 ? 6.0 : 4.0) * (row_entries * (sizeof(double) + sizeof(int)) + sizeof(int)) + backend_host_memory(nx, ny, use_gpu);
    double mem_avail = available_memory();
    if (mem_avail >= 0. && mem_needed > mem_avail)
    {
        printf("** Error: a %d x %d grid needs about %.1f GB of memory; about %.1f GB is available **\n",
               nx, ny, mem_needed / 1E9, mem_avail / 1E9);
        return 1;
    }

    int periodic = flow != CASE_CAVITY;
    if (periodic && poisson_type != 3)
    {
        printf("** Error: the periodic cases need the FFT Poisson solver (poisson_type = 3) **\n");
        return 1;
    }
    if (!periodic && forcing_given)
    {
        printf("** Error: the forcing and drag options apply to the periodic cases **\n");
        return 1;
    }
    if (!periodic && spectrum_interval > 0)
    {
        printf("** Error: --spectrum-interval applies to the periodic cases **\n");
        return 1;
    }
    if (!periodic && advection == 1)
    {
        printf("** Error: --advection skew applies to the periodic cases **\n");
        return 1;
    }
    // Turbulence cases default to the skew-symmetric form: in under-resolved
    // turbulence the advective form makes enstrophy at the grid cutoff
    if (advection < 0) advection = flow == CASE_FORCED || flow == CASE_DECAYING;
    if (spectrum_interval < 0) spectrum_interval = output_interval;
    // The cases' defaults for what was not given
    if (drag < 0.) drag = flow == CASE_FORCED ? 0.1 : 0.;
    if (kolmogorov_amp < 0.) kolmogorov_amp = flow == CASE_KOLMOGOROV ? 1. : 0.;
    if (forcing_rate < 0.) forcing_rate = flow == CASE_FORCED ? 0.1 : 0.;
    if (periodic && (poisson_order != 2 || wall_closure != 0 || velocity_order != 2))
    {
        printf("** Error: --poisson-order, --wall-closure and --velocity-order apply to walls; "
               "the periodic cases have none **\n");
        return 1;
    }
    if (order == FD_COMPACT6 && !periodic)
    {
        printf("** Error: --order compact6 applies to the periodic cases **\n");
        return 1;
    }
    if (velocity_order == 4 && (nx < 10 || ny < 10))
    {
        printf("** Error: --velocity-order 4 needs at least 10 grid points in x and y **\n");
        return 1;
    }
    if (poisson_order == 4 && poisson_type != 3)
    {
        printf("** Error: --poisson-order 4 needs the FFT Poisson solver (poisson_type = 3) **\n");
        return 1;
    }

    // Grid spacing. With walls, nodes 0 and nx-1 lie on the walls, so nx
    // nodes span Lx with nx-1 intervals. On a periodic grid node nx would be
    // node 0 again, so nx nodes span Lx with nx intervals.
    double dx = (double)Lx / (periodic ? nx : nx - 1);
    double dy = (double)Ly / (periodic ? ny : ny - 1);

    double beta = sor_beta(nx, ny, dx, dy); // optimal SOR parameter

    static const char *case_names[] = {"lid-driven cavity", "Taylor-Green vortex (periodic)",
                                       "double shear layer (periodic)", "Kolmogorov flow (periodic)",
                                       "randomly forced flow (periodic)", "decaying turbulence (periodic)"};
    printf("Case: %s | Re: %g\n", case_names[flow], Re);
    if (drag > 0. || kolmogorov_amp != 0. || forcing_rate > 0.)
        printf("Forcing: drag %g | Kolmogorov A %g, k = 2 pi x %d | random rate %g, |k|/2pi in %g +- %g, seed %llu\n",
               drag, kolmogorov_amp, kolmogorov_n, forcing_rate, forcing_k, forcing_width, seed);
    if (periodic) printf("Nonlinear term: %s\n", advection == 1 ? "skew-symmetric" : "advective");
    if (hyperviscosity > 0. || hypodrag > 0.)
        printf("Damping: hyperviscosity %g, order %d | hypodrag %g\n", hyperviscosity, hyper_order, hypodrag);
    printf("Grid: %d x %d | dt: %lf | tf: %lf | derivatives: %s\n", nx, ny, dt, tf,
           order == FD_COMPACT6 ? "compact, order 6" : order == 4 ? "order 4"
                                                   : order == 2   ? "order 2"
                                                                  : "order 6");
#ifdef _OPENMP
    default_threads();
    printf("OpenMP threads: %d\n", omp_get_max_threads());
#endif
    if (!periodic) printf("Poisson SOR parameter: %lf\n", beta);

    // Boundary conditions (Dirichlet): wall velocities on the left (1), right
    // (2), bottom (3) and top (4) walls; the top wall is the moving lid
    double ui = 0., vi = 0.;
    double u1 = 0., u2 = 0., u3 = 0., u4 = 1.;
    double v1 = 0., v2 = 0., v3 = 0., v4 = 0.;
    wall_bc bc = {{u1, u2, u3, u4}, {v1, v2, v3, v4}};

    // Build sparse 1D operators then free them after Kronecker
    smtrx sd_x = periodic ? SDiff1_periodic(nx, order, dx) : SDiff1(nx, order, dx);
    smtrx sd_y = periodic ? SDiff1_periodic(ny, order, dy) : SDiff1(ny, order, dy);
    smtrx sd_x2 = periodic ? SDiff2_periodic(nx, order, dx) : SDiff2(nx, order, dx);
    smtrx sd_y2 = periodic ? SDiff2_periodic(ny, order, dy) : SDiff2(ny, order, dy);
    smtrx sIx = seye(nx);
    smtrx sIy = seye(ny);

    // Stability checks, before the 2D operators are built so that a run that
    // cannot work fails at once. The velocity scale is the fastest wall, or
    // for the periodic cases 1 (the largest initial speed of the decaying
    // cases), or the speed the forcing drives: the laminar Kolmogorov speed
    // A / (nu k^2 + drag), and three times the r.m.s. speed sqrt(eps / drag)
    // where random forcing balances drag. Turbulence can exceed these; the
    // check guards against the clear cases only. Decaying turbulence starts
    // with r.m.s. speed 1 and peaks of a few times that.
    double u_max = periodic ? (flow == CASE_DECAYING ? 4.0 : 1.0) : 0.;
    if (kolmogorov_amp != 0.)
    {
        double kk = 2.0 * CASE_PI * kolmogorov_n;
        u_max = fmax(u_max, fabs(kolmogorov_amp) / (kk * kk / Re + drag));
    }
    if (forcing_rate > 0. && drag > 0.) u_max = fmax(u_max, 3.0 * sqrt(forcing_rate / drag));
    for (i = 0; i < 4 && !periodic; i++)
    {
        if (fabs(bc.u[i]) > u_max) u_max = fabs(bc.u[i]);
        if (fabs(bc.v[i]) > u_max) u_max = fabs(bc.v[i]);
    }
    forcing_config damping = {0};
    damping.drag = drag;
    damping.hyperviscosity = hyperviscosity;
    damping.hyper_order = hyper_order;
    damping.hypodrag = hypodrag;
    dt_limits lim = time_step_limits_forced(&sd_x2, &sd_y2, fmin(dx, dy), Re, u_max, max_co, time_scheme, &damping);
    if (dt > lim.accept)
    {
        printf("** Error: dt = %g is too large; use --dt %.3g or less. Limits: Courant number <= %g "
               "gives dt <= %.3g, the stability limit of the %s scheme for the viscous and damping "
               "terms on this grid gives dt <= %.3g",
               dt, round_down_3(lim.suggest), max_co, round_down_3(lim.courant),
               time_scheme == 1 ? "Euler" : "RK4", round_down_3(lim.viscous));
        if (time_scheme == 1)
            printf(", forward Euler's limit for centered advection gives dt <= %.3g",
                   round_down_3(lim.advection));
        printf(" **\n");
        exit(1);
    }
    // Forward Euler's limit for centered advection is conservative here, so it
    // is a warning, not an error
    if (dt > lim.advection)
        printf("** Warning: dt = %g is above 2/(Re u^2) = %.3g, the stability limit of forward Euler "
               "for centered advection at the wall speed. It assumes that speed everywhere and is "
               "conservative for the cavity (at Re = 1000 runs stayed stable up to about 3x it), so "
               "the run goes ahead, but it may diverge **\n",
               dt, round_down_3(lim.advection));

    // Sparse 2D operators: DX = I_y x d_x,  DY = d_y x I_x
    smtrx DX = skronecker(sIy, sd_x);
    smtrx DY = skronecker(sd_y, sIx);
    smtrx DX2 = skronecker(sIy, sd_x2);
    smtrx DY2 = skronecker(sd_y2, sIx);
    smtrx DXv = {0}, DYv = {0};
    if (velocity_order == 4)
    {
        smtrx vx = SDiff1_wall4(nx, order, dx), vy = SDiff1_wall4(ny, order, dy);
        DXv = skronecker(sIy, vx);
        DYv = skronecker(vy, sIx);
        freesm(vx);
        freesm(vy);
    }

    freesm(sd_x);
    freesm(sd_y);
    freesm(sd_x2);
    freesm(sd_y2);
    freesm(sIx);
    freesm(sIy);

    // Everything that defines the run, for whichever backend runs it
    solver_config cfg;
    cfg.nx = nx;
    cfg.ny = ny;
    cfg.dx = dx;
    cfg.dy = dy;
    cfg.Re = Re;
    cfg.dt = dt;
    cfg.time_scheme = time_scheme;
    cfg.poisson_type = poisson_type;
    cfg.poisson_max_it = poisson_max_it;
    cfg.poisson_order = poisson_order;
    cfg.wall_closure = wall_closure;
    cfg.poisson_tol = poisson_tol;
    cfg.beta = beta;
    cfg.periodic = periodic;
    cfg.advection = advection;
    cfg.bc = bc;
    cfg.DX = &DX;
    cfg.DY = &DY;
    cfg.DX2 = &DX2;
    cfg.DY2 = &DY2;
    cfg.DXv = velocity_order == 4 ? &DXv : NULL;
    cfg.DYv = velocity_order == 4 ? &DYv : NULL;
    cfg.t0 = 0.0;
    cfg.vorticity_source = NULL;
    cfg.source_data = NULL;
    cfg.forcing = (forcing_config){0};
    cfg.forcing.drag = drag;
    cfg.forcing.kolmogorov_amp = kolmogorov_amp;
    cfg.forcing.kolmogorov_n = kolmogorov_n;
    cfg.forcing.random_rate = forcing_rate;
    cfg.forcing.random_kf = forcing_k;
    cfg.forcing.random_dk = forcing_width;
    cfg.forcing.random_seed = seed;
    cfg.forcing.hyperviscosity = hyperviscosity;
    cfg.forcing.hyper_order = hyper_order;
    cfg.forcing.hypodrag = hypodrag;

    int it_max = (int)((tf / dt) - 1);

    // Backend selection: GPU when built with CUDA=1 and a device is usable
    backend *solver = backend_create(&cfg, use_gpu);

    // Initial condition, on the solver's host fields. Fields are ny rows (y)
    // of nx values (x).
    mtrx *u, *v, *w;
    backend_fields(solver, &u, &v, &w);
    if (flow == CASE_CAVITY)
        for (i = 1; i < ny - 1; i++)
            for (j = 1; j < nx - 1; j++)
            {
                MAt(*u, i, j) = ui;
                MAt(*v, i, j) = vi;
            }
    else if (flow == CASE_DECAYING)
        random_initial_field(*w, *u, *v, dx, dy, peak_k, 0.5, seed);
    else
        for (i = 0; i < ny; i++)
            for (j = 0; j < nx; j++)
            {
                if (flow == CASE_TAYLOR_GREEN)
                    taylor_green(j * dx, i * dy, 0.0, Re, &MAt(*w, i, j), &MAt(*u, i, j), &MAt(*v, i, j));
                else if (flow == CASE_SHEAR_LAYER)
                    shear_layer(j * dx, i * dy, &MAt(*w, i, j), &MAt(*u, i, j), &MAt(*v, i, j));
                else
                {
                    // Kolmogorov: a small random vorticity perturbation, from
                    // which the instability grows; forced: rest
                    MAt(*w, i, j) = 0.0;
                    if (flow == CASE_KOLMOGOROV)
                    {
                        seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
                        MAt(*w, i, j) = 1E-3 * ((double)(seed >> 11) / 9007199254740992.0 - 0.5);
                    }
                    MAt(*u, i, j) = MAt(*v, i, j) = 0.0;
                }
            }
    backend_set_fields(solver, u, v, w);
    printf("Backend: %s\n", backend_name(solver));
    if (backend_device(solver)) printf("CUDA device: %s\n", backend_device(solver));

    // Energy, enstrophy and palinstrophy every integrals_interval steps, and
    // on periodic grids the spectra with every VTK frame
    FILE *integrals = NULL;
    flow_integrals fi;
    if (integrals_interval > 0)
    {
        if (!(integrals = fopen("./output/integrals.csv", "w")))
        {
            printf("\nError while opening file\n");
            exit(1);
        }
        fi = backend_integrals(solver);
        fprintf(integrals, "step,t,E,Z,P,I,I_disc\n0,0,%.17g,%.17g,%.17g,%.17g,%.17g\n", fi.E, fi.Z, fi.P, fi.I,
                fi.I_disc);
    }
    spectra *spec = periodic && spectrum_interval > 0 ? spectra_setup(&cfg) : NULL;
    double *spec_out = spec ? (double *)malloc((size_t)SPECTRA_COLUMNS * spectra_bins(spec) * sizeof(double)) : NULL;
    if (spec && !spec_out)
    {
        printf("** Error: insufficient memory **\n");
        exit(1);
    }

    struct timespec t_start, t_end;
    clock_gettime(CLOCK_MONOTONIC, &t_start);

    // Main time loop
    for (t = 0; t <= it_max; t++)
    {
        double cmax, cmin;

        // Boundary conditions, time advancement + Poisson solve
        backend_step(solver);

        // Continuity check: du/dx + dv/dy ~ 0
        backend_continuity(solver, &cmax, &cmin);

        printf("Iteration: %d | Time: %.4lf | Progress: %.2lf%%\n",
               t, (double)t * dt, it_max > 0 ? (double)100 * t / it_max : 100.);
        printf("Continuity max: %E | min: %E\n", cmax, cmin);

        if (integrals && (t + 1) % integrals_interval == 0)
        {
            fi = backend_integrals(solver);
            fprintf(integrals, "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n", t + 1, (double)(t + 1) * dt, fi.E, fi.Z,
                    fi.P, fi.I, fi.I_disc);
        }

        if (output_interval > 0 && t % output_interval == 0)
        {
            backend_fields(solver, NULL, NULL, &w);
            printvtk(*w, "vorticity", dx, dy);
        }
        if (spec && t % spectrum_interval == 0)
        {
            backend_spectra(solver, spec, spec_out);
            spectra_write_frame(spec, spec_out, (double)(t + 1) * dt);
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t_end);
    double elapsed = (double)(t_end.tv_sec - t_start.tv_sec) + 1E-9 * (double)(t_end.tv_nsec - t_start.tv_nsec);

    if (integrals) fclose(integrals);
    spectra_free(spec);
    free(spec_out);
    backend_fields(solver, &u, &v, &w);

    if (flow == CASE_CAVITY)
    {
        // Re-apply wall BCs before sampling centerline. This writes to the
        // read view without backend_set_fields(), which is fine only because
        // the solver is not stepped again.
        apply_wall_bc(*u, *v, &bc);

        // Write centerline profiles and compare against Ghia et al. (1982)
        print_centerline(*u, *v, nx, ny, dx, dy);
    }
    else if (flow == CASE_TAYLOR_GREEN)
    {
        double t_end = (double)(it_max + 1) * dt, err = 0.0, peak = 0.0;
        for (i = 0; i < ny; i++)
            for (j = 0; j < nx; j++)
            {
                double we, ue, ve;
                taylor_green(j * dx, i * dy, t_end, Re, &we, &ue, &ve);
                err = fmax(err, fabs(MAt(*w, i, j) - we));
                peak = fmax(peak, fabs(we));
            }
        printf("Taylor-Green vortex at t = %g: max |w - exact| = %E, relative to max |w| %E\n",
               t_end, err, err / peak);
    }

    printf("Wall-clock time: %.3lf s total | %.4lf ms per step (%d steps, %s)\n",
           elapsed, 1E3 * elapsed / (it_max + 1), it_max + 1, backend_name(solver));

    backend_free(solver);
    freesm(DX);
    freesm(DY);
    freesm(DX2);
    freesm(DY2);
    if (velocity_order == 4)
    {
        freesm(DXv);
        freesm(DYv);
    }

    printf("Simulation complete!\n");
    return 0;
}
