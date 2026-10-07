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

static void usage(const char *prog)
{
    printf("Usage: %s [options]\n", prog);
    printf("  --nx N               grid points in x\n");
    printf("  --ny N               grid points in y\n");
    printf("  --n N                grid points in both x and y\n");
    printf("  --dt DT              time step\n");
    printf("  --tf TF              final time\n");
    printf("  --output-interval N  write VTK every N iterations (0 = never)\n");
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

    // Command-line overrides
    static struct option long_opts[] = {
        {"n", required_argument, 0, 'n'},
        {"nx", required_argument, 0, 'x'},
        {"ny", required_argument, 0, 'y'},
        {"dt", required_argument, 0, 'd'},
        {"tf", required_argument, 0, 'f'},
        {"output-interval", required_argument, 0, 'o'},
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
    double mem_needed = (double)nx * ny * 4.0 * (7.0 * (sizeof(double) + sizeof(int)) + sizeof(int)) + backend_host_memory(nx, ny, use_gpu);
    double mem_avail = available_memory();
    if (mem_avail >= 0. && mem_needed > mem_avail)
    {
        printf("** Error: a %d x %d grid needs about %.1f GB of memory; about %.1f GB is available **\n",
               nx, ny, mem_needed / 1E9, mem_avail / 1E9);
        return 1;
    }

    // Grid spacing: nodes 0 and nx-1 lie on the walls, so nx nodes span Lx
    // with nx-1 intervals
    double dx = (double)Lx / (nx - 1);
    double dy = (double)Ly / (ny - 1);

    double beta = sor_beta(nx, ny, dx, dy); // optimal SOR parameter

    printf("Grid: %d x %d | dt: %lf | tf: %lf\n", nx, ny, dt, tf);
#ifdef _OPENMP
    default_threads();
    printf("OpenMP threads: %d\n", omp_get_max_threads());
#endif
    printf("Poisson SOR parameter: %lf\n", beta);

    // Boundary conditions (Dirichlet): wall velocities on the left (1), right
    // (2), bottom (3) and top (4) walls; the top wall is the moving lid
    double ui = 0., vi = 0.;
    double u1 = 0., u2 = 0., u3 = 0., u4 = 1.;
    double v1 = 0., v2 = 0., v3 = 0., v4 = 0.;
    wall_bc bc = {{u1, u2, u3, u4}, {v1, v2, v3, v4}};

    // Build sparse 1D operators then free them after Kronecker
    smtrx sd_x = SDiff1(nx, order, dx);
    smtrx sd_y = SDiff1(ny, order, dy);
    smtrx sd_x2 = SDiff2(nx, order, dx);
    smtrx sd_y2 = SDiff2(ny, order, dy);
    smtrx sIx = seye(nx);
    smtrx sIy = seye(ny);

    // Stability checks, before the 2D operators are built so that a run that
    // cannot work fails at once. The fastest wall is the velocity scale.
    double u_max = 0.;
    for (i = 0; i < 4; i++)
    {
        if (fabs(bc.u[i]) > u_max) u_max = fabs(bc.u[i]);
        if (fabs(bc.v[i]) > u_max) u_max = fabs(bc.v[i]);
    }
    dt_limits lim = time_step_limits(&sd_x2, &sd_y2, fmin(dx, dy), Re, u_max, max_co, time_scheme);
    if (dt > lim.accept)
    {
        printf("** Error: dt = %g is too large; use --dt %.3g or less. Limits: Courant number <= %g "
               "gives dt <= %.3g, the viscous stability limit of the %s scheme for this grid, Re and "
               "order gives dt <= %.3g",
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
    cfg.poisson_tol = poisson_tol;
    cfg.beta = beta;
    cfg.bc = bc;
    cfg.DX = &DX;
    cfg.DY = &DY;
    cfg.DX2 = &DX2;
    cfg.DY2 = &DY2;
    cfg.t0 = 0.0;
    cfg.vorticity_source = NULL;
    cfg.source_data = NULL;

    int it_max = (int)((tf / dt) - 1);

    // Backend selection: GPU when built with CUDA=1 and a device is usable
    backend *solver = backend_create(&cfg, use_gpu);

    // Initial condition, on the solver's host fields. Fields are ny rows (y)
    // of nx values (x).
    mtrx *u, *v, *w;
    backend_fields(solver, &u, &v, &w);
    for (i = 1; i < ny - 1; i++)
        for (j = 1; j < nx - 1; j++)
        {
            MAt(*u, i, j) = ui;
            MAt(*v, i, j) = vi;
        }
    backend_set_fields(solver, u, v, w);
    printf("Backend: %s\n", backend_name(solver));
    if (backend_device(solver)) printf("CUDA device: %s\n", backend_device(solver));

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

        if (output_interval > 0 && t % output_interval == 0)
        {
            backend_fields(solver, NULL, NULL, &w);
            printvtk(*w, "vorticity", dx, dy);
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t_end);
    double elapsed = (double)(t_end.tv_sec - t_start.tv_sec) + 1E-9 * (double)(t_end.tv_nsec - t_start.tv_nsec);

    backend_fields(solver, &u, &v, &w);

    // Re-apply wall BCs before sampling centerline. This writes to the read
    // view without backend_set_fields(), which is fine only because the
    // solver is not stepped again.
    apply_wall_bc(*u, *v, &bc);

    // Write centerline profiles and compare against Ghia et al. (1982)
    print_centerline(*u, *v, nx, ny, dx, dy);

    printf("Wall-clock time: %.3lf s total | %.4lf ms per step (%d steps, %s)\n",
           elapsed, 1E3 * elapsed / (it_max + 1), it_max + 1, backend_name(solver));

    backend_free(solver);
    freesm(DX);
    freesm(DY);
    freesm(DX2);
    freesm(DY2);

    printf("Simulation complete!\n");
    return 0;
}
