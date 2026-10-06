#define _GNU_SOURCE // sched_getaffinity
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <limits.h>
#include <getopt.h>
#include <sched.h>
#include "linearalg.h"
#include "finitediff.h"
#include "utils.h"
#include "poisson.h"
#include "fluiddyn.h"
#ifdef USE_CUDA
#include "cudasolver.h"
#endif
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

#ifdef _OPENMP
// Number of distinct physical cores among the CPUs this process may run on:
// the (package, core) pairs of the CPUs in its affinity mask, read from Linux
// sysfs. This handles hybrid CPUs, whose cores have different numbers of
// hardware threads, and affinity masks that cover part of the machine.
// Returns 0 if it cannot be determined.
static int physical_cores(void)
{
    cpu_set_t set;
    int cpu, count = 0, n = 0;
    long seen[CPU_SETSIZE];

    if (sched_getaffinity(0, sizeof(set), &set) != 0) return 0;
    for (cpu = 0; cpu < CPU_SETSIZE; cpu++)
    {
        char path[128];
        long core = -1, package = -1, key;
        int k, found = 0;
        FILE *f;

        if (!CPU_ISSET(cpu, &set)) continue;
        snprintf(path, sizeof(path), "/sys/devices/system/cpu/cpu%d/topology/core_id", cpu);
        if ((f = fopen(path, "r"))) { if (fscanf(f, "%ld", &core) != 1) core = -1; fclose(f); }
        snprintf(path, sizeof(path), "/sys/devices/system/cpu/cpu%d/topology/physical_package_id", cpu);
        if ((f = fopen(path, "r"))) { if (fscanf(f, "%ld", &package) != 1) package = -1; fclose(f); }
        if (core < 0 || package < 0) return 0;

        key = package * 1000000 + core;
        for (k = 0; k < n; k++)
            if (seen[k] == key) { found = 1; break; }
        if (!found) { seen[n++] = key; count++; }
    }
    return count;
}

// Unless OMP_NUM_THREADS says otherwise, use one thread per physical core the
// process may run on, instead of the OpenMP default of one per logical CPU.
// The sparse products and transforms are limited by memory bandwidth, so a
// second hardware thread per core adds little, and every extra thread is one
// more to wait for at each of the ~70 barriers per step; on a busy machine
// that makes surplus threads very costly.
static void default_threads(void)
{
    int cores = physical_cores();
    if (!getenv("OMP_NUM_THREADS") && cores > 0 && cores < omp_get_max_threads())
        omp_set_num_threads(cores);
}
#endif

static void usage(const char *prog)
{
    printf("Usage: %s [options]\n", prog);
    printf("  --n N                grid points per side (the grid is N x N)\n");
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
    int poisson_type  = 3; // 1=Gauss-Seidel  2=SOR  3=FFT (direct, exact)
    int time_scheme   = 2; // 1=Euler  2=RK4
    int use_gpu       = 1; // only meaningful when built with CUDA=1

    // Command-line overrides
    static struct option long_opts[] = {
        {"n",               required_argument, 0, 'n'},
        {"dt",              required_argument, 0, 'd'},
        {"tf",              required_argument, 0, 'f'},
        {"output-interval", required_argument, 0, 'o'},
        {"cpu",             no_argument,       0, 'c'},
        {"help",            no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };
    int opt, opt_index;
    while ((opt = getopt_long(argc, argv, "", long_opts, &opt_index)) != -1)
    {
        int ok = 1;
        switch (opt)
        {
        case 'n': ok = parse_int(optarg, &nx); ny = nx; break;
        case 'd': ok = parse_double(optarg, &dt); break;
        case 'f': ok = parse_double(optarg, &tf); break;
        case 'o': ok = parse_int(optarg, &output_interval); break;
        case 'c': use_gpu = 0; break;
        case 'h': usage(argv[0]); return 0;
        default:  usage(argv[0]); return 1;
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

    // Host memory estimate: 21 arrays of nx*ny doubles (fields, derivatives,
    // RK4 stages, Poisson right-hand side, FFT buffer) and four CSR operators
    // with at most 7 non-zeros per row, compared with the memory available
    // now. A run that would not fit stops here with a message instead of being
    // killed by the kernel once the pages are touched. Other processes can
    // still take memory after this check, so it is a guard against the clear
    // cases only.
    double mem_needed = (double)nx * ny * (21.0 * sizeof(double)
                      + 4.0 * (7.0 * (sizeof(double) + sizeof(int)) + sizeof(int)));
    double mem_avail  = available_memory();
    if (mem_avail >= 0. && mem_needed > mem_avail)
    {
        printf("** Error: a %d x %d grid needs about %.1f GB of memory; about %.1f GB is available **\n",
               nx, ny, mem_needed / 1E9, mem_avail / 1E9);
        return 1;
    }

    // Spectral radius of Jacobi on the (nx-2) x (ny-2) interior nodes
    double rho  = 0.5 * (cos(PI / (nx - 1)) + cos(PI / (ny - 1)));
    double beta  = 2.0 / (1.0 + sqrt(1.0 - rho * rho));  // optimal SOR parameter

    printf("Grid: %d x %d | dt: %lf | tf: %lf\n", nx, ny, dt, tf);
#ifdef _OPENMP
    default_threads();
    printf("OpenMP threads: %d\n", omp_get_max_threads());
#endif
    printf("Poisson SOR parameter: %lf\n", beta);


    // Boundary conditions (Dirichlet)
    double ui = 0., vi = 0.;
    double u1 = 0., u2 = 0., u3 = 0., u4 = 1.;
    double v1 = 0., v2 = 0., v3 = 0., v4 = 0.;
    wall_bc bc = {{u1, u2, u3, u4}, {v1, v2, v3, v4}};

    // Grid spacing: nodes 0 and nx-1 lie on the walls, so nx nodes span Lx
    // with nx-1 intervals
    double dx = (double)Lx / (nx - 1);
    double dy = (double)Ly / (ny - 1);

    // Build sparse 1D operators then free them after Kronecker
    smtrx sd_x  = SDiff1(nx, order, dx);
    smtrx sd_y  = SDiff1(ny, order, dy);
    smtrx sd_x2 = SDiff2(nx, order, dx);
    smtrx sd_y2 = SDiff2(ny, order, dy);
    smtrx sIx   = seye(nx);
    smtrx sIy   = seye(ny);

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
    smtrx DX  = skronecker(sIy,   sd_x);
    smtrx DY  = skronecker(sd_y,  sIx);
    smtrx DX2 = skronecker(sIy,   sd_x2);
    smtrx DY2 = skronecker(sd_y2, sIx);

    freesm(sd_x); freesm(sd_y); freesm(sd_x2); freesm(sd_y2);
    freesm(sIx);  freesm(sIy);

    // Solver workspace (allocated once, reused every timestep)
    rk4_ctx ctx = rk4_alloc(nx, ny);
    ctx.DX = &DX; ctx.DY = &DY; ctx.DX2 = &DX2; ctx.DY2 = &DY2;
    ctx.Re = Re; ctx.dx = dx; ctx.dy = dy;
    ctx.poisson_type = poisson_type;
    ctx.poisson_max_it = poisson_max_it; ctx.poisson_tol = poisson_tol;
    ctx.beta = beta;

    int it_max = (int)((tf / dt) - 1);

    // Dense field matrices
    mtrx u   = initm(nx, ny);
    mtrx v   = initm(nx, ny);
    mtrx w   = initm(nx, ny);

    // Continuity check workspace — pre-allocated once, reused every iteration
    mtrx dudx   = initm(nx, ny);
    mtrx dvdy   = initm(nx, ny);
    mtrx check_continuity = initm(nx, ny);

    // Initial condition
    for (i = 1; i < nx - 1; i++)
        for (j = 1; j < ny - 1; j++)
        {
            MAt(u, i, j) = ui;
            MAt(v, i, j) = vi;
        }

    // Backend selection: GPU when built with CUDA=1 and a device is usable
    int on_gpu = 0;
#ifdef USE_CUDA
    gpu_solver *gpu = NULL;
    if (use_gpu)
    {
        gpu = gpu_init(&ctx, dt, time_scheme, &bc);
        if (gpu)
            gpu_set_fields(gpu, &u, &v, &w);
        else
            printf("No usable CUDA device - falling back to the CPU\n");
    }
    on_gpu = (gpu != NULL);
#else
    (void)use_gpu;
#endif
    printf("Backend: %s\n", on_gpu ? "CUDA" : "CPU");
#ifdef USE_CUDA
    if (on_gpu) printf("CUDA device: %s\n", gpu_device_name());
#endif

    if (!on_gpu && poisson_type == 3) fft_setup(nx, ny);

    struct timespec t_start, t_end;
    clock_gettime(CLOCK_MONOTONIC, &t_start);

    // Main time loop
    for (t = 0; t <= it_max; t++)
    {
        double cmax, cmin;

#ifdef USE_CUDA
        if (on_gpu)
        {
            gpu_step(gpu);
            gpu_continuity(gpu, &cmax, &cmin);
        }
        else
#endif
        {
            // Boundary conditions, time advancement + Poisson solve
            step(w, u, v, dt, time_scheme, &bc, &ctx);

            // Continuity check: du/dx + dv/dy ~ 0
            spmv(DX, u.M, dudx.M);
            spmv(DY, v.M, dvdy.M);

            // reuse check_continuity storage
            for (i = 0; i < nx; i++)
                for (j = 0; j < ny; j++)
                    MAt(check_continuity, i, j) = MAt(dudx, i, j) + MAt(dvdy, i, j);

            cmax = maxel(check_continuity);
            cmin = minel(check_continuity);
        }

        printf("Iteration: %d | Time: %.4lf | Progress: %.2lf%%\n",
               t, (double)t * dt, it_max > 0 ? (double)100 * t / it_max : 100.);
        printf("Continuity max: %E | min: %E\n", cmax, cmin);

        if (output_interval > 0 && t % output_interval == 0)
        {
#ifdef USE_CUDA
            if (on_gpu) gpu_get_fields(gpu, NULL, NULL, &w);
#endif
            printvtk(w, "vorticity");
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t_end);
    double elapsed = (double)(t_end.tv_sec - t_start.tv_sec)
                   + 1E-9 * (double)(t_end.tv_nsec - t_start.tv_nsec);

#ifdef USE_CUDA
    if (on_gpu) gpu_get_fields(gpu, &u, &v, &w);
#endif

    // Re-apply wall BCs before sampling centerline
    apply_wall_bc(u, v, &bc);

    // Write centerline profiles and compare against Ghia et al. (1982)
    print_centerline(u, v, nx, ny, dx, dy);

    // Free dense fields
    freem(&u);
    freem(&v);
    freem(&w);

    // Free continuity workspace
    freem(&dudx);   freem(&dvdy);
    freem(&check_continuity);


    // Free sparse operators
    freesm(DX); freesm(DY); freesm(DX2); freesm(DY2);

    rk4_free(&ctx);
    if (!on_gpu && poisson_type == 3) fft_cleanup();
#ifdef USE_CUDA
    if (on_gpu) gpu_free(gpu);
#endif

    printf("Wall-clock time: %.3lf s total | %.4lf ms per step (%d steps, %s)\n",
           elapsed, 1E3 * elapsed / (it_max + 1), it_max + 1, on_gpu ? "CUDA" : "CPU");
    printf("Simulation complete!\n");
    return 0;
}
