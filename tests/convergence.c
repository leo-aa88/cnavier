// Spatial convergence study against the manufactured solutions in mms.h, with
// walls and on a doubly periodic grid.
//
//   make convergence                 grids 17 ... 257
//   ./convergence_study --max-n 513  up to 513 (several minutes)
//
// Every run uses RK4 and the FFT Poisson solver, starts from the exact state
// at t = 0 and is compared with the exact solution at t = 0.25. The time step
// is small enough that the time error is negligible: halving it changes no
// digit printed. The observed order between two grids is
// log2(e_coarse / e_fine), since each grid halves the spacing.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mms.h"

#define T_FINAL 0.25
#define RE      100.0
#define DT      2.5E-3 // capped by mms_run() at half the stability limit

// One error norm and its observed order against the previous grid
static void cell(double e, double prev)
{
    if (prev > 0.0 && e > 0.0)
        printf(" | %.2e  %5.2f", e, log2(prev / e));
    else
        printf(" | %.2e       ", e);
}

// Refine an (n-1)*sx x (n-1)*sy cell grid on an Lx x Ly domain, from n = 17
// nodes along the shorter side up to max_n nodes along the longer one
static void study(const char *title, int order, double Lx, double Ly, int sx, int sy, int max_n,
                  int poisson_order, int wall_closure, int velocity_order)
{
    int n, s = sx > sy ? sx : sy;
    mms_errors prev;

    printf("\n%s\n", title);
    printf("%-9s | %-15s | %-15s | %-15s | %-15s | %-15s | %-15s\n", "grid", "psi max   order",
           "u max     order", "v max     order", "w max     order", "w int.    order", "w rms     order");
    memset(&prev, 0, sizeof(prev));
    for (n = 17; (n - 1) * s + 1 <= max_n; n = 2 * n - 1)
    {
        int nx = (n - 1) * sx + 1, ny = (n - 1) * sy + 1;
        mms_errors e = mms_run_closures(nx, ny, Lx, Ly, RE, order, poisson_order, wall_closure, velocity_order, DT, 0.0, T_FINAL);
        char grid[32];

        snprintf(grid, sizeof(grid), "%dx%d", nx, ny);
        printf("%-9s", grid);
        cell(e.psi.max, prev.psi.max);
        cell(e.u.max, prev.u.max);
        cell(e.v.max, prev.v.max);
        cell(e.w.max, prev.w.max);
        cell(e.w_interior.max, prev.w_interior.max);
        cell(e.w.rms, prev.w.rms);
        printf("\n");
        fflush(stdout);
        prev = e;
    }
}

// The periodic solution on n x n grids, n = 16, 32, ..., up to max_n
static void study_periodic(int order, int max_n)
{
    int n;
    mms_errors prev;

    printf("\nDoubly periodic unit square, derivative order %d\n", order);
    printf("%-9s | %-15s | %-15s | %-15s | %-15s | %-15s\n", "grid", "psi max   order", "u max     order",
           "v max     order", "w max     order", "w rms     order");
    memset(&prev, 0, sizeof(prev));
    for (n = 16; n <= max_n; n *= 2)
    {
        mms_errors e = mms_run_periodic(n, n, 1.0, 1.0, RE, order, 2, DT, 0.0, T_FINAL);
        char grid[32];

        snprintf(grid, sizeof(grid), "%dx%d", n, n);
        printf("%-9s", grid);
        cell(e.psi.max, prev.psi.max);
        cell(e.u.max, prev.u.max);
        cell(e.v.max, prev.v.max);
        cell(e.w.max, prev.w.max);
        cell(e.w.rms, prev.w.rms);
        printf("\n");
        fflush(stdout);
        prev = e;
    }
}

int main(int argc, char **argv)
{
    int max_n = 257, order;

    if (argc == 3 && strcmp(argv[1], "--max-n") == 0)
    {
        char *end;
        long n = strtol(argv[2], &end, 10);
        if (*end != '\0' || n < 17 || n > 4097)
        {
            printf("** Error: --max-n must be an integer from 17 to 4097 **\n");
            return 1;
        }
        max_n = (int)n;
    }
    else if (argc != 1)
    {
        printf("Usage: %s [--max-n N]\n", argv[0]);
        return 1;
    }

    printf("Manufactured solution psi = g(t) sin^2(pi x/Lx) sin^2(pi y/Ly), Re = %g, t = %g,\n"
           "RK4 + FFT Poisson solver. Errors are max |numerical - exact| over the nodes\n"
           "(w max includes the wall nodes; w interior excludes them), with the observed\n"
           "order against the previous grid.\n",
           RE, T_FINAL);
    for (order = 2; order <= 6; order += 2)
    {
        char title[96];
        snprintf(title, sizeof(title), "Unit square, derivative order %d", order);
        study(title, order, 1.0, 1.0, 1, 1, max_n, 2, 0, 2);
    }
    study("2 x 1 domain, nx - 1 = 2 (ny - 1), derivative order 6", 6, 2.0, 1.0, 2, 1, max_n, 2, 0, 2);
    study("Unit square, nx - 1 = 2 (ny - 1) (dx = dy / 2), derivative order 6", 6, 1.0, 1.0, 2, 1, max_n, 2, 0, 2);

    printf("\nHigher order with walls (issues #26, #29): the compact Poisson operator\n"
           "(--poisson-order 4) and the third-order wall closure from psi\n"
           "(--wall-closure briley), alone and together, and with fourth-order velocity\n"
           "rows next to the walls (--velocity-order 4).\n");
    study("Unit square, order 6, compact Poisson operator", 6, 1.0, 1.0, 1, 1, max_n, 4, 0, 2);
    study("Unit square, order 6, wall closure from psi", 6, 1.0, 1.0, 1, 1, max_n, 2, 1, 2);
    study("Unit square, order 6, compact Poisson operator and wall closure from psi", 6, 1.0, 1.0, 1, 1, max_n, 4, 1, 2);
    study("Unit square, order 6, both, and fourth-order velocity rows next to the walls", 6, 1.0, 1.0, 1, 1, max_n, 4,
          1, 4);
    study("2 x 1 domain, nx - 1 = 2 (ny - 1), order 6, all three", 6, 2.0, 1.0, 2, 1, max_n, 4, 1, 4);

    printf("\nWithout walls: the periodic solution, a sum of three Fourier modes, on a\n"
           "doubly periodic grid; the Poisson operator is DX2 + DY2, of the same order\n"
           "as the derivatives.\n");
    for (order = 2; order <= 6; order += 2)
        study_periodic(order, max_n);
    return 0;
}
