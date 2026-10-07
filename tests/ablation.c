// Ablation of the spatial order (issue #25): which part of the discretization
// holds the global order of the manufactured-solution study to two?
//
//   make ablation                    grids 17 ... 129
//   ./ablation_study --max-n 257     one grid further (several minutes)
//
// Each variant replaces parts of the RK4 step with the exact solution at each
// stage's time (see mms_run_ablated() in mms.c):
//   computed           the solver as it is
//   exact wall w       wall vorticity exact; Poisson solve and D_x, D_y as usual
//   exact psi          psi exact; u = D_y psi, v = -D_x psi, wall w computed
//   exact psi, wall w  both
//   exact u, v         velocity exact at every node; wall w computed
//   exact u, v, wall w velocity and wall vorticity exact: what is left is the
//                      transport derivatives of w (with their near-wall rows)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mms.h"

static const struct
{
    const char *name;
    int flags;
} variants[] = {
    {"computed", 0},
    {"exact wall w", MMS_EXACT_WALL_W},
    {"exact psi", MMS_EXACT_PSI},
    {"exact psi, wall w", MMS_EXACT_PSI | MMS_EXACT_WALL_W},
    {"exact u, v", MMS_EXACT_VELOCITY},
    {"exact u, v, wall w", MMS_EXACT_VELOCITY | MMS_EXACT_WALL_W},
};

int main(int argc, char **argv)
{
    int max_n = 129, order, v, n;

    if (argc == 3 && strcmp(argv[1], "--max-n") == 0)
    {
        char *end;
        long m = strtol(argv[2], &end, 10);
        if (*end != '\0' || m < 33 || m > 4097)
        {
            printf("** Error: --max-n must be an integer from 33 to 4097 **\n");
            return 1;
        }
        max_n = (int)m;
    }
    else if (argc != 1)
    {
        printf("Usage: %s [--max-n N]\n", argv[0]);
        return 1;
    }

    printf("Manufactured solution of mms.h on the unit square, Re = 100, t = 0.25, RK4 + FFT.\n"
           "Max error of w over all nodes on each grid, with the observed order against the\n"
           "previous grid; then the observed order of u on the last two grids ('-' when u is\n"
           "exact).\n");
    for (order = 2; order <= 6; order += 2)
    {
        printf("\nDerivative order %d\n%-19s", order, "variant");
        for (n = 17; n <= max_n; n = 2 * n - 1)
            printf(" | w at %-3d  order", n);
        printf(" | u order\n");
        for (v = 0; v < (int)(sizeof(variants) / sizeof(variants[0])); v++)
        {
            double prev_w = 0.0, prev_u = 0.0, order_u = 0.0;

            printf("%-19s", variants[v].name);
            for (n = 17; n <= max_n; n = 2 * n - 1)
            {
                mms_errors e = mms_run_ablated(n, n, 1.0, 1.0, 100.0, order, variants[v].flags, 2.5E-3, 0.0, 0.25);
                if (prev_w > 0.0)
                    printf(" | %.2e  %5.2f", e.w.max, log2(prev_w / e.w.max));
                else
                    printf(" | %.2e       ", e.w.max);
                if (prev_u > 0.0 && e.u.max > 0.0) order_u = log2(prev_u / e.u.max);
                prev_w = e.w.max;
                prev_u = e.u.max;
                fflush(stdout);
            }
            if (variants[v].flags & MMS_EXACT_VELOCITY)
                printf(" |       -\n");
            else
                printf(" | %7.2f\n", order_u);
        }
    }
    return 0;
}
