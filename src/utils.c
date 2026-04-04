#include <stdlib.h>
#include <stdio.h>
#include "linearalg.h"
#include "utils.h"

double randdouble(double min, double max)
{
    double range = (max - min);
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

void printvtk(mtrx A, char *title)
{
    int i, j;
    char c[320];
    static int count = 0;
    char name[64];
    FILE *pf;

    if (A.M == NULL)
    {
        printf("\n** Error: Aborting program **\n");
        exit(1);
    }
    if ((A.m < 1) || (A.n < 1))
    {
        printf("\n** Error: Invalid parameter **\n");
        exit(1);
    }

    snprintf(name, sizeof(name), "./output/%s-1-%d.vtk", title, count);

    if ((pf = fopen(name, "a")) == NULL)
    {
        printf("\nError while opening file\n");
        exit(1);
    }

    printf("%s\n", name);

    fprintf(pf, "# vtk DataFile Version 2.0\n"); // vtk file headers
    fprintf(pf, "test\n");
    fprintf(pf, "ASCII\n");
    fprintf(pf, "DATASET STRUCTURED_POINTS\n");
    fprintf(pf, "DIMENSIONS %d %d 1\n", A.m, A.n);
    fprintf(pf, "ORIGIN 0 0 0\n");
    fprintf(pf, "SPACING 1 1 1\n");
    fprintf(pf, "POINT_DATA %d\n", A.m * A.n);
    fprintf(pf, "SCALARS values float\n");
    fprintf(pf, "LOOKUP_TABLE default");

    for (i = 0; i < A.m; i++)
    {
        fprintf(pf, "\n");
        for (j = 0; j < A.n; j++)
        {
            if ((j == 0))
            {
                sprintf(c, "%.6lf", MAt(A, i, j));
                fprintf(pf, "%s", c);
            }
            else
            {
                sprintf(c, " %.6lf", MAt(A, i, j));
                fprintf(pf, "%s", c);
            }
        }
    }
    fclose(pf);
    count++;
}

void print_centerline(mtrx u, mtrx v, int nx, int ny, double dx, double dy)
{
    int i, j;
    FILE *f;

    // Ghia et al. (1982), Table 1 — Re=1000
    // u-velocity along vertical centerline x=0.5, y in [0,1]
    static const double ghia_y[]  = {0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719,
                                      0.2813, 0.4531, 0.5000, 0.6172, 0.7344, 0.8516,
                                      0.9531, 0.9609, 0.9688, 0.9766, 1.0000};
    static const double ghia_u[]  = {0.0000,-0.3869,-0.4441,-0.5000,-0.6117,-0.7026,
                                     -0.7014,-0.6504,-0.6174,-0.5155,-0.3732,-0.1877,
                                      0.0570, 0.1867, 0.3322, 0.4664, 1.0000};

    // Ghia et al. (1982), Table 2 — Re=1000
    // v-velocity along horizontal centerline y=0.5, x in [0,1]
    static const double ghia_x[]  = {0.0000, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563,
                                      0.2266, 0.2344, 0.5000, 0.8047, 0.8594, 0.9063,
                                      0.9453, 0.9531, 0.9609, 0.9688, 1.0000};
    static const double ghia_v[]  = {0.0000, 0.1812, 0.2060, 0.2274, 0.2572, 0.3273,
                                      0.3769, 0.3780, 0.0258,-0.4297,-0.5155,-0.5512,
                                     -0.5622,-0.5589,-0.5494,-0.5345, 0.0000};
    int n_ghia = 17;

    // --- u along vertical centerline (i = nx/2) ---
    int ci = nx / 2;
    f = fopen("./output/centerline_u.csv", "w");
    if (!f) { printf("Error opening centerline_u.csv\n"); return; }
    fprintf(f, "y,u_cnavier,y_ghia,u_ghia\n");
    for (j = 0; j < ny; j++)
    {
        double y = (j + 0.5) * dy;
        double u_val = MAt(u, ci, j);
        if (j < n_ghia)
            fprintf(f, "%.6f,%.6f,%.6f,%.6f\n", y, u_val, ghia_y[j], ghia_u[j]);
        else
            fprintf(f, "%.6f,%.6f,,\n", y, u_val);
    }
    fclose(f);

    // --- v along horizontal centerline (j = ny/2) ---
    int cj = ny / 2;
    f = fopen("./output/centerline_v.csv", "w");
    if (!f) { printf("Error opening centerline_v.csv\n"); return; }
    fprintf(f, "x,v_cnavier,x_ghia,v_ghia\n");
    for (i = 0; i < nx; i++)
    {
        double x = (i + 0.5) * dx;
        double v_val = MAt(v, i, cj);
        if (i < n_ghia)
            fprintf(f, "%.6f,%.6f,%.6f,%.6f\n", x, v_val, ghia_x[i], ghia_v[i]);
        else
            fprintf(f, "%.6f,%.6f,,\n", x, v_val);
    }
    fclose(f);

    printf("Centerline profiles written to output/centerline_u.csv and output/centerline_v.csv\n");
}
