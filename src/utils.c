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
    fprintf(pf, "LOOKUP_TABLE default\n");

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