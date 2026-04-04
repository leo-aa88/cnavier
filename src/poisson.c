#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linearalg.h"
#include "poisson.h"

double error(mtrx u1, mtrx u2)
{
    double e = 0;
    int i, j;

    for (i = 0; i < u1.m; i++)
        for (j = 0; j < u1.n; j++)
            e += sqrt(pow(u2.M[i][j] - u1.M[i][j], 2));
    return e;
}

// Gauss-Seidel Poisson solver.
// u and u0 are pre-allocated by the caller; u holds the result on return.
void poisson(mtrx f, mtrx u, mtrx u0, double dx, double dy, int itmax, double tol)
{
    int i, j, k;
    int nx = f.m, ny = f.n;
    double e;
    double dx2 = dx * dx, dy2 = dy * dy;
    double denom = 2.0 * (dx2 + dy2);

    zerosm(u);

    for (k = 0; k < itmax; k++)
    {
        mtrxcpy(u0, u);
        for (i = 1; i < nx - 1; i++)
            for (j = 1; j < ny - 1; j++)
                u.M[i][j] = (dy2 * (u.M[i+1][j] + u.M[i-1][j])
                            + dx2 * (u.M[i][j+1] + u.M[i][j-1])
                            - dx2 * dy2 * f.M[i][j]) / denom;

        e = error(u, u0);
        if (e < tol)
        {
            printf("Poisson solved in %d iterations - RSS error: %E\n", k, e);
            return;
        }
    }
    printf("Error: max iterations reached for Poisson solver.\n");
    exit(1);
}

// SOR Poisson solver.
// u and u0 are pre-allocated by the caller; u holds the result on return.
void poisson_SOR(mtrx f, mtrx u, mtrx u0, double dx, double dy, int itmax, double tol, double beta)
{
    int i, j, k;
    int nx = f.m, ny = f.n;
    double e;
    double dx2 = dx * dx, dy2 = dy * dy;
    double denom = 2.0 * (dx2 + dy2);

    zerosm(u);

    for (k = 0; k < itmax; k++)
    {
        mtrxcpy(u0, u);
        for (i = 1; i < nx - 1; i++)
            for (j = 1; j < ny - 1; j++)
                u.M[i][j] = beta  * (dy2 * (u.M[i+1][j] + u.M[i-1][j])
                                    + dx2 * (u.M[i][j+1] + u.M[i][j-1])
                                    - dx2 * dy2 * f.M[i][j]) / denom
                           + (1.0 - beta) * u0.M[i][j];

        e = error(u, u0);
        if (e < tol)
        {
            printf("Poisson SOR solved in %d iterations - RSS error: %E\n", k, e);
            return;
        }
    }
    printf("Error: max iterations reached for Poisson SOR solver.\n");
    exit(1);
}
