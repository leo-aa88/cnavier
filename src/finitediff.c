#include <stdio.h>
#include <stdlib.h>
#include "linearalg.h"
#include "finitediff.h"

// Computes finite-difference matrices for the first derivative
mtrx Diff1(int n, int o, double dx)
{
    int i;
    mtrx D;
    D = initm(n, n);

    if (o == 2) // second order
    {
        D.M[0][0] = -1. / dx;
        D.M[0][1] = 1. / dx;
        for (i = 1; i < (n - 1); i++)
        {
            D.M[i][i - 1] = -0.5 / dx;
            D.M[i][i] = 0. / dx;
            D.M[i][i + 1] = 0.5 / dx;
        }
        D.M[n - 1][n - 1] = D.M[0][1];
        D.M[n - 1][n - 2] = D.M[0][0];
        return D;
    }
    else if (o == 4) // Fourth-order
    {
        D.M[0][0] = (double)-1 / dx;
        D.M[0][1] = (double)1 / dx;
        D.M[1][0] = (double)-0.5 / dx;
        D.M[1][1] = (double)0 / dx;
        D.M[1][2] = (double)0.5 / dx;
        for (i = 2; i < (n - 2); i++)
        {
            D.M[i][i - 2] =  1.0 / 12.0 / dx;
            D.M[i][i - 1] = -2.0 /  3.0 / dx;
            D.M[i][i] = 0;
            D.M[i][i + 1] =  2.0 /  3.0 / dx;
            D.M[i][i + 2] = -1.0 / 12.0 / dx;
        }
        D.M[n - 1][n - 1] = D.M[0][1];
        D.M[n - 1][n - 2] = D.M[0][0];
        D.M[n - 2][n - 1] = D.M[1][2];
        D.M[n - 2][n - 2] = D.M[1][1];
        D.M[n - 2][n - 3] = D.M[1][0];
        return D;
    }
    else if (o == 6) // Sixth-order
    {
        D.M[0][0] = (double)-1 / dx;
        D.M[0][1] = (double)1 / dx;
        D.M[1][0] = (double)-0.5 / dx;
        D.M[1][1] = (double)0 / dx;
        D.M[1][2] = (double)0.5 / dx;
        D.M[2][0] =  1.0 / 12.0 / dx;
        D.M[2][1] = -2.0 /  3.0 / dx;
        D.M[2][2] =  0.0;
        D.M[2][3] =  2.0 /  3.0 / dx;
        D.M[2][4] = -1.0 / 12.0 / dx;
        for (i = 3; i < (n - 3); i++)
        {
            D.M[i][i - 3] = -1.0 / 60.0 / dx;
            D.M[i][i - 2] =  3.0 / 20.0 / dx;
            D.M[i][i - 1] = -3.0 /  4.0 / dx;
            D.M[i][i] = 0.0;
            D.M[i][i + 1] =  3.0 /  4.0 / dx;
            D.M[i][i + 2] = -3.0 / 20.0 / dx;
            D.M[i][i + 3] =  1.0 / 60.0 / dx;
        }
        D.M[n - 1][n - 1] = D.M[0][1];
        D.M[n - 1][n - 2] = D.M[0][0];
        D.M[n - 2][n - 1] = D.M[1][2];
        D.M[n - 2][n - 2] = D.M[1][1];
        D.M[n - 2][n - 3] = D.M[1][0];
        D.M[n - 3][n - 1] = D.M[2][4];
        D.M[n - 3][n - 2] = D.M[2][3];
        D.M[n - 3][n - 3] = D.M[2][2];
        D.M[n - 3][n - 4] = D.M[2][1];
        D.M[n - 3][n - 5] = D.M[2][0];
        return D;
    }
    else
    {
        printf("** Error: valid orders are 2, 4 or 6 **\n");
        exit(1);
    }
}

// Computes finite - difference matrices for the second derivative
mtrx Diff2(int n, int o, double dx)
{
    int i;
    mtrx D = initm(n, n);

    if (o == 2) // Second-order
    {
        D.M[0][0] = 2. / (dx * dx); // Forward scheme (second-order)
        D.M[0][1] = -5. / (dx * dx);
        D.M[0][2] = 4. / (dx * dx);
        D.M[0][3] = -1. / (dx * dx);
        for (i = 1; i < (n - 1); i++)
        {
            D.M[i][i - 1] = 1. / (dx * dx);
            D.M[i][i] = -2. / (dx * dx);
            D.M[i][i + 1] = 1. / (dx * dx);
        }
        D.M[n - 1][n - 1] = D.M[0][0];
        D.M[n - 1][n - 2] = D.M[0][1];
        D.M[n - 1][n - 3] = D.M[0][2];
        D.M[n - 1][n - 4] = D.M[0][3];
        return D;
    }
    else if (o == 4) // Fourth-order
    {
        D.M[0][0] = (double)2 / (dx * dx); // ForwarD.M scheme (second-order)
        D.M[0][1] = (double)-5 / (dx * dx);
        D.M[0][2] = (double)4 / (dx * dx);
        D.M[0][3] = (double)-1 / (dx * dx);
        D.M[1][0] = (double)1 / (dx * dx); // Central scheme (second-order)
        D.M[1][1] = (double)-2 / (dx * dx);
        D.M[1][2] = (double)1 / (dx * dx);
        for (i = 2; i < (n - 2); i++)
        {
            D.M[i][i - 2] =  -1.0 / 12.0 / (dx * dx);
            D.M[i][i - 1] =   4.0 /  3.0 / (dx * dx);
            D.M[i][i] =      -5.0 /  2.0 / (dx * dx);
            D.M[i][i + 1] =   4.0 /  3.0 / (dx * dx);
            D.M[i][i + 2] =  -1.0 / 12.0 / (dx * dx);
        }
        D.M[n - 1][n - 1] = D.M[0][0];
        D.M[n - 1][n - 2] = D.M[0][1];
        D.M[n - 1][n - 3] = D.M[0][2];
        D.M[n - 1][n - 4] = D.M[0][3];
        D.M[n - 2][n - 1] = D.M[1][0];
        D.M[n - 2][n - 2] = D.M[1][1];
        D.M[n - 2][n - 3] = D.M[1][2];
        return D;
    }
    else if (o == 6) // Sixth-order
    {
        D.M[0][0] = (double)2 / (dx * dx); // Forward-scheme (second-order)
        D.M[0][1] = (double)-5 / (dx * dx);
        D.M[0][2] = (double)4 / (dx * dx);
        D.M[0][3] = (double)-1 / (dx * dx);
        D.M[1][0] = (double)1 / (dx * dx); // Central-scheme (second-order)
        D.M[1][1] = (double)-2 / (dx * dx);
        D.M[1][2] = (double)1 / (dx * dx);
        D.M[2][0] =  -1.0 / 12.0 / (dx * dx); // Central-scheme (fourth-order)
        D.M[2][1] =   4.0 /  3.0 / (dx * dx);
        D.M[2][2] =  -5.0 /  2.0 / (dx * dx);
        D.M[2][3] =   4.0 /  3.0 / (dx * dx);
        D.M[2][4] =  -1.0 / 12.0 / (dx * dx);
        for (i = 3; i < (n - 3); i++)
        {
            D.M[i][i - 3] =   1.0 / 90.0 / (dx * dx);
            D.M[i][i - 2] =  -3.0 / 20.0 / (dx * dx);
            D.M[i][i - 1] =   3.0 /  2.0 / (dx * dx);
            D.M[i][i] =     -49.0 / 18.0 / (dx * dx);
            D.M[i][i + 1] =   3.0 /  2.0 / (dx * dx);
            D.M[i][i + 2] =  -3.0 / 20.0 / (dx * dx);
            D.M[i][i + 3] =   1.0 / 90.0 / (dx * dx);
        }
        D.M[n - 1][n - 1] = D.M[0][0];
        D.M[n - 1][n - 2] = D.M[0][1];
        D.M[n - 1][n - 3] = D.M[0][2];
        D.M[n - 1][n - 4] = D.M[0][3];
        D.M[n - 2][n - 1] = D.M[1][0];
        D.M[n - 2][n - 2] = D.M[1][1];
        D.M[n - 2][n - 3] = D.M[1][2];
        D.M[n - 3][n - 1] = D.M[2][0];
        D.M[n - 3][n - 2] = D.M[2][1];
        D.M[n - 3][n - 3] = D.M[2][2];
        D.M[n - 3][n - 4] = D.M[2][3];
        D.M[n - 3][n - 5] = D.M[2][4];
        return D;
    }
    else
    {
        printf("** Error: valid orders are 2, 4 or 6 **\n");
        exit(1);
    }
}
