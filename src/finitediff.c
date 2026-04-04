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
        MAt(D, 0, 0) = -1. / dx;
        MAt(D, 0, 1) = 1. / dx;
        for (i = 1; i < (n - 1); i++)
        {
            MAt(D, i, i - 1) = -0.5 / dx;
            MAt(D, i, i) = 0. / dx;
            MAt(D, i, i + 1) = 0.5 / dx;
        }
        MAt(D, n - 1, n - 1) = MAt(D, 0, 1);
        MAt(D, n - 1, n - 2) = MAt(D, 0, 0);
        return D;
    }
    else if (o == 4) // Fourth-order
    {
        MAt(D, 0, 0) = (double)-1 / dx;
        MAt(D, 0, 1) = (double)1 / dx;
        MAt(D, 1, 0) = (double)-0.5 / dx;
        MAt(D, 1, 1) = (double)0 / dx;
        MAt(D, 1, 2) = (double)0.5 / dx;
        for (i = 2; i < (n - 2); i++)
        {
            MAt(D, i, i - 2) =  1.0 / 12.0 / dx;
            MAt(D, i, i - 1) = -2.0 /  3.0 / dx;
            MAt(D, i, i) = 0;
            MAt(D, i, i + 1) =  2.0 /  3.0 / dx;
            MAt(D, i, i + 2) = -1.0 / 12.0 / dx;
        }
        MAt(D, n - 1, n - 1) = MAt(D, 0, 1);
        MAt(D, n - 1, n - 2) = MAt(D, 0, 0);
        MAt(D, n - 2, n - 1) = MAt(D, 1, 2);
        MAt(D, n - 2, n - 2) = MAt(D, 1, 1);
        MAt(D, n - 2, n - 3) = MAt(D, 1, 0);
        return D;
    }
    else if (o == 6) // Sixth-order
    {
        MAt(D, 0, 0) = (double)-1 / dx;
        MAt(D, 0, 1) = (double)1 / dx;
        MAt(D, 1, 0) = (double)-0.5 / dx;
        MAt(D, 1, 1) = (double)0 / dx;
        MAt(D, 1, 2) = (double)0.5 / dx;
        MAt(D, 2, 0) =  1.0 / 12.0 / dx;
        MAt(D, 2, 1) = -2.0 /  3.0 / dx;
        MAt(D, 2, 2) =  0.0;
        MAt(D, 2, 3) =  2.0 /  3.0 / dx;
        MAt(D, 2, 4) = -1.0 / 12.0 / dx;
        for (i = 3; i < (n - 3); i++)
        {
            MAt(D, i, i - 3) = -1.0 / 60.0 / dx;
            MAt(D, i, i - 2) =  3.0 / 20.0 / dx;
            MAt(D, i, i - 1) = -3.0 /  4.0 / dx;
            MAt(D, i, i) = 0.0;
            MAt(D, i, i + 1) =  3.0 /  4.0 / dx;
            MAt(D, i, i + 2) = -3.0 / 20.0 / dx;
            MAt(D, i, i + 3) =  1.0 / 60.0 / dx;
        }
        MAt(D, n - 1, n - 1) = MAt(D, 0, 1);
        MAt(D, n - 1, n - 2) = MAt(D, 0, 0);
        MAt(D, n - 2, n - 1) = MAt(D, 1, 2);
        MAt(D, n - 2, n - 2) = MAt(D, 1, 1);
        MAt(D, n - 2, n - 3) = MAt(D, 1, 0);
        MAt(D, n - 3, n - 1) = MAt(D, 2, 4);
        MAt(D, n - 3, n - 2) = MAt(D, 2, 3);
        MAt(D, n - 3, n - 3) = MAt(D, 2, 2);
        MAt(D, n - 3, n - 4) = MAt(D, 2, 1);
        MAt(D, n - 3, n - 5) = MAt(D, 2, 0);
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
        MAt(D, 0, 0) = 2. / (dx * dx); // Forward scheme (second-order)
        MAt(D, 0, 1) = -5. / (dx * dx);
        MAt(D, 0, 2) = 4. / (dx * dx);
        MAt(D, 0, 3) = -1. / (dx * dx);
        for (i = 1; i < (n - 1); i++)
        {
            MAt(D, i, i - 1) = 1. / (dx * dx);
            MAt(D, i, i) = -2. / (dx * dx);
            MAt(D, i, i + 1) = 1. / (dx * dx);
        }
        MAt(D, n - 1, n - 1) = MAt(D, 0, 0);
        MAt(D, n - 1, n - 2) = MAt(D, 0, 1);
        MAt(D, n - 1, n - 3) = MAt(D, 0, 2);
        MAt(D, n - 1, n - 4) = MAt(D, 0, 3);
        return D;
    }
    else if (o == 4) // Fourth-order
    {
        MAt(D, 0, 0) = (double)2 / (dx * dx); // ForwarD.M scheme (second-order)
        MAt(D, 0, 1) = (double)-5 / (dx * dx);
        MAt(D, 0, 2) = (double)4 / (dx * dx);
        MAt(D, 0, 3) = (double)-1 / (dx * dx);
        MAt(D, 1, 0) = (double)1 / (dx * dx); // Central scheme (second-order)
        MAt(D, 1, 1) = (double)-2 / (dx * dx);
        MAt(D, 1, 2) = (double)1 / (dx * dx);
        for (i = 2; i < (n - 2); i++)
        {
            MAt(D, i, i - 2) =  -1.0 / 12.0 / (dx * dx);
            MAt(D, i, i - 1) =   4.0 /  3.0 / (dx * dx);
            MAt(D, i, i) =      -5.0 /  2.0 / (dx * dx);
            MAt(D, i, i + 1) =   4.0 /  3.0 / (dx * dx);
            MAt(D, i, i + 2) =  -1.0 / 12.0 / (dx * dx);
        }
        MAt(D, n - 1, n - 1) = MAt(D, 0, 0);
        MAt(D, n - 1, n - 2) = MAt(D, 0, 1);
        MAt(D, n - 1, n - 3) = MAt(D, 0, 2);
        MAt(D, n - 1, n - 4) = MAt(D, 0, 3);
        MAt(D, n - 2, n - 1) = MAt(D, 1, 0);
        MAt(D, n - 2, n - 2) = MAt(D, 1, 1);
        MAt(D, n - 2, n - 3) = MAt(D, 1, 2);
        return D;
    }
    else if (o == 6) // Sixth-order
    {
        MAt(D, 0, 0) = (double)2 / (dx * dx); // Forward-scheme (second-order)
        MAt(D, 0, 1) = (double)-5 / (dx * dx);
        MAt(D, 0, 2) = (double)4 / (dx * dx);
        MAt(D, 0, 3) = (double)-1 / (dx * dx);
        MAt(D, 1, 0) = (double)1 / (dx * dx); // Central-scheme (second-order)
        MAt(D, 1, 1) = (double)-2 / (dx * dx);
        MAt(D, 1, 2) = (double)1 / (dx * dx);
        MAt(D, 2, 0) =  -1.0 / 12.0 / (dx * dx); // Central-scheme (fourth-order)
        MAt(D, 2, 1) =   4.0 /  3.0 / (dx * dx);
        MAt(D, 2, 2) =  -5.0 /  2.0 / (dx * dx);
        MAt(D, 2, 3) =   4.0 /  3.0 / (dx * dx);
        MAt(D, 2, 4) =  -1.0 / 12.0 / (dx * dx);
        for (i = 3; i < (n - 3); i++)
        {
            MAt(D, i, i - 3) =   1.0 / 90.0 / (dx * dx);
            MAt(D, i, i - 2) =  -3.0 / 20.0 / (dx * dx);
            MAt(D, i, i - 1) =   3.0 /  2.0 / (dx * dx);
            MAt(D, i, i) =     -49.0 / 18.0 / (dx * dx);
            MAt(D, i, i + 1) =   3.0 /  2.0 / (dx * dx);
            MAt(D, i, i + 2) =  -3.0 / 20.0 / (dx * dx);
            MAt(D, i, i + 3) =   1.0 / 90.0 / (dx * dx);
        }
        MAt(D, n - 1, n - 1) = MAt(D, 0, 0);
        MAt(D, n - 1, n - 2) = MAt(D, 0, 1);
        MAt(D, n - 1, n - 3) = MAt(D, 0, 2);
        MAt(D, n - 1, n - 4) = MAt(D, 0, 3);
        MAt(D, n - 2, n - 1) = MAt(D, 1, 0);
        MAt(D, n - 2, n - 2) = MAt(D, 1, 1);
        MAt(D, n - 2, n - 3) = MAt(D, 1, 2);
        MAt(D, n - 3, n - 1) = MAt(D, 2, 0);
        MAt(D, n - 3, n - 2) = MAt(D, 2, 1);
        MAt(D, n - 3, n - 3) = MAt(D, 2, 2);
        MAt(D, n - 3, n - 4) = MAt(D, 2, 3);
        MAt(D, n - 3, n - 5) = MAt(D, 2, 4);
        return D;
    }
    else
    {
        printf("** Error: valid orders are 2, 4 or 6 **\n");
        exit(1);
    }
}
// Convert a dense finite-difference matrix to CSR sparse format.
// This is called once at startup, so the O(n^2) scan is acceptable.
static smtrx dense_to_csr(mtrx D)
{
    int i, j, nnz = 0;

    // Count non-zeros
    for (i = 0; i < D.m; i++)
        for (j = 0; j < D.n; j++)
            if (MAt(D, i, j) != 0.0)
                nnz++;

    smtrx S = initsm(D.m, D.n, nnz);
    int pos = 0;
    for (i = 0; i < D.m; i++)
    {
        S.row_ptr[i] = pos;
        for (j = 0; j < D.n; j++)
        {
            if (MAt(D, i, j) != 0.0)
            {
                S.values[pos]  = MAt(D, i, j);
                S.col_idx[pos] = j;
                pos++;
            }
        }
    }
    S.row_ptr[D.m] = pos;
    return S;
}

smtrx SDiff1(int n, int o, double dx)
{
    mtrx  D = Diff1(n, o, dx);
    smtrx S = dense_to_csr(D);
    freem(&D);
    return S;
}

smtrx SDiff2(int n, int o, double dx)
{
    mtrx  D = Diff2(n, o, dx);
    smtrx S = dense_to_csr(D);
    freem(&D);
    return S;
}
