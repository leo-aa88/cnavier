#include <stdio.h>
#include <stdlib.h>
#include "linearalg.h"
#include "finitediff.h"

// The operators are assembled as a list of (row, column, value) entries and
// converted to CSR, so building one costs O(n) memory instead of a dense n x n
// matrix. Entries are written in the same way the dense matrix used to be:
// a later write to the same position replaces the earlier one, and zeros are
// not stored.
typedef struct
{
    int    *row, *col;
    double *val;
    int     count, cap;
} op_builder;

static void op_set(op_builder *D, int i, int j, double v)
{
    if (D->count == D->cap)
    {
        D->cap = D->cap ? 2 * D->cap : 64;
        D->row = (int *)realloc(D->row, D->cap * sizeof(int));
        D->col = (int *)realloc(D->col, D->cap * sizeof(int));
        D->val = (double *)realloc(D->val, D->cap * sizeof(double));
        if (!D->row || !D->col || !D->val)
        {
            printf("** Error: insufficient memory for finite-difference operator **\n");
            exit(1);
        }
    }
    D->row[D->count] = i;
    D->col[D->count] = j;
    D->val[D->count] = v;
    D->count++;
}

// Value most recently written at (i, j), or 0. Only used for the few boundary
// rows that are copied from the first rows.
static double op_get(const op_builder *D, int i, int j)
{
    int k;
    for (k = D->count - 1; k >= 0; k--)
        if (D->row[k] == i && D->col[k] == j)
            return D->val[k];
    return 0.0;
}

// Convert the entry list to CSR with columns in ascending order
static smtrx op_to_csr(op_builder *D, int n)
{
    int i, k, a, b, pos;
    int *count = (int *)calloc(n + 1, sizeof(int));
    int *order = (int *)malloc((D->count > 0 ? D->count : 1) * sizeof(int));
    if (!count || !order)
    {
        printf("** Error: insufficient memory for finite-difference operator **\n");
        exit(1);
    }

    // Stable bucket sort of the entries by row
    for (k = 0; k < D->count; k++)
        count[D->row[k] + 1]++;
    for (i = 0; i < n; i++)
        count[i + 1] += count[i];
    for (k = 0; k < D->count; k++)
        order[count[D->row[k]]++] = k;
    for (i = n; i > 0; i--)
        count[i] = count[i - 1];
    count[0] = 0;

    // Within each row: stable insertion sort by column, so that of two writes
    // to the same position the later one ends up last
    for (i = 0; i < n; i++)
        for (a = count[i] + 1; a < count[i + 1]; a++)
        {
            int key = order[a];
            for (b = a - 1; b >= count[i] && D->col[order[b]] > D->col[key]; b--)
                order[b + 1] = order[b];
            order[b + 1] = key;
        }

    smtrx S = initsm(n, n, D->count > 0 ? D->count : 1);
    pos = 0;
    for (i = 0; i < n; i++)
    {
        S.row_ptr[i] = pos;
        for (a = count[i]; a < count[i + 1]; a++)
        {
            k = order[a];
            // Skip a write that a later one to the same position replaces
            if (a + 1 < count[i + 1] && D->col[order[a + 1]] == D->col[k])
                continue;
            if (D->val[k] != 0.0)
            {
                S.values[pos]  = D->val[k];
                S.col_idx[pos] = D->col[k];
                pos++;
            }
        }
    }
    S.row_ptr[n] = pos;
    S.nnz = pos;

    free(count);
    free(order);
    free(D->row);
    free(D->col);
    free(D->val);
    return S;
}

static mtrx csr_to_dense(smtrx S)
{
    int i, k;
    mtrx D = initm(S.m, S.n);
    for (i = 0; i < S.m; i++)
        for (k = S.row_ptr[i]; k < S.row_ptr[i + 1]; k++)
            MAt(D, i, S.col_idx[k]) = S.values[k];
    return D;
}

// Computes finite-difference matrices for the first derivative
static void build_diff1(op_builder *D, int n, int o, double dx)
{
    int i;

    if (o == 2) // second order
    {
        op_set(D, 0, 0, -1. / dx);
        op_set(D, 0, 1, 1. / dx);
        for (i = 1; i < (n - 1); i++)
        {
            op_set(D, i, i - 1, -0.5 / dx);
            op_set(D, i, i, 0. / dx);
            op_set(D, i, i + 1, 0.5 / dx);
        }
        op_set(D, n - 1, n - 1, op_get(D, 0, 1));
        op_set(D, n - 1, n - 2, op_get(D, 0, 0));
        return;
    }
    else if (o == 4) // Fourth-order
    {
        op_set(D, 0, 0, (double)-1 / dx);
        op_set(D, 0, 1, (double)1 / dx);
        op_set(D, 1, 0, (double)-0.5 / dx);
        op_set(D, 1, 1, (double)0 / dx);
        op_set(D, 1, 2, (double)0.5 / dx);
        for (i = 2; i < (n - 2); i++)
        {
            op_set(D, i, i - 2,  1.0 / 12.0 / dx);
            op_set(D, i, i - 1, -2.0 /  3.0 / dx);
            op_set(D, i, i, 0);
            op_set(D, i, i + 1,  2.0 /  3.0 / dx);
            op_set(D, i, i + 2, -1.0 / 12.0 / dx);
        }
        op_set(D, n - 1, n - 1, op_get(D, 0, 1));
        op_set(D, n - 1, n - 2, op_get(D, 0, 0));
        op_set(D, n - 2, n - 1, op_get(D, 1, 2));
        op_set(D, n - 2, n - 2, op_get(D, 1, 1));
        op_set(D, n - 2, n - 3, op_get(D, 1, 0));
        return;
    }
    else if (o == 6) // Sixth-order
    {
        op_set(D, 0, 0, (double)-1 / dx);
        op_set(D, 0, 1, (double)1 / dx);
        op_set(D, 1, 0, (double)-0.5 / dx);
        op_set(D, 1, 1, (double)0 / dx);
        op_set(D, 1, 2, (double)0.5 / dx);
        op_set(D, 2, 0,  1.0 / 12.0 / dx);
        op_set(D, 2, 1, -2.0 /  3.0 / dx);
        op_set(D, 2, 2,  0.0);
        op_set(D, 2, 3,  2.0 /  3.0 / dx);
        op_set(D, 2, 4, -1.0 / 12.0 / dx);
        for (i = 3; i < (n - 3); i++)
        {
            op_set(D, i, i - 3, -1.0 / 60.0 / dx);
            op_set(D, i, i - 2,  3.0 / 20.0 / dx);
            op_set(D, i, i - 1, -3.0 /  4.0 / dx);
            op_set(D, i, i, 0.0);
            op_set(D, i, i + 1,  3.0 /  4.0 / dx);
            op_set(D, i, i + 2, -3.0 / 20.0 / dx);
            op_set(D, i, i + 3,  1.0 / 60.0 / dx);
        }
        op_set(D, n - 1, n - 1, op_get(D, 0, 1));
        op_set(D, n - 1, n - 2, op_get(D, 0, 0));
        op_set(D, n - 2, n - 1, op_get(D, 1, 2));
        op_set(D, n - 2, n - 2, op_get(D, 1, 1));
        op_set(D, n - 2, n - 3, op_get(D, 1, 0));
        op_set(D, n - 3, n - 1, op_get(D, 2, 4));
        op_set(D, n - 3, n - 2, op_get(D, 2, 3));
        op_set(D, n - 3, n - 3, op_get(D, 2, 2));
        op_set(D, n - 3, n - 4, op_get(D, 2, 1));
        op_set(D, n - 3, n - 5, op_get(D, 2, 0));
        return;
    }
    else
    {
        printf("** Error: valid orders are 2, 4 or 6 **\n");
        exit(1);
    }
}

// Computes finite - difference matrices for the second derivative
static void build_diff2(op_builder *D, int n, int o, double dx)
{
    int i;

    if (o == 2) // Second-order
    {
        op_set(D, 0, 0, 2. / (dx * dx)); // Forward scheme (second-order)
        op_set(D, 0, 1, -5. / (dx * dx));
        op_set(D, 0, 2, 4. / (dx * dx));
        op_set(D, 0, 3, -1. / (dx * dx));
        for (i = 1; i < (n - 1); i++)
        {
            op_set(D, i, i - 1, 1. / (dx * dx));
            op_set(D, i, i, -2. / (dx * dx));
            op_set(D, i, i + 1, 1. / (dx * dx));
        }
        op_set(D, n - 1, n - 1, op_get(D, 0, 0));
        op_set(D, n - 1, n - 2, op_get(D, 0, 1));
        op_set(D, n - 1, n - 3, op_get(D, 0, 2));
        op_set(D, n - 1, n - 4, op_get(D, 0, 3));
        return;
    }
    else if (o == 4) // Fourth-order
    {
        op_set(D, 0, 0, (double)2 / (dx * dx)); // ForwarD.M scheme (second-order)
        op_set(D, 0, 1, (double)-5 / (dx * dx));
        op_set(D, 0, 2, (double)4 / (dx * dx));
        op_set(D, 0, 3, (double)-1 / (dx * dx));
        op_set(D, 1, 0, (double)1 / (dx * dx)); // Central scheme (second-order)
        op_set(D, 1, 1, (double)-2 / (dx * dx));
        op_set(D, 1, 2, (double)1 / (dx * dx));
        for (i = 2; i < (n - 2); i++)
        {
            op_set(D, i, i - 2,  -1.0 / 12.0 / (dx * dx));
            op_set(D, i, i - 1,   4.0 /  3.0 / (dx * dx));
            op_set(D, i, i,      -5.0 /  2.0 / (dx * dx));
            op_set(D, i, i + 1,   4.0 /  3.0 / (dx * dx));
            op_set(D, i, i + 2,  -1.0 / 12.0 / (dx * dx));
        }
        op_set(D, n - 1, n - 1, op_get(D, 0, 0));
        op_set(D, n - 1, n - 2, op_get(D, 0, 1));
        op_set(D, n - 1, n - 3, op_get(D, 0, 2));
        op_set(D, n - 1, n - 4, op_get(D, 0, 3));
        op_set(D, n - 2, n - 1, op_get(D, 1, 0));
        op_set(D, n - 2, n - 2, op_get(D, 1, 1));
        op_set(D, n - 2, n - 3, op_get(D, 1, 2));
        return;
    }
    else if (o == 6) // Sixth-order
    {
        op_set(D, 0, 0, (double)2 / (dx * dx)); // Forward-scheme (second-order)
        op_set(D, 0, 1, (double)-5 / (dx * dx));
        op_set(D, 0, 2, (double)4 / (dx * dx));
        op_set(D, 0, 3, (double)-1 / (dx * dx));
        op_set(D, 1, 0, (double)1 / (dx * dx)); // Central-scheme (second-order)
        op_set(D, 1, 1, (double)-2 / (dx * dx));
        op_set(D, 1, 2, (double)1 / (dx * dx));
        op_set(D, 2, 0,  -1.0 / 12.0 / (dx * dx)); // Central-scheme (fourth-order)
        op_set(D, 2, 1,   4.0 /  3.0 / (dx * dx));
        op_set(D, 2, 2,  -5.0 /  2.0 / (dx * dx));
        op_set(D, 2, 3,   4.0 /  3.0 / (dx * dx));
        op_set(D, 2, 4,  -1.0 / 12.0 / (dx * dx));
        for (i = 3; i < (n - 3); i++)
        {
            op_set(D, i, i - 3,   1.0 / 90.0 / (dx * dx));
            op_set(D, i, i - 2,  -3.0 / 20.0 / (dx * dx));
            op_set(D, i, i - 1,   3.0 /  2.0 / (dx * dx));
            op_set(D, i, i,     -49.0 / 18.0 / (dx * dx));
            op_set(D, i, i + 1,   3.0 /  2.0 / (dx * dx));
            op_set(D, i, i + 2,  -3.0 / 20.0 / (dx * dx));
            op_set(D, i, i + 3,   1.0 / 90.0 / (dx * dx));
        }
        op_set(D, n - 1, n - 1, op_get(D, 0, 0));
        op_set(D, n - 1, n - 2, op_get(D, 0, 1));
        op_set(D, n - 1, n - 3, op_get(D, 0, 2));
        op_set(D, n - 1, n - 4, op_get(D, 0, 3));
        op_set(D, n - 2, n - 1, op_get(D, 1, 0));
        op_set(D, n - 2, n - 2, op_get(D, 1, 1));
        op_set(D, n - 2, n - 3, op_get(D, 1, 2));
        op_set(D, n - 3, n - 1, op_get(D, 2, 0));
        op_set(D, n - 3, n - 2, op_get(D, 2, 1));
        op_set(D, n - 3, n - 3, op_get(D, 2, 2));
        op_set(D, n - 3, n - 4, op_get(D, 2, 3));
        op_set(D, n - 3, n - 5, op_get(D, 2, 4));
        return;
    }
    else
    {
        printf("** Error: valid orders are 2, 4 or 6 **\n");
        exit(1);
    }
}
smtrx SDiff1(int n, int o, double dx)
{
    op_builder D = {0};
    build_diff1(&D, n, o, dx);
    return op_to_csr(&D, n);
}

smtrx SDiff2(int n, int o, double dx)
{
    op_builder D = {0};
    build_diff2(&D, n, o, dx);
    return op_to_csr(&D, n);
}

mtrx Diff1(int n, int o, double dx)
{
    smtrx S = SDiff1(n, o, dx);
    mtrx  D = csr_to_dense(S);
    freesm(S);
    return D;
}

mtrx Diff2(int n, int o, double dx)
{
    smtrx S = SDiff2(n, o, dx);
    mtrx  D = csr_to_dense(S);
    freesm(S);
    return D;
}
