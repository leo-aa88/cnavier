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
    int *row, *col;
    double *val;
    int count, cap;
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
    int *order = (int *)calloc(D->count > 0 ? D->count : 1, sizeof(int));
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
                S.values[pos] = D->val[k];
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
            op_set(D, i, i - 2, 1.0 / 12.0 / dx);
            op_set(D, i, i - 1, -2.0 / 3.0 / dx);
            op_set(D, i, i, 0);
            op_set(D, i, i + 1, 2.0 / 3.0 / dx);
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
        op_set(D, 2, 0, 1.0 / 12.0 / dx);
        op_set(D, 2, 1, -2.0 / 3.0 / dx);
        op_set(D, 2, 2, 0.0);
        op_set(D, 2, 3, 2.0 / 3.0 / dx);
        op_set(D, 2, 4, -1.0 / 12.0 / dx);
        for (i = 3; i < (n - 3); i++)
        {
            op_set(D, i, i - 3, -1.0 / 60.0 / dx);
            op_set(D, i, i - 2, 3.0 / 20.0 / dx);
            op_set(D, i, i - 1, -3.0 / 4.0 / dx);
            op_set(D, i, i, 0.0);
            op_set(D, i, i + 1, 3.0 / 4.0 / dx);
            op_set(D, i, i + 2, -3.0 / 20.0 / dx);
            op_set(D, i, i + 3, 1.0 / 60.0 / dx);
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
            op_set(D, i, i - 2, -1.0 / 12.0 / (dx * dx));
            op_set(D, i, i - 1, 4.0 / 3.0 / (dx * dx));
            op_set(D, i, i, -5.0 / 2.0 / (dx * dx));
            op_set(D, i, i + 1, 4.0 / 3.0 / (dx * dx));
            op_set(D, i, i + 2, -1.0 / 12.0 / (dx * dx));
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
        op_set(D, 2, 0, -1.0 / 12.0 / (dx * dx)); // Central-scheme (fourth-order)
        op_set(D, 2, 1, 4.0 / 3.0 / (dx * dx));
        op_set(D, 2, 2, -5.0 / 2.0 / (dx * dx));
        op_set(D, 2, 3, 4.0 / 3.0 / (dx * dx));
        op_set(D, 2, 4, -1.0 / 12.0 / (dx * dx));
        for (i = 3; i < (n - 3); i++)
        {
            op_set(D, i, i - 3, 1.0 / 90.0 / (dx * dx));
            op_set(D, i, i - 2, -3.0 / 20.0 / (dx * dx));
            op_set(D, i, i - 1, 3.0 / 2.0 / (dx * dx));
            op_set(D, i, i, -49.0 / 18.0 / (dx * dx));
            op_set(D, i, i + 1, 3.0 / 2.0 / (dx * dx));
            op_set(D, i, i + 2, -3.0 / 20.0 / (dx * dx));
            op_set(D, i, i + 3, 1.0 / 90.0 / (dx * dx));
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

// Centered stencils, coefficients of f[i-h] ... f[i+h] (before dividing by
// dx or dx^2). They are the interior rows of SDiff1/SDiff2.
static const double d1_coef[3][7] = {
    {-1. / 2., 0., 1. / 2.},
    {1. / 12., -2. / 3., 0., 2. / 3., -1. / 12.},
    {-1. / 60., 3. / 20., -3. / 4., 0., 3. / 4., -3. / 20., 1. / 60.},
};
static const double d2_coef[3][7] = {
    {1., -2., 1.},
    {-1. / 12., 4. / 3., -5. / 2., 4. / 3., -1. / 12.},
    {1. / 90., -3. / 20., 3. / 2., -49. / 18., 3. / 2., -3. / 20., 1. / 90.},
};

// The centered stencil of order o on every row, wrapping around the ends
static smtrx periodic_op(int n, int o, double scale, const double coef[3][7])
{
    op_builder D = {0};
    int i, m, h = o / 2;

    if (o != 2 && o != 4 && o != 6)
    {
        printf("** Error: valid orders are 2, 4 or 6 **\n");
        exit(1);
    }
    if (n < 2 * h + 1)
    {
        printf("** Error: a periodic order-%d stencil needs at least %d points **\n", o, 2 * h + 1);
        exit(1);
    }
    for (i = 0; i < n; i++)
        for (m = -h; m <= h; m++)
            if (coef[h - 1][m + h] != 0.0)
                op_set(&D, i, ((i + m) % n + n) % n, coef[h - 1][m + h] / scale);
    return op_to_csr(&D, n);
}

smtrx SDiff1_periodic(int n, int o, double dx)
{
    return periodic_op(n, o, dx, d1_coef);
}

smtrx SDiff2_periodic(int n, int o, double dx)
{
    return periodic_op(n, o, dx * dx, d2_coef);
}

// SDiff1 with fourth-order one-sided rows at the ends and next to them (orders
// 4 and 6; order 2 is returned unchanged). The deeper rows are SDiff1's: the
// next one is already fourth-order centered.
//   row 0:   (-25 f0 + 48 f1 - 36 f2 + 16 f3 - 3 f4) / (12 dx)
//   row 1:   (-3 f0 - 10 f1 + 18 f2 - 6 f3 + f4) / (12 dx)
// and the mirror images, with the opposite sign, at the other end
smtrx SDiff1_wall4(int n, int o, double dx)
{
    static const double row0[5] = {-25., 48., -36., 16., -3.}, row1[5] = {-3., -10., 18., -6., 1.};
    op_builder D = {0};
    int k;

    build_diff1(&D, n, o, dx);
    if (o >= 4)
    {
        if (n < 10)
        {
            printf("** Error: the fourth-order wall rows need at least 10 points **\n");
            exit(1);
        }
        // A later write to the same position replaces the earlier one, so
        // first zero every entry of the four rows, then write the new ones
        for (k = 0; k < 7; k++)
        {
            op_set(&D, 0, k, 0.0);
            op_set(&D, 1, k, 0.0);
            op_set(&D, n - 1, n - 1 - k, 0.0);
            op_set(&D, n - 2, n - 1 - k, 0.0);
        }
        for (k = 0; k < 5; k++)
        {
            op_set(&D, 0, k, row0[k] / (12.0 * dx));
            op_set(&D, 1, k, row1[k] / (12.0 * dx));
            op_set(&D, n - 1, n - 1 - k, -row0[k] / (12.0 * dx));
            op_set(&D, n - 2, n - 1 - k, -row1[k] / (12.0 * dx));
        }
    }
    return op_to_csr(&D, n);
}
