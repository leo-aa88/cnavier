#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "linearalg.h"

void zerosm(mtrx A)
{
    memset(A.M, 0, (size_t)A.m * A.n * sizeof(double));
}

// Single contiguous allocation: no pointer-of-pointers, cache-friendly
double *allocm(int m, int n)
{
    double *A;
    if ((m < 1) || (n < 1))
    {
        printf("** Error: invalid parameter **\n");
        exit(1);
    }
    A = (double *)malloc((size_t)m * n * sizeof(double));
    if (A == NULL)
    {
        printf("** Error: insufficient memory **\n");
        exit(1);
    }
    return A;
}

// Takes pointer-to-mtrx so it can NULL the field after freeing
void freem(mtrx *A)
{
    if (A == NULL || A->M == NULL) return;
    free(A->M);
    A->M = NULL;
}

double *readm(char *filename, int *m, int *n)
{
    int i, j;
    FILE *f;
    double *A;
    f = fopen(filename, "r");
    fscanf(f, "%d", m);
    fscanf(f, "%d", n);
    A = allocm(*m, *n);
    for (i = 0; i < *m; i++)
        for (j = 0; j < *n; j++)
            fscanf(f, "%lf", &A[i * (*n) + j]);
    fclose(f);
    return A;
}

void printm(mtrx A)
{
    int i, j;
    printf("\n");
    for (i = 0; i < A.m; i++)
    {
        printf("[");
        for (j = 0; j < A.n; j++)
        {
            if (j == A.n - 1)
                printf(" %.4lf ", MAt(A, i, j));
            else
                printf(" %.4lf", MAt(A, i, j));
        }
        printf("]\n");
    }
}

void zerosv(vec v)
{
    memset(v.v, 0, v.n * sizeof(double));
}

double *allocv(int n)
{
    double *v;
    if (n < 1)
    {
        printf("** Error: invalid parameter **\n");
        exit(1);
    }
    v = (double *)malloc(n * sizeof(double));
    if (v == NULL)
    {
        printf("** Error: insufficient memory **\n");
        exit(1);
    }
    return v;
}

double *freev(vec v)
{
    if (v.v == NULL) return NULL;
    free(v.v);
    return NULL;
}

double *readv(char *filename, int *n)
{
    int i;
    FILE *f;
    double *v;
    f = fopen(filename, "r");
    fscanf(f, "%d", n);
    v = allocv(*n);
    for (i = 0; i < *n; i++)
        fscanf(f, "%lf", &v[i]);
    fclose(f);
    return v;
}

void printv(vec v)
{
    int i;
    printf("\n[ ");
    for (i = 0; i < v.n; i++)
    {
        if ((i == 0) || (i == (v.n - 1)))
            printf("%.4lf", v.v[i]);
        else
            printf(" %.4lf ", v.v[i]);
    }
    printf(" ]\n");
}

mtrx mtrxmul(mtrx A, mtrx B)
{
    int i, j, k;
    mtrx C;

    if (A.n != B.m)
    {
        printf("** Error: incompatible dimensions for matrix multiply **\n");
        printf("A columns: %d, B rows: %d\n", A.n, B.m);
        exit(1);
    }

    C = initm(A.m, B.n);
    for (i = 0; i < C.m; i++)
        for (j = 0; j < C.n; j++)
        {
            double sum = 0.0;
            for (k = 0; k < A.n; k++)
                sum += MAt(A, i, k) * MAt(B, k, j);
            MAt(C, i, j) = sum;
        }
    return C;
}

vec gaussian(mtrx A, vec b)
{
    int i, j, m, n, k;
    // Augmented matrix as a flat array
    int cols = b.n + 1;
    double *a = (double *)malloc((size_t)b.n * cols * sizeof(double));
    if (!a) { printf("** Error: insufficient memory **\n"); exit(1); }

    for (i = 0; i < b.n; i++)
    {
        for (j = 0; j < b.n; j++)
            a[i * cols + j] = MAt(A, i, j);
        a[i * cols + b.n] = b.v[i];
    }
    m = b.n;
    n = cols;

    vec x;
    x.v = allocv(m);
    x.n = m;

    for (i = 0; i < m - 1; i++)
    {
        // Partial pivoting
        for (k = i + 1; k < m; k++)
        {
            if (fabs(a[i * n + i]) < fabs(a[k * n + i]))
            {
                for (j = 0; j < n; j++)
                {
                    double temp    = a[i * n + j];
                    a[i * n + j]   = a[k * n + j];
                    a[k * n + j]   = temp;
                }
            }
        }
        // Elimination
        for (k = i + 1; k < m; k++)
        {
            double term = a[k * n + i] / a[i * n + i];
            for (j = 0; j < n; j++)
                a[k * n + j] -= term * a[i * n + j];
        }
    }
    // Back-substitution
    for (i = m - 1; i >= 0; i--)
    {
        x.v[i] = a[i * n + (n - 1)];
        for (j = i + 1; j < m; j++)
            x.v[i] -= a[i * n + j] * x.v[j];
        x.v[i] /= a[i * n + i];
    }
    free(a);
    return x;
}

mtrx kronecker(mtrx A, mtrx B)
{
    int i, j, ia, ib;
    int m = A.m * B.m;
    int n = A.n * B.n;
    mtrx C = initm(m, n);

    for (ia = 0; ia < A.m; ia++)
        for (ib = 0; ib < B.m; ib++)
            for (i = 0; i < A.n; i++)  // col block in A
                for (j = 0; j < B.n; j++)
                    MAt(C, ia * B.m + ib, i * B.n + j) =
                        MAt(A, ia, i) * MAt(B, ib, j);
    return C;
}

mtrx reshape(mtrx A, int m, int n)
{
    if (A.m * A.n != m * n)
    {
        printf("** Error: reshape element count mismatch: %d vs %d **\n",
               A.m * A.n, m * n);
        exit(1);
    }
    mtrx B = initm(m, n);
    memcpy(B.M, A.M, (size_t)m * n * sizeof(double));
    return B;
}

mtrx eye(int n)
{
    int i;
    mtrx A = initm(n, n);
    for (i = 0; i < n; i++)
        MAt(A, i, i) = 1.0;
    return A;
}

mtrx initm(int m, int n)
{
    mtrx A;
    A.M = allocm(m, n);
    A.m = m;
    A.n = n;
    zerosm(A);
    return A;
}

void invsig(mtrx A)
{
    int i, total = A.m * A.n;
    for (i = 0; i < total; i++)
        A.M[i] = -A.M[i];
}

double maxel(mtrx A)
{
    int i, total = A.m * A.n;
    double max = -__DBL_MAX__;
    for (i = 0; i < total; i++)
        if (A.M[i] > max) max = A.M[i];
    return max;
}

double minel(mtrx A)
{
    int i, total = A.m * A.n;
    double min = __DBL_MAX__;
    for (i = 0; i < total; i++)
        if (A.M[i] < min) min = A.M[i];
    return min;
}

void mtrxcpy(mtrx A, mtrx B)
{
    memcpy(A.M, B.M, (size_t)A.m * A.n * sizeof(double));
}

// ---------------------------------------------------------------------------
// Sparse matrix (CSR) operations
// ---------------------------------------------------------------------------

smtrx initsm(int m, int n, int nnz)
{
    smtrx A;
    A.m = m; A.n = n; A.nnz = nnz;
    A.values  = (double *)malloc(nnz * sizeof(double));
    A.col_idx = (int *)   malloc(nnz * sizeof(int));
    A.row_ptr = (int *)   malloc((m + 1) * sizeof(int));
    if (!A.values || !A.col_idx || !A.row_ptr)
    {
        printf("** Error: insufficient memory for sparse matrix **\n");
        exit(1);
    }
    return A;
}

void freesm(smtrx A)
{
    free(A.values);
    free(A.col_idx);
    free(A.row_ptr);
}

void spmv(smtrx A, double *x, double *y)
{
    int i;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (i = 0; i < A.m; i++)
    {
        double sum = 0.0;
        int k;
        for (k = A.row_ptr[i]; k < A.row_ptr[i + 1]; k++)
            sum += A.values[k] * x[A.col_idx[k]];
        y[i] = sum;
    }
}

smtrx seye(int n)
{
    int i;
    smtrx I = initsm(n, n, n);
    for (i = 0; i < n; i++)
    {
        I.values[i]  = 1.0;
        I.col_idx[i] = i;
        I.row_ptr[i] = i;
    }
    I.row_ptr[n] = n;
    return I;
}

smtrx skronecker(smtrx A, smtrx B)
{
    int ia, ib, ka, kb;
    smtrx C = initsm(A.m * B.m, A.n * B.n, A.nnz * B.nnz);
    int pos = 0, row = 0;

    for (ia = 0; ia < A.m; ia++)
        for (ib = 0; ib < B.m; ib++)
        {
            C.row_ptr[row] = pos;
            for (ka = A.row_ptr[ia]; ka < A.row_ptr[ia + 1]; ka++)
                for (kb = B.row_ptr[ib]; kb < B.row_ptr[ib + 1]; kb++)
                {
                    C.values[pos]  = A.values[ka] * B.values[kb];
                    C.col_idx[pos] = A.col_idx[ka] * B.n + B.col_idx[kb];
                    pos++;
                }
            row++;
        }
    C.row_ptr[row] = pos;
    return C;
}

void flatten(mtrx A, double *out, int nx, int ny)
{
    // With flat storage A.M is already row-major — just memcpy
    (void)nx; (void)ny; // nx*ny == A.m*A.n by contract
    memcpy(out, A.M, (size_t)A.m * A.n * sizeof(double));
}

void unflatten(double *in, mtrx A, int nx, int ny)
{
    (void)nx; (void)ny;
    memcpy(A.M, in, (size_t)A.m * A.n * sizeof(double));
}
