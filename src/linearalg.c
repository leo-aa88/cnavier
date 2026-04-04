#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linearalg.h"

void zerosm(mtrx A)
{
    int i, j;
    for (i = 0; i < A.m; i++)
    {
        for (j = 0; j < A.n; j++)
        {
            A.M[i][j] = 0;
        }
    }
}

double **allocm(int m, int n)
{
    int i;
    double **A;

    if ((m < 1) || (n < 1))
    {
        printf("** Error: invalid parameter **\n");
        exit(1);
    }

    A = (double **)malloc(m * sizeof(double *));

    if (A == NULL)
    {
        printf("** Error: insufficient memory **");
        exit(1);
    }

    for (i = 0; i < m; i++)
    {
        A[i] = (double *)malloc(n * sizeof(double));
        if (A[i] == NULL)
        {
            printf("** Error: insufficient memory **");
            exit(1);
        }
    }
    return (A);
}

double **freem(mtrx A)
{
    int i;
    if (A.M == NULL)
        return (NULL);
    if ((A.m < 1) || (A.n < 1))
    {
        printf("** Error: invalid parameter **\n");
        exit(1);
    }
    for (i = 0; i < A.m; i++)
        free(A.M[i]);
    free(A.M);
    //printf("Liberated successfully\n");
    return (NULL);
}

double **readm(char *filename, int *m, int *n)
{
    int i, j;
    FILE *f;
    double **A;
    f = fopen(filename, "r"); // opens file
    fscanf(f, "%d", m);       // read row size of matrix
    fscanf(f, "%d", n);       // read col size of matrix
    A = allocm(*m, *n);       // allocate memory for matrix
    for (i = 0; i < *m; i++)
    {
        for (j = 0; j < *n; j++)
        {
            fscanf(f, "%lf", &A[i][j]);
        }
    }
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
            {
                printf(" %.4lf ", A.M[i][j]);
            }
            else
            {
                printf(" %.4lf", A.M[i][j]);
            }
        }
        if (i == A.n - 1)
        {
            printf("]");
        }
        else
        {
            printf("]");
        }
        printf("\n");
    }
}

void zerosv(vec v)
{
    int i;
    for (i = 0; i < v.n; i++)
    {
        v.v[i] = 0;
    }
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
        printf("** Error: insufficient memory **");
        exit(1);
    }
    return (v);
}

double *freev(vec v)
{
    if (v.v == NULL)
        return (NULL);
    if (v.n < 1)
    {
        printf("** Error: invalid parameter **\n");
        exit(1);
    }
    free(v.v);
    return (NULL);
}

double *readv(char *filename, int *n)
{
    int i;
    FILE *f;
    double *v;
    f = fopen(filename, "r"); // opens file
    fscanf(f, "%d", n);       // read size of vector
    v = allocv(*n);           // allocate memory for matrix
    for (i = 0; i < *n; i++)
    {
        fscanf(f, "%lf", &v[i]);
    }
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
        {
            printf("%.4lf", v.v[i]);
        }
        else
        {
            printf(" %.4lf ", v.v[i]);
        }
    }
    printf(" ]\n");
}

mtrx mtrxmul(mtrx A, mtrx B)
{
    mtrx C;
    int i, j, k;

    if (A.n != B.m)
    {
        printf("** Error: the first matrix number of columns must be equal to the second matrix number of rows **\n");
        printf("Columns of first matrix: %d\n", A.n);
        printf("Rows of second matrix: %d\n", B.m);
        exit(1);
    }

    C = initm(A.m, B.n);

    for (i = 0; i < C.m; i++)
    {
        for (j = 0; j < C.n; j++)
        {
            C.M[i][j] = 0;
            for (k = 0; k < B.m; k++)
            {
                C.M[i][j] = C.M[i][j] + A.M[i][k] * B.M[k][j];
            }
        }
    }
    return C;
}

vec gaussian(mtrx A, vec b)
{
    int i, j, m, n, k;

    // Augmented matrix
    double **a;
    a = allocm(b.n, b.n + 1);
    for (i = 0; i < b.n; i++)
    {
        for (j = 0; j < b.n; j++)
        {
            a[i][j] = A.M[i][j];
        }
    }
    for (i = 0; i < b.n; i++)
    {
        a[i][b.n] = b.v[i];
    }
    m = b.n;
    n = b.n + 1;

    vec x;
    x.v = allocv(n - 1);
    x.n = n - 1;

    for (i = 0; i < m - 1; i++)
    {
        // Partial Pivoting
        for (k = i + 1; k < m; k++)
        {
            // If diagonal element(absolute vallue) is smaller than any of the terms below it
            if (fabs(a[i][i]) < fabs(a[k][i]))
            {
                // Swap the rows
                for (j = 0; j < n; j++)
                {
                    double temp;
                    temp = a[i][j];
                    a[i][j] = a[k][j];
                    a[k][j] = temp;
                }
            }
        }
        // Begin Gauss Elimination
        for (k = i + 1; k < m; k++)
        {
            double term = a[k][i] / a[i][i];
            for (j = 0; j < n; j++)
            {
                a[k][j] = a[k][j] - term * a[i][j];
            }
        }
    }
    // Begin Back-substitution
    for (i = m - 1; i >= 0; i--)
    {
        x.v[i] = a[i][n - 1];
        for (j = i + 1; j < n - 1; j++)
        {
            x.v[i] = x.v[i] - a[i][j] * x.v[j];
        }
        x.v[i] = x.v[i] / a[i][i];
    }
    return x;
}

mtrx kronecker(mtrx A, mtrx B)
{
    int i, j;
    int n = A.n * B.n;
    mtrx C;
    C.M = allocm(n, n);
    C.m = n;
    C.n = n;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            C.M[i][j] = A.M[i / A.n][j / A.n] * B.M[i % B.n][j % B.n];
        }
    }
    return C;
}

mtrx reshape(mtrx A, int m, int n)
{
    mtrx B;
    int i, j, k, l;

    if ((A.m * A.n) != (m * n))
    {
        printf("** Error: the reshaped matrix must have the same number of elements **\n");
        printf("Number of elements of input matrix: %d\n", A.m * A.n);
        printf("Number of elements of output matrix: %d\n", m * n);
        exit(1);
    }

    B = initm(m, n);

    k = 0;
    l = 0;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            B.M[i][j] = A.M[k][l];
            if (l < (A.n - 1))
            {
                l++;
            }
            else
            {
                k++;
                l = 0;
            }
        }
    }
    return B;
}

mtrx eye(int n)
{
    int i, j;
    mtrx A;
    A.M = allocm(n, n);
    A.m = n;
    A.n = n;
    zerosm(A);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i == j)
            {
                A.M[i][j] = 1;
            }
        }
    }
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
    int i, j;
    for (i = 0; i < A.m; i++)
    {
        for (j = 0; j < A.n; j++)
        {
            A.M[i][j] = -A.M[i][j];
        }
    }
}

double maxel(mtrx A)
{
    int i, j;
    double max_element = -__DBL_MAX__;

    for (i = 0; i < A.m; i++)
    {
        for (j = 0; j < A.n; j++)
        {
            if (A.M[i][j] > max_element)
            {
                max_element = A.M[i][j];
            }
        }
    }
    return max_element;
}

double minel(mtrx A)
{
    int i, j;
    double min_element = __DBL_MAX__;

    for (i = 0; i < A.m; i++)
    {
        for (j = 0; j < A.n; j++)
        {
            if (A.M[i][j] < min_element)
            {
                min_element = A.M[i][j];
            }
        }
    }
    return min_element;
}

void mtrxcpy(mtrx A, mtrx B)
{
    int i, j;

    for (i = 0; i < A.m; i++)
    {
        for (j = 0; j < A.n; j++)
        {
            A.M[i][j] = B.M[i][j];
        }
    }
}
// ---------------------------------------------------------------------------
// Sparse matrix (CSR) operations
// ---------------------------------------------------------------------------

smtrx initsm(int m, int n, int nnz)
{
    smtrx A;
    A.m = m;
    A.n = n;
    A.nnz = nnz;
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

// y = A * x   (y pre-allocated by caller, fully overwritten)
void spmv(smtrx A, double *x, double *y)
{
    int i, k;
    for (i = 0; i < A.m; i++)
    {
        double sum = 0.0;
        for (k = A.row_ptr[i]; k < A.row_ptr[i + 1]; k++)
            sum += A.values[k] * x[A.col_idx[k]];
        y[i] = sum;
    }
}

// Sparse n x n identity matrix
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

// Sparse Kronecker product: C = A x B
// C is (A.m*B.m) x (A.n*B.n) with A.nnz*B.nnz non-zeros.
smtrx skronecker(smtrx A, smtrx B)
{
    int ia, ib, ka, kb;
    int total_nnz = A.nnz * B.nnz;
    int m = A.m * B.m;
    int n = A.n * B.n;
    smtrx C = initsm(m, n, total_nnz);

    int pos = 0;
    int row = 0;

    for (ia = 0; ia < A.m; ia++)
    {
        for (ib = 0; ib < B.m; ib++)
        {
            C.row_ptr[row] = pos;
            for (ka = A.row_ptr[ia]; ka < A.row_ptr[ia + 1]; ka++)
            {
                for (kb = B.row_ptr[ib]; kb < B.row_ptr[ib + 1]; kb++)
                {
                    C.values[pos]  = A.values[ka] * B.values[kb];
                    C.col_idx[pos] = A.col_idx[ka] * B.n + B.col_idx[kb];
                    pos++;
                }
            }
            row++;
        }
    }
    C.row_ptr[row] = pos;
    return C;
}

// Flatten dense (nx x ny) matrix into pre-allocated flat array (row-major)
void flatten(mtrx A, double *out, int nx, int ny)
{
    int i, j;
    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
            out[i * ny + j] = A.M[i][j];
}

// Write flat array back into pre-allocated dense matrix (row-major)
void unflatten(double *in, mtrx A, int nx, int ny)
{
    int i, j;
    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
            A.M[i][j] = in[i * ny + j];
}
