// Linear algebra library

#ifndef LINEARALG_H_INCLUDED
#define LINEARALG_H_INCLUDED

typedef struct matrix
{
    double **M; // matrix double pointer
    int m;      // number of rows
    int n;      // number of columns
} mtrx;

typedef struct vector
{
    double *v; // vector pointer
    int n;     // number of elements
} vec;

void zerosm(mtrx A);                            // Initialize matrix with zeros
double **allocm(int m, int n);                  // Allocate matrix with size m x n
double **freem(mtrx A);                         // Free memory for matrix
double **readm(char *filename, int *m, int *n); // Read matrix
void printm(mtrx A);                            // Print matrix
void zerosv(vec v);                             // Initialize vector with zeros
double *allocv(int n);                          // Allocate vector with size n
double *freev(vec v);                           // Free memory for vector
double *readv(char *filename, int *n);          // Read vector
void printv(vec v);                             // Print vector;
mtrx mtrxmul(mtrx A, mtrx B);                   // Matrix multiplication
vec gaussian(mtrx A, vec b);                    // Gaussian elimination
mtrx kronecker(mtrx A, mtrx B);                 // Kronecker matrix product
mtrx reshape(mtrx A, int m, int n);             // Reshape matrix
mtrx eye(int n);                                // Generates identiy matrix
mtrx initm(int m, int n);                       // Allocates memory and return matrix
void invsig(mtrx A);                            // Inverts signal of all matrix elements
double maxel(mtrx A);                           // Returns max element from matrix
double minel(mtrx A);                           // Returns min element from matrix
void mtrxcpy(mtrx A, mtrx B);                   // Copy matrix B into A

// Sparse matrix in Compressed Sparse Row (CSR) format
typedef struct sparse_matrix
{
    double *values;  // non-zero values        [nnz]
    int    *col_idx; // column index per value [nnz]
    int    *row_ptr; // row i spans values[row_ptr[i]..row_ptr[i+1]-1] [m+1]
    int     nnz;     // number of non-zero entries
    int     m;       // number of rows
    int     n;       // number of columns
} smtrx;

smtrx initsm(int m, int n, int nnz);       // Allocate sparse matrix (nnz slots)
void  freesm(smtrx A);                     // Free sparse matrix memory
void  spmv(smtrx A, double *x, double *y); // y = A*x  (y must be pre-allocated, size m)
smtrx skronecker(smtrx A, smtrx B);        // Sparse Kronecker product A ⊗ B
smtrx seye(int n);                         // Sparse n×n identity matrix
// Flatten row-major nx*ny dense matrix into pre-allocated flat array
void flatten(mtrx A, double *out, int nx, int ny);
// Write flat array back into pre-allocated dense matrix (row-major)
void unflatten(double *in, mtrx A, int nx, int ny);

#endif // LINEARALG_H_INCLUDED
