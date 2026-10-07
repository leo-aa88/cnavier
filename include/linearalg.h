// Linear algebra library

#ifndef LINEARALG_H_INCLUDED
#define LINEARALG_H_INCLUDED

typedef struct matrix
{
    double *M; // flat row-major storage: element [i][j] -> M[i*n + j]
    int m;     // number of rows
    int n;     // number of columns
} mtrx;

// Element access macro — use instead of A.M[i][j]
#define MAt(A, i, j) ((A).M[(i) * (A).n + (j)])

typedef struct vector
{
    double *v; // vector pointer
    int n;     // number of elements
} vec;

void zerosm(mtrx A);                // Initialize matrix with zeros
double *allocm(int m, int n);       // Allocate flat matrix buffer m x n
void freem(mtrx *A);                // Free matrix memory and NULL the pointer
void printm(mtrx A);                // Print matrix
void zerosv(vec v);                 // Initialize vector with zeros
double *allocv(int n);              // Allocate vector with size n
double *freev(vec v);               // Free memory for vector
void printv(vec v);                 // Print vector
mtrx mtrxmul(mtrx A, mtrx B);       // Matrix multiplication
vec gaussian(mtrx A, vec b);        // Gaussian elimination
mtrx kronecker(mtrx A, mtrx B);     // Kronecker matrix product
mtrx reshape(mtrx A, int m, int n); // Reshape matrix
mtrx eye(int n);                    // Generates identity matrix
mtrx initm(int m, int n);           // Allocate and zero-initialise matrix
void invsig(mtrx A);                // Negate all matrix elements
double maxel(mtrx A);               // Return max element
double minel(mtrx A);               // Return min element
void mtrxcpy(mtrx A, mtrx B);       // Copy B into A
void negcpy(mtrx A, mtrx B);        // A = -B (A may be B)

// OpenMP builds run a loop serially when the grid has fewer than this many
// points (about 45x45): below it, waking the threads costs more than the loop
// saves. The comparison is on grid points, whatever the loop iterates over
// (rows, single points, or batches of transforms).
#ifndef OMP_MIN_WORK
#define OMP_MIN_WORK 2048
#endif

// Sparse matrix in Compressed Sparse Row (CSR) format
typedef struct sparse_matrix
{
    double *values; // non-zero values        [nnz]
    int *col_idx;   // column index per value [nnz]
    int *row_ptr;   // row i spans values[row_ptr[i]..row_ptr[i+1]-1] [m+1]
    int nnz;        // number of non-zero entries
    int m;          // number of rows
    int n;          // number of columns
} smtrx;

smtrx initsm(int m, int n, int nnz);      // Allocate sparse matrix (nnz slots)
void freesm(smtrx A);                     // Free sparse matrix memory
void spmv(smtrx A, double *x, double *y); // y = A*x  (y must be pre-allocated, size m)
smtrx skronecker(smtrx A, smtrx B);       // Sparse Kronecker product A x B
smtrx seye(int n);                        // Sparse n x n identity matrix
// Flatten row-major nx*ny dense matrix into pre-allocated flat array
void flatten(mtrx A, double *out, int nx, int ny);
// Write flat array back into pre-allocated dense matrix (row-major)
void unflatten(double *in, mtrx A, int nx, int ny);

#endif // LINEARALG_H_INCLUDED
