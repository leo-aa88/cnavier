// Poisson solver library

#ifndef POISSON_H_INCLUDED
#define POISSON_H_INCLUDED

#include "linearalg.h"
#include <fftw3.h>

#define PI 3.14159265359

double error(mtrx u1, mtrx u2);

// SOR relaxation parameter for the 5-point Laplacian on nx x ny nodes with
// spacing dx, dy (ny rows, nx columns), from the spectral radius of Jacobi on
// the interior: rho = (dy^2 cos(pi/(nx-1)) + dx^2 cos(pi/(ny-1))) / (dx^2 + dy^2),
// beta = 2 / (1 + sqrt(1 - rho^2)).
double sor_beta(int nx, int ny, double dx, double dy);

// Solvers write the result into the pre-allocated matrix u.
// u0 is a same-sized scratch buffer, also pre-allocated by the caller.
// Both must be initm(f.m, f.n) before the first call.
// They return the number of iterations taken.
int poisson(mtrx f, mtrx u, mtrx u0, double dx, double dy, int itmax, double tol);
int poisson_SOR(mtrx f, mtrx u, mtrx u0, double dx, double dy, int itmax, double tol, double beta);

// FFT-based direct Poisson solver (exact, O(n² log n), no iteration needed).
// Uses a DST-I (sine transform) of the interior nodes, so u = 0 on the wall
// nodes, as in the iterative solvers. Result written into pre-allocated
// matrix u. No scratch buffer needed.
void poisson_FFT(mtrx f, mtrx u, double dx, double dy);

// Call once at program start to pre-plan FFTW transforms for an nx*ny grid.
// Call fft_cleanup() at program end.
void fft_setup(int nx, int ny);
void fft_cleanup(void);

#endif // POISSON_H_INCLUDED

