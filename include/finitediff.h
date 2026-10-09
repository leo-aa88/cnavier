// Finite difference library

#ifndef FINITEDIFF_H_INCLUDED
#define FINITEDIFF_H_INCLUDED

#include "linearalg.h"

smtrx SDiff1(int n, int o, double dx); // Sparse finite-difference matrix, first derivative
smtrx SDiff2(int n, int o, double dx); // Sparse finite-difference matrix, second derivative

// SDiff1 with fourth-order one-sided rows at both ends and next to them, for
// recovering the velocity from the stream function (orders 4 and 6; order 2
// is SDiff1). Needs n >= 10.
smtrx SDiff1_wall4(int n, int o, double dx);

// The same derivatives on a periodic grid of n points x_i = i*dx, i = 0..n-1
// (x_n is x_0): every row is the centered stencil of order o, wrapped around.
// o = FD_COMPACT6: Lele's sixth-order tridiagonal compact schemes instead,
// A f' = B f, applied as the circulant matrix A^-1 B (periodic grids only)
#define FD_COMPACT6 106
// o = FD_SPECTRAL: pseudospectral differentiation, symbols i k and -k^2. The
// matrices are dense; the solver applies them in Fourier space
// (solver_config.fourier), and they exist for their symbols
#define FD_SPECTRAL 199
smtrx SDiff1_periodic(int n, int o, double dx);
smtrx SDiff2_periodic(int n, int o, double dx);

// Symbol of the compact schemes at theta = k h, in units of 1/h (deriv 1) or
// 1/h^2 (deriv 2): k* h for the first derivative (whose symbol is i k*), and
// the real symbol -(k* h)^2 of the second
double compact6_symbol(int deriv, double theta);

#endif // FINITEDIFF_H_INCLUDED
