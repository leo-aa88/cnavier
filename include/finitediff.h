// Finite difference library

#ifndef FINITEDIFF_H_INCLUDED
#define FINITEDIFF_H_INCLUDED

#include "linearalg.h"

smtrx SDiff1(int n, int o, double dx); // Sparse finite-difference matrix, first derivative
smtrx SDiff2(int n, int o, double dx); // Sparse finite-difference matrix, second derivative

// The same derivatives on a periodic grid of n points x_i = i*dx, i = 0..n-1
// (x_n is x_0): every row is the centered stencil of order o, wrapped around
smtrx SDiff1_periodic(int n, int o, double dx);
smtrx SDiff2_periodic(int n, int o, double dx);

#endif // FINITEDIFF_H_INCLUDED
