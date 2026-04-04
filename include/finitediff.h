// Finite difference library

#ifndef FINITEDIFF_H_INCLUDED
#define FINITEDIFF_H_INCLUDED

#include "linearalg.h"

mtrx  Diff1(int n, int o, double dx);  // Dense finite-difference matrix, first derivative
smtrx SDiff1(int n, int o, double dx); // Sparse finite-difference matrix, first derivative
mtrx  Diff2(int n, int o, double dx);  // Dense finite-difference matrix, second derivative
smtrx SDiff2(int n, int o, double dx); // Sparse finite-difference matrix, second derivative

#endif // FINITEDIFF_H_INCLUDED
