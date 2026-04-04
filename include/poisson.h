// Poisson solver library

#ifndef POISSON_H_INCLUDED
#define POISSON_H_INCLUDED

#include "linearalg.h"

#define PI 3.14159265359

double error(mtrx u1, mtrx u2);

// Solvers write the result into the pre-allocated matrix u.
// u0 is a same-sized scratch buffer, also pre-allocated by the caller.
// Both must be initm(f.m, f.n) before the first call.
void poisson(mtrx f, mtrx u, mtrx u0, double dx, double dy, int itmax, double tol);
void poisson_SOR(mtrx f, mtrx u, mtrx u0, double dx, double dy, int itmax, double tol, double beta);

#endif // POISSON_H_INCLUDED
