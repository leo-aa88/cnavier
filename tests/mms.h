// Manufactured solution for verifying the spatial discretization
//
//   psi(x, y, t) = g(t) sin^2(pi x / Lx) sin^2(pi y / Ly)
//
// psi = 0 and u = v = 0 on every wall, so it satisfies the solver's
// stationary no-slip walls and its psi = 0 Poisson condition. The body force
// that makes it an exact solution of
//
//   dw/dt = -u dw/dx - v dw/dy + (1/Re) lap(w) + f,   w = -lap(psi),
//   u = dpsi/dy, v = -dpsi/dx
//
// is f = dw/dt + u dw/dx + v dw/dy - (1/Re) lap(w), in closed form.

#ifndef MMS_H_INCLUDED
#define MMS_H_INCLUDED

#include "linearalg.h"

typedef struct
{
    double Lx, Ly; // domain size
    double Re;     // Reynolds number
    double dx, dy; // grid spacing: node (i, j) is at x = j*dx, y = i*dy
} mms_case;

// Exact w, u, v and psi at time t on the case's grid. NULL arguments are skipped.
void mms_exact(const mms_case *c, double t, mtrx *w, mtrx *u, mtrx *v, mtrx *psi);

// The body force at time t, in the form solver_config.forcing expects;
// data is a const mms_case *
void mms_forcing(double t, mtrx f, void *data);

// Errors of a run against the exact solution: the largest and the root-mean-
// square difference over the nodes. w is measured on all nodes and,
// separately, on the wall and on the interior nodes; u and v with the wall
// velocities (zero) imposed.
typedef struct
{
    double max, rms;
} mms_norm;

typedef struct
{
    mms_norm psi, u, v, w, w_wall, w_interior;
    int steps;
} mms_errors;

// Run the solver from the exact state at t = 0 to t = T on an nx x ny grid
// with the given interior derivative order, time scheme (1 Euler, 2 RK4) and
// Poisson solver (1-3), and compare with the exact solution at the final
// time. The step is at most dt and at most half the viscous stability limit,
// shortened so that a whole number of steps reaches T.
mms_errors mms_run(int nx, int ny, double Lx, double Ly, double Re, int order, int time_scheme,
                   int poisson_type, double dt, double T);

#endif // MMS_H_INCLUDED
