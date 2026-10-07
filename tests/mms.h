// Manufactured solution for verifying the spatial discretization
//
//   psi(x, y, t) = g(t) sin^2(pi x / Lx) sin^2(pi y / Ly)
//
// psi = 0 and u = v = 0 on every wall, so it satisfies the solver's
// stationary no-slip walls and its psi = 0 Poisson condition. The vorticity source
// that makes it an exact solution of
//
//   dw/dt = -u dw/dx - v dw/dy + (1/Re) lap(w) + f,   w = -lap(psi),
//   u = dpsi/dy, v = -dpsi/dx
//
// is f = dw/dt + u dw/dx + v dw/dy - (1/Re) lap(w), in closed form.
//
// On a doubly periodic grid the solution is instead a sum of three Fourier
// modes of different wavenumbers (so that the nonlinear term does not vanish):
//
//   psi = g(t) sum_m A_m cos(2 pi (a_m x / Lx + b_m y / Ly) + phi_m)

#ifndef MMS_H_INCLUDED
#define MMS_H_INCLUDED

#include "linearalg.h"

typedef struct
{
    double Lx, Ly; // domain size
    double Re;     // Reynolds number
    double dx, dy; // grid spacing: node (i, j) is at x = j*dx, y = i*dy
    int periodic;  // 0: the wall solution, 1: the periodic one
} mms_case;

// Exact w, u, v and psi at time t on the case's grid. NULL arguments are skipped.
void mms_exact(const mms_case *c, double t, mtrx *w, mtrx *u, mtrx *v, mtrx *psi);

// The vorticity source at time t, in the form solver_config.vorticity_source expects;
// data is a const mms_case *
void mms_source(double t, mtrx f, void *data);

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

// Run the solver from the exact state at t = t0 to t = t0 + T on an nx x ny grid
// with the given interior derivative order, time scheme (1 Euler, 2 RK4) and
// Poisson solver (1-3), and compare with the exact solution at the final
// time. The step is at most dt and at most half the viscous stability limit,
// shortened so that a whole number of steps reaches T.
mms_errors mms_run(int nx, int ny, double Lx, double Ly, double Re, int order, int time_scheme,
                   int poisson_type, double dt, double t0, double T);

// mms_run() with RK4 and the FFT solver, with a choice of Poisson operator
// (poisson_order 2 or 4) and wall-vorticity closure (0: D_x v - D_y u, 1:
// third order from psi)
mms_errors mms_run_closures(int nx, int ny, double Lx, double Ly, double Re, int order, int poisson_order,
                            int wall_closure, double dt, double t0, double T);

// mms_run() on a doubly periodic nx x ny grid (dx = Lx/nx) with the periodic
// solution, RK4 or Euler, and the periodic FFT Poisson solver
mms_errors mms_run_periodic(int nx, int ny, double Lx, double Ly, double Re, int order, int time_scheme,
                            double dt, double t0, double T);

// Ablation of the discretization (issue #25): parts of the RK4 step replaced
// by the exact solution at each stage's time. Flags may be combined.
#define MMS_EXACT_WALL_W 1   // exact wall vorticity instead of D_x v - D_y u
#define MMS_EXACT_PSI 2      // exact psi instead of the Poisson solve; u, v from D_y psi, -D_x psi
#define MMS_EXACT_VELOCITY 4 // exact u, v at every node (no psi at all)

// mms_run() with RK4 and the FFT solver, done by a copy of step() in which
// the parts named by `flags` are replaced with the exact solution. With
// flags = 0 it computes exactly what mms_run() does, to the last bit.
mms_errors mms_run_ablated(int nx, int ny, double Lx, double Ly, double Re, int order, int flags,
                           double dt, double t0, double T);

#endif // MMS_H_INCLUDED
