// Fluid dynamics library

#ifndef FLUIDDYN_H_INCLUDED
#define FLUIDDYN_H_INCLUDED

#include "linearalg.h"

void euler(mtrx w, mtrx dwdx, mtrx dwdy, mtrx d2wdx2, mtrx d2wdy2, mtrx u, mtrx v, double Re, double dt); // Euler time-advancement


// RK4 context — holds all workspace needed to evaluate the vorticity RHS
// at each stage without allocating inside the loop.
typedef struct {
    smtrx *DX, *DY, *DX2, *DY2; // sparse derivative operators
    mtrx   dwdx, dwdy;           // first derivatives of w
    mtrx   d2wdx2, d2wdy2;       // second derivatives of w
    mtrx   dpsidx, dpsidy;       // stream-function derivatives
    mtrx   psi, psi_scratch;     // Poisson solution and scratch
    mtrx   k1, k2, k3, k4;      // RK4 stage increments
    mtrx   w_tmp;                // temporary w for intermediate stages
    mtrx   rhs;                  // Poisson right-hand side, -w
    int    nx, ny;
    int    poisson_type;
    int    poisson_max_it;
    double poisson_tol, beta, dx, dy, Re;
} rk4_ctx;

// Allocate all RK4 workspace for an nx*ny grid. Only square grids are
// supported: exits with an error unless nx == ny.
rk4_ctx rk4_alloc(int nx, int ny);
// Free all RK4 workspace
void rk4_free(rk4_ctx *ctx);

// Evaluate vorticity RHS: dw/dt = -u*dw/dx - v*dw/dy + (1/Re)*(d2w/dx2 + d2w/dy2)
// Writes result into out. Updates u and v via Poisson solve for the given w.
void dwdt(mtrx w, mtrx u, mtrx v, mtrx out, rk4_ctx *ctx);

// RK4 time advancement — advances w, u, v by dt
void rk4(mtrx w, mtrx u, mtrx v, double dt, rk4_ctx *ctx);

// Largest stable time step of the explicit scheme (time_scheme 1=Euler,
// 2=RK4) for the viscous term (1/Re) * (d_xx + d_yy), from the 1D
// second-derivative operators in x and y.
double max_stable_dt(const smtrx *dxx, const smtrx *dyy, double Re, int time_scheme);

// The stability limit of forward Euler with centered advection at speeds up
// to u_max: dt <= 2 nu / u_max^2 = 2 / (Re u_max^2). It assumes u_max
// everywhere, so it is conservative for the cavity (at Re = 1000 runs stayed
// stable up to about 3x it). Returns HUGE_VAL if u_max is 0.
double euler_advection_dt(double Re, double u_max);

// The time-step limits main() checks before a run. h is the smaller grid
// spacing, u_max the fastest wall speed, max_co the largest Courant number.
typedef struct
{
    double courant;   // max_co * h / u_max (HUGE_VAL if the walls are at rest)
    double viscous;   // max_stable_dt()
    double advection; // euler_advection_dt() for Euler, HUGE_VAL for RK4
    double accept;    // largest dt that is run: min(courant, viscous)
    double suggest;   // dt to suggest: min of all three. The advection limit
                      // is conservative, so a dt above it only gets a warning,
                      // but a suggestion should not be one that is warned about
} dt_limits;
dt_limits time_step_limits(const smtrx *dxx, const smtrx *dyy, double h, double Re, double u_max,
                           double max_co, int time_scheme);

// x rounded down to three significant digits, for suggesting a time step: the
// value printed with %.3g is then never above x.
double round_down_3(double x);

// Dirichlet wall velocities, indexed by wall:
// 0 -> j=0, 1 -> j=ny-1, 2 -> i=0, 3 -> i=nx-1
typedef struct {
    double u[4];
    double v[4];
} wall_bc;

// Impose the wall velocities on u and v
void apply_wall_bc(mtrx u, mtrx v, const wall_bc *bc);

// One full timestep: wall BCs, vorticity BCs, time advancement
// (time_scheme 1=Euler, 2=RK4), Poisson solve and velocity recovery.
void step(mtrx w, mtrx u, mtrx v, double dt, int time_scheme, const wall_bc *bc, rk4_ctx *ctx);

#endif // FLUIDDYN_H_INCLUDED
