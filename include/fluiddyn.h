// Fluid dynamics library

#ifndef FLUIDDYN_H_INCLUDED
#define FLUIDDYN_H_INCLUDED

#include "linearalg.h"

void euler(mtrx w, mtrx dwdx, mtrx dwdy, mtrx d2wdx2, mtrx d2wdy2, mtrx u, mtrx v, double Re, double dt); // Euler time-advancement
mtrx continuity(mtrx dudx, mtrx dvdy);                                                                    // computes continuity equation
mtrx vorticity(mtrx dvdx, mtrx dudy);                                                                     // computes vorticity
mtrx pressure(void);                                                                                      // computes pressure


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
    double *flat_w, *flat_psi, *flat_tmp; // flat work buffers
    int    nx, ny;
    int    poisson_type;
    int    poisson_max_it;
    double poisson_tol, beta, dx, dy, Re;
} rk4_ctx;

// Allocate all RK4 workspace for an nx*ny grid
rk4_ctx rk4_alloc(int nx, int ny);
// Free all RK4 workspace
void rk4_free(rk4_ctx *ctx);

// Evaluate vorticity RHS: dw/dt = -u*dw/dx - v*dw/dy + (1/Re)*(d2w/dx2 + d2w/dy2)
// Writes result into out. Updates u and v via Poisson solve for the given w.
void dwdt(mtrx w, mtrx u, mtrx v, mtrx out, rk4_ctx *ctx);

// RK4 time advancement — advances w, u, v by dt
void rk4(mtrx w, mtrx u, mtrx v, double dt, rk4_ctx *ctx);

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
