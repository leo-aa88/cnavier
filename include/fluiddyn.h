// Fluid dynamics library

#ifndef FLUIDDYN_H_INCLUDED
#define FLUIDDYN_H_INCLUDED

#include "linearalg.h"

// Forward Euler update of w with the given derivatives and velocity; f is the
// vorticity source at the start of the step, or NULL for none
void euler(mtrx w, mtrx dwdx, mtrx dwdy, mtrx d2wdx2, mtrx d2wdy2, mtrx u, mtrx v, const double *f,
           double Re, double dt);

// Dirichlet wall velocities, indexed by wall:
// 0 -> left   (x = 0,  column j = 0)
// 1 -> right  (x = Lx, column j = nx-1)
// 2 -> bottom (y = 0,  row i = 0)
// 3 -> top    (y = Ly, row i = ny-1)   the lid of the default cavity
// At the corners the left and right walls take precedence.
typedef struct
{
    double u[4];
    double v[4];
} wall_bc;

// Everything that defines a run, shared by the CPU and GPU backends. Every
// solver takes a copy when it is created, so later changes to the caller's
// struct have no effect on it. The operators are owned by the caller and must
// outlive the solver.
typedef struct
{
    int nx, ny;                       // grid points in x and y
    double dx, dy;                    // grid spacing
    double Re;                        // Reynolds number
    double dt;                        // time step
    double t0;                        // time at the start of the run
    int time_scheme;                  // 1=Euler, 2=RK4
    int poisson_type;                 // 1=Gauss-Seidel, 2=SOR, 3=FFT
    int poisson_max_it;               // iteration limit of Gauss-Seidel/SOR
    int poisson_order;                // FFT solver with walls: 2 (5-point) or 4 (compact 9-point)
    int wall_closure;                 // wall vorticity: 0 = D_x v - D_y u, 1 = third-order formula from psi
    double poisson_tol, beta;         // their tolerance and SOR parameter
    int periodic;                     // 0: four walls with velocities bc; 1: doubly periodic
    wall_bc bc;                       // wall velocities (walls only)
    const smtrx *DX, *DY, *DX2, *DY2; // sparse derivative operators

    // Optional source term f in the vorticity equation,
    // dw/dt = -u.grad(w) + (1/Re) lap(w) + f (the curl of a body force in the
    // momentum equation), used to verify the solver against a manufactured
    // solution. vorticity_source(t, f, source_data) fills the ny x nx field f
    // with f(x, y, t); set it to NULL for none. The solver calls it once per
    // right-hand-side evaluation, at t = t0 + (steps taken) * dt plus the
    // stage offset; the GPU solver calls it on the host and copies the field
    // to the device, so a run with a source is slow there.
    void (*vorticity_source)(double t, mtrx f, void *data);
    void *source_data;
} solver_config;

// CPU workspace: everything needed to advance a timestep without allocating
// inside the loop. Fields are stored as ny rows (y) of nx values (x): element
// (i, j) is at y = i*dy, x = j*dx.
typedef struct
{
    solver_config cfg;                // copy taken by rk4_alloc(); do not change
    mtrx dwdx, dwdy;                  // first derivatives of w
    mtrx d2wdx2, d2wdy2;              // second derivatives of w
    mtrx psi, psi_scratch;            // Poisson solution and scratch
    mtrx k1, k2, k3, k4;              // RK4 stage increments
    mtrx w_tmp;                       // temporary w for intermediate stages
    mtrx rhs;                         // Poisson right-hand side, -w
    struct fft_solver *fft;           // FFT Poisson solver (poisson_type 3, walls), else NULL
    struct periodic_solver *periodic; // periodic Poisson solver (cfg.periodic), else NULL
    mtrx source;                      // vorticity source of the current stage (cfg.vorticity_source only)
    long steps;                       // steps taken; the time is cfg.t0 + steps * cfg.dt
} rk4_ctx;

// Allocate the CPU workspace for the run described by cfg (copied)
rk4_ctx rk4_alloc(const solver_config *cfg);
// Free the CPU workspace
void rk4_free(rk4_ctx *ctx);

// Evaluate the vorticity RHS dw/dt = -u*dw/dx - v*dw/dy + (1/Re)*(d2w/dx2 + d2w/dy2) + f
// into out, for the interior of w at time t (f only with cfg.vorticity_source). On the way it changes all three inputs:
// u and v become the velocity of w's interior (Poisson solve) with the wall
// velocities imposed, and the boundary entries of w become the wall
// vorticity of that velocity.
void dwdt(mtrx w, mtrx u, mtrx v, mtrx out, double t, rk4_ctx *ctx);

// RK4 time advancement: advances the interior of w by cfg->dt and leaves u, v
// as the velocity of the new interior (without the wall velocities; step()
// completes the state)
void rk4(mtrx w, mtrx u, mtrx v, rk4_ctx *ctx);

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

// Impose the wall velocities on u and v
void apply_wall_bc(mtrx u, mtrx v, const wall_bc *bc);

// One full timestep on the CPU: wall BCs, vorticity BCs, time advancement
// (cfg->time_scheme), Poisson solve and velocity recovery. On return u and v
// are the velocity of the new w (from the stream function, walls included;
// apply_wall_bc() gives the wall velocities), and the boundary of w is the
// wall vorticity of that velocity with the wall velocities imposed, so w is
// consistent with u and v, walls included.
void step(mtrx w, mtrx u, mtrx v, rk4_ctx *ctx);

#endif // FLUIDDYN_H_INCLUDED
