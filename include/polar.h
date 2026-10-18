// Polar (log-polar) grid around a circular cylinder: the vorticity-
// streamfunction equations in r = R e^xi, theta, with uniform xi and theta.
// The map is conformal, so with h = r the scale factor
//
//   psi_xixi + psi_thth = -r^2 w                       (Poisson)
//   u_r = (1/r) psi_th,  u_th = -(1/r) psi_xi          (velocity)
//   dw/dt = -(u_th w_th + u_r w_xi) / r + (w_xixi + w_thth) / (Re r^2)
//
// Fields keep the solver's layout, ny rows of nx values: column j is
// theta_j = j * dtheta (periodic, nx = ntheta points), row i is xi_i = i * dxi
// (row 0 on the cylinder, row ny - 1 the outer boundary). The solver's u is
// u_theta and its v is u_r, so the right-hand side is the Cartesian one with
// the factors 1/r and 1/r^2 on each row. (theta, xi) is left-handed where
// (r, theta) is right-handed: the velocity is u = -(1/r) DY psi,
// v = (1/r) DX psi, the opposite signs of the Cartesian u = DY psi,
// v = -DX psi.
//
// theta = 0 points downstream (the wake axis); x = r cos theta, y = r sin theta.

#ifndef POLAR_H_INCLUDED
#define POLAR_H_INCLUDED

#include "linearalg.h"

// Outer boundary condition. psi is always Dirichlet there (the Poisson
// solver needs one type for every theta), from psi_out; w is set as follows.
enum polar_outer_bc
{
    POLAR_DIRICHLET,  // psi = psi_out, w = w_out, both set by the caller
    POLAR_FREESTREAM, // psi = r_inf sin(theta); w = 0 on inflow, zero normal gradient on outflow
    POLAR_FARFIELD    // psi from the Schlichting wake and a mass correction (polar_farfield_set()); w as freestream
};

typedef struct polar_grid
{
    int ntheta, nr;        // grid points in theta (nx) and xi (ny)
    double R, r_inf;       // cylinder radius and outer radius
    double dtheta, dxi;    // 2 pi / ntheta and ln(r_inf / R) / (nr - 1)
    double *r;             // r_i = R e^(i dxi), [nr]
    double *inv_r;         // 1 / r_i
    double *inv_r2;        // 1 / r_i^2
    double *cos_t, *sin_t; // cos and sin of theta_j, [ntheta]
    int outer_bc;          // enum polar_outer_bc
    double wall_ut;        // tangential velocity of the cylinder surface (0: at rest)
    double *psi_out;       // psi on the outer boundary, [ntheta]
    double *w_out;         // w on the outer boundary (POLAR_DIRICHLET), [ntheta]
    double cd;             // C_D of the far-field condition (POLAR_FARFIELD)
    double Re;             // Reynolds number of the far-field condition
    // Changed by every change of psi_out, so that a solver does not reuse a
    // stream function solved for the old boundary values
    long version;
} polar_grid;

// Grid of ntheta x nr points between R and r_inf, outer condition freestream.
// nr <= 0: square cells, dxi as close to dtheta as a whole number of cells
// allows. Free with polar_grid_free().
polar_grid polar_grid_make(int ntheta, int nr, double R, double r_inf);
void polar_grid_free(polar_grid *g);
// The number of points in xi for square cells
int polar_square_nr(int ntheta, double R, double r_inf);

// Outer boundary values. polar_freestream_set(): psi_out = r_inf sin(theta),
// w_out = 0. polar_farfield_set(): psi_out = r_inf sin(theta) plus the stream
// function of the far wake of a body with drag coefficient cd at Reynolds
// number Re: the Schlichting deficit
//   d = cd sqrt(Re) / (4 sqrt(pi)) exp(-Re r (1 - cos theta) / 2) / sqrt(r)
// in u_x, with u_y = -(sin theta / 2) d (Gautier et al. 2013), and a uniform
// radial outflow that returns the flux the deficit removes (Hasan et al.
// 2005), so that psi_out is single-valued. Both bump g->version; call them
// between steps only.
void polar_freestream_set(polar_grid *g);
void polar_farfield_set(polar_grid *g, double cd, double Re);

// Exact solution with w = 0: potential flow past the cylinder inside the
// circle r_inf with psi = r_inf sin(theta) there, psi = A (r - R^2/r) sin(theta)
// with A = r_inf / (r_inf - R^2 / r_inf)
double polar_potential_psi(const polar_grid *g, double r, double theta);

// Poisson solver: psi_xixi + psi_thth = f, psi = 0 on row 0 and psi = psi_out
// on row nr - 1, by a real FFT in theta and a tridiagonal solve in xi for each
// mode. In theta the operator is DX2 (its symbols, from the first row of the
// 2-D operator), in xi the 3-point one (order 2) or the compact fourth-order
// one, delta^2 / (1 + dxi^2 delta^2 / 12) (order 4). Order 4 reads f on the
// boundary rows; they are extrapolated from the interior, so f is read on the
// interior rows only.
typedef struct polar_poisson polar_poisson;
polar_poisson *polar_poisson_setup(const polar_grid *g, const smtrx *DX2, int order);
void polar_poisson_free(polar_poisson *s);
// psi for w (nr x ntheta): f = -r^2 w on the interior rows
void polar_poisson_solve(polar_poisson *s, const double *w, const double *psi_out, double *psi);
// The same with f given (interior rows read)
void polar_poisson_solve_f(polar_poisson *s, const double *f, const double *psi_out, double *psi);

// u = u_theta, v = u_r from psi with the first-derivative operators
void polar_velocity(const polar_grid *g, const smtrx *DX, const smtrx *DY, const double *psi, double *u, double *v);

// Velocity on the cylinder (u = wall_ut, v = 0 on row 0)
void polar_apply_bc(const polar_grid *g, mtrx u, mtrx v);

// Wall vorticity on the cylinder (row 0) and w on the outer boundary (row
// nr - 1). closure 1: third order from psi (Briley), w = -(1/R^2) psi_xixi
// with psi_xi = -R wall_ut; closure 0: w = (1/r^2) (DY(r u) - DX(r v)) with
// the one-sided rows of DY. The outer row follows g->outer_bc.
void polar_boundary_vorticity(const polar_grid *g, int closure, const smtrx *DX, const smtrx *DY, mtrx w,
                              const mtrx *psi, mtrx u, mtrx v);

// du_r/dr + u_r/r + (1/r) du_th/dth = (1/r^2) (DY(r v) + DX(r u)) into div;
// s1, s2: scratch, all ntheta * nr values
void polar_divergence(const polar_grid *g, const smtrx *DX, const smtrx *DY, const double *u, const double *v,
                      double *div, double *s1, double *s2);

// ---------------------------------------------------------------------------
// Diagnostics of the flow past the cylinder (D = 2R = 1, U = 1 in the units
// the drag coefficient and lengths are given in)
// ---------------------------------------------------------------------------

typedef struct
{
    double cd, cdp, cdf; // drag coefficient: total, pressure and friction parts
    double cl;           // lift coefficient
    double dpdt_mean;    // mean of dp/dtheta on the wall; 0 for a consistent pressure
} polar_forces;

// Forces from the wall vorticity: the friction from w on row 0, the pressure
// from dp/dtheta = (1/Re) dw/dxi on the wall (theta-momentum on a no-slip
// wall), integrated in theta by FFT. DY gives dw/dxi (its one-sided row 0).
// p_wall (ntheta values, or NULL) receives the wall pressure with zero mean.
// Lengths are made dimensionless with D = 2R and the speed with U = 1.
polar_forces polar_wall_forces(const polar_grid *g, const smtrx *DY, const double *w, double Re, double *p_wall);

// Separation angle in degrees from the front stagnation point (180 - the
// angle from the rear where w on the wall changes sign on 0 < theta < pi), or
// -1 if the wall vorticity does not change sign (no separation)
double polar_separation_angle(const polar_grid *g, const double *w);

// Length of the recirculation bubble from the rear stagnation point: where
// u_r (v) on theta = 0 turns from negative to positive, minus R. 0 if there is
// no reverse flow.
double polar_wake_length(const polar_grid *g, const double *v);

// Centre of the upper recirculation bubble: the minimum of psi on
// 0 < theta < pi/2 inside the bubble, refined by a quadratic fit. a: its
// distance downstream of the rear stagnation point; b: the distance between
// the two centres (2 y). Returns 0 if psi has no interior minimum below 0.
int polar_vortex_centre(const polar_grid *g, const double *psi, double *a, double *b);

// Fields at a point (x, y) inside the grid by Lagrange interpolation of degree
// 3 in xi and theta: the Cartesian velocity and the vorticity; any output may
// be NULL. Returns 0 if (x, y) is outside the grid.
int polar_interpolate(const polar_grid *g, const double *u, const double *v, const double *w, double x, double y,
                      double *ux, double *uy, double *wv);

// max |w(theta) + w(-theta)| over the grid: the flow is symmetric about the
// wake axis when w is odd in theta
double polar_asymmetry(const polar_grid *g, const double *w);

// Write a VTK structured grid output/<title>-1-<n>.vtk of the field f in the
// Cartesian plane, with column 0 repeated at theta = 2 pi so that the annulus
// is closed; with u and v (both or neither) also the Cartesian velocity
void printvtk_polar(const polar_grid *g, const mtrx *f, const mtrx *u, const mtrx *v, const char *title);

#endif // POLAR_H_INCLUDED
