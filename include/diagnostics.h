// Flow diagnostics: domain integrals, and on periodic grids spectra and
// spectral fluxes

#ifndef DIAGNOSTICS_H_INCLUDED
#define DIAGNOSTICS_H_INCLUDED

#include "linearalg.h"
#include "fluiddyn.h"

// Domain means of the kinetic energy E = 1/2 (u^2 + v^2), the enstrophy
// Z = 1/2 w^2 and the palinstrophy P = 1/2 |grad w|^2 (grad by DX, DY). On a
// periodic grid the mean is over the nodes. With walls it is the trapezoidal
// mean over the domain, and u, v on the wall nodes are the wall velocities.
typedef struct
{
    double E, Z, P;
    // Energy input of the built-in forcing (forcing.h): the work <u A sin(k y)>
    // of the Kolmogorov force (the continuum expression), plus eps, the
    // expected input of the random kicks. A single kick also changes E by
    // <u . du>, which averages to zero but not on any one step, so a
    // step-by-step energy budget holds only on average.
    double I;
    // The same with the Kolmogorov part replaced by the exact rate at which
    // the source f_K changes the discrete energy, (k k1 / Q) <u A sin(k y)>,
    // with k1 and Q the symbols of DY and -DY2 at k (kolmogorov_factor());
    // on a periodic grid the discrete energy budget closes with it exactly
    double I_disc;
} flow_integrals;

// k k1 / Q for the Kolmogorov mode: the ratio of the exact rate at which f_K
// changes the discrete energy to the work <u A sin(k y)>. f_K is a single
// mode, so psi_f = (A k / Q) cos(k y) up to sign and DY gives it a velocity
// (A k k1 / Q) sin(k y), with k1 the first-derivative symbol. 1 + O(h^p) for a
// resolved mode; 1 with walls, where it is not defined.
double kolmogorov_factor(const solver_config *cfg);

// wx, wy: scratch of nx*ny values each
flow_integrals compute_integrals(const solver_config *cfg, mtrx u, mtrx v, mtrx w, double *wx, double *wy);

// Spectra on a periodic grid, summed over shells of |k| (not averaged): bins
// of width dk = 2 pi / max(Lx, Ly), bin b holding b dk - dk/2 <= |k| < b dk + dk/2:
//   E(k)  energy, 1/2 (|u^|^2 + |v^|^2), summing to E over the bins
//   Z(k)  enstrophy, 1/2 |w^|^2, summing to Z
//   PE(k) energy flux through k, -sum over bins <= k of (A/Q) Re(conj(psi^) N^)
//   PZ(k) enstrophy flux through k, -sum over bins <= k of Re(conj(w^) N^)
// with N = -(u DX w + v DY w) the solver's nonlinear term, psi the solution of
// (DX2 + DY2) psi = -w, A = |symbol of DX|^2 + |symbol of DY|^2 and Q = -(symbol
// of DX2 + DY2). A mode's energy is 1/2 A |psi^|^2 and w^ = Q psi^, so these
// are exactly the rates at which N changes E and Z: summed over all shells
// they equal 1/2 d<u^2 + v^2>/dt and 1/2 d<w^2>/dt from N, to round-off. A
// positive flux carries energy or enstrophy to larger k. For a spectral
// density, divide E(k) and Z(k) by dk.
typedef struct spectra spectra;

spectra *spectra_setup(const solver_config *cfg);
int spectra_bins(const spectra *s);
double spectra_dk(const spectra *s);
// Arrays of spectra_bins() values; any may be NULL
void spectra_compute(spectra *s, mtrx u, mtrx v, mtrx w, double *E, double *Z, double *PE, double *PZ);
// Dissipation of energy and enstrophy in each shell, as the discrete
// equations have it. A damping term of symbol -sigma removes
// sigma (A/Q^2) |w^|^2 of energy and sigma |w^|^2 of enstrophy per mode:
//   DE, DZ  small scales: viscosity nu (DX2 + DY2) w, sigma = nu Q, plus the
//           hyperviscosity, nu_h Q^p. For viscosity alone,
//           DE(k) = nu sum over the shell of (A/Q) |w^|^2   (continuum: 2 nu Z(k))
//           DZ(k) = nu sum over the shell of Q |w^|^2       (continuum: 2 nu P(k))
//   FE, FZ  large scales: drag, sigma = alpha (FE = 2 alpha E(k)), plus the
//           hypodrag, alpha_h / Q
// so dE/dt = sum of the nonlinear transfers - sum (DE + FE) + the input of
// the forcing (flow_integrals.I_disc) exactly, where 2 nu Z holds only to the
// order of the scheme; the difference grows where A/Q departs from 1, near the
// grid cutoff. Arrays of spectra_bins() values; any may be NULL.
void spectra_dissipation(spectra *s, mtrx w, double *DE, double *DZ, double *FE, double *FZ);
// Write output/spectrum-1-<n>.csv (k, E, Z, PE, PZ, DE, DZ, FE, FZ) for time t
void spectra_write(spectra *s, mtrx u, mtrx v, mtrx w, double t);
void spectra_free(spectra *s);

#endif // DIAGNOSTICS_H_INCLUDED
