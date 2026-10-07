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
} flow_integrals;

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
// Write output/spectrum-1-<n>.csv (k, E, Z, PE, PZ) for time t
void spectra_write(spectra *s, mtrx u, mtrx v, mtrx w, double t);
void spectra_free(spectra *s);

#endif // DIAGNOSTICS_H_INCLUDED
