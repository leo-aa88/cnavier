// Periodic grids: the derivative operators by their symbols, and applied in
// Fourier space (solver_config.fourier)

#ifndef FOURIER_H_INCLUDED
#define FOURIER_H_INCLUDED

#include "linearalg.h"

struct solver_config;

// The symbols of the periodic operators: on exp(i k x) at wavenumber index k
// (0..n-1), the first derivative gives (d1_re + i d1_im) exp(i k x) and the
// second d2 exp(i k x). Centred stencils give d1_re = 0 (to round-off) and
// d2 <= 0. With dealias 1, the 2/3 rule: modes with |k| > n/3 along an axis
// are cut (the masks).
typedef struct
{
    int nx, ny;
    double *d1x_re, *d1x_im, *d2x; // nx values
    double *d1y_re, *d1y_im, *d2y; // ny values
    double *maskx, *masky;         // 1 or 0 (all 1 without dealiasing)
} periodic_symbols;

// The symbols of cfg's operators: from the 1-D operators cfg->D1x, D1y, D2x,
// D2y when set, else from the first rows of the 2-D DX, DY, DX2, DY2. Exits
// if neither is there. Free with periodic_symbols_free().
void periodic_symbols_of(const struct solver_config *cfg, periodic_symbols *s);
void periodic_symbols_free(periodic_symbols *s);

// Workspace for applying the operators in Fourier space (FFTW): every
// operation is exact for the circulant operators, so the results equal those
// of the sparse products to round-off
typedef struct fourier_ops fourier_ops;
fourier_ops *fourier_setup(const struct solver_config *cfg);
void fourier_free(fourier_ops *f);
const periodic_symbols *fourier_symbols(const fourier_ops *f);

// wx = DX w, wy = DY w, lap = (DX2 + DY2) w; NULL outputs are skipped
void fourier_derivatives(fourier_ops *f, const double *w, double *wx, double *wy, double *lap);
// (DX2 + DY2) psi = -w with psi of zero mean, then u = DY psi, v = -DX psi
// (psi, u, v may be NULL to skip)
void fourier_poisson(fourier_ops *f, const double *w, double *psi, double *u, double *v);
// u = DY psi, v = -DX psi
void fourier_velocity(fourier_ops *f, const double *psi, double *u, double *v);
// out = DX a + DY b
void fourier_divergence(fourier_ops *f, const double *a, const double *b, double *out);
// out = (-(DX2 + DY2))^p w
void fourier_power(fourier_ops *f, const double *w, int p, double *out);
// w = w without the modes the dealiasing masks cut (nothing without dealiasing)
void fourier_filter(fourier_ops *f, double *w);
// With 3/2 padding (solver_config.dealias 2): out = -(u DX w + v DY w), the
// velocity from w, with the products formed on a grid 3/2 as fine and
// truncated back, so that every mode of the grid (except the Nyquist modes,
// which are left out) is free of aliasing: the Fourier-Galerkin nonlinear term
void fourier_nonlinear_padded(fourier_ops *f, const double *w, double *out);

#endif // FOURIER_H_INCLUDED
