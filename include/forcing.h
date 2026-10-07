// Built-in forcing and drag: linear drag, Kolmogorov forcing and
// narrow-band random forcing, for the CPU and GPU solvers alike

#ifndef FORCING_H_INCLUDED
#define FORCING_H_INCLUDED

#include "linearalg.h"

// dw/dt = ... - drag w + f_K, and before every step a random kick w += dw:
//   drag              linear (Ekman) drag alpha; 0: none
//   kolmogorov_amp    A of the body force A sin(k y) in x, whose curl is the
//   kolmogorov_n      vorticity source f_K = -A k cos(k y), k = 2 pi n / Ly; 0: none
//   random_rate       energy injected per unit time by the kicks, eps; 0: none
//   random_kf/dk      the kicks act on the wavevectors with
//                     |kf - |k|/dk0| <= dk, dk0 = 2 pi / max(Lx, Ly)
//   random_seed       seed of the kicks' random phases
// The random kicks need a periodic grid. Each kick has random phases and
// injects exactly eps dt of the solver's discrete kinetic energy
// 1/2 <u^2 + v^2>; with independent phases from step to step this is
// white-in-time forcing at rate eps.
typedef struct
{
    double drag;
    double kolmogorov_amp;
    int kolmogorov_n;
    double random_rate;
    double random_kf, random_dk;
    unsigned long long random_seed;
} forcing_config;

// The Kolmogorov source f_K = -A k cos(k y) at the ny rows y = i dy, with
// k = 2 pi n / Ly; NULL without Kolmogorov forcing
double *kolmogorov_rows(const forcing_config *f, int ny, double dy, double Ly);

// The random kicks: their modes, fixed amplitudes and the current phases
typedef struct
{
    int modes;                // number of wavevectors in the shell (half plane)
    double *kx, *ky;          // wavevectors
    double *amp;              // amplitude of each mode's vorticity cosine
    double *phase;            // phases of the current kick
    unsigned long long state; // random number generator
} random_forcing;

// NULL without random forcing. nx, ny, dx, dy: the periodic grid; DX, DY,
// DX2, DY2: the solver's operators, whose symbols set the amplitudes so that
// each kick carries eps dt of their discrete energy.
random_forcing *random_forcing_setup(const forcing_config *f, int nx, int ny, double dx, double dy, double dt,
                                     const smtrx *DX, const smtrx *DY, const smtrx *DX2, const smtrx *DY2);
void random_forcing_free(random_forcing *rf);

// Draw the phases of the next kick
void random_forcing_draw(random_forcing *rf);

// w += the current kick, sum over modes of amp cos(kx x + ky y + phase)
void random_forcing_add(const random_forcing *rf, mtrx w, double dx, double dy);

#endif // FORCING_H_INCLUDED
