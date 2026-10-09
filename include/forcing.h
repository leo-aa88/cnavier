// Built-in forcing and drag: linear drag, Kolmogorov forcing and
// narrow-band random forcing, for the CPU and GPU solvers alike

#ifndef FORCING_H_INCLUDED
#define FORCING_H_INCLUDED

#include "linearalg.h"
#include "fourier.h"

// dw/dt = ... - drag w + f_K - nu_h (-L)^p w - alpha_h psi, and before every
// step a random kick w += dw:
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
    // Small- and large-scale dissipation for turbulence runs:
    //   hyperviscosity    nu_h of -nu_h (-L)^p w, L = DX2 + DY2; 0: none
    //   hyper_order       p >= 2
    //   hypodrag          alpha_h of -alpha_h psi = -alpha_h (-L)^-1 w; 0: none
    double hyperviscosity;
    int hyper_order;
    double hypodrag;
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
    // amp cos(kx x) and amp sin(kx x) of each mode at the nx columns
    // (cx[m nx + j]), and cos(ky y + phase), sin(ky y + phase) at the ny rows
    // (cy[m ny + i]), so that a kick costs no trigonometry per node:
    // amp cos(kx x + ky y + phase) = cx cy - sx sy
    int nx, ny;
    double *cx, *sx, *cy, *sy;
} random_forcing;

// NULL without random forcing. nx, ny, dx, dy: the periodic grid; sym: the
// symbols of the solver's operators (periodic_symbols_of()), which set the
// amplitudes so that each kick carries eps dt of their discrete energy.
random_forcing *random_forcing_setup(const forcing_config *f, int nx, int ny, double dx, double dy, double dt,
                                     const periodic_symbols *sym);
void random_forcing_free(random_forcing *rf);

// Draw the phases of the next kick
void random_forcing_draw(random_forcing *rf);

// w += the current kick, sum over modes of amp cos(kx x + ky y + phase), on
// the grid of random_forcing_setup() (dy: its row spacing)
void random_forcing_add(random_forcing *rf, mtrx w, double dy);

// Initial condition for decaying turbulence on the periodic grid of nx x ny
// nodes spacing dx, dy: random phases, and modal amplitudes whose shell
// envelope is the energy spectrum E(k) ~ (k/k0)^4 exp(-2 (k/k0)^2), which
// peaks at |k| = k0 dk0 (dk0 as for the kicks); the shell sums of the lattice
// modes fluctuate about it with the number of modes per shell. Scaled to the
// energy 1/2 <u^2 + v^2> = energy. w, u and v from
// the stream function with exact (spectral) derivatives; no Nyquist modes.
// The phase of each wavevector depends on the seed and the wavevector only,
// so grids of the same domain that resolve the spectrum get the same field.
void random_initial_field(mtrx w, mtrx u, mtrx v, double dx, double dy, double k0, double energy,
                          unsigned long long seed);

#endif // FORCING_H_INCLUDED
