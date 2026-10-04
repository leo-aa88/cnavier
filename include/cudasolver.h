// CUDA backend (optional, built with `make CUDA=1`)

#ifndef CUDASOLVER_H_INCLUDED
#define CUDASOLVER_H_INCLUDED

#include "linearalg.h"
#include "fluiddyn.h"

#ifdef __cplusplus
extern "C" {
#endif

// Opaque solver state: all fields, operators and workspace live on the device.
typedef struct gpu_solver gpu_solver;

// Upload the operators and settings held by ctx and allocate the device workspace.
// Fields start at zero.
// Returns NULL if no CUDA device can be initialised, so the caller can fall
// back to the CPU path. Any failure after that point (out of device memory, a
// non-square grid, an unknown Poisson solver type) is fatal: it prints an
// error and exits.
gpu_solver *gpu_init(const rk4_ctx *ctx, double dt, int time_scheme, const wall_bc *bc);
void gpu_free(gpu_solver *g);

// Name of the device used by gpu_init
const char *gpu_device_name(void);

// One full timestep on the device — same sequence as step() in fluiddyn.h
void gpu_step(gpu_solver *g);

// Max and min of du/dx + dv/dy for the current velocity field
void gpu_continuity(gpu_solver *g, double *cmax, double *cmin);

// Copy fields between host and device. NULL arguments are skipped.
void gpu_set_fields(gpu_solver *g, const mtrx *u, const mtrx *v, const mtrx *w);
void gpu_get_fields(gpu_solver *g, mtrx *u, mtrx *v, mtrx *w);

// Building blocks on host arrays of size nx*ny, exposed for the tests.
// op: 0=DX, 1=DY, 2=DX2, 3=DY2
void gpu_spmv(gpu_solver *g, int op, const double *x, double *y);
// Solve nabla^2 psi = -w with the configured Poisson solver
void gpu_poisson(gpu_solver *g, const double *w, double *psi);

#ifdef __cplusplus
}
#endif

#endif // CUDASOLVER_H_INCLUDED
