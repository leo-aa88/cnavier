// Solver backends: one interface over the CPU solver (fluiddyn.c) and, in
// builds with CUDA=1, the GPU solver (cudasolver.cu)

#ifndef BACKEND_H_INCLUDED
#define BACKEND_H_INCLUDED

#include "linearalg.h"
#include "fluiddyn.h"

typedef struct backend backend;

// Set up a solver for the run described by cfg. With prefer_gpu set, a CUDA
// build uses the GPU if one can be initialised and otherwise says so and
// falls back to the CPU; without CUDA the CPU is always used. Fields start at
// zero. cfg is copied: changing the caller's struct later has no effect. The
// operators it points to must outlive the solver. Exits with an error if the
// chosen backend's host arrays do not fit in the memory available.
backend *backend_create(const solver_config *cfg, int prefer_gpu);

// Host memory, in bytes, of the arrays backend_create() allocates for an
// nx x ny grid: CPU workspace and fields, or on the GPU only host copies of
// the fields. prefer_gpu is ignored in builds without CUDA.
double backend_host_memory(int nx, int ny, int prefer_gpu);
void backend_free(backend *b);

// "CPU" or "CUDA"
const char *backend_name(const backend *b);
// Name of the GPU in use, or NULL on the CPU
const char *backend_device(const backend *b);

// One full timestep
void backend_step(backend *b);

// Max and min of du/dx + dv/dy for the current velocity field
void backend_continuity(backend *b, double *cmax, double *cmin);

// Host copies of the current fields; NULL arguments are skipped. On the CPU
// these are the solver's own arrays, on the GPU arrays filled from the device
// by this call. They stay allocated until backend_free() and are current
// until the next backend_step().
void backend_fields(backend *b, mtrx **u, mtrx **v, mtrx **w);

// Copy fields into and out of the solver. NULL arguments are skipped. The
// arrays from backend_fields() may be passed back to backend_set_fields().
void backend_set_fields(backend *b, const mtrx *u, const mtrx *v, const mtrx *w);
void backend_get_fields(backend *b, mtrx *u, mtrx *v, mtrx *w);

#endif // BACKEND_H_INCLUDED
