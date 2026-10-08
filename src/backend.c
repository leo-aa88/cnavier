#include <stdio.h>
#include <stdlib.h>
#include "linearalg.h"
#include "fluiddyn.h"
#include "poisson.h"
#include "backend.h"
#include "utils.h"
#include "diagnostics.h"
#ifdef USE_CUDA
#include "cudasolver.h"
#endif

struct backend
{
    solver_config cfg;

#ifdef USE_CUDA
    gpu_solver *gpu; // NULL when running on the CPU
#endif

    // Fields on the host: the CPU solver's own, or the GPU solver's host copies
    mtrx u, v, w;

    // CPU solver only: workspace, and scratch for the continuity check
    rk4_ctx ws;
    mtrx dudx, dvdy;
};

static int on_gpu(const backend *b)
{
#ifdef USE_CUDA
    return b->gpu != NULL;
#else
    (void)b;
    return 0;
#endif
}

// Arrays of nx*ny doubles: CPU workspace (14 and the FFT buffer) plus fields
// (3) and continuity scratch (2); on the GPU, host copies of the fields (3)
#define CPU_ARRAYS 20
#define GPU_ARRAYS 3

double backend_host_memory(int nx, int ny, int prefer_gpu)
{
#ifndef USE_CUDA
    prefer_gpu = 0;
#endif
    return (double)nx * ny * sizeof(double) * (prefer_gpu ? GPU_ARRAYS : CPU_ARRAYS);
}

// Stop with a message if `bytes` more cannot be allocated now
static void require_memory(double bytes, const char *what)
{
    double avail = available_memory();
    if (avail >= 0 && bytes > avail) // -1: unknown
    {
        printf("** Error: the %s needs about %.1f GB more memory; about %.1f GB is available **\n",
               what, bytes / 1E9, avail / 1E9);
        exit(1);
    }
}

backend *backend_create(const solver_config *cfg, int prefer_gpu)
{
    backend *b = (backend *)calloc(1, sizeof(backend));
    if (!b)
    {
        printf("** Error: insufficient memory **\n");
        exit(1);
    }
    b->cfg = *cfg;

#ifdef USE_CUDA
    if (prefer_gpu)
    {
        b->gpu = gpu_init(&b->cfg);
        if (b->gpu)
        {
            require_memory(backend_host_memory(cfg->nx, cfg->ny, 1), "GPU backend");
            b->u = initm(cfg->ny, cfg->nx);
            b->v = initm(cfg->ny, cfg->nx);
            b->w = initm(cfg->ny, cfg->nx);
            return b;
        }
        printf("No usable CUDA device - falling back to the CPU\n");
    }
#else
    (void)prefer_gpu;
#endif

    require_memory(backend_host_memory(cfg->nx, cfg->ny, 0), "CPU backend");
    b->ws = rk4_alloc(&b->cfg);
    b->u = initm(cfg->ny, cfg->nx);
    b->v = initm(cfg->ny, cfg->nx);
    b->w = initm(cfg->ny, cfg->nx);
    b->dudx = initm(cfg->ny, cfg->nx);
    b->dvdy = initm(cfg->ny, cfg->nx);
    return b;
}

void backend_free(backend *b)
{
    if (!b) return;
#ifdef USE_CUDA
    if (b->gpu)
    {
        gpu_free(b->gpu);
        freem(&b->u);
        freem(&b->v);
        freem(&b->w);
        free(b);
        return;
    }
#endif
    rk4_free(&b->ws);
    freem(&b->u);
    freem(&b->v);
    freem(&b->w);
    freem(&b->dudx);
    freem(&b->dvdy);
    free(b);
}

const char *backend_name(const backend *b)
{
    return on_gpu(b) ? "CUDA" : "CPU";
}

const char *backend_device(const backend *b)
{
#ifdef USE_CUDA
    if (b->gpu) return gpu_device_name();
#endif
    (void)b;
    return NULL;
}

void backend_step(backend *b)
{
#ifdef USE_CUDA
    if (b->gpu)
    {
        gpu_step(b->gpu);
        return;
    }
#endif
    step(b->w, b->u, b->v, &b->ws);
}

void backend_continuity(backend *b, double *cmax, double *cmin)
{
    int k, n = b->cfg.nx * b->cfg.ny;

#ifdef USE_CUDA
    if (b->gpu)
    {
        gpu_continuity(b->gpu, cmax, cmin);
        return;
    }
#endif
    spmv(b->cfg.DXv ? *b->cfg.DXv : *b->cfg.DX, b->u.M, b->dudx.M);
    spmv(b->cfg.DYv ? *b->cfg.DYv : *b->cfg.DY, b->v.M, b->dvdy.M);
    *cmax = -__DBL_MAX__;
    *cmin = __DBL_MAX__;
    for (k = 0; k < n; k++)
    {
        double c = b->dudx.M[k] + b->dvdy.M[k];
        if (c > *cmax) *cmax = c;
        if (c < *cmin) *cmin = c;
    }
}

flow_integrals backend_integrals(backend *b)
{
#ifdef USE_CUDA
    if (b->gpu)
    {
        flow_integrals r;
        gpu_integrals(b->gpu, &r.E, &r.Z, &r.P, &r.I);
        r.I_disc = kolmogorov_factor(&b->cfg) * (r.I - b->cfg.forcing.random_rate) + b->cfg.forcing.random_rate;
        return r;
    }
#endif
    return compute_integrals(&b->cfg, b->u, b->v, b->w, b->dudx.M, b->dvdy.M);
}

void backend_fields(backend *b, mtrx **u, mtrx **v, mtrx **w)
{
#ifdef USE_CUDA
    if (b->gpu)
        gpu_get_fields(b->gpu, u ? &b->u : NULL, v ? &b->v : NULL, w ? &b->w : NULL);
#endif
    if (u) *u = &b->u;
    if (v) *v = &b->v;
    if (w) *w = &b->w;
}

void backend_set_fields(backend *b, const mtrx *u, const mtrx *v, const mtrx *w)
{
#ifdef USE_CUDA
    if (b->gpu)
    {
        gpu_set_fields(b->gpu, u, v, w);
        return;
    }
#endif
    // Arrays from backend_fields() are already the solver's
    if (u && u->M != b->u.M) mtrxcpy(b->u, *u);
    if (v && v->M != b->v.M) mtrxcpy(b->v, *v);
    if (w && w->M != b->w.M) mtrxcpy(b->w, *w);
}

void backend_get_fields(backend *b, mtrx *u, mtrx *v, mtrx *w)
{
#ifdef USE_CUDA
    if (b->gpu)
    {
        gpu_get_fields(b->gpu, u, v, w);
        return;
    }
#endif
    if (u) mtrxcpy(*u, b->u);
    if (v) mtrxcpy(*v, b->v);
    if (w) mtrxcpy(*w, b->w);
}
