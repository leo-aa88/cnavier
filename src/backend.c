#include <stdio.h>
#include <stdlib.h>
#include "linearalg.h"
#include "fluiddyn.h"
#include "poisson.h"
#include "backend.h"
#ifdef USE_CUDA
#include "cudasolver.h"
#endif

struct backend
{
    solver_config cfg;

#ifdef USE_CUDA
    gpu_solver *gpu; // NULL when running on the CPU
#endif

    // CPU solver: fields, workspace, and scratch for the continuity check
    rk4_ctx ws;
    mtrx    u, v, w;
    mtrx    dudx, dvdy;
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
            return b;
        printf("No usable CUDA device - falling back to the CPU\n");
    }
#else
    (void)prefer_gpu;
#endif

    b->ws   = rk4_alloc(&b->cfg);
    b->u    = initm(cfg->ny, cfg->nx);
    b->v    = initm(cfg->ny, cfg->nx);
    b->w    = initm(cfg->ny, cfg->nx);
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
        free(b);
        return;
    }
#endif
    rk4_free(&b->ws);
    freem(&b->u); freem(&b->v); freem(&b->w);
    freem(&b->dudx); freem(&b->dvdy);
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
    spmv(*b->cfg.DX, b->u.M, b->dudx.M);
    spmv(*b->cfg.DY, b->v.M, b->dvdy.M);
    *cmax = -__DBL_MAX__;
    *cmin = __DBL_MAX__;
    for (k = 0; k < n; k++)
    {
        double c = b->dudx.M[k] + b->dvdy.M[k];
        if (c > *cmax) *cmax = c;
        if (c < *cmin) *cmin = c;
    }
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
    if (u) mtrxcpy(b->u, *u);
    if (v) mtrxcpy(b->v, *v);
    if (w) mtrxcpy(b->w, *w);
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
