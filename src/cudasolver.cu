// CUDA backend.
//
// Every field (u, v, w, psi, RK4 stages) and the four CSR derivative operators
// live on the device, and a timestep is a fixed sequence of kernels that
// mirrors step() in fluiddyn.c. Data returns to the host only for output and
// for the scalar diagnostics.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <cuda_runtime.h>
#include <cufft.h>

extern "C"
{
#include "linearalg.h"
#include "fluiddyn.h"
#include "poisson.h"
}
#include "cudasolver.h"

#define BLOCK      256 // threads per block (power of two, required by reduce_kernel)
#define RED_BLOCKS 256 // max blocks used by a reduction
#define RB_BATCH   32  // red-black sweeps queued between looks at the convergence state

#define CUDA_CHECK(call)                                          \
    do                                                            \
    {                                                             \
        cudaError_t err_ = (call);                                \
        if (err_ != cudaSuccess)                                  \
        {                                                         \
            printf("** CUDA error: %s (%s:%d) **\n",              \
                   cudaGetErrorString(err_), __FILE__, __LINE__); \
            exit(1);                                              \
        }                                                         \
    } while (0)

#define CUFFT_CHECK(call)                             \
    do                                                \
    {                                                 \
        cufftResult res_ = (call);                    \
        if (res_ != CUFFT_SUCCESS)                    \
        {                                             \
            printf("** cuFFT error: %d (%s:%d) **\n", \
                   (int)res_, __FILE__, __LINE__);    \
            exit(1);                                  \
        }                                             \
    } while (0)

// Launch a kernel with one thread per element of an n-element array
#define LAUNCH(kernel, n, ...)                                     \
    do                                                             \
    {                                                              \
        kernel<<<((n) + BLOCK - 1) / BLOCK, BLOCK>>>(__VA_ARGS__); \
        CUDA_CHECK(cudaGetLastError());                            \
    } while (0)

// Device copy of a CSR matrix
typedef struct
{
    double *values;
    int *col_idx;
    int *row_ptr;
} csr_dev;

struct gpu_solver
{
    solver_config cfg; // copy taken by gpu_init(); later changes to the caller's have no effect
    int nx, ny, n;     // cfg.nx, cfg.ny and their product, for brevity

    csr_dev DX, DY, DX2, DY2;

    double *u, *v, *w, *psi;
    double *k1, *k2, *k3, *k4; // RK4 stage increments
    double *w_tmp;             // temporary w for intermediate stages
    double *scratch;           // continuity field / Poisson work array
    double *partial;           // per-block results of a reduction

    // Convergence state of the iterative Poisson solvers, kept on the device
    // so that a batch of sweeps runs without waiting for the host
    struct rb_state *rb;

    // FFT Poisson solver (poisson_type 3)
    cufftHandle plan_rows;       // batched real FFTs of the odd extensions of the rows
    cufftHandle plan_cols;       // ... and of the columns
    double *ext;                 // odd extensions, (ny-2) x 2(nx-1) or (nx-2) x 2(ny-1)
    cufftDoubleComplex *spec;    // their spectra, (ny-2) x nx or (nx-2) x ny
    double *lambda_i, *lambda_j; // eigenvalues of the 1D second differences
};

// ---------------------------------------------------------------------------
// Kernels
// ---------------------------------------------------------------------------

// Row r of A*x
__device__ inline double csr_row(csr_dev A, const double *x, int r)
{
    double sum = 0.0;
    for (int k = A.row_ptr[r]; k < A.row_ptr[r + 1]; k++)
        sum += A.values[k] * x[A.col_idx[k]];
    return sum;
}

// y = A*x
__global__ void spmv_kernel(csr_dev A, const double *x, double *y, int n)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < n)
        y[k] = csr_row(A, x, k);
}

// Dirichlet wall velocities. The j walls win at the corners, as on the CPU.
__global__ void wall_bc_kernel(double *u, double *v, wall_bc bc, int nx, int ny)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= nx * ny) return;
    int i = k / nx, j = k % nx;
    int wall;

    if (j == 0)
        wall = 0;
    else if (j == nx - 1)
        wall = 1;
    else if (i == 0)
        wall = 2;
    else if (i == ny - 1)
        wall = 3;
    else
        return;

    u[k] = bc.u[wall];
    v[k] = bc.v[wall];
}

// Vorticity BCs: w = dv/dx - du/dy evaluated at boundaries
__global__ void vorticity_bc_kernel(csr_dev DX, csr_dev DY, const double *u, const double *v,
                                    double *w, int nx, int ny)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= nx * ny) return;
    int i = k / nx, j = k % nx;

    if (i == 0 || i == ny - 1 || j == 0 || j == nx - 1)
        w[k] = csr_row(DX, v, k) - csr_row(DY, u, k);
}

// out = -u*(dw/dx) - v*(dw/dy) + (1/Re)*(d2w/dx2 + d2w/dy2)
__global__ void rhs_kernel(csr_dev DX, csr_dev DY, csr_dev DX2, csr_dev DY2,
                           const double *w, const double *u, const double *v,
                           double Re, double *out, int n)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < n)
        out[k] = -u[k] * csr_row(DX, w, k) - v[k] * csr_row(DY, w, k) + (1.0 / Re) * (csr_row(DX2, w, k) + csr_row(DY2, w, k));
}

// u = dpsi/dy, v = -dpsi/dx
__global__ void velocity_kernel(csr_dev DX, csr_dev DY, const double *psi,
                                double *u, double *v, int n)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < n)
    {
        u[k] = csr_row(DY, psi, k);
        v[k] = -csr_row(DX, psi, k);
    }
}

// out = du/dx + dv/dy
__global__ void continuity_kernel(csr_dev DX, csr_dev DY, const double *u, const double *v,
                                  double *out, int n)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < n)
        out[k] = csr_row(DX, u, k) + csr_row(DY, v, k);
}

// out = x + a*y  (out may be x)
__global__ void axpy_kernel(const double *x, double a, const double *y, double *out, int n)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < n)
        out[k] = x[k] + a * y[k];
}

// w += c*(k1 + 2*k2 + 2*k3 + k4)
__global__ void rk4_combine_kernel(double *w, double c, const double *k1, const double *k2,
                                   const double *k3, const double *k4, int n)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < n)
        w[k] += c * (k1[k] + 2.0 * k2[k] + 2.0 * k3[k] + k4[k]);
}

// Device-side state of a red-black solve
struct rb_state
{
    int done;   // set once the stopping rule is met; later sweeps do nothing
    int iter;   // sweep at which it was met
    double err; // sum of |change| in that sweep
};

// One colour of a red-black sweep for nabla^2 psi = -w; beta = 1 is Gauss-Seidel.
// Points of one colour only read the other colour, so the update is safe in
// parallel. delta receives |change| at every updated point.
__global__ void redblack_kernel(double *psi, const double *w, double *delta, int nx, int ny,
                                double dx2, double dy2, double beta, int colour,
                                const struct rb_state *rb)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= nx * ny || rb->done) return;
    int i = k / nx, j = k % nx;

    if (i < 1 || i >= ny - 1 || j < 1 || j >= nx - 1 || ((i + j) & 1) != colour)
        return;

    double denom = 2.0 * (dx2 + dy2);
    double old = psi[k];
    double upd = beta * (dx2 * (psi[k + nx] + psi[k - nx]) // y-neighbours
                         + dy2 * (psi[k + 1] + psi[k - 1]) // x-neighbours
                         + dx2 * dy2 * w[k]) /
                     denom +
                 (1.0 - beta) * old;
    psi[k] = upd;
    delta[k] = fabs(upd - old);
}

// The 2D DST-I of the interior — the transform FFTW calls RODFT00, which cuFFT
// does not provide — is done as two batches of 1D transforms, first along the
// rows (x), then along the columns (y). The DST-I of x[1..m] is -Im of the
// real FFT of its odd extension [0, x1..xm, 0, -xm..-x1] of length 2(m+1).

// Odd extension of every interior row: row r is the extension of
// src(r+1, 1..nx-2), of length 2(nx-1)
__global__ void extend_rows_kernel(const double *src, double *ext, int nx, int ny)
{
    int len = 2 * (nx - 1);
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= (ny - 2) * len) return;
    int r = k / len, q = k % len;

    if (q == 0 || q == nx - 1)
        ext[k] = 0.0;
    else
        ext[k] = q < nx - 1 ? src[(r + 1) * nx + q] : -src[(r + 1) * nx + (len - q)];
}

// Odd extension of every column of the row transforms: column c is the
// extension of the DST along x of rows 0..ny-3 at mode c, of length 2(ny-1).
// Mode c of row r is -Im(spec_rows[r][c + 1]).
__global__ void extend_cols_kernel(const cufftDoubleComplex *spec_rows, double *ext, int nx, int ny)
{
    int len = 2 * (ny - 1);
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= (nx - 2) * len) return;
    int c = k / len, p = k % len;

    if (p == 0 || p == ny - 1)
    {
        ext[k] = 0.0;
        return;
    }
    int r = p < ny - 1 ? p - 1 : len - p - 1;
    double mode = -spec_rows[r * nx + c + 1].y;
    ext[k] = p < ny - 1 ? mode : -mode;
}

// Interior node (i, j) corresponds to sine mode (i-1, j-1): mode i-1 of the
// column transform of column j-1. The column spectra have ny entries each.
__device__ inline double dst_mode(const cufftDoubleComplex *spec_cols, int i, int j, int ny)
{
    return -spec_cols[(j - 1) * ny + i].y;
}

// Divide each sine mode by its eigenvalue of the 2D Laplacian. Wall nodes are
// set to 0.
__global__ void spectral_divide_kernel(const cufftDoubleComplex *spec, double *out,
                                       const double *lambda_i, const double *lambda_j,
                                       int nx, int ny)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= nx * ny) return;
    int i = k / nx, j = k % nx;

    if (i == 0 || i == ny - 1 || j == 0 || j == nx - 1)
        out[k] = 0.0;
    else
        out[k] = dst_mode(spec, i, j, ny) / (lambda_i[i - 1] + lambda_j[j - 1]);
}

// Scale each sine mode. Wall nodes are set to 0.
__global__ void spectral_scale_kernel(const cufftDoubleComplex *spec, double *out,
                                      double scale, int nx, int ny)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= nx * ny) return;
    int i = k / nx, j = k % nx;

    if (i == 0 || i == ny - 1 || j == 0 || j == nx - 1)
        out[k] = 0.0;
    else
        out[k] = dst_mode(spec, i, j, ny) * scale;
}

enum
{
    RED_SUM,
    RED_MAX,
    RED_MIN
};

__host__ __device__ inline double red_identity(int op)
{
    return op == RED_SUM ? 0.0 : (op == RED_MAX ? -DBL_MAX : DBL_MAX);
}

__host__ __device__ inline double red_combine(double a, double b, int op)
{
    return op == RED_SUM ? a + b : (op == RED_MAX ? fmax(a, b) : fmin(a, b));
}

// Reduce x[0..n-1] to one value per block; the host combines the blocks.
__global__ void reduce_kernel(const double *x, int n, int op, double *partial)
{
    __shared__ double s[BLOCK];
    double acc = red_identity(op);

    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < n; k += gridDim.x * blockDim.x)
        acc = red_combine(acc, x[k], op);
    s[threadIdx.x] = acc;
    __syncthreads();

    for (int half = blockDim.x / 2; half > 0; half >>= 1)
    {
        if (threadIdx.x < half)
            s[threadIdx.x] = red_combine(s[threadIdx.x], s[threadIdx.x + half], op);
        __syncthreads();
    }
    if (threadIdx.x == 0)
        partial[blockIdx.x] = s[0];
}

// Finish the sum of |change| for sweep `iter` and apply the stopping rule of
// the CPU solvers. The block results are added in order by one thread, as
// reduce() does on the host, so the sum and the stopping sweep are the same.
__global__ void rb_check_kernel(const double *partial, int blocks, double tol, int iter,
                                struct rb_state *rb)
{
    if (blockIdx.x != 0 || threadIdx.x != 0 || rb->done) return;
    double acc = red_identity(RED_SUM);
    for (int b = 0; b < blocks; b++)
        acc = red_combine(acc, partial[b], RED_SUM);
    if (acc < tol)
    {
        rb->done = 1;
        rb->iter = iter;
        rb->err = acc;
    }
}

// ---------------------------------------------------------------------------
// Host-side helpers
// ---------------------------------------------------------------------------

static double *dev_alloc(size_t count)
{
    double *p;
    CUDA_CHECK(cudaMalloc((void **)&p, count * sizeof(double)));
    CUDA_CHECK(cudaMemset(p, 0, count * sizeof(double)));
    return p;
}

static csr_dev csr_upload(const smtrx *A)
{
    csr_dev d;
    int nnz = A->row_ptr[A->m];

    CUDA_CHECK(cudaMalloc((void **)&d.values, nnz * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void **)&d.col_idx, nnz * sizeof(int)));
    CUDA_CHECK(cudaMalloc((void **)&d.row_ptr, (A->m + 1) * sizeof(int)));
    CUDA_CHECK(cudaMemcpy(d.values, A->values, nnz * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d.col_idx, A->col_idx, nnz * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d.row_ptr, A->row_ptr, (A->m + 1) * sizeof(int), cudaMemcpyHostToDevice));
    return d;
}

static void csr_free(csr_dev d)
{
    cudaFree(d.values);
    cudaFree(d.col_idx);
    cudaFree(d.row_ptr);
}

static double reduce(gpu_solver *g, const double *x, int n, int op)
{
    double host[RED_BLOCKS];
    int b, blocks = (n + BLOCK - 1) / BLOCK;
    if (blocks > RED_BLOCKS) blocks = RED_BLOCKS;

    reduce_kernel<<<blocks, BLOCK>>>(x, n, op, g->partial);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpy(host, g->partial, blocks * sizeof(double), cudaMemcpyDeviceToHost));

    double acc = red_identity(op);
    for (b = 0; b < blocks; b++)
        acc = red_combine(acc, host[b], op);
    return acc;
}

// Gauss-Seidel / SOR in red-black ordering, with the stopping rule of the CPU solvers
static void poisson_redblack(gpu_solver *g, const double *w)
{
    int k, n = g->n;
    double dx2 = g->cfg.dx * g->cfg.dx, dy2 = g->cfg.dy * g->cfg.dy;
    double beta = (g->cfg.poisson_type == 2) ? g->cfg.beta : 1.0;
    const char *name = (g->cfg.poisson_type == 2) ? "Poisson SOR" : "Poisson";

    struct rb_state state;
    int blocks = (n + BLOCK - 1) / BLOCK;
    if (blocks > RED_BLOCKS) blocks = RED_BLOCKS;

    CUDA_CHECK(cudaMemset(g->psi, 0, n * sizeof(double)));
    CUDA_CHECK(cudaMemset(g->scratch, 0, n * sizeof(double)));
    CUDA_CHECK(cudaMemset(g->rb, 0, sizeof(struct rb_state)));

    // Sweeps are queued in batches and the host looks at the result once per
    // batch. Sweeps after the one that meets the stopping rule return at
    // once, so the answer is the same as checking after every sweep.
    for (k = 0; k < g->cfg.poisson_max_it; k += RB_BATCH)
    {
        int s, last = k + RB_BATCH < g->cfg.poisson_max_it ? k + RB_BATCH : g->cfg.poisson_max_it;
        for (s = k; s < last; s++)
        {
            LAUNCH(redblack_kernel, n, g->psi, w, g->scratch, g->nx, g->ny, dx2, dy2, beta, 0, g->rb);
            LAUNCH(redblack_kernel, n, g->psi, w, g->scratch, g->nx, g->ny, dx2, dy2, beta, 1, g->rb);
            reduce_kernel<<<blocks, BLOCK>>>(g->scratch, n, RED_SUM, g->partial);
            CUDA_CHECK(cudaGetLastError());
            rb_check_kernel<<<1, 1>>>(g->partial, blocks, g->cfg.poisson_tol, s, g->rb);
            CUDA_CHECK(cudaGetLastError());
        }
        CUDA_CHECK(cudaMemcpy(&state, g->rb, sizeof(state), cudaMemcpyDeviceToHost));
        if (state.done)
        {
            printf("%s solved in %d iterations - RSS error: %E\n", name, state.iter, state.err);
            return;
        }
    }
    printf("Error: max iterations reached for %s solver.\n", name);
    exit(1);
}

// 2D DST-I of the interior of src, left in g->spec as column spectra
static void dst2d(gpu_solver *g, const double *src)
{
    int nx = g->nx, ny = g->ny;

    LAUNCH(extend_rows_kernel, (ny - 2) * 2 * (nx - 1), src, g->ext, nx, ny);
    CUFFT_CHECK(cufftExecD2Z(g->plan_rows, g->ext, g->spec));
    LAUNCH(extend_cols_kernel, (nx - 2) * 2 * (ny - 1), g->spec, g->ext, nx, ny);
    CUFFT_CHECK(cufftExecD2Z(g->plan_cols, g->ext, g->spec));
}

// Direct solve: DST-I, divide by the eigenvalues, DST-I again, normalise.
// The sign of the right-hand side (-w) is folded into the final scale.
static void poisson_fft(gpu_solver *g, const double *w)
{
    int nx = g->nx, ny = g->ny, n = g->n;
    double inv_norm = 1.0 / (4.0 * (double)(nx - 1) * (double)(ny - 1));

    dst2d(g, w);
    LAUNCH(spectral_divide_kernel, n, g->spec, g->scratch, g->lambda_i, g->lambda_j, nx, ny);
    dst2d(g, g->scratch);
    LAUNCH(spectral_scale_kernel, n, g->spec, g->psi, -inv_norm, nx, ny);
}

// Solve nabla^2 psi = -w
static void solve_poisson(gpu_solver *g, const double *w)
{
    if (g->cfg.poisson_type == 3)
        poisson_fft(g, w);
    else
        poisson_redblack(g, w);
}

// Solve nabla^2 psi = -w and recover u = dpsi/dy, v = -dpsi/dx
static void velocity_from_vorticity(gpu_solver *g, const double *w)
{
    solve_poisson(g, w);
    LAUNCH(velocity_kernel, g->n, g->DX, g->DY, g->psi, g->u, g->v, g->n);
}

// Evaluate dw/dt into out and update u, v consistent with w
static void dwdt(gpu_solver *g, double *w, double *out)
{
    // Velocity of this stage, then the wall vorticity that goes with it, as
    // in dwdt() in fluiddyn.c
    velocity_from_vorticity(g, w);
    LAUNCH(wall_bc_kernel, g->n, g->u, g->v, g->cfg.bc, g->nx, g->ny);
    LAUNCH(vorticity_bc_kernel, g->n, g->DX, g->DY, g->u, g->v, w, g->nx, g->ny);
    LAUNCH(rhs_kernel, g->n, g->DX, g->DY, g->DX2, g->DY2, w, g->u, g->v, g->cfg.Re, out, g->n);
}

// ---------------------------------------------------------------------------
// Public interface
// ---------------------------------------------------------------------------

gpu_solver *gpu_init(const solver_config *cfg)
{
    int i, count = 0;
    int nx = cfg->nx, ny = cfg->ny, n = nx * ny;

    // Probe for a device. cudaFree(0) forces the context to be created, so a
    // device that is present but cannot be used is also reported here.
    if (cudaGetDeviceCount(&count) != cudaSuccess || count < 1)
        return NULL;
    if (cudaSetDevice(0) != cudaSuccess || cudaFree(0) != cudaSuccess)
        return NULL;

    if (cfg->poisson_type < 1 || cfg->poisson_type > 3)
    {
        printf("** Error: valid Poisson solver types are 1, 2 or 3 **\n");
        exit(1);
    }

    gpu_solver *g = (gpu_solver *)calloc(1, sizeof(gpu_solver));
    if (!g)
    {
        printf("** Error: insufficient memory **\n");
        exit(1);
    }

    g->nx = nx;
    g->ny = ny;
    g->n = n;
    g->cfg = *cfg;

    g->DX = csr_upload(cfg->DX);
    g->DY = csr_upload(cfg->DY);
    g->DX2 = csr_upload(cfg->DX2);
    g->DY2 = csr_upload(cfg->DY2);

    g->u = dev_alloc(n);
    g->v = dev_alloc(n);
    g->w = dev_alloc(n);
    g->psi = dev_alloc(n);
    g->k1 = dev_alloc(n);
    g->k2 = dev_alloc(n);
    g->k3 = dev_alloc(n);
    g->k4 = dev_alloc(n);
    g->w_tmp = dev_alloc(n);
    g->scratch = dev_alloc(n);
    g->partial = dev_alloc(RED_BLOCKS);
    CUDA_CHECK(cudaMalloc((void **)&g->rb, sizeof(struct rb_state)));

    if (g->cfg.poisson_type == 3)
    {
        int len_x = 2 * (nx - 1), len_y = 2 * (ny - 1); // lengths of the odd extensions
        size_t ext_rows = (size_t)(ny - 2) * len_x, ext_cols = (size_t)(nx - 2) * len_y;
        size_t spec_rows = (size_t)(ny - 2) * nx, spec_cols = (size_t)(nx - 2) * ny;
        double *lambda = (double *)malloc((size_t)(nx > ny ? nx : ny) * sizeof(double));
        if (!lambda)
        {
            printf("** Error: insufficient memory **\n");
            exit(1);
        }

        CUFFT_CHECK(cufftPlanMany(&g->plan_rows, 1, &len_x, NULL, 1, len_x,
                                  NULL, 1, len_x / 2 + 1, CUFFT_D2Z, ny - 2));
        CUFFT_CHECK(cufftPlanMany(&g->plan_cols, 1, &len_y, NULL, 1, len_y,
                                  NULL, 1, len_y / 2 + 1, CUFFT_D2Z, nx - 2));
        g->ext = dev_alloc(ext_rows > ext_cols ? ext_rows : ext_cols);
        CUDA_CHECK(cudaMalloc((void **)&g->spec, (spec_rows > spec_cols ? spec_rows : spec_cols) * sizeof(cufftDoubleComplex)));

        // Eigenvalues of the 2D Laplacian under the DST-I of the interior, as
        // in poisson_FFT(); rows (index i) are y, columns (index j) are x:
        //   λ_ij = (2*cos(π*(i+1)/(ny-1)) - 2) / dy²
        //         + (2*cos(π*(j+1)/(nx-1)) - 2) / dx²
        g->lambda_i = dev_alloc(ny - 2);
        g->lambda_j = dev_alloc(nx - 2);
        for (i = 0; i < ny - 2; i++)
            lambda[i] = (2.0 * cos(PI * (i + 1) / (double)(ny - 1)) - 2.0) / (g->cfg.dy * g->cfg.dy);
        CUDA_CHECK(cudaMemcpy(g->lambda_i, lambda, (ny - 2) * sizeof(double), cudaMemcpyHostToDevice));
        for (i = 0; i < nx - 2; i++)
            lambda[i] = (2.0 * cos(PI * (i + 1) / (double)(nx - 1)) - 2.0) / (g->cfg.dx * g->cfg.dx);
        CUDA_CHECK(cudaMemcpy(g->lambda_j, lambda, (nx - 2) * sizeof(double), cudaMemcpyHostToDevice));
        free(lambda);
    }
    return g;
}

void gpu_free(gpu_solver *g)
{
    if (!g) return;

    csr_free(g->DX);
    csr_free(g->DY);
    csr_free(g->DX2);
    csr_free(g->DY2);
    cudaFree(g->u);
    cudaFree(g->v);
    cudaFree(g->w);
    cudaFree(g->psi);
    cudaFree(g->k1);
    cudaFree(g->k2);
    cudaFree(g->k3);
    cudaFree(g->k4);
    cudaFree(g->w_tmp);
    cudaFree(g->scratch);
    cudaFree(g->partial);
    cudaFree(g->rb);
    if (g->cfg.poisson_type == 3)
    {
        cufftDestroy(g->plan_rows);
        cufftDestroy(g->plan_cols);
        cudaFree(g->ext);
        cudaFree(g->spec);
        cudaFree(g->lambda_i);
        cudaFree(g->lambda_j);
    }
    free(g);
}

const char *gpu_device_name(void)
{
    static cudaDeviceProp prop;

    if (cudaGetDeviceProperties(&prop, 0) != cudaSuccess)
        return "unknown";
    return prop.name;
}

void gpu_step(gpu_solver *g)
{
    int n = g->n;
    double dt = g->cfg.dt;

    // Boundary conditions
    LAUNCH(wall_bc_kernel, n, g->u, g->v, g->cfg.bc, g->nx, g->ny);
    LAUNCH(vorticity_bc_kernel, n, g->DX, g->DY, g->u, g->v, g->w, g->nx, g->ny);

    if (g->cfg.time_scheme == 1)
    {
        // Euler: single RHS evaluation, then one Poisson solve
        LAUNCH(rhs_kernel, n, g->DX, g->DY, g->DX2, g->DY2, g->w, g->u, g->v, g->cfg.Re, g->k1, n);
        LAUNCH(axpy_kernel, n, g->w, dt, g->k1, g->w, n);
        velocity_from_vorticity(g, g->w);
    }
    else
    {
        // Classical RK4: w_{n+1} = w_n + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
        dwdt(g, g->w, g->k1);
        LAUNCH(axpy_kernel, n, g->w, 0.5 * dt, g->k1, g->w_tmp, n);
        dwdt(g, g->w_tmp, g->k2);
        LAUNCH(axpy_kernel, n, g->w, 0.5 * dt, g->k2, g->w_tmp, n);
        dwdt(g, g->w_tmp, g->k3);
        LAUNCH(axpy_kernel, n, g->w, dt, g->k3, g->w_tmp, n);
        dwdt(g, g->w_tmp, g->k4);
        LAUNCH(rk4_combine_kernel, n, g->w, dt / 6.0, g->k1, g->k2, g->k3, g->k4, n);

        // Final Poisson solve so u, v are consistent with w_{n+1}
        velocity_from_vorticity(g, g->w);
    }
}

void gpu_continuity(gpu_solver *g, double *cmax, double *cmin)
{
    LAUNCH(continuity_kernel, g->n, g->DX, g->DY, g->u, g->v, g->scratch, g->n);
    *cmax = reduce(g, g->scratch, g->n, RED_MAX);
    *cmin = reduce(g, g->scratch, g->n, RED_MIN);
}

void gpu_set_fields(gpu_solver *g, const mtrx *u, const mtrx *v, const mtrx *w)
{
    size_t bytes = g->n * sizeof(double);

    if (u) CUDA_CHECK(cudaMemcpy(g->u, u->M, bytes, cudaMemcpyHostToDevice));
    if (v) CUDA_CHECK(cudaMemcpy(g->v, v->M, bytes, cudaMemcpyHostToDevice));
    if (w) CUDA_CHECK(cudaMemcpy(g->w, w->M, bytes, cudaMemcpyHostToDevice));
}

void gpu_get_fields(gpu_solver *g, mtrx *u, mtrx *v, mtrx *w)
{
    size_t bytes = g->n * sizeof(double);

    if (u) CUDA_CHECK(cudaMemcpy(u->M, g->u, bytes, cudaMemcpyDeviceToHost));
    if (v) CUDA_CHECK(cudaMemcpy(v->M, g->v, bytes, cudaMemcpyDeviceToHost));
    if (w) CUDA_CHECK(cudaMemcpy(w->M, g->w, bytes, cudaMemcpyDeviceToHost));
}

void gpu_spmv(gpu_solver *g, int op, const double *x, double *y)
{
    csr_dev A = (op == 0) ? g->DX : (op == 1) ? g->DY
                                : (op == 2)   ? g->DX2
                                              : g->DY2;
    size_t bytes = g->n * sizeof(double);

    CUDA_CHECK(cudaMemcpy(g->w_tmp, x, bytes, cudaMemcpyHostToDevice));
    LAUNCH(spmv_kernel, g->n, A, g->w_tmp, g->scratch, g->n);
    CUDA_CHECK(cudaMemcpy(y, g->scratch, bytes, cudaMemcpyDeviceToHost));
}

void gpu_poisson(gpu_solver *g, const double *w, double *psi)
{
    size_t bytes = g->n * sizeof(double);

    CUDA_CHECK(cudaMemcpy(g->w_tmp, w, bytes, cudaMemcpyHostToDevice));
    solve_poisson(g, g->w_tmp);
    CUDA_CHECK(cudaMemcpy(psi, g->psi, bytes, cudaMemcpyDeviceToHost));
}
