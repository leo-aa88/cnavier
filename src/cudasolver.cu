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
#include "diagnostics.h"
#include "fourier.h"
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
    csr_dev DXv, DYv; // velocity operators: cfg.DXv, cfg.DYv, or copies of DX, DY

    double *u, *v, *w, *psi;
    double *k1, *k2, *k3, *k4;                     // RK4 stage increments
    double *w_tmp;                                 // temporary w for intermediate stages
    double *scratch;                               // continuity field / Poisson work array
    double *partial;                               // per-block results of a reduction
    double *source;                                // vorticity source of the current stage (cfg.vorticity_source only)
    double *kolmogorov;                            // Kolmogorov source of each row (cfg.forcing), else NULL
    double *hyp1, *hyp2;                           // scratch for (-L)^p w (cfg.forcing.hyperviscosity only)
    double *uw, *vw;                               // u w and v w (cfg.advection 1 only)
    random_forcing *kicks;                         // random forcing: modes on the host, else NULL
    double *kick_ky, *kick_phase;                  // ... and on the device, with the tables
    double *kick_cx, *kick_sx, *kick_cy, *kick_sy; // of random_forcing
    mtrx source_host;                              // ... and the host field cfg.vorticity_source fills
    long steps;                                    // steps taken; the time is cfg.t0 + steps * cfg.dt
    int psi_of_w;                                  // psi was solved for the current interior of w
    long solves;                                   // Poisson solves done

    // Convergence state of the iterative Poisson solvers, kept on the device
    // so that a batch of sweeps runs without waiting for the host
    struct rb_state *rb;

    // FFT Poisson solver (poisson_type 3)
    cufftHandle plan_rows;       // batched real FFTs of the odd extensions of the rows
    cufftHandle plan_cols;       // ... and of the columns
    double *ext;                 // odd extensions, (ny-2) x 2(nx-1) or (nx-2) x 2(ny-1)
    cufftDoubleComplex *spec;    // their spectra, (ny-2) x nx or (nx-2) x ny
    double *lambda_i, *lambda_j; // eigenvalues of the 1D second differences

    // Periodic Poisson solver (cfg.periodic): 2D real FFT of the whole field
    cufftHandle plan_r2c, plan_c2r;
    double *plam_x, *plam_y; // eigenvalues of DX2 (nx/2+1 of them) and DY2 (ny)

    // The operators in Fourier space (cfg.fourier): the symbols on the device
    // (half spectrum in x: kx values; ny in y), spectra and work arrays
    struct fsym
    {
        double *d1x_re, *d1x_im, *d2x, *maskx;
        double *d1y_re, *d1y_im, *d2y, *masky;
    } fs;
    cufftDoubleComplex *f_hat, *f_hat2, *f_work; // kx * ny each
    double *f_wx, *f_wy, *f_lap, *f_tmp;         // n each

    // Spectra on the device (gpu_spectra()), set up on first use for the
    // spectra object with spectra_id() sp_id (0: none yet)
    unsigned long sp_id;
    int sp_bins;
    int *sp_bin_ptr, *sp_modes;            // the modes of each shell, in index order
    double *sp_weight, *sp_lap, *sp_ratio; // spectra_tables()
    cufftDoubleComplex *sp_hat;            // spectra of u, v, w and N
    double *sp_nl, *sp_out;                // nonlinear term; shell sums
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

// Third-order wall vorticity from psi (wall_closure 1), as
// set_wall_vorticity_psi() in fluiddyn.c: corners take the bottom/top value
__device__ inline double briley(double p0, double p1, double p2, double p3, double U, double h)
{
    return (85.0 * p0 - 108.0 * p1 + 27.0 * p2 - 4.0 * p3) / (18.0 * h * h) + 11.0 * U / (3.0 * h);
}

__global__ void vorticity_bc_psi_kernel(const double *psi, double *w, wall_bc bc, int nx, int ny,
                                        double dx, double dy)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= nx * ny) return;
    int i = k / nx, j = k % nx;

    if (i == 0)
        w[k] = briley(psi[j], psi[nx + j], psi[2 * nx + j], psi[3 * nx + j], bc.u[2], dy);
    else if (i == ny - 1)
        w[k] = briley(psi[k], psi[k - nx], psi[k - 2 * nx], psi[k - 3 * nx], -bc.u[3], dy);
    else if (j == 0)
        w[k] = briley(psi[k], psi[k + 1], psi[k + 2], psi[k + 3], -bc.v[0], dx);
    else if (j == nx - 1)
        w[k] = briley(psi[k], psi[k - 1], psi[k - 2], psi[k - 3], bc.v[1], dx);
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
                                       double cross, int nx, int ny)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= nx * ny) return;
    int i = k / nx, j = k % nx;

    if (i == 0 || i == ny - 1 || j == 0 || j == nx - 1)
        out[k] = 0.0;
    else
        out[k] = dst_mode(spec, i, j, ny) /
                 (lambda_i[i - 1] + lambda_j[j - 1] + cross * lambda_i[i - 1] * lambda_j[j - 1]);
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
    csr_dev d = {NULL, NULL, NULL};
    if (!A) return d; // not used: Fourier operators
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

// w at interior node (i, j), or at a wall node the quadratic extrapolation of
// the first three interior nodes along the wall normal, as
// f_or_extrapolated() in poisson.c
__device__ inline double w_or_extrapolated(const double *w, int i, int j, int nx, int ny)
{
    if (j == 0) return 3.0 * w[i * nx + 1] - 3.0 * w[i * nx + 2] + w[i * nx + 3];
    if (j == nx - 1) return 3.0 * w[i * nx + nx - 2] - 3.0 * w[i * nx + nx - 3] + w[i * nx + nx - 4];
    if (i == 0) return 3.0 * w[nx + j] - 3.0 * w[2 * nx + j] + w[3 * nx + j];
    if (i == ny - 1) return 3.0 * w[(ny - 2) * nx + j] - 3.0 * w[(ny - 3) * nx + j] + w[(ny - 4) * nx + j];
    return w[i * nx + j];
}

// Right-hand side of the compact operator (poisson_order 4) on the interior:
// w + (w_E + w_W + w_N + w_S - 4 w)/12
__global__ void compact_rhs_kernel(const double *w, double *out, int nx, int ny)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= nx * ny) return;
    int i = k / nx, j = k % nx;

    if (i == 0 || i == ny - 1 || j == 0 || j == nx - 1)
        out[k] = 0.0;
    else
        out[k] = w[k] + (w_or_extrapolated(w, i, j - 1, nx, ny) + w_or_extrapolated(w, i, j + 1, nx, ny) +
                         w_or_extrapolated(w, i - 1, j, nx, ny) + w_or_extrapolated(w, i + 1, j, nx, ny) -
                         4.0 * w[k]) /
                            12.0;
}

// Direct solve: DST-I, divide by the eigenvalues, DST-I again, normalise.
// The sign of the right-hand side (-w) is folded into the final scale.
static void poisson_fft(gpu_solver *g, const double *w)
{
    int nx = g->nx, ny = g->ny, n = g->n;
    double inv_norm = 1.0 / (4.0 * (double)(nx - 1) * (double)(ny - 1));
    double cross = 0.0;

    if (g->cfg.poisson_order == 4)
    {
        // The corrected right-hand side goes to scratch; dst2d() has read it
        // before spectral_divide_kernel writes scratch again
        cross = (g->cfg.dx * g->cfg.dx + g->cfg.dy * g->cfg.dy) / 12.0;
        LAUNCH(compact_rhs_kernel, n, w, g->scratch, nx, ny);
        w = g->scratch;
    }
    dst2d(g, w);
    LAUNCH(spectral_divide_kernel, n, g->spec, g->scratch, g->lambda_i, g->lambda_j, cross, nx, ny);
    dst2d(g, g->scratch);
    LAUNCH(spectral_scale_kernel, n, g->spec, g->psi, -inv_norm, nx, ny);
}

// Solve nabla^2 psi = -w
// out = -x
__global__ void negate_kernel(const double *x, double *out, int n)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < n)
        out[k] = -x[k];
}

// Divide the spectrum by the eigenvalues of DX2 + DY2 and by nx*ny (cuFFT
// does not normalise); the mean (mode 0) is set to zero
__global__ void periodic_divide_kernel(cufftDoubleComplex *spec, const double *lx, const double *ly,
                                       int kx, int ny, double inv_n)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < kx * ny)
    {
        int i = k / kx, j = k % kx;
        double scale = k == 0 ? 0.0 : inv_n / (lx[j] + ly[i]);
        spec[k].x *= scale;
        spec[k].y *= scale;
    }
}

// (DX2 + DY2) psi = -w on the periodic grid, as poisson_periodic() does
static void poisson_periodic_gpu(gpu_solver *g, const double *w)
{
    int kx = g->nx / 2 + 1;
    LAUNCH(negate_kernel, g->n, w, g->scratch, g->n);
    CUFFT_CHECK(cufftExecD2Z(g->plan_r2c, g->scratch, g->spec));
    LAUNCH(periodic_divide_kernel, kx * g->ny, g->spec, g->plam_x, g->plam_y, kx, g->ny, 1.0 / g->n);
    CUFFT_CHECK(cufftExecZ2D(g->plan_c2r, g->spec, g->psi));
}

// ---------------------------------------------------------------------------
// Operators in Fourier space (cfg.fourier), as fourier.c
// ---------------------------------------------------------------------------

enum
{
    F_DX,
    F_MINUS_DX,
    F_DY,
    F_LAP,
    F_PSI, // -(DX2 + DY2)^-1, zero for the mean
    F_FILTER,
    F_POWER // (-(DX2 + DY2))^p, p = op - F_POWER
};

// dst = scale * multiplier(op) * src, mode by mode of the half spectrum
__global__ void fourier_mul_kernel(const cufftDoubleComplex *src, cufftDoubleComplex *dst, int op,
                                   gpu_solver::fsym s, double scale, int kx, int ny)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= kx * ny) return;
    int i = k / kx, j = k % kx;
    double mask = s.maskx[j] * s.masky[i], mr, mi = 0.0;
    switch (op)
    {
    case F_DX:
        mr = s.d1x_re[j] * mask;
        mi = s.d1x_im[j] * mask;
        break;
    case F_MINUS_DX:
        mr = -s.d1x_re[j] * mask;
        mi = -s.d1x_im[j] * mask;
        break;
    case F_DY:
        mr = s.d1y_re[i] * mask;
        mi = s.d1y_im[i] * mask;
        break;
    case F_LAP:
        mr = (s.d2x[j] + s.d2y[i]) * mask;
        break;
    case F_PSI:
        mr = k == 0 ? 0.0 : -mask / (s.d2x[j] + s.d2y[i]);
        break;
    case F_FILTER:
        mr = mask;
        break;
    default:
    {
        double q = -(s.d2x[j] + s.d2y[i]);
        mr = mask;
        for (int r = F_POWER; r < op; r++)
            mr *= q;
        break;
    }
    }
    mr *= scale;
    mi *= scale;
    cufftDoubleComplex a = src[k];
    dst[k].x = mr * a.x - mi * a.y;
    dst[k].y = mr * a.y + mi * a.x;
}

// dst = scale * (DX a + DY b) from the spectra of a and b
__global__ void fourier_div_kernel(const cufftDoubleComplex *a, const cufftDoubleComplex *b,
                                   cufftDoubleComplex *dst, gpu_solver::fsym s, double scale, int kx, int ny)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= kx * ny) return;
    int i = k / kx, j = k % kx;
    double mask = s.maskx[j] * s.masky[i] * scale;
    double xr = s.d1x_re[j] * mask, xi = s.d1x_im[j] * mask, yr = s.d1y_re[i] * mask, yi = s.d1y_im[i] * mask;
    dst[k].x = xr * a[k].x - xi * a[k].y + yr * b[k].x - yi * b[k].y;
    dst[k].y = xr * a[k].y + xi * a[k].x + yr * b[k].y + yi * b[k].x;
}

// out = -(u wx + v wy) + lap / Re, from derivative arrays
__global__ void rhs_arrays_kernel(const double *wx, const double *wy, const double *lap, const double *u,
                                  const double *v, double Re, double *out, int n)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < n)
        out[k] = -u[k] * wx[k] - v[k] * wy[k] + (1.0 / Re) * lap[k];
}

// out = -(u wx + v wy), the advective nonlinear term from arrays
__global__ void nonlinear_arrays_kernel(const double *wx, const double *wy, const double *u, const double *v,
                                        double *out, int n)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < n)
        out[k] = -(u[k] * wx[k] + v[k] * wy[k]);
}

// out += 1/2 (u wx + v wy) - 1/2 div, div = DX(u w) + DY(v w): the
// skew-symmetric correction from arrays
__global__ void skew_arrays_kernel(const double *u, const double *v, const double *wx, const double *wy,
                                   const double *div, double *out, int n)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < n)
        out[k] += 0.5 * (u[k] * wx[k] + v[k] * wy[k]) - 0.5 * div[k];
}

// The spectrum of x into hat (x is copied first: cuFFT may overwrite the input)
static void f_forward(gpu_solver *g, const double *x, cufftDoubleComplex *hat)
{
    CUDA_CHECK(cudaMemcpy(g->f_tmp, x, g->n * sizeof(double), cudaMemcpyDeviceToDevice));
    CUFFT_CHECK(cufftExecD2Z(g->plan_r2c, g->f_tmp, hat));
}

// out = multiplier(op) applied to the spectrum src, back in physical space
static void f_inverse(gpu_solver *g, int op, const cufftDoubleComplex *src, double *out)
{
    int kx = g->nx / 2 + 1;
    LAUNCH(fourier_mul_kernel, kx * g->ny, src, g->f_work, op, g->fs, 1.0 / g->n, kx, g->ny);
    CUFFT_CHECK(cufftExecZ2D(g->plan_c2r, g->f_work, out));
}

// g->f_wx = DX w, g->f_wy = DY w, and with lap g->f_lap = (DX2 + DY2) w
static void f_derivatives(gpu_solver *g, const double *w, int lap)
{
    f_forward(g, w, g->f_hat);
    f_inverse(g, F_DX, g->f_hat, g->f_wx);
    f_inverse(g, F_DY, g->f_hat, g->f_wy);
    if (lap) f_inverse(g, F_LAP, g->f_hat, g->f_lap);
}

// out = DX a + DY b (out may be a or b)
static void f_divergence(gpu_solver *g, const double *a, const double *b, double *out)
{
    int kx = g->nx / 2 + 1;
    f_forward(g, a, g->f_hat);
    f_forward(g, b, g->f_hat2);
    LAUNCH(fourier_div_kernel, kx * g->ny, g->f_hat, g->f_hat2, g->f_work, g->fs, 1.0 / g->n, kx, g->ny);
    CUFFT_CHECK(cufftExecZ2D(g->plan_c2r, g->f_work, out));
}

static void solve_poisson(gpu_solver *g, const double *w)
{
    if (g->cfg.periodic)
        poisson_periodic_gpu(g, w);
    else if (g->cfg.poisson_type == 3)
        poisson_fft(g, w);
    else
        poisson_redblack(g, w);
}

// Solve nabla^2 psi = -w and recover u = dpsi/dy, v = -dpsi/dx. The solve is
// skipped when psi is already that of w: RK4's first stage would otherwise
// repeat the last solve of the step before. psi_of_w is set by a solve for
// g->w and cleared by everything that changes the interior of g->w (the wall
// entries, which the FFT solvers do not read, may change). The iterative
// solvers start from the previous psi and are never skipped, as on the CPU.
static void velocity_from_vorticity(gpu_solver *g, const double *w)
{
    if (g->cfg.fourier)
    {
        // The solve and the velocity in one pass through Fourier space;
        // f_hat2 holds the spectrum of psi
        int kx = g->nx / 2 + 1;
        if (w == g->w && g->psi_of_w)
            f_forward(g, g->psi, g->f_hat2);
        else
        {
            f_forward(g, w, g->f_hat);
            LAUNCH(fourier_mul_kernel, kx * g->ny, g->f_hat, g->f_hat2, F_PSI, g->fs, 1.0, kx, g->ny);
            f_inverse(g, F_FILTER, g->f_hat2, g->psi);
            g->solves++;
            g->psi_of_w = w == g->w;
        }
        f_inverse(g, F_DY, g->f_hat2, g->u);
        f_inverse(g, F_MINUS_DX, g->f_hat2, g->v);
        return;
    }
    if (!(w == g->w && g->psi_of_w && (g->cfg.periodic || g->cfg.poisson_type == 3)))
    {
        solve_poisson(g, w);
        g->solves++;
        g->psi_of_w = w == g->w;
    }
    LAUNCH(velocity_kernel, g->n, g->DXv, g->DYv, g->psi, g->u, g->v, g->n);
}

// out += f(t), the vorticity source, when cfg.vorticity_source is set. The source is a host
// function, so the field is filled on the host and copied over.
static void add_vorticity_source(gpu_solver *g, double t, double *out)
{
    if (!g->cfg.vorticity_source) return;
    g->cfg.vorticity_source(t, g->source_host, g->cfg.source_data);
    CUDA_CHECK(cudaMemcpy(g->source, g->source_host.M, g->n * sizeof(double), cudaMemcpyHostToDevice));
    LAUNCH(axpy_kernel, g->n, out, 1.0, g->source, out, g->n);
}

// Wall vorticity of w from u, v (wall_closure 0) or from psi (wall_closure 1)
static void wall_vorticity(gpu_solver *g, const double *u, const double *v, double *w)
{
    if (g->cfg.wall_closure == 1)
        LAUNCH(vorticity_bc_psi_kernel, g->n, g->psi, w, g->cfg.bc, g->nx, g->ny, g->cfg.dx, g->cfg.dy);
    else
        LAUNCH(vorticity_bc_kernel, g->n, g->DX, g->DY, u, v, w, g->nx, g->ny);
}

// out += -drag w + the Kolmogorov source of the row, as add_forcing_terms()
__global__ void forcing_terms_kernel(double *out, const double *w, double drag, const double *kolmogorov, int nx,
                                     int n)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < n)
        out[k] += -drag * w[k] + (kolmogorov ? kolmogorov[k / nx] : 0.0);
}

// The row tables of the kick: cos and sin of ky y + phase, entry m ny + i
__global__ void kick_rows_kernel(double *cy, double *sy, const double *ky, const double *phase, int ny, double dy,
                                 int n)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= n) return;
    int m = k / ny, i = k % ny;
    sincos(ky[m] * (i * dy) + phase[m], &sy[k], &cy[k]);
}

// w += sum over the modes of amp cos(kx x + ky y + phase) = cx cy - sx sy,
// as random_forcing_add()
__global__ void kick_kernel(double *w, const double *cx, const double *sx, const double *cy, const double *sy,
                            int modes, int nx, int ny, int n)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= n) return;
    int i = k / nx, j = k % nx;
    double s = 0.0;
    for (int m = 0; m < modes; m++)
        s += cx[m * nx + j] * cy[m * ny + i] - sx[m * nx + j] * sy[m * ny + i];
    w[k] += s;
}

// uw = u w, vw = v w
__global__ void products_kernel(const double *u, const double *v, const double *w, double *uw, double *vw, int n)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < n)
    {
        uw[k] = u[k] * w[k];
        vw[k] = v[k] * w[k];
    }
}

// out += 1/2 (u DX w + v DY w) - 1/2 (DX(u w) + DY(v w)): the skew-symmetric
// nonlinear term in place of the advective one, as skew_correction()
__global__ void skew_kernel(csr_dev DX, csr_dev DY, const double *u, const double *v, const double *w,
                            const double *uw, const double *vw, double *out, int n)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < n)
        out[k] += 0.5 * (u[k] * csr_row(DX, w, k) + v[k] * csr_row(DY, w, k)) -
                  0.5 * (csr_row(DX, uw, k) + csr_row(DY, vw, k));
}

static void add_skew_correction(gpu_solver *g, const double *w, double *out)
{
    if (g->cfg.advection != 1) return;
    if (g->cfg.fourier)
    {
        // Uses g->f_wx, g->f_wy of w, which the caller has computed
        LAUNCH(products_kernel, g->n, g->u, g->v, w, g->uw, g->vw, g->n);
        f_divergence(g, g->uw, g->vw, g->uw);
        LAUNCH(skew_arrays_kernel, g->n, g->u, g->v, g->f_wx, g->f_wy, g->uw, out, g->n);
        return;
    }
    LAUNCH(products_kernel, g->n, g->u, g->v, w, g->uw, g->vw, g->n);
    LAUNCH(skew_kernel, g->n, g->DX, g->DY, g->u, g->v, w, g->uw, g->vw, out, g->n);
}

// spec *= -scale Q^p, with Q = -(lx + ly) the eigenvalue of -(DX2 + DY2):
// (-L)^p in spectral space, times -scale
__global__ void hyper_multiply_kernel(cufftDoubleComplex *spec, const double *lx, const double *ly, int p,
                                      double scale, int kx, int ny)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < kx * ny)
    {
        double q = -(lx[k % kx] + ly[k / kx]), f = -scale;
        for (int i = 0; i < p; i++)
            f *= q;
        spec[k].x *= f;
        spec[k].y *= f;
    }
}

// out = -(u DX w + v DY w), the advective nonlinear term
__global__ void nonlinear_kernel(csr_dev DX, csr_dev DY, const double *w, const double *u, const double *v,
                                 double *out, int n)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < n)
        out[k] = -(u[k] * csr_row(DX, w, k) + v[k] * csr_row(DY, w, k));
}

#define SP_THREADS 128 // threads per shell in spectra_kernel (power of two)

// The shell sums of spectra_all() for shell blockIdx.x, into out[q * bins + b]
// for q = E, Z, TE, TZ, DE, DZ, FE, FZ (TE, TZ: the nonlinear transfer into
// the shell, which the host turns into fluxes). Each block sums its shell's
// modes and reduces them in a fixed order, so the result is deterministic.
__global__ void spectra_kernel(const cufftDoubleComplex *uh, const cufftDoubleComplex *vh,
                               const cufftDoubleComplex *wh, const cufftDoubleComplex *nh, const int *bin_ptr,
                               const int *modes, const double *weight, const double *lap, const double *ratio,
                               double inv, double nu, double drag, double nu_h, int p, double alpha_h, int bins,
                               double *out)
{
    __shared__ double sh[SPECTRA_COLUMNS][SP_THREADS];
    double acc[SPECTRA_COLUMNS] = {0.0};
    int b = blockIdx.x, t = threadIdx.x;

    for (int e = bin_ptr[b] + t; e < bin_ptr[b + 1]; e += SP_THREADS)
    {
        int k = modes[e];
        double c = weight[k] * inv, wr = wh[k].x, wi = wh[k].y, w2 = wr * wr + wi * wi;
        acc[0] += 0.5 * c * (uh[k].x * uh[k].x + uh[k].y * uh[k].y + vh[k].x * vh[k].x + vh[k].y * vh[k].y);
        acc[1] += 0.5 * c * w2;
        acc[3] += c * (wr * nh[k].x + wi * nh[k].y);
        if (k == 0)
        {
            // The mean vorticity carries no energy and feels only the drag
            acc[7] += drag * w2 * inv;
            continue;
        }
        // psi^ = -w^ / lap; energy transfer (A/Q) Re(psi^* N^); a damping of
        // symbol -sigma removes sigma (A/Q^2) |w^|^2 and sigma |w^|^2
        double Q = -lap[k], pr = wr / Q, pi = wi / Q, small = nu * Q, large = drag + alpha_h / Q;
        acc[2] += c * ratio[k] * (pr * nh[k].x + pi * nh[k].y);
        if (nu_h > 0.0)
        {
            double qp = 1.0;
            for (int i = 0; i < p; i++)
                qp *= Q;
            small += nu_h * qp;
        }
        acc[4] += small * c * ratio[k] / Q * w2;
        acc[5] += small * c * w2;
        acc[6] += large * c * ratio[k] / Q * w2;
        acc[7] += large * c * w2;
    }
    for (int q = 0; q < SPECTRA_COLUMNS; q++)
        sh[q][t] = acc[q];
    __syncthreads();
    for (int half = SP_THREADS / 2; half > 0; half /= 2)
    {
        if (t < half)
            for (int q = 0; q < SPECTRA_COLUMNS; q++)
                sh[q][t] += sh[q][t + half];
        __syncthreads();
    }
    if (t < SPECTRA_COLUMNS) out[t * bins + b] = sh[t][0];
}

// out = -(DX2 + DY2) in
__global__ void neg_laplacian_kernel(csr_dev DX2, csr_dev DY2, const double *in, double *out, int n)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < n)
        out[k] = -(csr_row(DX2, in, k) + csr_row(DY2, in, k));
}

// out += -drag w + the Kolmogorov source - nu_h (-L)^p w - alpha_h psi, as
// add_forcing_terms() in fluiddyn.c; g->psi must be that of w
static void add_forcing_terms(gpu_solver *g, double *out, const double *w)
{
    const forcing_config *fc = &g->cfg.forcing;
    if (fc->hyperviscosity > 0.0 && g->cfg.fourier)
    {
        // With the other Fourier operators, masked like them (fourier_power())
        int kx = g->nx / 2 + 1;
        f_forward(g, w, g->f_hat);
        LAUNCH(fourier_mul_kernel, kx * g->ny, g->f_hat, g->f_work, F_POWER + fc->hyper_order, g->fs,
               -fc->hyperviscosity / g->n, kx, g->ny);
        CUFFT_CHECK(cufftExecZ2D(g->plan_c2r, g->f_work, g->hyp2));
        LAUNCH(axpy_kernel, g->n, out, 1.0, g->hyp2, out, g->n);
    }
    else if (fc->hyperviscosity > 0.0 && g->cfg.periodic)
    {
        // DX2 + DY2 is circulant, so (-L)^p is diagonal in Fourier space with
        // eigenvalue Q^p: one transform each way instead of 2p sparse products
        // (cuFFT may overwrite the input of a multidimensional transform, so
        // w is copied first)
        int kx = g->nx / 2 + 1;
        CUDA_CHECK(cudaMemcpy(g->hyp1, w, g->n * sizeof(double), cudaMemcpyDeviceToDevice));
        CUFFT_CHECK(cufftExecD2Z(g->plan_r2c, g->hyp1, g->spec));
        LAUNCH(hyper_multiply_kernel, kx * g->ny, g->spec, g->plam_x, g->plam_y, fc->hyper_order,
               fc->hyperviscosity / g->n, kx, g->ny);
        CUFFT_CHECK(cufftExecZ2D(g->plan_c2r, g->spec, g->hyp2));
        LAUNCH(axpy_kernel, g->n, out, 1.0, g->hyp2, out, g->n);
    }
    else if (fc->hyperviscosity > 0.0)
    {
        const double *a = w;
        double *bufs[2] = {g->hyp1, g->hyp2};
        for (int i = 0; i < fc->hyper_order; i++)
        {
            LAUNCH(neg_laplacian_kernel, g->n, g->DX2, g->DY2, a, bufs[i % 2], g->n);
            a = bufs[i % 2];
        }
        LAUNCH(axpy_kernel, g->n, out, -fc->hyperviscosity, a, out, g->n);
    }
    if (fc->hypodrag > 0.0) LAUNCH(axpy_kernel, g->n, out, -fc->hypodrag, g->psi, out, g->n);
    if (fc->drag == 0.0 && !g->kolmogorov) return;
    LAUNCH(forcing_terms_kernel, g->n, out, w, g->cfg.forcing.drag, g->kolmogorov, g->nx, g->n);
}

// out = -(u DX w + v DY w) + (DX2 + DY2) w / Re, the transport right-hand
// side, by the sparse operators or in Fourier space (then leaving DX w, DY w
// in g->f_wx, g->f_wy for the skew-symmetric correction)
static void rhs(gpu_solver *g, const double *w, double *out)
{
    if (g->cfg.fourier)
    {
        f_derivatives(g, w, 1);
        LAUNCH(rhs_arrays_kernel, g->n, g->f_wx, g->f_wy, g->f_lap, g->u, g->v, g->cfg.Re, out, g->n);
    }
    else
        LAUNCH(rhs_kernel, g->n, g->DX, g->DY, g->DX2, g->DY2, w, g->u, g->v, g->cfg.Re, out, g->n);
}

// Evaluate dw/dt at time t into out and update u, v consistent with w
static void dwdt(gpu_solver *g, double *w, double *out, double t)
{
    // Velocity of this stage, then the wall vorticity that goes with it, as
    // in dwdt() in fluiddyn.c
    velocity_from_vorticity(g, w);
    if (!g->cfg.periodic)
    {
        LAUNCH(wall_bc_kernel, g->n, g->u, g->v, g->cfg.bc, g->nx, g->ny);
        wall_vorticity(g, g->u, g->v, w);
    }
    rhs(g, w, out);
    add_skew_correction(g, w, out);
    add_vorticity_source(g, t, out);
    add_forcing_terms(g, out, w);
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
    g->DXv = cfg->DXv ? csr_upload(cfg->DXv) : g->DX;
    g->DYv = cfg->DYv ? csr_upload(cfg->DYv) : g->DY;

    g->u = dev_alloc(n);
    g->v = dev_alloc(n);
    g->w = dev_alloc(n);
    g->psi = dev_alloc(n);
    g->k1 = dev_alloc(n);
    g->k2 = dev_alloc(n);
    g->k3 = dev_alloc(n);
    g->k4 = dev_alloc(n);
    g->w_tmp = dev_alloc(n);
    {
        double *rows = kolmogorov_rows(&cfg->forcing, ny, cfg->dy, cfg->dy * (cfg->periodic ? ny : ny - 1));
        if (rows)
        {
            g->kolmogorov = dev_alloc(ny);
            CUDA_CHECK(cudaMemcpy(g->kolmogorov, rows, ny * sizeof(double), cudaMemcpyHostToDevice));
            free(rows);
        }
    }
    if (cfg->advection == 1)
    {
        if (!cfg->periodic)
        {
            printf("** Error: the skew-symmetric nonlinear term needs a periodic grid **\n");
            exit(1);
        }
        g->uw = dev_alloc(n);
        g->vw = dev_alloc(n);
    }
    if (cfg->forcing.hyperviscosity > 0.0)
    {
        if (cfg->forcing.hyper_order < 2)
        {
            printf("** Error: the hyperviscosity order must be at least 2 **\n");
            exit(1);
        }
        g->hyp1 = dev_alloc(n);
        g->hyp2 = dev_alloc(n);
    }
    if (cfg->forcing.random_rate > 0.0 && !cfg->periodic)
    {
        printf("** Error: random forcing needs a periodic grid **\n");
        exit(1);
    }
    g->kicks = NULL;
    if (cfg->forcing.random_rate > 0.0)
    {
        periodic_symbols sym;
        periodic_symbols_of(cfg, &sym);
        g->kicks = random_forcing_setup(&cfg->forcing, nx, ny, cfg->dx, cfg->dy, cfg->dt, &sym);
        periodic_symbols_free(&sym);
    }
    if (g->kicks)
    {
        int modes = g->kicks->modes;
        g->kick_ky = dev_alloc(modes);
        g->kick_phase = dev_alloc(modes);
        g->kick_cx = dev_alloc(modes * nx);
        g->kick_sx = dev_alloc(modes * nx);
        g->kick_cy = dev_alloc(modes * ny);
        g->kick_sy = dev_alloc(modes * ny);
        CUDA_CHECK(cudaMemcpy(g->kick_ky, g->kicks->ky, modes * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(g->kick_cx, g->kicks->cx, (size_t)modes * nx * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(g->kick_sx, g->kicks->sx, (size_t)modes * nx * sizeof(double), cudaMemcpyHostToDevice));
    }
    if (cfg->vorticity_source)
    {
        g->source = dev_alloc(n);
        g->source_host = initm(ny, nx);
    }
    g->scratch = dev_alloc(n);
    g->partial = dev_alloc(RED_BLOCKS);
    CUDA_CHECK(cudaMalloc((void **)&g->rb, sizeof(struct rb_state)));

    if (cfg->poisson_order == 4 && cfg->poisson_type != 3 && !cfg->periodic)
    {
        printf("** Error: the fourth-order Poisson operator needs the FFT solver (poisson_type 3) **\n");
        exit(1);
    }
    if (cfg->periodic)
    {
        int kx = nx / 2 + 1;
        double *lx = (double *)malloc((size_t)kx * sizeof(double)), *ly = (double *)malloc((size_t)ny * sizeof(double));
        if (cfg->poisson_type != 3)
        {
            printf("** Error: periodic boundaries need the FFT Poisson solver (poisson_type 3) **\n");
            exit(1);
        }
        if (!lx || !ly)
        {
            printf("** Error: insufficient memory **\n");
            exit(1);
        }
        CUFFT_CHECK(cufftPlan2d(&g->plan_r2c, ny, nx, CUFFT_D2Z));
        CUFFT_CHECK(cufftPlan2d(&g->plan_c2r, ny, nx, CUFFT_Z2D));
        CUDA_CHECK(cudaMalloc((void **)&g->spec, (size_t)kx * ny * sizeof(cufftDoubleComplex)));
        periodic_symbols sym;
        periodic_symbols_of(cfg, &sym);
        for (i = 0; i < kx; i++)
            lx[i] = sym.d2x[i];
        for (i = 0; i < ny; i++)
            ly[i] = sym.d2y[i];
        if (cfg->fourier)
        {
            // The symbols on the device: the half spectrum in x, all of y
            double **dst[8] = {&g->fs.d1x_re, &g->fs.d1x_im, &g->fs.d2x, &g->fs.maskx,
                               &g->fs.d1y_re, &g->fs.d1y_im, &g->fs.d2y, &g->fs.masky};
            const double *src[8] = {sym.d1x_re, sym.d1x_im, sym.d2x, sym.maskx, sym.d1y_re, sym.d1y_im, sym.d2y, sym.masky};
            for (int q = 0; q < 8; q++)
            {
                int m = q < 4 ? kx : ny;
                *dst[q] = dev_alloc(m);
                CUDA_CHECK(cudaMemcpy(*dst[q], src[q], m * sizeof(double), cudaMemcpyHostToDevice));
            }
            CUDA_CHECK(cudaMalloc((void **)&g->f_hat, (size_t)kx * ny * sizeof(cufftDoubleComplex)));
            CUDA_CHECK(cudaMalloc((void **)&g->f_hat2, (size_t)kx * ny * sizeof(cufftDoubleComplex)));
            CUDA_CHECK(cudaMalloc((void **)&g->f_work, (size_t)kx * ny * sizeof(cufftDoubleComplex)));
            g->f_wx = dev_alloc(n);
            g->f_wy = dev_alloc(n);
            g->f_lap = dev_alloc(n);
            g->f_tmp = dev_alloc(n);
        }
        periodic_symbols_free(&sym);
        g->plam_x = dev_alloc(kx);
        g->plam_y = dev_alloc(ny);
        CUDA_CHECK(cudaMemcpy(g->plam_x, lx, kx * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(g->plam_y, ly, ny * sizeof(double), cudaMemcpyHostToDevice));
        free(lx);
        free(ly);
    }
    else if (g->cfg.poisson_type == 3)
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
    if (g->cfg.DXv) csr_free(g->DXv);
    if (g->cfg.DYv) csr_free(g->DYv);
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
    cudaFree(g->kolmogorov);
    cudaFree(g->hyp1);
    cudaFree(g->hyp2);
    cudaFree(g->uw);
    cudaFree(g->vw);
    if (g->kicks)
    {
        cudaFree(g->kick_ky);
        cudaFree(g->kick_phase);
        cudaFree(g->kick_cx);
        cudaFree(g->kick_sx);
        cudaFree(g->kick_cy);
        cudaFree(g->kick_sy);
        random_forcing_free(g->kicks);
    }
    if (g->source)
    {
        cudaFree(g->source);
        freem(&g->source_host);
    }
    if (g->cfg.periodic)
    {
        cufftDestroy(g->plan_r2c);
        cufftDestroy(g->plan_c2r);
        cudaFree(g->spec);
        cudaFree(g->plam_x);
        cudaFree(g->plam_y);
        double *fourier_arrays[12] = {g->fs.d1x_re, g->fs.d1x_im, g->fs.d2x, g->fs.maskx, g->fs.d1y_re, g->fs.d1y_im,
                                      g->fs.d2y, g->fs.masky, g->f_wx, g->f_wy, g->f_lap, g->f_tmp};
        for (int q = 0; q < 12; q++)
            cudaFree(fourier_arrays[q]);
        cudaFree(g->f_hat);
        cudaFree(g->f_hat2);
        cudaFree(g->f_work);
        cudaFree(g->sp_bin_ptr);
        cudaFree(g->sp_modes);
        cudaFree(g->sp_weight);
        cudaFree(g->sp_lap);
        cudaFree(g->sp_ratio);
        cudaFree(g->sp_hat);
        cudaFree(g->sp_nl);
        cudaFree(g->sp_out);
    }
    else if (g->cfg.poisson_type == 3)
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

// With dealiasing, w without the modes the 2/3 rule cuts
static void dealias(gpu_solver *g)
{
    if (!g->cfg.fourier || !g->cfg.dealias) return;
    f_forward(g, g->w, g->f_hat);
    f_inverse(g, F_FILTER, g->f_hat, g->w);
}

void gpu_step(gpu_solver *g)
{
    int n = g->n;
    double dt = g->cfg.dt, t = g->cfg.t0 + (double)g->steps * dt;

    // Boundary conditions (walls only). The formula from psi needs psi,
    // which a first step does not have yet.
    if (!g->cfg.periodic && g->cfg.wall_closure == 1 && g->steps == 0)
        velocity_from_vorticity(g, g->w);
    if (!g->cfg.periodic)
    {
        LAUNCH(wall_bc_kernel, n, g->u, g->v, g->cfg.bc, g->nx, g->ny);
        wall_vorticity(g, g->u, g->v, g->w);
    }

    // Random forcing: phases drawn on the host by the same generator as the
    // CPU solver's, the kick built on the device; as in step()
    if (g->kicks)
    {
        random_forcing_draw(g->kicks);
        CUDA_CHECK(cudaMemcpy(g->kick_phase, g->kicks->phase, g->kicks->modes * sizeof(double),
                              cudaMemcpyHostToDevice));
        LAUNCH(kick_rows_kernel, g->kicks->modes * g->ny, g->kick_cy, g->kick_sy, g->kick_ky, g->kick_phase, g->ny,
               g->cfg.dy, g->kicks->modes * g->ny);
        LAUNCH(kick_kernel, n, g->w, g->kick_cx, g->kick_sx, g->kick_cy, g->kick_sy, g->kicks->modes, g->nx, g->ny, n);
        g->psi_of_w = 0;
        if (g->cfg.time_scheme == 1) velocity_from_vorticity(g, g->w);
    }
    // hypodrag needs the psi of w, which a first step does not have yet
    if (g->cfg.time_scheme == 1 && g->cfg.forcing.hypodrag > 0.0 && g->steps == 0 && !g->kicks)
        velocity_from_vorticity(g, g->w);

    if (g->cfg.time_scheme == 1)
    {
        // Euler: single RHS evaluation, then one Poisson solve
        rhs(g, g->w, g->k1);
        add_skew_correction(g, g->w, g->k1);
        add_vorticity_source(g, t, g->k1);
        add_forcing_terms(g, g->k1, g->w);
        LAUNCH(axpy_kernel, n, g->w, dt, g->k1, g->w, n);
        g->psi_of_w = 0;
        dealias(g);
        velocity_from_vorticity(g, g->w);
    }
    else
    {
        // Classical RK4: w_{n+1} = w_n + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
        dwdt(g, g->w, g->k1, t);
        LAUNCH(axpy_kernel, n, g->w, 0.5 * dt, g->k1, g->w_tmp, n);
        dwdt(g, g->w_tmp, g->k2, t + 0.5 * dt);
        LAUNCH(axpy_kernel, n, g->w, 0.5 * dt, g->k2, g->w_tmp, n);
        dwdt(g, g->w_tmp, g->k3, t + 0.5 * dt);
        LAUNCH(axpy_kernel, n, g->w, dt, g->k3, g->w_tmp, n);
        dwdt(g, g->w_tmp, g->k4, t + dt);
        LAUNCH(rk4_combine_kernel, n, g->w, dt / 6.0, g->k1, g->k2, g->k3, g->k4, n);
        g->psi_of_w = 0;
        dealias(g);

        // Final Poisson solve so u, v are consistent with w_{n+1}
        velocity_from_vorticity(g, g->w);
    }

    // The wall vorticity of the new velocity in place of the wall entries the
    // update advanced, with the wall velocities imposed on copies of u and v,
    // as at the end of step()
    if (!g->cfg.periodic)
    {
        CUDA_CHECK(cudaMemcpy(g->k1, g->u, n * sizeof(double), cudaMemcpyDeviceToDevice));
        CUDA_CHECK(cudaMemcpy(g->k2, g->v, n * sizeof(double), cudaMemcpyDeviceToDevice));
        LAUNCH(wall_bc_kernel, n, g->k1, g->k2, g->cfg.bc, g->nx, g->ny);
        wall_vorticity(g, g->k1, g->k2, g->w);
    }
    g->steps++;
}

// Weighted integrand of compute_integrals() at each node: which = 0 energy
// (wall velocities on the wall nodes), 1 enstrophy, 2 palinstrophy, 3 the
// Kolmogorov work u sin(k y)
__global__ void integrand_kernel(csr_dev DX, csr_dev DY, const double *dwx, const double *dwy, const double *u,
                                 const double *v, const double *w, wall_bc bc, int periodic, int which, int nx, int ny,
                                 double dy, double kk, double *out)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= nx * ny) return;
    int i = k / nx, j = k % nx, wall = -1;
    double weight = 1.0, val;

    if (!periodic)
    {
        weight = ((i == 0 || i == ny - 1) ? 0.5 : 1.0) * ((j == 0 || j == nx - 1) ? 0.5 : 1.0);
        wall = j == 0 ? 0 : j == nx - 1 ? 1
                        : i == 0        ? 2
                        : i == ny - 1   ? 3
                                        : -1;
    }
    if (which == 0)
    {
        double uu = wall < 0 ? u[k] : bc.u[wall], vv = wall < 0 ? v[k] : bc.v[wall];
        val = uu * uu + vv * vv;
    }
    else if (which == 1)
        val = w[k] * w[k];
    else if (which == 3)
        val = (wall < 0 ? u[k] : bc.u[wall]) * sin(kk * (i * dy));
    else
    {
        // DX w, DY w: given (Fourier operators) or from the sparse ones
        double wx = dwx ? dwx[k] : csr_row(DX, w, k), wy = dwy ? dwy[k] : csr_row(DY, w, k);
        val = wx * wx + wy * wy;
    }
    out[k] = weight * val;
}

void gpu_integrals(gpu_solver *g, double *E, double *Z, double *P, double *I)
{
    double *res[3] = {E, Z, P};
    double norm = g->cfg.periodic ? (double)g->nx * g->ny : (double)(g->nx - 1) * (g->ny - 1);
    double Ly = g->cfg.dy * (g->cfg.periodic ? g->ny : g->ny - 1), kk = 2.0 * PI * g->cfg.forcing.kolmogorov_n / Ly;
    int q;

    if (g->cfg.fourier) f_derivatives(g, g->w, 0);
    const double *dwx = g->cfg.fourier ? g->f_wx : NULL, *dwy = g->cfg.fourier ? g->f_wy : NULL;
    for (q = 0; q < 3; q++)
    {
        LAUNCH(integrand_kernel, g->n, g->DX, g->DY, dwx, dwy, g->u, g->v, g->w, g->cfg.bc, g->cfg.periodic, q, g->nx,
               g->ny, g->cfg.dy, kk, g->scratch);
        *res[q] = 0.5 * reduce(g, g->scratch, g->n, RED_SUM) / norm;
    }
    *I = g->cfg.forcing.random_rate;
    if (g->cfg.forcing.kolmogorov_amp != 0.0)
    {
        LAUNCH(integrand_kernel, g->n, g->DX, g->DY, dwx, dwy, g->u, g->v, g->w, g->cfg.bc, g->cfg.periodic, 3, g->nx,
               g->ny, g->cfg.dy, kk, g->scratch);
        *I += g->cfg.forcing.kolmogorov_amp * reduce(g, g->scratch, g->n, RED_SUM) / norm;
    }
}

// Upload the mode tables of s, with the modes of each shell listed together
static void spectra_upload(gpu_solver *g, const spectra *s)
{
    int kx = g->nx / 2 + 1, modes = kx * g->ny, bins = spectra_bins(s), b, k;
    const int *bin;
    const double *weight, *lap, *ratio;
    int *ptr = (int *)calloc(bins + 1, sizeof(int)), *list = (int *)malloc(modes * sizeof(int));
    int *fill = (int *)malloc(bins * sizeof(int));

    if (!ptr || !list || !fill)
    {
        printf("** Error: insufficient memory **\n");
        exit(1);
    }
    spectra_tables(s, &bin, &weight, &lap, &ratio);
    for (k = 0; k < modes; k++)
        ptr[bin[k] + 1]++;
    for (b = 0; b < bins; b++)
    {
        ptr[b + 1] += ptr[b];
        fill[b] = ptr[b];
    }
    for (k = 0; k < modes; k++)
        list[fill[bin[k]]++] = k;
    if (!g->sp_hat)
    {
        CUDA_CHECK(cudaMalloc((void **)&g->sp_hat, 4 * (size_t)modes * sizeof(cufftDoubleComplex)));
        g->sp_nl = dev_alloc(g->n);
        g->sp_weight = dev_alloc(modes);
        g->sp_lap = dev_alloc(modes);
        g->sp_ratio = dev_alloc(modes);
        CUDA_CHECK(cudaMalloc((void **)&g->sp_modes, modes * sizeof(int)));
    }
    cudaFree(g->sp_bin_ptr);
    cudaFree(g->sp_out);
    CUDA_CHECK(cudaMalloc((void **)&g->sp_bin_ptr, (bins + 1) * sizeof(int)));
    g->sp_out = dev_alloc(SPECTRA_COLUMNS * bins);
    CUDA_CHECK(cudaMemcpy(g->sp_bin_ptr, ptr, (bins + 1) * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(g->sp_modes, list, modes * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(g->sp_weight, weight, modes * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(g->sp_lap, lap, modes * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(g->sp_ratio, ratio, modes * sizeof(double), cudaMemcpyHostToDevice));
    g->sp_id = spectra_id(s);
    g->sp_bins = bins;
    free(ptr);
    free(list);
    free(fill);
}

void gpu_spectra(gpu_solver *g, const spectra *s, double *out)
{
    int kx = g->nx / 2 + 1, modes = kx * g->ny, B, b, q;
    const forcing_config *fc = &g->cfg.forcing;
    double inv = 1.0 / ((double)g->n * g->n);
    const double *fields[3] = {g->u, g->v, g->w};

    if (!g->cfg.periodic)
    {
        printf("** Error: spectra need a periodic grid **\n");
        exit(1);
    }
    if (g->sp_id != spectra_id(s)) spectra_upload(g, s);
    B = g->sp_bins;

    // The nonlinear term of the current fields, in the solver's form
    if (g->cfg.fourier)
    {
        f_derivatives(g, g->w, 0);
        LAUNCH(nonlinear_arrays_kernel, g->n, g->f_wx, g->f_wy, g->u, g->v, g->sp_nl, g->n);
    }
    else
        LAUNCH(nonlinear_kernel, g->n, g->DX, g->DY, g->w, g->u, g->v, g->sp_nl, g->n);
    add_skew_correction(g, g->w, g->sp_nl);

    // Spectra of u, v, w (copied first: cuFFT may overwrite the input) and N
    for (q = 0; q < 3; q++)
    {
        CUDA_CHECK(cudaMemcpy(g->scratch, fields[q], g->n * sizeof(double), cudaMemcpyDeviceToDevice));
        CUFFT_CHECK(cufftExecD2Z(g->plan_r2c, g->scratch, g->sp_hat + (size_t)q * modes));
    }
    CUFFT_CHECK(cufftExecD2Z(g->plan_r2c, g->sp_nl, g->sp_hat + 3 * (size_t)modes));

    spectra_kernel<<<B, SP_THREADS>>>(g->sp_hat, g->sp_hat + modes, g->sp_hat + 2 * (size_t)modes,
                                      g->sp_hat + 3 * (size_t)modes, g->sp_bin_ptr, g->sp_modes, g->sp_weight,
                                      g->sp_lap, g->sp_ratio, inv, 1.0 / g->cfg.Re, fc->drag, fc->hyperviscosity,
                                      fc->hyper_order, fc->hypodrag, B, g->sp_out);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpy(out, g->sp_out, SPECTRA_COLUMNS * B * sizeof(double), cudaMemcpyDeviceToHost));

    // Fluxes: what the nonlinear term removes from all shells up to k
    for (b = 0; b < B; b++)
    {
        out[2 * B + b] = (b ? out[2 * B + b - 1] : 0.0) - out[2 * B + b];
        out[3 * B + b] = (b ? out[3 * B + b - 1] : 0.0) - out[3 * B + b];
    }
}

void gpu_continuity(gpu_solver *g, double *cmax, double *cmin)
{
    if (g->cfg.fourier)
        f_divergence(g, g->u, g->v, g->scratch);
    else
        LAUNCH(continuity_kernel, g->n, g->DXv, g->DYv, g->u, g->v, g->scratch, g->n);
    *cmax = reduce(g, g->scratch, g->n, RED_MAX);
    *cmin = reduce(g, g->scratch, g->n, RED_MIN);
}

void gpu_set_fields(gpu_solver *g, const mtrx *u, const mtrx *v, const mtrx *w)
{
    size_t bytes = g->n * sizeof(double);

    if (u) CUDA_CHECK(cudaMemcpy(g->u, u->M, bytes, cudaMemcpyHostToDevice));
    if (v) CUDA_CHECK(cudaMemcpy(g->v, v->M, bytes, cudaMemcpyHostToDevice));
    if (w)
    {
        CUDA_CHECK(cudaMemcpy(g->w, w->M, bytes, cudaMemcpyHostToDevice));
        g->psi_of_w = 0;
    }
}

long gpu_poisson_solves(const gpu_solver *g)
{
    return g->solves;
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
