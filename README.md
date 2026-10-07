# cnavier

Explicit incompressible Navier-Stokes solver written in C.

Solves the 2D lid-driven cavity problem using the **vorticity-streamfunction formulation**, which eliminates pressure from the equations and enforces incompressibility exactly.

## Results

![Re=1000 lid-driven cavity flow](re1000.gif)

*Vorticity field for the lid-driven cavity at Re=1000.*

## Method

The governing equations are:

```
∂ω/∂t + u·∂ω/∂x + v·∂ω/∂y = (1/Re) ∇²ω          (vorticity transport)
∇²ψ = −ω                                        (Poisson equation)
u = ∂ψ/∂y,  v = −∂ψ/∂x                          (velocity recovery)
ω = ∂v/∂x − ∂u/∂y                               (vorticity definition)
```

At each timestep:
1. Compute vorticity boundary conditions from current velocity field
2. Evaluate spatial derivatives of ω using finite differences
3. Advance ω in time (Euler or RK4)
4. Solve the Poisson equation for ψ
5. Recover u and v from ψ

## Features

- **Spatial discretisation**: finite differences of selectable order (2nd, 4th, or 6th)
- **Time integration**: explicit Euler or classical RK4. The wall vorticity is updated once per step, which limits both to first order in time overall, and RK4 is currently no more accurate than Euler at the same step (measured; see [Tests](#tests))
- **Poisson solver**: three options — Gauss-Seidel, SOR, or FFTW3-based direct DST-I solver (default)
- **Sparse operators**: 2D derivative operators built as CSR sparse matrices via Kronecker products, replacing dense O(n³) matrix-vector multiplies with O(7n) SpMV
- **GPU acceleration**: optional CUDA backend that runs the whole time loop on an NVIDIA GPU (see [CUDA](#cuda-gpu))
- **Output**: VTK files for visualisation in ParaView

## Dependencies

- `gcc`
- `libfftw3` — required for the FFT Poisson solver
- **OpenMP** — optional; used only if you build with `make OPENMP=1` (see below)
- **CUDA toolkit** (`nvcc`, cuFFT) and an NVIDIA GPU — optional; used only if you build with `make CUDA=1` (see [CUDA](#cuda-gpu))

**Ubuntu/Debian**
```bash
sudo apt install libfftw3-dev
```

For parallel builds, install a full GCC toolchain (OpenMP comes with GCC as **libgomp**; nothing extra beyond the compiler is usually required):
```bash
sudo apt install build-essential
```

**Arch Linux**
```bash
sudo pacman -S fftw
```

GCC on Arch already includes OpenMP support for `make OPENMP=1`.

**macOS**
```bash
brew install fftw
```

Apple’s default compiler is often Clang without OpenMP enabled for `-fopenmp`. To use OpenMP, install GCC from Homebrew and point `CC` at it (the exact version suffix may vary):
```bash
brew install gcc
CC=gcc-14 make OPENMP=1
```
Alternatively, install `libomp` and use Clang with the appropriate `-fopenmp` / `-lomp` flags if you maintain a custom toolchain.

### OpenMP (quick facts)

You typically **do not** install a package named “openmp” alone. **GCC** implements OpenMP and links the runtime (**libgomp** on Linux). If `gcc -fopenmp` works on your system, `make OPENMP=1` should work.

## Build

```bash
make
```

Optional **OpenMP** (shared-memory parallelism in the hot loops — SpMV, Euler/RK4 RHS, Poisson):

```bash
make OPENMP=1
```

The default `make` build is unchanged (no OpenMP). Results are bit-identical to the serial build with any number of threads (with Gauss-Seidel/SOR, identical to other OpenMP runs: the parallel sweep order differs from the serial one).

**Threads.** Without `OMP_NUM_THREADS`, the solver uses one thread per physical core it may run on (the distinct cores among the CPUs in its affinity mask, read from Linux sysfs; elsewhere the OpenMP default applies) and prints the count at start-up. Set `OMP_NUM_THREADS` to override. If a thread placement is set (`OMP_PROC_BIND`, `OMP_PLACES` or `GOMP_CPU_AFFINITY`), the OpenMP default is kept, which is one thread per logical CPU; combine the placement with `OMP_NUM_THREADS` to choose the count, e.g. `OMP_PLACES=cores OMP_NUM_THREADS=8`.

Time per step on an i7-12650H under WSL2 (which presents it as 8 cores with 2 hardware threads each), RK4 + FFT, best of 2–3 runs:

| Threads | 129×129 | 1025×1025 |
|---|---|---|
| serial build | 3.66 ms | 326 ms |
| 1 | 3.55 ms | 330 ms |
| 2 | 2.06 ms | 213 ms |
| 4 | 1.57 ms | 179 ms |
| 8 (default here) | 1.55 ms | 157 ms |
| 16 | 2.16 ms | 191 ms |

The 8- and 16-thread rows come from one interleaved comparison. The speed-up levels off at about 2× from 4 threads on: the sparse derivatives and the transforms are limited by memory bandwidth rather than by arithmetic, so a second hardware thread per core does not help, and 16 threads were 25–40% slower than 8.

**On a busy machine, use fewer threads.** A step has about 70 parallel loops, each ending at a barrier where all threads wait for the slowest. When other programs occupy the cores, a thread that is descheduled holds everyone up: with other jobs running, 16 threads took 577 ms per step at 1025×1025 and the 8-thread default slowed down severalfold too. On a shared machine set `OMP_NUM_THREADS` below the number of free cores, or `OMP_WAIT_POLICY=passive` so that waiting threads sleep instead of spinning.

Loops over grids of fewer than 2048 points (`OMP_MIN_WORK` in `linearalg.h`) stay serial, since below that waking the threads costs more than it saves.

Optional **CUDA** (GPU) build — see [CUDA](#cuda-gpu) for details:

```bash
make CUDA=1
```

Switching between `make`, `make OPENMP=1` and `make CUDA=1` rebuilds everything automatically; no `make clean` is needed in between.

This produces the `cnavier` binary. Output VTK files are written to `output/` — create it first:

```bash
mkdir -p output
./cnavier
```

To clean build artifacts:
```bash
make clean
```

## CUDA (GPU)

`make CUDA=1` builds the same `cnavier` binary with a GPU backend. It is entirely optional: the default build does not need CUDA, and a CUDA build started on a machine without a usable GPU prints a notice and runs on the CPU.

**Requirements**: an NVIDIA GPU with a working driver, and the CUDA toolkit (`nvcc` and cuFFT). On Ubuntu/Debian the distribution package is enough:

```bash
sudo apt install nvidia-cuda-toolkit
```

NVIDIA's own packages work too. `nvcc` is taken from `PATH`, falling back to `/usr/local/cuda/bin/nvcc`; set `NVCC` or `CUDA_HOME` to point elsewhere:

```bash
make CUDA=1 CUDA_HOME=/opt/cuda
```

**Usage**: a CUDA build uses the GPU by default. Pass `--cpu` to run the CPU path from the same binary, which is handy for comparing the two:

```bash
./cnavier          # GPU
./cnavier --cpu    # CPU
```

**What runs on the GPU**: everything inside the time loop. The fields (`u`, `v`, `ω`, `ψ`, RK4 stages) and the four CSR derivative operators are uploaded once and stay on the device; boundary conditions, derivatives (CSR SpMV kernels), Euler/RK4 and all three Poisson solvers run as kernels. Data returns to the host only for the VTK files, the per-iteration continuity min/max and the final centerline profiles. All arithmetic is double precision, as on the CPU.

- *FFT solver*: cuFFT has no sine transform, so the DST-I is computed with batches of 1D real-to-complex FFTs of odd extensions, first along the rows and then along the columns, as on the CPU.
- *Gauss-Seidel / SOR*: red-black ordering, the same one the OpenMP build uses, with the same stopping rule. The convergence test runs on the device and the host looks at it once every 32 sweeps; later sweeps in the batch do nothing once it is met, so the result is the same as checking after every sweep.

**Agreement with the CPU**

- *FFT solver*: on the default case (Re=100, 64×64, RK4 + FFT, 6000 steps) the GPU and CPU runs write identical centerline profiles and byte-identical VTK files.
- *Gauss-Seidel / SOR*: the answer depends on which CPU build you compare with. The default tolerance (`poisson_tol = 1E-3`) stops the iteration well before convergence, so the sweep order shows in the result. Against the `OPENMP=1` build, which also sweeps red-black, the GPU takes the same number of iterations and agrees to round-off. Against the default serial build, which sweeps lexicographically, it does not: with SOR on the default case the iteration counts differ (67 against 90 on the first solve) and the centerline velocities differ by about 2e-5. The two only meet when the tolerance is tight enough for the iteration to converge.

`make CUDA=1 test` checks every GPU building block and full timesteps against the CPU (see [Tests](#tests)).

### Performance

Time per timestep, RK4 + FFT, Re=100, `dt = 10/n²`, VTK output off. Measured on an Intel i7-12650H and an NVIDIA GeForce RTX 3050 Laptop GPU (4 GB) under WSL2 (Ubuntu 22.04, gcc 11.4, CUDA 11.5):

| Grid | CPU, `make` | OpenMP, 8 threads | CUDA | CUDA vs CPU |
|---|---|---|---|---|
| 65×65 | 0.85 ms | 0.44 ms | 0.94 ms | 0.9× |
| 129×129 | 3.66 ms | 1.46 ms | 1.05 ms | 3.5× |
| 257×257 | 17.3 ms | 6.26 ms | 3.62 ms | 4.8× |
| 513×513 | 85.9 ms | 40.9 ms | 14.8 ms | 5.8× |
| 1025×1025 | 326 ms | 171 ms | 59.6 ms | 5.5× |

These are single runs on a laptop; repeat runs usually vary by 10–15%, occasionally by more. All builds use the Makefile's default `-O2`. Each row is a run such as:

```bash
./cnavier --n 513 --dt 3.80e-5 --tf 7.62e-3 --output-interval 0
```

Things to keep in mind:

- The GPU pays off from roughly 129×129 upwards. On small grids the fixed cost of launching kernels dominates and the CPU, especially with OpenMP, is faster.
- Much of the GPU time goes to double-precision FFTs. Consumer GeForce cards are much slower in double than in single precision, so expect larger gains on workstation/datacenter GPUs.
- Pick grid sizes where `n − 1` has only small prime factors (65, 129, 257, 513, 1025, ...). The sine transform of the interior works on length `2(n−1)`, and awkward lengths are slow on both backends: 1024×1024 takes 77 ms per step on the GPU, against 60 ms for 1025×1025.
- The table is for the FFT solver only. Gauss-Seidel and SOR are on the GPU so that every solver option works there, not because they are fast: on the default 64×64 case SOR takes about 30 ms per step on the GPU, against about 1 ms for the FFT solver.

## Configuration

All parameters are set at the top of `src/main.c`:

### Physical
| Parameter | Default | Description |
|---|---|---|
| `Re` | `100` | Reynolds number |
| `Lx`, `Ly` | `1` | Domain size |

### Numerical
| Parameter | Default | Description |
|---|---|---|
| `nx`, `ny` | `64` | Grid points in x and y |
| `dt` | `0.005` | Time step |
| `tf` | `30` | Final time |
| `order` | `6` | Finite difference order (2, 4, or 6) |
| `time_scheme` | `2` | `1` = Euler, `2` = RK4 |
| `poisson_type` | `3` | `1` = Gauss-Seidel, `2` = SOR, `3` = FFT |
| `output_interval` | `20` | Write VTK every N iterations |

### Command-line options
A few numerical parameters can be overridden without recompiling; anything not given keeps its value from `src/main.c`:

| Option | Description |
|---|---|
| `--nx N`, `--ny N` | Grid points in x and y (8 to 16384 each) |
| `--n N` | Same number of grid points in x and y |
| `--dt DT` | Time step |
| `--tf TF` | Final time |
| `--output-interval N` | Write VTK every N iterations (`0` disables VTK output) |
| `--cpu` | CUDA builds only: run on the CPU instead of the GPU |
| `--help` | Show the option list |

Both time schemes are explicit, so `dt` has to shrink with the grid spacing. Before anything is computed, `dt` is checked against two limits, and the run stops if it is above either. The message gives both limits and a `dt` that passes (the smaller limit, rounded down):

- the Courant number `u dt/dx` must not exceed 1 (`u` is the fastest wall);
- the viscous stability limit of the chosen scheme. It is computed from the actual second-derivative operator, so it follows the grid, `Re` and the finite-difference order; it scales with `dx²`.

For Euler there is a third limit, `dt ≤ 2/(Re·u²)` (that is `2ν/u²`), the stability limit of forward Euler for centered advection. It assumes the wall speed everywhere and is conservative for the cavity (at Re=1000 runs stayed stable up to about 3× it), so exceeding it only prints a warning. When an Euler run is refused, the suggested `dt` respects this limit as well, so following the suggestion does not lead to the warning.

At the defaults (64×64, Re=100, 6th order) the limits are `dt ≤ 0.00580` for RK4 and `dt ≤ 0.00416` for Euler. The benchmarks above use `dt = 10/n²`.

```bash
./cnavier --n 127 --dt 6.2e-4 --tf 20
```

Values that cannot be used (not a number, a grid outside 8–16384, `tf/dt` below one step or beyond the integer range) are rejected with an error before anything is computed or written.

The wall-clock time of the time loop is printed at the end of every run.

### Grid
The grid is nodal: node `j` sits at `x = j·Lx/(nx−1)` and node `i` at `y = i·Ly/(ny−1)`, so the first and last row and column of nodes lie on the walls. `nx` and `ny` are independent. Fields are stored as `ny` rows of `nx` values (`x` varies fastest), which is also the layout of the VTK files. Wall velocities are imposed on those nodes, ψ = 0 there for every Poisson solver, and the centerline CSVs are written at the node coordinates (interpolated onto `x = Lx/2` or `y = Ly/2` when no node lies on the centerline, i.e. for an even number of nodes).

### Boundary conditions
The default case is the **lid-driven cavity**: the top wall moves at u=1, all other walls are stationary no-slip. Boundary conditions are set via `u1`–`u4` and `v1`–`v4` in `main.c`: the left, right, bottom and top walls, in that order (`u4` is the lid).

## Output

VTK files are written to `output/` and can be opened in [ParaView](https://www.paraview.org/). Only the vorticity field is exported; other fields can be written by adding `printvtk` calls in `main.c`. Each field is its own numbered series (`vorticity-1-0.vtk`, `vorticity-1-1.vtk`, ...). The first time a run writes a series, it deletes the files of that series an earlier run left in `output/`, so a series always comes from one run; other files in `output/` are left alone.

At the end of a run the velocity profiles along the two centerlines are written to `output/centerline_u_sim.csv` and `output/centerline_v_sim.csv`, next to the Ghia et al. (1982) reference data in `centerline_*_ghia.csv`. The `_sim` files are results and are not tracked by git.

## Tests

```bash
make test
```

builds and runs `test_cnavier`:

- **Building blocks**: sparse operations against dense ones; every row of the finite-difference operators (boundary rows included) exact on the polynomials its stencil is built for, for orders 2, 4 and 6, and rows written out; the derivative operators acting along the right axis on non-square grids; wall velocities, including the corners; the VTK and centerline writers.
- **Poisson solvers**: the FFT solver against an exact eigenmode, at sizes that exercise every batching case of the transform; FFT, SOR and Gauss-Seidel giving the same answer once converged; residuals of the iterative solvers; SOR iteration counts on anisotropic grids.
- **Time stepping**: short cavity runs with both schemes on square and non-square grids (finite, divergence-free); the stability limit (0.95× runs, 1.05× diverges) and the suggested `dt`; the observed order in `dt` on the interior nodes. Euler is first order. RK4 is also first order overall, not fourth, and its error is about Euler's: the wall vorticity is computed once per step from the velocities at its start and not updated between the stages. The test checks that RK4 stays at least first order.
- **Backends**: driving the solver through `backend.c` gives exactly what calling it directly gives, and every solver ignores changes to the caller's configuration after it is created.
- **OpenMP** (with `OPENMP=1`): results on 1 and on 4 threads are bitwise identical above the size where loops go parallel.

```bash
make CUDA=1 test
```

additionally checks the GPU backend against the CPU: SpMV for all four operators, each Poisson solver, the continuity diagnostic, and full timesteps for the Euler/RK4 and FFT/SOR/Gauss-Seidel combinations on square and non-square grids. These comparisons are what keeps the GPU code in step with the CPU code, so run them on a machine with a GPU after changing either.

```bash
make OPENMP=1 CUDA=1 test
```

adds the comparison of the iterative solvers at the default tolerance, which is only meaningful when the CPU sweeps in the same red-black order as the GPU.

The exit status tells the three outcomes apart:

| Status | Meaning |
|---|---|
| `0` | every check ran and passed |
| `1` | at least one check failed |
| `77` | nothing failed, but the GPU tests could not run because no usable CUDA device was found |

A CUDA build tested on a machine without a GPU therefore does **not** count as a pass: the summary line says how many GPU tests were skipped and `make` reports an error.

### Other checks

| Command | What it checks |
|---|---|
| `make test-cli` | Invalid command lines (NaN, garbage after numbers, grids out of range, unstable `dt`, ...) are refused with exit status 1 and an error message, and write nothing. A crash counts as a failure |
| `make regression` | 100 steps of the default case against `tests/reference/`, and the steady default case within 0.004 (`u`) and 0.010 (`v`) of the Ghia et al. (1982) data stored in `tests/reference/`. Non-finite values fail; a self-test checks that |
| `make test-asan` | The test suite and a short run under AddressSanitizer and UndefinedBehaviorSanitizer |
| `make valgrind` | The test suite and a short run under valgrind; any leak, including memory still reachable at exit, fails |
| `make format-check` | Formatting against `.clang-format` (`make format` applies it) |
| `make cppcheck`, `make tidy` | Static analysis with cppcheck and clang-tidy (checks in `.clang-tidy`), including the code behind `#ifdef USE_CUDA` and `_OPENMP`. Neither reads `cudasolver.cu` |
| `make WERROR=1` | Build with `-Wextra -Werror`; with `CUDA=1`, nvcc and its host compiler treat warnings as errors too |

The results in `tests/reference/` are deterministic, so `make regression` fails only when the numerics change. If they change on purpose, regenerate them:

```bash
mkdir -p /tmp/ref/output && (cd /tmp/ref && $OLDPWD/cnavier --tf 0.5 --output-interval 0)
cp /tmp/ref/output/centerline_u_sim.csv tests/reference/centerline_u_short.csv
cp /tmp/ref/output/centerline_v_sim.csv tests/reference/centerline_v_short.csv
```

### Continuous integration

GitHub Actions (`.github/workflows/ci.yml`) runs, in order: formatting; cppcheck, clang-tidy and `-Wextra -Werror` builds with gcc and clang; serial, OpenMP and CUDA builds; then the test suite, the command-line and regression tests for the serial and OpenMP builds, and the sanitizer and valgrind checks. The CUDA build's test program must report its GPU tests as skipped (exit status 77).

The hosted runners have no GPU, so the GPU-vs-CPU checks only run where someone runs `make CUDA=1 test` on a machine with one. The pull request template asks for that whenever the numerics change.

## Project structure

```
cnavier/
├── .github/workflows/
│   └── ci.yml          # Formatting, static analysis, builds, tests, memory checks
├── src/
│   ├── main.c          # Configuration, command line, time loop and output
│   ├── backend.c       # One interface over the CPU and CUDA solvers
│   ├── linearalg.c     # Dense and sparse (CSR) linear algebra
│   ├── finitediff.c    # Finite difference operators (dense + sparse)
│   ├── fluiddyn.c      # CPU timestep: Euler/RK4, wall and vorticity BCs, stability limit
│   ├── poisson.c       # Gauss-Seidel, SOR, and FFT Poisson solvers
│   ├── cudasolver.cu   # CUDA backend (built only with CUDA=1)
│   └── utils.c         # VTK output, random utilities
├── include/
│   ├── linearalg.h
│   ├── finitediff.h
│   ├── fluiddyn.h      # solver_config (shared by both backends) and the CPU solver
│   ├── backend.h
│   ├── poisson.h
│   ├── cudasolver.h
│   └── utils.h
├── tests/
│   ├── test_solver.c   # Unit and solver tests, GPU-vs-CPU checks
│   ├── cli.sh          # Command-line tests
│   ├── regression.sh   # Stored-result and Ghia et al. checks
│   └── reference/      # Stored results and the Ghia et al. data
├── output/             # VTK output files
├── Re1000_cavity_flow_example.png
├── Re1000_cavity_flow_example.mp4
├── LICENSE
├── Makefile
└── README.md
```

## Notes on the FFT Poisson solver

The default `poisson_type = 3` uses FFTW3's `RODFT00` plan (DST-I) of the interior nodes to solve the Poisson equation exactly in O(n² log n). It solves the same discrete problem as Gauss-Seidel and SOR (the 5-point Laplacian with ψ = 0 on the wall nodes), so once those have converged all three give the same ψ. This is orders of magnitude faster than the iterative solvers at high Reynolds numbers where many iterations are required for convergence. FFTW plans are computed once at startup via `fft_setup()` and reused every timestep.

With RK4 (`time_scheme = 2`), the Poisson equation is solved once per RK4 stage (5 solves per timestep total). Since each solve is O(n² log n), this remains fast.
