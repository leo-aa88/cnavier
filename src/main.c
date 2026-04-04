#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "linearalg.h"
#include "finitediff.h"
#include "utils.h"
#include "poisson.h"
#include "fluiddyn.h"

int main(int argc, char *argv[])
{
    int i, j, t;

    // Physical parameters
    double Re = 1000.;
    int Lx = 1;
    int Ly = 1;

    // Numerical parameters
    int nx = 64;
    int ny = 64;
    double dt = 0.005;
    double tf = 20;
    double max_co = 1.;
    int order = 6;
    int poisson_max_it = 10000;
    double poisson_tol = 1E-3;
    int output_interval = 10;
    int poisson_type = 3; // 1=Gauss-Seidel  2=SOR  3=FFT (direct, exact)
    double rho  = 0.5 * (cos(PI / nx) + cos(PI / ny)); // spectral radius of Gauss-Seidel
    double beta  = 2.0 / (1.0 + sqrt(1.0 - rho * rho));  // optimal SOR parameter

    printf("Poisson SOR parameter: %lf\n", beta);

    if (poisson_type == 3) fft_setup(nx, ny);

    // Boundary conditions (Dirichlet)
    double ui = 0., vi = 0.;
    double u1 = 0., u2 = 0., u3 = 0., u4 = 1.;
    double v1 = 0., v2 = 0., v3 = 0., v4 = 0.;

    // Cell sizes
    double dx = (double)Lx / nx;
    double dy = (double)Ly / ny;

    // Build sparse 1D operators then free them after Kronecker
    smtrx sd_x  = SDiff1(nx, order, dx);
    smtrx sd_y  = SDiff1(ny, order, dy);
    smtrx sd_x2 = SDiff2(nx, order, dx);
    smtrx sd_y2 = SDiff2(ny, order, dy);
    smtrx sIx   = seye(nx);
    smtrx sIy   = seye(ny);

    // Sparse 2D operators: DX = I_y x d_x,  DY = d_y x I_x
    smtrx DX  = skronecker(sIy,   sd_x);
    smtrx DY  = skronecker(sd_y,  sIx);
    smtrx DX2 = skronecker(sIy,   sd_x2);
    smtrx DY2 = skronecker(sd_y2, sIx);

    freesm(sd_x); freesm(sd_y); freesm(sd_x2); freesm(sd_y2);
    freesm(sIx);  freesm(sIy);

    int N = nx * ny;

    // Courant check (u4 is the only non-zero BC velocity)
    int it_max = (int)((tf / dt) - 1);
    double r1 = u4 * dt / dx;
    double r2 = u4 * dt / dy;
    if ((r1 > max_co) || (r2 > max_co))
    {
        printf("Unstable Solution! r1=%lf r2=%lf\n", r1, r2);
        exit(1);
    }

    // Dense field matrices
    mtrx u   = initm(nx, ny);
    mtrx v   = initm(nx, ny);
    mtrx w   = initm(nx, ny);
    mtrx psi         = initm(nx, ny);
    mtrx psi_scratch = initm(nx, ny); // scratch buffer for Poisson solver

    // Derivative matrices — pre-allocated once, reused every iteration
    mtrx dwdx   = initm(nx, ny);
    mtrx dwdy   = initm(nx, ny);
    mtrx d2wdx2 = initm(nx, ny);
    mtrx d2wdy2 = initm(nx, ny);
    mtrx dpsidx = initm(nx, ny);
    mtrx dpsidy = initm(nx, ny);
    mtrx dudy   = initm(nx, ny);
    mtrx dvdx   = initm(nx, ny);
    mtrx dudx   = initm(nx, ny);
    mtrx dvdy   = initm(nx, ny);
    mtrx check_continuity = initm(nx, ny);

    // Flat work buffers — pre-allocated once, no malloc/free in the time loop
    double *flat_u   = (double *)malloc(N * sizeof(double));
    double *flat_v   = (double *)malloc(N * sizeof(double));
    double *flat_w   = (double *)malloc(N * sizeof(double));
    double *flat_psi = (double *)malloc(N * sizeof(double));
    double *flat_tmp = (double *)malloc(N * sizeof(double));
    if (!flat_u || !flat_v || !flat_w || !flat_psi || !flat_tmp)
    {
        printf("** Error: insufficient memory **\n");
        exit(1);
    }

    // Initial condition
    for (i = 1; i < nx - 1; i++)
        for (j = 1; j < ny - 1; j++)
        {
            MAt(u, i, j) = ui;
            MAt(v, i, j) = vi;
        }

    // Main time loop
    for (t = 0; t <= it_max; t++)
    {
        // Boundary conditions
        for (j = 0; j < ny; j++)
        {
            MAt(u, 0, j)    = u3;  MAt(v, 0, j)    = v3;
            MAt(u, nx-1, j) = u4;  MAt(v, nx-1, j) = v4;
        }
        for (i = 0; i < nx; i++)
        {
            MAt(u, i, 0)    = u1;  MAt(v, i, 0)    = v1;
            MAt(u, i, ny-1) = u2;  MAt(v, i, ny-1) = v2;
        }

        // Vorticity BCs: w = dv/dx - du/dy evaluated at boundaries
        flatten(u, flat_u, nx, ny);
        flatten(v, flat_v, nx, ny);
        spmv(DY, flat_u, flat_tmp);  unflatten(flat_tmp, dudy, nx, ny);
        spmv(DX, flat_v, flat_tmp);  unflatten(flat_tmp, dvdx, nx, ny);

        for (j = 0; j < ny; j++)
        {
            MAt(w, 0, j)    = MAt(dvdx, 0, j)    - MAt(dudy, 0, j);
            MAt(w, nx-1, j) = MAt(dvdx, nx-1, j) - MAt(dudy, nx-1, j);
        }
        for (i = 0; i < nx; i++)
        {
            MAt(w, i, 0)    = MAt(dvdx, i, 0)    - MAt(dudy, i, 0);
            MAt(w, i, ny-1) = MAt(dvdx, i, ny-1) - MAt(dudy, i, ny-1);
        }

        // Vorticity derivatives via SpMV (no reshape, no malloc)
        flatten(w, flat_w, nx, ny);
        spmv(DX,  flat_w, flat_tmp);  unflatten(flat_tmp, dwdx,   nx, ny);
        spmv(DY,  flat_w, flat_tmp);  unflatten(flat_tmp, dwdy,   nx, ny);
        spmv(DX2, flat_w, flat_tmp);  unflatten(flat_tmp, d2wdx2, nx, ny);
        spmv(DY2, flat_w, flat_tmp);  unflatten(flat_tmp, d2wdy2, nx, ny);

        // Time advancement (Euler)
        euler(w, dwdx, dwdy, d2wdx2, d2wdy2, u, v, Re, dt);

        // Poisson solve: nabla^2 psi = -w

        invsig(w);
        if (poisson_type == 1)
            poisson(w, psi, psi_scratch, dx, dy, poisson_max_it, poisson_tol);
        else if (poisson_type == 2)
            poisson_SOR(w, psi, psi_scratch, dx, dy, poisson_max_it, poisson_tol, beta);
        else
            poisson_FFT(w, psi, dx, dy);
        invsig(w);

        // Velocities from stream function: u = dpsi/dy, v = -dpsi/dx
        flatten(psi, flat_psi, nx, ny);
        spmv(DY, flat_psi, flat_tmp);  unflatten(flat_tmp, dpsidy, nx, ny);
        spmv(DX, flat_psi, flat_tmp);  unflatten(flat_tmp, dpsidx, nx, ny);

        mtrxcpy(u, dpsidy);
        invsig(dpsidx);
        mtrxcpy(v, dpsidx);

        // Continuity check: du/dx + dv/dy ~ 0
        flatten(u, flat_u, nx, ny);
        flatten(v, flat_v, nx, ny);
        spmv(DX, flat_u, flat_tmp);  unflatten(flat_tmp, dudx, nx, ny);
        spmv(DY, flat_v, flat_tmp);  unflatten(flat_tmp, dvdy, nx, ny);

        // reuse check_continuity storage
        for (i = 0; i < nx; i++)
            for (j = 0; j < ny; j++)
                MAt(check_continuity, i, j) = MAt(dudx, i, j) + MAt(dvdy, i, j);

        printf("Iteration: %d | Time: %.4lf | Progress: %.2lf%%\n",
               t, (double)t * dt, (double)100 * t / it_max);
        printf("Continuity max: %E | min: %E\n",
               maxel(check_continuity), minel(check_continuity));

        if (t % output_interval == 0)
            printvtk(w, "vorticity");
    }

    // Free dense fields
    freem(&u);
    freem(&v);
    freem(&w);
    freem(&psi);

    // Free derivative matrices
    freem(&dwdx);   freem(&dwdy);
    freem(&d2wdx2); freem(&d2wdy2);
    freem(&dpsidx); freem(&dpsidy);
    freem(&dudy);   freem(&dvdx);
    freem(&dudx);   freem(&dvdy);
    freem(&check_continuity);

    // Free flat buffers
    free(flat_u); free(flat_v); free(flat_w);
    free(flat_psi); free(flat_tmp);

    // Free sparse operators
    freesm(DX); freesm(DY); freesm(DX2); freesm(DY2);

    if (poisson_type == 3) fft_cleanup();

    printf("Simulation complete!\n");
    return 0;
}
