// Flow past a circular cylinder on a polar grid (--geometry polar): the run
// and its diagnostics, compared with the reference solution of Gautier,
// Biau and Lamballais (2013) at Re = 40

#ifndef CYLINDER_H_INCLUDED
#define CYLINDER_H_INCLUDED

typedef struct
{
    int ntheta, nr;      // grid points in theta and xi (nr <= 0: square cells)
    double r_inf;        // outer radius, in diameters D = 1 (cylinder radius 1/2)
    double Re;           // U D / nu
    double dt;           // time step; <= 0: chosen from the stability limits
    double tf;           // final time
    int order;           // finite differences: 2, 4 or 6
    int poisson_order;   // xi operator of the Poisson solver: 2 or 4 (compact)
    int wall_closure;    // 0: velocity, 1: Briley's formula from psi
    int time_scheme;     // 1 Euler, 2 RK4
    int outer_bc;        // POLAR_FREESTREAM or POLAR_FARFIELD (polar.h)
    double steady_tol;   // stop when max |dw/dt| falls below it (0: never)
    int diag_interval;   // steps between diagnostics (history, steady check)
    int output_interval; // steps between VTK frames (0: the final state only)
    int use_gpu;         // asked for the GPU (the polar grid runs on the CPU)
} cylinder_options;

// Run the case and write output/cylinder_*.csv and VTK frames; returns the
// exit status
int cylinder_run(const cylinder_options *o);

#endif // CYLINDER_H_INCLUDED
