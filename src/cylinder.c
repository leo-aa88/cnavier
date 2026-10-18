// Flow past a circular cylinder on a polar grid (cylinder.h)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "cylinder.h"
#include "polar.h"
#include "linearalg.h"
#include "finitediff.h"
#include "fluiddyn.h"
#include "backend.h"
#include "poisson.h"

// The reference solution at Re = 40 (Gautier et al. 2013, Tables 1 and 2 and
// Appendix A; tests/reference/cylinder_re40.csv)
#define REF_CD      1.4931
#define REF_LW      2.2360
#define REF_THETA_S 126.3945
#define REF_A       0.71
#define REF_B       0.59
#define REF_UX33    1.0670111347099145
#define REF_UY33    0.015134054077428072

// Radius of the cylinder: D = 1
#define CYL_R 0.5

static FILE *open_output(const char *name)
{
    FILE *f = fopen(name, "w");
    if (!f)
    {
        printf("\nError while opening file\n");
        exit(1);
    }
    return f;
}

// Largest |a - b| over the interior rows (1..nr-2)
static double interior_change(const mtrx *a, const mtrx *b)
{
    int k;
    double e = 0.0;
    for (k = a->n; k < (a->m - 1) * a->n; k++)
        e = fmax(e, fabs(a->M[k] - b->M[k]));
    return e;
}

int cylinder_run(const cylinder_options *o)
{
    int i, j, t, nt = o->ntheta;
    double R = CYL_R, dt = o->dt;
    polar_grid g = polar_grid_make(nt, o->nr, R, o->r_inf);
    int nr = g.nr;

    if (o->outer_bc == POLAR_FARFIELD)
        polar_farfield_set(&g, 1.5, o->Re); // C_D is updated as the run finds it
    else
        polar_freestream_set(&g);

    // 1-D operators: theta periodic, xi with one-sided rows at both ends
    smtrx d_t = SDiff1_periodic(nt, o->order, g.dtheta), d_t2 = SDiff2_periodic(nt, o->order, g.dtheta);
    smtrx d_x = SDiff1(nr, o->order, g.dxi), d_x2 = SDiff2(nr, o->order, g.dxi);

    // Time step: the viscous term (w_xixi + w_thth) / (Re r^2) is stiffest on
    // the cylinder, so the limits are those of a uniform grid of spacing
    // dtheta, dxi at Reynolds number Re R^2; the Courant number uses the
    // smallest cell, R min(dxi, dtheta), and twice the free stream (the peak
    // of potential flow)
    double h = R * fmin(g.dxi, g.dtheta);
    dt_limits lim = time_step_limits(&d_t2, &d_x2, h, o->Re * R * R, 2.0, 1.0, o->time_scheme);
    if (dt <= 0.0) dt = round_down_3(0.9 * lim.suggest);
    if (dt > lim.accept)
    {
        printf("** Error: dt = %g is too large; use --dt %.3g or less (Courant limit %.3g, stability limit of the "
               "viscous term on the cylinder %.3g) **\n",
               dt, round_down_3(lim.suggest), round_down_3(lim.courant), round_down_3(lim.viscous));
        polar_grid_free(&g);
        freesm(d_t);
        freesm(d_t2);
        freesm(d_x);
        freesm(d_x2);
        return 1;
    }

    smtrx It = seye(nt), Ix = seye(nr);
    smtrx DX = skronecker(Ix, d_t), DY = skronecker(d_x, It);
    smtrx DX2 = skronecker(Ix, d_t2), DY2 = skronecker(d_x2, It);
    freesm(It);
    freesm(Ix);

    solver_config cfg = {0};
    cfg.nx = nt;
    cfg.ny = nr;
    cfg.dx = g.dtheta;
    cfg.dy = g.dxi;
    cfg.Re = o->Re;
    cfg.dt = dt;
    cfg.time_scheme = o->time_scheme;
    cfg.poisson_type = 3;
    cfg.poisson_order = o->poisson_order;
    cfg.wall_closure = o->wall_closure;
    cfg.DX = &DX;
    cfg.DY = &DY;
    cfg.DX2 = &DX2;
    cfg.DY2 = &DY2;
    cfg.polar = &g;

    printf("Case: flow past a cylinder (polar grid) | Re: %g\n", o->Re);
    printf("Grid: %d (theta) x %d (xi), r from %g to %g | cells: dtheta %.4g, dxi %.4g, on the wall %.4g x %.4g\n", nt,
           nr, R, g.r_inf, g.dtheta, g.dxi, R * g.dtheta, R * (exp(g.dxi) - 1.0));
    printf("Derivatives: order %d | Poisson: order %d in xi | wall vorticity: %s | outer boundary: %s\n", o->order,
           o->poisson_order, o->wall_closure ? "Briley" : "velocity",
           o->outer_bc == POLAR_FARFIELD ? "far field (Schlichting wake + mass correction)" : "free stream");
    printf("dt: %g (limits: Courant %.3g, viscous %.3g) | tf: %g | %s\n", dt, lim.courant, lim.viscous, o->tf,
           o->time_scheme == 1 ? "Euler" : "RK4");

    backend *solver = backend_create(&cfg, o->use_gpu);
    printf("Backend: %s\n", backend_name(solver));

    // Impulsive start from rest: w = 0, so the first solve gives potential
    // flow and the wall closure the vortex sheet on the cylinder
    mtrx *u, *v, *w;
    mtrx w_prev = initm(nr, nt);
    backend_fields(solver, &u, &v, &w);

    FILE *hist = open_output("./output/cylinder_history.csv");
    fprintf(hist, "step,t,cd,cdp,cdf,cl,wake_length,separation_angle,dwdt_max,cd_boundary\n");

    int it_max = (int)(o->tf / dt), steps = 0, steady = 0;
    double residual = HUGE_VAL;
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (t = 1; t <= it_max; t++)
    {
        backend_step(solver);
        steps = t;
        if (o->output_interval > 0 && t % o->output_interval == 0)
        {
            backend_fields(solver, NULL, NULL, &w);
            printvtk_polar(&g, w, NULL, NULL, "vorticity");
        }
        if (t % o->diag_interval != 0 && t != it_max) continue;

        backend_fields(solver, &u, &v, &w);
        int k, n = nt * nr, finite = 1;
        for (k = 0; k < n && finite; k++)
            finite = isfinite(w->M[k]);
        if (!finite)
        {
            printf("** Error: the solution is not finite at step %d (t = %g); reduce dt **\n", t, t * dt);
            fclose(hist);
            freem(&w_prev);
            backend_free(solver);
            polar_grid_free(&g);
            return 1;
        }
        residual = interior_change(w, &w_prev) / (o->diag_interval * dt);
        mtrxcpy(w_prev, *w);
        polar_forces f = polar_wall_forces(&g, &DY, w->M, o->Re, NULL);
        double lw = polar_wake_length(&g, v->M), ts = polar_separation_angle(&g, w->M);
        fprintf(hist, "%d,%.10g,%.12g,%.12g,%.12g,%.6e,%.10g,%.10g,%.6e,%.10g\n", t, t * dt, f.cd, f.cdp, f.cdf, f.cl,
                lw, ts, residual, g.cd);
        fflush(hist);
        printf("Step %d | t = %.3f | C_D = %.6f (pressure %.6f, friction %.6f) | L_w/D = %.4f | theta_s = %.3f | "
               "max |dw/dt| = %.3e\n",
               t, t * dt, f.cd, f.cdp, f.cdf, lw, ts, residual);

        // The far-field condition carries the drag of the body: follow the
        // run's C_D, under-relaxed, once the impulsive start has passed
        if (o->outer_bc == POLAR_FARFIELD && t * dt >= 5.0 && f.cd > 0.0 && f.cd < 10.0)
            polar_farfield_set(&g, g.cd + 0.5 * (f.cd - g.cd), o->Re);

        if (o->steady_tol > 0.0 && residual < o->steady_tol && t * dt >= 10.0)
        {
            steady = 1;
            break;
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    fclose(hist);
    double elapsed = (double)(t1.tv_sec - t0.tv_sec) + 1E-9 * (double)(t1.tv_nsec - t0.tv_nsec);

    // Final state: velocity with the wall values, psi solved again for the
    // diagnostics.
    // Continuity of the velocity of the Poisson solution, before the wall
    // velocity is imposed (as the Cartesian run measures it)
    double cmax, cmin;
    backend_continuity(solver, &cmax, &cmin);
    backend_fields(solver, &u, &v, &w);
    polar_apply_bc(&g, *u, *v);
    mtrx psi = initm(nr, nt);
    polar_poisson *ps = polar_poisson_setup(&g, &DX2, o->poisson_order);
    polar_poisson_solve(ps, w->M, g.psi_out, psi.M);
    polar_poisson_free(ps);

    double *p_wall = (double *)malloc((size_t)nt * sizeof(double));
    if (!p_wall)
    {
        printf("** Error: insufficient memory **\n");
        exit(1);
    }
    polar_forces f = polar_wall_forces(&g, &DY, w->M, o->Re, p_wall);
    double lw = polar_wake_length(&g, v->M), ts = polar_separation_angle(&g, w->M), a = 0.0, b = 0.0, ux, uy, w33;
    int centre = polar_vortex_centre(&g, psi.M, &a, &b);
    int inside = polar_interpolate(&g, u->M, v->M, w->M, 3.0, 3.0, &ux, &uy, &w33);

    FILE *fs = open_output("./output/cylinder_surface.csv");
    fprintf(fs, "theta_deg,w_wall,p_wall\n");
    for (j = 0; j < nt; j++)
        fprintf(fs, "%.10g,%.12g,%.12g\n", j * g.dtheta * 180.0 / PI, MAt(*w, 0, j), p_wall[j]);
    fclose(fs);
    FILE *fc = open_output("./output/cylinder_centreline.csv");
    fprintf(fc, "x,u_x\n");
    if (nt % 2 == 0)
        for (i = nr - 1; i >= 0; i--)
            fprintf(fc, "%.12g,%.12g\n", -g.r[i], -MAt(*v, i, nt / 2)); // upstream, theta = pi
    for (i = 0; i < nr; i++)
        fprintf(fc, "%.12g,%.12g\n", g.r[i], MAt(*v, i, 0)); // the wake, theta = 0
    fclose(fc);
    printvtk_polar(&g, w, u, v, "vorticity");
    printvtk_polar(&g, &psi, NULL, NULL, "streamfunction");

    printf("\n%s at t = %g after %d steps (max |dw/dt| = %.3e)\n", steady ? "Steady" : "Stopped", steps * dt, steps,
           residual);
    printf("Continuity max: %E | min: %E | asymmetry max |w(theta) + w(-theta)| = %.3e\n", cmax, cmin,
           polar_asymmetry(&g, w->M));
    printf("Mean dp/dtheta on the wall: %.3e (0 for a consistent pressure)\n", f.dpdt_mean);
    printf("Comparison with Gautier et al. (2013), Re = 40 (their 200 x 1024 grid, r_inf = 40):\n");
    printf("  C_D       %10.5f   reference %8.4f   (%+.2f%%); pressure %.5f, friction %.5f, C_L %.2e\n", f.cd, REF_CD,
           100.0 * (f.cd - REF_CD) / REF_CD, f.cdp, f.cdf, f.cl);
    printf("  L_w/D     %10.5f   reference %8.4f   (%+.2f%%)\n", lw, REF_LW, 100.0 * (lw - REF_LW) / REF_LW);
    printf("  theta_s   %10.4f   reference %8.4f   (%+.3f deg)\n", ts, REF_THETA_S, ts - REF_THETA_S);
    if (centre)
        printf("  a/D, b/D  %6.4f %6.4f  reference %.2f %.2f\n", a, b, REF_A, REF_B);
    else
        printf("  a/D, b/D  no recirculation\n");
    if (inside)
        printf("  at (3, 3) u_x %.6f (reference %.6f, %+.1e), u_y %.6f (reference %.6f, %+.1e), w %.2e\n", ux,
               REF_UX33, ux - REF_UX33, uy, REF_UY33, uy - REF_UY33, w33);
    printf("Wall-clock time: %.3lf s total | %.4lf ms per step (%d steps, %s)\n", elapsed,
           steps ? 1E3 * elapsed / steps : 0.0, steps, backend_name(solver));

    free(p_wall);
    freem(&psi);
    freem(&w_prev);
    backend_free(solver);
    freesm(DX);
    freesm(DY);
    freesm(DX2);
    freesm(DY2);
    freesm(d_t);
    freesm(d_t2);
    freesm(d_x);
    freesm(d_x2);
    polar_grid_free(&g);
    printf("Simulation complete!\n");
    return 0;
}
