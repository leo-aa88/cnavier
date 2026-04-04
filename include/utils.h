// Utilities library

#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include "linearalg.h"

double randdouble(double min, double max); // generates random double number between min and max
void printvtk(mtrx A, char *title);        // prints matrix A to a vtk file

// Write centerline velocity profiles to CSV for validation against Ghia et al. (1982).
// u_centerline: u-velocity along vertical centerline (x=0.5), sampled at each j
// v_centerline: v-velocity along horizontal centerline (y=0.5), sampled at each i
void print_centerline(mtrx u, mtrx v, int nx, int ny, double dx, double dy);

#endif
