// Utilities library

#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include "linearalg.h"

double randdouble(double min, double max); // generates random double number between min and max
// Write A (ny rows of nx values) to output/<title>-1-<n>.vtk as a structured
// grid with spacing dx, dy
void printvtk(mtrx A, char *title, double dx, double dy);

// Memory the process can allocate now, in bytes: the kernel's MemAvailable,
// or the room left under a cgroup memory limit if that is smaller. 0 means
// none (a cgroup at or above its limit); -1 means it cannot be determined.
double available_memory(void);

// Write centerline velocity profiles to CSV for validation against Ghia et al. (1982).
// u along the vertical centerline (x = Lx/2) at every node y = i*dy, and v along
// the horizontal centerline (y = Ly/2) at every node x = j*dx. Values are
// interpolated linearly onto the centerline when no node lies on it.
void print_centerline(mtrx u, mtrx v, int nx, int ny, double dx, double dy);

#endif
