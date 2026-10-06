// Default number of OpenMP threads

#ifndef THREADS_H_INCLUDED
#define THREADS_H_INCLUDED

// Number of distinct physical cores among the CPUs the calling thread may run
// on, from Linux sysfs. 0 if it cannot be determined, or without OpenMP.
int physical_cores(void);

// Unless OMP_NUM_THREADS or an OpenMP thread placement is set, use one thread
// per physical core (see threads.c). Does nothing without OpenMP.
void default_threads(void);

#endif
