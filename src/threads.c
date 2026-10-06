#ifdef __linux__
#define _GNU_SOURCE // NOLINT(bugprone-reserved-identifier): needed for sched_getaffinity
#include <sched.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "threads.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(_OPENMP) && defined(__linux__)
// Number of distinct physical cores among the CPUs this process may run on.
// CPUs that are hardware threads of one core share the same sysfs list of
// sibling CPUs (topology/core_cpus_list, thread_siblings_list on kernels
// before 5.7), so this counts distinct lists. That works on every
// architecture and for hybrid CPUs, whose cores have different numbers of
// hardware threads. The mask is the calling thread's, which is the process's
// as long as the OpenMP runtime has not bound the initial thread; see
// default_threads(). Returns 0 if it cannot be determined.
int physical_cores(void)
{
    cpu_set_t set;
    int cpu, k, count = 0;
    char(*seen)[64];

    if (sched_getaffinity(0, sizeof(set), &set) != 0) return 0;
    seen = malloc(sizeof(*seen) * CPU_SETSIZE);
    if (!seen) return 0;
    for (cpu = 0; cpu < CPU_SETSIZE; cpu++)
    {
        char path[128], list[64] = "";
        FILE *f;

        if (!CPU_ISSET(cpu, &set)) continue;
        snprintf(path, sizeof(path), "/sys/devices/system/cpu/cpu%d/topology/core_cpus_list", cpu);
        if (!(f = fopen(path, "r")))
        {
            snprintf(path, sizeof(path), "/sys/devices/system/cpu/cpu%d/topology/thread_siblings_list", cpu);
            f = fopen(path, "r");
        }
        if (!f || !fgets(list, sizeof(list), f) || list[0] == '\0')
        {
            if (f) fclose(f);
            count = 0;
            break;
        }
        fclose(f);
        for (k = 0; k < count; k++)
            if (strcmp(seen[k], list) == 0) break;
        if (k == count) snprintf(seen[count++], sizeof(seen[0]), "%s", list);
    }
    free(seen);
    return count;
}

#else
// Elsewhere the topology is not read, and the OpenMP default applies
int physical_cores(void) { return 0; }
#endif

#ifdef _OPENMP
// Unless OMP_NUM_THREADS says otherwise, use one thread per physical core the
// process may run on, instead of the OpenMP default of one per logical CPU.
// The sparse products and transforms are limited by memory bandwidth, so a
// second hardware thread per core adds little, and every extra thread is one
// more to wait for at each of the ~70 barriers per step; on a busy machine
// that makes surplus threads very costly.
// If thread placement is configured (OMP_PROC_BIND, OMP_PLACES or
// GOMP_CPU_AFFINITY), the OpenMP default is left alone: the user has chosen
// the placement, and the runtime has already bound this thread to one place,
// so its affinity mask no longer describes the process.
void default_threads(void)
{
    int cores;
    if (getenv("OMP_NUM_THREADS") || omp_get_proc_bind() != omp_proc_bind_false) return;
    cores = physical_cores();
    if (cores > 0 && cores < omp_get_max_threads()) omp_set_num_threads(cores);
}
#else
void default_threads(void) {}
#endif
