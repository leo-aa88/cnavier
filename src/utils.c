#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <dirent.h>
#include <unistd.h>
#include "linearalg.h"
#include "utils.h"

// Read the first number on the line of `path` that starts with `key` (the
// whole file if key is ""). Returns -1 if there is none.
static double read_number(const char *path, const char *key)
{
    char line[256];
    double val = -1.0;
    FILE *f = fopen(path, "r");

    if (!f) return -1.0;
    while (fgets(line, sizeof(line), f))
        if (strncmp(line, key, strlen(key)) == 0)
        {
            char *end;
            double v = strtod(line + strlen(key), &end);
            if (end != line + strlen(key)) val = v;
            break;
        }
    fclose(f);
    return val;
}

// Room left under the memory limit of cgroup directory `dir`, in bytes.
// Returns 0 if the cgroup has no limit. Page cache on the file LRU lists
// (active_file + inactive_file) counts as free: the kernel reclaims both
// before it fails an allocation. Shared memory and tmpfs are charged as file
// memory too, but they sit on the anon lists and are not reclaimable without
// swap; the two LRU counters leave them out, so they stay counted as used.
static int cgroup_room(const char *dir, int v2, double *room)
{
    char path[4300];
    double limit, used, active, inactive;

    snprintf(path, sizeof(path), "%s/%s", dir, v2 ? "memory.max" : "memory.limit_in_bytes");
    limit = read_number(path, ""); // v2 "max" reads as -1, i.e. no limit
    if (limit <= 0 || limit >= 1E18) return 0;
    snprintf(path, sizeof(path), "%s/%s", dir, v2 ? "memory.current" : "memory.usage_in_bytes");
    used = read_number(path, "");
    snprintf(path, sizeof(path), "%s/memory.stat", dir);
    active   = read_number(path, v2 ? "active_file " : "total_active_file ");
    inactive = read_number(path, v2 ? "inactive_file " : "total_inactive_file ");
    if (used < 0) used = 0.0;
    if (active > 0 && active < used) used -= active;
    if (inactive > 0 && inactive < used) used -= inactive;
    *room = limit > used ? limit - used : 0.0;
    return 1;
}

// The process's cgroup, relative to the hierarchy root, from
// /proc/self/cgroup. Returns 1 for the cgroup v1 memory controller's line,
// 2 for the cgroup v2 "0::" line, 0 if there is neither. The v1 line wins:
// on hybrid systems the v2 hierarchy has no memory controller.
static int own_cgroup(char *cg, size_t size)
{
    char line[4096];
    int version = 0;
    FILE *f = fopen("/proc/self/cgroup", "r");

    if (!f) return 0;
    while (version != 1 && fgets(line, sizeof(line), f))
    {
        char *ctrl = strchr(line, ':');
        char *path = ctrl ? strchr(ctrl + 1, ':') : NULL;
        int found = 0;

        if (!path) continue;
        *ctrl++ = '\0';
        *path++ = '\0';
        path[strcspn(path, "\n")] = '\0';
        if (strcmp(line, "0") == 0 && *ctrl == '\0') found = 2;
        else
            for (char *c = strtok(ctrl, ","); c; c = strtok(NULL, ","))
                if (strcmp(c, "memory") == 0) found = 1;
        if (found)
        {
            version = found;
            snprintf(cg, size, "%s", path);
        }
    }
    fclose(f);
    return version;
}

double available_memory(void)
{
    char cg[4096], dir[4200];
    double avail = -1.0;
    int version;

    // What the kernel thinks can be allocated without swapping
    double kb = read_number("/proc/meminfo", "MemAvailable:");
    if (kb > 0) avail = kb * 1024.0;
    else
    {
        long pages = sysconf(_SC_PHYS_PAGES), size = sysconf(_SC_PAGESIZE);
        if (pages > 0 && size > 0) avail = (double)pages * (double)size;
    }

    // Memory limits on the process's cgroup and on each of its ancestors
    // (containers, systemd scopes, batch jobs): the tightest one wins. Inside
    // a container with its own cgroup namespace the path is "/".
    version = own_cgroup(cg, sizeof(cg));
    if (version == 0) return avail;
    for (;;)
    {
        char *slash;
        double room;

        snprintf(dir, sizeof(dir), "/sys/fs/cgroup%s%s", version == 1 ? "/memory" : "",
                 strcmp(cg, "/") == 0 ? "" : cg);
        if (cgroup_room(dir, version == 2, &room) && (avail < 0 || room < avail)) avail = room;
        slash = strrchr(cg, '/');
        if (!slash || strcmp(cg, "/") == 0) break;
        if (slash == cg) cg[1] = '\0'; // last step: the root itself
        else *slash = '\0';
    }
    return avail;
}

double randdouble(double min, double max)
{
    double range = (max - min);
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

// Delete output/<title>-1-<number>.vtk, the series an earlier run left behind.
// A shorter run would otherwise leave the earlier run's later frames in place,
// and ParaView would show them as part of this run's series.
static void remove_old_series(const char *title)
{
    char prefix[64], path[512];
    size_t len;
    struct dirent *e;
    DIR *d = opendir("./output");

    if (!d) return;
    snprintf(prefix, sizeof(prefix), "%s-1-", title);
    len = strlen(prefix);
    while ((e = readdir(d)) != NULL)
    {
        const char *rest;
        size_t digits;
        if (strncmp(e->d_name, prefix, len) != 0) continue;  // also skips names shorter than prefix
        rest = e->d_name + len;
        digits = strspn(rest, "0123456789");
        if (digits > 0 && strcmp(rest + digits, ".vtk") == 0)
        {
            snprintf(path, sizeof(path), "./output/%s", e->d_name);
            unlink(path);
        }
    }
    closedir(d);
}

// Next frame number of the series `title`. The first call for a title starts
// its series at 0 and deletes what an earlier run left of it.
static int next_frame(const char *title)
{
    static char titles[16][64];
    static int frames[16], n = 0;
    int k;

    for (k = 0; k < n; k++)
        if (strcmp(titles[k], title) == 0) return frames[k]++;
    if (n == 16)
    {
        printf("\n** Error: more than 16 VTK series **\n");
        exit(1);
    }
    snprintf(titles[n], sizeof(titles[n]), "%s", title);
    frames[n] = 1;
    n++;
    remove_old_series(title);
    return 0;
}

void printvtk(mtrx A, char *title)
{
    int i, j, count;
    char c[320];
    char name[64];
    FILE *pf;

    if (A.M == NULL)
    {
        printf("\n** Error: Aborting program **\n");
        exit(1);
    }
    if ((A.m < 1) || (A.n < 1))
    {
        printf("\n** Error: Invalid parameter **\n");
        exit(1);
    }

    count = next_frame(title);
    snprintf(name, sizeof(name), "./output/%s-1-%d.vtk", title, count);

    if ((pf = fopen(name, "w")) == NULL)
    {
        printf("\nError while opening file\n");
        exit(1);
    }

    printf("%s\n", name);

    fprintf(pf, "# vtk DataFile Version 2.0\n"); // vtk file headers
    fprintf(pf, "test\n");
    fprintf(pf, "ASCII\n");
    fprintf(pf, "DATASET STRUCTURED_POINTS\n");
    fprintf(pf, "DIMENSIONS %d %d 1\n", A.m, A.n);
    fprintf(pf, "ORIGIN 0 0 0\n");
    fprintf(pf, "SPACING 1 1 1\n");
    fprintf(pf, "POINT_DATA %d\n", A.m * A.n);
    fprintf(pf, "SCALARS values float\n");
    fprintf(pf, "LOOKUP_TABLE default");

    for (i = 0; i < A.m; i++)
    {
        fprintf(pf, "\n");
        for (j = 0; j < A.n; j++)
        {
            if ((j == 0))
            {
                sprintf(c, "%.6lf", MAt(A, i, j));
                fprintf(pf, "%s", c);
            }
            else
            {
                sprintf(c, " %.6lf", MAt(A, i, j));
                fprintf(pf, "%s", c);
            }
        }
    }
    fclose(pf);
}
void print_centerline(mtrx u, mtrx v, int nx, int ny, double dx, double dy)
{
    int i, j;
    FILE *f;

    // Ghia et al. (1982), Table 1 — Re=100
    // u-velocity along vertical centerline x=0.5, y in [0,1]
    static const double ghia_y[]  = {0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719,
                                      0.2813, 0.4531, 0.5000, 0.6172, 0.7344, 0.8516,
                                      0.9531, 0.9609, 0.9688, 0.9766, 1.0000};
    static const double ghia_u[]  = {0.0000,-0.0372,-0.0419,-0.0477,-0.0643,-0.1015,
                                     -0.1566,-0.2109,-0.2058,-0.1364, 0.0033, 0.2315,
                                      0.6872, 0.7372, 0.7887, 0.8412, 1.0000};

    // Ghia et al. (1982), Table 2 — Re=100
    // v-velocity along horizontal centerline y=0.5, x in [0,1]
    static const double ghia_x[]  = {0.0000, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563,
                                      0.2266, 0.2344, 0.5000, 0.8047, 0.8594, 0.9063,
                                      0.9453, 0.9531, 0.9609, 0.9688, 1.0000};
    static const double ghia_v[]  = {0.0000, 0.0923, 0.1009, 0.1089, 0.1232, 0.1608,
                                      0.1751, 0.1753, 0.0545,-0.2453,-0.2245,-0.1691,
                                     -0.1031,-0.0886,-0.0739,-0.0591, 0.0000};
    int n_ghia = 17;

    // --- u along vertical centerline (i = nx/2): simulation data ---
    int ci = ny / 2;
    f = fopen("./output/centerline_u_sim.csv", "w");
    if (!f) { printf("Error opening centerline_u_sim.csv\n"); return; }
    fprintf(f, "y,u\n");
    for (i = 0; i < nx; i++)
        fprintf(f, "%.6f,%.6f\n", (i + 0.5) * dy, MAt(u, i, ci));
    fclose(f);

    // --- u along vertical centerline: Ghia et al. (1982) reference ---
    f = fopen("./output/centerline_u_ghia.csv", "w");
    if (!f) { printf("Error opening centerline_u_ghia.csv\n"); return; }
    fprintf(f, "y,u\n");
    for (j = 0; j < n_ghia; j++)
        fprintf(f, "%.6f,%.6f\n", ghia_y[j], ghia_u[j]);
    fclose(f);

    // --- v along horizontal centerline (j = ny/2): simulation data ---
    int cj = nx / 2;
    f = fopen("./output/centerline_v_sim.csv", "w");
    if (!f) { printf("Error opening centerline_v_sim.csv\n"); return; }
    fprintf(f, "x,v\n");
    for (j = 0; j < ny; j++)
        fprintf(f, "%.6f,%.6f\n", (j + 0.5) * dx, MAt(v, cj, j));
    fclose(f);

    // --- v along horizontal centerline: Ghia et al. (1982) reference ---
    f = fopen("./output/centerline_v_ghia.csv", "w");
    if (!f) { printf("Error opening centerline_v_ghia.csv\n"); return; }
    fprintf(f, "x,v\n");
    for (j = 0; j < n_ghia; j++)
        fprintf(f, "%.6f,%.6f\n", ghia_x[j], ghia_v[j]);
    fclose(f);

    printf("Centerline profiles written to output/centerline_u_sim.csv, centerline_u_ghia.csv,\n");
    printf("                             output/centerline_v_sim.csv, centerline_v_ghia.csv\n");
}
