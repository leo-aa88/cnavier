#include <stdio.h>
#include <stdlib.h>
#include "fluiddyn.h"
#include "linearalg.h"

void euler(mtrx w, mtrx dwdx, mtrx dwdy, mtrx d2wdx2, mtrx d2wdy2, mtrx u, mtrx v, double Re, double dt)
{
    int i, j;

    for (i = 0; i < w.m; i++)
    {
        for (j = 0; j < w.n; j++)
        {
            MAt(w, i, j) = (-MAt(u, i, j) * MAt(dwdx, i, j) - MAt(v, i, j) * MAt(dwdy, i, j) + (1. / Re) * (MAt(d2wdx2, i, j) + MAt(d2wdy2, i, j))) * dt + MAt(w, i, j);
        }
    }
}

mtrx continuity(mtrx dudx, mtrx dvdy)
{
    int i, j;
    mtrx temp;
    temp = initm(dudx.m, dudx.n);

    for (i = 0; i < temp.m; i++)
    {
        for (j = 0; j < temp.n; j++)
        {
            MAt(temp, i, j) = MAt(dudx, i, j) + MAt(dvdy, i, j);
        }
    }
    return temp;
}

mtrx vorticity(mtrx dudy, mtrx dvdx)
{
    int i, j;
    mtrx temp;
    temp = initm(dudy.m, dudy.n);

    for (i = 0; i < temp.m; i++)
    {
        for (j = 0; j < temp.n; j++)
        {
            MAt(temp, i, j) = MAt(dvdx, i, j) - MAt(dudy, i, j);
        }
    }
    return temp;
}
