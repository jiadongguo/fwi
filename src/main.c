/*
seismic wavefield modeling with effective absorbing layer
*/
#include "cstd.h"
#include "waveutils.h"
#include <mpi/mpi.h>
float laplace(int n1, int n2, int i1, int i2, float *curr, float d1, float d2);
void eal_init(acpar par, float alpha_, int mode);
void eal_apply(acpar par, float *pre, float *curr, float *next);
void eal_close();
void fdfor(acpar par, float *pre, float *curr, float *next)
{
    float *vv;
    int nz, nx, nzb, nxb, top, lft;
    float dz, dx, dt;
    nz = par->nz;
    nx = par->nx;
    nzb = par->nzb;
    nxb = par->nxb;
    top = par->top;
    lft = par->lft;
    dz = par->dz;
    dx = par->dx;
    dt = par->dt;
    vv = par->vv;
    /*only for inner grid*/
    for (int ix = lft; ix < lft + nx; ix++)
    {
        for (int iz = top; iz < top + nz; iz++)
        {
            next[ix * nzb + iz] = laplace(nzb, nxb, iz, ix, curr, dz, dx) * (vv[ix * nzb + iz] * dt) * (vv[ix * nzb + iz] * dt) + 2 * curr[ix * nzb + iz] - pre[ix * nzb + iz];
        }
    }
    eal_apply(par, pre, curr, next);
}
int main(int argc, char **argv)
{
    initargs(argc, argv);
    char *shots; /*observed seismogram*/
    if (!getparstring("shots", &shots))
        err("need shots for observed seismogram");
    FILE *fshots = fopen(shots, "rb");
    
    MPI_Finalize();
    return 0;
}