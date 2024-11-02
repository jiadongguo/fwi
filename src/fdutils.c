#include "fdutils.h"
/* Laplace filter */
void eal_apply(acpar par, float *pre, float *curr, float *next, float *vv_);

float laplace(float *p, int n1, int n2, int i1, int i2, float d1, float d2)
{
    static float c0 = -5. / 2, c1 = 4. / 3, c2 = -1. / 12;
    float fdx[5], fdz[5];
    float lapz, lapx;
    d1 = 1. / pow(d1, 2);
    d2 = 1. / pow(d2, 2);
    for (int i = 0; i < 5; i++)
    {
        fdx[i] = ((i2 + i - 2 < 0) || (i2 + i - 2 >= n2)) ? 0 : p[(i2 + i - 2) * n1 + i1];
        fdz[i] = ((i1 + i - 2 < 0) || (i1 + i - 2 >= n1)) ? 0 : p[i2 * n1 + i1 + i - 2];
    }
    lapz = (c2 * (fdz[0] + fdz[4]) + c1 * (fdz[1] + fdz[3]) + c0 * fdz[2]) * d1;
    lapx = (c2 * (fdx[0] + fdx[4]) + c1 * (fdx[1] + fdx[3]) + c0 * fdx[2]) * d2;
    return lapx + lapz;
}
void step_forward(acpar par, float *pre, float *curr, float *next, float *vv)
{
    float dz = par->dz, dx = par->dx, dt = par->dt;
    int nz = par->nz, nx = par->nx, nzb = par->nzb, nxb = par->nxb, top = par->top, lft = par->lft;
    /*only for inner grid*/
    for (int ix = lft; ix < lft + nx; ix++)
    {
        for (int iz = top; iz < top + nz; iz++)
        {

            next[ix * nzb + iz] = laplace(curr, nzb, nxb, iz, ix, dz, dx) * (vv[ix * nzb + iz] * dt) * (vv[ix * nzb + iz] * dt) + 2 * curr[ix * nzb + iz] - pre[ix * nzb + iz];
        }
    }
    eal_apply(par, pre, curr, next, vv);
}
void step_backward(acpar par, float *pre, float *curr, float *next, float *vv)
{
    float dz = par->dz, dx = par->dx, dt = par->dt;
    int nz = par->nz, nx = par->nx, nzb = par->nzb, nxb = par->nxb, top = par->top, lft = par->lft;
    /*only for inner grid*/
    for (int ix = lft; ix < lft + nx; ix++)
    {
        for (int iz = top; iz < top + nz; iz++)
        {

            next[ix * nzb + iz] = laplace(curr, nzb, nxb, iz, ix, dz, dx) * (vv[ix * nzb + iz] * dt) * (vv[ix * nzb + iz] * dt) + 2 * curr[ix * nzb + iz] - pre[ix * nzb + iz];
        }
    }
    eal_apply(par, pre, curr, next, vv);
}