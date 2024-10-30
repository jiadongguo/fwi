#include "waveutils.h"
#include "cstd.h"
static float *d;    /* absorbing coefficient */
static float alpha; /*theoretical reflection coefficient*/
static int mode;
void eal_init(acpar par, float alpha_, int mode_, float *vv_)
{
    alpha = 1. / alpha_;
    mode = mode_;
    int nz, nx, nzb, nxb, nzxb, lft, top, bot, rht;
    nz = par->nz, nx = par->nx, nzb = par->nzb, nxb = par->nxb, nzxb = par->nzxb, lft = par->lft, top = par->top, bot = par->bot, rht = par->rht;
    float *vv = vv_;
    float refl = log(alpha);
    d = alloc1float(nzxb);
    float dz = par->dz, dx = par->dx, thick, damp;
    /*inner*/
    for (int iz = top; iz < top + nz; iz++)
    {
        for (int ix = lft; ix < nx + lft; ix++)
        {
            d[ix * nzb + iz] = 0;
        }
    }
    if (mode == 0)
    {
        /*top*/
        damp = top * dz;
        for (int iz = 0; iz < top; iz++)
        {
            for (int ix = lft; ix < nx + lft; ix++)
            {
                thick = (top - iz) * dz;
                d[ix * nzb + iz] = refl * 1.5 * vv[ix * nzb + iz] / damp * pow(thick / damp, 2);
            }
        }
        /*bottom*/
        damp = bot * dz;
        for (int iz = top + nz; iz < nzb; iz++)
        {
            for (int ix = lft; ix < nx + lft; ix++)
            {
                thick = (iz + 1 - top - nz) * dz;
                d[ix * nzb + iz] = refl * 1.5 * vv[ix * nzb + iz] / damp * pow(thick / damp, 2);
            }
        }
        /*left*/
        damp = lft * dx;
        for (int iz = top; iz < top + nz; iz++)
        {
            for (int ix = 0; ix < lft; ix++)
            {
                thick = (lft - ix) * dx;
                d[ix * nzb + iz] = refl * 1.5 * vv[ix * nzb + iz] / damp * pow(thick / damp, 2);
            }
        }
        /*right*/
        damp = rht * dx;
        for (int iz = top; iz < top + nz; iz++)
        {
            for (int ix = lft + nx; ix < nxb; ix++)
            {
                thick = (ix + 1 - lft - nx) * dx;
                d[ix * nzb + iz] = refl * 1.5 * vv[ix * nzb + iz] / damp * pow(thick / damp, 2);
            }
        }
        /*left top corner*/
        damp = hypot(top * dz, lft * dx);
        for (int ix = 0; ix < lft; ix++)
        {
            for (int iz = 0; iz < top; iz++)
            {
                thick = hypot((top - iz) * dz, (lft - ix) * dx);
                d[ix * nzb + iz] = refl * 1.5 * vv[ix * nzb + iz] / damp * pow(thick / damp, 2);
            }
        }
        /*right top corner*/
        damp = hypot(top * dz, rht * dx);
        for (int ix = nx + lft; ix < nxb; ix++)
        {
            for (int iz = 0; iz < top; iz++)
            {
                thick = hypot((top - iz) * dz, (ix + 1 - nx - lft) * dx);
                d[ix * nzb + iz] = refl * 1.5 * vv[ix * nzb + iz] / damp * pow(thick / damp, 2);
            }
        }
        /*left bottom corner*/
        damp = hypot(bot * dz, lft * dx);
        for (int ix = 0; ix < lft; ix++)
        {
            for (int iz = top + nz; iz < nzb; iz++)
            {
                thick = hypot((iz + 1 - nz - top) * dz, (lft - ix) * dx);
                d[ix * nzb + iz] = refl * 1.5 * vv[ix * nzb + iz] / damp * pow(thick / damp, 2);
            }
        }
        /*right bottom corner*/
        damp = hypot(bot * dz, rht * dx);
        for (int ix = nx + lft; ix < nxb; ix++)
        {
            for (int iz = top + nz; iz < nzb; iz++)
            {
                thick = hypot((iz + 1 - nz - top) * dz, (ix + 1 - nx - lft) * dx);
                d[ix * nzb + iz] = refl * 1.5 * vv[ix * nzb + iz] / damp * pow(thick / damp, 2);
            }
        }
    }
    else
    {
        float coeff = 2, l2 = log(2);
        /*top*/
        damp = top * dz;
        for (int iz = 0; iz < top; iz++)
        {
            for (int ix = lft; ix < nx + lft; ix++)
            {
                thick = (top - iz) * dz;
                d[ix * nzb + iz] = refl * 1.5 * vv[ix * nzb + iz] / damp * (exp(l2 * pow(thick / damp, coeff)) - 1);
            }
        }
        /*bottom*/
        damp = bot * dz;
        for (int iz = top + nz; iz < nzb; iz++)
        {
            for (int ix = lft; ix < nx + lft; ix++)
            {
                thick = (iz + 1 - top - nz) * dz;
                d[ix * nzb + iz] = refl * 1.5 * vv[ix * nzb + iz] / damp * (exp(l2 * pow(thick / damp, coeff)) - 1);
            }
        }
        /*left*/
        damp = lft * dx;
        for (int iz = top; iz < top + nz; iz++)
        {
            for (int ix = 0; ix < lft; ix++)
            {
                thick = (lft - ix) * dx;
                d[ix * nzb + iz] = refl * 1.5 * vv[ix * nzb + iz] / damp * (exp(l2 * pow(thick / damp, coeff)) - 1);
            }
        }
        /*right*/
        damp = rht * dx;
        for (int iz = top; iz < top + nz; iz++)
        {
            for (int ix = lft + nx; ix < nxb; ix++)
            {
                thick = (ix + 1 - lft - nx) * dx;
                d[ix * nzb + iz] = refl * 1.5 * vv[ix * nzb + iz] / damp * (exp(l2 * pow(thick / damp, coeff)) - 1);
            }
        }
        /*left top corner*/
        damp = hypot(top * dz, lft * dx);
        for (int ix = 0; ix < lft; ix++)
        {
            for (int iz = 0; iz < top; iz++)
            {
                thick = hypot((top - iz) * dz, (lft - ix) * dx);
                d[ix * nzb + iz] = refl * 1.5 * vv[ix * nzb + iz] / damp * (exp(l2 * pow(thick / damp, coeff)) - 1);
            }
        }
        /*right top corner*/
        damp = hypot(top * dz, rht * dx);
        for (int ix = nx + lft; ix < nxb; ix++)
        {
            for (int iz = 0; iz < top; iz++)
            {
                thick = hypot((top - iz) * dz, (ix + 1 - nx - lft) * dx);
                d[ix * nzb + iz] = refl * 1.5 * vv[ix * nzb + iz] / damp * (exp(l2 * pow(thick / damp, coeff)) - 1);
            }
        }
        /*left bottom corner*/
        damp = hypot(bot * dz, lft * dx);
        for (int ix = 0; ix < lft; ix++)
        {
            for (int iz = top + nz; iz < nzb; iz++)
            {
                thick = hypot((iz + 1 - nz - top) * dz, (lft - ix) * dx);
                d[ix * nzb + iz] = refl * 1.5 * vv[ix * nzb + iz] / damp * (exp(l2 * pow(thick / damp, coeff)) - 1);
            }
        }
        /*right bottom corner*/
        damp = hypot(bot * dz, rht * dx);
        for (int ix = nx + lft; ix < nxb; ix++)
        {
            for (int iz = top + nz; iz < nzb; iz++)
            {
                thick = hypot((iz + 1 - nz - top) * dz, (ix + 1 - nx - lft) * dx);
                d[ix * nzb + iz] = refl * 1.5 * vv[ix * nzb + iz] / damp * (exp(l2 * pow(thick / damp, coeff)) - 1);
            }
        }
    }
}

void eal_apply(acpar par, float *pre, float *curr, float *next, float *vv_)
{
    int nz, nx, nzb, nxb, lft, top;
    nz = par->nz, nx = par->nx, nzb = par->nzb, nxb = par->nxb, lft = par->lft, top = par->top;
    float *vv = vv_, dt = par->dt, dz = par->dz, dx = par->dx, tmp, lap;
    for (int ix = 0; ix < nxb; ix++)
    {
        for (int iz = 0; iz < nzb; iz++)
        {
            if (ix >= lft && ix < nx + lft && iz >= top && iz < top + nz)
                continue;
            /* absorbing area */
            tmp = 1. / (d[ix * nzb + iz] * dt + 1);
            lap = laplace(nzb, nxb, iz, ix, curr, dz, dx);
            next[ix * nzb + iz] = (d[ix * nzb + iz] * dt - 1) * tmp * pre[ix * nzb + iz] + (2 - pow(d[ix * nzb + iz] * dt, 2)) * tmp * curr[ix * nzb + iz] + pow(vv[ix * nzb + iz] * dt, 2) * tmp * lap;
        }
    }
}
void eal_close()
{
    free1(d);
}