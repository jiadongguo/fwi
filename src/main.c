/*
a simple program for fwi
*/

#include "cstd.h"
#include "waveutils.h"

/* ======================================================================================

laplace : laplace opeprator
forward : Seismic waves travel forward
backward : Seismic waves travel backward
obj : objective function of fwi
record :
gradient : gradient of current model

===========================================================================================*/

void eal_init(acpar par, float alpha_, int mode);
void eal_apply(acpar par, float *pre, float *curr, float *next);
void eal_close();
float laplace(int n1, int n2, int i1, int i2, float *curr, float d1, float d2)
{
    static float c0 = -5. / 2, c1 = 4. / 3, c2 = -1. / 12;
    float fdx[5], fdz[5];
    float lapz, lapx;
    d1 = 1. / pow(d1, 2);
    d2 = 1. / pow(d2, 2);
    for (int i = 0; i < 5; i++)
    {
        fdx[i] = ((i2 + i - 2 < 0) || (i2 + i - 2 >= n2)) ? 0 : curr[(i2 + i - 2) * n1 + i1];
        fdz[i] = ((i1 + i - 2 < 0) || (i1 + i - 2 >= n1)) ? 0 : curr[i2 * n1 + i1 + i - 2];
    }
    lapz = (c2 * (fdz[0] + fdz[4]) + c1 * (fdz[1] + fdz[3]) + c0 * fdz[2]) * d1;
    lapx = (c2 * (fdx[0] + fdx[4]) + c1 * (fdx[1] + fdx[3]) + c0 * fdx[2]) * d2;
    return lapx + lapz;
}

void fun_forward(acpar par, float *pre, float *curr, float *next)
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
void fun_backward(acpar par, float *pre, float *curr, float *next)
{
    float *vv;
    int nz, nx, nzb, nxb, top, lft;
    int nr = par->nr;
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
float fun_obj(acpar par, float *wt)
{
    float obj = 0;
}

void fun_gradient(float *g /* backward wavefield */, float *lap, int it, acpar par, float *grad)
{
    int nx = par->nx, nz = par->nz;
    int nzb = par->nzb, lft = par->lft, top = par->top;
    float *v = par->v;
    for (int ix = 0; ix < nx; ix++)
    {
        for (int iz = 0; iz < nz; iz++)
        {
            grad[ix * nz + iz] += g[(ix + lft) * nzb + iz + top] * lap[it * nz * nx + ix * nz + iz] / pow(v[ix * nz + iz], 2);
        }
    }
}
void fun_record(acpar par, float *curr, float *rcd)
{
    int nr = par->nr, nzb = par->nzb;
    int rx, rz;
    for (int ir = 0; ir < nr; ir++)
    {
        rx = par->rx + par->lft + par->jsx * ir;
        rz = par->rz + par->top + par->jsz * ir;
        rcd[ir] = curr[rx * nzb + rz];
    }
}
float fun_alphatest(int n, acpar par, float *grad)
{
    float *v = par->v;
    float maxv = 0, maxg = 0;
    for (int i = 0; i < n; i++)
    {
        maxv = MAX(maxv, fabs(v[i]));
        maxg = MAX(maxg, fabs(grad[i]));
    }
    return 0.01 * maxv / maxg;
}
float fun_update(acpar par, float alpha, float *grad)
{
    int nzx = par->nzx;
    int nz = par->nz, nx = par->nx, lft = par->lft, rht = par->rht, top = par->top, bot = par->bot;
    float *v = par->v, *vv = par->vv;
    for (int i = 0; i < nzx; i++)
    {
        v[i] = v[i] - alpha * grad[i];
    }
    pad2(v, vv, nz, nx, lft, rht, top, bot);
}
int main(int argc, char **argv)
{
    initargs(argc, argv);
    int nter;                       /* total number of iterations */
    char *fwt, *fvel, *shots, *out; /*observed seismogram*/
    int nz, nx, nt, top, bot, lft, rht;
    int ns, sz, sx, jsx, jsz, rz, rx, jrx, jrz, nr;
    float dz, dx, dt;
    float *wt /* wavelet */, *vel /* initial velocity model */;
    float *vtmp, *vtmppad;
    bool verb;
    if (!getparbool("verb", &verb))
        verb = false;
    if (!getparint("nter", &nter))
    {
        nter = 10;
    }
    if (!getparint("n1", &nz))
        err("need nz");
    if (!getparint("n2", &nx))
        err("need nx");
    if (!getparint("nt", &nt))
        err("need nt");
    if (!getparint("top", &top))
        err("need top");
    if (!getparint("bot", &bot))
        err("need bot");
    if (!getparint("lft", &lft))
        err("need lft");
    if (!getparint("rht", &rht))
        err("need rht");
    if (!getparint("ns", &ns))
        err("need ns");
    if (!getparint("nr", &nr))
        err("need nr");
    if (!getparint("sz", &sz))
        err("need sz");
    if (!getparint("sx", &sx))
        err("need sx");
    if (!getparint("jsx", &jsx))
        err("need jsx");
    if (!getparint("jsz", &jsz))
        err("need jsz");
    if (!getparint("rz", &rz))
        err("need rz");
    if (!getparint("rx", &rx))
        err("need rx");
    if (!getparint("jrx", &jrx))
        err("need jrx");
    if (!getparint("jrz", &jrz))
        err("need jrz");
    if (!getparfloat("dt", &dt))
        err("need dt");
    if (!getparfloat("d2", &dx))
        err("need dx");
    if (!getparfloat("d1", &dz))
        err("need dz");
    if (!getparstring("shots", &shots))
        err("need shots for observed seismogram");
    if (!getparstring("vpfile", &fvel))
        err("need vp");
    if (!getparstring("wtfile", &fwt))
        err("need fwt");
    if (!getparstring("out", &out))
        err("need out");
    if ((sx + (ns - 1) * jsx >= nx) || (sz + (ns - 1) * jsz >= nz))
        err("source position outer bound");
    if ((rx + (nr - 1) * jrx >= nx) || (rz + (nr - 1) * jrz >= nz))
        err("receiver position outer bound");
    /* ===================================================================================================== */
    acpar par = creat_acpar(nz, nx, dz, dx, top, bot, lft, rht, nt, dt, ns, sz, sx, jsx, jsz, nr, rz, rx, jrx, jrz, vel);

    FILE *fshots = fopen(shots, "rb");
    FILE *fout = fopen(out, "wb");
    int nzxb = par->nzxb, nzb = par->nzb, nxb = par->nxb, nzx = par->nzx;
    float *p0 = alloc1float(par->nzxb);
    float *p1 = alloc1float(par->nzxb);
    float *p2 = alloc1float(par->nzxb);
    float *grad = alloc1float(par->nzx);
    float *tmp, *dobs, *dtobs, *lap, *derr, *dcal;
    wt = alloc1float(nt);
    dtobs = alloc1float(nt * nr);
    dobs = alloc1float(nt * nr);
    dcal = alloc1float(nt * nr);
    derr = alloc1float(nt * nr);
    lap = alloc1float(nzx);
    vtmp = alloc1float(nzx);
    memcpy(vtmp, vel, sizeof(float) * nzx);
    acpar partmp = creat_acpar(nz, nx, dz, dx, top, bot, lft, rht, nt, dt, ns, sz, sx, jsx, jsz, nr, rz, rx, jrx, jrz, vtmp); /* test model */
    eal_init(par, 1e-3, 0);
    for (int iter = 0; iter < nter; iter++)
    {
        float obj = 0;
        memset(grad, 0, sizeof(float) * nzxb);
        /*read observed seismogram*/
        fseek(fshots, 0, SEEK_SET);
        /* forward modeling */
        for (int is = 0; is < ns; is++)
        {
            fread(dtobs, sizeof(float) * nr * nt, 1, fshots);
            transp(nt, nr, dtobs, dobs);
            memset(p0, 0, sizeof(float) * nzxb);
            memset(p1, 0, sizeof(float) * nzxb);
            memset(p2, 0, sizeof(float) * nzxb);
            int sx = lft + jsx * is, sz = top + jsz * is;
            sx = sx * nzb + sz;
            for (int it = 0; it < nt; it++)
            {
                fun_forward(par, p0, p1, p2);
                tmp = p0, p0 = p1, p1 = p2, p2 = tmp;
                p1[sx] += wt[it];
                for (int ix = 0; ix < nx; ix++)
                {
                    for (int iz = 0; iz < nz; iz++)
                    {
                        lap[it * nz * nx + ix * nz + iz] = laplace(nzb, nxb, iz + top, ix + lft, p1, dz, dx);
                    }
                }
                fun_record(par, p1, dcal + it * nr);
                /* calculate wavefield residuals and  objective value */
                obj += fun_obj(dcal, dobs, derr + it * nr, nr);
            }
            /* calculate gradient */
            memset(p0, 0, sizeof(float) * nzxb);
            memset(p1, 0, sizeof(float) * nzxb);
            memset(p2, 0, sizeof(float) * nzxb);
            for (int it = nt - 1; it >= 0; it--)
            {
                fun_backward(par, p0, p1, p2);
                tmp = p0, p0 = p1, p1 = p2, p2 = tmp;
                for (int ir = 0; ir < nr; ir++)
                {
                    p1[(rx + lft + jrx * ir) * nzb + top + jrz * ir + rz] += derr[ir * nt + it];
                }
                fun_gradient(p1, lap, it, par, grad);
            }
        }
        /* update velocity model */
        /* test model */
        float alphatest = fun_alphatest(nzx, par, grad);
        memcpy(vtmp, vel, sizeof(float) * nzx);
        fun_update(partmp, alphatest, grad);
        pad2(vtmp, partmp->vv, nz, nx, lft, rht, top, bot);
    }
    eal_close();
    return 0;
}