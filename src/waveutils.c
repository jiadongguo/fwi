#include "waveutils.h"

acpar creat_acpar(const int nz, const int nx,
                  const float dz, const float dx,
                  const int top, const int bot, const int lft, const int rht,
                  const int nt, const float dt,
                  int ns, const int sz, const int sx, const int jsx, const int jsz,
                  int nr, const int rz, const int rx, const int jrx, const int jrz)
{
    acpar ac = alloc1(1, sizeof(*ac));
    /*computation area dimension*/
    ac->nz = nz, ac->nx = nx, ac->nzx = nz * nx;
    ac->dz = dz, ac->dx = dx;
    /*boundary*/
    ac->lft = lft, ac->rht = rht, ac->top = top, ac->bot = bot, ac->nzb = nz + top + bot, ac->nxb = nx + lft + rht, ac->nzxb = ac->nzb * ac->nxb;
    /*recorde time and wavelet*/
    ac->nt = nt, ac->dt = dt;
    /*source*/
    ac->sz = sz, ac->sx = sx, ac->jsx = jsx, ac->jsz = jsz, ac->ns = ns;
    /*receiver*/
    ac->rz = rz, ac->rx = rx, ac->jrx = jrx, ac->jrz = jrz, ac->nr = nr;
    return ac;
}
void record(acpar par, float *p, float *rcd)
{
    int nzb = par->nzb;
    int nr = par->nr;
    int top = par->top, lft = par->lft;
    int jrx = par->jrx, jrz = par->jrz;
    int rx = par->rx + lft, rz = par->rz + top;
    for (int ir = 0; ir < nr; ir++)
    {
        rcd[ir] = p[(rx + ir * jrx) * nzb + (rz + ir * jrz)];
    }
}
void add_src(acpar par, float *p, float wt, int is)
{
    int nzb = par->nzb;
    int top = par->top, lft = par->lft;
    int jsx = par->jsx, jsz = par->jsz;
    int sx = par->sx + lft, sz = par->sz + top;
    sx += is * jsx;
    sz += is * jsz;

    p[sx * nzb + sz] += wt;
}