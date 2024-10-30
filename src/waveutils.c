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