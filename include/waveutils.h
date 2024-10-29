#ifndef _waveutils_h_
#define _waveutils_h_
#include "cstd.h"
typedef struct acpar
{
    /*computation area dimension*/
    int nz, nx;
    float dx, dz;
    int nzx;
    /*boundary*/
    int lft, rht, top, bot;
    int nzb, nxb, nzxb;
    /*recorde time and wavelet*/
    int nt;
    float dt;
    /*velocity model*/
    float *v, *vv;
    /*source*/
    int sz, sx /*start form 0*/, jsx, jsz, ns;
    /*receiver*/
    int rz, rx /*start form 0*/, jrx, jrz, nr;
} *acpar;
acpar creat_acpar(const int nz, const int nx,
                  const float dz, const float dx,
                  const int top, const int bot, const int lft, const int rht,
                  const int nt, const float dt,
                  int ns, const int sz, const int sx, const int jsx, const int jsz,
                  int nr, const int rz, const int rx, const int jrx, const int jrz,
                  float *vel);
#endif