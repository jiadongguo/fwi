#ifndef _opt_h_
#define _opt_h_
#include "cstd.h"
void opt_init(bool verb_);
float line_search(int n, float *x, float *g, float *d, float fk, opt_fg fg);
void descent(int n, float *g, float *d);
void cg(int n, float *g0, float *g1, float *d0, float *d1, int mode);
#endif
