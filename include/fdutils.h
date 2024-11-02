#ifndef _fdutils_h_
#define _fdutils_h_
#include "cstd.h"
#include "waveutils.h"
float laplace(float *p, int n1, int n2, int i1, int i2, float d1, float d2);
void step_backward(acpar par, float *pre, float *curr, float *next, float *vv);
void step_forward(acpar par, float *pre, float *curr, float *next, float *vv);
#endif