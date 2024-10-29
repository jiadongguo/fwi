#ifndef _cstd_h_
#define _cstd_h_
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdarg.h>
#include <errno.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include <limits.h>
#include <cblas.h>

#ifndef PI
#define PI 3.1415927
#endif

#ifndef PI2
#define PI2 9.8696044
#endif

#ifndef PI_1
#define PI_1 0.3183099
#endif

#ifndef PI_2
#define PI_2 0.1013212
#endif

#ifndef TOL
#define TOL 1e-12
#endif

#ifndef EPS
#define EPS FLT_EPSILON
#endif

#ifndef NULL
#define NULL ((void *)0)
#endif

#ifndef SGN
#define SGN(x) ((x) < 0 ? -1.0 : 1.0)
#endif

#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif

#ifndef MAX
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#endif

#ifndef NINT
#define NINT(x) ((int)((x) > 0.0 ? (x) + 0.5 : (x) - 0.5))
#endif
/* complex data type */
typedef float cpx[2];
typedef double zpx[2];
/* interface of optimization algorithm*/
typedef float (*opt_fg)(float *x, float *g);
/* allocate and free multi-dimensional arrays */
void *alloc1(size_t n1, size_t size);
void *realloc1(void *v, size_t n1, size_t size);
void **alloc2(size_t n1, size_t n2, size_t size);
void ***alloc3(size_t n1, size_t n2, size_t n3, size_t size);
void ****alloc4(size_t n1, size_t n2, size_t n3, size_t n4, size_t size);
void *****alloc5(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t size);
void ******alloc6(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t n6,
                  size_t size);
void free1(void *p);
void free2(void **p);
void free3(void ***p);
void free4(void ****p);
void free5(void *****p);
void free6(void ******p);

int *alloc1int(size_t n1);
int *realloc1int(int *v, size_t n1);
int **alloc2int(size_t n1, size_t n2);
int ***alloc3int(size_t n1, size_t n2, size_t n3);
int ****alloc4int(size_t n1, size_t n2, size_t n3, size_t n4);
int *****alloc5int(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5);
void free1int(int *p);
void free2int(int **p);
void free3int(int ***p);
void free4int(int ****p);
void free5int(int *****p);

float *alloc1float(size_t n1);
float *realloc1float(float *v, size_t n1);
float **alloc2float(size_t n1, size_t n2);
float ***alloc3float(size_t n1, size_t n2, size_t n3);
float ****alloc4float(size_t n1, size_t n2, size_t n3, size_t n4);
float *****alloc5float(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5);
float ******alloc6float(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t n6);
void free1float(float *p);
void free2float(float **p);
void free3float(float ***p);
void free4float(float ****p);
void free5float(float *****p);
void free6float(float ******p);

double *alloc1double(size_t n1);
double *realloc1double(double *v, size_t n1);
double **alloc2double(size_t n1, size_t n2);
double ***alloc3double(size_t n1, size_t n2, size_t n3);
void free1double(double *p);
void free2double(double **p);
void free3double(double ***p);

cpx *alloc1complexf(size_t n1);
cpx *realloc1complexf(cpx *v, size_t n1);
cpx **alloc2complexf(size_t n1, size_t n2);
cpx ***alloc3complexf(size_t n1, size_t n2, size_t n3);
cpx ****alloc4complexf(size_t n1, size_t n2, size_t n3, size_t n4);
void free1complexf(cpx *p);
void free2complexf(cpx **p);
void free3complexf(cpx ***p);
void free4complexf(cpx ****p);

zpx *alloc1complex(size_t n1);
zpx *realloc1complex(zpx *v, size_t n1);
zpx **alloc2complex(size_t n1, size_t n2);
zpx ***alloc3complex(size_t n1, size_t n2, size_t n3);
zpx ****alloc4complex(size_t n1, size_t n2, size_t n3, size_t n4);
void free1complex(zpx *p);
void free2complex(zpx **p);
void free3complex(zpx ***p);
void free4complex(zpx ****p);

char *alloc1char(size_t n1);
char *realloc1char(char *v, size_t n1);
void free1char(char *p);
/* getpar parameter parsing */
void initargs(int argc, char **argv);
bool getparbool(char *name, bool *p);
int getparint(char *name, int *p);
int getparuint(char *name, unsigned int *p);
int getparshort(char *name, short *p);
int getparushort(char *name, unsigned short *p);
int getparlong(char *name, long *p);
int getparulong(char *name, unsigned long *p);
int getparfloat(char *name, float *p);
int getpardouble(char *name, double *p);
int getparstring(char *name, char **p);
int getparstringarray(char *name, char **p);
int getnparint(int n, char *name, int *p);
int getnparuint(int n, char *name, unsigned int *p);
int getnparshort(int n, char *name, short *p);
int getnparushort(int n, char *name, unsigned short *p);
int getnparlong(int n, char *name, long *p);
int getnparulong(int n, char *name, unsigned long *p);
int getnparfloat(int n, char *name, float *p);
int getnpardouble(int n, char *name, double *p);
int getnparstring(int n, char *name, char **p);
int getnparstringarray(int n, char *name, char **p);
int getnpar(int n, char *name, char *type, void *ptr);
int countparname(char *name);
int countparval(char *name);
int countnparval(int n, char *name);
void checkpars(void);

/* string to numeric conversion with error checking */
short eatoh(char *s);
unsigned short eatou(char *s);
int eatoi(char *s);
unsigned int eatop(char *s);
long eatol(char *s);
unsigned long eatov(char *s);
float eatof(char *s);
double eatod(char *s);

/* string to numeric conversion with error checking */
short eatoh(char *s);
unsigned short eatou(char *s);
int eatoi(char *s);
unsigned int eatop(char *s);
long eatol(char *s);
unsigned long eatov(char *s);
float eatof(char *s);
double eatod(char *s);

void err(char *fmt, ...);
void warn(char *fmt, ...);

/* Claerbout's box-triangle smoothing adapted for 2D */
void triangle_smoothing(float **mod, int n1, int n2, int r1, int r2, int repeat);

/* pading for 2D model */
void pad2(const float *x, float *xx, const int nz, const int nx, const int lft, const int rht, const int top, const int bot);

/* blas function used openblas*/
/* blas2 */
void sscale(int n, float alpha, float *a, float *b);
float sdot(int n, float *a, float *b);
double ddot(int n, double *a, double *b);
void saxpy(int n, float alpha, float *a, float *b);
/* blas3 */
void sgemm(int m, int n, int k, float *A, float *B, float *C);
void dgemm(int m, int n, int k, double *A, double *B, double *C);
#endif
