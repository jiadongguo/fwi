#include "cstd.h"

/* allocate a 1-d array */
void *alloc1(size_t n1, size_t size)
{
    void *p;

    if ((p = malloc(n1 * size)) == NULL)
        return NULL;
    return p;
}

/* re-allocate a 1-d array */
void *realloc1(void *v, size_t n1, size_t size)
{
    void *p;

    if ((p = realloc(v, n1 * size)) == NULL)
        return NULL;
    return p;
}

/* free a 1-d array */
void free1(void *p)
{
    free(p);
}

/* allocate a 2-d array */
void **alloc2(size_t n1, size_t n2, size_t size)
{
    size_t i2;
    void **p;

    if ((p = (void **)malloc(n2 * sizeof(void *))) == NULL)
        return NULL;
    if ((p[0] = (void *)malloc(n2 * n1 * size)) == NULL)
    {
        free(p);
        return NULL;
    }
    for (i2 = 0; i2 < n2; i2++)
        p[i2] = (char *)p[0] + size * n1 * i2;
    return p;
}

/* free a 2-d array */
void free2(void **p)
{
    free(p[0]);
    free(p);
}

/* allocate a 3-d array */
void ***alloc3(size_t n1, size_t n2, size_t n3, size_t size)
{
    size_t i3, i2;
    void ***p;

    if ((p = (void ***)malloc(n3 * sizeof(void **))) == NULL)
        return NULL;
    if ((p[0] = (void **)malloc(n3 * n2 * sizeof(void *))) == NULL)
    {
        free(p);
        return NULL;
    }
    if ((p[0][0] = (void *)malloc(n3 * n2 * n1 * size)) == NULL)
    {
        free(p[0]);
        free(p);
        return NULL;
    }

    for (i3 = 0; i3 < n3; i3++)
    {
        p[i3] = p[0] + n2 * i3;
        for (i2 = 0; i2 < n2; i2++)
            p[i3][i2] = (char *)p[0][0] + size * n1 * (i2 + n2 * i3);
    }
    return p;
}

/* free a 3-d array */
void free3(void ***p)
{
    free(p[0][0]);
    free(p[0]);
    free(p);
}

/* allocate a 4-d array */
void ****alloc4(size_t n1, size_t n2, size_t n3, size_t n4, size_t size)
{
    size_t i4, i3, i2;
    void ****p;

    if ((p = (void ****)malloc(n4 * sizeof(void ***))) == NULL)
        return NULL;
    if ((p[0] = (void ***)malloc(n4 * n3 * sizeof(void **))) == NULL)
    {
        free(p);
        return NULL;
    }
    if ((p[0][0] = (void **)malloc(n4 * n3 * n2 * sizeof(void *))) == NULL)
    {
        free(p[0]);
        free(p);
        return NULL;
    }
    if ((p[0][0][0] = (void *)malloc(n4 * n3 * n2 * n1 * size)) == NULL)
    {
        free(p[0][0]);
        free(p[0]);
        free(p);
        return NULL;
    }
    for (i4 = 0; i4 < n4; i4++)
    {
        p[i4] = p[0] + i4 * n3;
        for (i3 = 0; i3 < n3; i3++)
        {
            p[i4][i3] = p[0][0] + n2 * (i3 + n3 * i4);
            for (i2 = 0; i2 < n2; i2++)
                p[i4][i3][i2] = (char *)p[0][0][0] +
                                size * n1 * (i2 + n2 * (i3 + n3 * i4));
        }
    }
    return p;
}

/* free a 4-d array */
void free4(void ****p)
{
    free(p[0][0][0]);
    free(p[0][0]);
    free(p[0]);
    free(p);
}

/* The following two functions were added by Zhaobo Meng, Jan. 1997*/
/* allocate a 5-d array */
void *****alloc5(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t size)
{
    size_t i5, i4, i3, i2;
    void *****p;

    if ((p = (void *****)malloc(n5 * sizeof(void ****))) == NULL)
        return NULL;
    if ((p[0] = (void ****)malloc(n5 * n4 * sizeof(void ***))) == NULL)
    {
        free(p);
        return NULL;
    }
    if ((p[0][0] = (void ***)malloc(n5 * n4 * n3 * sizeof(void **))) == NULL)
    {
        free(p[0]);
        free(p);
        return NULL;
    }
    if ((p[0][0][0] = (void **)malloc(n5 * n4 * n3 * n2 * sizeof(void *))) == NULL)
    {
        free(p[0][0]);
        free(p[0]);
        free(p);
        return NULL;
    }
    if ((p[0][0][0][0] = (void *)malloc(n5 * n4 * n3 * n2 * n1 * size)) == NULL)
    {
        free(p[0][0][0]);
        free(p[0][0]);
        free(p[0]);
        free(p);
        return NULL;
    }
    for (i5 = 0; i5 < n5; i5++)
    {
        p[i5] = p[0] + n4 * i5;
        for (i4 = 0; i4 < n4; i4++)
        {
            p[i5][i4] = p[0][0] + n3 * (i4 + n4 * i5);
            for (i3 = 0; i3 < n3; i3++)
            {
                p[i5][i4][i3] = p[0][0][0] + n2 * (i3 + n3 * (i4 + n4 * i5));
                for (i2 = 0; i2 < n2; i2++)
                    p[i5][i4][i3][i2] = (char *)p[0][0][0][0] +
                                        size * n1 * (i2 + n2 * (i3 + n3 * (i4 + n4 * i5)));
            }
        }
    }
    return p;
}

/* free a 5-d array */
void free5(void *****p)
{
    free(p[0][0][0][0]);
    free(p[0][0][0]);
    free(p[0][0]);
    free(p[0]);
    free(p);
}

/* The following two functions were added by Zhaobo Meng, Jan. 1997*/
/* allocate a 6-d array */
void ******alloc6(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t n6,
                  size_t size)
{
    size_t i6, i5, i4, i3, i2;
    void ******p;

    if ((p = (void ******)malloc(n6 * sizeof(void *****))) == NULL)
        return NULL;

    if ((p[0] = (void *****)malloc(n6 * n5 * sizeof(void ****))) == NULL)
    {
        free(p);
        return NULL;
    }

    if ((p[0][0] = (void ****)malloc(n6 * n5 * n4 * sizeof(void ***))) == NULL)
    {
        free(p[0]);
        free(p);
        return NULL;
    }
    if ((p[0][0][0] = (void ***)malloc(n6 * n5 * n4 * n3 * sizeof(void **))) == NULL)
    {
        free(p[0][0]);
        free(p[0]);
        free(p);
        return NULL;
    }
    if ((p[0][0][0][0] = (void **)malloc(n6 * n5 * n4 * n3 * n2 * sizeof(void *))) == NULL)
    {
        free(p[0][0][0]);
        free(p[0][0]);
        free(p[0]);
        free(p);
        return NULL;
    }
    if ((p[0][0][0][0][0] = (void *)malloc(n6 * n5 * n4 * n3 * n2 * n1 * size)) == NULL)
    {
        free(p[0][0][0][0]);
        free(p[0][0][0]);
        free(p[0][0]);
        free(p[0]);
        free(p);
        return NULL;
    }

    for (i6 = 0; i6 < n6; i6++)
    {
        p[i6] = p[0] + n5 * i6;
        for (i5 = 0; i5 < n5; i5++)
        {
            p[i6][i5] = p[0][0] + n4 * (i5 + n5 * i6);
            for (i4 = 0; i4 < n4; i4++)
            {
                p[i6][i5][i4] = p[0][0][0] + n3 * (i4 + n4 * (i5 + n5 * i6));
                for (i3 = 0; i3 < n3; i3++)
                {
                    p[i6][i5][i4][i3] = p[0][0][0][0] + n2 * (i3 + n3 * (i4 + n4 * (i5 + n5 * i6)));
                    for (i2 = 0; i2 < n2; i2++)
                        p[i6][i5][i4][i3][i2] =
                            (char *)p[0][0][0][0][0] +
                            size * n1 * (i2 + n2 * (i3 + n3 * (i4 + n4 * (i5 + n5 * i6))));
                }
            }
        }
    }
    return p;
}

/* free a 6-d array */
void free6(void ******p)
{
    free(p[0][0][0][0][0]);
    free(p[0][0][0][0]);
    free(p[0][0][0]);
    free(p[0][0]);
    free(p[0]);
    free(p);
}

/* allocate a 1-d array of ints */
int *alloc1int(size_t n1)
{
    return (int *)alloc1(n1, sizeof(int));
}

/* re-allocate a 1-d array of ints */
int *realloc1int(int *v, size_t n1)
{
    return (int *)realloc1(v, n1, sizeof(int));
}

/* free a 1-d array of ints */
void free1int(int *p)
{
    free1(p);
}

/* allocate a 2-d array of ints */
int **alloc2int(size_t n1, size_t n2)
{
    return (int **)alloc2(n1, n2, sizeof(int));
}

/* free a 2-d array of ints */
void free2int(int **p)
{
    free2((void **)p);
}

/* allocate a 3-d array of ints */
int ***alloc3int(size_t n1, size_t n2, size_t n3)
{
    return (int ***)alloc3(n1, n2, n3, sizeof(int));
}

/* free a 3-d array of ints */
void free3int(int ***p)
{
    free3((void ***)p);
}

/* allocate a 4-d array of ints, added by Zhaobo Meng, 1997 */
int ****alloc4int(size_t n1, size_t n2, size_t n3, size_t n4)
{
    return (int ****)alloc4(n1, n2, n3, n4, sizeof(int));
}

/* free a 4-d array of ints, added by Zhaobo Meng, 1997 */
void free4int(int ****p)
{
    free4((void ****)p);
}

/* allocate a 5-d array of ints, added by Zhaobo Meng, 1997 */
int *****alloc5int(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5)
{
    return (int *****)alloc5(n1, n2, n3, n4, n5, sizeof(int));
}

/* free a 5-d array of ints, added by Zhaobo Meng, 1997 */
void free5int(int *****p)
{
    free5((void *****)p);
}

/* allocate a 1-d array of floats */
float *alloc1float(size_t n1)
{
    return (float *)alloc1(n1, sizeof(float));
}

/* re-allocate a 1-d array of floats */
float *realloc1float(float *v, size_t n1)
{
    return (float *)realloc1(v, n1, sizeof(float));
}

/* free a 1-d array of floats */
void free1float(float *p)
{
    free1(p);
}

/* allocate a 2-d array of floats */
float **alloc2float(size_t n1, size_t n2)
{
    return (float **)alloc2(n1, n2, sizeof(float));
}

/* free a 2-d array of floats */
void free2float(float **p)
{
    free2((void **)p);
}

/* allocate a 3-d array of floats */
float ***alloc3float(size_t n1, size_t n2, size_t n3)
{
    return (float ***)alloc3(n1, n2, n3, sizeof(float));
}

/* free a 3-d array of floats */
void free3float(float ***p)
{
    free3((void ***)p);
}

/* allocate a 4-d array of floats, added by Zhaobo Meng, 1997 */
float ****alloc4float(size_t n1, size_t n2, size_t n3, size_t n4)
{
    return (float ****)alloc4(n1, n2, n3, n4, sizeof(float));
}

/* free a 4-d array of floats, added by Zhaobo Meng, 1997 */
void free4float(float ****p)
{
    free4((void ****)p);
}

/* allocate a 5-d array of floats, added by Zhaobo Meng, 1997 */
float *****alloc5float(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5)
{
    return (float *****)alloc5(n1, n2, n3, n4, n5, sizeof(float));
}

/* free a 5-d array of floats, added by Zhaobo Meng, 1997 */
void free5float(float *****p)
{
    free5((void *****)p);
}

/* allocate a 6-d array of floats, added by Zhaobo Meng, 1997 */
float ******alloc6float(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t n6)
{
    return (float ******)alloc6(n1, n2, n3, n4, n5, n6, sizeof(float));
}

/* free a 6-d array of floats, added by Zhaobo Meng, 1997 */
void free6float(float ******p)
{
    free6((void ******)p);
}

/* allocate a 1-d array of doubles */
double *alloc1double(size_t n1)
{
    return (double *)alloc1(n1, sizeof(double));
}

/* re-allocate a 1-d array of doubles */
double *realloc1double(double *v, size_t n1)
{
    return (double *)realloc1(v, n1, sizeof(double));
}

/* free a 1-d array of doubles */
void free1double(double *p)
{
    free1(p);
}

/* allocate a 2-d array of doubles */
double **alloc2double(size_t n1, size_t n2)
{
    return (double **)alloc2(n1, n2, sizeof(double));
}

/* free a 2-d array of doubles */
void free2double(double **p)
{
    free2((void **)p);
}

/* allocate a 3-d array of doubles */
double ***alloc3double(size_t n1, size_t n2, size_t n3)
{
    return (double ***)alloc3(n1, n2, n3, sizeof(double));
}

/* free a 3-d array of doubles */
void free3double(double ***p)
{
    free3((void ***)p);
}

/* allocate a 1-d array of complexs */
cpx *alloc1complexf(size_t n1)
{
    return (cpx *)alloc1(n1, sizeof(cpx));
}

/* re-allocate a 1-d array of complexs */
cpx *realloc1complexf(cpx *v, size_t n1)
{
    return (cpx *)realloc1(v, n1, sizeof(cpx));
}

/* free a 1-d array of complexs */
void free1complexf(cpx *p)
{
    free1(p);
}

/* allocate a 2-d array of complexs */
cpx **alloc2complexf(size_t n1, size_t n2)
{
    return (cpx **)alloc2(n1, n2, sizeof(cpx));
}

/* free a 2-d array of complexs */
void free2complexf(cpx **p)
{
    free2((void **)p);
}

/* allocate a 3-d array of complexs */
cpx ***alloc3complexf(size_t n1, size_t n2, size_t n3)
{
    return (cpx ***)alloc3(n1, n2, n3, sizeof(cpx));
}

/* free a 3-d array of complexs */
void free3complexf(cpx ***p)
{
    free3((void ***)p);
}

/* allocate a 4-d array of complexs */
cpx ****alloc4complexf(size_t n1, size_t n2, size_t n3, size_t n4)
{
    return (cpx ****)alloc4(n1, n2, n3, n4, sizeof(cpx));
}

/* free a 4-d array of complexs */
void free4complexf(cpx ****p)
{
    free4((void ****)p);
}

/* allocate a 1-d array of complexs */
zpx *alloc1complex(size_t n1)
{
    return (zpx *)alloc1(n1, sizeof(zpx));
}

/* re-allocate a 1-d array of complexs */
zpx *realloc1complex(zpx *v, size_t n1)
{
    return (zpx *)realloc1(v, n1, sizeof(zpx));
}

/* free a 1-d array of complexs */
void free1complex(zpx *p)
{
    free1(p);
}

/* allocate a 2-d array of complexs */
zpx **alloc2complex(size_t n1, size_t n2)
{
    return (zpx **)alloc2(n1, n2, sizeof(zpx));
}

/* free a 2-d array of complexs */
void free2complex(zpx **p)
{
    free2((void **)p);
}

/* allocate a 3-d array of complexs */
zpx ***alloc3complex(size_t n1, size_t n2, size_t n3)
{
    return (zpx ***)alloc3(n1, n2, n3, sizeof(zpx));
}

/* free a 3-d array of complexs */
void free3complex(zpx ***p)
{
    free3((void ***)p);
}

/* allocate a 4-d array of complexs */
zpx ****alloc4complex(size_t n1, size_t n2, size_t n3, size_t n4)
{
    return (zpx ****)alloc4(n1, n2, n3, n4, sizeof(zpx));
}

/* free a 4-d array of complexs */
void free4complex(zpx ****p)
{
    free4((void ****)p);
}

//============================================================
/* allocate a 1-d array of chars */
char *alloc1char(size_t n1)
{
    return (char *)alloc1(n1, sizeof(char));
}

/* re-allocate a 1-d array of chars */
char *realloc1char(char *v, size_t n1)
{
    return (char *)realloc1(v, n1, sizeof(char));
}

/* free a 1-d array of chars */
void free1char(char *p)
{
    free1(p);
}

/*
initargs 	        Makes command line args available to subroutines (re-entrant).
                    Every par program starts with this call!
getparbool		    get bool
getparint		    get integers
getparuint		    get unsigned integers
getparshort		    get short integers
getparushort		get unsigned short integers
getparlong		    get long integers
getparulong		    get unsigned long integers
getparfloat		    get float
getpardouble		get double
getparstring		get a single string
getparstringarray	get string array (fields delimited by commas)
getpar			    get parameter by type
getnparint		    get n'th occurrence of integer
getnparuint		    get n'th occurrence of unsigned int
getnparshort		get n'th occurrence of short integer
getnparushort		get n'th occurrence of unsigned short int
getnparlong		    get n'th occurrence of long integer
getnparulong		get n'th occurrence of unsigned long int
getnparfloat		get n'th occurrence of float
getnpardouble		get n'th occurrence of double
getnparstring		get n'th occurrence of string
getnparstringarray	get n'th occurrence of string array
getnpar			    get n'th occurrence by type
countparname		return the number of times a parameter names is used
countparval		    return the number of values in the last occurrence
                    of a parameter
countnparval		return the number of values in the n'th occurrence
                    of a parameter
getPar			    Promax compatible version of getpar
checkpars()		    check the argument list for typos
******************************************************************************
Notes:
Here are some usage examples:

    ... if integer n not specified, then default to zero.
    if (!getparint("n", &n)) n = 0;

    ... if array of floats vx is specified, then
    if (nx=countparval("vx")) {
        ... allocate space for array
        vx = alloc1float(nx);
        ... and get the floats
        getparfloat("vx",vx);
    }

The command line for the above examples might look like:
    progname n=35 vx=3.21,4,9.5
    Every par program starts with this call!

More examples are provided in the DTEST code at the end of this file.

The functions: eatoh, eatou, eatol, eatov, eatoi, eatop used
below are versions of atoi that check for overflow.  The source
file for these functions is atopkge.c.

******************************************************************************
*/
/* parameter table */
typedef struct
{
    char *name;     /* external name of parameter	*/
    char *asciival; /* ascii value of parameter	*/
} pointer_table;
/* global variables declared and used internally */
static pointer_table *argtbl; /* parameter table		*/
static int nargs;             /* number of args that parse	*/
static int tabled = false;    /* true when parameters tabled 	*/
static size_t targc;          /* total number of args		*/
static char **targv;          /* pointer to arg strings	*/
static char *argstr;          /* storage for command line	*/
/* functions declared and used internally */
static int getparindex(int n, char *name);
static void getparinit(void);
static void tabulate(size_t argc, char **argv);
static char *getpfname(void);
static int ccount(char c, char *s);
/*--------------------------------------------------------------------*\
  These variables are used by checkpars() to warn of parameter typos.
  par= does not use getpar() so we need to store that as the first
  parameter name.  lheaders= is buried in fgettr() so we intialize
  that also
  \*--------------------------------------------------------------------*/
#define PAR_NAMES_MAX 512

static char *par_names[PAR_NAMES_MAX];
static int par_count = 0;
static int parcheck = 0;

int xargc;
char **xargv;
/* make command line args available to subroutines -- re-entrant version */
void initargs(int argc, char **argv)
{
    memset(par_names, 0, sizeof(par_names));
    par_names[0] = "par";
    par_names[1] = "lheader";
    par_count = 2;

    xargc = argc;
    xargv = argv;
    if (tabled == true)
    {
        free(argstr);
        free(targv);
        free(argtbl);
    }
    tabled = false;
    return;
}
void strchop(char *s, char *t)
/***********************************************************************
strchop - chop off the tail end of a string "s" after a "," returning
      the front part of "s" as "t".
      ************************************************************************
Notes:
Based on strcpy in Kernighan and Ritchie's C [ANSI C] book, p. 106.
************************************************************************
Author: CWP: John Stockwell and Jack K. Cohen, July 1995
***********************************************************************/
{

    while ((*s != ',') && (*s != '\0'))
    {
        *t++ = *s++;
    }
    *t = '\0';
}

/* functions to get values for the last occurrence of a parameter name */
bool getparbool(char *name, bool *ptr)
{
    return getnpar(0, name, "b", ptr);
}
int getparint(char *name, int *ptr)
{
    return getnpar(0, name, "i", ptr);
}
int getparuint(char *name, unsigned int *ptr)
{
    return getnpar(0, name, "p", ptr);
}
int getparshort(char *name, short *ptr)
{
    return getnpar(0, name, "h", ptr);
}
int getparushort(char *name, unsigned short *ptr)
{
    return getnpar(0, name, "u", ptr);
}
int getparlong(char *name, long *ptr)
{
    return getnpar(0, name, "l", ptr);
}
int getparulong(char *name, unsigned long *ptr)
{
    return getnpar(0, name, "v", ptr);
}
int getparfloat(char *name, float *ptr)
{
    return getnpar(0, name, "f", ptr);
}
int getpardouble(char *name, double *ptr)
{
    return getnpar(0, name, "d", ptr);
}
int getparstring(char *name, char **ptr)
{
    return getnpar(0, name, "s", ptr);
}
int getparstringarray(char *name, char **ptr)
{
    return getnpar(0, name, "a", ptr);
}
int getpar(char *name, char *type, void *ptr)
{
    return getnpar(0, name, type, ptr);
}

/* functions to get values for the n'th occurrence of a parameter name */
int getnparint(int n, char *name, int *ptr)
{
    return getnpar(n, name, "i", ptr);
}
int getnparuint(int n, char *name, unsigned int *ptr)
{
    return getnpar(n, name, "p", ptr);
}
int getnparshort(int n, char *name, short *ptr)
{
    return getnpar(n, name, "h", ptr);
}
int getnparushort(int n, char *name, unsigned short *ptr)
{
    return getnpar(n, name, "u", ptr);
}
int getnparlong(int n, char *name, long *ptr)
{
    return getnpar(n, name, "l", ptr);
}
int getnparulong(int n, char *name, unsigned long *ptr)
{
    return getnpar(n, name, "v", ptr);
}
int getnparfloat(int n, char *name, float *ptr)
{
    return getnpar(n, name, "f", ptr);
}
int getnpardouble(int n, char *name, double *ptr)
{
    return getnpar(n, name, "d", ptr);
}
int getnparstring(int n, char *name, char **ptr)
{
    return getnpar(n, name, "s", ptr);
}
int getnparstringarray(int n, char *name, char **ptr)
{
    return getnpar(n, name, "a", ptr);
}
int getnpar(int n, char *name, char *type, void *ptr)
{
    int i;      /* index of name in symbol table	*/
    int j;      /* index for par_names[]		*/
    int nval;   /* number of parameter values found	*/
    char *aval; /* ascii field of symbol		*/

    /*--------------------------------------------------------------------*\
      getpar gets called in loops reading traces in some programs.  So
      check for having seen this name before. Also make sure we don't
      walk off the end of the table.
      \*--------------------------------------------------------------------*/

    if (parcheck && strcmp("lheader", name))
    {
        fprintf(stderr, "getpar() call after checkpars(): %s\n", name);
    }

    for (j = 0; j < par_count; j++)
    {
        if (!strcmp(par_names[j], name))
        {
            break;
        }
    }

    if (j >= par_count && par_count < PAR_NAMES_MAX)
    {
        par_names[par_count++] = name;
    }

    if (par_count == PAR_NAMES_MAX)
    {
        fprintf(stderr, " %s exceeded PAR_NAMES_MAX %d \n", xargv[0], PAR_NAMES_MAX);
    }

    if (xargc == 1)
        return 0;
    if (!tabled)
        getparinit();         /* Tabulate command line and parfile */
    i = getparindex(n, name); /* Get parameter index */
    if (i < 0)
        return 0; /* Not there */

    if (0 == ptr)
    {
        fprintf(stderr, "%s: getnpar called with 0 pointer, type = %s\n", __FILE__, type);
    }

    /*
     * handle string type as a special case, since a string
     * may contain commas.
     */
    if (type[0] == 's')
    {
        *((char **)ptr) = argtbl[i].asciival;
        return 1;
    }

    /* convert vector of ascii values to numeric values */
    for (nval = 0, aval = argtbl[i].asciival; *aval; nval++)
    {
        switch (type[0])
        {
        case 'b':
            *(bool *)ptr = eatob(aval);
            ptr = (bool *)ptr + 1;
            break;
        case 'i':
            *(int *)ptr = eatoi(aval);
            ptr = (int *)ptr + 1;
            break;
        case 'p':
            *(unsigned int *)ptr = eatop(aval);
            ptr = (unsigned int *)ptr + 1;
            break;
        case 'h':
            *(short *)ptr = eatoh(aval);
            ptr = (short *)ptr + 1;
            break;
        case 'u':
            *(unsigned short *)ptr = eatou(aval);
            ptr = (unsigned short *)ptr + 1;
            break;
        case 'l':
            *(long *)ptr = eatol(aval);
            ptr = (long *)ptr + 1;
            break;
        case 'v':
            *(unsigned long *)ptr = eatov(aval);
            ptr = (unsigned long *)ptr + 1;
            break;
        case 'f':
            *(float *)ptr = eatof(aval);
            ptr = (float *)ptr + 1;
            break;
        case 'd':
            *(double *)ptr = eatod(aval);
            ptr = (double *)ptr + 1;
            break;
        case 'a':
        {
            char *tmpstr = "";
            tmpstr = alloc1(strlen(aval) + 1, 1);

            strchop(aval, tmpstr);
            *(char **)ptr = tmpstr;
            ptr = (char **)ptr + 1;
        }
        break;
        default:
            fprintf(stderr, "%s: invalid parameter type = %s",
                    __FILE__, type);
        }
        while (*aval++ != ',')
        {
            if (!*aval)
                break;
        }
    }
    return nval;
}

void checkpars(void)
{

    int i;
    int j;
    char buf[256];

    if (getparint("verbose", &i) && i == 1)
    {

#ifdef SUXDR
        fprintf(stderr, "Using Big Endian SU data format w/ XDR.\n");
#else
        fprintf(stderr, "Using native byte order SU data format w/o XDR.\n");
#endif
    }

    for (j = 1; j < xargc; j++)
    {

        for (i = 0; i < par_count; i++)
        {
            sprintf(buf, "%s=", par_names[i]);

            if (!strncmp(buf, xargv[j], strlen(buf)))
            {
                break;
            }
        }
        if (i == par_count && strchr(xargv[j], '='))
        {
            fprintf(stderr, "Unknown %s argument %s\n", xargv[0], xargv[j]);
        }
    }

    parcheck = 1;
}

/* return number of occurrences of parameter name */
int countparname(char *name)
{
    int i, nname;

    if (xargc == 1)
        return 0;
    if (!tabled)
        getparinit();
    for (i = 0, nname = 0; i < nargs; ++i)
        if (!strcmp(name, argtbl[i].name))
            ++nname;
    return nname;
}

/* return number of values in n'th occurrence of parameter name */
int countnparval(int n, char *name)
{
    int i;

    if (xargc == 1)
        return 0;
    if (!tabled)
        getparinit();
    i = getparindex(n, name);
    if (i >= 0)
        return ccount(',', argtbl[i].asciival) + 1;
    else
        return 0;
}

/* return number of values in last occurrence of parameter name */
int countparval(char *name)
{
    return countnparval(0, name);
}

/*
 * Return the index of the n'th occurrence of a parameter name,
 * except if n==0, return the index of the last occurrence.
 * Return -1 if the specified occurrence does not exist.
 */
static int getparindex(int n, char *name)
{
    int i;
    if (n == 0)
    {
        for (i = nargs - 1; i >= 0; --i)
            if (!strcmp(name, argtbl[i].name))
                break;
        return i;
    }
    else
    {
        for (i = 0; i < nargs; ++i)
            if (!strcmp(name, argtbl[i].name))
                if (--n == 0)
                    break;
        if (i < nargs)
            return i;
        else
            return -1;
    }
}

/* Initialize getpar */
static void getparinit(void)
{
    static char *pfname; /* name of parameter file		*/
    FILE *pffd = NULL;   /* file id of parameter file		*/
    int pflen;           /* length of parameter file in bytes	*/
    int parfile;         /* parfile existence flag		*/
    int argstrlen = 0;
    char *pargstr; /* storage for parameter file args	*/
    int nread = 0; /* bytes fread				*/
    int i, j;      /* counters				*/
    int start = true;
    int debug = false;
    int valencete = false;

    tabled = true; /* remember table is built		*/

    /* Check if xargc was initiated */

    if (!xargc)
        fprintf(stderr, "%s: xargc=%d -- not initiated in main\n", __FILE__, xargc);

    /* Space needed for command lines */

    for (i = 1, argstrlen = 0; i < xargc; i++)
    {
        argstrlen += strlen(xargv[i]) + 1;
    }

    /* Get parfile name if there is one */

    if ((pfname = getpfname()))
    {
        parfile = true;
    }
    else
    {
        parfile = false;
    }

    if (parfile)
    {
        pffd = fopen(pfname, "r");

        /* Get the length */
        fseek(pffd, 0, SEEK_END);

        pflen = ftell(pffd);

        rewind(pffd);
        argstrlen += pflen;
    }
    else
    {
        pflen = 0;
    }

    /*--------------------------------------------------------------------*\
      Allocate space for command line and parameter file. The pointer
      table could be as large as the string buffer, but no larger.

      The parser logic has been completely rewritten to prevent bad
      input from crashing the program.

      Reginald H. Beardsley			    rhb@acm.org
      \*--------------------------------------------------------------------*/

    argstr = (char *)alloc1(argstrlen + 1, 1);
    targv = (char **)alloc1((argstrlen + 1) / 4, sizeof(char *));

    if (parfile)
    {
        /* Read the parfile */

        nread = fread(argstr, 1, pflen, pffd);
        if (nread != pflen)
        {
            fprintf(stderr, "%s: fread only %d bytes out of %d from %s\n",
                    __FILE__, nread, pflen, pfname);
        }
        fclose(pffd);
    }

    /* force input to valid 7 bit ASCII */

    for (i = 0; i < nread; i++)
    {
        argstr[i] &= 0x7F;
    }

    /* tokenize the input */

    j = 0;

    for (i = 0; i < nread; i++)
    {

        /* look for start of token */

        if (start)
        {

            /* getpars.c:475: warning: subscript has type `char' */
            if (isgraph((int)argstr[i]))
            {
                targv[j] = &(argstr[i]);
                start = !start;
                j++;
            }
            else
            {
                argstr[i] = 0;
            }

            /* terminate token */

            /* getpars.c:487: warning: subscript has type `char' */
        }
        else if (!valencete && isspace((int)argstr[i]))
        {
            argstr[i] = 0;
            start = !start;
        }

        /* toggle valencete semaphore */

        if (argstr[i] == '\'' || argstr[i] == '\"')
        {
            valencete = !valencete;
        }
    }

    /* display all tokens */

    if (debug)
    {

        i = 0;
        while (i < j && targv[i] != 0)
        {
            if (strlen(targv[i]))
            {
                fprintf(stderr, "%d -> %s\n", i, targv[i]);
            }
            i++;
        }
    }

    /* discard non-parameter tokens */

    i = 0;
    targc = 0;
    while (i < j && targv[i] != 0)
    {
        if (strchr(targv[i], '='))
        {
            targv[targc] = targv[i];
            targc++;
        }
        i++;
    }

    /* Copy command line arguments */

    for (j = 1, pargstr = argstr + pflen + 1; j < xargc; j++)
    {
        strcpy(pargstr, xargv[j]);
        targv[targc++] = pargstr;
        pargstr += strlen(xargv[j]) + 1;
    }

    /* Allocate space for the pointer table */

    argtbl = (pointer_table *)alloc1(targc, sizeof(pointer_table));

    /* Tabulate targv */

    tabulate(targc, targv);

    return;
}
#define PFNAME "par="
/* Get name of parameter file */
static char *getpfname(void)
{
    int i;
    size_t pfnamelen;

    pfnamelen = strlen(PFNAME);
    for (i = xargc - 1; i > 0; i--)
    {
        if (!strncmp(PFNAME, xargv[i], pfnamelen) && strlen(xargv[i]) != pfnamelen)
        {
            return xargv[i] + pfnamelen;
        }
    }
    return NULL;
}

#define iswhite(c) ((c) == ' ' || (c) == '\t' || (c) == '\n')

/* Install symbol table */
static void tabulate(size_t argc, char **argv)
{
    int i;
    char *eqptr;
    int debug = false;

    for (i = 0, nargs = 0; i < argc; i++)
    {
        eqptr = strchr(argv[i], '=');
        if (eqptr)
        {
            argtbl[nargs].name = argv[i];
            argtbl[nargs].asciival = eqptr + 1;
            *eqptr = (char)0;

            /* Debugging dump */
            if (debug)
            {
                fprintf(stderr,
                        "argtbl[%d]: name=%s asciival=%s\n",
                        nargs, argtbl[nargs].name, argtbl[nargs].asciival);
            }
            nargs++;
        }
    }
    return;
}

/* Count characters in a string */
static int ccount(char c, char *s)
{
    int i, count;
    for (i = 0, count = 0; s[i] != 0; i++)
        if (s[i] == c)
            count++;
    return count;
}
/* eatoh - convert string s to bool {true:false} */
bool eatob(char *s)
{
    bool b = false;
    char *p[] = {"YES", "yes", "Y", "y"};
    for (int i = 0; i < 4; i++)
    {
        if (strcmp(p[i], s) == 0)
        {
            b = true;
            break;
        }
    }
    return b;
}
/* eatoh - convert string s to short integer {SHRT_MIN:SHRT_MAX} */
short eatoh(char *s)
{
    long n = strtol(s, NULL, 10);
    return (short)n;
}
/* eatou - convert string s to unsigned short integer {0:USHRT_MAX} */
unsigned short eatou(char *s)
{
    unsigned long n = strtoul(s, NULL, 10);
    return (unsigned short)n;
}
/* eatoi - convert string s to integer {INT_MIN:INT_MAX} */
int eatoi(char *s)
{
    long n = strtol(s, NULL, 10);
    return (int)n;
}
/* eatop - convert string s to unsigned integer {0:UINT_MAX} */
unsigned int eatop(char *s)
{
    unsigned long n = strtoul(s, NULL, 10);
    return (unsigned int)n;
}
/* eatol - convert string s to long integer {LONG_MIN:LONG_MAX} */
long eatol(char *s)
{
    long n = strtol(s, NULL, 10);
    return n;
}
/* eatov - convert string s to unsigned long {0:ULONG_MAX} */
unsigned long eatov(char *s)
{
    unsigned long n = strtoul(s, NULL, 10);
    return n;
}
/* eatof - convert string s to float {-FLT_MAX:FLT_MAX} */
float eatof(char *s)
{
    float x = strtod(s, NULL);
    return (float)x;
}
/* eatod - convert string s to double {-DBL_MAX:DBL_MAX} */
double eatod(char *s)
{
    double x = strtod(s, NULL);
    return x;
}

/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */
/* errpkge.c
   err	 print warning on application program error and die
   warn print warning on application program error
   Examples:
   err("Cannot divide %f by %f", x, y);
   warn("fmax = %f exceeds half nyquist= %f", fmax, 0.25/dt);

   if (NULL == (fp = fopen(xargv[1], "r")))
   err("can't open %s", xargv[1]);
   ...
   if (-1 == close(fd))
   err("close failed");
*/
void err(char *fmt, ...)
{
    va_list args;

    if (EOF == fflush(stdout))
    {
        fprintf(stderr, "\nerr: fflush failed on stdout");
    }
    fprintf(stderr, "\n%s: ", xargv[0]);
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fprintf(stderr, "\n");
    exit(EXIT_FAILURE);
}

void warn(char *fmt, ...)
{
    va_list args;

    if (EOF == fflush(stdout))
    {
        fprintf(stderr, "\nwarn: fflush failed on stdout");
    }
    fprintf(stderr, "\n%s: ", xargv[0]);
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fprintf(stderr, "\n");
    return;
}

/* pading for 2D model */
void pad2(const float *x, float *xx, const int nz, const int nx, const int lft, const int rht, const int top, const int bot)
{
    int nxb = nx + lft + rht;
    int nzb = nz + top + bot;
    int i, j;
    for (int ix = lft; ix < nx + lft; ix++)
    {
        /*inner*/
        for (int iz = top; iz < nz + top; iz++)
        {
            i = ix * nzb + iz;
            j = (ix - lft) * nz + (iz - top);
            xx[i] = x[j];
        }
        /*top*/
        j = (ix - lft) * nz;
        for (int iz = 0; iz < top; iz++)
        {
            i = ix * nzb + iz;
            xx[i] = x[j];
        }
        /*bot*/
        j = (ix - lft) * nz + (nz - 1);
        for (int iz = nz + top; iz < nzb; iz++)
        {
            i = ix * nzb + iz;
            xx[i] = x[j];
        }
    }
    for (int iz = 0; iz < nzb; iz++)
    {
        /*left*/
        j = lft * nzb + iz;
        for (int ix = 0; ix < lft; ix++)
        {
            i = ix * nzb + iz;
            xx[i] = xx[j];
        }
        /*right*/
        j = (lft + nx - 1) * nzb + iz;
        for (int ix = lft + nx; ix < nxb; ix++)
        {
            i = ix * nzb + iz;
            xx[i] = xx[j];
        }
    }
}
/* blas function used openblas*/
/* blas2 */
void sscale(int n, float alpha, float *a, float *b)
{
    memcpy(b, a, sizeof(float) * n);
    cblas_sscal(n, alpha, b, 1);
}
float sdot(int n, float *a, float *b)
{
    return cblas_sdot(n, a, 1, b, 1);
}
double ddot(int n, double *a, double *b)
{
    return cblas_ddot(n, a, 1, b, 1);
}
void saxpy(int n, float alpha, float *a, float *b)
{
    cblas_saxpy(n, alpha, a, 1, b, 1);
}
/* blas3 */
void sgemm(int m, int n, int k, float *A, float *B, float *C)
{
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1, A, m, B, k, 0, C, m);
}
void dgemm(int m, int n, int k, double *A, double *B, double *C)
{
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1, A, m, B, k, 0, C, m);
}