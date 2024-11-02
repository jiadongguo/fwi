#include "cstd.h"

static float *pp, *qq, *bb1, *bb2, *yy;
static int np, nq;

static void triangle_init(int nbox, int nd)
{
    np = nbox + nd - 1;
    nq = nbox + np - 1;
    pp = alloc1float(np);
    qq = alloc1float(nq);
    bb1 = alloc1float(np);
    bb2 = alloc1float(nq);
    yy = alloc1float(nd);
}

static void triangle_close()
{
    free1float(pp);
    free1float(qq);
    free1float(bb1);
    free1float(bb2);
    free1float(yy);
}

static void boxconv_lop(int nbox, int nx, float *xx, float *yy, float *bb)
{
    int ny;
    int i;

    ny = nx + nbox - 1;
    for (i = 0; i < ny; i++)
        bb[i] = 0.;
    bb[0] = xx[0];
    for (i = 1; i < nx; i++)
        bb[i] = bb[i - 1] + xx[i]; // make B(z) = X(z)/(1-z)
    for (i = nx; i < ny; i++)
        bb[i] = bb[i - 1];
    for (i = 0; i < nbox; i++)
        yy[i] = bb[i];
    for (i = nbox; i < ny; i++)
        yy[i] = bb[i] - bb[i - nbox]; // make Y(z) = B(z)*(1-z**nbox)
    for (i = 0; i < ny; i++)
        yy[i] /= nbox;
}

static void triangle_lop(int nbox, int nd, float *xx)
{
    int i;

    boxconv_lop(nbox, nd, xx, pp, bb1);
    boxconv_lop(nbox, np, pp, qq, bb2);
    for (i = 0; i < nd; i++)
        yy[i] = qq[i + nbox - 1];
    for (i = 0; i < nbox - 1; i++)
        yy[i] += qq[nbox - 2 - i]; // fold back near end
    for (i = 0; i < nbox - 1; i++)
        yy[nd - i - 1] += qq[nd + (nbox - 1) + i]; // fold back far end
    memcpy(xx, yy, nd * sizeof(float));
}

void triangle_smoothing(float *mod, int n1, int n2, int r1, int r2, int repeat)
{
    int i, i1, i2;
    float *tmp;

    tmp = alloc1float(n2);

    for (i = 0; i < repeat; i++)
    {
        /*-----------------------------------------*/
        triangle_init(r1, n1);
        for (i2 = 0; i2 < n2; i2++)
            triangle_lop(r1, n1, mod + i2 * n1);
        triangle_close();

        /*-----------------------------------------*/
        triangle_init(r2, n2);

        for (i1 = 0; i1 < n1; i1++)
        {
            for (i2 = 0; i2 < n2; i2++)
                tmp[i2] = mod[i2 * n1 + i1];
            triangle_lop(r2, n2, tmp);
            for (i2 = 0; i2 < n2; i2++)
                mod[i2 * n1 + i1] = tmp[i2];
        }
        triangle_close();
    }
    free1float(tmp);
}
int main(int argc, char **argv)
{
    initargs(argc, argv);
    int nz, nx;
    int r1, r2;
    float *v;
    char *s, *out;
    int repeat;
    if (!getparint("repeat", &repeat))
    {
        repeat = 1;
    }
    if (!getparint("n1", &nz))
    {
        err("need nz");
    }
    if (!getparint("n2", &nx))
    {
        err("need nx");
    }
    if (!getparint("r1", &r1))
    {
        err("need r1");
    }
    if (!getparint("r2", &r2))
    {
        err("need r2");
    }
    if (!getparstring("vpfile", &s))
    {
        err("need vpfile");
    }
    if (!getparstring("out", &out))
    {
        err("need output file");
    }
    v = alloc1float(nz * nx);
    {
        FILE *fd = fopen(s, "rb");
        if (fd == NULL)
        {
            err("can't open velocity file");
        }
        fread(v, sizeof(float) * nz * nx, 1, fd);
        fclose(fd);
    }
    triangle_smoothing(v, nz, nx, r1, r2, repeat);
    {
        FILE *fd = fopen(out, "wb");
        if (fd == NULL)
        {
            err("can't open output file");
        }
        fwrite(v, sizeof(float) * nz * nx, 1, fd);
        fclose(fd);
    }
}