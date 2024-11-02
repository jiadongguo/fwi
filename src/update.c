#include "cstd.h"
int main(int argc, char **argv)
{
    int n;
    float alpha, *dir, *old;
    char *fin, *fout, *fdir, *falpha;
    bool neg;
    initargs(argc, argv);
    if (!getparstring("falpha", &falpha))
        err("need falpha");
    if (!getparstring("in", &fin))
        err("need fin");
    if (!getparstring("out", &fout))
        err("need fout");
    if (!getparstring("fdir", &fdir))
        err("need fdir");
    if (!getparint("n", &n))
        err("need n");
    if (!getparbool("neg", &neg))
        neg = true; /*true: need opposite direction*/
    dir = alloc1float(n);
    old = alloc1float(n);
    {
        FILE *fd = fopen(fin, "rb");
        if (fd == NULL)
            err("can't open fin");
        fread(old, sizeof(float) * n, 1, fd);
        fclose(fd);
    }
    {
        FILE *fd = fopen(fdir, "rb");
        if (fd == NULL)
            err("can't open fdir");
        fread(dir, sizeof(float) * n, 1, fd);
        fclose(fd);
    }

    {
        FILE *fd = fopen(falpha, "rb");
        if (fd == NULL)
            err("can't open falpha");
        fread(&alpha, sizeof(float), 1, fd);
        fclose(fd);
    }
    for (int i = 0; i < n; i++)
    {
        if (neg)
        {
            old[i] -= alpha * dir[i];
        }
        else
        {
            old[i] += alpha * dir[i];
        }
    }
    {
        FILE *fd = fopen(fdir, "wb");
        if (fd == NULL)
            err("can't open fout");
        fwrite(old, sizeof(float) * n, 1, fd);
        fclose(fd);
    }
    free1float(dir);
    free1float(old);
}