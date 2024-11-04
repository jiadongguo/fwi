#include "cstd.h"
int main(int argc, char **argv)
{
    initargs(argc, argv);
    int nz, nx;
    float *v;
    if (!getparint("n1", &nz))
        err("need nz");
    if (!getparint("n2", &nx))
        err("need nx");
    v = alloc1float(nz * nx);
    for (int ix = 0; ix < nx; ix++)
    {
        for (int iz = 0; iz < nz; iz++)
        {

            if (iz > nz / 2)
            {
                v[ix * nz + iz] = 5000;
            }
            else
            {
                v[ix * nz + iz] = 4000;
            }
        }
    }
    char *fvel;
    if (!getparstring("vpfile", &fvel))
        err("need vpfile");
    FILE *fd = fopen(fvel, "wb");
    fwrite(v, sizeof(float) * nz * nx, 1, fd);
    fclose(fd);
    free1float(v);
}