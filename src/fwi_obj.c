#include "cstd.h"
int main(int argc, char **argv)
{
    initargs(argc, argv);
    int ns, nt, ng;
    char *out, *dcal, *dobs;
    FILE *fobs, *fcal, *ferr;
    float *pcal, *pobs, obj;
    if (!getparint("ns", &ns))
        err("need ns");
    if (!getparint("nt", &nt))
        err("need nt");
    if (!getparint("ng", &ng))
        err("need ng");
    if (!getparstring("out", &out))
        err("need out");
    if (!getparstring("fcal", &dcal))
        err("need fcal");
    if (!getparstring("fobs", &dobs))
        err("need fobs");
    fobs = fopen(dobs, "rb");
    fcal = fopen(dcal, "rb");
    ferr = fopen(out, "w+");
    pobs = alloc1float(nt * ng);
    pcal = alloc1float(nt * ng);
    obj = 0;
    for (int is = 0; is < ns; is++)
    {
        fread(pcal, sizeof(float) * ng * nt, 1, fcal);
        fread(pobs, sizeof(float) * ng * nt, 1, fobs);
        for (int i = 0; i < nt * ng; i++)
        {
            obj += pow(pobs[i] - pcal[i], 2) / 2;
        }
    }
    fprintf(ferr, "%g\n", obj);
    fclose(ferr);
    fclose(fobs);
    fclose(fcal);
    free1float(pcal);
    free1float(pobs);
    return 0;
}