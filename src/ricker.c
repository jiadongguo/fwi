#include "cstd.h"
void ricker(const float dt, const float fm, const int nt, float amp, const float t0, float *wt)
{
    float t;
    for (int i = 0; i < nt; i++)
    {
        t = i * dt - t0;
        t = PI2 * pow(t * fm, 2);
        wt[i] = amp * (1. - 2. * t) * exp(-t);
    }
}

int main(int argc, char **argv)
{
    initargs(argc, argv);
    int nt;
    float dt, fm, t0, amp;
    char *out;
    float *wt;
    if (!getparint("n1", &nt))
        err("need n1");
    if (!getparfloat("d1", &dt))
        err("need d1");
    if (!getparfloat("amp", &amp))
        amp = 10;
    if (!getparfloat("fm", &fm))
        err("need fm");
    if (!getparfloat("t0", &t0))
        t0 = 1. / fm;
    if (!getparstring("out", &out))
        err("need out");
    FILE *fd = output(out);
    wt = alloc1float(nt);
    ricker(dt, fm, nt, amp, t0, wt);
    fwrite(wt, sizeof(float) * nt, 1, fd);
    fclose(fd);
    free1float(wt);
}