/*
seismic wavefield modeling with effective absorbing layer
*/
#include "cstd.h"
#include "waveutils.h"
#include "fdutils.h"
#include <mpi/mpi.h>

void eal_init(acpar par, float alpha_, int mode_, float *vv_);
void eal_apply(acpar par, float *pre, float *curr, float *next, float *vv_);
void eal_close();

int main(int argc, char **argv)
{
    initargs(argc, argv);
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int nz, nx, nt, top, bot, lft, rht;
    int ns, sz, sx, jsx, jsz, rz, rx, jrx, jrz, nr;
    float dz, dx, dt;
    char *fwt, *fvel, *out;
    float *wt, *vel;
    int mode;
    if (!getparint("mode", &mode))
    {
        mode = 0;
        warn("set mode=0");
    }
    if (!getparint("n1", &nz))
        err("need nz");
    if (!getparint("n2", &nx))
        err("need nx");
    if (!getparint("nt", &nt))
        err("need nt");
    if (!getparint("top", &top))
        err("need top");
    if (!getparint("bot", &bot))
        err("need bot");
    if (!getparint("lft", &lft))
        err("need lft");
    if (!getparint("rht", &rht))
        err("need rht");
    if (!getparint("ns", &ns))
        err("need ns");
    if (!getparint("nr", &nr))
        err("need nr");
    if (!getparint("sz", &sz))
        err("need sz");
    if (!getparint("sx", &sx))
        err("need sx");
    if (!getparint("jsx", &jsx))
        err("need jsx");
    if (!getparint("jsz", &jsz))
        err("need jsz");
    if (!getparint("rz", &rz))
        err("need rz");
    if (!getparint("rx", &rx))
        err("need rx");
    if (!getparint("jrx", &jrx))
        err("need jrx");
    if (!getparint("jrz", &jrz))
        err("need jrz");
    if (!getparfloat("dt", &dt))
        err("need dt");
    if (!getparfloat("d2", &dx))
        err("need dx");
    if (!getparfloat("d1", &dz))
        err("need dz");
    if (!getparstring("vpfile", &fvel))
        err("need vp");
    if (!getparstring("wtfile", &fwt))
        err("need fwt");
    if (!getparstring("out", &out))
        err("need out");
    if ((sx + (ns - 1) * jsx >= nx) || (sz + (ns - 1) * jsz >= nz))
        err("source position outer bound");
    if ((rx + (nr - 1) * jrx >= nx) || (rz + (nr - 1) * jrz >= nz))
        err("receiver position outer bound");
    /*--------------------------------------------------------------------------------------------------------------*/

    vel = alloc1float(nz * nx);
    wt = alloc1float(nt);
    /*-----------------------------------------------read velocity model--------------------------------------------*/
    {
        FILE *fd = input(fvel);
        fread(vel, sizeof(float) * nz * nx, 1, fd);
        fclose(fd);
    }
    /*-----------------------------------------------read wavelet--------------------------------------------*/
    {
        FILE *fd = input(fwt);
        fread(wt, sizeof(float) * nt, 1, fd);
        fclose(fd);
    }
    acpar par = creat_acpar(nz, nx, dz, dx, top, bot, lft, rht, nt, dt, ns, sz, sx, jsx, jsz, nr, rz, rx, jrx, jrz);
    int ns0 = ns;
    int nzb, nxb, nzxb;
    float *pre, *curr, *next, *tmp;
    float *rcd, *seis;
    FILE *fd;
    nzb = par->nzb, nxb = par->nxb, nzxb = nzb * nxb;

    pre = alloc1float(nzxb);
    curr = alloc1float(nzxb);
    next = alloc1float(nzxb);
    rcd = alloc1float(nr * nt);

    if (ns % size != 0)
    {
        ns += size - ns % size;
    }

    if (rank == 0)
    {
        fd = output(out);
        seis = alloc1float(size * nr * nt);
    }

    float *vpd = alloc1float(nzb * nxb);
    pad2(vel, vpd, nz, nx, lft, rht, top, bot);
    eal_init(par, 1e-4, mode, vpd);
    MPI_Barrier(MPI_COMM_WORLD);
    clock_t start = clock();
    for (int is = rank; is < ns; is += size)
    {
        if (is < ns0)
        {
            memset(pre, 0, sizeof(float) * nzxb);
            memset(curr, 0, sizeof(float) * nzxb);
            memset(next, 0, sizeof(float) * nzxb);
            for (int it = 0; it < nt; it++)
            {
                if (it % 100 == 0)
                {
                    printf("forward modeling is=%d/%d,it=%d/%d\n", is + 1, ns0, it + 1, nt);
                }
                step_forward(par, pre, curr, next, vpd);
                tmp = pre, pre = curr, curr = next, next = tmp;
                add_src(par, curr, wt[it], is);
                record(par, curr, rcd + it * nr);
            }
            transp(nr, nt, rcd, NULL);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gather(rcd, nr * nt, MPI_FLOAT, seis, nr * nt, MPI_FLOAT, 0, MPI_COMM_WORLD);
        if (rank == 0)
        {
            if (is + size >= ns0)
            {
                fwrite(seis, sizeof(float), (ns0 - is) * nr * nt, fd);
            }
            else
            {
                fwrite(seis, sizeof(float), nr * nt * size, fd);
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    clock_t finish = clock();
    if (rank == 0)
    {
        printf("Forward modeling costs time %.4fs\n", 1. * (finish - start) / CLOCKS_PER_SEC);
        free1float(seis);
        fclose(fd);
    }
    eal_close();
    free1float(rcd);
    free1float(pre);
    free1float(curr);
    free1float(next);
    free1float(wt);
    free1float(vel);
    free(par);
    MPI_Finalize();
    return 0;
}