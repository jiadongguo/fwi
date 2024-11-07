#include "cstd.h"
#include "waveutils.h"
#include <mpi/mpi.h>
void eal_init(acpar par, float alpha_, int mode_, float *vv);
void eal_apply(acpar par, float *pre, float *curr, float *next, float *vv);
void eal_close();
float laplace(int n1, int n2, int i1, int i2, float *curr, float d1, float d2);
void bound(float *v, int n, float max, float min)
{
    for (int i = 0; i < n; i++)
    {
        v[i] = MAX(v[i], min);
        v[i] = MIN(v[i], max);
    }
}
void fdfor(acpar par, float *pre, float *curr, float *next, float *vv, float *lap)
{
    int nz, nx, nzb, nxb, top, lft;
    float dz, dx, dt;
    float tmp;
    nz = par->nz;
    nx = par->nx;
    nzb = par->nzb;
    nxb = par->nxb;
    top = par->top;
    lft = par->lft;
    dz = par->dz;
    dx = par->dx;
    dt = par->dt;
    /*only for inner grid*/
    for (int ix = lft; ix < lft + nx; ix++)
    {
        for (int iz = top; iz < top + nz; iz++)
        {
            tmp = laplace(nzb, nxb, iz, ix, curr, dz, dx);
            if (lap != NULL)
            {
                lap[(ix - lft) * nz + (iz - top)] = tmp;
            }
            next[ix * nzb + iz] = tmp * (vv[ix * nzb + iz] * dt) * (vv[ix * nzb + iz] * dt) + 2 * curr[ix * nzb + iz] - pre[ix * nzb + iz];
        }
    }
    eal_apply(par, pre, curr, next, vv);
}
float findalpha(acpar par, float vmax, float *grad, float refl)
{
    float gmax = 0;
    int nzx = par->nzx;
    for (int ix = 0; ix < nzx; ix++)
    {
        gmax = MAX(gmax, grad[ix]);
    }
    return refl * vmax / (gmax + TOL);
}
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    initargs(argc, argv);
    bool verb, ill;
    int nz, nx, nt, top, bot, lft, rht;
    int ns, sz, sx, jsx, jsz, rz, rx, jrx, jrz, nr;
    float dz, dx, dt;
    char *fwt, *fvel, *out, *fobs, *fobj;
    float *wt, *vel;
    int mode;
    int niter;
    float max, min;
    if (!getparbool("verb", &verb))
        verb = true;
    if (!getparbool("illum", &ill))
        ill = true;
    if (!getparint("mode", &mode))
        mode = 0;
    if (!getparint("niter", &niter))
        niter = 10;
    if (!getparstring("fobs", &fobs))
        err("need fobs");
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
    if (!getparfloat("min", &min))
        err("need min");
    if (!getparfloat("max", &max))
        err("need max");
    if (!getparstring("vpfile", &fvel))
        err("need vp");
    if (!getparstring("wtfile", &fwt))
        err("need fwt");
    if (!getparstring("fobj", &fobj))
        err("need fobj");
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
        FILE *fd = fopen(fvel, "rb");
        fread(vel, sizeof(float) * nz * nx, 1, fd);
        fclose(fd);
    }
    /*-----------------------------------------------read wavelet--------------------------------------------*/
    {
        FILE *fd = fopen(fwt, "rb");
        fread(wt, sizeof(float) * nt, 1, fd);
        fclose(fd);
    }
    acpar par = creat_acpar(nz, nx, dz, dx, top, bot, lft, rht, nt, dt, ns, sz, sx, jsx, jsz, nr, rz, rx, jrx, jrz);
    int ns0 = ns;
    int nzb, nxb, nzxb, nzx;
    float *pre, *curr, *next, *tmp;
    float *dcal, *dobs;
    float *illum;
    nzb = par->nzb, nxb = par->nxb, nzxb = nzb * nxb, nzx = nz * nx;
    if (ill)
        illum = alloc1float(nzx);
    pre = alloc1float(nzxb);
    curr = alloc1float(nzxb);
    next = alloc1float(nzxb);
    dcal = alloc1float(nr * nt);
    dobs = alloc1float(nr * nt);
    if (ns % size != 0)
    {
        ns += size - ns % size;
    }
    int is_x, is_z, ir_x, ir_z;
    float *vv = alloc1float(nzxb);
    float *grad = alloc1float(nzx);
    float *gradient;
    float *lap = alloc1float(nzx * nt);
    float alpha = 0.1;
    float objsum, obj;
    float vmax = 0;
    FILE *Fobj;
    if (rank == 0)
    {
        Fobj = fopen(fobj, "w");
        gradient = alloc1float(nzx);
        for (int ix = 0; ix < nzx; ix++)
        {
            vmax = MAX(vmax, vel[ix]);
        }
    }
    clock_t start, finish;
    eal_init(par, 1e-4, mode, vv);
    MPI_File fh;
    MPI_Status status;
    MPI_Barrier(MPI_COMM_WORLD);
    for (int iter = 0; iter < niter; iter++)
    {
        start = clock();
        obj = 0;
        objsum = 0;
        pad2(vel, vv, nz, nx, lft, rht, top, bot);
        memset(grad, 0, sizeof(float) * nzx);
        if (rank == 0)
        {
            memset(gradient, 0, sizeof(float) * nzx);
        }
        eal_init(par, 1e-4, 0, vv);
        for (int is = rank; is < ns; is += size)
        {
            if (is < ns0)
            {
                memset(pre, 0, sizeof(float) * nzxb);
                memset(curr, 0, sizeof(float) * nzxb);
                memset(next, 0, sizeof(float) * nzxb);
                is_x = lft + sx + (is - 1) * jsx, is_z = top + sz + (is - 1) * jsz;
                MPI_File_open(MPI_COMM_WORLD, fobs, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
                MPI_File_read_at(fh, is * nt * nr * sizeof(float), dobs, nt * nr, MPI_FLOAT, &status);
                for (int it = 0; it < nt; it++)
                {
                    // if (verb && it % 100 == 0)
                    // {
                    //     printf("forward modeling is=%d/%d,it=%d/%d\n", is + 1, ns0, it + 1, nt);
                    // }
                    fdfor(par, pre, curr, next, vv, lap + it * nzx);
                    if (ill)
                    {
                        for (int ix = 0; ix < nx; ix++)
                        {
                            for (int iz = 0; iz < nz; iz++)
                            {
                                illum[ix * nz + iz] = pow(curr[(ix + lft) * nzb + iz + top], 2);
                            }
                        }
                    }
                    tmp = pre, pre = curr, curr = next, next = tmp;
                    curr[is_x * nzb + is_z] += wt[it];
                    for (int ir = 0; ir < nr; ir++)
                    {
                        ir_x = lft + rx + (ir - 1) * jrx, ir_z = top + rz + (ir - 1) * jrz;
                        dcal[ir * nt + it] = curr[ir_x * nzb + ir_z];
                        obj += pow(dcal[ir * nt + it] - dobs[ir * nt + it], 2) / 2;
                    }
                }
                MPI_File_close(&fh);
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Reduce(&obj, &objsum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
                memset(pre, 0, sizeof(float) * nzxb);
                memset(curr, 0, sizeof(float) * nzxb);
                memset(next, 0, sizeof(float) * nzxb);
                for (int it = nt - 1; it >= 0; it--)
                {
                    fdfor(par, pre, curr, next, vv, NULL);
                    tmp = pre, pre = curr, curr = next, next = tmp;
                    for (int ir = 0; ir < nr; ir++)
                    {
                        ir_x = lft + rx + (ir - 1) * jrx, ir_z = top + rz + (ir - 1) * jrz;
                        curr[ir_x * nzb + ir_z] += dcal[ir * nt + it] - dobs[ir * nt + it];
                    }
                    for (int ix = lft; ix < nx + lft; ix++)
                    {
                        for (int iz = top + 10; iz < nz + top; iz++)
                        {
                            if (ill)
                                grad[(ix - lft) * nz + (iz - top)] += 2 * curr[ix * nzb + iz] * lap[it * nz * nx + (ix - lft) * nz + (iz - top)] / vv[ix * nzb + iz] / (illum[(ix - lft) * nz + (iz - top)] + EPS);
                            else
                                grad[(ix - lft) * nz + (iz - top)] += 2 * curr[ix * nzb + iz] * lap[it * nz * nx + (ix - lft) * nz + (iz - top)] / vv[ix * nzb + iz];
                        }
                    }
                }
            }
        }
        eal_close();
        if (rank == 0)
        {
            fprintf(Fobj, "\t%d\t%g\n", iter + 1, objsum);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        /*gradient*/
        MPI_Reduce(grad, gradient, nz * nx, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        /*alpha*/
        if (rank == 0)
        {
            if (iter < 4)
                alpha = findalpha(par, vmax, gradient, 0.01);
            else if (iter < 100)
            {
                alpha = findalpha(par, vmax, gradient, 0.001);
            }
            else
            {
                alpha = findalpha(par, vmax, gradient, 0.0001);
            }
            alpha = 1e-2;
            for (int ix = 0; ix < nzx; ix++)
            {
                vel[ix] -= alpha * gradient[ix];
            }
            bound(vel, nzx, max, min);
        }
        MPI_Bcast(vel, nzx, MPI_FLOAT, 0, MPI_COMM_WORLD);
        finish = clock();
        if (rank == 0)
        {
            warn("iter%d finish,obj=%g,cost time %g", iter + 1, objsum, 1. * (finish - start) / CLOCKS_PER_SEC);
        }
    }
    if (rank == 0)
    {
        FILE *fd = fopen(out, "wb");
        fwrite(vel, sizeof(float), nzx, fd);
        fclose(fd);
        fclose(Fobj);
    }
    MPI_Finalize();
    return 0;
}
