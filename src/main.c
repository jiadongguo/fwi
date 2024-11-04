/*
a simple program for fwi
*/

#include "cstd.h"
#include "waveutils.h"

/* ======================================================================================

fun_fun_laplace : fun_laplace opeprator
fun_forward : Seismic waves travel forward
fun_backward : Seismic waves travel backward
fun_obj : objective function of fwi
fun_record :
fun_gradient : gradient of current model

===========================================================================================*/

void eal_init(acpar par, float alpha_, int mode);
void eal_apply(acpar par, float *pre, float *curr, float *next);
void eal_close();
float fun_laplace(int n1, int n2, int i1, int i2, float *curr, float d1, float d2)
{
    static float c0 = -5. / 2, c1 = 4. / 3, c2 = -1. / 12;
    float fdx[5], fdz[5];
    float lapz, lapx;
    d1 = 1. / pow(d1, 2);
    d2 = 1. / pow(d2, 2);
    for (int i = 0; i < 5; i++)
    {
        fdx[i] = ((i2 + i - 2 < 0) || (i2 + i - 2 >= n2)) ? 0 : curr[(i2 + i - 2) * n1 + i1];
        fdz[i] = ((i1 + i - 2 < 0) || (i1 + i - 2 >= n1)) ? 0 : curr[i2 * n1 + i1 + i - 2];
    }
    lapz = (c2 * (fdz[0] + fdz[4]) + c1 * (fdz[1] + fdz[3]) + c0 * fdz[2]) * d1;
    lapx = (c2 * (fdx[0] + fdx[4]) + c1 * (fdx[1] + fdx[3]) + c0 * fdx[2]) * d2;
    return lapx + lapz;
}

void fun_forward(acpar par, float *pre, float *curr, float *next, float *lap)
{
    float *vv;
    int nz, nx, nzb, nxb, top, lft;
    float dz, dx, dt;
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
    float tmp;
    for (int ix = lft; ix < lft + nx; ix++)
    {
        for (int iz = top; iz < top + nz; iz++)
        {
            tmp = fun_laplace(nzb, nxb, iz, ix, curr, dz, dx);
            if (lap != NULL)
            {
                lap[(ix - lft) * nz + (iz - top)] = tmp;
            }
            next[ix * nzb + iz] = tmp * (vv[ix * nzb + iz] * dt) * (vv[ix * nzb + iz] * dt) + 2 * curr[ix * nzb + iz] - pre[ix * nzb + iz];
        }
    }
    eal_apply(par, pre, curr, next);
}
void fun_backward(acpar par, float *pre, float *curr, float *next)
{
    float *vv;
    int nz, nx, nzb, nxb, top, lft;
    float dz, dx, dt;
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
            next[ix * nzb + iz] = fun_laplace(nzb, nxb, iz, ix, curr, dz, dx) * (vv[ix * nzb + iz] * dt) * (vv[ix * nzb + iz] * dt) + 2 * curr[ix * nzb + iz] - pre[ix * nzb + iz];
        }
    }
    eal_apply(par, pre, curr, next);
}
void fun_record(acpar par, float *curr, float *rcd)
{
    int nr = par->nr, nzb = par->nzb;
    int rx, rz;
    for (int ir = 0; ir < nr; ir++)
    {
        rx = par->rx + par->lft + par->jrx * ir;
        rz = par->rz + par->top + par->jrz * ir;
        rcd[ir] = curr[rx * nzb + rz];
    }
}
void fun_addsrc(bool adj, acpar par, float wt, float *p, int is /* is or ir */)
{
    int nzb = par->nzb;
    int lft = par->lft, top = par->top;
    int jsx, jsz, sx, sz;
    if (adj)
    {
        sx = par->rx;
        sz = par->rz;
        jsx = par->jrx;
        jsz = par->jrz;
    }
    else
    {
        sx = par->sx;
        sz = par->sz;
        jsx = par->jsx;
        jsz = par->jsz;
    }
    sx = sx + lft + is * jsx;
    sz = sz + top + is * jsz;
    p[sx * nzb + sz] += wt;
}
float fun_alphatest(acpar par, float *grad)
{
    int n = par->nzx;
    float maxv = 0, maxg = 0;
    for (int i = 0; i < n; i++)
    {
        maxv = MAX(maxv, fabs(v[i]));
        maxg = MAX(maxg, fabs(grad[i]));
    }
    return 0.01 * maxv / maxg;
}

int main(int argc, char **argv)
{
    initargs(argc, argv);
    int nter;                              /* total number of iterations */
    char *fwt, *fvel, *shots, *out, *fobj; /*observed seismogram*/
    int nz, nx, nt, top, bot, lft, rht;
    int ns, sz, sx, jsx, jsz, rz, rx, jrx, jrz, nr;
    float dz, dx, dt;
    float *wt /* wavelet */, *vel /* initial velocity model */;
    float *vtmp;
    bool verb;
    int cut;
    if (!getparbool("verb", &verb))
        verb = false;
    if (!getparint("nter", &nter))
    {
        nter = 10;
    }
    if (!getparint("n1", &nz))
        err("need nz");
    if (!getparint("cut", &cut))
    {
        cut = 10;
        warn("set cut=10");
    }
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
    if (!getparstring("obj", &fobj))
    {
        err("need obj");
    }
    if (!getparstring("shots", &shots))
        err("need shots for observed seismogram");
    if (!getparstring("vpfile", &fvel))
        err("need fvel");
    if (!getparstring("wtfile", &fwt))
        err("need fwt");
    if (!getparstring("out", &out))
        err("need out");
    if ((sx + (ns - 1) * jsx >= nx) || (sz + (ns - 1) * jsz >= nz))
        err("source position outer bound");
    if ((rx + (nr - 1) * jrx >= nx) || (rz + (nr - 1) * jrz >= nz))
        err("receiver position outer bound");
    /* ===================================================================================================== */

    float alpha;
    float *tmp, *dobs, *lap, *dcal, *grad, *dcaltmp;
    wt = alloc1float(nt);
    {
        FILE *fd = fopen(fwt, "rb");
        if (fd == NULL)
            err("can't open wavelet file");
        fread(wt, sizeof(float) * nt, 1, fd);
        fclose(fd);
    }
    vel = alloc1float(nz * nx);
    {
        FILE *fd = fopen(fvel, "rb");
        if (fd == NULL)
            err("can't open initial velocity file");
        fread(vel, sizeof(float) * nz * nx, 1, fd);
        fclose(fd);
    }
    vtmp = alloc1float(nz * nx);
    memcpy(vtmp, vel, sizeof(float) * nz * nx);
    dobs = alloc1float(nt * nr * ns);
    {
        FILE *fd = fopen(shots, "rb");
        if (fd == NULL)
            err("can't open shots file");
        fread(dobs, sizeof(float) * nr * nt * ns, 1, fd);
        fclose(fd);
    }
    /*========================================transpose for observed seismogram =================================================*/
    {
        float *dtobs = alloc1float(nt * nr);
        for (int is = 0; is < ns; is++)
        {
            transp(nt, nr, &dobs[is * nr * nt], dtobs);
            memcpy(dtobs, &dobs[is * nr * nt], sizeof(float) * nr * nt);
        }
        free1(dtobs);
    }
    /*=========================================================================================*/
    dcal = alloc1float(nt * nr);
    dcaltmp = alloc1float(nr * nt);

    /* ======================================================observation system definition======================================================================================================= */
    acpar partmp = creat_acpar(nz, nx, dz, dx, top, bot, lft, rht, nt, dt, ns, sz, sx, jsx, jsz, nr, rz, rx, jrx, jrz, vtmp); /* test model */
    acpar par = creat_acpar(nz, nx, dz, dx, top, bot, lft, rht, nt, dt, ns, sz, sx, jsx, jsz, nr, rz, rx, jrx, jrz, vel);
    int nzxb = par->nzxb, nzb = par->nzb, nxb = par->nxb, nzx = par->nzx;
    lap = alloc1float(nzx * nt);
    float *p0, *p1, *p2, *illum;
    p0 = alloc1float(nzxb);
    p1 = alloc1float(nzxb);
    p2 = alloc1float(nzxb);
    grad = alloc1float(nzx);
    illum = alloc1float(nzx);
    clock_t start, finish;
    double num1, num2;
    /*====================================================================================================================*/
    FILE *Fobj = fopen(fobj, "w");
    for (int iter = 0; iter < nter; iter++)
    {
        if (verb)
        {
            start = clock();
        }
        float obj = 0;
        memset(grad, 0, sizeof(float) * nzx);
        memset(illum, 0, sizeof(float) * nzx);
        eal_init(par, 1e-3, 0);
        for (int is = 0; is < ns; is++)
        {
            memset(p0, 0, sizeof(float) * nzxb);
            memset(p1, 0, sizeof(float) * nzxb);
            memset(p2, 0, sizeof(float) * nzxb);
            /* ====================================forward modeling ===========================*/

            fun_forward(par, p0, p1, p2, NULL);
            fun_addsrc(false, par, wt[0], p1, is);
            fun_record(par, p1, dcal + is * nr * nt + 0 * nr);
            tmp = p0, p0 = p1, p1 = p2, p2 = tmp;
            /*==============================illum===============================================*/
            for (int ix = 0; ix < nx; ix++)
            {
                for (int iz = 0; iz < nz; iz++)
                {
                    illum[ix * nz + iz] += pow(p1[(ix + lft) * nzb + iz + top], 2);
                }
            }

            /*==================================================================================*/
            for (int it = 1; it < nt; it++)
            {
                if (verb)
                {
                    warn("forward start iter=%d/%d,is=%d/%d,it=%d/%d", iter + 1, nter, is + 1, ns, it, nt);
                }
                fun_forward(par, p0, p1, p2, lap + (it - 1) * nz * nx);
                tmp = p0, p0 = p1, p1 = p2, p2 = tmp;
                fun_addsrc(false, par, wt[it], p1, is);
                fun_record(par, p1, dcal + is * nr * nt + it * nr);
                for (int ig = 0; ig < nr; ig++)
                {
                    obj += pow(dcal[is * nr * nt + it * nr + ig] - dobs[is * nr * nt + it * nr + ig], 2);
                }
                for (int ix = 0; ix < nx; ix++)
                {
                    for (int iz = 0; iz < nz; iz++)
                    {
                        illum[ix * nz + iz] += pow(p1[(ix + lft) * nzb + iz + top], 2);
                    }
                }
            }

            for (int ix = 0; ix < nx; ix++)
            {
                for (int iz = 0; iz < nz; iz++)
                {
                    lap[(nt - 1) * nz * nx + ix * nz + iz] = fun_laplace(nzb, nxb, ix + lft, iz + top, p1, dz, dx);
                }
            }
            /* ==========================================================================================*/
            /* ===========================================backward imaging ===============================================*/
            memset(p0, 0, sizeof(float) * nzxb);
            memset(p1, 0, sizeof(float) * nzxb);
            memset(p2, 0, sizeof(float) * nzxb);
            for (int it = nt - 1; it >= 0; it--)
            {

                if (verb)
                {
                    warn("backward start iter=%d/%d,is=%d/%d,it=%d/%d", iter + 1, nter, is + 1, ns, it, nt);
                }
                fun_backward(par, p0, p1, p2);
                tmp = p0, p0 = p1, p1 = p2, p2 = tmp;
                for (int ir = 0; ir < nr; ir++)
                {
                    fun_addsrc(true, par, dcal[is * nr * nt + it * nr + ir] - dobs[is * nr * nt + it * nr + ir], p1, ir);
                }
                for (int ix = 0; ix < nx; ix++)
                {
                    for (int iz = cut; iz < nz; iz++)
                    {
                        grad[ix * nz + iz] += 2 * lap[it * nz * nx + ix * nz + iz] * p1[(ix + lft) * nzb + iz + top] / vel[ix * nz + iz] / (illum[ix * nz + iz] + EPS);
                    }
                }
            }
        }
        eal_close();
        /* ==========================================================================================*/

        /*====================================== test model ===========================================*/
        /* test step length */
        float alpha_test = fun_alphatest(par, grad);
        /* forward modeling for tmp model */
        for (int ix = 0; ix < nzx; ix++)
        {
            partmp->v[ix] = par->v[ix] - alpha_test * grad[ix];
        }
        pad2(partmp->v, partmp->vv, nz, nx, lft, rht, top, bot);
        eal_init(partmp, 1e-3, 0);
        num1 = 0;
        num2 = 0;
        for (int is = 0; is < ns; is++)
        {
            memset(p0, 0, sizeof(float) * nzxb);
            memset(p1, 0, sizeof(float) * nzxb);
            memset(p2, 0, sizeof(float) * nzxb);
            for (int it = 0; it < nt; it++)
            {
                fun_forward(partmp, p0, p1, p2, NULL);
                tmp = p0, p0 = p1, p1 = p2, p2 = tmp;
                fun_addsrc(false, partmp, wt[it], p1, is);
                fun_record(partmp, p1, dcaltmp);
                for (int ir = 0; ir < nr; ir++)
                {
                    num1 += (dcaltmp[ir] - dcal[is * nr * nt + it * nr + ir]) * (dcal[is * nr * nt + it * nr + ir] - dobs[is * nr * nt + it * nr + ir]);
                    num2 += dcaltmp[ir] * dcaltmp[ir];
                }
            }
        }
        eal_close();
        /*==============================================================================================================*/
        alpha = -alpha_test * num1 / num2;
        for (int ix = 0; ix < nzx; ix++)
        {
            par->v[ix] -= alpha * grad[ix];
        }
        pad2(par->v, par->vv, nz, nx, lft, rht, top, bot);
        if (verb)
        {
            finish = clock();
            warn("costtime=%.4f,iter=%d,obj=%g,alpha=%g", 1. * (finish - start) / CLOCKS_PER_SEC, iter + 1, obj, alpha);
        }
        fprintf(Fobj, "%g\n", obj);
        fflush(Fobj);
    }
    {
        FILE *fd = fopen(out, "wb");
        if (fd == NULL)
            err("can't open output file");
        fwrite(par->v, sizeof(float) * nzx, 1, fd);
        fclose(fd);
    }
    fclose(Fobj);
    return 0;
}