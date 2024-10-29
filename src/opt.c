#include "opt.h"
static float c1 = 1e-4, c2 = 0.9;
static bool verb;
void opt_init(bool verb_)
{
    verb = verb_;
}

float line_search(int n,    // dimension of x
                  float *x, // input vector x
                  float *g, // gradient of misfit function
                  float *d, // descent direction
                  float fk, // value of misfit function
                  opt_fg fg // subroutine to evaluation function and gradient
)
/* bisection line search  based on Wolfe condition */
{
    int j;
    float gxd, c1_gxd, c2_gxd, fcost, fxx, alpha1, alpha2;
    float *xk, alpha;
    static float infinity = 1e10;
    static int nls = 10;
    alpha = 1.;
    alpha1 = 0;
    alpha2 = infinity;

    xk = alloc1float(n);
    memcpy(xk, x, n * sizeof(float));

    gxd = sdot(n, g, d);
    c1_gxd = c1 * gxd;
    c2_gxd = c2 * gxd;
    for (int ils = 0; ils < nls; ils++)
    {
        for (j = 0; j < n; j++)
            x[j] = xk[j] + alpha * d[j];

        fcost = fg(x, g);

        gxd = sdot(n, g, d);
        fxx = fk + alpha * c1_gxd;

        /* check Wolfe condition for current step length, see Nocedal book, eqn 3.6 */
        /* sufﬁcient decrease condition */
        if (fcost > fxx)
        {
            if (verb)
                printf("Wolfe condition 1 fails: insufficient misfit decrease!\n");
            alpha2 = alpha;
            alpha = 0.5 * (alpha1 + alpha2); // shrink the search interval
        }
        /* check Wolfe condition for current step length, see Nocedal book, eqn 3.6 */
        /* curvature condition */
        else if (gxd < c2_gxd)
        {
            if (verb)
                printf("Wolfe condition 2 fails: stepsize is too small!\n");
            alpha1 = alpha;
            if (alpha2 < infinity)
                alpha = 0.5 * (alpha1 + alpha2); // shrink the search interval
            else
                alpha *= 5.; // extend the search interval if alpha2==0
        }
        else
        {
            free1float(xk);
            return alpha;
        }
        if (verb)
        {
            printf("#line search %d, alp1=%f alp2=%f alp=%f\n", ils, alpha1, alpha2, alpha);
            printf("------------------------------------------------\n");
        }
    }
    free1float(xk);
    return 1e-3;
}
/* gradient descent */
void descent(int n, float *g, float *d)
{
    sscale(n, -1, g, d);
}
void cg(int n,     // dimension of x
        float *g0, // previous gradient of misfit function
        float *g1, // current gradient of misfit function
        float *d0, // previous update direction
        float *d1, // current update direction
        int mode   // beta select mode
)
{

    for (int i = 0; i < n; i++)
    {
        d1[i] = -g1[i];
    }
    if (d0 == NULL || g0 == NULL)
    {
        return;
    }
    float beta;
    if (mode == 1) // FR
    {
        float gdg1, gdg0;
        gdg1 = sdot(n, g1, g1);
        gdg0 = sdot(n, g0, g0);
        beta = gdg1 / gdg0;
    }
    else if (mode == 2) // PR
    {
        float gdg1, gdg0, gdg;
        gdg1 = sdot(n, g1, g1);
        gdg0 = sdot(n, g0, g0);
        gdg = sdot(n, g0, g1);
        beta = (gdg1 - gdg) / gdg0;
        beta = MAX(beta, 0);
    }
    else if (mode == 3) // HS
    {
        float gdg1, pdg0, pdg1, gdg;
        gdg1 = sdot(n, g1, g1);
        gdg = sdot(n, g0, g1);
        pdg1 = sdot(n, g1, d0);
        pdg0 = sdot(n, g0, d0);
        beta = (gdg1 - gdg) / (pdg1 - pdg0);
    }
    else // FR-PR
    {
        float betaPR, betaFR;
        float gdg1, gdg0, gdg;
        gdg1 = sdot(n, g1, g1);
        gdg0 = sdot(n, g0, g0);
        gdg = sdot(n, g0, g1);
        betaPR = (gdg1 - gdg) / gdg0;
        betaFR = gdg1 / gdg0;
        if (betaPR < -betaFR)
        {
            beta = -betaFR;
        }
        else if (betaPR > betaFR)
        {
            beta = betaFR;
        }
        else
        {
            beta = betaPR;
        }
    }
    saxpy(n, beta, d0, d1);
}