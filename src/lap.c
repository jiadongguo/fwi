/* finite difference laplace */
#include <math.h>
float laplace(int n1, int n2, int i1, int i2, float *curr, float d1, float d2)
{
    static float c0 = -5. / 2, c1 = 4. / 3, c2 = -1. / 12;
    float fdx[5], fdz[5];
    float lapz, lapx;
    d1 = 1./pow(d1, 2);
    d2 = 1./pow(d2, 2);
    for (int i = 0; i < 5; i++)
    {
        fdx[i] = ((i2 + i - 2 < 0) || (i2 + i - 2 >= n2)) ? 0 : curr[(i2 + i - 2) * n1 + i1];
        fdz[i] = ((i1 + i - 2 < 0) || (i1 + i - 2 >= n1)) ? 0 : curr[i2 * n1 + i1 + i - 2];
    }
    lapz = (c2 * (fdz[0] + fdz[4]) + c1 * (fdz[1] + fdz[3]) + c0 * fdz[2]) * d1;
    lapx = (c2 * (fdx[0] + fdx[4]) + c1 * (fdx[1] + fdx[3]) + c0 * fdx[2]) * d2;
    return lapx + lapz;
}