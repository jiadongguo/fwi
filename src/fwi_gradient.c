#include "cstd.h"
#include "waveutils.h"
#include <mpi/mpi.h>
int main(int argc, char **argv)
{
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    initargs(argc, argv);
    int nz, nx, nt, top, bot, lft, rht;
    int ns, sz, sx, jsx, jsz, rz, rx, jrx, jrz, nr;
    float dz, dx, dt;
    char *fwt, *fvel, *out;
    float *wt, *vel;
    MPI_Finalize();
    return 0;
}