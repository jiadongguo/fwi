#include "cstd.h"
#include "waveutils.h"
#include <mpi/mpi.h>
int main(int argc, char **argv)
{
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    initargs(argc, argv);
    MPI_Finalize();
    return 0;
}