#include <mpi.h>
#include <stdio.h>
#include <unistd.h>

int main(int argc, char *argv[])
{
    int rank, size, len;
    char hostname[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Get_processor_name(hostname, &len);

    printf("Hostname: %s, Rank: %d, Size: %d, CPU cores: %ld\n", hostname, rank, size, sysconf(_SC_NPROCESSORS_ONLN));

    MPI_Finalize();
    return 0;
}
