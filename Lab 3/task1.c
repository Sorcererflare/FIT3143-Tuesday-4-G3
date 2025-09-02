#include <stdio.h>
#include <mpi.h>


int main(int argc, char* argv[]){

    // *Vars
    int numTasks, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf("Hello World from process %d of %d\n",
    rank, numTasks);

    MPI_Finalize();

    return 0;

}
