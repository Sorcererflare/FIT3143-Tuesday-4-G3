#include <stdio.h>
#include <mpi.h>
// By: Abdullah (33861641) (16.7%), Hayden (33861889) (16.7%), Samuel (33114110) (16.7%), 
// Hesamreza (33861544) (16.7%), Peter (33143722) (16.7%), Khang (33048258) (16.7%)

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
