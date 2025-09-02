#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

// By: Abdullah (33861641) (16.7%), Hayden (33861889) (16.7%), Samuel (33114110) (16.7%), 
// Hesamreza (33861544) (16.7%), Peter (33143722) (16.7%), Khang (33048258) (16.7%)

int main(int argc, char **argv) 
{
    int rank, size;
    int value;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) 
    {
        
        do 
        {
            
            printf("Process 0: Enter a non-negative integer: ");
            
            fflush(stdout);
            
            scanf("%d", &value);

            for (int i = 1; i < size; i++) 
            {
                MPI_Send(&value, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }

            printf("Process %d: Sent %d to all processes.\n", rank, value);
            
            fflush(stdout);

        } while (value >= 0);

    } else 
    {
        do 
        {
            
            MPI_Recv(&value, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            printf("Process %d: Received value %d.\n", rank, value);
            
            fflush(stdout);
            
        } while (value >= 0);
    }

    MPI_Finalize();
    return 0;
}
