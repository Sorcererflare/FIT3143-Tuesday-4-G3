// BCast is a lot simpler than Send, can send/receive through BCast and can 
// send to all processes in a single statement

// By: Abdullah (33861641) (16.7%), Hayden (33861889) (16.7%), Samuel (33114110) (16.7%), 
// Hesamreza (33861544) (16.7%), Peter (33143722) (16.7%), Khang (33048258) (16.7%)

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
    int rank, size;
    int value = -1;

    MPI_Init(&argc, &argv);                 // Start MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   // Get process rank
    MPI_Comm_size(MPI_COMM_WORLD, &size);   // Get total number of processes

    while (1) {
        if (rank == 0) {
            // Root process reads input
            printf("Enter an integer (negative to quit): ");
            fflush(stdout);
            scanf("%d", &value);
        }
        
        // Broadcast value from rank 0 to ALL processes (including rank 0 itself)
        MPI_Bcast(&value, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (errCode != MPI_SUCCESS)
        {
            printf("Rank: %d, encountered error: %d", rank, errCode);
            fflush(stdout);
            MPI_Abort(MPI_COMM_WORLD, errCode);
        }

        if (value < 0) break;
        

        // Print rank and received value
        printf("Process %d: Sent %d to all processes.\n", rank, value);
        fflush(stdout);
    }

    MPI_Finalize(); // Shut down MPI
    return 0;
}
