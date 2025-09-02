// BCast is a lot simpler than Send, can send/receive through BCast and can 
// send to all processes in a single statement

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
            if (value < 0) break;
        }
        
        // Broadcast value from rank 0 to ALL processes (including rank 0 itself)
        MPI_Bcast(&value, 1, MPI_INT, 0, MPI_COMM_WORLD);

        

        // Print rank and received value
        printf("Process %d: Sent %d to all processes.\n", rank, value);
        fflush(stdout);
    }

    MPI_Finalize(); // Shut down MPI
    return 0;
}