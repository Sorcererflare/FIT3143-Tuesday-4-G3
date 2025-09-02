// By: Abdullah (33861641) (16.7%), Hayden (33861889) (16.7%), Samuel (33114110) (16.7%), 
// Hesamreza (33861544) (16.7%), Peter (33143722) (16.7%), Khang (33048258) (16.7%)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char* argv[]) {
    int rank, size;
    long long N;
    double local_sum = 0.0;
    double total_sum = 0.0;
    double piVal;
    double start_time, end_time;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        printf("Enter the value of N: ");
        scanf("%lld", &N);
    }

    MPI_Bcast(&N, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    start_time = MPI_Wtime();

    long long chunk_size = N / size;
    long long start_index = rank * chunk_size;
    long long end_index = (rank == size - 1) ? N : start_index + chunk_size;

    for (long long i = start_index; i < end_index; i++) {
        double x = (2.0 * i + 1.0) / (2.0 * N);
        local_sum += 4.0 / (1.0 + x * x);
    }

    MPI_Reduce(&local_sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        piVal = total_sum / N;
        end_time = MPI_Wtime();
        printf("Approximated Pi value: %.15f\n", piVal);
        printf("Time taken: %.6f seconds\n", end_time - start_time);
    }

    MPI_Finalize();
    return 0;
}
