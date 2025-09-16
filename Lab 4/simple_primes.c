#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define MAX_PRIMES 500000
//By: Abdullah (33861641) (16.7%), Hayden (33861889) (16.7%), Samuel (33114110) (16.7%),  Hesamreza (33861544) (16.7%), Peter (33143722) (16.7%), Khang (33048258) (16.7%)

int primality_check(int n) {
    if (n == 2) return 1;
    if (n == 0 || n == 1 || n % 2 == 0) return 0;
    int square_root = (int)sqrt(n);
    for (int i = 3; i <= square_root; i += 2) {
        if (n % i == 0) return 0;
    }
    return 1;
}

int compare(const void *a, const void *b) {
    return (*(int*)a - *(int*)b);
}

int main(int argc, char *argv[]) {
    int rank, size, n;
    double start_time, end_time, compute_start, compute_end;
    double comm_time = 0.0, compute_time = 0.0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    start_time = MPI_Wtime();


    if (rank == 0) {
        if (argc != 2) {
            printf("Usage: mpirun -np <processes> %s <n>\n", argv[0]);
            return 1;
        }
        n = atoi(argv[1]);
        if (n <= 2) {
            printf("Error: n must be greater than 2\n");
            return 1;
        }
    }


    double comm_start = MPI_Wtime();
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    comm_time += MPI_Wtime() - comm_start;


    int local_primes[MAX_PRIMES];
    int local_count = 0;


    compute_start = MPI_Wtime();


    if (rank == 0 && 2 < n) {
        local_primes[local_count] = 2;
        local_count++;
    }


    int start_odd = 3 + 2 * rank;
    for (int num = start_odd; num < n; num += 2 * size) {
        if (primality_check(num)) {
            if (local_count < MAX_PRIMES) {
                local_primes[local_count] = num;
                local_count++;
            }
        }
    }

    compute_end = MPI_Wtime();
    compute_time = compute_end - compute_start;


    qsort(local_primes, local_count, sizeof(int), compare);


    comm_start = MPI_Wtime();
    int all_counts[size];

    MPI_Gather(&local_count, 1, MPI_INT, all_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);


    int total_primes = 0;
    int displacements[size];
    if (rank == 0) {
        for (int i = 0; i < size; i++) {
            displacements[i] = total_primes;
            total_primes += all_counts[i];
        }
    }


    int all_primes[total_primes > 0 ? total_primes : 1];

    MPI_Gatherv(local_primes, local_count, MPI_INT,
                all_primes, all_counts, displacements, MPI_INT,
                0, MPI_COMM_WORLD);

    comm_time += MPI_Wtime() - comm_start;


    if (rank == 0) {

        qsort(all_primes, total_primes, sizeof(int), compare);


        FILE *file = fopen("primes.txt", "w");
        if (file != NULL) {
            for (int i = 0; i < total_primes; i++) {
                fprintf(file, "%d ", all_primes[i]);
            }
            fclose(file);
        }
    }

    end_time = MPI_Wtime();


    if (rank == 0) {
        double total_time = end_time - start_time;
        printf("Total primes found: %d\n", total_primes);
        printf("Total execution time: %.4f seconds\n", total_time);
    }

    MPI_Finalize();
    return 0;
}
