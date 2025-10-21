#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <math.h>
#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define ROOT_RANK 0
#define NUM_READ_THREADS 8
#define OMP_THREADS_PER_PROCESS 2

static int write_matrix_to_file_ull(unsigned long long* matrix, char* filename, int row, int col){
    FILE *fp = fopen(filename, "wb");
    if (fp == NULL){
        return 0;
    }
    fwrite(&row, sizeof(int), 1, fp);
    fwrite(&col, sizeof(int), 1, fp);
    for (int i = 0; i < row; i++) {
        fwrite(&matrix[i * col], sizeof(unsigned long long), col, fp);
    }
    fclose(fp);
    return 1;
}

static int is_perfect_square(int p){
    int q = (int)floor(sqrt((double)p));
    return q*q == p;
}

static void extract_tile_int(const int *matrix, int rows, int cols,
                             int startRow, int startCol, int blockSize,
                             int *tileOut){
    for (int r = 0; r < blockSize; r++){
        const int *src = matrix + ((startRow + r) * cols) + startCol;
        int *dst = tileOut + r * blockSize;
        for (int c = 0; c < blockSize; c++){
            dst[c] = src[c];
        }
    }
}

static void place_tile_ull(unsigned long long *matrix, int rows, int cols,
                           int startRow, int startCol, int blockSize,
                           const unsigned long long *tile){
    for (int r = 0; r < blockSize; r++){
        unsigned long long *dst = matrix + ((startRow + r) * cols) + startCol;
        const unsigned long long *src = tile + r * blockSize;
        for (int c = 0; c < blockSize; c++){
            dst[c] = src[c];
        }
    }
}

int main(int argc, char *argv[]){
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 4){
        if (rank == ROOT_RANK){
            printf("Usage: %s <MatrixA_File> <MatrixB_File> <MatrixC_Output_File>\n", argv[0]);
            printf("Example: %s MA_1024x1024.bin MB_1024x1024.bin MC_1024x1024.bin\n", argv[0]);
        }
        MPI_Finalize();
        return 0;
    }

    char *filenameA = argv[1];
    char *filenameB = argv[2];
    char *filenameC = argv[3];

    // Basic validation for Fox: size must be a perfect square, and matrices square with N % q == 0
    if (!is_perfect_square(size)){
        if (rank == ROOT_RANK){
            printf("Error: Number of MPI processes (%d) must be a perfect square for Fox.\n", size);
        }
        MPI_Abort(MPI_COMM_WORLD, 20);
        return 0;
    }
    int q = (int)round(sqrt((double)size));

    // Timers
    double t_total_start = MPI_Wtime();
    double t_read_start = 0.0, t_read_end = 0.0, t_dist_start = 0.0, t_dist_end = 0.0;
    double t_compute_start = 0.0, t_compute_end = 0.0, t_gather_start = 0.0, t_gather_end = 0.0;

    int N = 0; // dimension (square)
    int rowA=0,colA=0,rowB=0,colB=0;
    int *A_full = NULL;
    int *B_full = NULL;

    if (rank == ROOT_RANK){
        t_read_start = MPI_Wtime();
        int fdA = open(filenameA, O_RDONLY);
        if (fdA < 0){
            printf("Error opening %s\n", filenameA);
            MPI_Abort(MPI_COMM_WORLD, 21);
            return 0;
        }
        ssize_t nread = read(fdA, &rowA, sizeof(int));
        if (nread != (ssize_t)sizeof(int)){
            printf("Error reading rowA header from %s\n", filenameA);
            close(fdA);
            MPI_Abort(MPI_COMM_WORLD, 22);
            return 0;
        }
        nread = read(fdA, &colA, sizeof(int));
        if (nread != (ssize_t)sizeof(int)){
            printf("Error reading colA header from %s\n", filenameA);
            close(fdA);
            MPI_Abort(MPI_COMM_WORLD, 23);
            return 0;
        }
        int elemsA = rowA * colA;
        int offset = 2 * (int)sizeof(int);
        int chunkA = (NUM_READ_THREADS > 0) ? (elemsA / NUM_READ_THREADS) : elemsA;
        A_full = (int*)malloc((size_t)elemsA * sizeof(int));
        if (A_full == NULL){
            printf("Error: malloc A_full\n");
            close(fdA);
            MPI_Abort(MPI_COMM_WORLD, 24);
            return 0;
        }
        #pragma omp parallel num_threads(NUM_READ_THREADS)
        {
            int tid = 0;
            #ifdef _OPENMP
            tid = omp_get_thread_num();
            #endif
            int start = tid * chunkA;
            int end = (tid == NUM_READ_THREADS - 1) ? elemsA : start + chunkA;
            int count = end - start;
            if (count > 0){
                int totalOffset = offset + start * (int)sizeof(int);
                int *dest = A_full + start;
                pread(fdA, dest, (size_t)count * sizeof(int), totalOffset);
            }
        }
        close(fdA);

        int fdB = open(filenameB, O_RDONLY);
        if (fdB < 0){
            printf("Error opening %s\n", filenameB);
            MPI_Abort(MPI_COMM_WORLD, 25);
            return 0;
        }
        nread = read(fdB, &rowB, sizeof(int));
        if (nread != (ssize_t)sizeof(int)){
            printf("Error reading rowB header from %s\n", filenameB);
            close(fdB);
            MPI_Abort(MPI_COMM_WORLD, 26);
            return 0;
        }
        nread = read(fdB, &colB, sizeof(int));
        if (nread != (ssize_t)sizeof(int)){
            printf("Error reading colB header from %s\n", filenameB);
            close(fdB);
            MPI_Abort(MPI_COMM_WORLD, 27);
            return 0;
        }
        int elemsB = rowB * colB;
        int chunkB = (NUM_READ_THREADS > 0) ? (elemsB / NUM_READ_THREADS) : elemsB;
        B_full = (int*)malloc((size_t)elemsB * sizeof(int));
        if (B_full == NULL){
            printf("Error: malloc B_full\n");
            close(fdB);
            MPI_Abort(MPI_COMM_WORLD, 28);
            return 0;
        }
        #pragma omp parallel num_threads(NUM_READ_THREADS)
        {
            int tid = 0;
            #ifdef _OPENMP
            tid = omp_get_thread_num();
            #endif
            int start = tid * chunkB;
            int end = (tid == NUM_READ_THREADS - 1) ? elemsB : start + chunkB;
            int count = end - start;
            if (count > 0){
                int totalOffset = offset + start * (int)sizeof(int);
                int *dest = B_full + start;
                pread(fdB, dest, (size_t)count * sizeof(int), totalOffset);
            }
        }
        close(fdB);
        t_read_end = MPI_Wtime();

        if (rowA != colA || rowB != colB || colA != rowB){
            printf("Error: Fox assumes square matrices with compatible dims (A NxN, B NxN).\n");
            free(A_full); free(B_full);
            MPI_Abort(MPI_COMM_WORLD, 29);
            return 0;
        }
        N = rowA;
        if (N % q != 0){
            printf("Error: Matrix size N=%d must be divisible by q=%d (sqrt(p)).\n", N, q);
            free(A_full); free(B_full);
            MPI_Abort(MPI_COMM_WORLD, 30);
            return 0;
        }
    }

    // Broadcast N to all ranks
    MPI_Bcast(&N, 1, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);
    if (N == 0){
        MPI_Finalize();
        return 0;
    }

    int blockSize = N / q;
    int myRow = rank / q;
    int myCol = rank % q;
    MPI_Comm rowComm, colComm;
    MPI_Comm_split(MPI_COMM_WORLD, myRow, myCol, &rowComm);
    MPI_Comm_split(MPI_COMM_WORLD, myCol, myRow, &colComm);

    // Allocate local tiles
    int *A_local = (int*)malloc((size_t)blockSize * (size_t)blockSize * sizeof(int));
    int *B_local = (int*)malloc((size_t)blockSize * (size_t)blockSize * sizeof(int));
    unsigned long long *C_local = (unsigned long long*)calloc((size_t)blockSize * (size_t)blockSize, sizeof(unsigned long long));
    if (A_local == NULL || B_local == NULL || C_local == NULL){
        if (A_local) free(A_local);
        if (B_local) free(B_local);
        if (C_local) free(C_local);
        if (rank == ROOT_RANK){ if (A_full) free(A_full); if (B_full) free(B_full); }
        MPI_Abort(MPI_COMM_WORLD, 31);
        return 0;
    }

    // Root distributes initial tiles Aij and Bij
    t_dist_start = MPI_Wtime();
    if (rank == ROOT_RANK){
        for (int pr = 0; pr < size; pr++){
            int r = pr / q;
            int c = pr % q;
            int startRow = r * blockSize;
            int startCol = c * blockSize;

            if (pr == ROOT_RANK){
                extract_tile_int(A_full, N, N, startRow, startCol, blockSize, A_local);
                extract_tile_int(B_full, N, N, startRow, startCol, blockSize, B_local);
            } else {
                int *tmpA = (int*)malloc((size_t)blockSize * (size_t)blockSize * sizeof(int));
                int *tmpB = (int*)malloc((size_t)blockSize * (size_t)blockSize * sizeof(int));
                if (!tmpA || !tmpB){
                    if (tmpA) free(tmpA);
                    if (tmpB) free(tmpB);
                    free(A_local); free(B_local); free(C_local);
                    free(A_full); free(B_full);
                    MPI_Abort(MPI_COMM_WORLD, 32);
                    return 0;
                }
                extract_tile_int(A_full, N, N, startRow, startCol, blockSize, tmpA);
                extract_tile_int(B_full, N, N, startRow, startCol, blockSize, tmpB);
                MPI_Send(tmpA, blockSize * blockSize, MPI_INT, pr, 100, MPI_COMM_WORLD);
                MPI_Send(tmpB, blockSize * blockSize, MPI_INT, pr, 101, MPI_COMM_WORLD);
                free(tmpA); free(tmpB);
            }
        }
        free(A_full); free(B_full);
    } else {
        MPI_Status st;
        MPI_Recv(A_local, blockSize * blockSize, MPI_INT, ROOT_RANK, 100, MPI_COMM_WORLD, &st);
        MPI_Recv(B_local, blockSize * blockSize, MPI_INT, ROOT_RANK, 101, MPI_COMM_WORLD, &st);
    }
    t_dist_end = MPI_Wtime();

    // Fox iterations
    t_compute_start = MPI_Wtime();
    int *A_bcast = (int*)malloc((size_t)blockSize * (size_t)blockSize * sizeof(int));
    if (A_bcast == NULL){
        free(A_local); free(B_local); free(C_local);
        MPI_Abort(MPI_COMM_WORLD, 33);
        return 0;
    }
    int bcastRoot;

    for (int k = 0; k < q; k++){
        bcastRoot = (myRow + k) % q;
        if (myCol == bcastRoot){
            for (int i = 0; i < blockSize * blockSize; i++) A_bcast[i] = A_local[i];
        }
        MPI_Bcast(A_bcast, blockSize * blockSize, MPI_INT, bcastRoot, rowComm);

        #ifdef _OPENMP
        omp_set_num_threads(OMP_THREADS_PER_PROCESS);
        #endif
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < blockSize; i++){
            for (int j = 0; j < blockSize; j++){
                unsigned long long sum = 0ULL;
                for (int t = 0; t < blockSize; t++){
                    sum += (unsigned long long)A_bcast[i * blockSize + t] * (unsigned long long)B_local[t * blockSize + j];
                }
                C_local[i * blockSize + j] += sum;
            }
        }

        // Shift B up by 1 along the column communicator (circular)
        int src = (myRow + 1) % q; // receive from below
        int dst = (myRow - 1 + q) % q; // send to above
        MPI_Status st;
        MPI_Sendrecv_replace(B_local, blockSize * blockSize, MPI_INT,
                             /*dest*/ dst, 200,
                             /*src*/ src, 200,
                             colComm, &st);
    }
    t_compute_end = MPI_Wtime();

    // Gather C tiles on root
    t_gather_start = MPI_Wtime();
    if (rank == ROOT_RANK){
        unsigned long long *C_full = (unsigned long long*)calloc((size_t)N * (size_t)N, sizeof(unsigned long long));
        if (!C_full){
            free(A_local); free(B_local); free(C_local); free(A_bcast);
            MPI_Abort(MPI_COMM_WORLD, 34);
            return 0;
        }
        for (int pr = 0; pr < size; pr++){
            int r = pr / q;
            int c = pr % q;
            int startRow = r * blockSize;
            int startCol = c * blockSize;
            if (pr == ROOT_RANK){
                place_tile_ull(C_full, N, N, startRow, startCol, blockSize, C_local);
            } else {
                unsigned long long *tmp = (unsigned long long*)malloc((size_t)blockSize * (size_t)blockSize * sizeof(unsigned long long));
                if (!tmp){
                    free(C_full); free(A_local); free(B_local); free(C_local); free(A_bcast);
                    MPI_Abort(MPI_COMM_WORLD, 35);
                    return 0;
                }
                MPI_Status st;
                MPI_Recv(tmp, blockSize * blockSize, MPI_UNSIGNED_LONG_LONG, pr, 300, MPI_COMM_WORLD, &st);
                place_tile_ull(C_full, N, N, startRow, startCol, blockSize, tmp);
                free(tmp);
            }
        }
        // Write output
        write_matrix_to_file_ull(C_full, filenameC, N, N);
        free(C_full);
    } else {
        MPI_Send(C_local, blockSize * blockSize, MPI_UNSIGNED_LONG_LONG, ROOT_RANK, 300, MPI_COMM_WORLD);
    }
    t_gather_end = MPI_Wtime();

    if (rank == ROOT_RANK){
        double t_total_end = MPI_Wtime();
        printf("Timing (s): read=%.6f distribute=%.6f compute=%.6f gather=%.6f total=%.6f\n",
               (t_read_end - t_read_start), (t_dist_end - t_dist_start),
               (t_compute_end - t_compute_start), (t_gather_end - t_gather_start),
               (t_total_end - t_total_start));
    }

    free(A_local); free(B_local); free(C_local); free(A_bcast);
    MPI_Comm_free(&rowComm);
    MPI_Comm_free(&colComm);
    MPI_Finalize();
    return 0;
}


