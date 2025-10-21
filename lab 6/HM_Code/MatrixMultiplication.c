#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <math.h>
#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define NUM_THREADS 8
#define THREADS_PER_PROCESS 2

void transpose(int *pInMatrix, int rows, int cols, int *pTransposedMatrix);

int write_matrix_to_file_ull(unsigned long long* matrix, char* filename, int row, int col){
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

static void factor_workers_into_grid(int workers, int maxRows, int maxCols, int *outRows, int *outCols){
    if (workers <= 0){
        *outRows = 0; *outCols = 0; return;
    }
    int r = (int)floor(sqrt((double)workers));
    while (r > 1 && (workers % r != 0)) {
        r--;
    }
    int c = workers / r;
    if (r > maxRows) r = maxRows;
    if (c > maxCols) c = maxCols;
    if (r * c > workers){
        while (r * c > workers && c > 1) c--;
        if (r * c > workers && r > 1) r--;
    }
    if (r < 1) r = 1;
    if (c < 1) c = 1;
    *outRows = r;
    *outCols = c;
}

int main(int argc, char *argv[]) {
    int rank, size;
    
    //Initialise MPI Environment 
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //Determine number of tiles
    int numTiles = size - 1; //Leave root process as designator and controller

    if (rank == 0){
        //Controller Process
        int tileRows = 0, tileCols = 0;

        if (argc != 4) {
        printf("Usage: %s <MatrixA_File> <MatrixB_File> <MatrixC_Output_File>\n", argv[0]);
        printf("Example: %s MA_6x7.bin MB_7x6.bin MC_6x6.bin\n", argv[0]);
        return 0;
        }

        char* filenameA = argv[1];
        char* filenameB = argv[2];
        char* filenameC = argv[3];

        // Timers
        double t_total_start = MPI_Wtime();
        double t_read_start = 0.0, t_read_end = 0.0, t_dist_start = 0.0, t_dist_end = 0.0, t_gather_start = 0.0, t_gather_end = 0.0;

        //Reading matrices in parallel
        int rowA, colA, rowB, colB;

        int i, j, k;

        t_read_start = MPI_Wtime();
        int fdA = open(filenameA, O_RDONLY);
        if (fdA < 0){
            printf("error opening %s\n",filenameA);
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 0;
        }

        ssize_t nread = read(fdA, &rowA, sizeof(int));
        if (nread != (ssize_t)sizeof(int)){
            printf("Error reading header rowA from %s\n", filenameA);
            close(fdA);
            MPI_Abort(MPI_COMM_WORLD, 9);
            return 0;
        }
        nread = read(fdA, &colA, sizeof(int));
        if (nread != (ssize_t)sizeof(int)){
            printf("Error reading header colA from %s\n", filenameA);
            close(fdA);
            MPI_Abort(MPI_COMM_WORLD, 10);
            return 0;
        }

        int elemsA = rowA * colA;
        int offset = 2*sizeof(int);
        int chunkA = elemsA / NUM_THREADS;

        printf("%d\n", elemsA);

        //Allocate Memory to store matrix A
        int * A = (int*)malloc(elemsA*sizeof(int));

        #pragma omp parallel num_threads(NUM_THREADS)
        {
            int tid = omp_get_thread_num();
            int start = tid*chunkA;
            int end = (tid==NUM_THREADS - 1) ? elemsA : start+chunkA;
            int count = end-start;

            int totalOffset = offset + start*sizeof(int);
            int* dest = (int*) A + start;

            pread(fdA, dest, count * sizeof(int), totalOffset);
        }

        close(fdA);

        printf("Succesfully Read Matrix A\n");

        int fdB = open(filenameB, O_RDONLY);
        if (fdB < 0){
            printf("error opening %s\n",filenameB);
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 0;
        }

        nread = read(fdB, &rowB, sizeof(int));
        if (nread != (ssize_t)sizeof(int)){
            printf("Error reading header rowB from %s\n", filenameB);
            close(fdB);
            MPI_Abort(MPI_COMM_WORLD, 11);
            return 0;
        }
        nread = read(fdB, &colB, sizeof(int));
        if (nread != (ssize_t)sizeof(int)){
            printf("Error reading header colB from %s\n", filenameB);
            close(fdB);
            MPI_Abort(MPI_COMM_WORLD, 12);
            return 0;
        }

        int elemsB = rowB * colB;
        int chunkB = elemsB / NUM_THREADS;

        //Allocate Memory to store matrix B
        int * B = (int*)malloc(elemsB*sizeof(int));

        #pragma omp parallel num_threads(NUM_THREADS)
        {
            int tid = omp_get_thread_num();
            int start = tid*chunkB;
            int end = (tid==NUM_THREADS - 1) ? elemsB : start+chunkB;
            int count = end-start;

            int totalOffset = offset + start*sizeof(int);
            int* dest = (int*) B + start;

            pread(fdB, dest, count * sizeof(int), totalOffset);
        }

        close(fdB);

        printf("Succesfully Read Matrix B\n");
        t_read_end = MPI_Wtime();

        if (colA != rowB){
            printf("Matrix dimension mismatch: A(%d x %d) cannot multiply B(%d x %d)\n", rowA, colA, rowB, colB);
            free(A);
            free(B);
            MPI_Abort(MPI_COMM_WORLD, 2);
            return 0;
        }

        // Transpose B for contiguous column access during distribution/compute
        int *BT = (int*)malloc((rowB * colB) * sizeof(int));
        if (BT == NULL){
            printf("Error: Memory allocation failed for transposed B.\n");
            free(A);
            free(B);
            MPI_Abort(MPI_COMM_WORLD, 3);
            return 0;
        }
        transpose(B, rowB, colB, BT); // BT has shape (colB x rowB)

        int rowC = rowA, colC = colB;

        // Decide tile grid based on number of workers and matrix dims
        factor_workers_into_grid(numTiles, rowC, colC, &tileRows, &tileCols);
        int activeWorkers = tileRows * tileCols;
        if (activeWorkers == 0){
            // Single-process fallback (no workers)
            unsigned long long *C = (unsigned long long*)calloc((size_t)rowC * (size_t)colC, sizeof(unsigned long long));
            if (C == NULL){
                printf("Error: Memory allocation failed for C.\n");
                free(A); free(B); free(BT);
                MPI_Abort(MPI_COMM_WORLD, 4);
                return 0;
            }
            omp_set_num_threads(NUM_THREADS);
            #pragma omp parallel for schedule(static)
            for (int r = 0; r < rowC; r++){
                for (int c = 0; c < colC; c++){
                    unsigned long long sum = 0ULL;
                    for (int k2 = 0; k2 < colA; k2++){
                        sum += (unsigned long long)A[r * colA + k2] * (unsigned long long)BT[c * rowB + k2];
                    }
                    C[r * colC + c] = sum;
                }
            }
            write_matrix_to_file_ull(C, filenameC, rowC, colC);
            free(A); free(B); free(BT); free(C);
            MPI_Finalize();
            return 0;
        }

        // Precompute row and column partitions
        int baseRows = rowC / tileRows;
        int remRows = rowC % tileRows;
        int baseCols = colC / tileCols;
        int remCols = colC % tileCols;

        // Send tasks to all workers (even those beyond activeWorkers get zero sizes)
        t_dist_start = MPI_Wtime();
        for (int rankId = 1; rankId < size; rankId++){
            int tileIndex = rankId - 1;
            int tr = (tileIndex / tileCols);
            int tc = (tileIndex % tileCols);
            int rowCount = 0, colCount = 0, rowStart = 0, colStart = 0;
            if (tileIndex < activeWorkers){
                rowCount = baseRows + (tr < remRows ? 1 : 0);
                colCount = baseCols + (tc < remCols ? 1 : 0);
                rowStart = tr * baseRows + (tr < remRows ? tr : remRows);
                colStart = tc * baseCols + (tc < remCols ? tc : remCols);
            }

            int header[5];
            header[0] = rowStart;
            header[1] = rowCount;
            header[2] = colStart;
            header[3] = colCount;
            header[4] = colA; // shared inner dimension
            MPI_Send(header, 5, MPI_INT, rankId, 0, MPI_COMM_WORLD);

            if (rowCount > 0 && colCount > 0){
                // Send A sub-block (rowCount x colA), contiguous
                int *A_sub = A + (rowStart * colA);
                MPI_Send(A_sub, rowCount * colA, MPI_INT, rankId, 1, MPI_COMM_WORLD);

                // Send BT sub-block (colCount x colA), contiguous rows from BT
                int *BT_sub = BT + (colStart * rowB);
                MPI_Send(BT_sub, colCount * colA, MPI_INT, rankId, 2, MPI_COMM_WORLD);
            }
        }

        t_dist_end = MPI_Wtime();

        // Prepare C and receive results
        unsigned long long *C = (unsigned long long*)calloc((size_t)rowC * (size_t)colC, sizeof(unsigned long long));
        if (C == NULL){
            printf("Error: Memory allocation failed for C.\n");
            free(A); free(B); free(BT);
            MPI_Abort(MPI_COMM_WORLD, 5);
            return 0;
        }

        t_gather_start = MPI_Wtime();
        for (int rankId = 1; rankId < size; rankId++){
            int header[4];
            MPI_Status st;
            // Receive position and size info back: rowStart, rowCount, colStart, colCount
            MPI_Recv(header, 4, MPI_INT, rankId, 3, MPI_COMM_WORLD, &st);
            int rowStart = header[0];
            int rowCount = header[1];
            int colStart = header[2];
            int colCount = header[3];
            if (rowCount == 0 || colCount == 0){
                continue;
            }
            // Receive the tile data
            unsigned long long *tileBuf = (unsigned long long*)malloc((size_t)rowCount * (size_t)colCount * sizeof(unsigned long long));
            if (tileBuf == NULL){
                printf("Error: Memory allocation failed for tile buffer.\n");
                free(A); free(B); free(BT); free(C);
                MPI_Abort(MPI_COMM_WORLD, 6);
                return 0;
            }
            MPI_Recv(tileBuf, rowCount * colCount, MPI_UNSIGNED_LONG_LONG, rankId, 4, MPI_COMM_WORLD, &st);
            // Place into C
            for (int rr = 0; rr < rowCount; rr++){
                unsigned long long *dst = C + ((rowStart + rr) * colC) + colStart;
                unsigned long long *src = tileBuf + (rr * colCount);
                for (int cc = 0; cc < colCount; cc++){
                    dst[cc] = src[cc];
                }
            }
            free(tileBuf);
        }

        t_gather_end = MPI_Wtime();
        // Write result
        write_matrix_to_file_ull(C, filenameC, rowC, colC);

        double t_total_end = MPI_Wtime();
        printf("Timing (s): read=%.6f, distribute=%.6f, gather=%.6f, total=%.6f\n",
               (t_read_end - t_read_start), (t_dist_end - t_dist_start), (t_gather_end - t_gather_start), (t_total_end - t_total_start));

        free(A);
        free(B);
        free(BT);
        free(C);
        MPI_Finalize();
        return 0;
    }else{
        //Worker Processes: receive header, then sub-blocks, compute with OMP, send back
        int header[5];
        MPI_Status st;
        MPI_Recv(header, 5, MPI_INT, 0, 0, MPI_COMM_WORLD, &st);
        int rowStart = header[0];
        int rowCount = header[1];
        int colStart = header[2];
        int colCount = header[3];
        int innerDim = header[4];

        if (rowCount <= 0 || colCount <= 0){
            int ack[4] = {0,0,0,0};
            MPI_Send(ack, 4, MPI_INT, 0, 3, MPI_COMM_WORLD);
            MPI_Finalize();
            return 0;
        }

        int *A_sub = (int*)malloc((size_t)rowCount * (size_t)innerDim * sizeof(int));
        int *BT_sub = (int*)malloc((size_t)colCount * (size_t)innerDim * sizeof(int));
        if (A_sub == NULL || BT_sub == NULL){
            if (A_sub) free(A_sub);
            if (BT_sub) free(BT_sub);
            MPI_Abort(MPI_COMM_WORLD, 7);
            return 0;
        }
        MPI_Recv(A_sub, rowCount * innerDim, MPI_INT, 0, 1, MPI_COMM_WORLD, &st);
        MPI_Recv(BT_sub, colCount * innerDim, MPI_INT, 0, 2, MPI_COMM_WORLD, &st);

        unsigned long long *C_sub = (unsigned long long*)calloc((size_t)rowCount * (size_t)colCount, sizeof(unsigned long long));
        if (C_sub == NULL){
            free(A_sub); free(BT_sub);
            MPI_Abort(MPI_COMM_WORLD, 8);
            return 0;
        }

        omp_set_num_threads(THREADS_PER_PROCESS);
        #pragma omp parallel for collapse(2) schedule(static)
        for (int r = 0; r < rowCount; r++){
            for (int c = 0; c < colCount; c++){
                unsigned long long sum = 0ULL;
                int *aRow = A_sub + (r * innerDim);
                int *btRow = BT_sub + (c * innerDim);
                for (int k = 0; k < innerDim; k++){
                    sum += (unsigned long long)aRow[k] * (unsigned long long)btRow[k];
                }
                C_sub[r * colCount + c] = sum;
            }
        }

        int retHeader[4] = {rowStart, rowCount, colStart, colCount};
        MPI_Send(retHeader, 4, MPI_INT, 0, 3, MPI_COMM_WORLD);
        MPI_Send(C_sub, rowCount * colCount, MPI_UNSIGNED_LONG_LONG, 0, 4, MPI_COMM_WORLD);

        free(A_sub);
        free(BT_sub);
        free(C_sub);
        MPI_Finalize();
        return 0;
    }
    
}

void transpose(int *pInMatrix, int rows, int cols, int *pTransposedMatrix)
{
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            pTransposedMatrix[(j * rows) + i] = pInMatrix[(i * cols) + j];
        }
    }
}