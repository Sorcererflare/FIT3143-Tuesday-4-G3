#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#define NUM_THREADS 8
#define THREADS_PER_PROCESS 2

void transpose(int *pInMatrix, int rows, int cols, int *pTransposedMatrix);

int write_matrix_to_file(int* matrix, char* filename, int row, int col){
    FILE *fp = fopen(filename, "wb");
    if (fp == NULL){
        return 0;
    }
    fwrite(&row, sizeof(int), 1, fp); 
    fwrite(&col, sizeof(int), 1, fp);
    int i;

    for (i = 0; i < row; i++) {
        fwrite(&matrix[i * col], sizeof(int), col, fp);
    }
    fclose(fp);
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
        int tileRows, tileCols;
        
        //Assume a square number of worker processes
        tileRows = (int)sqrt(numTiles);
        tileCols = tileRows;

        if (argc != 4) {
        printf("Usage: %s <MatrixA_File> <MatrixB_File> <MatrixC_Output_File>\n", argv[0]);
        printf("Example: %s MA_6x7.bin MB_7x6.bin MC_6x6.bin\n", argv[0]);
        return 0;
        }

        char* filenameA = argv[1];
        char* filenameB = argv[2];
        char* filenameC = argv[3];

        //Reading matrices in parallel
        int rowA, colA, rowB, colB;

        int i, j, k;

        int fdA = open(filenameA, O_RDONLY);
        if (fdA < 0){
            printf("error opening %s\n",filenameA);
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 0;
        }

        read(fdA, &rowA, sizeof(int));
        read(fdA, &colA, sizeof(int));

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

        //Test reading by writing file back
        //write_matrix_to_file(A, "A_writeback.bin", rowA, colA);

        printf("Succesfully Read Matrix A\n");

        int fdB = open(filenameB, O_RDONLY);
        if (fdB < 0){
            printf("error opening %s\n",filenameA);
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 0;
        }

        read(fdA, &rowB, sizeof(int));
        read(fdA, &colB, sizeof(int));

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

        int rowC = rowA, colC = colB;
        
        int rowsPerProcess = (int) rowC/tileRows;
        int colsPerProcess = (int) colC/tileCols;

        int tileRowNum, tileColNum;

        for (i = 1; i<size; i++){
            tileRowNum = (i-1)/(tileCols);
            tileColNum = (i-1)%tileCols;
            printf("Process %d will take on tile (%d, %d)\n", i, tileRowNum, tileColNum);
        }
        MPI_Finalize();
        return 0;
    }else{
        //Worker Processes go here
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