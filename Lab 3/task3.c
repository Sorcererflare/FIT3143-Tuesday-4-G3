/*
FIT3143 Lab3 Task 3 Group 3
// By: Abdullah (33861641) (16.7%), Hayden (33861889) (16.7%), Samuel (33114110) (16.7%), 
// Hesamreza (33861544) (16.7%), Peter (33143722) (16.7%), Khang (33048258) (16.7%)

*/


#include <stdio.h>
#include <mpi.h>
struct valuestruct {
    int a;
    double b;
};


int main(int argc, char** argv){

    //Declare Parameters
    int source = 0; //Source Process to broadcast from
    int size = 1; // Number of items to broadcast
    struct valuestruct values;
    int myrank;
    MPI_Datatype Valuetype; // Declare derived datatype
    MPI_Datatype type[2] = { MPI_INT, MPI_DOUBLE }; // Sub-types
    int blocklen[2] = { 1, 1}; //Number of each sub-type
    MPI_Aint disp[2]; //Displacements to each block
    MPI_Init(&argc, &argv);

    //Initialise MPI environment
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Get_address(&values.a, &disp[0]);
    MPI_Get_address(&values.b, &disp[1]);
    //Make relative
    disp[1]=disp[1]-disp[0];
    disp[0]=0;

    // Create MPI struct
    MPI_Type_create_struct(2, blocklen, disp, type, &Valuetype);
    MPI_Type_commit(&Valuetype);
    
    do{
        if (myrank == 0){
            printf("Enter an round number (>0) & a real number: ");
            fflush(stdout);
            scanf("%d%lf", &values.a, &values.b);
        }
        //Broadcast item
        int errCode;
        errCode = MPI_Bcast(&values, size, Valuetype, source, MPI_COMM_WORLD);
        if (errCode != MPI_SUCCESS){
            printf("Rank: %d, encountered error: %d", myrank, errCode);
            fflush(stdout);
            MPI_Abort(MPI_COMM_WORLD, errCode);
        }
        //Print
        printf("Rank: %d. values.a = %d. values.b = %lf\n",
        myrank, values.a, values.b);
        fflush(stdout);
    }while(values.a > 0);
    /* Clean up the type */
    MPI_Type_free(&Valuetype);
    MPI_Finalize();
    return 0;
}
