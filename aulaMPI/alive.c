#include <mpi.h>
#include <stdio.h>
#include <string.h>

#define STD_TAG 0
int main(int argc, char *argv[])
{
    int i, my_rank, n_procs; char msg[100]; MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    if(my_rank != 0){
        sprintf(msg,"I'm alive!");
        // Proc send message
        MPI_Send(msg, strlen(msg) + 1, MPI_CHAR, 0, STD_TAG, MPI_COMM_WORLD);
        // Proc wait for response
        MPI_Recv(msg, 100, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        printf("Proc %d: %s \n", status.MPI_SOURCE, msg);
    }else{
        for (i = 1; i < n_procs; i++)
        {
            // wait for msg from i
            MPI_Recv(msg, 100, MPI_CHAR, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            printf("Proc %d: %s \n", status.MPI_SOURCE, msg);
            sprintf(msg,"Thank you proc %d",i);
            // send response
            MPI_Send(msg, strlen(msg) + 1, MPI_CHAR, i, STD_TAG, MPI_COMM_WORLD);
        }
        
    }
    MPI_Finalize();
}
