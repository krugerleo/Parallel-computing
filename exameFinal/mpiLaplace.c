#define NPES 4
#define NC 5000       /* Number of Cols        sk*/
#define NR 5000       /* Number of Rows        */
#define NRL NR / NPES /* Number of Rows per PE */
#define DOWN 100      /* Tag for messages down */
#define UP 101        /* Tag for messages up   */
#define ROOT 0        /* The root PE           */
#define MAX(x, y) (((x) > (y)) ? x : y)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

void initialize(float **t);
void set_bcs(float **t, int mype, int npes);

int main(int argc, char **argv)
{

    int npes;  /* Number of PEs */
    int mype;  /* My PE number  */
    int stat;  /* Error Status  */
    int niter; /* iter counter  */

    // float      t[NRL+2][NC+2], told[NRL+2][NC+2];
    float **t, **told;
    float dt;  /* Delta t       */
    float dtg; /* Delta t global*/

    int i, j, iter;

    int my_rank;
    int num_procs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);   // grab this process's rank
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs); // grab the total num of processes

    t = calloc(NR, sizeof(float *));
    t[0] = calloc(NR * NC, sizeof(float));
    for (i = 0; i <= NRL + 1; i++) /* Copy the values into told */
        t[i] = &t[0][i * NC];

    told = calloc(NR, sizeof(float *));
    told[0] = calloc(NR * NC, sizeof(float));
    for (i = 0; i <= NRL + 1; i++) /* Copy the values into told */
        told[i] = &told[0][i * NC];

    initialize(t); /* Give initial guesss of 0. */

    set_bcs(t, mype, npes); /* Set the Boundary values   */

    for (i = 0; i <= NRL + 1; i++) /* Copy the values into told */
        for (j = 0; j <= NC + 1; j++)
            told[i][j] = t[i][j];

    /*-------------------------------------------------*/
    /* Do Computation on Sub-grid for Niter iterations */
    /*-------------------------------------------------*/

    niter = 1000; // executa 10x
    int chunk_size = NRL / num_procs;
    int start_index = chunk_size * my_rank + 1;
    int finish_index = chunk_size * (my_rank + 1);
    float matrixTold[chunk_size];

    for (iter = 1; iter <= niter; iter++)
    {

        MPI_Scatter(told[0], chunk_size, MPI_FLOAT, &matrixTold, chunk_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

        for (i = 1; i <= chunk_size; i++)
            for (j = 1; j <= NC; j++)
                t[i][j] = 0.25 * (told[i + 1][j] + told[i - 1][j] + told[i][j + 1] + told[i][j - 1]);

        dt = 0.;

        for (i = 1; i <= chunk_size; i++) /* Copy for next iteration  */
            for (j = 1; j <= NC; j++)
            {
                dt = MAX(abs(t[i][j] - told[i][j]), dt);
                told[i][j] = t[i][j];
            }

        MPI_Gather(told[0], chunk_size, MPI_FLOAT, &matrixTold, chunk_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
        /*------------------------*/
        /* Print some test values */
        /*------------------------*/
        printf("ITEER");
        if ((iter % 100) == 0)
        {
            if (mype == 0)
                printf("Iter = %4d: PE = %d: t[10][10] = %20.8f\n",
                       iter, mype, t[10][10]);
        }

    } /* End of iteration */
    MPI_Finalize();
} /* End of Program */

/*-----------------------------------------------------*/
/* Initialize all the values to 0. as a starting value */
/*-----------------------------------------------------*/

void initialize(float **t)
{

    int i, j, iter;

    for (i = 0; i <= NRL + 1; i++) /* Initialize */
        for (j = 0; j <= NC + 1; j++)
            t[i][j] = 0.0;
}

/*----------------------------------------------------------------*/
/* Set the values at the boundary.  Values at the boundary do not */
/* Change through out the execution of the program                */
/*----------------------------------------------------------------*/

void set_bcs(float **t, int mype, int npes)
{

    int i, j;

    for (i = 0; i <= NRL + 1; i++)
    { /* Set Left and Right boundary */
        t[i][0] = 100.0;
        t[i][NC + 1] = 100.0;
    }

    if (mype == 0) /* Top boundary */
        for (j = 0; j <= NC + 1; j++)
            t[0][j] = 100.0;

    if (mype == npes - 1) /* Bottom boundary */
        for (j = 0; j <= NC + 1; j++)
            t[NRL + 1][j] = 100.0;
}
