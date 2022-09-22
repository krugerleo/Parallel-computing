#define NPES 4
#define NC 5000       /* Number of Cols        */
#define NR 5000       /* Number of Rows        */
#define NRL NR / NPES /* Number of Rows per PE */
#define DOWN 100      /* Tag for messages down */
#define UP 101        /* Tag for messages up   */
#define ROOT 0        /* The root PE           */
#define MAX(x, y) (((x) > (y)) ? x : y)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

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

    clock_t begin = clock();

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

    niter = 1000;

    for (iter = 1; iter <= niter; iter++)
    {

        for (i = 1; i <= NRL; i++)
            for (j = 1; j <= NC; j++)
                t[i][j] = 0.25 * (told[i + 1][j] + told[i - 1][j] +
                                  told[i][j + 1] + told[i][j - 1]);
        dt = 0.;

        for (i = 1; i <= NRL; i++) /* Copy for next iteration  */
            for (j = 1; j <= NC; j++)
            {
                dt = MAX(abs(t[i][j] - told[i][j]), dt);
                told[i][j] = t[i][j];
            }

        /*------------------------*/
        /* Print some test values */
        /*------------------------*/

        if ((iter % 100) == 0)
        {
            if (mype == 0)
                printf("Iter = %4d: PE = %d: t[10][10] = %20.8f\n",
                       iter, mype, t[10][10]);
        }

    } /* End of iteration */
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("%f", time_spent);

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
