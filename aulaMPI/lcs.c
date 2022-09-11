#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "mpi.h"

#ifndef max
#define max(a, b) (((a) > (b)) ? (a) : (b))
#endif

void print_matrix(short *x, int row, int col)
{
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            printf("%d ", x[(i * (col + 1)) + j]);
        }
        printf("\n");
    }
}
int get_index_of_character(char *str, char x, int len)
{
    for (int i = 0; i < len; i++)
    {
        if (str[i] == x)
        {
            return i;
        }
    }
    return -1; // not found the character x in str
}
void calc_P_matrix_v1(short *P, char *b, int len_b, char *c, int len_c, int myrank, int chunk_size)
{
    char receive_array_for_scatter_c[chunk_size];
    short receive_array_for_scatter_p[chunk_size * (len_b + 1)];
    // Scatter the char array chunks by sending each process a particular chunk
    MPI_Scatter(c, chunk_size, MPI_CHAR, &receive_array_for_scatter_c, chunk_size, MPI_CHAR, 0, MPI_COMM_WORLD);
    // Scatter the char array chunks by sending each process a particular chunk
    MPI_Scatter(P, chunk_size * (len_b + 1), MPI_SHORT, &receive_array_for_scatter_p, chunk_size * (len_b + 1), MPI_SHORT, 0, MPI_COMM_WORLD);
    // Broadcast the whole b  array to everybody
    MPI_Bcast(b, len_b, MPI_CHAR, 0, MPI_COMM_WORLD);

    for (int i = 0; i < chunk_size; i++)
    {
        for (int j = 2; j < len_b + 1; j++)
        {
            if (b[j - 2] == receive_array_for_scatter_c[i]) // j-2 as b we assume here that b has a empty character in the beginning
            {
                receive_array_for_scatter_p[(i * (len_b + 1)) + j] = j - 1;
            }
            else
            {
                receive_array_for_scatter_p[(i * (len_b + 1)) + j] = receive_array_for_scatter_p[(i * (len_b + 1)) + j - 1];
            }
        }
    }

    // now gather all the calculated values of P matrix in process 0
    MPI_Gather(receive_array_for_scatter_p, chunk_size * (len_b + 1), MPI_SHORT, P, chunk_size * (len_b + 1), MPI_SHORT, 0, MPI_COMM_WORLD);
}

int lcs_yang_v1(short **DP, short *P, char *A, char *B, char *C, int m, int n, int u, int myrank, int chunk_size)
{

    MPI_Bcast(P, (u * (n + 1)), MPI_SHORT, 0, MPI_COMM_WORLD);
    for (int i = 1; i < m + 1; i++)
    {

        int c_i = get_index_of_character(C, A[i - 1], u);
        short dp_i_receive[chunk_size];
        // Broadcast the  whole B  array to everybody
        MPI_Scatter(DP[i], chunk_size, MPI_SHORT, &dp_i_receive, chunk_size, MPI_SHORT, 0, MPI_COMM_WORLD);
        int start_id = (myrank * chunk_size);
        int end_id = (myrank * chunk_size) + chunk_size;
        for (int j = start_id; j < end_id; j++) // if myrank=0 then j=start_id+1 else j=start_id
        {
            if (j == start_id && myrank == 0)
                j = j + 1;
            if (A[i - 1] == B[j - 1])
            {
                dp_i_receive[j - start_id] = DP[i - 1][j - 1] + 1;
            }
            else if (P[(c_i * (n + 1)) + j] == 0)
            {
                dp_i_receive[j - start_id] = max(DP[i - 1][j], 0);
            }
            else
            {
                dp_i_receive[j - start_id] = max(DP[i - 1][j], DP[i - 1][P[(c_i * (n + 1)) + j] - 1] + 1);
            }
        }
        // now gather all the calculated values of P matrix in process 0
        MPI_Allgather(dp_i_receive, chunk_size, MPI_SHORT, DP[i], chunk_size, MPI_SHORT, MPI_COMM_WORLD);
    }

    return DP[m - 1][n - 1];
}
void calPMatrixMPI(short *pMatrix, char *seqB, int sizeB, char *uniqAB, int sizeUniqAB, int myRank, int chunkSize)
{
    char uniqABBufferToReceiveFromScatter[chunkSize];
    short pMatrixBufferToReceiveFromScatter[chunkSize * (sizeB + 1)];

    // Scatter the char array chunks by sending each process a particular chunk
    MPI_Scatter(uniqAB, chunkSize, MPI_CHAR, &uniqABBufferToReceiveFromScatter, chunkSize, MPI_CHAR, 0, MPI_COMM_WORLD);
    // Scatter the char array chunks by sending each process a particular chunk
    MPI_Scatter(pMatrix, chunkSize * (sizeB + 1), MPI_SHORT, &pMatrixBufferToReceiveFromScatter, chunkSize * (sizeB + 1), MPI_SHORT, 0, MPI_COMM_WORLD);
    // Broadcast the whole b  array to everybody
    MPI_Bcast(seqB, sizeB, MPI_CHAR, 0, MPI_COMM_WORLD);

    for (int i = 0; i < chunkSize; i++)
    {
        for (int j = 2; j < sizeB + 1; j++)
        {
            if (seqB[j - 2] == uniqABBufferToReceiveFromScatter[i]) // j-2 as b we assume here that b has a empty character in the beginning
            {
                pMatrixBufferToReceiveFromScatter[(i * (sizeB + 1)) + j] = j - 1;
            }
            else
            {
                pMatrixBufferToReceiveFromScatter[(i * (sizeB + 1)) + j] = pMatrixBufferToReceiveFromScatter[(i * (sizeB + 1)) + j - 1];
            }
        }
    }

    // now gather all the calculated values of P matrix in process 0
    MPI_Gather(pMatrixBufferToReceiveFromScatter, chunkSize * (sizeB + 1), MPI_SHORT, pMatrix, chunkSize * (sizeB + 1), MPI_SHORT, 0, MPI_COMM_WORLD);
}

int lcsMPI(short **DP, short *P, char *A, char *B, char *C, int m, int n, int u, int myrank, int chunk_size)
{

    MPI_Bcast(P, (u * (n + 1)), MPI_SHORT, 0, MPI_COMM_WORLD);
    for (int i = 1; i < m + 1; i++)
    {

        // Broadcast the c_i  array to everybody
        // MPI_Bcast(A_i, 1, MPI_CHAR, 0, MPI_COMM_WORLD);
        // Broadcast the  whole B  array to everybody

        // Scatter the char array A chunks by sending each process a particular chunk
        // MPI_Scatter(A, chunk_size, MPI_CHAR,&scatter_receive_a,chunk_size,MPI_CHAR, 0, MPI_COMM_WORLD);

        int c_i = get_index_of_character(C, A[i - 1], u);
        short dp_i_receive[chunk_size];
        // Broadcast the  whole B  array to everybody
        MPI_Scatter(DP[i], chunk_size, MPI_SHORT, &dp_i_receive, chunk_size, MPI_SHORT, 0, MPI_COMM_WORLD);
        int start_id = (myrank * chunk_size);
        int end_id = (myrank * chunk_size) + chunk_size;
        for (int j = start_id; j < end_id; j++) // if myrank=0 then j=start_id+1 else j=start_id
        {
            if (j == start_id && myrank == 0)
                j = j + 1;
            if (A[i - 1] == B[j - 1])
            {
                dp_i_receive[j - start_id] = DP[i - 1][j - 1] + 1;
            }
            else if (P[(c_i * (n + 1)) + j] == 0)
            {
                dp_i_receive[j - start_id] = max(DP[i - 1][j], 0);
            }
            else
            {
                dp_i_receive[j - start_id] = max(DP[i - 1][j], DP[i - 1][P[(c_i * (n + 1)) + j] - 1] + 1);
            }
        }
        // now gather all the calculated values of P matrix in process 0
        MPI_Allgather(dp_i_receive, chunk_size, MPI_SHORT, DP[i], chunk_size, MPI_SHORT, MPI_COMM_WORLD);
    }
    return DP[m][n];
}

char *read_seq(char *fname)
{
    // file pointer
    FILE *fseq = NULL;
    // sequence size
    long size = 0;
    // sequence pointer
    char *seq = NULL;
    // sequence index
    int i = 0;

    // open file
    fseq = fopen(fname, "rt");
    if (fseq == NULL)
    {
        printf("Error reading file %s\n", fname);
        exit(1);
    }

    // find out sequence size to allocate memory afterwards
    fseek(fseq, 0L, SEEK_END);
    size = ftell(fseq);
    rewind(fseq);

    // allocate memory (sequence)
    seq = (char *)calloc(size + 1, sizeof(char));
    if (seq == NULL)
    {
        printf("Erro allocating memory for sequence %s.\n", fname);
        exit(1);
    }

    // read sequence from file
    while (!feof(fseq))
    {
        seq[i] = fgetc(fseq);
        if ((seq[i] != '\n') && (seq[i] != EOF))
            i++;
    }
    // insert string terminator
    seq[i] = '\0';

    // close file
    fclose(fseq);

    // return sequence pointer
    return seq;
}

short **allocateScoreMatrix(int cols, int rows)
{
    // Allocate memory for LCS score matrix
    short **scoreMatrix = (short **)malloc((rows + 1) * sizeof(short *));
    for (int i = 0; i < (rows + 1); i++)
        scoreMatrix[i] = (short *)calloc((cols + 1), sizeof(short));
    return scoreMatrix;
}

void initScoreMatrix(short **scoreMatrix, int sizeA, int sizeB)
{
    int i, j;
    // Fill first line of LCS score matrix with zeroes
    for (j = 0; j < (sizeA + 1); j++)
        scoreMatrix[0][j] = 0;

    // Do the same for the first collumn
    for (i = 1; i < (sizeB + 1); i++)
        scoreMatrix[i][0] = 0;
}
void clearScoreMatrix(short **scoreMatrix, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            scoreMatrix[i][j] = 0;
        }
    }
}
int LCS(short **scoreMatrix, int sizeA, int sizeB, char *seqA, char *seqB)
{
    int i, j;
    for (i = 1; i < sizeB + 1; i++)
    {
        for (j = 1; j < sizeA + 1; j++)
        {
            if (seqA[j - 1] == seqB[i - 1])
            {
                /* if elements in both sequences match,
                 the corresponding score will be the score from
                 previous elements + 1*/
                scoreMatrix[i][j] = scoreMatrix[i - 1][j - 1] + 1;
            }
            else
            {
                /* else, pick the maximum value (score) from left and upper elements*/
                scoreMatrix[i][j] = max(scoreMatrix[i - 1][j], scoreMatrix[i][j - 1]);
            }
        }
    }
    return scoreMatrix[sizeB][sizeA];
}
void printMatrix(char *seqA, char *seqB, short **scoreMatrix, int sizeA, int sizeB)
{
    int i, j;

    // print header
    printf("Score Matrix:\n");
    printf("========================================\n");

    // print LCS score matrix allong with sequences

    printf("    ");
    printf("%5c   ", ' ');

    for (j = 0; j < sizeA; j++)
        printf("%5c   ", seqA[j]);
    printf("\n");
    for (i = 0; i < sizeB + 1; i++)
    {
        if (i == 0)
            printf("    ");
        else
            printf("%c   ", seqB[i - 1]);
        for (j = 0; j < sizeA + 1; j++)
        {
            printf("%5d   ", scoreMatrix[i][j]);
        }
        printf("\n");
    }
    printf("========================================\n");
}

void freeScoreMatrix(short **scoreMatrix, int sizeB)
{
    int i;
    for (i = 0; i < (sizeB + 1); i++)
        free(scoreMatrix[i]);
    free(scoreMatrix);
}

int main(int argc, char **argv)
{
    if (argc <= 3)
    {
        printf("Error: No input files specified! Please specify the input files, and run again!\n  Example: a.out fileA fileB alphabetAB");
        return 0;
    }
    // sequence pointers for both sequences

    int myRank;
    int numProcs;
    int pMatrixChunkSize, scoreMatrixChunkSize; // chunk size for P matrix and Score matrix

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);   // get process's rank
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs); // get the total num of processes

    seqA = read_seq(argv[1]);   // rows
    seqB = read_seq(argv[2]);   // colunms
    uniqAB = read_seq(argv[3]); // alphabet

    // find out sizes
    sizeA = strlen(seqA);
    sizeB = strlen(seqB);
    sizeUniqAB = strlen(uniqAB);

    // allocate LCS score matrix
    short **scoreMatrix = allocateScoreMatrix(sizeA, sizeB);
    // initialize LCS score matrix
    initScoreMatrix(scoreMatrix, sizeA, sizeB);

    double start_time = MPI_Wtime();

    // fill up the rest of the matrix and return final score (element locate at the last line and collumn)
    short score = LCS(scoreMatrix, sizeA, sizeB, seqA, seqB);

    double stop_time = MPI_Wtime();

    if (myRank == 0)
    {
        // Print Iput Sizes
        printf("(%d,%d,%d);", sizeA, sizeB, sizeUniqAB);
        // Print Serial LCS Score
        printf("%d;", score);
        // Print Serial LCS Time
        printf("%f;", stop_time - start_time);
    }

    clearScoreMatrix(scoreMatrix, sizeB + 1, sizeA + 1);

    // -----------------------------------------------------------------------------------------------------------------------------------

    short *pMatrix = (short *)malloc((sizeUniqAB * (sizeB + 1)) * sizeof(short));

    pMatrixChunkSize = (sizeUniqAB / numProcs);
    scoreMatrixChunkSize = ((sizeB + 1) / numProcs);

    start_time = MPI_Wtime();

    calc_P_matrix_v1(pMatrix, seqB, sizeB, uniqAB, sizeUniqAB, myRank, pMatrixChunkSize);
    // score = lcs_yang_v1(scoreMatrix, pMatrix, seqA, seqB, uniqAB, sizeA, sizeB, sizeUniqAB, myRank, scoreMatrixChunkSize);

    stop_time = MPI_Wtime();
    // Print Parallel LCS Score
    if (myRank == 0)
    {
        printf("%d;", score);
        // Print Serial LCS Time
        printf("%lf;\n", stop_time - start_time);
    }

    // free score matrix
    freeScoreMatrix(scoreMatrix, sizeB);
    free(pMatrix);

    MPI_Finalize();

    return EXIT_SUCCESS;
}
