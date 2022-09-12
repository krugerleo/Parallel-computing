#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "mpi.h"

#ifndef max
#define max(a, b) (((a) > (b)) ? (a) : (b))
#endif

int get_index_of_character(char *str, char x, int len);
void print_matrix(short **x, int row, int col);
void calcPMatrix(short *P, char *b, int len_b, char *c, int len_c, int myrank, int chunk_size);
int lcsMPI(short **scoreMatrix, short *P, char *A, char *B, char *C, int m, int n, int u, int myrank, int chunk_size);
int lcs(short **scoreMatrix, char *A, char *B, int m, int n);

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

void calcPMatrix(short *P, char *b, int len_b, char *c, int len_c, int myrank, int chunk_size)
{
    char bufferToReceiveAlphabet[chunk_size];
    short bufferToReceivePmatrix[chunk_size * (len_b + 1)];
    // Send chunk
    MPI_Scatter(c, chunk_size, MPI_CHAR, &bufferToReceiveAlphabet, chunk_size, MPI_CHAR, 0, MPI_COMM_WORLD);
    // Send chunk
    MPI_Scatter(P, chunk_size * (len_b + 1), MPI_SHORT, &bufferToReceivePmatrix, chunk_size * (len_b + 1), MPI_SHORT, 0, MPI_COMM_WORLD);
    // Broadcast B
    // max throughput
    MPI_Bcast(b, len_b, MPI_CHAR, 0, MPI_COMM_WORLD);

    for (int i = 0; i < chunk_size; i++)
    {
        for (int j = 2; j < len_b + 1; j++)
        {
            if (b[j - 2] == bufferToReceiveAlphabet[i]) // j-2 as b we assume here that b has a empty character in the beginning
            {
                bufferToReceivePmatrix[(i * (len_b + 1)) + j] = j - 1;
            }
            else
            {
                bufferToReceivePmatrix[(i * (len_b + 1)) + j] = bufferToReceivePmatrix[(i * (len_b + 1)) + j - 1];
            }
        }
    }

    // Juntar contrario do scatter
    MPI_Gather(bufferToReceivePmatrix, chunk_size * (len_b + 1), MPI_SHORT, P, chunk_size * (len_b + 1), MPI_SHORT, 0, MPI_COMM_WORLD);
}

int lcsMPI(short **scoreMatrix, short *P, char *A, char *B, char *C, int m, int n, int u, int myrank, int chunk_size)
{

    MPI_Bcast(A, m, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(P, (u * (n + 1)), MPI_SHORT, 0, MPI_COMM_WORLD);
    for (int i = 1; i < m + 1; i++)
    {
        int c_i = get_index_of_character(C, A[i - 1], u);

        short scoreMatrix_i_receive[chunk_size];

        MPI_Scatter(scoreMatrix[i], chunk_size, MPI_SHORT, &scoreMatrix_i_receive, chunk_size, MPI_SHORT, 0, MPI_COMM_WORLD);
        int start_id = (myrank * chunk_size);
        int end_id = (myrank * chunk_size) + chunk_size;
        for (int j = start_id; j < end_id; j++) // if myrank=0 then j=start_id+1 else j=start_id
        {
            if (j == start_id && myrank == 0)
                j = j + 1;
            if (A[i - 1] == B[j - 1])
            {
                scoreMatrix_i_receive[j - start_id] = scoreMatrix[i - 1][j - 1] + 1;
            }
            else if (P[(c_i * (n + 1)) + j] == 0)
            {
                scoreMatrix_i_receive[j - start_id] = max(scoreMatrix[i - 1][j], 0);
            }
            else
            {
                scoreMatrix_i_receive[j - start_id] = max(scoreMatrix[i - 1][j], scoreMatrix[i - 1][P[(c_i * (n + 1)) + j] - 1] + 1);
            }
        }
        // Junta em todos os processos
        MPI_Allgather(scoreMatrix_i_receive, chunk_size, MPI_SHORT, scoreMatrix[i], chunk_size, MPI_SHORT, MPI_COMM_WORLD);
    }
    return scoreMatrix[m - 1][n - 1];
}

void print_matrix(short **x, int row, int col)
{
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            printf("%d ", x[i][j]);
        }
        printf("\n");
    }
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
    int i;
    // Allocate memory for LCS score matrix
    short **scoreMatrix = (short **)malloc((rows + 1) * sizeof(short *));
    for (i = 0; i < (rows + 1); i++)
        scoreMatrix[i] = (short *)malloc((cols + 1) * sizeof(short));
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
void printMatrix(char *seqA, char *seqB, short **scoreMatrix, int sizeA,
                 int sizeB)
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
    int my_rank;
    int num_procs;
    int chunk_size_p, chunk_size_scoreMatrix; // chunk_size for P matrix and scoreMatrix matrix
    short score;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);   // grab this process's rank
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs); // grab the total num of processes

    // sequence pointers for both sequences
    char *seqA, *seqB, *uniqAB;

    // sizes of both sequences
    int sizeA, sizeB, sizeUniqAB;

    seqA = read_seq(argv[1]);   // rows
    seqB = read_seq(argv[2]);   // colunms
    uniqAB = read_seq(argv[3]); // alphabet

    // read both sequences

    // find out sizes
    sizeA = strlen(seqA);
    sizeB = strlen(seqB);
    sizeUniqAB = strlen(uniqAB);

    chunk_size_p = (sizeUniqAB / num_procs);
    chunk_size_scoreMatrix = ((sizeB + 1) / num_procs);
    // allocate LCS score matrix
    short **scoreMatrix = allocateScoreMatrix(sizeA, sizeB);

    // initialize LCS score matrix
    initScoreMatrix(scoreMatrix, sizeA, sizeB);
    double start_time, stop_time;

    if (my_rank == 0)
    {
        start_time = MPI_Wtime();
        score = LCS(scoreMatrix, sizeA, sizeB, seqA, seqB);
        stop_time = MPI_Wtime();

        // Print Iput Sizes
        printf("(%d,%d,%d);", sizeB, sizeA, sizeUniqAB);
        // Print Serial LCS Score
        printf("%d;", score);
        // Print Serial LCS Time
        printf("%f;", stop_time - start_time);
    }

    clearScoreMatrix(scoreMatrix, sizeB + 1, sizeA + 1);

    // -----------------------------------------------------------------------------------------------------------------------------------

    short *pMatrix = (short *)malloc((sizeUniqAB * (sizeB + 1)) * sizeof(short));

    // clearScoreMatrix(pMatrix, sizeUniqAB, sizeA)

    start_time = MPI_Wtime();

    calcPMatrix(pMatrix, seqB, sizeB, uniqAB, sizeUniqAB, my_rank, chunk_size_p);
    score = lcsMPI(scoreMatrix, pMatrix, seqA, seqB, uniqAB, sizeA, sizeB, sizeUniqAB, my_rank, chunk_size_scoreMatrix);

    stop_time = MPI_Wtime();

    if (my_rank == 0)
    {
        // Print Parallel LCS Score
        printf("%d;", score);
        // Print Serial LCS Time
        printf("%lf;", stop_time - start_time);
        printf("%d;", num_procs);
    }

    // free score matrix
    free(scoreMatrix);
    free(pMatrix);
    // Shutdown MPI (important - don't forget!)
    MPI_Finalize();
    return EXIT_SUCCESS;
}
