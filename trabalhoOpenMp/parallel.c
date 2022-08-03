#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>     
#include "omp.h"  

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

typedef unsigned short mtype;

/* Read sequence from a file to a char vector.
 Filename is passed as parameter */
void print_matrix(mtype **x, int row, int col)
{
    for(int i=0;i<row;i++)
    {
        for(int j=0;j<col;j++)
        {
            printf("%d ",x[i][j]);
        }
        printf("\n");
    }
}
int get_index_of_character(char *str,char x, int len)
{
    for(int i=0;i<len;i++)
    {
        if(str[i]== x)
        {
            return i;
        }
    }
    return -1;//not found the character x in str
}
mtype lcs_yang_v1(mtype **DP, mtype **P, char *A, char *B, char *C, int m, int n, int u)
{
    for(int i=1;i<m+1;i++)
    {
        int c_i = get_index_of_character(C,A[i-1],u);
        #pragma omp parallel for schedule(static) num_threads(1)
        for(int j=0;j<n+1;j++)
        {
            if(A[i-1]==B[j-1])
            {
                DP[i][j] = DP[i-1][j-1] + 1;
            }
            else if(P[c_i][j]==0)
            {
                DP[i][j] = max(DP[i-1][j], 0);
            }
            else
            {
                DP[i][j] = max(DP[i-1][j], DP[i-1][P[c_i][j]-1] + 1);
            }
        }
    }
    return DP[m][n];
}
void calc_P_matrix_v1(mtype **P, char *b, int len_b, char *c, int len_c)
{
    #pragma omp parallel for
    for(int i=0;i<len_c;i++)
    {
        for(int j=2;j<len_b+1;j++)
        {
            if(b[j-2]==c[i]) //j-2 as b we assume here that b has a empty character in the beginning
            {
                P[i][j] = j-1;
            }
            else
            {
                P[i][j] = P[i][j-1];
            }
        }
    }
}
char* read_seq(char *fname) {
	//file pointer
	FILE *fseq = NULL;
	//sequence size
	long size = 0;
	//sequence pointer
	char *seq = NULL;
	//sequence index
	int i = 0;

	//open file
	fseq = fopen(fname, "rt");
	if (fseq == NULL ) {
		printf("Error reading file %s\n", fname);
		exit(1);
	}

	//find out sequence size to allocate memory afterwards
	fseek(fseq, 0L, SEEK_END);
	size = ftell(fseq);
	rewind(fseq);

	//allocate memory (sequence)
	seq = (char *) calloc(size + 1, sizeof(char));
	if (seq == NULL ) {
		printf("Erro allocating memory for sequence %s.\n", fname);
		exit(1);
	}

	//read sequence from file
	while (!feof(fseq)) {
		seq[i] = fgetc(fseq);
		if ((seq[i] != '\n') && (seq[i] != EOF))
			i++;
	}
	//insert string terminator
	seq[i] = '\0';

	//close file
	fclose(fseq);

	//return sequence pointer
	return seq;
}

mtype ** allocateScoreMatrix(int cols, int rows) {
	int i;
	//Allocate memory for LCS score matrix
	mtype ** scoreMatrix = (mtype **) malloc((rows + 1) * sizeof(mtype *));
	for (i = 0; i < (rows + 1); i++)
		scoreMatrix[i] = (mtype *) malloc((cols + 1) * sizeof(mtype));
	return scoreMatrix;
}

void initScoreMatrix(mtype ** scoreMatrix, int sizeA, int sizeB) {
	int i, j;
	//Fill first line of LCS score matrix with zeroes
	for (j = 0; j < (sizeA + 1); j++)
		scoreMatrix[0][j] = 0;

	//Do the same for the first collumn
	for (i = 1; i < (sizeB + 1); i++)
		scoreMatrix[i][0] = 0;
}
void clearScoreMatrix(mtype ** scoreMatrix, int rows, int cols){
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            scoreMatrix[i][j] = 0;
        }
    }
}
int LCS(mtype ** scoreMatrix, int sizeA, int sizeB, char * seqA, char *seqB) {
	int i, j;
	for (i = 1; i < sizeB + 1; i++) {
		for (j = 1; j < sizeA + 1; j++) {
			if (seqA[j - 1] == seqB[i - 1]) {
				/* if elements in both sequences match,
				 the corresponding score will be the score from
				 previous elements + 1*/
				scoreMatrix[i][j] = scoreMatrix[i - 1][j - 1] + 1;
			} else {
				/* else, pick the maximum value (score) from left and upper elements*/
				scoreMatrix[i][j] =max(scoreMatrix[i-1][j], scoreMatrix[i][j-1]);
			}
		}
	}
	return scoreMatrix[sizeB][sizeA];
}
void printMatrix(char * seqA, char * seqB, mtype ** scoreMatrix, int sizeA,
		int sizeB) {
	int i, j;

	//print header
	printf("Score Matrix:\n");
	printf("========================================\n");

	//print LCS score matrix allong with sequences

	printf("    ");
	printf("%5c   ", ' ');

	for (j = 0; j < sizeA; j++)
		printf("%5c   ", seqA[j]);
	printf("\n");
	for (i = 0; i < sizeB + 1; i++) {
		if (i == 0)
			printf("    ");
		else
			printf("%c   ", seqB[i - 1]);
		for (j = 0; j < sizeA + 1; j++) {
			printf("%5d   ", scoreMatrix[i][j]);
		}
		printf("\n");
	}
	printf("========================================\n");
}

void freeScoreMatrix(mtype **scoreMatrix, int sizeB) {
	int i;
	for (i = 0; i < (sizeB + 1); i++)
		free(scoreMatrix[i]);
	free(scoreMatrix);
}

int main(int argc, char ** argv) {
	double time_spent = 0.0;

	// sequence pointers for both sequences
	char *seqA, *seqB ,*uniqAB;

	// sizes of both sequences
	int sizeA, sizeB, sizeUniqAB;

	//read both sequences
	seqA = read_seq("./entradas/fileA.in"); // rows
	seqB = read_seq("./entradas/fileB.in"); // colunms
    uniqAB = read_seq("entradas/uniqAB.in");

	//find out sizes
	sizeA = strlen(seqA);
	sizeB = strlen(seqB);
    sizeUniqAB = strlen(uniqAB);

	// allocate LCS score matrix
	mtype ** scoreMatrix = allocateScoreMatrix(sizeA, sizeB);

    clock_t begin = clock();
	//initialize LCS score matrix
	initScoreMatrix(scoreMatrix, sizeA, sizeB);

	//fill up the rest of the matrix and return final score (element locate at the last line and collumn)
	mtype score = LCS(scoreMatrix, sizeA, sizeB, seqA, seqB);

	/* if you wish to see the entire score matrix,
	 for debug purposes, define DEBUGMATRIX. */
#ifdef DEBUGMATRIX
	printMatrix(seqA, seqB, scoreMatrix, sizeA, sizeB);
#endif
    
    clock_t end = clock();
    // calculate elapsed time by finding difference (end - begin) and
    // dividing the difference by CLOCKS_PER_SEC to convert to seconds
    time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
    //print score
    printf("\nScore: %d\n", score);    
    printf("The elapsed time for normal LCS is %f seconds\n", time_spent);

	
	
    clearScoreMatrix(scoreMatrix, sizeB+1, sizeA+1);

    // -----------------------------------------------------------------------------------------------------------------------------------

    mtype ** pMatrix = allocateScoreMatrix(sizeA, sizeUniqAB-1);

    mtype ** DP_Results = (mtype **)malloc((sizeB+1) * sizeof(mtype *));
    for(int k=0;k<sizeB+1;k++)
    {
        DP_Results[k] = (mtype *)calloc((sizeA+1), sizeof(mtype));
    }
    mtype ** P_Matrix = (mtype **)malloc(sizeUniqAB * sizeof(mtype *));
    for(int k=0;k<sizeUniqAB;k++)
    {
        P_Matrix[k] = (mtype *)calloc((sizeA+1), sizeof(mtype));
    }

    clearScoreMatrix(pMatrix,sizeUniqAB,sizeA);

    double start_time = omp_get_wtime();
    calc_P_matrix_v1(P_Matrix,seqA,sizeA,uniqAB,sizeUniqAB);
    printf("lcs is: %d\n",lcs_yang_v1(DP_Results,P_Matrix,seqB,seqA,uniqAB,sizeB,sizeA,sizeUniqAB));

    // calc_P_matrix_v1(pMatrix,seqA,sizeA,uniqAB,sizeUniqAB);
    // score = lcs_yang_v1(scoreMatrix,pMatrix,seqB,seqA,uniqAB,sizeB,sizeA,sizeUniqAB);
    printf("\nScore: %d\n", score);
 	
    

    //print score
	printf("\nParallel Score: %d\n", score);

 	double  stop_time = omp_get_wtime();
 
    printf("\ntime taken by parallel algorithm is: %lf\n",stop_time-start_time);
    
	//free score matrix
	freeScoreMatrix(scoreMatrix, sizeB);

	return EXIT_SUCCESS;
}
