#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

typedef unsigned short mtype;

/* Read sequence from a file to a char vector.
 Filename is passed as parameter */
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

mtype ** allocateScoreMatrix(int sizeA, int sizeB) {
	int i;
	//Allocate memory for LCS score matrix
	mtype ** scoreMatrix = (mtype **) malloc((sizeA + 1) * sizeof(mtype *));
	for (i = 0; i < (sizeA + 1); i++)
		scoreMatrix[i] = (mtype *) malloc((sizeB + 1) * sizeof(mtype));
	return scoreMatrix;
}
mtype ** allocatePMatrix(int sizeA, int sizeB) {
	//Allocate memory for LCS score matrix
	mtype ** scoreMatrix = (mtype **) malloc((sizeA) * sizeof(mtype *));
	for (int i = 0; i < (sizeA); i++)
		scoreMatrix[i] = (mtype *) malloc((sizeB + 1) * sizeof(mtype));
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
void initPMatrix(mtype **P, char *b, int len_b, char *c, int len_c) {
    #pragma omp parallel for
    for(int i=0;i<len_c;i++)
    {
        for(int j=1;j<len_b+1;j++)
        {
            if(b[j-1]==c[i])
            {
                P[i][j] = j;
            }
            else
            {
                P[i][j] = P[i][j-1];
            }
        }
    }
}
void clearScoreMatrix(mtype ** scoreMatrix, int sizeA, int sizeB){
    for (int i = 1; i < sizeB + 1; i++) {
        for (int j = 1; j < sizeA + 1; j++) {
            scoreMatrix[i][j] = 0;
        }
    }
}
int parallelLCS(mtype ** scoreMatrix, mtype ** pMatrix, int sizeA, int sizeB, int sizeUniqAB, char * seqA, char *seqB, char *uniqAB) {
    for(int i=1;i<sizeB+1;i++)
    {
        int c_i = get_index_of_character(uniqAB,seqB[i-1],sizeUniqAB);
        #pragma omp parallel for schedule(static)
	for(int j=0;j<sizeA+1;j++)
        {
            if(seqA[i-1]==seqB[j-1])
            {
                scoreMatrix[i][j] = scoreMatrix[i-1][j-1] + 1;
            }
            else if(pMatrix[c_i][j]==0)
            {
                scoreMatrix[i][j] = max(scoreMatrix[i-1][j], 0);
            }
            else
            {
                scoreMatrix[i][j] = max(scoreMatrix[i-1][j], scoreMatrix[i-1][pMatrix[c_i][j]-1] + 1);
            }
        }
    }
    return scoreMatrix[sizeB][sizeA];
}
int LCS(mtype ** scoreMatrix, int sizeA, int sizeB, char * seqA, char *seqB) {
	int i, j;

    {
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
	double start_time,stop_time;

	// sequence pointers for both sequences
	char *seqA, *seqB, *uniqAB;

	// sizes of both sequences
	int sizeA, sizeB, sizeUniqAB;

	//read both sequences
	seqA = read_seq("entradas/fileA.in");
	seqB = read_seq("entradas/fileB.in");
    uniqAB = read_seq("entradas/uniqAB.in");

	//find out sizes
	sizeA = strlen(seqA); // colunms
	sizeB = strlen(seqB); // rows
    sizeUniqAB = strlen(uniqAB);

	// allocate LCS score matrix
	mtype ** scoreMatrix = allocateScoreMatrix(sizeA, sizeB);
    mtype ** pMatrix = allocatePMatrix(sizeUniqAB, sizeA);
	start_time = omp_get_wtime();
	//initialize LCS score matrix
	initScoreMatrix(scoreMatrix, sizeA, sizeB);
    
	//fill up the rest of the matrix and return final score (element locate at the last line and collumn)
	mtype score = LCS(scoreMatrix, sizeA, sizeB, seqA, seqB);
	
    stop_time = omp_get_wtime();

    printf("time taken Paralel algorithm is: %lf",stop_time-start_time);
	printf("\nScore: %d\n", score);

	/* if you wish to see the entire score matrix,
	 for debug purposes, define DEBUGMATRIX. */
#ifdef DEBUGMATRIX
	printMatrix(seqA, seqB, scoreMatrix, sizeA, sizeB);
#endif
    clearScoreMatrix(scoreMatrix, sizeA, sizeB);

	start_time = omp_get_wtime();

    initPMatrix(pMatrix,seqA,sizeA,uniqAB, sizeUniqAB);
    mtype parallelScore = parallelLCS(scoreMatrix,pMatrix,sizeA,sizeB,sizeUniqAB,seqA,seqB,uniqAB);

    stop_time = omp_get_wtime();
    printf("time taken Paralel algorithm is: %lf",stop_time-start_time);
	//print score
    printf("\nScore: %d\n", parallelScore);
	
	//free score matrix
	freeScoreMatrix(scoreMatrix, sizeB);
    

	return EXIT_SUCCESS;
}
