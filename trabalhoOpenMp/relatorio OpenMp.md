# Relatorio Trabalho 1 OpenMP
Aluno: Leonardo Bueno Nogueira Kruger
GRR20180130
## Introdução
## Funcionamento Core LCS
Algoritmo LCS (Longest Common Subsequence)
Algoritmo utilizado para encontrar a maior subsequencia presente em duas sequencias (Strings), uma subsequencia é caracteriza como uma sequencia que aparece na mesma ordem relativa mas não necessariamente continua.

O Algoritmo trabalha em cima de uma matriz de tamanho sizeA x sizeB onde sizeA é o tamanho da string A e sizeB é o tamanho da string B, com a primeira linha e coluna inicializadas em 0.
```C
for (j = 0; j < (sizeA + 1); j++)
		scoreMatrix[0][j] = 0;

for (i = 1; i < (sizeB + 1); i++)
    scoreMatrix[i][0] = 0;
```
A partir dessa estrutura ocorre o LCS sobre a matriz, percorrendo a matriz caso encontre um 'match' de character (char de seqA e seqB correspondem) é pego o valor na diagonal anterior e somado 1 para o tamanho da subsequencia, caso contrario é maior valor do campo superior(cima) ou anterior(esquerda), ao final teremos o valor da maior subsequencia.
```C
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
```
## Estratégia de paralelização
## Historico 
Após conhecimento sobre Algoritmo LCS (Longest Common Subsequence) é facil perceber uma dependencia de dados entre o elemento atual e anteriores para o funcionamento do algoritmo, logo a paralelização a seguir não demonstrou ganhos significativos de desempenho:
```C
#pragma omp parallel shared(scoreMatrix) private(i, k)
    {
        #pragma omp for
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
```


## Informações 