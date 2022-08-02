#include <stdio.h>
#include <omp.h>
static long num_steps = 100000;
double step;

int main () {
    int i; double pi = 0.0;
    step = 1.0/(double) num_steps;
    double result[4];

    #pragma omp parallel num_threads(4)
    {
        int numeroThreads = omp_get_num_threads();
        int ID = omp_get_thread_num();
        int chunk = num_steps/numeroThreads;
        double x,sum = 0.0;
        
        for (int i = chunk*ID; i < chunk*ID + chunk-1; i++)
        {
            x = (i + 0.5) * step; // Largura do retângulo
            sum = sum + 4.0 / (1.0 + x*x); // Sum += Área do retângulo
        }
        result[ID] = sum;
    }

    pi = (result[0] + result[1] + result[2] + result[3]) * step;
    printf("PI: (%f)\n",pi);
}   