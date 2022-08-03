all: main
main: parallel.c
		gcc -O3 -Wall -pedantic -pthread -fopenmp parallel.c -o parallel
clean:
		rm -rf *.o