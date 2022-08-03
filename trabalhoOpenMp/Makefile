all: parallel
parallel: parallel.o 
		gcc -o parallel parallel.o  -fopenmp

main.o: parallel.c
		gcc -o parallel.o -c parallel.c  -fopenmp
clean:
		rm -rf *.o
mrproper: clean
		rm -rf teste
# gcc -O3 -Wall -pedantic -pthread -fopenmp parallel.c -o parallel