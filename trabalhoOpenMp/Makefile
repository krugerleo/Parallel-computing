all: parallel
parallel: parallel.o 
		gcc -o parallel parallel.o -O3 -fopenmp

main.o: parallel.c
		gcc -o parallel.o -c parallel.c -O3 -fopenmp
clean:
		rm -rf *.o
mrproper: clean
		rm -rf teste