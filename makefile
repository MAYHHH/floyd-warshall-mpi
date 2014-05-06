CC= mpicc

CFLAGS= -Wall -g -fopenmp

install: 
	$(CC) $(CFLAGS) -o  floydWarshallMPI floydWarshallMPI.c -lm
	
clean: 	
	rm -f *.o palmieri floydWarshallMPI