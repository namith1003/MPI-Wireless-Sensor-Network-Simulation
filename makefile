CC = mpicc
CFLAGS = -Wall
OMPFLAGS = -fopenmp

all: ev ev

ev: ev.c
	$(CC) -o ev ev.c $(CFLAGS) $(OMPFLAGS)


run: ev 
	./mpirun -np 26 ./ev 5 5

clean:
	rm -f ass1 ass1p
