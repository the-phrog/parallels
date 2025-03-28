CC=g++
FLAGS=-fopenmp -o
SER=serial
PAR=parallel
PARF=parallel_for

all:
	$(CC) $(SER).cpp $(FLAGS) $(SER)

par:
	$(CC) $(PAR).cpp $(FLAGS) $(PAR)

parf:
	$(CC) $(PARF).cpp $(FLAGS) $(PARF)

clean:
	rm -f $(SER)
	rm -f $(PAR)
	rm -f $(PARF)
