CC = g++
CFLAGS = -Wall -O3 --std=c++11

ex3 : ex3.o ../Random/random.o
	$(CC) ../Random/random.o ex3.o -o ex3
ex3.o : ex3.cpp ex3.h
	$(CC) -c ex3.cpp -o ex3.o $(CFLAGS)
../Random/random.o : ../Random/random.cpp ../Random/random.h
	$(CC) -c ../Random/random.cpp -o ../Random/random.o $(CFLAGS)
clean :
	rm *.o ex3
