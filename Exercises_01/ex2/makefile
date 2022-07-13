CC = g++
CFLAGS = -Wall -O3 --std=c++11

ex2 : ex2.o ../Random/random.o
	$(CC) ../Random/random.o ex2.o -o ex2
ex2.o : ex2.cpp ex2.h
	$(CC) -c ex2.cpp -o ex2.o $(CFLAGS)
../Random/random.o : ../Random/random.cpp ../Random/random.h
	$(CC) -c ../Random/random.cpp -o ../Random/random.o $(CFLAGS)
clean :
	rm *.o ex2
