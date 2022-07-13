CC = g++
CFLAGS = -Wall -O3 --std=c++11

integralMonteCarlo : integralMonteCarlo.o random.o
	$(CC) random.o integralMonteCarlo.o -o integralMonteCarlo
integralMonteCarlo.o : integralMonteCarlo.cpp
	$(CC) -c integralMonteCarlo.cpp -o integralMonteCarlo.o $(CFLAGS)

random.o : ../../Random/random.cpp ../../Random/random.h
	$(CC) -c ../../Random/random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o ex1 ex2
