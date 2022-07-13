CC = g++
CFLAGS = -Wall -O3 --std=c++11

ex2.exe : ex2.o random.o
	$(CC) random.o ex2.o -o ex2.exe
ex2.o : ex2.cpp ex2.h
	$(CC) -c ex2.cpp -o ex2.o $(CFLAGS)
random.o : ../Random/random.cpp ../Random/random.h
	$(CC) -c ../Random/random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o ex2.exe output_*
