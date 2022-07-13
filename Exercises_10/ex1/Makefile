CC = g++
MPI = mpicxx
CFLAGS = -Wall -O3 --std=c++11

ex1.exe : ex1.o random.o
	$(MPI) random.o ex1.o -o ex1.exe $(CFLAGS)
ex1.o : ex1.cpp ex1.h ../Random/random.h
	$(MPI) -c ex1.cpp -o ex1.o $(CFLAGS)
random.o : ../Random/random.cpp ../Random/random.h
	$(MPI) -c ../Random/random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o ex1.exe
