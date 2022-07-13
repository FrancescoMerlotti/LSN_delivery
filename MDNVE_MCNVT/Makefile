CC = g++
CFLAGS = -Wall -O3 --std=c++11

NVE_NVT.exe : MD_MC.o random.o
	$(CC) random.o MD_MC.o -o NVE_NVT.exe
MD_MC.o : MD_MC.cpp MD_MC.h
	$(CC) -c MD_MC.cpp -o MD_MC.o $(CFLAGS)
random.o : ./Random/random.cpp ./Random/random.h
	$(CC) -c ./Random/random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o NVE_NVT.exe
