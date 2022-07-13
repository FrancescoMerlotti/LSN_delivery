#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include "mpi.h"

using namespace std;

int Pbc(int, int);

int main(int argc, char ** argv) {

	// Parallelization
	int size, rank;
	const int nconf = 9, nmigr = 4, ncities = 10;
	int cost = pow(2, nmigr - 1);
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int ** pop = new int*[nconf];
	for(int ic = 0; ic < nconf; ic++)
		pop[ic] = new int[ncities];
	
	for(int i = 0; i < nconf; i++)
		for(int j = 0; j < ncities; j++)
			pop[i][j] = (i + 1) * 10 + j;

	if(size % 2 != 0) {
		cout << endl << "Processes need to be even" << endl;
		return -3;
	}

	MPI_Status stat0, stat1;
	// MPI_Request req;
	if(rank % 2 == 1) {
		MPI_Send(pop[0], (ncities) * nmigr + cost, MPI_INTEGER, Pbc(rank+1, size), rank, MPI_COMM_WORLD);
		MPI_Recv(pop[nconf - nmigr], (ncities) * nmigr + cost, MPI_INTEGER, Pbc(rank-1, size), Pbc(rank-1, size), MPI_COMM_WORLD, &stat0);
	}
	else if(rank % 2 == 0) {
		MPI_Recv(pop[nconf - nmigr], (ncities) * nmigr + cost, MPI_INTEGER, Pbc(rank-1, size), Pbc(rank-1, size), MPI_COMM_WORLD, &stat1);
		MPI_Send(pop[0], (ncities) * nmigr + cost, MPI_INTEGER, Pbc(rank+1, size), rank, MPI_COMM_WORLD);
		
	}

	if(rank == 0) {
		cout << endl << "Population in node " << rank << " after migration:" << endl << endl;
		for(int i = 0; i < nconf; i++) {
			cout << setw(10) << i + 1;
			for(int j = 0; j < ncities; j++)
				cout << setw(10) << pop[i][j];
			cout << endl;
		}
	}

	MPI_Finalize();
	
	return 0;
}

// Periodic boundary condition on indeces
int Pbc(int rank, int size) {
	if(rank >= size)
		rank -= size;
	if(rank < 0)
		rank += size;
	return rank;
}