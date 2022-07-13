#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include "mpi.h"
#include "ex1.h"

using namespace std;

int main(int argc, char** argv) {

	// Parallelization
	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double t_beg = MPI_Wtime();

	if(argc != 2) {
		cout << endl << "Correct usage: mpiexec -np <nproc> ./ex1.exe <nseed>" << endl;
		return -1;
	}

	if(size % 2 != 0) {
		cout << endl << "Processes need to be even" << endl;
		return -3;
	}

	ofstream Output;
	// Input
	Input(rank, atoi(argv[1]));
	// Simulation performing
	Loading();
	Configurate();
	Order();
	best = population[0];
	// Display();
	for(int istep = 0; istep < nsteps; istep++) {
		// Crossing over between the selected chromosomes (steady state)
		Selection();
		population[ncities - 1] = CrossingOver(in1, in2);
		// Ordering the new population
		Order();		
		// Mutations on the selected chromosome
		// Permutation
		if(rnd.Rannyu() <= 0.1)
			Permutation(int(0.5 * ncities * rnd.Rannyu()));
		// Permutation Cluster
		if(rnd.Rannyu() <= 0.1)
			PermutationCluster(int(0.5 * ncities * rnd.Rannyu()));
		// Reverse
		if(rnd.Rannyu() <= 0.1)
			Reverse(int(0.5 * ncities * rnd.Rannyu()));
		// Queue
		if(rnd.Rannyu() <= 0.1)
			Queue(int(0.5 * ncities * rnd.Rannyu()));
		// Tracing the best solution
		IsBest();
		// Migration
		if(istep % ntime == 0) {
			MPI_Status stat0, stat1;
			// MPI_Request req;
			LoadPopulation();
			if(rank % 2 == 0) {
				for(int im = 0; im < nmigr; im++) {
					MPI_Send(pop[im], ncities, MPI_INTEGER, Pbc(rank+1, size), im * 10 + rank, MPI_COMM_WORLD);
					MPI_Recv(pop[nconf - nmigr + im], ncities, MPI_INTEGER, Pbc(rank-1, size), im * 10 + Pbc(rank-1, size), MPI_COMM_WORLD, &stat1);
				}
			}
			else if(rank % 2 == 1) {
				for(int im = 0; im < nmigr; im++) {
					MPI_Recv(pop[nconf - nmigr + im], ncities, MPI_INTEGER, Pbc(rank-1, size), im * 10 + Pbc(rank-1, size), MPI_COMM_WORLD, &stat0);
					MPI_Send(pop[im], ncities, MPI_INTEGER, Pbc(rank+1, size), im * 10 + rank, MPI_COMM_WORLD);
				}
			}
			SavePopulation();
		}
		// Printing on file
		PrintFile(istep, rank);
	}
	// Printing the best path
	Output.open("./output/bestPath/output_bestPath" + to_string(rank) + ".dat");
	Output << "Best path (Length = " << setw(8) << PathLength(population[0]) << "):" << endl;
	for(int ic = 0; ic < ncities; ic++)
		Output << setw(4) << population[0][ic] << endl;
	Output.close();

	double t_end = MPI_Wtime();
	cout << endl << "Execution time for node " << rank << ": " << t_end - t_beg << endl;

	MPI_Finalize();

	return 0;
}

// Input
void Input(int rank, int nseed) {
	// Read seed for random numbers
	ifstream Primes, Seed, Setup;
	int p1, p2;
	Setup.open("setup.dat");
	Setup >> ncities >> nconf >> nsteps >> nmigr >> ntime;
	Setup.close();
	pop = new int*[nconf];
	for(int ic = 0; ic < nconf; ic++)
		pop[ic] = new int[ncities];
	Primes.open("../Random/Primes");
	for(int i = 0; i < rank * nseed; i++)
    	Primes >> p1 >> p2;
	Primes >> p1 >> p2;
	Primes.close();
	Seed.open("../Random/Seed");
	Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	Seed.close();
	// Random setup
	rnd.SetRandom(seed, p1, p2);
	cout << endl << "Node: " << rank << ", Performing a " << nsteps << " steps genetic simulation with " << ncities;
	cout << " cities, using a population of " << nconf << " elements" << endl;
}

// Check
bool Check(vector<int> p) {
	bool flag = true;
	int ic = 1;
	vector<int> appo(ncities-1, 0);
	if(p[0] != 1)
		flag = false;
	while(flag and ic < ncities){
		appo[p[ic]-2] += 1;
		ic++;
	}
	if(any_of(appo.begin(), appo.end(), [](int n){return n != 1;}))
		flag = false;
	return flag;
}

// Cities input
void Loading() {
	double ax, ay;
	// int al;
	string as, an;
	ifstream Cities;
	cities = new City[ncities];
	Cities.open("american_capitals.dat");
	if(Cities.fail()) {
		cout << "Error in american_capitals.dat opening!" << endl << endl;
		return;
	}
	// Cities loading
	for(int ic = 0; ic < ncities; ic++) {
		Cities >> as >> an >> ax >> ay;
		cities[ic].SetLabel(ic+1);
		cities[ic].SetX(ax);
		cities[ic].SetY(ay);
	}
	Cities.close();
	// cout << endl << "Cities successfully loaded" << endl;
}

// Configurate
void Configurate() {
	int num;
	// Vector initialization
	vector<int> ac;
	for(int ic = 0; ic <= ncities - 2; ic++)
		ac.push_back(ic+2);
	// Building initial configurations
	for(int iconf = 0; iconf < nconf; iconf++) {
		vector<int> appo {1};
		for(int ic = 0; ic < ncities-1; ic++) {
			int o = int((ncities - ic - 1) * rnd.Rannyu());
			appo.push_back(ac[o]);
			num = ac[o];
			ac[o] = ac[ncities - ic - 2];
			ac[ncities - ic - 2] = num;
		}
		if(Check(appo))
			population.push_back(appo);
		for(int ic = 0; ic <= ncities - 2; ic++)
			ac[ic] = ic+2 ;
	}
	// cout << endl << "Configuration successfully completed" << endl;
}

// Population display
void Display() {
	cout << endl;
	for(int i = 0; i < nconf; i++) {
		cout << setw(4) << i + 1 << " - ";
		for(int j = 0; j < ncities; j++) {
			cout << setw(3) << population[i][j] << " ";
		}
		cout << "- Distance = " << setw(8) << PathLength(population[i]) << endl;
	}
	cout << endl;
}

// Distance between cities (L1)
double Distance(City& a, City &b) {
	return sqrt(pow(a.GetX() - b.GetX(), 2) + pow(a.GetY() - b.GetY(), 2)); 
}

// Path length
double PathLength(vector<int> p) {
	double appo = Distance(cities[p[ncities-1]-1], cities[p[0]-1]);
	for(int ic = 0; ic < ncities-1; ic++)
		appo += Distance(cities[p[ic]-1], cities[p[ic+1]-1]);
	return appo;
}

// Swap rule
bool SwapRule(vector<int> a, vector<int> b) {
	return (PathLength(a) < PathLength(b));
}

// Sort the paths in a population (increasing length of the path)
void Order() {
	sort(population.begin(), population.end(), SwapRule);
}

// Permutation between two cities in a path
void Permutation(int index) {
	int temp;
	int ic = int(rnd.Rannyu(1., double(ncities))), jc = int(rnd.Rannyu(1., double(ncities)));
	while(ic == jc)
		jc = int(rnd.Rannyu(1., double(ncities)));
	temp = population[index][ic];
	population[index][ic] = population[index][jc];
	population[index][jc] = temp;
}

// Permutation between two clusters of cities in a path
void PermutationCluster(int index) {
	int len = int(rnd.Rannyu(0., double((ncities - 1) / 2))) + 1;
	int pos = int(rnd.Rannyu(0., double(ncities - 2 * len))) + 1;
	int temp;

	for(int ic = 0; ic < len; ic++) {
		temp = population[index][pos + ic];
		population[index][pos + ic] = population[index][pos + len + ic];
		population[index][pos + len + ic] = temp;
	}
}

// Reverse of a cluster of cities in a path
void Reverse(int index) {
	int len = int(rnd.Rannyu(0., double(ncities - 1))) + 1;
	int pos = int(rnd.Rannyu(0., double(ncities - len))) + 1;
	int temp;
	for(int ic = 0; ic < len / 2; ic++) {
		temp = population[index][pos+ic];
		population[index][pos+ic] = population[index][pos+len-1-ic];
		population[index][pos+len-1-ic] = temp;
	}

}

// Queue a cluster of cities in a path
void Queue(int index) {
	vector<int> cluster, tail;
	int len = int(rnd.Rannyu(0., double(ncities - 1))) + 1;
	int pos = int(rnd.Rannyu(0., double(ncities - len - 1))) + 1;

	for(int ic = pos; ic < pos + len; ic++)
		cluster.push_back(population[index][ic]);
	for(int ic = pos+len; ic < ncities; ic++)
		tail.push_back(population[index][ic]);

	for(int ic = 0; ic < ncities - pos - len; ic++)
		population[index][pos + ic] = tail[ic];
	for(int ic = 0; ic < len; ic++)
		population[index][ncities - len + ic] = cluster[ic];
}

// Crossover
vector<int> CrossingOver(int i1, int i2) {

	int pos = int(rnd.Rannyu(1., double(ncities)));
	int cont = 0;
	vector<int> a1, a2, g1 = population[i1], g2 = population[i2];

	// Inserting elements of g2 in g1, checking if elements of the 2nd part of g2 are present in the 1st part of g1
	for(int ic = pos; ic < ncities; ic++) {
		bool flag = true;
		for(int jc = 0; flag and jc < pos; jc++) {
			if(population[i1][jc] == population[i2][ic]) {
				flag = false;
				a2.push_back(population[i2][ic]);
				cont += 1;
			}
		}
		if(flag)
			g1[ic- cont] = population[i2][ic];
	}
	cont = 0;
	// Inserting elements of g1 not present in g2
	for(int ic = pos; ic < ncities; ic++) {
		bool flag = true;
		for(int jc = 0; flag and jc < pos; jc++) {
			if(population[i2][jc] == population[i1][ic]) {
				flag = false;
				a1.push_back(population[i1][ic]);
				cont += 1;
			}
		}
		if(flag)
			g2[ic - cont] = population[i1][ic];	
	}

	// Inserting remainig elements for g1
	if(a1.size() != 0)
		for(int ic = 0; ic < a1.size(); ic++) 
			g1[ncities - 1 - ic] = a1[ic];
	// Inserting remainig elements for g2
	if(a2.size() != 0)
		for(int ic = 0; ic < a2.size(); ic++) 
			g2[ncities - 1 - ic] = a2[ic];
	
	// Double crossing over
	// if(rnd.Rannyu() <= (1. / 16.))		
	
	if(SwapRule(g1, g2))
		return g1;
	else
		return g2;
}

// Selection
void Selection() {
	in1 = int(ncities * pow(rnd.Rannyu(), 1.5));
	in2 = int(ncities * pow(rnd.Rannyu(), 1.5));
	while(in1 == in2)
		in2 = int(ncities * rnd.Rannyu());
}

// Print data on file
void PrintFile(int istep, int rank) {
	ofstream Best, Mean;
	double appo = 0;
	double half = nconf / 2;
	// Order
	Order();
	// Estimating mean path length of the best half of the population
	for(int ip = 0; ip < half; ip++)
		appo += PathLength(population[ip]);
	// Selection output mode
	if(istep == 0) {
		Best.open("./output/best/output_best" + to_string(rank) + ".dat");
		Mean.open("./output/mean/output_mean" + to_string(rank) + ".dat");
	}
	else {
		Best.open("./output/best/output_best" + to_string(rank) + ".dat", ios::app);
		Mean.open("./output/mean/output_mean" + to_string(rank) + ".dat", ios::app);
	}
	// Printing the best path of the current population
	Best << setw(5) << istep + 1 << setw(10) << PathLength(population[0]) << endl;
	// Printing the mean path length of the best half of the population
	Mean << setw(5) << istep + 1 << setw(10) << appo / half << endl;
	// Closing stream
	Best.close();
	Mean.close();
}

// Best path
bool IsBest() {
	Order();
	if(PathLength(population[0]) < PathLength(best)) {
		best = population[0];
		return true;
	}
	else 
		return false;
}

// Load population from vector on c-array
void LoadPopulation() {
	for(int ic = 0; ic < nconf; ic++)
		for(int jc = 0; jc < ncities; jc++)
			pop[ic][jc] = population[ic][jc];
}

// Save population from c-array on vector
void SavePopulation() {
	for(int ic = 0; ic < nconf; ic++)
		for(int jc = 0; jc < ncities; jc++)
			population[ic][jc] = pop[ic][jc];
}

// Periodic boundary condition on indeces
int Pbc(int rank, int size) {
	if(rank >= size)
		rank -= size;
	if(rank < 0)
		rank += size;
	return rank;
}

// City no arg constructor
City::City() {
	m_x = 0.;
	m_y = 0.;
	m_label = 0;
}

// City position and label constructor
City::City(double x, double y, int label) {
	m_x = x;
	m_y = y;
	m_label = label;
}

// City copy constructor
City::City(City& c) {
	m_x = c.GetX();
	m_y = c.GetY();
	m_label = c.GetLabel();
}

// City get x
double City::GetX() const{
	return m_x;
}

// City get y
double City::GetY() const{
	return m_y;
}

// City get label
int City::GetLabel() const{
	return m_label;
}

// City set x
void City::SetX(double x) {
	m_x = x;
}

// City set y
void City::SetY(double y) {
	m_y = y;
}

// City set label
void City::SetLabel(int label) {
	m_label = label;
}