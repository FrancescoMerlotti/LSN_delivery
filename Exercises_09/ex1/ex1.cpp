#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>
#include "ex1.h"

using namespace std;

int main(int argc, char** argv) {
	time_t t_start = time(NULL), t_end;
	// In line arguments check
	if(argc != 2){
		cout << endl << "Right usage: ./ex1.exe <nseed>" << endl << endl;
		return -1;
	}
	// Input
	Input(atoi(argv[1]));
	// Generation
	if(gen){
		Generation(atoi(argv[3]));
		return -2;
	}
	// Simulation performing
	Loading();
	Configurate();
	Order();
	best = population[0];
	// Display();
	for(int istep = 0; istep < nsteps; istep++) {
		// Crossing over between the selected chromosomes
		if(rnd.Rannyu() <= 0.8) {
			Selection();
			// population[nconf - 1] = CrossingOver(in1, in2);
			CrossingOver2(in1, in2);
		}
		// Ordering the new population
		Order();
		// Mutations on the selected chromosome
		// Permutation
		if(rnd.Rannyu() <= 0.1)
			Permutation(int(0.5 * nconf * rnd.Rannyu()));
		// Permutation Cluster
		if(rnd.Rannyu() <= 0.1)
			PermutationCluster(int(0.5 * nconf * rnd.Rannyu()));
		// Reverse
		if(rnd.Rannyu() <= 0.1)
			Reverse(int(0.5 * nconf * rnd.Rannyu()));
		// Queue
		if(rnd.Rannyu() <= 0.1)
			Queue(int(0.5 * nconf * rnd.Rannyu()));
		// Tracing the best solution
		IsBest();
		// Printing on file
		PrintFile(istep);
	}
	// Printing the best path
	PrintBest();
	// Printing on screen the length of the best path
	cout << "Best path length = " << PathLength(best) << endl << endl;
	// Estimating execution time
	t_end = time(NULL);
	cout << "Execution time: " << difftime(t_end, t_start) << " seconds" << endl << endl;
	return 0;
}

// Input
void Input(int nseed) {
	// Read seed for random numbers
	ifstream Primes, Seed, Setup;
	int p1, p2;
	Primes.open("../Random/primes32001.dat");
	for(int i = 1; i <= nseed; i++)
    	Primes >> p1 >> p2;
	Primes >> p1 >> p2;
	Primes.close();
	Seed.open("../Random/seed.in");
	Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	Seed.close();
	// Random setup
	rnd.SetRandom(seed, p1, p2);
	Setup.open("setup.dat");
	Setup >> ncities >> nconf >> nsteps >> gen >> shape;
	Setup.close();
	cout << endl << "Performing a " << nsteps << " steps genetic simulation with " << ncities;
	if(shape)
		cout << " cities in a Square, ";
	else
		cout << " cities on a Circle, ";
	cout << "using a population of " << nconf << " elements" << endl << endl;
	
}

// Generation
void Generation(int conf) {
	// conf = 0: circle, conf = 1: square
	ofstream Cities;
	const int wd = 12;
	cities = new City[ncities];
	Cities.open("./output/cities.dat");
	// Cities generation and output
	for(int ic = 0; ic < ncities; ic++) {
		cities[ic].SetLabel(ic+1);
		// Circle
		if(conf == 0) {
			double theta = rnd.Rannyu(0., 2*M_PI);
			cities[ic].SetX(cos(theta));
			cities[ic].SetY(sin(theta));
		}
		// Square
		if(conf == 1) {
			cities[ic].SetX(rnd.Rannyu(-1., 1.));
			cities[ic].SetY(rnd.Rannyu(-1., 1.));
		}
		// Output on file
		Cities << setw(4) << cities[ic].GetLabel() << setw(wd) << cities[ic].GetX() << setw(wd) << cities[ic].GetY() << endl;
	}
	Cities.close();
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
	int al;
	ifstream Cities;
	cities = new City[ncities];
	Cities.open("./output/cities.dat");
	for(int ic = 0; ic < ncities; ic++) {
		Cities >> al >> ax >> ay;
		cities[ic].SetLabel(al);
		cities[ic].SetX(ax);
		cities[ic].SetY(ay);
	}
	Cities.close();
	cout << "Cities successfully loaded" << endl << endl;
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
	cout << "Configuration successfully completed" << endl << endl;
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

// Crossing over
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

	if(SwapRule(g1, g2))
		return g1;
	else
		return g2;
}

// Crossing over
void CrossingOver2(int i1, int i2) {

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

	population[nconf - 2] = g1;
	population[nconf - 1] = g2;
}

// Selection
void Selection() {
	in1 = int(nconf * pow(rnd.Rannyu(), 2));
	in2 = int(nconf * pow(rnd.Rannyu(), 2));
	while(in1 == in2)
			in2 = int(nconf * pow(rnd.Rannyu(), 2));
}

// Print data on file
void PrintFile(int istep) {
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
		Best.open("./output/output_best.dat");
		Mean.open("./output/output_mean.dat");
	}
	else {
		Best.open("./output/output_best.dat", ios::app);
		Mean.open("./output/output_mean.dat", ios::app);
	}
	// Printing the best path of the current population
	Best << setw(5) << istep + 1 << setw(10) << PathLength(population[0]) << endl;
	// Printing the mean path length of the best half of the population
	Mean << setw(5) << istep + 1 << setw(10) << appo / half << endl;
	// Closing stream
	Best.close();
	Mean.close();
}

// Print best on file
void PrintBest() {
	ofstream Output;
	Output.open("./output/output_bestPath.dat");
	Output << "Best path (Length = " << setw(8) << PathLength(population[0]) << "):" << endl;
	for(int ic = 0; ic < ncities; ic++)
		Output << setw(4) << population[0][ic] << endl;
	Output.close();
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
