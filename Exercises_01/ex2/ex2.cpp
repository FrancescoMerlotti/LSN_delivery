#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include "../Random/random.h"
#include "ex2.h"

using namespace std;

int main(int argc, char** argv) {

	if (argc != 2) {
		cout << endl << "Correct usage: ./ex2 seed" << endl << endl;
		return -1;
	}
	// Randon Number Generator Setup
	Input(atoi(argv[1]));
	// Standard Dice
	Standard();
	// Exponential Dice
	Exp();
	// Lorentzian Dice
	Lorentz();
	return 0;
}

// Input seutp
void Input(int nseed) {
	ifstream Primes, Seed;
	// Prime loading
	Primes.open("../Random/Primes");
	for(int i = 0; i < nseed; i++)
		Primes >> p1 >> p2;		
	Primes.close();
	// Seed loading
	Seed.open("../Random/Seed");
	Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	Seed.close();
	// Generator setup
	generator.SetRandom(seed, p1, p2);
}

// Standard dice
void Standard() {
	ofstream Output;
	Output.open("output_standard.dat");
	for (int iiter = 0; iiter < niter; iiter++) {
		vector<double> avg(4, 0);
		Output << setw(wd) << iiter;
		for (int i = 0; i < 4; i++) {
			double appo = 0.;
			for (int ithrow = 0; ithrow < nthrows[i]; ithrow++)
				appo += int(6 * generator.Rannyu()) + 1;
			avg[i] = appo / nthrows[i];
			Output << setw(wd) << avg[i];
		}
		Output << endl;
	}
	Output.close();
}

// Exponential dice
void Exp() {
	ofstream Output;
	Output.open("output_exponential.dat");
	for (int iiter = 0; iiter < niter; iiter++) {
		vector<double> avg(4, 0);
		Output << setw(wd) << iiter;
		for (int i = 0; i < 4; i++) {
			double appo = 0.;
			for (int ithrow = 0; ithrow < nthrows[i]; ithrow++)
				appo += generator.Exponential(1.);
			avg[i] = appo / nthrows[i];
			Output << setw(wd) << avg[i];
		}
		Output << endl;
	}
	Output.close();
}

// Lorentzian dice
void Lorentz() {
	ofstream Output;
	Output.open("output_lorentzian.dat");
	for (int iiter = 0; iiter < niter; iiter++) {
		vector<double> avg(4, 0);
		Output << setw(wd) << iiter;
		for (int i = 0; i < 4; i++) {
			double appo = 0.;
			for (int ithrow = 0; ithrow < nthrows[i]; ithrow++)
				appo += generator.Lorentzian(0., 1.);
			avg[i] = appo / nthrows[i];
			Output << setw(wd) << avg[i];
		}
		Output << endl;
	}
	Output.close();
}