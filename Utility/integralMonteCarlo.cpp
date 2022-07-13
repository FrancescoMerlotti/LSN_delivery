#include <iostream>
#include <fstream>
#include <cmath>
#include "../../Random/random.h"

using namespace std;

double integrand(double x) {
	return (M_PI / 2.) * cos(M_PI * x / 2.);
}

double pdf(double x) {
	return 1. - M_PI * (x - 0.5) / 2.;
}

double pdf2(double x) {
	return 1.5 * (1. - pow(x, 2));
}

int main(int argc, char** argv) {

	if (argc != 2) {
		cout << endl << "Usage: ./ex1 seed" << endl << endl;
		return -1;
	}

	// Randon Number Generator Setup
	Random generator;
	int seed[4]{ 0, 0, 0, 0 };
	int p1 = 0, p2 = 0;
	ifstream primes("../../Random/primes32001.in");
	ifstream input("../../Random/seed.in");
	string property;

	if (primes.is_open())
		for (int i = 0; i < atoi(argv[1]); i++)
			primes >> p1 >> p2;
	else
		cerr << "PROBLEM: Unable to open Primes" << endl;
	primes.close();

	if (input.is_open()) {
		while (!input.eof()) {
			input >> property;
			if (property == "RANDOMSEED") {
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				generator.SetRandom(seed, p1, p2);
			}
		}
	}
	else
		cerr << "PROBLEM: Unable to open seed.in" << endl;
	input.close();

	//----------------------------------------------------

	double app1 = 0., app2 = 0., app3 = 0., app4 = 0.;
	const int nthrows = 1e9;

	for (int i = 0; i < nthrows; i++) {
		double x1 = generator.Rannyu();
		// double x2 = generator.Rannyu();
		app1 += integrand(x1);
		app2 += pow(integrand(x1), 2);
		app3 += pow(integrand(x1), 2) / pdf(x1);
		app4 += pow(integrand(x1), 2) / pdf2(x1);
	}

	cout << endl << "Integral = " << app1 / nthrows << endl;
	cout << endl << "Sigma_UN = " << sqrt(app2 / nthrows - pow(app1 / nthrows, 2)) << endl;
	cout << endl << "Sigma_IS = " << sqrt(app3 / nthrows - pow(app1 / nthrows, 2)) << endl;
	cout << endl << "Sigma_AR = " << sqrt(app4 / nthrows - pow(app1 / nthrows, 2)) << endl;

	return 0;
}