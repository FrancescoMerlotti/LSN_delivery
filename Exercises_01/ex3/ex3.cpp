#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include "../Random/random.h"
#include "ex3.h"

using namespace std;

int main(int argc, char** argv) {
	// Arguments check
	if (argc != 2) {
		cout << endl << "Correct usage: ./ex3 seed" << endl << endl;
		return -1;
	}
	// Setup
	Input(atoi(argv[1]));
	
	for (int iblk = 0; iblk < nblk; iblk++) {
		Reset(iblk);
		for(int istep = 0; istep < nsteps; istep++) {
			nhit = 0;
			for (int ithrow = 0; ithrow < nthrows; ithrow++)
				Needle();
			Measure();
		}
		Average(iblk);
	}
	return 0;
}

// Input setup
void Input(int nseed) {
	ifstream Primes, Seed, Setup;
	// Prime loading
	Primes.open("../Random/Primes");
	for(int i = 0; i < nseed; i++)
		Primes >> p1 >> p2;		
	Primes.close();
	// Seed loading
	Seed.open("../Random/Seed");
	Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	Seed.close();
	Setup.open("setup.in");
	Setup >> nblk >> nsteps >> nthrows;
	// Generator setup
	generator.SetRandom(seed, p1, p2);
}

// Reset
void Reset(int iblk) {
	if(iblk == 0) {
		glob_avg = 0.;
		glob_avg2 = 0.;
	}
	blk_avg = 0.;
	blk_norm = 0;
}

// Needle generation
void Needle() {
	// Direction of the needle as an angle
	ux = generator.Rannyu(-1., 1.);
	uy = generator.Rannyu();
	while (pow(ux, 2) + pow(uy, 2) >= 1) {
		ux = generator.Rannyu(-1., 1.);
		uy = generator.Rannyu();
	}
	theta = atan(uy / ux);
	// center of the needle
	xcenter = generator.Rannyu(1., 2.);
	// supposing vertical line on the integers
	if(Hit(xcenter - cos(theta) * len / 2, xcenter + cos(theta) * len / 2))
		nhit++;
}

// Measures
void Measure() {
	blk_avg += constant * double(nthrows) / nhit;
	blk_norm++;
}

// Hit
bool Hit(double x1, double x2) {
	if (int(x1) != int(x2))
		return true;
	if (x1 == int(x1))
		return true;
	if (x2 == int(x2))
		return true;
	else
		return false;
}

// Averages
void Average(int iblk) {
	ofstream Pi;
	int wd = 16;
	if(iblk == 0) {
		Pi.open("output_pi.dat");
	}
	else {
		Pi.open("output_pi.dat", ios::app);
	}

	blk_avg /= blk_norm;
	glob_avg += blk_avg;
	glob_avg2 +=  blk_avg * blk_avg;
	Pi << setw(wd) << iblk << setw(wd) << blk_avg << setw(wd) << glob_avg / (iblk + 1) << setw(wd) << StatUncertainty(glob_avg, glob_avg2, iblk) << endl;

}

// Error
double StatUncertainty(double sum, double sum2, int iblk) {
   	return sqrt(fabs(sum2 / double(iblk + 1) - pow(sum / double(iblk + 1), 2)) / double(iblk + 1));
}
