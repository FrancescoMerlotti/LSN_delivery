#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include "../Random/random.h"
#include "ex2.h"

using namespace std;

int main(int argc, char** argv) {
	// Arguments check
	if (argc != 2) {
		cout << endl << "Usage: ./ex2 seed" << endl << endl;
		return -1;
	}
	// Setup
	Input(atoi(argv[1]));
	// Random walk
	for (int istep = 0; istep < nsteps; istep++) {
		for (int iblk = 0; iblk < nblk; iblk++) {
			Reset(iblk);
			for (int i = 0; i < len; i++) {
				Measure(i + iblk * len);
				Accumulate();	
			}
			Average(iblk);
		}
		Print(istep, flag);
		Move(flag);
	}

	return 0;
}

// Input seutp
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
	// Setup loading
	Setup.open("setup.in");
	Setup >> flag;
	Setup.close();
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

// Measures
void Measure(int index) {
	est_pos = MagnitudeSquared(walker[index]);
}

// Accumulate
void Accumulate() {
	blk_avg += est_pos;
	blk_norm++;
}

// Averages
void Average(int iblk) {
	blk_avg /= blk_norm;
	glob_avg += blk_avg;
	glob_avg2 +=  blk_avg * blk_avg;	
}

// Move
void Move(int flag) {
	// Random walk in a cubic lattice
	if(flag == 0) {
		for(int ithrow = 0; ithrow < nthrows; ithrow++) {
			int dir = int(generator.Rannyu(0., 6.));
			if(dir % 2 == 0)
				walker[ithrow][dir / 2] += a;
			if(dir % 2 == 1)
				walker[ithrow][dir / 2] -= a;
		}
	}
	// Random walk in a 
	if(flag == 1) {
		for(int ithrow = 0; ithrow < nthrows; ithrow++) {
			double theta = generator.Rannyu(0., M_PI);
			double phi = generator.Rannyu(0., 2 * M_PI);
			walker[ithrow][0] += a * sin(theta) * cos(phi);
			walker[ithrow][1] += a * sin(theta) * sin(phi);
			walker[ithrow][2] += a * cos(theta);;
		}
	}
}

// Print
void Print(int istep, int flag) {
	ofstream Output;
	const int wd = 15;
	vector<string> filenames {"lattice", "continuum"};
	if(istep == 0)
		Output.open("output_" + filenames[flag] + ".dat");
	else
		Output.open("output_" + filenames[flag] + ".dat", ios::app);
	Output << setw(wd) << istep << setw(wd) <<  sqrt(glob_avg / nblk) << setw(wd) << StatUncertainty(glob_avg, glob_avg2, nblk - 1) << endl;
}

// Error
double StatUncertainty(double sum, double sum2, int iblk) {
   	return sqrt(fabs(sum2 / double(iblk + 1) - pow(sum / double(iblk + 1), 2)) / double(iblk + 1));
}

// Magnitude Squared
double MagnitudeSquared(vector<double>& position) {
	if (position.size() != 3)
		return -1;
	return pow(position[0], 2) + pow(position[1], 2) + pow(position[2], 2);
}