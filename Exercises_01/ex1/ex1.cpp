#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include "../Random/random.h"
#include "ex1.h"

using namespace std;

int main(int argc, char** argv) {

	if (argc != 2) {
		cout << endl << "Correct usage: ./ex1 nseed" << endl << endl;
		return -1;
	}
	// Setup
	Input(atoi(argv[1]));
	// Data generation
	for (int idat = 0; idat < nblk * nsteps; idat++) {
		data.push_back(generator.Rannyu());
	}
	// Block average (measures production)
	for (int iblk = 0; iblk < nblk; iblk++) {
		Reset(iblk);
		Measure(iblk);
		Average(iblk);
	}

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

// Reset
void Reset(int iblk) {
	if(iblk == 0) {
		glob_avg = 0.;
		glob_avg2 = 0.;
		glob_var = 0.;
		glob_var2 = 0.;
	}
	blk_avg = 0.;
	blk_var = 0.;
	vector<int> reset(100, 0);
	cont = reset;
	blk_norm = 0;
}

// Measures
void Measure(int iblk) {
	for(int istep = 0; istep < nsteps; istep++) {
		blk_avg += data[istep + iblk * nsteps];
		blk_var += pow(data[istep + iblk * nsteps] - 0.5, 2);
		cont[int(nint * data[istep + iblk * nsteps])]++;
		blk_norm++;
	}
}

// Averages
void Average(int iblk) {
	ofstream Mean, Variance, Chi;
	int wd = 12;
	if(iblk == 0) {
		Mean.open("output_mean.dat");
		Variance.open("output_variance.dat");
		Chi.open("output_chi.dat");
	}
	else {
		Mean.open("output_mean.dat", ios::app);
		Variance.open("output_variance.dat", ios::app);
		Chi.open("output_chi.dat", ios::app);
	}

	est_avg = blk_avg / blk_norm;
	glob_avg += est_avg;
	glob_avg2 +=  est_avg * est_avg;
	Mean << setw(wd) << iblk << setw(wd) << est_avg << setw(wd) << glob_avg / (iblk + 1) << setw(wd) << StatUncertainty(glob_avg, glob_avg2, iblk) << endl;	

	est_var = blk_var / blk_norm;
	glob_var += est_var;
	glob_var2 += est_var * est_var;
	Variance << setw(wd) << iblk << setw(wd) << est_var << setw(wd) << glob_var / (iblk + 1) << setw(wd) << StatUncertainty(glob_var, glob_var2, iblk) << endl;

	for(int iint = 0; iint < nint; iint++)
		est_chi += pow(cont[iint] - nsteps / nint, 2);
	est_chi /= (nsteps / nint);
	Chi << setw(wd) << iblk << setw(wd) << est_chi << endl;
}

// Error
double StatUncertainty(double sum, double sum2, int iblk) {
   	return sqrt(fabs(sum2 / double(iblk + 1) - pow(sum / double(iblk + 1), 2)) / double(iblk + 1));
}
