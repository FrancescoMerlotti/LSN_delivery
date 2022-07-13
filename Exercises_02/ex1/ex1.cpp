#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include "../Random/random.h"
#include "ex1.h"

using namespace std;

int main(int argc, char** argv) {
	// Arguments check
	if (argc != 2) {
		cout << endl << "Usage: ./ex1 seed" << endl << endl;
		return -1;
	}
	// Setup
	Input(atoi(argv[1]));
	for(int flag = 0; flag < 3; flag++) {
		for(int iblk = 0; iblk < nblk; iblk++) {
		Reset(iblk);
		for(int istep = 0; istep < nsteps; istep++) {
			Measure(flag);
			Accumulate();
		}		
		Average(iblk, flag);
		}
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
	}
	blk_avg = 0.;
	blk_norm = 0;
}

// Measures
void Measure(int flag) {
	est_int = 0.;
	// Uniform sampling
	if(flag == 0) {
		for(int ithrow = 0; ithrow < nthrows; ithrow++)
			est_int += Integrand(generator.Rannyu());
	}
	// Importance sampling
	if(flag == 1) {
		for(int ithrow = 0; ithrow < nthrows; ithrow++) {
			double x = Cumulative(generator.Rannyu());
			est_int += Integrand(x) / Pdf(x);
		}			
	}
	// Importance sampling & accept-reject
	if(flag == 2) {
		for(int ithrow = 0; ithrow < nthrows; ithrow++) {
			double x = generator.Rannyu();
			double y = generator.Rannyu(0., 1.5);
			while (y > PdfAR(x)) {
				x = generator.Rannyu();
				y = generator.Rannyu(0., 1.5);
			}				
			est_int += Integrand(x) / PdfAR(x);
		}			
	}
	est_int /= nthrows;
}

// Accumulate
void Accumulate() {
	blk_avg += est_int;
	blk_norm++;
}

// Averages
void Average(int iblk, int flag) {
	int wd = 15;
	ofstream Output;
	vector<string> filenames {"standard", "importance", "acc_rej"};
	if(iblk == 0)
		Output.open("output_" + filenames[flag] + ".dat");
	else
		Output.open("output_" + filenames[flag] + ".dat", ios::app);
	blk_avg /= blk_norm;
	glob_avg += blk_avg;
	glob_avg2 +=  blk_avg * blk_avg;
	Output << setw(wd) << iblk << setw(wd) << blk_avg << setw(wd) << glob_avg / (iblk + 1) << setw(wd) << StatUncertainty(glob_avg, glob_avg2, iblk) << endl;	
}

// Error
double StatUncertainty(double sum, double sum2, int iblk) {
   	return sqrt(fabs(sum2 / double(iblk + 1) - pow(sum / double(iblk + 1), 2)) / double(iblk + 1));
}

// Integrand
double Integrand(double x) {
	return (M_PI / 2.) * cos(M_PI * x / 2.);
}

// Probability distribution function
double Pdf(double x) {
	return 1. - M_PI * (x - 0.5) / 2.;
}

// Cumulative function
double Cumulative(double x) {
	double a = M_PI / 4.;
	double b = (1 + M_PI / 4.);
	return (b - sqrt(pow(b, 2) - 4. * a * x)) / (2. * a);
}

// Probability distribution function for accept-reject
double PdfAR(double x) {
	return 1.5 * (1. - pow(x, 2));
}