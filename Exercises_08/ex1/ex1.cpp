#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "ex1.h"

using namespace std;

int main(int argc, char** argv){
	if(argc != 2){
		cout << "Right usage: ./ex1.exe <nseed>" << endl << endl;
		return -1;
	}
	Input(atoi(argv[1]));
	for(int iblk = 1; iblk <= nblk; iblk++) {
    	Reset(iblk);
    	for(int i = 1; i <= nstep; i++) {
			Move();
			Accumulate();
			Print(iblk, i);
		}
		Averages(iblk);
	}
	return 0;
}

// Input
void Input(int nseed) {
	// Read seed for random numbers
	ifstream Primes, Seed, Setup;
	Primes.open("../Random/Primes");
	int p1, p2;
	for(int i = 0; i < nseed; i++)
    Primes >> p1 >> p2 ;
	Primes.close();
	Seed.open("../Random/Seed");
	Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	Seed.close();
	// RNG setup
	rnd.SetRandom(seed, p1, p2);
	Setup.open("setup.dat");
	// Data blocking setup
	Setup >> nblk >> nstep;
	// Metropolis setup
	Setup >> delta;
	// Psi setup
	Setup >> mu >> sigma;
	// Starting point
	Setup >> x;
	Setup.close();
	x_old = x;
	cout << endl << "Psi parameters:" << endl << "Mu = " << mu << ", Sigma = " << sigma << endl;
}

// Reset
void Reset(int iblk) {
   // Reset the global average
	if(iblk == 1) {
        glob_av  = 0;
        glob_av2 = 0;
	}
	// reset the block average and block norm
	blk_avg = 0;
	blk_norm = 0;
	// reset the metropolis acceptance
	attempted = 0;
	accepted = 0;
}

// Move
void Move() {
	double p, psi2_old, psi2_new;
	x = x_old;
	// Old square module of psi
	psi2_old = pow(Psi(x), 2);
	// Move
	x += delta * rnd.Rannyu(-1., 1.);
    // New square module of psi
	psi2_new = pow(Psi(x), 2);
	//Metropolis test
	p = psi2_new / psi2_old;
	if(rnd.Rannyu() <= p) {
		// Update
		x_old = x;
		accepted += 1;
	}
	else {
		// Reweight
		x = x_old;
	}
	attempted += 1;
}

// Accumulate
void Accumulate() {
	blk_avg += Eloc(x);
	blk_norm += 1.;
}

// Averages
void Averages(int iblk) {
	ofstream Energy;
	const int wd = 15;

	est_energy = blk_avg / blk_norm;
	glob_av += est_energy;
	glob_av2 += est_energy * est_energy;
	err_energy = Error(glob_av, glob_av2, iblk);

	if(iblk == 1)
		Energy.open("output_energy_best.dat");
	else
		Energy.open("output_energy_best.dat", ios::app);
	Energy << setw(wd) << iblk << setw(wd) << glob_av / double(iblk) << setw(wd) << err_energy << endl;
	Energy.close();

	cout << endl << "Block Number = " << iblk << endl;
	cout << "Acceptance Rate = " << double(accepted) / attempted << endl;
}

// Print for Histogram
void Print(int i, int j) {
	ofstream Hist;
	if(i == 1 and j == 1)
		Hist.open("output_hist.dat");
	else
		Hist.open("output_hist.dat", ios::app);
	Hist << setw(10) << x_old << endl;
	Hist.close();
}

// Error
double Error(double sum, double sum2, int iblk) {
    return sqrt(fabs(sum2 / double(iblk) - pow(sum / double(iblk), 2)) / double(iblk));
}

// Wave function
double Psi(double x) {
	double sigma2 = sigma * sigma;
	return exp(-pow(x - mu, 2) / (2. * sigma2)) + exp(-pow(x + mu, 2) / (2. * sigma2));
}

// Potential
double V(double x) {
	return pow(x, 4) - 2.5 * pow(x, 2);
	// Harmonic oscillator 
	// return 0.5 * pow(x, 2);
}

// Local energy
double Eloc(double x) {
	double sigma2 = sigma * sigma;
	double mu1 = (x - mu) * (x - mu), mu2 = (x + mu) * (x + mu);
	return 0.5 * (1. - (mu1 * exp(- mu1 / (2. * sigma2)) + mu2 * exp(- mu2 / (2 * sigma2))) / (Psi(x) * sigma2) ) / sigma2 + V(x);
	// Harmonic oscillator 
	// return 0.5 * (1. - mu1 / sigma2) / sigma2 + V(x);
}