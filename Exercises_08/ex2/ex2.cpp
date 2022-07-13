#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include "../Random/random.h"
#include "./ex2.h"

using namespace std;

int main(int argc, char** argv){

	if(argc != 2){
		cout << "Right usage: ./ex1.exe <nseed>" << endl << endl;
		return -1;
	}
	// Setup the simulation
	Input(atoi(argv[1]));
	for(iann = 0; iann < 25; iann++) {
		// Fixing the temperature and beta
		temp = index[iann % 5] * pow(10, -iann / 5);
		nstep = nsteps[iann / 5];
		beta = 1. / temp;
		cout << endl << "Temperature = " << temp << ", Beta = " << beta << ", Number of steps = " << nstep << endl;
		accepted_param = attempted_param = 0;
		// Simulation at the fixed temperature
    	for(int istep = 1; istep <= nstep; istep++) {
			// Parameters move
			MoveParam();			
		}
		Print();
		cout << "Acceptance Parameters = " << double(accepted_param) / attempted_param << endl;
	}
	return 0;
}

// Input
void Input(int nseed) {
	// Read seed for random numbers
	ifstream Primes, Seed, Setup;
	Primes.open("../Random/Primes");
	if(Primes.fail()) {
		cout << "Error while opening Primes" << endl;
		return;
	}
	// Primes setup
	int p1, p2;
	for(int is = 0; is < nseed; is++)
    	Primes >> p1 >> p2;
	Primes >> p1 >> p2;
	Primes.close();
	// Seed setup
	Seed.open("../Random/Seed");
	if(Seed.fail()) {
		cout << "Error while opening seed.in" << endl;
		return;
	}
	Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	Seed.close();
	// RNG setup
	rnd.SetRandom(seed, p1, p2);
	Setup.open("setup.dat");
	if(!Setup.good()) {
		cout << "Error while opening setup.dat" << endl;
		return;
	}
	// Data blocking setup
	Setup >> nblk;
	cout << endl << "Performing a simulation with " << nblk << " blocks and " << nint << " steps" << endl;
	// Metropolis setup
	Setup >> dx >> dp[im] >> dp[is];
	cout << "dx = " << dx << ", dm = " << dp[im] << ", and ds = " << dp[is] << endl;
	// Psi setup
	Setup >> param_old[im] >> param_old[is];
	// Starting point
	Setup >> x_old;
	Setup.close();
	// Psi properties print
	cout << endl << "Psi parameters:" << endl << "Mu = " << param_old[im] << ", Sigma = " << param_old[is] << endl;
	cout << "The test wave function is exp(-(x-" << param_old[im] << ")^2 / (2*" << param_old[is] << "^2)) + exp(-(x+" << param_old[im] << ")^2 / (2*" << param_old[is] << "^2))" << endl;
	// Estimating initial energy
	for(int ip = 0; ip < npar; ip++)
		param[ip] = param_old[ip];
	energy_old = Energy();
	cout << endl << "Initial energy = " << energy_old << endl;
}

// Reset
void Reset(int iblk) {
	// Reset global average
	if(iblk == 1){
		glob_av  = 0;
       	glob_av2 = 0;
	}
	// Reset the block average and block norm
	blk_av = 0;
	blk_norm = 0;
	// Reset the metropolis acceptance
	attempted = 0;
	accepted = 0;
	return;
}

// Move
void Move() {
	double p, psi2_old, psi2_new;
	x = x_old;
	// Old square module of psi
	psi2_old = pow(Psi(x), 2);
	// Move
	x += dx * rnd.Rannyu(-0.5, 0.5);
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

// Move for parameters
void MoveParam() {
	double p, energy_new, bolztmann_new, bolztmann_old;
	for(int ip = 0; ip < npar; ip++)
		param[ip] = param_old[ip];
	// Old Boltzmann weight
	bolztmann_old = Boltzmann(energy_old);
	// Move
	for(int ip = 0; ip < npar; ip++)
		param[ip] += dp[ip] * rnd.Rannyu(-0.5, 0.5);
	// New energy
	energy_new = Energy();
	// New Boltzmann weight
	bolztmann_new = Boltzmann(energy_new);
	// Metropolis test
	p = bolztmann_new / bolztmann_old;
	attempted_param += 1;
	if(rnd.Rannyu() <= p) {
		for(int ip = 0; ip < npar; ip++)
			param_old[ip] = param[ip];
		accepted_param += 1;
		energy_old = energy_new;
	}
	else {
		for(int ip = 0; ip < npar; ip++)
			param[ip] = param_old[ip];
		energy_old = Energy();
	}
}

// Accumulate
void Accumulate() {
	blk_av += Eloc(x);
	blk_norm += 1.;
}

// Averages
void Averages(int iblk) {
	est_energy = blk_av / blk_norm;
	glob_av += est_energy;
	glob_av2 += est_energy * est_energy;
	err_energy = Error(glob_av, glob_av2, iblk);
}

// Print on file
void Print() {
	ofstream Energy, Mu, Sigma;
	const int wd = 15;
	// Selection the type of file writing
	if(iann == 0) {
		Energy.open("./output/output_energy.dat");
		Mu.open("./output/output_mu.dat");
		Sigma.open("./output/output_sigma.dat");
	}
	else {
		Energy.open("./output/output_energy.dat", ios::app);
		Mu.open("./output/output_mu.dat", ios::app);
		Sigma.open("./output/output_sigma.dat", ios::app);
	}
	// Output on file
	Energy << setw(wd) << iann << setw(wd) << temp << setw(wd) << glob_av / double(nblk) << setw(wd) << err_energy << endl;
	Mu << setw(wd) << iann << setw(wd) << temp << setw(wd) << param[0] << endl;
	Sigma << setw(wd) << iann << setw(wd) << temp << setw(wd) << param[1] << endl;
	// Closing streams
	Energy.close();	
	Mu.close();
	Sigma.close();
}

// Error
double Error(double sum, double sum2, int iblk) {
    return sqrt(fabs(sum2 / double(iblk) - pow(sum / double(iblk), 2)) / double(iblk));
}

// Wave function
double Psi(double x) {
	double sigma2 = param[is] * param[is];
	return exp(-pow(x - param[im], 2) / (2. * sigma2)) + exp(-pow(x + param[im], 2) / (2. * sigma2));
}

// Potential
double V(double x) {
	return pow(x, 4) - 2.5 * pow(x, 2);
	// Harmonic oscillator 
	// return 0.5 * pow(x, 2);
}

// Local energy
double Eloc(double x) {
	double sigma2 = param[is] * param[is];
	double mu1 = (x - param[im]) * (x - param[im]), mu2 = (x + param[im]) * (x + param[im]);
	return 0.5 * (1. - (mu1 * exp(- mu1 / (2. * sigma2)) + mu2 * exp(- mu2 / (2 * sigma2))) / (Psi(x) * sigma2) ) / sigma2 + V(x);
	// Harmonic oscillator 
	// return 0.5 * (1. - mu1 / sigma2) / sigma2 + V(x);
}

// Boltzmann
double Boltzmann(double energy) {
	return exp(- beta * energy);
}

// Energy computation
double Energy() {
	for(int iblk = 1; iblk <= nblk; iblk++) {
		Reset(iblk);
		for(int iint = 1; iint <= nint; iint++) {
			Move();
			Accumulate();
		}
		Averages(iblk);
	}
	return glob_av / double(nblk);
}