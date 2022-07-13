#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
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
	// Option price computation
	for (int iblk = 0; iblk < nblk; iblk++) {
		Reset(iblk);
		for (int i = 0; i < len; i++) {
			Measure();
			Accumulate();
		}
		Average(iblk);
	}

	/*
	// Computing the price for a European call option with 
	// Block averages
	for (int i = 0; i < nblocks; i++) {
		double appo = 0.;
		for (int j = 0; j < len; j++) {
			double s = s_0;
			for (int l = 0; l < nsteps; l++) {
				s = s * exp((r - 0.5 * pow(volatility, 2)) * dt + volatility * generator.Gauss(0., 1.) * sqrt(dt));
			}
			appo += exp(-r * dtime) * max(0., s - k);
		}
		avg_blocks[i] = appo / len;
		avg2_blocks[i] = pow(avg_blocks[i], 2);
	}

	// Progressive average and quadratic average
	for (int i = 0; i < nblocks; i++) {
		double appo = 0.;
		double app2 = 0.;
		for (int j = 0; j < i + 1; j++) {
			appo += avg_blocks[j];
			app2 += avg2_blocks[j];
		}
		avg_prog[i] = appo / (i + 1);
		avg2_prog[i] = app2 / (i + 1);
		sigma[i] = statUncertenty(avg_prog, avg2_prog, i);
	}

	// Print on file
	output.open("data1_call_disc.dat");
	for (int i = 0; i < nblocks; i++) {
		output << avg_prog[i] << " " << sigma[i] << endl;
	}
	output.close();

	//----------------------------------------------------

	// Vectors reinitialization
	for (int i = 0; i < nblocks; i++) {
		avg_blocks[i] = 0.;
		avg_prog[i] = 0.;
		avg2_blocks[i] = 0.;
		avg2_prog[i] = 0.;
		sigma[i] = 0.;
	}
	
	//----------------------------------------------------

	// Computing directly the price for a European put option
	// Block averages
	for (int i = 0; i < nblocks; i++) {
		double appo = 0.;
		for (int j = 0; j < len; j++) {
			appo += exp(-r * dtime) * max(0., k - s_0 * exp((r - 0.5 * pow(volatility, 2)) * dtime + volatility * generator.Gauss(0., dtime)));
		}
		avg_blocks[i] = appo / len;
		avg2_blocks[i] = pow(avg_blocks[i], 2);
	}

	// Progressive average and quadratic average
	for (int i = 0; i < nblocks; i++) {
		double appo = 0.;
		double app2 = 0.;
		for (int j = 0; j < i + 1; j++) {
			appo += avg_blocks[j];
			app2 += avg2_blocks[j];
		}
		avg_prog[i] = appo / (i + 1);
		avg2_prog[i] = app2 / (i + 1);
		sigma[i] = statUncertenty(avg_prog, avg2_prog, i);
	}

	// Print on file
	output.open("data1_put_dir.dat");
	for (int i = 0; i < nblocks; i++) {
		output << avg_prog[i] << " " << sigma[i] << endl;
	}
	output.close();

	//----------------------------------------------------

	// Vectors reinitialization
	for (int i = 0; i < nblocks; i++) {
		avg_blocks[i] = 0.;
		avg_prog[i] = 0.;
		avg2_blocks[i] = 0.;
		avg2_prog[i] = 0.;
		sigma[i] = 0.;
	}

	//----------------------------------------------------

	// Computing the price for a European put option with 
	// Block averages
	for (int i = 0; i < nblocks; i++) {
		double appo = 0.;
		for (int j = 0; j < len; j++) {
			double s = s_0;
			for (int l = 0; l < nsteps; l++) {
				s = s * exp((r - 0.5 * pow(volatility, 2)) * dt + volatility * generator.Gauss(0., 1.) * sqrt(dt));
			}
			appo += exp(-r * dtime) * max(0., k - s);
		}
		avg_blocks[i] = appo / len;
		avg2_blocks[i] = pow(avg_blocks[i], 2);
	}

	// Progressive average and quadratic average
	for (int i = 0; i < nblocks; i++) {
		double appo = 0.;
		double app2 = 0.;
		for (int j = 0; j < i + 1; j++) {
			appo += avg_blocks[j];
			app2 += avg2_blocks[j];
		}
		avg_prog[i] = appo / (i + 1);
		avg2_prog[i] = app2 / (i + 1);
		sigma[i] = statUncertenty(avg_prog, avg2_prog, i);
	}

	// Print on file
	output.open("data1_put_disc.dat");
	for (int i = 0; i < nblocks; i++) {
		output << avg_prog[i] << " " << sigma[i] << endl;
	}
	output.close();
	*/
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
	// Setup
	Setup.open("setup.in");
	Setup >> f1 >> f2 >> nblk >> nsteps >> nthrows;
	Setup.close();
	len = nthrows / nblk;
	dt = dtime / nsteps;
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
void Measure() {
	// Call option
	if(f1 == 0) {
		// Directly estimated
		if(f2 == 0) {
			est_opt = exp(-r * dtime) * max(0., s_0 * exp((r - 0.5 * pow(volatility, 2)) * dtime + volatility * generator.Gauss(0., dtime)) - k);
		}
		// Estimated with progressive sampling
		if(f2 == 1) {
			double s = s_0;
			for (int istep = 0; istep < nsteps; istep++) {
				s *= exp((r - 0.5 * pow(volatility, 2)) * dt + volatility * generator.Gauss(0., 1.) * sqrt(dt));
			}
			est_opt =  exp(-r * dtime) * max(0., s - k);
		}
	}
	// Put option
	if(f1 == 1) {
		// Directly estimated
		if(f2 == 0) {
			est_opt = exp(-r * dtime) * max(0., k - s_0 * exp((r - 0.5 * pow(volatility, 2)) * dtime + volatility * generator.Gauss(0., dtime)));
		}
		// Estimated with progressive sampling
		if(f2 == 1) {
			double s = s_0;
			for (int istep = 0; istep < nsteps; istep++) {
				s *= exp((r - 0.5 * pow(volatility, 2)) * dt + volatility * generator.Gauss(0., 1.) * sqrt(dt));
			}
			est_opt =  exp(-r * dtime) * max(0., k - s);
		}
	}
}

// Accumulate
void Accumulate() {
	blk_avg += est_opt;
	blk_norm++;
}

// Averages
void Average(int iblk) {
	int wd = 15;
	ofstream Output;
	vector<string> option {"call", "put"};
	vector<string> method {"dir", "smpl"};
	if(iblk == 0)
		Output.open("output_" + option[f1] + "_" + method[f2] + ".dat");
	else
		Output.open("output_" + option[f1] + "_" + method[f2] + ".dat", ios::app);
	blk_avg /= blk_norm;
	glob_avg += blk_avg;
	glob_avg2 +=  blk_avg * blk_avg;
	Output << setw(wd) << iblk << setw(wd) << blk_avg << setw(wd) << glob_avg / (iblk + 1) << setw(wd) << StatUncertainty(glob_avg, glob_avg2, iblk) << endl;
}

// Error
double StatUncertainty(double sum, double sum2, int iblk) {
   	return sqrt(fabs(sum2 / double(iblk + 1) - pow(sum / double(iblk + 1), 2)) / double(iblk + 1));
}