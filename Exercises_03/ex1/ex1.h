# ifndef __ex1_h__
#define __ex1_h__

using namespace std;

// Randon Number Generator Setup
int seed[4];
int p1 = 0, p2 = 0;
Random generator;

// Data blocking
int nblk, nsteps, nthrows, len;
double est_opt, blk_avg, glob_avg, glob_avg2;
int blk_norm;
int f1, f2;

// Parameters
const double dtime = 1., k = 100., r = 0.1, volatility = 0.25, s_0 = 100.;
double dt;

// Functions
void Input(int);
void Reset(int);
void Measure();
void Accumulate();
void Average(int);
double StatUncertainty(double, double, int);

#endif