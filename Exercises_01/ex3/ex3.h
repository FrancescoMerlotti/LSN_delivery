# ifndef __ex3_h__
#define __ex3_h__

using namespace std;

// Randon Number Generator Setup
int seed[4];
int p1 = 0, p2 = 0;
Random generator;

// Data blocking
int nblk, nsteps, nthrows;
double blk_avg, glob_avg, glob_avg2;
int blk_norm, nhit;

// Parameters
const double dist = 1., len = 0.9;
double ux, uy, xcenter, theta;
double constant = 2 * len / dist;
// Functions
void Input(int);
void Reset(int);
void Needle();
void Measure();
bool Hit(double, double);
void Average(int);
double StatUncertainty(double, double, int);

#endif