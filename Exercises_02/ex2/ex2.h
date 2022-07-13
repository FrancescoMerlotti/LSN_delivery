# ifndef __ex2_h__
#define __ex2_h__

using namespace std;

// Randon Number Generator Setup
int seed[4];
int p1 = 0, p2 = 0;
Random generator;

// Data blocking
const int nblk = 1e2, nsteps = 1e2, nthrows = 1e4;
int len = nthrows / nblk;
double est_pos, blk_avg, glob_avg, glob_avg2, sigma;
int blk_norm;
int flag;

// Parameters
const double a = 1.;
vector<vector<double>> walker(nthrows, vector<double>(3, 0));

// Functions
void Input(int);
void Reset(int);
void Measure(int);
void Accumulate();
void Average(int);
double StatUncertainty(double, double, int);
void Move(int);
void Print(int, int);
double MagnitudeSquared(vector<double>&);

#endif