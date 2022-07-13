# ifndef __ex1_h__
#define __ex1_h__

using namespace std;

// Randon Number Generator Setup
int seed[4];
int p1 = 0, p2 = 0;
Random generator;

// Data blocking
const int nblk = 1e2, nsteps = 1e4, nint = 1e2;
vector<double> data;
vector<int> cont(100, 0);
double est_avg, est_var, est_chi, blk_avg, blk_var, glob_avg, glob_avg2, glob_var, glob_var2;
int blk_norm;

// Functions
void Input(int);
void Reset(int);
void Measure(int);
void Average(int);
double StatUncertainty(double, double, int);

#endif