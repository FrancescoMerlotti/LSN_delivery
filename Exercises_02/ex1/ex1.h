# ifndef __ex1_h__
#define __ex1_h__

using namespace std;

// Randon Number Generator Setup
int seed[4];
int p1 = 0, p2 = 0;
Random generator;

// Data blocking
const int nblk = 1e2, nsteps = 1e3, nthrows = 1e4;
double est_int, blk_avg, glob_avg, glob_avg2;
int blk_norm;
int flag;

// Functions
void Input(int);
void Reset(int);
void Measure(int);
void Accumulate();
void Average(int, int);
double StatUncertainty(double, double, int);
double Integrand(double);
double Pdf(double);
double Cumulative(double);
double PdfAR(double);

#endif