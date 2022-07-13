# ifndef __ex2_h__
#define __ex2_h__

using namespace std;

// Randon Number Generator Setup
int seed[4];
int p1 = 0, p2 = 0;
Random generator;

// Data blocking
const int niter = 1e4;
int nthrows[4] = {1, 2, 10, 100};
const int wd = 20;

// Functions
void Input(int);
void Standard();
void Exp();
void Lorentz();

#endif