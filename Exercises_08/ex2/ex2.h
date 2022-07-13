#ifndef __ex2_h__
#define __ex2_h__

#include <vector>

// Random
int seed[4];
Random rnd;

// Psi parameters
const int npar = 2;
int im = 0, is = 1;
double param[npar], param_old[npar];
double x, x_old;

// Metropolis parameter
double dx;
double dp[npar];
int attempted, accepted;
int attempted_param, accepted_param;
double beta, temp;
int iann;

// Data blocking parameters
int nblk, nstep;
int nsteps[5] {500, 1000, 1500, 2000, 3000}; 
double index[5] {10., 8., 6., 4., 2.};
int nint = 5000;
double energy_old = 0.;
double est_energy, err_energy;
double blk_av, blk_norm, glob_av, glob_av2;

void Input(int);
void Reset(int);
void Move();
void MoveParam();
void Accumulate();
void Averages(int);
void Print();
double Error(double, double, int);
double Psi(double);
double V(double);
double Eloc(double);
double Boltzmann(double);
double Energy();

#endif