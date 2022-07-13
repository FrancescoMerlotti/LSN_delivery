#ifndef __ex1_h__
#define __ex1_h__

#include <iostream>
#include <cmath>
#include "../Random/random.h"

// Random
int seed[4];
Random rnd;

// Metropolis parameter
double delta;
int attempted, accepted;

// Data blocking parameters
int nblk, nstep, len;
double est_energy, blk_avg, blk_norm, glob_av, glob_av2, err_energy;

// Psi parameters
double x, x_old, mu, sigma;

void Input(int);
void Reset(int);
void Move();
void Accumulate();
void Averages(int);
void Print(int, int);
double Error(double, double, int);
double Psi(double);
double V(double);
double Eloc(double);



#endif