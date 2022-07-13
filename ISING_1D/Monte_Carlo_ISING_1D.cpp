/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main() { 
  for(t = 0.5; t < 2.; t += 0.015) {
    Input(); //Inizialization
    for(int iblk = 1; iblk <= nblk; ++iblk) { // Simulation
      Reset(iblk);   //Reset block averages
      for(int istep=1; istep <= nstep; ++istep) {
        Move(metro);
        Measure();
        Accumulate(); //Update block averages
      }
      Averages(iblk);   //Print results for current block
    }
  }
  // ConfFinal(); //Write final configuration
  return 0;
}

void Input() {
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

  //Read seed for random numbers
  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed, p1, p2);
  input.close();
  
  //Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp; // during this simulation it will use the one provided by the main
  beta = 1. / t;
  cout << "Temperature = " << t << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  if(metro==1)
    cout << "The program perform Metropolis moves" << endl;
  else
    cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


  //Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

  //initial configuration
  for (int i=0; i<nspin; ++i) {
    if(rnd.Rannyu() <= 0.5)
      s[i] = 1;
    else
      s[i] = -1;
  }
  
  //Evaluate energy etc. of the initial configuration
  Measure();

  //Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu] / double(nspin) << endl;
}

void Move(int metro) {
  int o, sm;
  double p, energy_old, energy_new;

  for(int i = 0; i < nspin; ++i) {
    //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = int(rnd.Rannyu() * nspin);
    
    if(metro==1) {
      sm = - s[o];
      energy_old = Boltzmann(s[o], o);
      energy_new = Boltzmann(sm, o);
      p = 1. / (1. + exp(beta * (energy_new - energy_old)));
      if(rnd.Rannyu() <= p) {
        accepted++;
        s[o] = sm;
      }
      //else { }
      attempted++;
    }
    else {
      // suggesting a move following the probabilty distribution to sample
      p = 1. / (1. + exp(-2. * beta * ((s[Pbc(o-1)] + s[Pbc(o+1)]) + h)));
      if(rnd.Rannyu() <= p) {
        s[o] = 1;
      }
      else {
        s[o] = -1;
      }
      accepted++;
      attempted++;
    }
  }
}

double Boltzmann(int sm, int ip) {
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure() {
  double u = 0., m = 0.;
  //cycle over spins
  for (int i = 0; i < nspin; i++) {
    u += -1. * J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
    m += s[i];
  }
  walker[iu] = u;
  walker[im] = m;
}

void Reset(int iblk) {
   if(iblk == 1) {
      for(int i=0; i < n_props; ++i) {
          glob_av[i] = 0;
          glob_av2[i] = 0;
       }
   }
   for(int i=0; i<n_props; ++i) {
      blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}

void Accumulate(void) {
  blk_av[iu] += walker[iu];
  blk_av[ic] += walker[iu] * walker[iu];
  blk_av[im] += walker[im];
  blk_av[ix] += walker[im] * walker[im];
  blk_norm += 1.;
}

void Averages(int iblk) {  
  ofstream Ene, Heat, Mag, Chi;
  const int wd = 12;
    
  cout << "Block number " << iblk << endl;
  cout << "Attempted moves " << attempted << endl;
  cout << "Acceptance rate " << accepted / attempted << endl << endl;

  stima_u = blk_av[iu] / blk_norm;
  stima_u2 = blk_av[ic] / blk_norm;
  stima_m = blk_av[im] / blk_norm;
  stima_m2 = blk_av[ix] / blk_norm;
  stima_c = beta * beta * (stima_u2 - stima_u * stima_u) / double(nspin);
  stima_x = beta * stima_m2 / double(nspin);
  stima_u /= double(nspin);
  stima_m /= double(nspin);

  glob_av[iu]  += stima_u;
  glob_av2[iu] += pow(stima_u, 2);
  err_u = Error(glob_av[iu], glob_av2[iu], iblk);

  glob_av[im]  += stima_m;
  glob_av2[im] += pow(stima_m, 2);
  err_m = Error(glob_av[im], glob_av2[im], iblk);

  glob_av[ic]  += stima_c;
  glob_av2[ic] += pow(stima_c, 2);
  err_c = Error(glob_av[ic], glob_av2[ic], iblk);

  glob_av[ix]  += stima_x;
  glob_av2[ix] += pow(stima_x, 2);
  err_x = Error(glob_av[ix], glob_av2[ix], iblk);

  if(iblk == 20){
    if(t == 0.5) {
      Ene.open("output_ene.dat");
      Mag.open("output_mag.dat");
      Heat.open("output_heat.dat");
      Chi.open("output_chi.dat");
    }
    else {
      Ene.open("output_ene.dat", ios::app);
      Mag.open("output_mag.dat", ios::app);
      Heat.open("output_heat.dat", ios::app);
      Chi.open("output_chi.dat", ios::app);
    }

    Ene << setw(wd) << t << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/double(iblk) << setw(wd) << err_u << endl;
    Mag << setw(wd) << t << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/double(iblk) << setw(wd) << err_m << endl;
    Heat << setw(wd) << t << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/double(iblk) << setw(wd) << err_c << endl;
    Chi << setw(wd) << t << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/double(iblk) << setw(wd) << err_x << endl;

    Ene.close();
    Mag.close();
    Heat.close();
    Chi.close();
  }

  cout << "----------------------------" << endl << endl;
}

void ConfFinal(void) {
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i) {
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk) {
    if(iblk==1)
      return 0.;
    else
      return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
