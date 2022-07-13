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
#include <string>
#include <cmath>
#include <iomanip>
#include "MD_MC.h"

using namespace std;

int main(int argc, char **argv) {
  ofstream out;
  if(argc != 2){
    cout << "Right usage: ./NVE_NVT.exe <seed>" << endl;
    return -1;
  }
  //Inizialization
  Input(atoi(argv[1]));
  //Equilibration
  for(int iequ = 1; iequ <= nequ; iequ++)
    Move();
  cout << "Equilibration completed successfully" << endl << endl;
  for(int iblk = 1; iblk <= nblk; iblk++) {
    //Reset block averages
    Reset(iblk); 
    for(int istep = 1; istep <= nstep; istep++) {
      Move();
      Measure();
      Accumulate();
    }
    //Print results for current block
    Averages(iblk);
  }
  //Write final configuration
  ConfFinal();
  return 0;
}

// Input
void Input(int nseed) {
  ifstream ReadInput, ReadConf, ReadVelocity, Primes, Seed;
  // Simulation header
  cout << "Classic Lennard-Jones fluid" << endl;
  cout << "MD(NVE) / MC(NVT) simulation" << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
  cout << "The program uses Lennard-Jones units" << endl;
  // Seed
  int p1, p2;
  Primes.open("./Random/Primes");
  for(int i = 0; i < nseed; i++)
    Primes >> p1 >> p2;
  Primes.close();
  // Input informations
  ReadInput.open("./input/input.in");
  // Type of Simulation (1 MC, 0 MD)
  ReadInput >> iNVET;
  // Restart the simulation from the last configuration
  ReadInput >> restart;
  if(restart)
    Seed.open("./Random/seed.out");
  else
    Seed.open("./Random/Seed");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  Seed.close();
  // Random number generator setup
  rnd.SetRandom(seed, p1, p2);
  // Temperature and estimating beta
  ReadInput >> temp;
  beta = 1. / temp;
  cout << "Temperature = " << temp << endl;
  // Number of particles
  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;
  // Numerical density
  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  // Kind of simulation
  if(iNVET) {
    simulation = "MC";
  }
  else{
    simulation = "MD";
  }
  // Phase
  if(rho == 1.1) {
    phase = "solid";
  }
  if(rho == 0.8) {
    phase = "liquid";
  }
  if(rho == 0.05) {
    phase = "gas"; 
  }
  // Box volume
  vol = double(npart) / rho;
  // Box side
  box = pow(vol, 1./3.);
  // max_len = sqrt(3.) * box / 2.;
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;
  // Cutoff distance
  ReadInput >> rcut;
  cout << "Cutoff of the interatomic potential = " << rcut << endl;
  rcut3 = pow(rcut, 3),
  rcut9 = pow(rcut, 9);
  vtail = 8. * M_PI * rho * (1. / (9. * rcut9) - 1. / (3. * rcut3));
  cout << "Potential tail correction = " << vtail << endl;
  ptail = 32. * M_PI * rho * (1. / (9. * rcut9) - 1. / (6. * rcut3)) / vol;
  cout << "Virial tail correction = " << ptail << endl << endl;
  // Time step ("big" for MC (fix the acceptance to 0.5), "small" for MD)
  ReadInput >> delta;
  //Number of blocks
  ReadInput >> nblk;
  //Number of stepsper block
  ReadInput >> nstep;
  // Close Input stream
  ReadInput.close();
  // Metropolis output
  if(iNVET) {
    cout << "The program performs Metropolis moves with uniform translations" << endl;
    cout << "Moves parameter = " << delta << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl << endl;
  }
  // Prepare arrays for measurements
  // Number of observables
  nprops = 5;
  iv = 0; // Potential energy
  it = 1; // Temperature
  ik = 2; // Kinetic energy
  ie = 3; // Total energy
  iw = 4; // Virial
  // Initial configuration
  cout << "Read initial configuration" << endl << endl;
  if(restart) {
    ReadConf.open("./configs/config.out");
    ReadVelocity.open("./configs/velocity.out");
    for (int i = 0; i < npart; ++i)
      ReadVelocity >> vx[i] >> vy[i] >> vz[i];
  }
  else {
    ReadConf.open("./configs/Config");
    cout << "Prepare velocities with center of mass velocity equal to zero " << endl;
    double sumv[3] = {0., 0., 0.};
    for (int i = 0; i < npart; ++i) {
      vx[i] = rnd.Gauss(0.,sqrt(temp));
      vy[i] = rnd.Gauss(0.,sqrt(temp));
      vz[i] = rnd.Gauss(0.,sqrt(temp));
      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
    for (int idim = 0; idim < 3; ++idim)
      sumv[idim] /= double(npart);
    double sumv2 = 0., fs;
    for (int i = 0; i < npart; ++i) {
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];
      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= double(npart);
    // Velocity scale factor 
    fs = sqrt(3 * temp / sumv2);
    cout << "velocity scale factor: " << fs << endl << endl;
    for (int i = 0; i < npart; ++i) {
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;
    }
  }
  // Reaad starting configuration
  for (int i = 0; i < npart; ++i) {
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = Pbc(x[i] * box);
    y[i] = Pbc(y[i] * box);
    z[i] = Pbc(z[i] * box);
  }
  ReadConf.close();

  for (int i = 0; i < npart; ++i) {
    if(iNVET) {
      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];
    }
    else {
      xold[i] = Pbc(x[i] - vx[i] * delta);
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    }
  }
  
  //Evaluate properties of the initial configuration
  Measure();

  //Print initial values for measured properties
  cout << "Initial potential energy = " << walker[iv] / (double)npart << endl;
  cout << "Initial temperature      = " << walker[it] << endl;
  cout << "Initial kinetic energy   = " << walker[ik] / (double)npart << endl;
  cout << "Initial total energy     = " << walker[ie] / (double)npart << endl;

  //Print initial values for measured properties
  cout << "----------------------------" << endl << endl;

  return;
}

// Move
void Move() {
  int o;
  double p, energy_old, energy_new;
  double xnew, ynew, znew;

  if(iNVET) {
    for(int i = 0; i < npart; ++i) {
      // Select randomly a particle
      o = int(rnd.Rannyu() * npart);
      // Old
      energy_old = Boltzmann(x[o], y[o], z[o], o);
      // New
      x[o] = Pbc(x[o] + delta * (rnd.Rannyu() - 0.5));
      y[o] = Pbc(y[o] + delta * (rnd.Rannyu() - 0.5));
      z[o] = Pbc(z[o] + delta * (rnd.Rannyu() - 0.5));
      energy_new = Boltzmann(x[o], y[o], z[o], o);
      // Metropolis test
      p = exp(beta * (energy_old - energy_new));
      if(rnd.Rannyu() <= p) {
        //Update
        xold[o] = x[o];
        yold[o] = y[o];
        zold[o] = z[o];
        accepted += 1.;
      }
      else {
        x[o] = xold[o];
        y[o] = yold[o];
        z[o] = zold[o];
      }
      attempted += 1.;
    }
  }
  else {
    double fx[m_part], fy[m_part], fz[m_part];
    //Force acting on i-particle
    for(int i = 0; i < npart; ++i) { 
      fx[i] = Force(i, 0);
      fy[i] = Force(i, 1);
      fz[i] = Force(i, 2);
    }

    //Verlet integration scheme
    for(int i = 0; i < npart; ++i){ 
      xnew = Pbc(2. * x[i] - xold[i] + fx[i] * pow(delta, 2));
      ynew = Pbc(2. * y[i] - yold[i] + fy[i] * pow(delta, 2));
      znew = Pbc(2. * z[i] - zold[i] + fz[i] * pow(delta, 2));

      vx[i] = Pbc(xnew - xold[i]) / (2. * delta);
      vy[i] = Pbc(ynew - yold[i]) / (2. * delta);
      vz[i] = Pbc(znew - zold[i]) / (2. * delta);

      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];

      x[i] = xnew;
      y[i] = ynew;
      z[i] = znew;

      accepted += 1.;
      attempted += 1.;
    }
  }
  return;
}

// Energy computation
double Boltzmann(double xx, double yy, double zz, int ip) {
  double ene = 0.;
  double dx, dy, dz, dr;
  for (int i=0; i<npart; ++i) {
    if(i != ip) {
      // distance ip-i in pbc
      dx = Pbc(xx - x[i]);
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);
      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);
      if(dr < rcut)
        ene += 1. / pow(dr, 12) - 1. / pow(dr, 6);
    }
  }
  return 4. * ene;
}

// Force computation
double Force(int id, int idir) {
  double f = 0.;
  double dvec[3], dr;

  for (int i = 0; i < npart; ++i) {
    if(i != id){
      dvec[0] = Pbc(x[id] - x[i]);
      dvec[1] = Pbc(y[id] - y[i]);
      dvec[2] = Pbc(z[id] - z[i]);
      dr = sqrt(dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2]);

      // -Grad_id V(r)
      if(dr < rcut) {
        f += dvec[idir] * (48. / pow(dr, 14) - 24. / pow(dr, 8));
      }
    }
  }
  return f;
}

// Measure
void Measure() {
  double v = 0., w = 0., kin = 0.;
  double vij, wij;
  double dx, dy, dz, dr;
  // Reset histogram
  for(int ibin = 0; ibin < nbins; ibin++) {
    histogram[ibin] = 0.;
  }
  // Cycle over pairs of particles
  for (int i = 0; i < npart - 1; ++i) {
    for (int j = i + 1; j < npart; ++j) {
      // Distance i-j in pbc
      dx = Pbc(x[i] - x[j]);
      dy = Pbc(y[i] - y[j]);
      dz = Pbc(z[i] - z[j]);
      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);
      // Radial distribution function
      histogram[int(2 * nbins * dr / box)] += 2.;
      // Test if dr is smaller of the cutoff radius to calculate properties
      if(dr < rcut) {
        double dr12 = pow(dr, 12), dr6 = pow(dr, 6);
        vij = 1. / dr12 - 1. / dr6;
        wij = 1. / dr12 - 0.5 / dr6;
        v += vij;
        w += wij;
      }
    }          
  }
  // This only applies to Molecular Dynamics
  for (int i = 0; i < npart; ++i)
    kin += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  
  walker[iv] = 4. * v;                            // Potential energy
  walker[ik] = kin;                               // Kinetic energy
  walker[it] = 2. * kin / (3. * double(npart));   // Temperature
  walker[ie] = 4. * v + kin;                      // Total energy
  walker[iw] = 48. * w;                           // Virial

  return;
}

// Reset
void Reset(int iblk) { 
  if(iblk == 1) {
    for(int i = 0; i < nprops; ++i) {
      glob_av[i] = 0.;
      glob_av2[i] = 0.;
    }
    for(int i = 0; i < nbins; i++) {
      hist_gl[i] = 0.;
      hist_gl2[i] = 0.;
    }
  }
  for(int i = 0; i < nprops; ++i) {
    blk_av[i] = 0.;
  }
  for(int i = 0; i < nbins; i++) {
    hist_av[i] = 0.;
  }
   
  blk_norm = 0.;
  attempted = 0.;
  accepted = 0.;
}

// Accumulate
void Accumulate(void) {
  for(int i = 0; i < nprops; ++i) {
    blk_av[i] += walker[i];
  }
  // g(r) histogram
  for(int ibin = 0; ibin < nbins; ibin++) {
    hist_av[ibin] += histogram[ibin];
  }
  blk_norm += 1.;
}

// Averages
void Averages(int iblk) {
  ofstream Histogram;
  ofstream Epot, Pres, Ekin, Etot, Temp;
  const int wd = 12;

  cout << "Block number " << iblk << endl;
  cout << "Acceptance rate " << accepted / attempted << endl << endl;

  if(iblk == 1) {
    Epot.open("./output/" + simulation + "/" + phase + "/output_epot.dat");
    Ekin.open("./output/" + simulation + "/" + phase + "/output_ekin.dat");
    Temp.open("./output/" + simulation + "/" + phase + "/output_temp.dat");
    Etot.open("./output/" + simulation + "/" + phase + "/output_etot.dat");
    Pres.open("./output/" + simulation + "/" + phase + "/output_pres.dat");
  }
  else {
    Epot.open("./output/" + simulation + "/" + phase + "/output_epot.dat", ios::app);
    Ekin.open("./output/" + simulation + "/" + phase + "/output_ekin.dat", ios::app);
    Temp.open("./output/" + simulation + "/" + phase + "/output_temp.dat", ios::app);
    Etot.open("./output/" + simulation + "/" + phase + "/output_etot.dat", ios::app);
    Pres.open("./output/" + simulation + "/" + phase + "/output_pres.dat", ios::app);
  }
  
  //Potential energy
  stima_pot = blk_av[iv] / (blk_norm * double(npart)) + vtail;
  glob_av[iv] += stima_pot;
  glob_av2[iv] += stima_pot * stima_pot;
  err_pot = Error(glob_av[iv], glob_av2[iv], iblk);
  
  //Kinetic energy
  stima_kin = blk_av[ik]/(blk_norm * double(npart));
  glob_av[ik] += stima_kin;
  glob_av2[ik] += stima_kin*stima_kin;
  err_kin=Error(glob_av[ik],glob_av2[ik],iblk);

  //Total energy
  stima_etot = blk_av[ie]/(blk_norm * double(npart));
  glob_av[ie] += stima_etot;
  glob_av2[ie] += stima_etot*stima_etot;
  err_etot=Error(glob_av[ie],glob_av2[ie],iblk);
  
  //Temperature
  stima_temp = blk_av[it] / blk_norm;
  glob_av[it] += stima_temp;
  glob_av2[it] += stima_temp * stima_temp;
  err_temp = Error(glob_av[it], glob_av2[it], iblk);
  
  //Pressure
  stima_pres = rho * stima_temp + blk_av[iw] / (3 * vol * blk_norm) + ptail;
  glob_av[iw] += stima_pres;
  glob_av2[iw] += stima_pres * stima_pres;
  err_pres = Error(glob_av[iw], glob_av2[iw], iblk);
  
  //Radial distribution function
  for(int ibin = 0; ibin < nbins; ibin++) {
    stima_gdir[ibin] = hist_av[ibin] / blk_norm;
    hist_gl[ibin] += stima_gdir[ibin];
    hist_gl2[ibin] += stima_gdir[ibin] * stima_gdir[ibin];
    err_gdir[ibin] = Error(hist_gl[ibin], hist_gl2[ibin], iblk);
  }
  
  //Potential energy per particle
  Epot << setw(wd) << iblk << setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;
  //Kinetic energy
  Ekin << setw(wd) << iblk << setw(wd) << stima_kin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_kin << endl;
  //Total energy
  Etot << setw(wd) << iblk << setw(wd) << stima_etot << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_etot << endl;
  //Temperature
  Temp << setw(wd) << iblk << setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;
  //Pressure
  Pres << setw(wd) << iblk << setw(wd) << stima_pres << setw(wd) << glob_av[iw]/(double)iblk << setw(wd) << err_pres << endl;
  //Radial distribution function
  if(iblk == nblk) {
    Histogram.open("./output/" + simulation + "/" + phase + "/output_histogram.dat");
    for(int ibin = 0; ibin < nbins; ibin++) {
      double r = 0.5 * double(ibin) * box / nbins;
      Histogram << setw(wd) << r << setw(wd) << hist_gl[ibin] / double(iblk) << setw(wd) << Error(hist_gl[ibin], hist_gl2[ibin], iblk) << endl;
    }
  }
     
  cout << "----------------------------" << endl << endl;

  Epot.close();
  Ekin.close();
  Etot.close();
  Temp.close();
  Pres.close();
}

// Print final configuration
void ConfFinal(void) {
  ofstream WriteConf, WriteVelocity, WriteSeed;
  cout << "Print final configuration to file config.out" << endl << endl;
  WriteConf.open("./configs/config.out");
  WriteVelocity.open("./configs/velocity.out");
  for (int i = 0; i < npart; ++i) {
    WriteConf << x[i] / box << "   " <<  y[i] / box << "   " << z[i] / box << endl;
    WriteVelocity << vx[i] << "   " <<  vy[i] << "   " << vz[i] << endl;
  }
  WriteConf.close();
  WriteVelocity.close();
  rnd.SaveSeed();
}

// Print configuration on file
void ConfXYZ(int nconf) {
  ofstream WriteXYZ;
  WriteXYZ.open("./frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i = 0; i < npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

// Periodic boundary conditions
double Pbc(double r) {
  return r - box * rint(r / box);
}

// Statistical error
double Error(double sum, double sum2, int iblk) {
  return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
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
