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
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   const int ncities = 34;
   int appo;
   vector<int> f(ncities, 0);
   ifstream Primes, Seed;
   ofstream Output;
   // Primes
   Primes.open("primes32001.dat");
   for(int i = 0; i < 99; i++)
      Primes >> p1 >> p2 ;
   Primes.close();
   // Seed
   Seed.open("seed.in");
   Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   Seed.close();
   // Set random
   rnd.SetRandom(seed,p1,p2);

   Output.open("output_selection.dat");

   for(int i = 0; i < 100; i++) {
      appo = int(ncities * pow(rnd.Rannyu(), 2));
      f[appo] += 1;
      Output << setw(6) << i+1 << setw(6) << appo + 1 << endl;
   }

   for(int i = 0; i < ncities; i++)
      cout << setw(6) << i+1 << setw(6) << f[i] << endl;

   Output.close();

   return 0;
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
