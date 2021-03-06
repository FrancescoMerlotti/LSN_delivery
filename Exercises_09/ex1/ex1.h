#ifndef __ex1_h__
#define __ex1_h__
#include <vector>
#include "../Random/random.h"

using namespace std;

// Random
int seed[4];
Random rnd;

// Classes
class City {
	public:
		City();
		City(double, double, int);
		City(City &);

		double GetX() const;
		double GetY() const;
		int GetLabel() const;
		void SetX(double);
		void SetY(double);
		void SetLabel(int);

	private:
		double m_x, m_y;
		int m_label;
};

// Simulation setup
int ncities, nconf, nsteps, gen, shape;
vector<vector<int>> population;
vector<int> best;
City * cities;

// Selection indexes
int in1, in2;

// Functions
void Input(int);
void Generation(int);
bool Check(vector<int>);
void Loading();
void Configurate();
void Display();
bool SwapRule(vector<int>, vector<int>);
void Order();
void Permutation(int);
void PermutationCluster(int);
void Reverse(int);
void Queue(int);
vector<int> CrossingOver(int, int);
void CrossingOver2(int, int);
double Distance(City &, City &);
double PathLength(vector<int>);
void Selection();
void PrintFile(int);
void PrintBest();
bool IsBest();

#endif