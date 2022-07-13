#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

int main(){

	ifstream input("output_epot.dat");
	ofstream output("autocorrelation.dat");
	vector<double> data, autocorrelation;
	double appo, app1 = 0., app2 = 0., app3 = 0., app4 = 0., app5 = 0., dt, den;
	const double tmax = 500000;
	const int width = 10;

	for(int i = 0; i < tmax; i++) {
		input >> appo >> appo;
		data.push_back(appo);
	}

	input.close();

	for(int i = 0; i < tmax; i++) {
		app4 += data[i] * data[i];
		app5 += data[i];
	}

	den = app4 / tmax - pow((app5 / tmax), 2);

	for(int i = 0; i < tmax; i++) {
		app1 = app2 = app3 = 0.;
		dt = tmax - i;
		for(int j = 0; j < dt; j++){
			app1 += data[j] * data[j + i];
			app2 += data[j + i];
			app3 += data[j];
		}
		
		autocorrelation.push_back((app1 / dt - app2 * app3 / (dt * dt)) / den);
		output << setw(width) << i << " " << setw(width) << autocorrelation[i] << endl;
	}

	return 0;
}