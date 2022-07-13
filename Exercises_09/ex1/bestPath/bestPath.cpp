#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

bool sorting(vector<double> a, vector<double> b) {
	return a[1] < b[1];
}

int main() {
	int wd = 10;
	double label, theta;
	ifstream Input;
	ofstream Output;
	vector<vector<double>> cities(34, {0, 0});
	// Input data
	Input.open("theta.dat");
	for(int i = 0; i < 34; i++) 
		Input >> cities[i][0] >> cities[i][1];
	Input.close();
	// Sort
	sort(cities.begin(), cities.end(), sorting);
	// Output
	Output.open("theta_sorted.dat");
	for(int i = 0; i < 34; i++)
		Output << setw(wd) << cities[i][0] << setw(wd) << cities[i][1] << endl;
	Output.close();
	return 0;
}