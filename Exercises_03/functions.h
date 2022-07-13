#ifndef __functions_h__
#define __functions_h__

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double statUncertenty(vector<double> avg_prog, vector<double> avg2_prog, int index) {
	double sigma = 0.;
	if(index == 0)
		return sigma;
	else {
		sigma = sqrt((avg2_prog[index] - pow(avg_prog[index], 2)) / index);
		return sigma;
	}
}

double statUncertenty(double avg1, double avg2, int index) {
	double sigma = 0.;
	if (index == 0)
		return sigma;
	else {
		sigma = sqrt((avg2 - pow(avg1, 2)) / index);
		return sigma;
	}
}

bool hit(double x1, double x2) {
	if (int(x1) != int(x2))
		return true;
	if (x1 == int(x1))
		return true;
	if (x2 == int(x2))
		return true;
	else
		return false;
}

double integrand(double x) {
	return (M_PI / 2.) * cos(M_PI * x / 2.);
}

double pdf(double x) {
	return 1. - M_PI * (x - 0.5) / 2.;
}

double pdf_ar(double x) {
	return 1.5 * (1. - pow(x, 2));
}

double value(double y) {
	double a = M_PI / 4.;
	double b = (1 + M_PI / 4.);
	return (b - sqrt(pow(b, 2) - 4. * a * y)) / (2. * a);
}

vector<double> move(double displacement, int direction) {
	vector<double> movement{ 0., 0., 0. };

	if (direction % 2 == 0)
		movement[direction / 2] += displacement;
	if (direction % 2 == 1)
		movement[direction / 2] -= displacement;

	return movement;
}

vector<double> move(double displacement, double theta, double phi) {
	vector<double> movement{ displacement * sin(theta) * cos(phi), displacement * sin(theta) * sin(phi), displacement * cos(theta) };
	return movement;
}

double magnitudeSquared(vector<double>& position) {
	if (position.size() != 3)
		return -2;
	return pow(position[0], 2) + pow(position[1], 2) + pow(position[2], 2);
}

#endif