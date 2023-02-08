#pragma once
#include<math.h>
#include "Matrix.h"

using namespace std;

double gauss(double x, double mean, double stdev)
{
	return (1 / (stdev * sqrt(2 * 3.1415926))) * exp(-0.5 * (x - mean) * (x - mean) / (stdev * stdev));
}

double gauss2d(double x, double y, double meanx, double meany, double stdevx, double stdevy, double theta)
{
	double a = pow(cos(theta), 2) / (2 * stdevx * stdevx) + pow(sin(theta), 2) / (2 * stdevy * stdevy);
	double b = sin(2 * theta) / 4 * (-1 / (stdevx * stdevx) + 1 / (stdevy * stdevy));
	double c = pow(sin(theta), 2) / (2 * stdevx * stdevx) + pow(cos(theta), 2) / (2 * stdevy * stdevy);
	return exp(-(a*(x-meanx)*(x-meanx) + 2*b*(x-meanx)*(y-meany) + c*(y-meany)*(y-meany)));
}

template<typename T>
int sign(T x)
{
	if (x > 0)
		return 1;
	else if (x == 0)
		return 0;
	else
		return -1;
}

template<typename T>
int clamp(T x, T d, T u)
{
	if (x > u)
		return u;
	else if (x < d)
		return d;
	else
		return x;
}

template<typename T>
int clamp_d(T x, T d)
{
	if (x < d)
		return d;
	else
		return x;
}

template<typename T>
int clamp_u(T x, T u)
{
	if (x > u)
		return u;
	else
		return x;
}


double sersik(double x, double y, double x0, double y0, double phi, double I, double R, double q, double n)
{
	double k = 1.9992 * n - 0.3271;

	x = x - x0;
	y = y - y0;

	int x1 = x;

	x = x * cos(phi) - y * sin(phi);
	y = y * cos(phi) + x1 * sin(phi);

	return I * exp(-k * (pow((sqrt(q * x * x + y * y / q)) / R, 1 / n) - 1));
}