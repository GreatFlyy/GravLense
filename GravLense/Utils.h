#pragma once
#include<math.h>

using namespace std;



double gauss(double x, double mean, double stdev);

double gauss2d(double x, double y, double meanx, double meany, double stdevx, double stdevy, double theta);

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


double sersik(double x, double y, double x0, double y0, double phi, double I, double R, double q, double n);