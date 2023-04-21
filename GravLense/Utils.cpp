#include "Utils.h"

double gauss(double x, double mean, double stdev)
{
	return (1 / (stdev * sqrt(2 * 3.1415926))) * exp(-0.5 * (x - mean) * (x - mean) / (stdev * stdev));
}

double gauss2d(double x, double y, double meanx, double meany, double stdevx, double stdevy, double theta)
{
	double a = pow(cos(theta), 2) / (2 * stdevx * stdevx) + pow(sin(theta), 2) / (2 * stdevy * stdevy);
	double b = sin(2 * theta) / 4 * (-1 / (stdevx * stdevx) + 1 / (stdevy * stdevy));
	double c = pow(sin(theta), 2) / (2 * stdevx * stdevx) + pow(cos(theta), 2) / (2 * stdevy * stdevy);
	return (1/(2 * 3.1415926 * stdevx * stdevy))*exp(-(a * (x - meanx) * (x - meanx) + 2 * b * (x - meanx) * (y - meany) + c * (y - meany) * (y - meany)));
}

double sersik(double x, double y, double x0, double y0, double phi, double I, double R, double q, double n)
{
	double k = 1.9992 * n - 0.3271;

	x = x - x0;
	y = y - y0;

	double x1 = x;

	x = x * cos(phi) - y * sin(phi);
	y = y * cos(phi) + x1 * sin(phi);

	return I * exp(-k * (pow((sqrt(q * x * x + y * y / q)) / R, 1 / n) - 1));
}

double ssersik(double x, double y, double x0, double y0, double phi, double I, double R, double q, double n, double r)
{
	vector<double> a(4, 0);
	a[0] = sersik(x + r, y, x0, y0, phi, I, R, q, n);
	a[1] = sersik(x - r, y, x0, y0, phi, I, R, q, n);
	a[2] = sersik(x, y + r, x0, y0, phi, I, R, q, n);
	a[3] = sersik(x, y - r, x0, y0, phi, I, R, q, n);
	int max = 0;
	int min = 1;
	for (int i = 0; i < 4; i++)
	{
		if (a[i] > a[max])
			max = i;
		if (a[i] < a[min])
			min = i;
	}
	double res = 0;
	for (int i = 0; i < 4; i++)
	{
		if (i != max && i != min)
		{
			res += a[i];
		}
	}
	res = res / 2;
	return res;
}