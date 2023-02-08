#include<iostream>
#include<math.h>
#include"Reader.h" //—амопальна€ запись в различные форматы(в том числе .bmp дл€ изображений)
#include"Matrix.h" //—амопальные матрицы, вектора, комплексные числа и кватернионы, авось понадоб€тс€(писал давно, много можно переписать покороче, но зачем))
#include<functional>
#include"Model.h"
#include"FFT.h" //—амопальный ‘урье

using namespace std;

double gauss(double x)
{
	return exp(-x * x);
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


//int main()
//{
//	string pathsource = "source.bmp";
//	string pathimage = "image.bmp"; //пути к файлам
//
//	read_bmp rbi(pathsource);
//	bmp source(rbi.read());
//	rbi.close(); //открытие и чтение картинки источника
//
//	const double c = 299792458;
//	const double G = 6.67430 * pow(10., -11);
//	const double PI = 3.141592654;
//
//	bmp image(256, 256, 24);
//
//	for (int j = 0; j < image.heightpx; j++)
//	{
//		for (int i = 0; i < image.widthpx; i++)
//		{
//			int x = i - 128;
//			int y = j - 128;
//
//			image.pixarr[j][i][1] = sersik(x, y, 0., 0., 0.5,  0.6, 256, 2, 4);
//		}
//	}
//
//	read_bmp rbo(pathimage);
//
//	rbo.print(image);
//
//	return 0;
//}



int main()
{
	string pathsource = "D:/source/repos/GravLense/GravLense/source.bmp";
	string pathmass = "D:/source/repos/GravLense/GravLense/lense.bmp";
	string pathresult = "D:/source/repos/GravLense/GravLense/image.bmp";

	double MassCoeff = (1. / 255.) * 0.5;
	int N = 128;

	vector<vector<double>> mass(N, vector<double>(N, 0));
	
	read_bmp rbm(pathmass);
	bmp mass_bmp = rbm.read();
	rbm.close();

	for (int i = 0; i < mass_bmp.heightpx; i++)
	{
		for (int j = 0; j < mass_bmp.widthpx; j++)
		{
			mass[i][j] = (mass_bmp.pixarr[i][j][0] + mass_bmp.pixarr[i][j][1] + mass_bmp.pixarr[i][j][2]) / 3 * MassCoeff;
		}
	}


	vector<vector<double>> source_g(N, vector<double>(N, 0));

	read_bmp rbs(pathsource);
	bmp source_bmp = rbs.read();
	rbs.close();

	for (int i = 0; i < mass_bmp.heightpx; i++)
	{
		for (int j = 0; j < mass_bmp.widthpx; j++)
		{
			source_g[i][j] = source_bmp.pixarr[i][j][1];
		}
	}

	vector<vector<double>> result(N, vector<double>(N, 0));

	cout << "making model" << endl;
	GeneralModel GM(mass, 2, 1, 1, 1);

	cout << "applying model" << endl;
	GM.apply(source_g, result);

	bmp output(source_bmp.widthpx, source_bmp.heightpx, source_bmp.bpp);

	for (int i = 0; i < mass_bmp.heightpx; i++)
	{
		for (int j = 0; j < mass_bmp.widthpx; j++)
		{
			output.pixarr[i][j][1] = result[i][j];
		}
	}

	read_bmp rbo(pathresult);
	rbo.print(output);
	rbo.close();



	return 0;
}