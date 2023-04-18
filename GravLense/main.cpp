#pragma once

#include<CCfits/CCfits>
#include<iostream>
#include<math.h>
#include"Reader.h" //—амопальна€ запись в различные форматы(в том числе .bmp дл€ изображений)
#include"Matrix.h" //—амопальные матрицы, вектора, комплексные числа и кватернионы, авось понадоб€тс€(писал давно, много можно переписать покороче, но зачем))
#include<functional>
#include"Model.h"
#include"FFT.h" //—амопальный ‘урье
#include "FitsUt.h"
#include<ctime>

using namespace std;




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


int notmain();



int main()
{
	//notmain();

	vector<vector<double>> data(128, vector<double>(128, 0));

	//string pathin = "C:/Users/GreatFly/Desktop/HST/MAST_2022-11-10T0750/HST/na1l91010/na1l91010_mos.fits";
	//readImageE(picture, pathin, 1);

	string pathin = "D:/INPUT.fits";
	readImageP(data, pathin);

	vector<double> x0(15, 1);
	x0[0] = 60;
	x0[1] = 50;
	x0[2] = 2;
	x0[3] = 30;
	x0[4] = 5;
	x0[5] = 0.7;
	x0[6] = 62;
	x0[7] = 54;
	x0[8] = 0.2;
	x0[9] = 50;
	x0[10] = 20;
	x0[11] = 0.5;
	x0[12] = 15;
	x0[13] = 0.4;
	x0[14] = 0.8;
	PModel M(x0, data);
	

	vector<double> delta(15, 0.00001);
	vector<double> _grad(15, 0);

	vector<double> x0prev(15, 0);
	vector<double> _gradprev(15, 0);

	long double xisq0 = 0;
	long double xisq1 = 0;

	double rate = 1e-13;
	double v = 0;

	cout << "x0:" << endl;
	for (int j = 0; j < M.n; j++)
	{
		cout << x0[j] << endl;
	}

	for (int k = 0; k < 1; k++)
	{
		M.Mstep(x0, x0prev, _grad,  delta, rate, v);

		xisq1 = M.xisq();
		cout << "xisq = " << xisq1 << "; deltaxisq = " << xisq1 - xisq0 << "; rate = " << rate << endl;

		cout << "x0:" << endl;
		for (int j = 0; j < M.n; j++)
		{
			cout << x0[j] << endl;
		}

		if (xisq1 - xisq0 > 0)
		{
			rate = rate / 1.3;
		}
		else
		{
			rate *= 1.01;
		}

		xisq0 = xisq1;
	}

	vector<vector<double>> box(M.N, vector<double>(M.N, 0));
	vector<vector<double>> picture(M.N, vector<double>(M.N, 0));

	for (int i = 0; i < M.N; i++)
	{
		for (int j = 0; j < M.N; j++)
		{
			box[i][j] = M.source.GetIb(i, j);
		}
	}
	M.MassModel.apply(box, picture);
	for (int i = 0; i < M.N; i++)
	{
		for (int j = 0; j < M.N; j++)
		{
			picture[i][j] += M.lense.GetI(i, j);
		}
	}



	string pathuuu = "D:/OUTPUT";

	srand(time(NULL)*(time(NULL) + 42732834) * (time(NULL) + 2354263));

	int ggg = 0;
	ggg = rand() % 100;

	pathuuu += to_string(ggg);
	pathuuu += ".fits";

	writeImage(picture, pathuuu);

	return 0;
}



int notmain()
{
	string pathsource = "D:/source/repos/GravLense/GravLense/source.bmp";
	string pathmass = "D:/source/repos/GravLense/GravLense/lense.bmp";
	string pathresult = "D:/source/repos/GravLense/GravLense/image.bmp";

	//double MassCoeff = (1. / 255.) * 0.5;
	double MassCoeff = 3.;
	int N = 128;

	vector<vector<double>> mass(N, vector<double>(N, 0));

	//read_bmp rbm(pathmass);
	//bmp mass_bmp = rbm.read();
	//rbm.close();


	GalaxyB SrcG(60, 50, 2, 30, 5, 0.7);

	GalaxyL LnsG(62, 54, 0.2, 50, 20, 0.5, 15, 0.4, 0.8);

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			mass[i][j] = LnsG.GetK(i, j);
		}
	}

	vector<vector<double>> source_g(N, vector<double>(N, 0));

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			source_g[i][j] = SrcG.GetIb(i, j);
		}
	}

	read_bmp rbs(pathsource);
	bmp source_bmp = rbs.read();
	rbs.close();

	/*for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			source_g[i][j] = source_bmp.pixarr[i][j][1];
		}
	}*/

	vector<vector<double>> result(N, vector<double>(N, 0));

	cout << "making model" << endl;
	GeneralModel GM(mass, 2, 1, 1, 1);

	cout << "applying model" << endl;
	GM.apply(source_g, result);

	//bmp output(source_bmp.widthpx, source_bmp.heightpx, source_bmp.bpp);

	//for (int i = 0; i < N; i++)
	//{
	//	for (int j = 0; j < N; j++)
	//	{
	//		output.pixarr[i][j][1] = result[i][j] + LnsG.GetI(i, j);
	//		//output.pixarr[i][j][2] = LnsG.GetI(i, j);
	//		//output.pixarr[i][j][0] = LnsG.GetK(i, j)*500;
	//	}
	//}

	//read_bmp rbo(pathresult);
	//rbo.print(output);
	//rbo.close();

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			result[i][j] = result[i][j] + LnsG.GetI(i, j);
			//output.pixarr[i][j][2] = LnsG.GetI(i, j);
			//output.pixarr[i][j][0] = LnsG.GetK(i, j)*500;
		}
	}

	string pathuuu = "D:/OUTPUTT";

	srand(time(NULL) * (time(NULL) + 42732834) * (time(NULL) + 2354263));

	int ggg = 0;
	ggg = rand() % 100;

	pathuuu += to_string(ggg);
	pathuuu += ".fits";

	writeImage(result, pathuuu);

	return 0;
}