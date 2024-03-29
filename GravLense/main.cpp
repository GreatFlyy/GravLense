#pragma once

#include<CCfits/CCfits>
#include<iostream>
#include<math.h>
#include"Reader.h" //����������� ������ � ��������� �������(� ��� ����� .bmp ��� �����������)
#include"Matrix.h" //����������� �������, �������, ����������� ����� � �����������, ����� �����������(����� �����, ����� ����� ���������� ��������, �� �����))
#include<functional>
#include"Model.h"
#include"FFT.h" //����������� �����
#include "FitsUt.h"
#include<ctime>

using namespace std;




//int main()
//{
//	string pathsource = "source.bmp";
//	string pathimage = "image.bmp"; //���� � ������
//
//	read_bmp rbi(pathsource);
//	bmp source(rbi.read());
//	rbi.close(); //�������� � ������ �������� ���������
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


int modeling();
int gradienting();



int main()
{
	modeling();
	return 0;
}



int modeling()
{
	string pathsource = "D:/source/repos/GravLense/GravLense/source.bmp";
	string pathmass = "D:/source/repos/GravLense/GravLense/lense.bmp";
	string pathresult = "D:/source/repos/GravLense/GravLense/image.bmp";

	//double MassCoeff = (1. / 255.) * 0.5;
	double MassCoeff = 3.;
	int N = 256;

	vector<vector<double>> mass(N, vector<double>(N, 0));

	//read_bmp rbm(pathmass);
	//bmp mass_bmp = rbm.read();
	//rbm.close();


	GalaxyB SrcG
	(
		129,		//x0
		125,		//y0
		2,		//theta
		1,		//R
		10,		//I
		1		//q
	);

	GalaxyL LnsG
	(
		128,		//x0
		128,		//y0
		0.2,	//theta
		15,		//R
		20,		//I
		1,	//qI
		30,		//s
		2500,	//M
		1		//qM
	);

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

	vector<vector<double>> massar(N, vector<double>(N, 0));

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			massar[i][j] = LnsG.GetK(i, j);
		}
	}

	string pathuuu = "D:/ProjectPics/OUTPUTT";

	srand(time(NULL) * (time(NULL) + 42732834) * (time(NULL) + 2354263));

	int ggg = 0;
	ggg = rand() % 100;

	pathuuu += to_string(ggg);
	pathuuu += ".fits";

	write2Image(result, massar, pathuuu);

	return 0;
}



int gradienting()
{
	vector<vector<double>> data(128, vector<double>(128, 0));

	//string pathin = "C:/Users/GreatFly/Desktop/HST/MAST_2022-11-10T0750/HST/na1l91010/na1l91010_mos.fits";
	//readImageE(picture, pathin, 1);

	string pathin = "D:/INPUT.fits";
	readImageP(data, pathin);

	vector<double> x0(15, 1);
	x0[0] = 68;
	x0[1] = 57;
	x0[2] = 5;
	x0[3] = 32;
	x0[4] = 7;
	x0[5] = 0.4;
	x0[6] = 58;
	x0[7] = 58;
	x0[8] = 0.4;
	x0[9] = 40;
	x0[10] = 11;
	x0[11] = 0.2;
	x0[12] = 17;
	x0[13] = 0.3;
	x0[14] = 1;
	PModel M(x0, data);


	vector<double> delta(15, 1e-14);
	vector<double> _grad(15, 0);

	vector<double> x0prev(15, 0);
	vector<double> _gradprev(15, 0);

	long double xisq0 = 0;
	long double xisq1 = 0;

	double rate = 1e-9;
	double v = 1e-6;

	cout << "x0:" << endl;
	for (int j = 0; j < M.n; j++)
	{
		cout << x0[j] << endl;
	}

	for (int k = 0; k < 1500; k++)
	{
		M.Mstep(x0, x0prev, _grad, delta, rate, v);

		xisq1 = M.xisq();
		cout << k << ". " << "xisq = " << xisq1 << "; deltaxisq = " << xisq1 - xisq0 << "; rate = " << rate << endl;

		cout << "x0:" << endl;
		for (int j = 0; j < M.n; j++)
		{
			cout << x0[j] << endl;
		}
		cout << "---------" << endl;
		if (xisq1 - xisq0 > 0)
		{
			//rate = rate / 1.3;
		}
		else
		{
			//rate *= 1.01;
		}

		xisq0 = xisq1;

		if (xisq1 < 10000)
		{
			break;
		}
	}

	//x0[0] = 59.9669;
	//x0[1] = 49.9575;
	//x0[2] = 2;
	//x0[3] = 30;
	//x0[4] = 5;
	//x0[5] = 0.7;
	//x0[6] = /*61.6503;*/ 62;
	//x0[7] = /*53.5208;*/ 54;
	//x0[8] = 0.2;
	//x0[9] = 50;
	//x0[10] = 20;
	//x0[11] = 0.5;
	//x0[12] = 15;
	//x0[13] = /*0.399999;*/ 0;
	//x0[14] = 0.8;

	//M.SetP(x0, true);
	//x0[14] = 0.9;
	//M.SetP(x0, true);
	cout << "itogovii xisq = " << M.xisq() << endl;



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

	srand(time(NULL) * (time(NULL) + 42732834) * (time(NULL) + 2354263));

	int ggg = 0;
	ggg = rand() % 100;

	pathuuu += to_string(ggg);
	pathuuu += ".fits";

	writeImage(picture, pathuuu);

	return 0;
}