#include<iostream>
#include<math.h>
#include"Reader.h" //����������� ������ � ��������� �������(� ��� ����� .bmp ��� �����������)
#include"Matrix.h" //����������� �������, �������, ����������� ����� � �����������, ����� �����������(����� �����, ����� ����� ���������� ��������, �� �����))
#include<functional>
#include"Model.h"
#include"FFT.h" //����������� �����
#include"Utils.h"

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



int main()
{
	string pathsource = "D:/source/repos/GravLense/GravLense/source.bmp";
	string pathmass = "D:/source/repos/GravLense/GravLense/lense.bmp";
	string pathresult = "D:/source/repos/GravLense/GravLense/image.bmp";

	//double MassCoeff = (1. / 255.) * 0.5;
	double MassCoeff = 0.;
	int N = 128;

	vector<vector<double>> mass(N, vector<double>(N, 0));
	
	//read_bmp rbm(pathmass);
	//bmp mass_bmp = rbm.read();
	//rbm.close();


	//GalaxyB SrcG(60, 50, 2, 30, 100, 0.7);





	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			mass[i][j] = gauss2d(i, j, 60, 70, 6, 1, 0.84) * MassCoeff;
		}
	}


	vector<vector<double>> source_g(N, vector<double>(N, 0));


	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			//source_g[i][j] = SrcG.GetI(i, j);
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

	bmp output(source_bmp.widthpx, source_bmp.heightpx, source_bmp.bpp);

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			output.pixarr[i][j][1] = result[i][j];
		}
	}

	read_bmp rbo(pathresult);
	rbo.print(output);
	rbo.close();



	return 0;
}