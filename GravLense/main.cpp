#include<iostream>
#include<math.h>
#include"Reader.h" //Самопальная запись в различные форматы(в том числе .bmp для изображений)
#include"Matrix.h" //Самопальные матрицы, вектора, комплексные числа и кватернионы, авось понадобятся(писал давно, много можно переписать покороче, но зачем))
#include<functional>
#include"Model.h"

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

int main()
{
	string pathsource = "source.bmp";
	string pathimage = "image.bmp"; //пути к файлам

	read_bmp rbi(pathsource);
	bmp source(rbi.read());
	rbi.close(); //открытие и чтение картинки источника

	const double c = 299792458;
	const double G = 6.67430 * pow(10., -11);
	const double PI = 3.141592654;
	
	AxleMassModel model1([](double x) {return 100000 * 0.8 * (1 - exp(-x * x / 0.8)); }, 64, 64, 1, 2, 1, 20);

	vector<vector<BYTE>> r(source.widthpx, vector<BYTE>(source.heightpx, 0));
	vector<vector<BYTE>> g(source.widthpx, vector<BYTE>(source.heightpx, 0));
	vector<vector<BYTE>> b(source.widthpx, vector<BYTE>(source.heightpx, 0));

	for (int i = 0; i < source.widthpx; i++)
	{
		for (int j = 0; j < source.heightpx; j++)
		{
			r[i][j] = (int)source.pixarr[i][j][0];
			g[i][j] = (int)source.pixarr[i][j][1];
			b[i][j] = (int)source.pixarr[i][j][2];
		}
	}

	r = model1.apply(r);
	g = model1.apply(g);
	b = model1.apply(b);

	

	bmp image(source.widthpx, source.heightpx, source.bpp);

	for (int i = 0; i < source.widthpx; i++)
	{
		for (int j = 0; j < source.heightpx; j++)
		{
			image.pixarr[i][j][0] = r[i][j];
			image.pixarr[i][j][1] = g[i][j];
			image.pixarr[i][j][2] = b[i][j];
		}
	}
	
	read_bmp rbo(pathimage);
	rbo.print(image);
	rbo.close();

	bmp MassDistr(source.widthpx, source.heightpx, source.bpp);

	



	return 0;
}