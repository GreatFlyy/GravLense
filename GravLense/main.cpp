#include<iostream>
#include<math.h>
#include"Reader.h" //Самопальная запись в различные форматы(в том числе .bmp для изображений)
#include"Matrix.h" //Самопальные матрицы, вектора, комплексные числа и кватернионы, авось понадобятся(писал давно, много можно переписать покороче, но зачем))



int main()
{
	string pathsource = "source.bmp"; //пути к файлам

	read_bmp rb(pathsource);
	bmp source(rb.read()); //открытие и чтение картинки источника

	const double c = 299792458;
	const double G = 6.67430 * pow(10., -11);

	/*double M = 1;
	double Ds = 2;
	double Dl = 1;
	double Dls = 1;
	
	double tethE = sqrt((4*G*M/(c*c))*(Dls/(Dl*Ds)))
	*/ //пока что без физики, тестим

	double tethE = 1;

	double scale = 40; //пикс. на радиан

	double rx = 0;
	double ry = 0;

	bmp image(source.widthpx, source.heightpx, source.bpp);

	double x = 0, y = 0;
	int x0 = 64, y0 = 64;

	for (int i = 0; i < image.widthpx; i++)
	{
		for (int j = 0; j < image.heightpx; j++)
		{
			x = i - x0;
			y = j - y0;
			rx = sqrt(x * x + y * y) / (scale * tethE);
			ry = rx + 1 / rx;

			x = -x * ry / rx;
			y = -y * ry / rx;


		}
	}





	return 0;
}