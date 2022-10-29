#include<iostream>
#include<math.h>
#include"Reader.h" //Самопальная запись в различные форматы(в том числе .bmp для изображений)
#include"Matrix.h" //Самопальные матрицы, вектора, комплексные числа и кватернионы, авось понадобятся(писал давно, много можно переписать покороче, но зачем))



int main()
{
	string pathsource = "source.bmp";
	string pathimage = "image.bmp"; //пути к файлам

	read_bmp rbi(pathsource);
	bmp source(rbi.read());
	rbi.close(); //открытие и чтение картинки источника

	const double c = 299792458;
	const double G = 6.67430 * pow(10., -11);

	/*double M = 1;
	double Ds = 2;
	double Dl = 1;
	double Dls = 1;
	
	double tethE = sqrt((4*G*M/(c*c))*(Dls/(Dl*Ds)))
	*/ //пока что без физики, тестим

	double tethE = 1;

	double scale = 20; //пикс. на радиан

	double rx = 0;
	double ry = 0;

	bmp image(source.widthpx, source.heightpx, source.bpp);

	double x = 0, y = 0;
	int x0 = 64, y0 = 64;
	vector<vector<BYTE>> vals(4, vector<BYTE>(3, 0));
	

	for (int i = 0; i < image.widthpx; i++)
	{
		for (int j = 0; j < image.heightpx; j++)
		{
			x = i - x0;
			y = j - y0;
			rx = sqrt(x * x + y * y) / (scale * tethE); //находим икс из формулы
			ry = rx - 1 / rx; //сама формула

			x = x * ry / rx; //нахождние координат на сурсе (не домножаю на коеффы, т.к. отношение)
			y = y * ry / rx;

			x += x0; //привод к кордам картинки
			y += y0;

			if (x >= image.widthpx - 1 || y >= image.heightpx - 1 || x <= 0 || y <= 0 || isnan(x) || isnan(y)) //проверка попадает ли наши корды в сурс
			{
				vals = vector<vector<BYTE>>(4, vector<BYTE>(3, 0)); //если в сурс не попал то очевидно 0
				for (int k = 0; k < 3; k++)
				{
					image.pixarr[i][j][k] = 0;
				}
			}
			else
			{
				vals[0] = source.pixarr[(int)x][(int)y]; //сохранение значений пикселей для интерполяции
				vals[1] = source.pixarr[(int)x][(int)y + 1];
				vals[2] = source.pixarr[(int)x + 1][(int)y];
				vals[3] = source.pixarr[(int)x + 1][(int)y + 1];


				x = x - (int)x;
				y = y - (int)y;

				for (int k = 0; k < 3; k++) //собственно интерполяция
				{
					image.pixarr[i][j][k] = vals[0][k] * (1 - x - y + x * y) + vals[1][k] * (y - x * y) + vals[2][k] * (x - x * y) + vals[3][k] * (x * y);
				}
			}
		}
	}

	read_bmp rbo(pathimage);
	rbo.print(image);
	rbo.close();





	return 0;
}