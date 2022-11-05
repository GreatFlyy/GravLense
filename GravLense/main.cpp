#include<iostream>
#include<math.h>
#include"Reader.h" //—амопальна€ запись в различные форматы(в том числе .bmp дл€ изображений)
#include"Matrix.h" //—амопальные матрицы, вектора, комплексные числа и кватернионы, авось понадоб€тс€(писал давно, много можно переписать покороче, но зачем))
#include<functional>

using namespace std;

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
	

	const double K1 = 3472.825799; //c^2/4PIG приведЄнное к ћпк, а то цифры не помещаюца в комп


	
	double Ds = 2;
	double Dl = 1;
	double Dls = 1; //в ћпк
	 //пока что без физики, тестим

	double Scr = K1 * Ds / (Dl * Dls);

	double csi0 = 1;


	function<double(double)> M;
	function<double(double)> m;

	M = [](double x)
	{
		return exp(-(x-1.1)*(x-1.1)/0.3) * 10000;
	};

	m = [&M, &csi0, &PI, &Scr](double x)
	{
		return M(x * csi0) / (PI * csi0 * csi0 * Scr);
	};


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
			rx = sqrt(x * x + y * y) / (scale * csi0); //находим икс из формулы
			ry = rx - m(rx) / rx; //сама формула

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
				vals[0] = source.pixarr[(int)x][(int)y]; //сохранение значений пикселей дл€ интерпол€ции
				vals[1] = source.pixarr[(int)x][(int)y + 1];
				vals[2] = source.pixarr[(int)x + 1][(int)y];
				vals[3] = source.pixarr[(int)x + 1][(int)y + 1];


				x = x - (int)x;
				y = y - (int)y;

				for (int k = 0; k < 3; k++) //собственно интерпол€ци€
				{
					image.pixarr[i][j][k] = vals[0][k] * (1 - x - y + x * y) + vals[1][k] * (y - x * y) + vals[2][k] * (x - x * y) + vals[3][k] * (x * y);
				}
			}
		}
	}

	read_bmp rbo(pathimage);
	rbo.print(image);
	rbo.close();

	bmp MassDistr(source.widthpx, source.heightpx, source.bpp);

	



	return 0;
}