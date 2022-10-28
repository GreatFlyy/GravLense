#include<iostream>
#include<math.h>
#include"Reader.h" //����������� ������ � ��������� �������(� ��� ����� .bmp ��� �����������)
#include"Matrix.h" //����������� �������, �������, ����������� ����� � �����������, ����� �����������(����� �����, ����� ����� ���������� ��������, �� �����))



int main()
{
	string pathsource = "source.bmp"; //���� � ������

	read_bmp rb(pathsource);
	bmp source(rb.read()); //�������� � ������ �������� ���������

	const double c = 299792458;
	const double G = 6.67430 * pow(10., -11);

	/*double M = 1;
	double Ds = 2;
	double Dl = 1;
	double Dls = 1;
	
	double tethE = sqrt((4*G*M/(c*c))*(Dls/(Dl*Ds)))
	*/ //���� ��� ��� ������, ������

	double tethE = 1;

	double scale = 40; //����. �� ������

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