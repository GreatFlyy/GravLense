#include<iostream>
#include<math.h>
#include"Reader.h" //����������� ������ � ��������� �������(� ��� ����� .bmp ��� �����������)
#include"Matrix.h" //����������� �������, �������, ����������� ����� � �����������, ����� �����������(����� �����, ����� ����� ���������� ��������, �� �����))



int main()
{
	string pathsource = "source.bmp";
	string pathimage = "image.bmp"; //���� � ������

	read_bmp rbi(pathsource);
	bmp source(rbi.read());
	rbi.close(); //�������� � ������ �������� ���������

	const double c = 299792458;
	const double G = 6.67430 * pow(10., -11);

	/*double M = 1;
	double Ds = 2;
	double Dl = 1;
	double Dls = 1;
	
	double tethE = sqrt((4*G*M/(c*c))*(Dls/(Dl*Ds)))
	*/ //���� ��� ��� ������, ������

	double tethE = 1;

	double scale = 20; //����. �� ������

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
			rx = sqrt(x * x + y * y) / (scale * tethE); //������� ��� �� �������
			ry = rx - 1 / rx; //���� �������

			x = x * ry / rx; //��������� ��������� �� ����� (�� �������� �� ������, �.�. ���������)
			y = y * ry / rx;

			x += x0; //������ � ������ ��������
			y += y0;

			if (x >= image.widthpx - 1 || y >= image.heightpx - 1 || x <= 0 || y <= 0 || isnan(x) || isnan(y)) //�������� �������� �� ���� ����� � ����
			{
				vals = vector<vector<BYTE>>(4, vector<BYTE>(3, 0)); //���� � ���� �� ����� �� �������� 0
				for (int k = 0; k < 3; k++)
				{
					image.pixarr[i][j][k] = 0;
				}
			}
			else
			{
				vals[0] = source.pixarr[(int)x][(int)y]; //���������� �������� �������� ��� ������������
				vals[1] = source.pixarr[(int)x][(int)y + 1];
				vals[2] = source.pixarr[(int)x + 1][(int)y];
				vals[3] = source.pixarr[(int)x + 1][(int)y + 1];


				x = x - (int)x;
				y = y - (int)y;

				for (int k = 0; k < 3; k++) //���������� ������������
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