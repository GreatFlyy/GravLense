#include "Model.h"

PointMassModel::PointMassModel(double M, int x0,int y0, double Dl, double Ds, double Dls, double scale)
{
	this->M = M;
	this->x0 = x0;
	this->y0 = y0;
	this->Dl = Dl;
	this->Ds = Ds;
	this->Dls = Dls;
	this->scale = scale;
	this->csi0 = sqrt(Dl/(Ds*Dls)/K1);
}

void PointMassModel::update()
{
	this->csi0 = sqrt(Dl / (Ds * Dls) / K1);
}

vector<vector<BYTE>> PointMassModel::apply(vector<vector<BYTE>>& source)
{
	vector<vector<BYTE>> image(source.size(), vector<BYTE>(source[0].size(), 0));
	for (int i = 0; i < image.size(); i++)
	{
		for (int j = 0; j < image[0].size(); j++)
		{
			double x = i - x0;
			double y = j - y0;
			double rx = sqrt(x * x + y * y) / (scale * csi0); //находим икс из формулы
			double ry = rx - 1 / rx; //сама формула

			x = x * ry / rx; //нахождние координат на сурсе (не домножаю на коеффы, т.к. отношение)
			y = y * ry / rx;

			x += x0; //привод к кордам картинки
			y += y0;

			vector<BYTE>vals(4,0);

			if (x >= image.size() - 1 || y >= image[0].size() - 1 || x <= 0 || y <= 0 || isnan((double)x) || isnan((double)y)) //проверка попадает ли наши корды в сурс
			{
				vals = vector<BYTE>(4, 0); //если в сурс не попал то очевидно 0
				
				image[i][j] = 0;					
			}
			else
			{
				vals[0] = source[(int)x][(int)y]; //сохранение значений пикселей для интерполяции
				vals[1] = source[(int)x][(int)y + 1];
				vals[2] = source[(int)x + 1][(int)y];
				vals[3] = source[(int)x + 1][(int)y + 1];


				x = x - (int)x;
				y = y - (int)y;

				//собственно интерполяция
				
				image[i][j] = vals[0] * (1 - x - y + x * y) + vals[1] * (y - x * y) + vals[2] * (x - x * y) + vals[3] * (x * y);
			}
		}
	}

	return image;
}




AxleMassModel::AxleMassModel(function<double(double)> M, double x0, double y0, double Dl, double Ds, double Dls, double scale)
{
	this->M = M;
	this->x0 = x0;
	this->y0 = y0;
	this->Dl = Dl;
	this->Ds = Ds;
	this->Dls = Dls;
	this->scale = scale;
	this->csi0 = 1;

	this->m = [&](double x)
	{
		return this->M(x * csi0) / (PI * csi0 * csi0 * K1 * Ds / (Dl * Dls));
	};
}

void AxleMassModel::update()
{
	this->m = [&](double x)
	{
		return M(x * csi0) / (PI * csi0 * csi0 * K1 * Ds / (Dl * Dls));
	};
}

vector<vector<BYTE>> AxleMassModel::apply(vector<vector<BYTE>>& source)
{
	vector<vector<BYTE>> image(source.size(), vector<BYTE>(source[0].size(), 0));
	for (int i = 0; i < image.size(); i++)
	{
		for (int j = 0; j < image[0].size(); j++)
		{
			double x = i - x0;
			double y = j - y0;
			double rx = sqrt(x * x + y * y) / (scale * csi0); //находим икс из формулы
			double ry = rx - m(rx) / rx; //сама формула

			x = x * ry / rx; //нахождние координат на сурсе (не домножаю на коеффы, т.к. отношение)
			y = y * ry / rx;

			x += x0; //привод к кордам картинки
			y += y0;

			vector<BYTE>vals(4, 0);

			if (x >= image.size() - 1 || y >= image[0].size() - 1 || x <= 0 || y <= 0 || isnan((double)x) || isnan((double)y)) //проверка попадает ли наши корды в сурс
			{
				vals = vector<BYTE>(4, 0); //если в сурс не попал то очевидно 0

				image[i][j] = 0;
			}
			else
			{
				vals[0] = source[(int)x][(int)y]; //сохранение значений пикселей для интерполяции
				vals[1] = source[(int)x][(int)y + 1];
				vals[2] = source[(int)x + 1][(int)y];
				vals[3] = source[(int)x + 1][(int)y + 1];


				x = x - (int)x;
				y = y - (int)y;

				//собственно интерполяция

				image[i][j] = vals[0] * (1 - x - y + x * y) + vals[1] * (y - x * y) + vals[2] * (x - x * y) + vals[3] * (x * y);
			}
		}
	}

	return image;
}




GeneralModel::GeneralModel(vector<vector<double>>& MassDistr, double Dl, double Ds, double Dls, double scale)
{
	this->Dl = Dl;
	this->Ds = Ds;
	this->Dls = Dls;
	this->scale = scale;
	this->csi0 = 1;
	
	this->N = MassDistr.size();

	this->MassDistr = vector<vector<Comp>>(N, vector<Comp>(N, 0));

	this->ConvolveMaskX = vector<vector<Comp>>(N, vector<Comp>(N, 0));
	this->ConvolveMaskY = vector<vector<Comp>>(N, vector<Comp>(N, 0));

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			this->MassDistr[i][j] = MassDistr[i][j];
			double x = i - N / 2 + 0.5;
			double y = j - N / 2 + 0.5;
			this->ConvolveMaskX[(i >= N / 2 ? i - N / 2 : i + N / 2)][(j >= N / 2 ? j - N / 2 : j + N / 2)] = x / (x * x + y * y);
			this->ConvolveMaskY[(i >= N / 2 ? i - N / 2 : i + N / 2)][(j >= N / 2 ? j - N / 2 : j + N / 2)] = y / (x * x + y * y);
		}
	}

	cout << "first fft" << endl;
	fft2(this->MassDistr, 1);
	cout << "second fft" << endl;
	fft2(this->ConvolveMaskX, 1);
	cout << "third fft" << endl;
	fft2(this->ConvolveMaskY, 1);

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			this->ConvolveMaskX[i][j] = this->ConvolveMaskX[i][j] * this->MassDistr[i][j];
			this->ConvolveMaskY[i][j] = this->ConvolveMaskY[i][j] * this->MassDistr[i][j];
		}
	}

	cout << "4 fft" << endl;
	fft2(this->ConvolveMaskX, -1);
	cout << "5 fft" << endl;
	fft2(this->ConvolveMaskY, -1);
}

void GeneralModel::setMassDistr(vector<vector<double>>& MassDistr)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			this->MassDistr[i][j] = MassDistr[i][j];
			double x = i - N / 2 + 0.5;
			double y = j - N / 2 + 0.5;
			this->ConvolveMaskX[(i >= N / 2 ? i - N / 2 : i + N / 2)][(j >= N / 2 ? j - N / 2 : j + N / 2)] = x / (x * x + y * y);
			this->ConvolveMaskY[(i >= N / 2 ? i - N / 2 : i + N / 2)][(j >= N / 2 ? j - N / 2 : j + N / 2)] = y / (x * x + y * y);
		}
	}

	fft2(this->MassDistr, 1);
	fft2(this->ConvolveMaskX, 1);
	fft2(this->ConvolveMaskY, 1);

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			this->ConvolveMaskX[i][j] = this->ConvolveMaskX[i][j] * this->MassDistr[i][j];
			this->ConvolveMaskY[i][j] = this->ConvolveMaskY[i][j] * this->MassDistr[i][j];
		}
	}

	fft2(this->ConvolveMaskX, -1);
	fft2(this->ConvolveMaskY, -1);
}

void GeneralModel::apply(const vector<vector<double>>& source, vector<vector<double>>& output)
{
	//vector<vector<double>> output = vector<vector<double>>(N, vector<double>(N, 0));
	if (output.size() < N || output[0].size() < N)
	{
		cout << "Model.cpp/230row/wrongsize" << endl;
		return;
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			double x = i;
			double y = j;
			x = x - this->ConvolveMaskX[i][j].r;
			y = y - this->ConvolveMaskY[i][j].r;

			vector<double> vals(4, 0);

			if (!(x > N - 2 || y > N - 2 || x < 0 || y < 0 || isnan((double)x) || isnan((double)y)))
			{
				vals[0] = source[(int)x][(int)y]; //сохранение значений пикселей для интерполяции
				vals[1] = source[(int)x][(int)y + 1];
				vals[2] = source[(int)x + 1][(int)y];
				vals[3] = source[(int)x + 1][(int)y + 1];

				x = x - (int)x;
				y = y - (int)y;

				//собственно интерполяция

				output[i][j] = vals[0] * (1 - x - y + x * y) + vals[1] * (y - x * y) + vals[2] * (x - x * y) + vals[3] * (x * y);
			}
		}
	}
}