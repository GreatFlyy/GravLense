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



GeneralModel::GeneralModel()
{
	this->Dl = 1;
	this->Ds = 2;
	this->Dls = 1;
	this->scale = 1;
	this->csi0 = 1;

	this->N = 0;

	this->MassDistr = vector<vector<Comp>>(N, vector<Comp>(N, 0)); //v pixelah na kartinke

	this->ConvolveMaskX = vector<vector<Comp>>(N, vector<Comp>(N, 0));
	this->ConvolveMaskY = vector<vector<Comp>>(N, vector<Comp>(N, 0));
}

GeneralModel::GeneralModel(vector<vector<double>>& MassDistr, double Dl, double Ds, double Dls, double scale)
{
	this->Dl = Dl;
	this->Ds = Ds;
	this->Dls = Dls;
	this->scale = scale;
	this->csi0 = 1;
	
	this->N = MassDistr.size();

	this->MassDistr = vector<vector<Comp>>(N, vector<Comp>(N, 0)); //v pixelah na kartinke

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

GeneralModel::GeneralModel(function<double(double, double)>& massfunc, int N, double Dl, double Ds, double Dls, double scale)
{
	this->Dl = Dl;
	this->Ds = Ds;
	this->Dls = Dls;
	this->scale = scale;
	this->csi0 = 1;

	this->N = N;

	this->MassDistr = vector<vector<Comp>>(N, vector<Comp>(N, 0)); //v pixelah na kartinke

	this->ConvolveMaskX = vector<vector<Comp>>(N, vector<Comp>(N, 0));
	this->ConvolveMaskY = vector<vector<Comp>>(N, vector<Comp>(N, 0));

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			this->MassDistr[i][j] = massfunc(i, j);
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

void GeneralModel::setMassDistr(function<double(double, double)>& massfunc)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			this->MassDistr[i][j] = massfunc(i, j);
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






















PModel::PModel(vector<double>& par, vector<vector<double>>& data)
{
	this->n = 15;
	this->nm = 3;

	if (par.size() != n)
	{
		cout << "PMODEL: vrong parameter size";
		throw;
	}

	
	this->source = GalaxyB(par[0], par[1], par[2], par[3], par[4], par[5]);
	this->lense = GalaxyL(par[6], par[7], par[8], par[9], par[10], par[11], par[12], par[13], par[14]);
	this->MassModel = GeneralModel(this->lense.GetK, data.size(), 1., 2., 1., 1.);
	this->data = data;
	this->N = data.size();
	this->sb = vector<vector<double>>(N, vector<double>(N, 0));
	this->rb = vector<vector<double>>(N, vector<double>(N, 0));
}

void PModel::SetP(vector<double>& par, bool needtochangemass)
{
	if (par.size() != n)
	{
		cout << "PMODEL: vrong parameter size";
		throw;
	}

	this->source.x0 = par[0];
	this->source.y0 = par[1];
	this->source.theta = par[2];
	this->source.R = par[3];
	this->source.I0 = par[4];
	this->source.q = par[5];

	this->lense.x0 = par[6];
	this->lense.y0 = par[7];
	this->lense.theta = par[8];
	this->lense.R = par[9];
	this->lense.I0 = par[10];
	this->lense.qI = par[11];
	this->lense.s = par[12];
	this->lense.M = par[13];
	this->lense.qM = par[14];

	if (needtochangemass)
	{
		this->MassModel.setMassDistr(this->lense.GetK);
	}
}

double PModel::xisq()
{
	double res = 0;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			sb[i][j] = source.GetI(i, j);
		}
	}
	this->MassModel.apply(sb, rb);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			rb[i][j] += lense.GetI(i, j);
			res += pow(data[i][j] - rb[i][j], 2.);
		}
	}
	return res;
}

void PModel::grad(vector<double>& x0, vector<double>& delta, vector<double>& res)
{
	if (x0.size() != n)
	{
		cout << "PModel grad wrong x0 size(458 row)" << endl;
		throw;
	}
	if (delta.size() != n)
	{
		cout << "PModel grad wrong delta size(463 row)" << endl;
		throw;
	}
	if (res.size() != n)
	{
		cout << "PModel grad wrong res size(468 row)" << endl;
		throw;
	}

	double f0 = this->xisq();
	double f1 = 0;
	bool ntcm = false;

	for (int i = 0; i < n; i++)
	{
		if (i != 0)
		{
			x0[i - 1] -= delta[i - 1];
		}
		x0[i] += delta[i];

		ntcm = (i >= (n - nm));

		this->SetP(x0, ntcm);
		f1 = this->xisq();
		res[i] = (f1 - f0) / delta[i];
	}
	x0[n - 1] -= delta[n - 1];
}

void PModel::step(vector<double>& x0, vector<double>& delta, vector<double>& _grad, double r)
{
	this->grad(x0, delta, _grad);
	for (int i = 0; i < n; i++)
	{
		x0[i] -= r * _grad[i];
	}
	this->SetP(x0, true);
}

void PModel::Mstep(vector<double>& x0, vector<double>& x0prev, vector<double>& _grad, vector<double>& delta, double r, double v)
{
	this->grad(x0, delta, _grad);
	for (int i = 0; i < n; i++)
	{
		x0[i] = x0[i] - r * _grad[i] + v * (x0[i] - x0prev[i]);
	}
	this->SetP(x0, true);
}