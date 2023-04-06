#pragma once
#include<vector>
#include<functional>
#include<math.h>
#include"Reader.h"
#include"FFT.h"
#include"Matrix.h"
//#include"Utils.h"

using namespace std;

const double c = 299792458;
const double G = 6.67430 * pow(10., -11);
const double PI = 3.141592654;

const double K1 = 3472.825799; //c^2/4PIG приведЄнное к ћпк, а то цифры не помещаюца в комп

class Model
{
public:
	double Dl;
	double Ds;
	double Dls;

	double scale;
	double csi0;

	//virtual void update() = 0;

	//virtual vector<vector<BYTE>> apply(vector<vector<BYTE>>& source) = 0;
};

class PointMassModel : public Model
{
public:
	double M;

	double x0;
	double y0;

	PointMassModel(double M, int x0, int y0, double Dl, double Ds, double Dls, double scale);

	void update();

	vector<vector<BYTE>> apply(vector<vector<BYTE>>& source);
};  

class AxleMassModel : public Model
{
public:
	function<double(double)> M;

	function<double(double)> m;

	double x0;
	double y0;

	AxleMassModel(function<double(double)> M, double x0, double y0, double Dl, double Ds, double Dls, double scale);

	void update();

	vector<vector<BYTE>> apply(vector<vector<BYTE>>& source);
};

class GeneralModel : public Model
{
public:
	GeneralModel(vector<vector<double>>& MassDistr, double Dl, double Ds, double Dls, double scale);

	void setMassDistr(vector<vector<double>>& MassDistr);

	void apply(const vector<vector<double>>& source, vector<vector<double>>& output);

private:
	vector<vector<Comp>> MassDistr;

	int N;

	vector<vector<Comp>> ConvolveMaskX;
	vector<vector<Comp>> ConvolveMaskY;

	//vector<vector<double>> 
};



class GalaxyB
{
public:
	double x0;
	double y0;
	double theta;
	double R;
	double I0;
	double q;

	GalaxyB()
	{
		x0 = 0;
		y0 = 0;
		theta = 0;
		R = 1;
		I0 = 1;
		q = 1;
	}

	GalaxyB(double x0, double y0, double theta, double R, double I0, double q)
	{
		this->x0 = x0;
		this->y0 = y0;
		this->theta = theta;
		this->R = R;
		this->I0 = I0;
		this->q = q;
	}

	double GetI(double x, double y)
	{
		return sersik(x, y, x0, y0, theta, I0, R, q, 4);
	};
};





class GalaxyL
{
public:
	double x0;
	double y0;
	double theta;

	double R;
	double I0;
	double q;

	double s;
	double M;

	GalaxyL()
	{
		x0 = 0;
		y0 = 0;
		theta = 0;
		R = 1; 
		I0 = 1;
		q = 1;
		s = 1;
		M = 1;
	}

	GalaxyL(double x0, double y0, double theta, double R, double I0, double q, double s, double M)
	{
		this->x0 = x0;
		this->y0 = y0;
		this->theta = theta;
		this->R = R;
		this->I0 = I0;
		this->q = q;
		this->s = s;
		this->M = M;
	}

	double GetI(double x, double y)
	{
		return sersik(x, y, x0, y0, theta, I0, R, q, 4);
	}

	double GetK(double x, double y)
	{
		return M * gauss2d(x, y, x0, y0, s / sqrt(q), s * sqrt(q), theta);
	}
};


