#pragma once
#include<vector>
#include<functional>
#include<math.h>
#include"Reader.h"
#include"FFT.h"
#include"Matrix.h"
#include"Utils.h"

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
	GeneralModel();
	GeneralModel(vector<vector<double>>& MassDistr, double Dl, double Ds, double Dls, double scale);
	GeneralModel(function<double(double, double)>& massfunc, int N, double Dl, double Ds, double Dls, double scale);

	void setMassDistr(vector<vector<double>>& MassDistr);
	void setMassDistr(function<double(double, double)>& massfunc);

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

	GalaxyB();

	GalaxyB(double _x0, double _y0, double _theta, double _R, double _I0, double _q);

	function<double(double, double)> GetIb;
};





class GalaxyL
{
public:
	double x0;
	double y0;
	double theta;

	double R;
	double I0;
	double qI;

	double s;
	double M;
	double qM;

	GalaxyL()
	{
		x0 = 0;
		y0 = 0;
		theta = 0;
		R = 1; 
		I0 = 1;
		qI = 1;
		s = 1;
		M = 1;
		qM = 1;

		this->GetI = [&](double x, double y)
		{
			return sersik(x, y, x0, y0, theta, I0, R, qI, 4);
		};

		this->GetK = [&](double x, double y)
		{
			return M * gauss2d(x, y, x0, y0, s / sqrt(qM), s * sqrt(qM), theta);
		};
	}

	GalaxyL(double x0, double y0, double theta, double R, double I0, double qI, double s, double M, double qM)
	{
		this->x0 = x0;
		this->y0 = y0;
		this->theta = theta;
		this->R = R;
		this->I0 = I0;
		this->qI = qI;
		this->s = s;
		this->M = M;
		this->qM = qM;

		this->GetI = [&](double x, double y)
		{
			return sersik(x, y, this->x0, this->y0, this->theta, this->I0, this->R, this->qI, 4);
		};

		this->GetK = [&](double x, double y)
		{
			return (this->M * gauss2d(x, y, this->x0, this->y0, this->s / sqrt(this->qM), this->s * sqrt(this->qM), this->theta));
		};
	}

	function<double(double, double)> GetI;

	function<double(double, double)> GetK;
};





class PModel
{
public:
	GalaxyB source;
	GalaxyL lense;
	GeneralModel MassModel;
	int n;
	int nm;
	vector<vector<double>> data;
	int N;

	PModel(vector<double>& par, vector<vector<double>>& data);

	void SetP(vector<double>& par, bool needtochangemass);

private:
	vector<vector<double>> sb;
	vector<vector<double>> rb;

public:
	double xisq();

	void grad(vector<double>& x0, vector<double>& delta, vector<double>& res);

	void step(vector<double>& x0, vector<double>& delta, vector<double>& _grad, double r);

	void BBstep(vector<double>& x0, vector<double>& x0prev, vector<double>& _grad, vector<double>& _gradprev, vector<double>& delta);
	
	void Mstep(vector<double>& x0, vector<double>& x0prev, vector<double>& _grad, vector<double>& delta, double r, double v);
};