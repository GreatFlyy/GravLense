#pragma once
#include<vector>
#include<functional>
#include<math.h>
#include"Reader.h"
#include"FFT.h"

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

	vector<vector<double>> apply(const vector<vector<double>>& source);

private:
	vector<vector<Comp>> MassDistr;

	int N;

	vector<vector<Comp>> ConvolveMaskX;
	vector<vector<Comp>> ConvolveMaskY;

	//vector<vector<double>> 
};




class SersicModel
{
public:
	
};

