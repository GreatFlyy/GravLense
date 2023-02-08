#pragma once
#include "Matrix.h"
#include <vector>
#include<math.h>

using namespace std;

template <typename T1, typename T2>
vector<T1> vec_dot(vector<T1> a, vector<T2> b)
{
	vector<T1> res;
	for (int i = 0; i < min(a.size(), b.size()); i++)
	{
		res.push_back(a[i] * b[i]);
	}
	return res;
}

std::vector<Comp> dft(vector<Comp> &x, int invornot);

void sort(vector<Comp>& x, int ssize);

void fft(vector<Comp>& x, int invornot);

void fft2(vector<vector<Comp>>& x, int invornot);

template<typename T>
void recomb(vector<vector<T>> &x)
{
	T box;
	for (int j = 0; j < x.size(); j++)
	{
		for (int i = 0; i < x[0].size() / 2; i++)
		{
			box = x[j][i];
			x[j][i] = x[(j >= x.size() / 2 ? j - x.size() / 2 : j + x.size() / 2)][(i >= x.size() / 2 ? i - x.size() / 2 : i + x.size() / 2)];
			x[(j >= x.size() / 2 ? j - x.size() / 2 : j + x.size() / 2)][(i >= x.size() / 2 ? i - x.size() / 2 : i + x.size() / 2)] = box;
		}
	}
}