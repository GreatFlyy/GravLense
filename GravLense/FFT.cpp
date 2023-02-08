#include "FFT.h"

double PI = 3.141592654;
Comp I(0, 1);

std::vector<Comp> dft(std::vector<Comp> &x, int invornot)
{
	Comp om = exp(I * invornot * (-2) * PI / x.size());
	std::vector<Comp> tf;
	tf.reserve((x.size() - 1) * (x.size() - 1) + 1);
	tf.push_back(1);
	for (int i = 1; i <= (x.size() - 1) * (x.size() - 1); ++i)
	{
		tf.push_back(tf[i - 1] * om);
	}
	std::vector<Comp> result;
	for (int k = 0; k < x.size(); ++k)
	{
		result.push_back(0);
		for (int j = 0; j < x.size(); ++j)
		{
			result[k] = result[k] + tf[k * j] * x[j];
		}
		if (invornot == -1)
		{
			result[k] = result[k] / (double)x.size();
		}
	}
	return result;
}

class Byte
{
public:
	char data;

	Byte(char data)
	{
		this->data = data;
	}

	bool get(short i)
	{
		return (bool)((((int)this->data)>>(i))%2);
	}

	void set(short i, bool data)
	{
		if (data && !this->get(i))
		{
			this->data += (char)(1 << i);
		}
		else if (!data && this->get(i))
		{
			this->data -= (char)(1 << i);
		}
	}
};

bool get_vec_byte(vector<Byte> &vec, int i)
{
	if (i < vec.size() * 8)
	{
		return vec[i / 8].get(i % 8);
	}
	else
	{
		return 0;
	}
}

void set_vec_byte(vector<Byte>& vec, int i, bool data)
{
	vec[i / 8].set(i % 8, data);
}

int move(int i, int N, int ssize)
{
	return (i % (N / ssize))* ssize + i / (N / ssize);
}

void sort(std::vector<Comp> &x, int ssize)
{
	/*int fnd = x.size() / ssize;

	vector<double> res(x.size(), 0);
	for (int i = 0; i < x.size(); i++)
	{
		res[(i % fnd) * ssize + i / fnd] = x[i];
	}

	x = res;*/

	vector<Byte> whbm(x.size()/8, (char)0);

	int k = 0;
	Comp box = 0;

	for (int i = 0; i < x.size(); i++)
	{
		if (!get_vec_byte(whbm, i))
		{
			k = i;
			box = x[k];
			while (!get_vec_byte(whbm, move(k, x.size(), ssize)))
			{				
				k = move(k, x.size(), ssize);
				swap(box, x[k]);
				set_vec_byte(whbm, k, true);
			}
		}
	}
}

void fft(std::vector<Comp>& x, int invornot)
{
	int ssize = 8;
	sort(x, ssize);
	int N = x.size();

	invornot = invornot <= -1 ? -1 : 1;

	vector<Comp> oms(ssize*ssize, 0);
	
	for (int k = 0; k < ssize; k++)
	{
		for (int j = k; j < ssize; j++)
		{
			oms[k*j] = exp(I * (-2) * invornot * PI * k * j / ssize);
		}
	}

	vector<Comp> box(ssize, 0);
	for (int i = 0; i < N / ssize; i++)
	{		
		for (int n = 0; n < ssize; n++)
		{
			box[n] = x[i * ssize + n];
		}

		for (int k = 0; k < ssize; k++)
		{
			x[i * ssize + k] = 0;
			for (int j = 0; j < ssize; j++)
			{
				x[i * ssize + k] = x[i * ssize + k] + box[j] * oms[k*j];
			}
			if (invornot == -1)
			{
				x[i * ssize + k] = x[i * ssize + k] / N;
			}
		}
	}

	oms.clear();

	

	for (int r = 0; r < ssize; r++)
	{
		box = vector<Comp>(N / ssize, 0);
		for (int n = 0; n < N / ssize; n++)
		{
			int k = n * ssize + r;
			for (int m = 0; m < N / ssize; m++)
			{
				box[n] = box[n] + exp(I * (-2) * invornot * PI * k * m / N) * x[m * ssize + r];
			}
		}

		for (int n = 0; n < N / ssize; n++)
		{
			x[n * ssize + r] = box[n];
		}
	}
}

void fft2(vector<vector<Comp>> & x, int invornot)
{
	for (int i = 0; i < x.size(); i++)
	{
		fft(x[i], invornot);
	}

	vector<Comp> box(x[0].size(), 0);

	for (int j = 0; j < x[0].size(); j++)
	{
		for (int i = 0; i < x.size(); i++)
		{
			box[i] = x[i][j];
		}
		fft(box, invornot);
		for (int i = 0; i < x.size(); i++)
		{
			x[i][j] = box[i];
		}
	}
}



