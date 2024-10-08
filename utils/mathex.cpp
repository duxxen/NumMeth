#include "mathex.h"

int randint(int min, int max)
{
	return min + rand() % (max - min + 1);
}

float randfloat(float min, float max)
{
	return min + (float)rand() / RAND_MAX * (max - min);
}

vect randintVect(size_t n, int min, int max)
{
	vect x(n);
	for (auto i = 0; i < n; i++)
		x[i] = randint(min, max);

	return x;
}

vect randfloatVect(size_t n, float min, float max)
{
	vect x(n);
	for (auto i = 0; i < n; i++)
		x[i] = randfloat(min, max);

	return x;
}

double norm2(const vect& v)
{
	double res = 0.0;
	for (auto& val : v)
		res += val * val;

	return sqrt(res);
}

double scalar(const vect& v1, const vect& v2)
{
	assert(v1.size() == v2.size());

	double scalar = 0;
	for (auto i = 0; i < v1.size(); i++)
		scalar += v1[i] * v2[i];

	return scalar;
}

vect operator +(const vect& v1, const vect& v2)
{
	assert(v1.size() == v2.size());

	vect res(v1.size());
	for (auto i = 0; i < v1.size(); i++)
		res[i] = v1[i] + v2[i];

	return res;
}

vect operator -(const vect& v1, const vect& v2)
{
	assert(v1.size() == v2.size());

	vect res(v1.size());
	for (auto i = 0; i < v1.size(); i++)
		res[i] = v1[i] - v2[i];

	return res;
}

vect operator *(double c, const vect& v)
{
	auto res = v;
	for (auto& val : res)
		val = c * val;

	return res;
}

vect operator *(const vect& v, double c)
{
	return c * v;
}

vect operator /(const vect& v, double c)
{
	auto res = v;
	for (auto& val : res)
		val = val / c;

	return res;
}

matr transposed(const matr& m)
{
	matr res(m.back().size(), vect(m.size()));
	for (auto i = 0; i < res.size(); i++)
		for (auto j = 0; j < m.size(); j++)
			res[i][j] = m[j][i];

	return res;
}

matr operator +(const matr& v1, const matr& v2)
{
	assert(v1.size() == v2.size());

	matr res(v1.size());
	for (auto i = 0; i < v1.size(); i++)
		res[i] = v1[i] + v2[i];

	return res;
}

matr operator -(const matr& v1, const matr& v2)
{
	assert(v1.size() == v2.size());

	matr res(v1.size());
	for (auto i = 0; i < v1.size(); i++)
		res[i] = v1[i] - v2[i];

	return res;
}

matr operator *(double c, const matr& m)
{
	auto res = m;
	for (auto& val : res)
		val = c * val;

	return res;
}

matr operator *(const matr& m, double c)
{
	return c * m;
}

matr operator /(const matr& m, double c)
{
	auto res = m;
	for (auto& val : res)
		val = val / c;

	return res;
}

matr operator*(const matr& m1, const matr& m2)
{
	auto l = m1.size();
	auto m = m1.back().size();
	auto n = m2.back().size();
	assert(m == m2.size());

	matr res(l, vect(n));
	for (auto i = 0; i < l; i++)
		for (auto j = 0; j < n; j++)
		{
			double sum = 0;
			for (auto k = 0; k < m; k++)
				sum += m1[i][k] * m2[k][j];
			res[i][j] = sum;
		}

	return res;
}

vect operator*(const matr& m, const vect& v)
{
	assert(m.back().size() == v.size());

	vect res(m.size());
	for (auto i = 0; i < res.size(); i++)
	{
		double sum = 0;
		for (auto k = 0; k < v.size(); k++)
			sum += m[i][k] * v[k];
		res[i] = sum;
	}

	return res;
}

vect operator*(const vect& v, const matr& m)
{
	assert(v.size() == m.size());

	vect res(m.back().size());
	for (auto i = 0; i < res.size(); i++)
	{
		double sum = 0;
		for (auto k = 0; k < v.size(); k++)
			sum += m[k][i] * v[k];
		res[i] = sum;
	}

	return res;
}