#include "mathex.h"

vect_t operator+(const vect_t& a, const vect_t& b)
{
	auto n = a.size();
	assert(n == b.size());

	vect_t result(n);
	for (size_t i = 0; i < n; i++)
		result[i] = a[i] + b[i];
	return result;
}

vect_t& operator+=(vect_t& rop, const vect_t& op)
{
	auto n = rop.size();
	assert(n == op.size());
	for (size_t i = 0; i < n; i++)
		rop[i] += op[i];
	return rop;
}

vect_t operator-(const vect_t& a, const vect_t& b)
{
	auto n = a.size();
	assert(n == b.size());

	vect_t result(n);
	for (size_t i = 0; i < n; i++)
		result[i] = a[i] - b[i];
	return result;
}

vect_t& operator-=(vect_t& rop, const vect_t& op)
{
	auto n = rop.size();
	assert(n == op.size());
	for (size_t i = 0; i < n; i++)
		rop[i] -= op[i];
	return rop;
}

vect_t dot(const vect_t& a, const vect_t& b)
{
	auto n = a.size();
	assert(n == b.size());

	vect_t result(n);
	for (size_t i = 0; i < n; i++)
		result[i] = a[i] * b[i];
	return result;
}

vect_t dot(const matr_t& a, const vect_t& b)
{
	auto m = a.size();
	auto n = a.back().size();
	assert(n == b.size());

	vect_t result(m);
	for (size_t i = 0; i < m; i++)
	{
		result[i] = 0.0;
		for (size_t j = 0; j < n; j++)
			result[i] += a[i][j] * b[j];
	}
	return result;
}

vect_t& dotadd(vect_t& rvct, const matr_t& mtr, const vect_t& vct)
{
	auto m = mtr.size();
	auto n = mtr.back().size();
	assert(n == vct.size() && m == rvct.size());

	for (size_t i = 0; i < m; i++)
	{
		for (size_t j = 0; j < n; j++)
			rvct[i] += mtr[i][j] * vct[j];
	}
	return rvct;
}

double norm2(const vect_t& vct)
{
	double sum = 0.0;
	for (auto& elm : vct)
		sum += pow(elm, 2);
	return sqrt(sum);
}