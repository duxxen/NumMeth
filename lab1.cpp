#include "lab1.h"

vect_t SeidelSolve(const matr_t& A, const vect_t& b, double eps)
{
	return SeidelSolve(A, b, vect_t(b.size(), 0), eps);
}

vect_t SeidelSolve(const matr_t& A, const vect_t& b, size_t maxit, double eps)
{
	return SeidelSolve(A, b, vect_t(b.size(), 0), maxit, eps);
}

vect_t SeidelSolve(const matr_t& A, const vect_t& b, const vect_t& x0, size_t maxit, double eps)
{
	size_t n = b.size();
	vect_t x(n);
	vect_t r(n);
	vect_t xp(n);
	bool converge = false;

	for (size_t i = 0; i < n; i++)
		xp[i] = x0[i];

	for (size_t it = 0; it < maxit; it++)
	{
		for (size_t i = 0; i < n; i++)
		{
			double sum = 0.0;
			for (size_t j = 0; j < i; j++)
				sum += A[i][j] * x[j];
			for (size_t j = i + 1; j < n; j++)
				sum += A[i][j] * xp[j];

			x[i] = (b[i] - sum) / A[i][i];
			r[i] = -b[i];
		}

		dotadd(r, A, x);
		if (norm2(r) <= eps)
			return x;

		for (size_t i = 0; i < n; i++)
			xp[i] = x[i];
	}
}

vect_t SeidelSolve(const matr_t& A, const vect_t& b, const vect_t& x0, double eps)
{
	size_t n = b.size();
	vect_t x(n);
	vect_t r(n);
	vect_t xp(n);
	bool converge = false;

	for (size_t i = 0; i < n; i++)
		xp[i] = x0[i];

	while (!converge)
	{
		for (size_t i = 0; i < n; i++)
		{
			double sum = 0.0;
			for (size_t j = 0; j < i; j++)
				sum += A[i][j] * x[j];
			for (size_t j = i + 1; j < n; j++)
				sum += A[i][j] * xp[j];

			x[i] = (b[i] - sum) / A[i][i];
			r[i] = -b[i];
		}
		dotadd(r, A, x);
		converge = norm2(r) <= eps;
		for (size_t i = 0; i < n; i++)
			xp[i] = x[i];
	}

	return x;
}