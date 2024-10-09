#include "lab2.h"
#include "../lab1/lab1.h"

double quad(const func& f, double a, double b, long double eps)
{
	double sum = 0.0;
	for (long double x = a; x < b; x += eps)
	{
		long double xn = x + eps;
		sum += f((x + xn) / 2.0) * (xn - x);
	}

	return sum;
}

double Simpson(const func& f, double a, double b, long double eps)
{
	double sum1 = 0.0;
	double sum2 = 0.0;
	for (auto x = a + eps; x < b; x += 2 * eps)
		sum1 += f(x);
	for (auto x = a + 2 * eps; x < b; x += 2 * eps)
		sum2 += f(x);

	return eps / 3.0 * (f(a) + f(b) + 4 * sum1 + 2 * sum2);
}

double a;
double b;
double h;
size_t n;

double p(double x)
{
	if (a >= x || x >= b)
		return 0.0;
	return -1.0;
}

double q(double x)
{
	if (a >= x || x >= b)
		return 0.0;
	return 1.0;
}

double f(double x)
{
	if (a >= x || x >= b)
		return 0.0;
	return pow(x - 1, 2) + 1;
}

pair<matr, vect> createSLE(bool useSimpson)
{
	auto integral = useSimpson ? Simpson : quad;
	double eps = Q_EPS;

	matr C(n, vect(n, 0));
	vect d(n, 0);

	for (int k = 0; k < n; k++)
	{
		double xi = a + k * h;
		double xp = xi - h;
		double xn = xi + h;

		d[k] = 1.0 / pow(h, 2) * (integral([=](double x) { return f(x) * (x - xp); }, xp, xi, eps) + integral([=](double x) { return f(x) * (xn - x); }, xi, xn, eps));

		if (k > 0)
			C[k][k - 1] = 1.0 / pow(h, 2) * integral([=](double x) { return q(x) * (x - xp) * (xi - x) - p(x); }, xp, xi, eps);
		if (k < n - 1)
			C[k][k + 1] = 1.0 / pow(h, 2) * integral([=](double x) { return q(x) * (x - xi) * (xn - x) - p(x); }, xi, xn, eps);
		
		C[k][k] = 1.0 / pow(h, 2) * (
			+ integral(p, xp, xn, eps) 
			+ integral([=](double x) { return q(x) * pow(x - xp, 2); }, xp, xi, eps) 
			+ integral([=](double x) { return q(x) * pow(xn - x, 2); }, xi, xn, eps));
	}

	return { C, d };
}

int lab2()
{
	a = 0;
	b = 2;
	n = 10;
	h = (b - a) / (n - 1);

	auto [C, d] = createSLE();
	auto [u, r, rn, it] = SeidelMethod(C, d, d, 1e-10);

	cout << "C = \n" << C << endl;
	cout << "d = " << d << endl;
	cout << "u = " << u << endl;
	cout << "rn = " << rn << endl;
	cout << "it = " << it << endl;

	std::ofstream fout("lab2/lab2_output.csv");
	if (!fout.is_open())
	{
		cout << "Error occuried while opening file!\n";
		return 1;
	}

	for (auto i = 0; i < u.size(); i++)
	{
		double x = a + i * h;
		fout << x << ";" << u[i] << endl;
	}
	fout.close();

	return 0;
}
