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

double Simpson(const func& f, double a, double b, long double h)
{
	double sum1 = 0.0;
	double sum2 = 0.0;
	for (auto x = a + h; x < b; x += 2 * h)
		sum1 += f(x);
	for (auto x = a + 2 * h; x < b; x += 2 * h)
		sum2 += f(x);

	return h / 3.0 * (f(a) + f(b) + 4 * sum1 + 2 * sum2);
}

pair<matr, vect> createSLE(const func& p, const func& q, const func& f, double h, double a, double b, long double eps, IntegralFunc ifunc)
{
	size_t n = (b - a) / h;
	matr C(n, vect(n, 0));
	vect d(n, 0);

	/*for (auto k = 0; k < n; k++)
	{
		double xi = a + k * h;
		double xp = xi - h;
		double xn = xi + h;

		if (k > 0)
			C[k][k - 1] = 1.0 / h - 1.0 / pow(h, 2) * (
				+ifunc([xp, xi, &q](double x) { return q(x) * (x - xi) * (x - xp); }, xp, xi, eps)
				+ifunc([xp, xi, &p](double x) { return p(x) * (x - xp); }, xi, xn, eps)
				);

		if (k < n - 1)
			C[k][k + 1] = 1.0 / h - 1.0 / pow(h, 2) * (
				+ifunc([xi, xn, &q](double x) { return q(x) * (x - xi) * (x - xn); }, xi, xn, eps)
				+ifunc([xi, xn, &p](double x) { return p(x) * (x - xn); }, xi, xn, eps)
				);

		C[k][k] = -2.0 / h + 1.0 / pow(h, 2) * (
			+ifunc([xp, xi, &p](double x) { return p(x) * (x - xp); }, xp, xi, eps)
			+ifunc([xp, xi, &q](double x) { return q(x) * pow(x - xp, 2); }, xp, xi, eps)
			+ifunc([xi, xn, &p](double x) { return p(x) * (x - xn); }, xp, xi, eps)
			+ifunc([xi, xn, &q](double x) { return q(x) * pow(x - xn, 2); }, xp, xi, eps)
		);

		d[k] = 1.0 / h * (
			+ifunc([xp, xi, &f](double x) { return f(x) * (x - xp); }, xp, xi, eps)
			-ifunc([xi, xn, &f](double x) { return f(x) * (x - xn); }, xi, xn, eps)
			);
	}*/

	for (auto k = 0; k < n; k++)
	{
		double xi = a + k * h;
		double xp = xi - h;
		double xn = xi + h;

		if (k > 0)
		C[k][k - 1] = 1.0 / pow(h, 2) * (
			+ifunc([xp, xi, &q](double x) { return q(x) * (xi - x) * (x - xp); }, xp, xi, eps)
			-ifunc(p, xp, xi, eps)
		);

		if (k < n - 1)
		C[k][k + 1] = 1.0 / pow(h, 2) * (
			+ifunc([xi, xn, &q](double x) { return q(x) * (x - xi) * (xn - x); }, xi, xn, eps)
			-ifunc(p, xi, xn, eps)
		);

		C[k][k] = (
			2.0 / pow(h, 2) * ifunc(p, xp, xn, eps) + 1.0 / pow(h, 2) * (
				+ifunc([xp, xi, &q](double x) { return q(x) * pow(x - xp, 2); }, xp, xi, eps)
				+ifunc([xi, xn, &q](double x) { return q(x) * pow(xn - x, 2); }, xi, xn, eps)
			));

		d[k] = 1.0 / h * (
			+ifunc([xp, xi, &f](double x) { return f(x) * (x - xp); }, xp, xi, eps)
			+ifunc([xi, xn, &f](double x) { return f(x) * (xn - x); }, xi, xn, eps)
		);
	}

	return { C, d };
}

int lab2()
{
	double C = 120;
	size_t N = 9;
	size_t K = 83;
	double a = (M_PI * (N + 10));
	double b = a + K / 50.0 + 2;
	double h = 1.0 / 10;

	cout << a << endl;
	cout << b << endl;

	func p = [N, K, C](double x) { return sin(N * pow(x, K % 4 + 2)) + cos(2*K)*x + C; };
	func q = [N, K](double x) { return pow(sin(N * x), 2) * exp(cos(K * x)) + pow(cos(K * x), 2) * pow(atan(x * sin(N * x)), 2) / (K + 1); };
	func f = [a, b, &p, &q](double x) { 
		if (x == a || x == b) 
			return 0.0;

		return q(x) * pow(sin(x), 2) - p(x) * (2 * pow(cos(x), 2) - 2 * pow(sin(x), 2));
	};

	auto [A, d] = createSLE(p, q, f, h, a, b);

	auto [u, r, rn, it] = SeidelMethod(A, d, d, 10e-10);

	cout << "u =\n";
	for (auto& v : u)
		cout << v << endl;
	cout << "it = " << it << endl;


	std::ofstream fout("lab2/lab2_output.csv");
	if (!fout.is_open())
	{
		cout << "Error occuried while opening file!\n";
		return 1;
	}

	for (auto i = 0; i < u.size(); i++)
	{
		fout << a + h * i << ";" << u[i] << endl;
	}

	fout.close();

	return 0;
}
