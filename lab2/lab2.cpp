#include "lab2.h"
#include "../lab1/lab1.h"

int K, N, C;
double a, b;
double cp, cq;

double fp(double x) 
{
	return sin(N * pow(x, (K % 4) + 2)) + cos(2 * K) * x + C;
}

double fq(double x) 
{
	return pow(sin(N * x), 2) * exp(cos(K * x)) + pow(cos(K * x), 2) * pow(atan(x * sin(N * x)), 2) / double(K + 1);
}

double dfp(double x)
{
	double s = (K % 4 + 2);
	return s * N * pow(x, s - 1) * cos(N * pow(x, s)) + cos(2 * K);
}

double u(double x)
{
	return pow(sin(M_PI * (x - b) / (b - a)), 2);
}

double du(double x)
{
	return M_PI * cos(M_PI * (x - b) / (b - a)) / (b - a);
}

double ddu(double x)
{
	return 2 * pow(M_PI, 2) * cos(2 * M_PI * (x - b) / (b - a)) / pow(b - a, 2);
}

double fc(double x)
{
	return -dfp(x) * du(x) - fp(x) * ddu(x) + fq(x) * u(x);
}

double ff(double x)
{
	return -cp * ddu(x) + cq * u(x);
}

double quad(func f, double a, double b, double q)
{
	double sum = 0.0;
	for (long double x = a; x <= b; x += q)
	{
		long double xn = x + q;
		sum += f((x + xn) / 2.0) * (xn - x);
	}
	return sum;
}

double simp(func f, double a, double b, double q)
{
	return 0;
}

pair<matr, vect> createCSLE(double h)
{
	size_t n = (b - a) / h;
	matr A1 = createThermalMatrix(n, 1);
	matr A2 = A1;
	for (auto& str : A2)
		for (auto& v : str)
			v *= v;

	cout << A1 << endl << A2 << endl;

	return std::make_pair(matr(), vect());
}

pair<matr, vect> createFSLE(double h)
{
	return std::make_pair(matr(), vect());
}

void lab2()
{
	N = 9;
	K = 83;
	a = M_PI * (N + 10);
	b = a + K / 50.0 + 2;
	double h = 0.1;
	createCSLE(h);
}
