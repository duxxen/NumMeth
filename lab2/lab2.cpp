#include "lab2.h"

double quad(const func& f, double a, double b, long double eps)
{
	return 0.0;
}

double simp(const func& f, double a, double b, long double eps)
{
	return 0.0;
}

double xi(size_t i, double a, double h)
{
	return a + i * h;
}

vect fv(const func& f, const vect& xv)
{
	vect res;
	for (auto& x : xv)
		res.push_back(f(x));
	return res;
}

matr createMatr(double p, double q, double a, double b, double h, const IntegralFunc& ifunc)
{
	size_t n = (b - a) / h;
	matr A1(n, vect(n, 0));
	matr A2(n, vect(n, 0));

	fill_diagonal(A1, 2);
	fill_diagonal(A1, -1, -1);
	fill_diagonal(A1, -1, 1);
	fill_diagonal(A2, 4);
	fill_diagonal(A2, 1, -1);
	fill_diagonal(A2, 1, 1);

	return p / pow(h, 2) * A1 + q / 6 * A2;
}

vect createRightVect(const func& f, double a, double b, double h, const IntegralFunc& ifunc)
{
	size_t n = (b - a) / h;
	vect d(n, 0);
	for (int k = 1; k <= n; k++)
	{
		auto xk = xi(k, a, h);
		auto xp = xi(k - 1, a, h);
		auto xn = xi(k + 1, a, h);

		d[k - 1] = (
			+ifunc([&](double x) { return f(x) * (x - xp); }, xp, xk, h / 10)
			+ ifunc([&](double x) {return f(x) * (xn - x); }, xk, xn, h / 10)
		);
	}

	return 1 / pow(h, 2) * d;
}

int printcsv(const string& path, const matr& mdata, const string& head, const string& del)
{
	ofstream fout(path);
	if (!fout.is_open()) return -1;

	fout << head << endl;
	for (auto& str : mdata)
	{
		for (auto& val : str)
			fout << val << del;
		fout << endl;
	}

	fout.close();
	return 0;
}

int printcsv(const string& path, const vect& x, const vect& y, const string& head, const string& del)
{
	assert(x.size() == y.size());

	ofstream fout(path);
	if (!fout.is_open()) return -1;

	fout << head << endl;
	for (auto i = 0; i < x.size(); i++)
		fout << x[i] << del << y[i] << endl;

	fout.close();
	return 0;
}

int lab2()
{
	// ---------------------------------------------------------------------- CONSTANTS

	double C = -38;
	double N = 9;
	double K = 83;
	double a = M_PI * (N + 10);
	double b = a + K / 50 + 2;
	double h = 0.0001;
	
	size_t n = (b - a) / h - 1;
	h = (b - a) / n;

	string path = "lab2/data/";

	// ---------------------------------------------------------------------- FUNCTION u(x)

	func u = [&](long double x) { return pow(sin(M_PI * (x - b) / (b - a)), 2); };
	func du = [&](long double x) { return M_PI * sin(2 * M_PI * (-x + b) / (b - a)) / (b - a); };
	func ddu = [&](long double x) { return 2 * pow(M_PI, 2) * cos(2 * M_PI * (x - b) / (b - a)) / pow(b - a, 2); };

	// ---------------------------------------------------------------------- |I|

	long double pc = K * exp(10.0 * (long double)N / K);
	func p = [&](long double x) { return sin((N * x) / M_PI) + x * cosl(2 * N) + C; };
	func dp = [&](long double x) { return cos((N * x) / M_PI) * N / M_PI; };

	long double qc = N * sin(pow(K, N)) + 2 * K;
	func q = [&](long double x) { return exp(cos((N * x) / M_PI)) + pow(cos((N * x) / 10), 2) + 1; };

	// ---------------------------------------------------------------------- |II|

	func fc = [&](long double x) { return -pc * ddu(x) + qc * u(x); };
	func f = [&](long double x) { return -dp(x) * du(x) - p(x) * ddu(x) + q(x) * u(x); };

	// ---------------------------------------------------------------------- PRINT DATA

	matr mdata(n, vect(8, 0));
	for (int i = 0; i < n; i++)
	{
		auto x = xi(i, a, h);
		mdata[i][0] = x;
		mdata[i][1] = u(x);
		mdata[i][2] = du(x);
		mdata[i][3] = ddu(x);
		mdata[i][4] = p(x);
		mdata[i][5] = q(x);
		mdata[i][6] = f(x);
		mdata[i][7] = fc(x);
	}
	printcsv(path + "init.csv", mdata, "x;u(x);du(x);ddu(x);p(x);q(x);f(x);fc(x)");

	return 0;
}
