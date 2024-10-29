#include "lab2.h"
#include "../lab1/lab1.h"

long double quad(const func& f, long double a, long double b, long double eps)
{
	size_t n = (b - a) / eps;
	long double sum = 0.0;

	auto start = system_clock::now();
	auto current = system_clock::now();
	size_t temp = 0;

	for (size_t i = 0; i < n; i++)
	{
		long double x = a + i * eps;
		long double xn = x + eps;
		sum += f((x + xn) / 2.0) * (xn - x);

		current = system_clock::now();
		if (duration_cast<seconds>(current - start).count() > 3)
		{
			start = current;
			double per = (double)i / n * 100.0;
			cout << format("\rCalculating Quad()... [{}/{}] ({:.2f}%)", i, n, per);
		}
	}
	

	return sum;
}

long double Simpson(const func& f, long double a, long double b, long double h)
{
	long double sum1 = 0.0;
	long double sum2 = 0.0;
	for (auto x = a + h; x < b; x += 2 * h)
		sum1 += f(x);
	for (auto x = a + 2 * h; x < b; x += 2 * h)
		sum2 += f(x);

	return h / 3.0 * (f(a) + f(b) + 4 * sum1 + 2 * sum2);
}

// functions p, q method
pair<matr, vect> createSLE(const func& p, const func& q, const func& f, long double h, long double a, double b, long double eps, IntegralFunc ifunc)
{
	size_t n = (b - a) / h;
	matr C(n, vect(n, 0));
	vect d(n, 0);

	for (auto k = 0; k < n; k++)
	{
		long double xi = a + (k + 1) * h;
		long double xp = xi - h;
		long double xn = xi + h;
		
		auto m = 1.0 / pow(h, 3); // коэффициент матриц 
		if (k > 0) // над главной 
		{
			auto ip = -ifunc(p, xp, xi, eps);
			auto iq = +ifunc([=](double x) { return q(x) * (x - xp) * (xi - x); }, xp, xi, eps);
			C[k][k - 1] = m * (ip + iq); 
		}

		if (k < n - 1) // под главной
		{
			auto ip = -ifunc(p, xi, xn, eps);
			auto iq = +ifunc([=](double x) { return q(x) * (x - xi) * (xn - x); }, xi, xn, eps);
			C[k][k + 1] = m * (ip + iq);
		}

		// главная диагональ
		auto ip = ifunc(p, xp, xn, eps);
		auto iq1 = ifunc([=](double x) { return q(x) * pow(x - xp, 2); }, xp, xi, eps);
		auto iq2 = ifunc([=](double x) { return q(x) * pow(xn - x, 2); }, xi, xn, eps);
		C[k][k] = m * (ip + iq1 + iq2);

		d[k] = (
			+ifunc([=](long double x) { return f(x) * (x - xp); }, xp, xi, eps)
			+ifunc([=](long double x) { return f(x) * (xn - x); }, xi, xn, eps)
		);
	}

	d = d * (1.0 / pow(h, 2));
	return { C, d };
}

// const p, q method
pair<matr, vect> createSLE(long double p, long double q, const func& f, long double h, long double a, double b, long double eps, IntegralFunc ifunc) 
{
	size_t n = (b - a) / h;
	matr C(n, vect(n, 0));
	vect d(n, 0);

	matr A1 = createThermalMatrix(n, 1);
	matr A2 = createThermalMatrix(n, 1);
	for (auto& s : A2)
		for (auto& v : s)
			v *= v;

	for (auto k = 0; k < n; k++)
	{
		long double xi = a + k * h;
		long double xp = xi - h;
		long double xn = xi + h;

		d[k] = (
			+ifunc([=](long double x) { return f(x) * (x - xp); }, xp, xi, eps)
			+ifunc([=](long double x) { return f(x) * (xn - x); }, xi, xn, eps)
		);
	}
	A1 = (p / pow(h, 2)) * A1;
	A2 = (q / 6.0) * A2;

	C = A1 + A2;
	d = d * (1.0 / pow(h, 2));
	return { C, d };
}

vect getFuncVect(const func& u, long double a, long double b, long double h)
{
	size_t n = (b - a) / h;
	vect uv(n, 0);
	for (int i = 0; i < n; i++)
	{
		auto x = a + i * h;
		uv[i] = u(x);
	}
	return uv;
}

int lab2()
{
	long double C = -38;
	size_t N = 9;
	size_t K = 83;
	long double a = (M_PI * (N + 10));
	long double b = a + K / 50.0 + 2;
	long double h = 1.0 / 10;
	size_t n = (b - a) / h;
	h = (b - a) / n;

	// ---------------------------------------------------------------------- |I|

	long double pc = K * exp(10.0 * (long double)N / K);
	func p = [=](long double x) { return sin((N * x) / M_PI) + x * cosl(2 * N) + C; };
	func dp = [=](long double x) { return cos((N * x) / M_PI) * N / M_PI; };

	long double qc = N * sin(pow(K, N)) + 2 * K;
	func q = [=](long double x) { return exp(cos((N * x) / M_PI)) + pow(cos((N * x) / 10), 2) + 1; };

	// ---------------------------------------------------------------------- |II|

	func u = [=](long double x) { return pow(sin(M_PI * (x - b) / (b - a)), 2); };
	func du = [=](long double x) { return M_PI * sin(2 * M_PI * (-x + b) / (b - a)) / (b - a); };
	func ddu = [=](long double x) { return 2 * pow(M_PI, 2) * cos(2 * M_PI * (x - b) / (b - a)) / pow(b - a, 2); };
	
	func fc = [=](long double x) { return -pc * ddu(x) + qc * u(x); };
	func f = [=](long double x) { return -dp(x) * du(x) - p(x) * ddu(x) + q(x) * u(x); };

	// ---------------------------------------------------------------------- |III|

	auto [A, d] = createSLE(p, q, f, h, a, b, h/10.0);
	cout << format("SLE created!\nsize = {}\n", A.size());
	auto [ux, r, rn, it] = SeidelMethod(A, d, vect(n, 0), 1e-10);
	cout << format("u(x) found!\nrn = {}\nit = {}\n", rn, it);
	auto uv = getFuncVect(u, a, b, h);

	// ---------------------------------------------------------------------- |IV|

	ofstream fout("lab2/lab2_s4_test.csv");
	fout << "x;u(x);ux(x)\n";
	for (int i = 0; i < n; i++)
	{
		long double x = a + i * h;
		fout << format("{};{};{}\n", x, uv[i], ux[i]);
	}

	// ---------------------------------------------------------------------- |V|

	/*size_t t = 100;
	matr um(t);
	for (int i = 0; i < t; i++)
	{
		long double h5 = 1.0 / (10 * (i + 1));
		size_t n5 = (b - a) / h5;
		h5 = (b - a) / n5;
		auto [A5, d5] = createSLE(pc, qc, fc, h5, a, b, h5 / 10.0);
		cout << format("[{}]: SLE created!\nsize = {}\n", i, A5.size());
		auto [ux5, r5, rn5, it5] = SeidelMethod(A5, d5, d5, 1e-10);
		cout << format("[{}]: u(x) found!\nrn = {}\nit = {}\n", i, rn5, it5);
		auto uv5 = getFuncVect(u, a, b, h5);
	}*/



	return 0;
}
