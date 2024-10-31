#include "lab2.h"

string GLOBAL_PATH;

double quad(const func& f, double a, double b, long double h)
{
	size_t n = (b - a) / h;
	double sum = 0;
	for (int i = 0; i <= n; i++)
		sum += f(a + i * h);
	return h * sum;
}

double simp(const func& f, double a, double b, long double h)
{
	double sum1 = 0;
	double sum2 = 0;
	for (auto x = a + h; x < b; x += 2 * h)
		sum1 += f(x);
	for (auto x = a + 2 * h; x < b; x += 2 * h)
		sum2 += f(x);

	return h / 3.0 * (f(a) + f(b) + 4 * sum1 + 2 * sum2);
}

vect createXVect(double a, double b, size_t n)
{
	double h = (b - a) / n;
	vect vx(n, 0);
	for (auto i = 0; i < n; i++)
		vx[i] = a + i * h;
	return vx;
}

vect createVectFunc(const func& f, const vect& vx)
{
	vect res;
	for (auto& x : vx)
		res.push_back(f(x));
	return res;
}

matr createMatr(double p, double q, double a, double b, double h)
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

	A1 = A1 * (p / pow(h, 2));
	A2 = A2 * (q / 6);

	return A1 + A2;
}

matr createMatr(const func& p, const func& q, double a, double b, double h, const IntegralFunc& ifunc)
{
	size_t n = (b - a) / h;
	matr A(n, vect(n, 0));

	double xk = a + h;
	double xp = xk - h;
	double xn = xk + h;

	double m3 = (1 / pow(h, 3));
	double m2 = (1 / pow(h, 2));

	double under = 1 * (-m3 * ifunc(p, xp, xk, h / 10) + m2 * ifunc([&](double x) { return q(x) * (x - xp) * (xk - x); }, xk, xn, h / 10));
	double above = 1 * (-m3 * ifunc(p, xk, xn, h / 10) + m2 * ifunc([&](double x) { return q(x) * (x - xk) * (xn - x); }, xk, xn, h / 10));
	double diagm = 1 * (m3 * ifunc(p, xp, xn, h / 10) + m2 * (ifunc([&](double x) { return q(x) * pow(x - xp, 2); }, xp, xk, h / 10) + ifunc([&](double x) { return q(x) * pow(xn - x, 2); }, xk, xn, h / 10)));

	fill_diagonal(A, under, -1);
	fill_diagonal(A, diagm, 0);
	fill_diagonal(A, above, 1);

	return A;
}

vect createRightVect(const func& f, double a, double b, double h, const IntegralFunc& ifunc)
{
	size_t n = (b - a) / h;
	vect d(n, 0);
	for (int k = 1; k <= n; k++)
	{
		auto xk = a + k * h;
		auto xp = xk - h;
		auto xn = xk + h;

		d[k - 1] = (
			+ifunc([&](double x) { return f(x) * (x - xp); }, xp, xk, h / 10)
			+ifunc([&](double x) { return f(x) * (xn - x); }, xk, xn, h / 10)
		);
	}

	return 1 / pow(h, 2) * d;
}

tuple<vect, double, size_t> SLESolver(const matr& a, const vect& b, const vect& x0, long double eps)
{
	auto n = a.size();
	vect x = x0;
	vect p;
	matr c = a;
	vect d = b;
	double en = 0;
	double rn = 0;

	double max = 0;
	for (int i = 0; i < n; i++)
	{
		double sum = 0;
		for (int j = 0; j < n; j++)
			sum += abs(a[i][j]);
		max = std::max(max, sum);
	}

	if (max >= 1)
	{
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				if (i != j) c[i][j] = -a[i][j] / a[i][i];
				else		c[i][j] = 0.0;
			}
			d[i] = b[i] / a[i][i];
		}
	}

	size_t m = 0;
	do
	{
		p = x;
		for (int i = 0; i < n; i++)
		{
			double sum = 0;
			for (int j = 0; j < i - 1; j++)
				sum += c[i][j] * x[j];
			for (int j = i; j < n; j++)
				sum += c[i][j] * p[j];

			x[i] = sum + d[i];
		}

		en = norm2(x - p);

		m++;
	} while (en > eps);

	return { x, en, m };
}

int printcsv(const string& path, const matr& mdata, const string& head, const string& del)
{
	ofstream fout(GLOBAL_PATH + path);
	if (!fout.is_open()) throw -1;

	if (!head.empty())	fout << head << endl;
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

	ofstream fout(GLOBAL_PATH + path);
	if (!fout.is_open()) throw -1;

	if (!head.empty())	fout << head << endl;
	for (auto i = 0; i < x.size(); i++)
		fout << x[i] << del << y[i] << endl;

	fout.close();
	return 0;
}

int lab2(const IntegralFunc& ifunc)
{
	// ------------------------- CONSTANTS INIT -----------------------------

	double S_EPS = 1e-10;
	double C = -38;
	double N = 9;
	double K = 83;
	double a = M_PI * (N + 10);
	double b = a + K / 50 + 2;
	double h = 0.1;

	size_t n = (b - a) / h - 1;
	h = (b - a) / n;

	// -------------------------- FUNCTION u(x) -----------------------------

	func u = [&](long double x) { return pow(sin(M_PI * (x - b) / (b - a)), 2); };
	func du = [&](long double x) { return M_PI * sin(2 * M_PI * (-x + b) / (b - a)) / (b - a); };
	func ddu = [&](long double x) { return 2 * pow(M_PI, 2) * cos(2 * M_PI * (x - b) / (b - a)) / pow(b - a, 2); };

	// ------------------------------- |I| ----------------------------------

	double pc = K * exp(10.0 * N / K);
	func pfc = [&](double x) { return pc; }; // const func

	func p = [&](double x) { return sin((N * x) / M_PI) + x * cosl(2 * N) + C; };
	func dp = [&](double x) { return cos((N * x) / M_PI) * N / M_PI; };

	double qc = N * sin(pow(K, N)) + 2 * K;
	func qfc = [&](double x) { return qc; }; // const func

	func q = [&](double x) { return exp(cos((N * x) / M_PI)) + pow(cos((N * x) / 10), 2) + 1; };

	// ------------------------------- |II| ----------------------------------

	func fc = [&](double x) { return -pc * ddu(x) + qc * u(x); };
	func f = [&](double x) { return -dp(x) * du(x) - p(x) * ddu(x) + q(x) * u(x); };

	// ---------------------------- PRINT DATA ------------------------------- 

	cout << "\n[INPUT DATA]\n";
	cout << format("\tN: {}\n\tK: {}\n\tC: {}\n\tp: {:.2f}\n\tq: {:.2f}\n", N, K, C, pc, qc);

	matr mdata(n, vect(8, 0));
	for (int i = 0; i < n; i++)
	{
		auto x = a + i * h;
		mdata[i][0] = x;		//			0:		x		
		mdata[i][1] = u(x);		//			1:		u(x)	
		mdata[i][2] = du(x);	//			2:		du(x)	
		mdata[i][3] = ddu(x);	//			3:		ddu(x)	
		mdata[i][4] = p(x);		//			4:		p(x)	
		mdata[i][5] = q(x);		//			5:		q(x)	
		mdata[i][6] = f(x);		//			6:		f(x)	
		mdata[i][7] = fc(x);	//			7:		fc(x)	
	}
	printcsv("init.csv", mdata, "x;u(x);du(x);ddu(x);p(x);q(x);f(x);fc(x)");
	auto mtdata = transposed(mdata);

	// ------------------------- CONST SOLUTION -----------------------------
	// ------------------------------ |IV| ----------------------------------

	cout << "\n[CONST SOLUTION]\n";
	auto mac = createMatr(pc, qc, a, b, h);
	auto vbc = createRightVect(fc, a, b, h, ifunc);
	cout << format("[IV]: SLE param:\n\tsize: {}\n", mac.size());
	
	auto [uxc, rnc, itc] = SLESolver(mac, vbc, mtdata[1], S_EPS);
	cout << format("[IV]: Solver param:\n\trn: {:e}\n\tit: {}\n", rnc, itc);

	printcsv("test_const_4.csv", mtdata[0], uxc, "x;uhc");

	// ------------------------------- |V| ----------------------------------

	size_t steps = 100;
	vect rh(steps, 0);
	vect vh(steps, 0);

	auto start = system_clock::now();
	auto temp = start;
	auto current = start;

	for (int i = 0; i < steps; i++)
	{
		h = 1.0 / (10 * (i + 1));
		n = (b - a) / h - 1;
		h = (b - a) / n;

		auto vx = createXVect(a, b, n);
		auto vu = createVectFunc(u, vx);

		mac = createMatr(pc, qc, a, b, h);
		vbc = createRightVect(fc, a, b, h, ifunc);
		auto [uxc1, rnc1, itc1] = SLESolver(mac, vbc, vu, 1e-10);

		vh[i] = h;
		rh[i] = sqrt(h) * norm2(vu - uxc1);

		current = system_clock::now();
		if (duration_cast<seconds>(current - temp).count() > 5)	// 10s timer
		{
			temp = current;
			auto elapsed = current - start;
			auto elapsedm = duration_cast<minutes>(elapsed);
			auto elapseds = duration_cast<seconds>(elapsed - elapsedm);

			auto per = (float)i / steps * 100.f;
			cout << format("[V]: Computing... [{}/{}] - ({:.2f}%) - {:2}:{:2}\n", i, steps, per, elapsedm.count(), elapseds.count());
		}
	}
	printcsv("residual_const_5.csv", vh, rh, "h;rhc");

	// ------------------------- FUNC SOLUTION ------------------------------
	// ------------------------------ |IV| ----------------------------------

	h = 0.1;
	n = (b - a) / h - 1;
	h = (b - a) / n;

	cout << "\n[FUNC SOLUTION]\n";
	auto ma = createMatr(p, q, a, b, h);
	auto vb = createRightVect(f, a, b, h, ifunc);
	cout << format("[IV]: SLE param:\n\tsize: {}\n", ma.size());

	auto [ux, rn, it] = SLESolver(ma, vb, mtdata[1], S_EPS);
	cout << format("[IV]: Solver param:\n\trn: {:e}\n\tit: {}\n", rn, it);

	printcsv("test_func_4.csv", mtdata[0], ux, "x;uhf");

	// ------------------------------- |V| ----------------------------------

	steps = 100;
	rh.resize(steps, n);
	vh.resize(steps, n);

	start = system_clock::now();
	temp = start;
	current = start;

	for (int i = 0; i < steps; i++)
	{
		h = 1.0 / (10 * (i + 1));
		n = (b - a) / h - 1;
		h = (b - a) / n;

		auto vx = createXVect(a, b, n);
		auto vu = createVectFunc(u, vx);

		mac = createMatr(pc, qc, a, b, h);
		vbc = createRightVect(fc, a, b, h, ifunc);
		auto [ux1, rn1, it1] = SLESolver(ma, vb, vu, 1e-10);

		vh[i] = h;
		rh[i] = sqrt(h) * norm2(vu - ux1);

		current = system_clock::now();
		if (duration_cast<seconds>(current - temp).count() > 5)	// 10s timer
		{
			temp = current;
			auto elapsed = current - start;
			auto elapsedm = duration_cast<minutes>(elapsed);
			auto elapseds = duration_cast<seconds>(elapsed - elapsedm);

			auto per = (float)i / steps * 100.f;
			cout << format("[V]: Computing... [{}/{}] - ({:.2f}%) - {:2}:{:2}\n", i, steps, per, elapsedm.count(), elapseds.count());
		}
	}
	printcsv("residual_func_5.csv", vh, rh, "h;rhf");

	return 0;
}

int lab2()
{
	GLOBAL_PATH = "lab2/data/quad/";
	cout << "[USING QUAD]";
	lab2(quad);

	GLOBAL_PATH = "lab2/data/simp/";
	cout << "\n[USING SIMP]";
	lab2(simp);

	

	 

	return 0;
}
