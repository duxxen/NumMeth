#include "lab1.h"

pair<matr, vect> toGaussMatrixRecurrent(const matr& A, const vect& b, size_t in)
{
	size_t n = A.size();
	assert(
		b.size() == n &&
		A.back().size() == n &&
		0 <= in && in < n
	);

	matr C = A;
	vect d = b;

	if (in == 0)
	{
		auto maxi = 0;
		for (auto i = 1; i < n; i++)
			if (C[i][0] > C[maxi][0])
				maxi = i;

		std::swap(C[maxi], C[0]);
		std::swap(d[maxi], d[0]);
	}

	for (auto i = in; i < n; i++)
	{
		auto temp = A[i][in];
		if (temp == 0)
			continue;

		for (auto j = in; j < n; j++)
			C[i][j] /= temp;
		d[i] /= temp;

		if (i > in)
		{
			C[i] = C[i] - C[in];
			d[i] = d[i] - d[in];
		}
	}

	if (in < n - 1)		
		return toGaussMatrixRecurrent(C, d, in + 1);
	return { C, d };
}

vect GaussMethod(const matr& A, const vect& b)
{
	size_t n = A.size();
	vect x(n, 0);
	auto [C, d] = toGaussMatrixRecurrent(A, b, 0);

	for (int i = n - 1; i >= 0; i--)
	{
		auto sum = d[i];
		for (auto j = i + 1; j < n; j++)
			sum -= C[i][j] * x[j];
		x[i] = sum;
	}

	return x;
}

pair<matr, vect> toSeidelMatrix(const matr& A, const vect& b)
{
	size_t n = A.size();
	assert(
		b.size() == n &&
		A.back().size() == n
	);

	matr C(n, vect(n, 0));
	vect d(n, 0);

	for (auto i = 0; i < n; i++)
	{
		d[i] = b[i] / A[i][i];
		for (auto j = 0; j < n; j++)
			if (i != j) C[i][j] = -A[i][j] / A[i][i];
	}

	return { C, d };
}

tuple<vect, vect, double, size_t> SeidelMethod(const matr& A, const vect& b, const vect& x0, long double eps, int q)
{
	size_t n = x0.size();
	size_t i = 0;
	vect x = x0;
	auto [C, d] = toSeidelMatrix(A, b);

	vect xp(n, 0);
	vect r(n, 0);
	double rn = 0;
	double en = 0;

	do
	{
		xp = x;
		x = C * xp + d;
		r = A * x - b;
		rn = norm2(r);
		en = norm2(x - xp);

		if (i == q)
			break;
		i++;
	} while (rn > eps);

	return { x, r, rn, i };
}

pair<vect, matr> basisExpansion(const vect& x)
{
	size_t n = x.size();
	double h = 1.0 / (n + 1);
	matr psi(n, vect(n, 0));
	matr phi(n, vect(n, 0));
	vect alpha(n, 0);

	for (auto k = 0; k < n; k++)
		for (auto j = 0; j < n; j++)
			psi[k][j] = sin(M_PI * (k + 1) * (j + 1) * h);

	for (auto k = 0; k < n; k++)
	{
		phi[k] = psi[k] / norm2(psi[k]);
		double sum = 0;
		for (auto i = 0; i < n; i++)
			sum += x[i] * phi[k][i];
		alpha[k] = sum;
	}

	return { alpha, phi };
}

matr createThermalMatrix(size_t n, double h)
{
	matr A(n, vect(n, 0));
	auto value = 1.0 / pow(h, 2);

	for (auto i = 0; i < n; i++)
	{
		A[i][i] = 2 * value;
		A[i][i < n - 1 ? i + 1 : i - 1] = -value;
		A[i][i > 0 ? i - 1 : i + 1] = -value;
	}

	return A;
}

int lab1()
{
	srand(time(0));

	size_t q2 = 20;
	size_t n2 = randint(4, 20);
	vect b2 = randintVect(n2, -5, 5);
	vect x02 = randintVect(n2, -5, 5);
	matr A2(n2, vect(n2, 0));

	for (int i = 0; i < n2; i++)
	{
		A2[i] = randintVect(n2, -3, 3);
		A2[i][i] = 10 * randint(1, 5);
	}

	cout << "------------------- 2) Seidel m-d realisation -------------------\n";
	cout << "Inputs:\n";
	cout << "A = \n" << A2 << endl;
	cout << "b = " << b2 << endl;
	cout << "x0 = " << x02 << endl;
	cout << "-------------------------\n";

	auto [x2, r2, rn2, i2] = SeidelMethod(A2, b2, x02);

	cout << "Outputs:\n";
	cout << "x = " << x2 << endl;
	cout << "r = " << r2 << endl;
	cout << "rn = " << rn2 << endl;
	cout << "i = " << i2 << endl;
	cout << "-------------------------\n";

	// ------------------------------------------------------------------

	size_t n3 = 4;
	vect x3 = randintVect(n3, -3, 3);

	cout << "------------------- 3) Vector expansion -------------------\n";
	cout << "Inputs:\n";
	cout << "x = " << x3 << endl;
	cout << "-------------------------\n";

	auto [alpha3, phi3] = basisExpansion(x3);
	vect sum3(4, 0);
	for (int i = 0; i < 4; i++)
		sum3 = sum3 + phi3[i] * alpha3[i];

	cout << "Outputs:\n";
	cout << "x = " << sum3 << endl;
	cout << "alpha = " << alpha3 << endl;
	cout << "phi = \n" << phi3 << endl;
	cout << "-------------------------\n";

	// ------------------------------------------------------------------

	size_t n4 = 16;
	matr A4 = createThermalMatrix(n4, 0.2);
	vect b4 = randintVect(n4, -50, 50);
	vect x04 = b4;

	cout << "------------------- 4) Seidel m-d with matrix-A -------------------\n";
	cout << "Inputs:\n";
	cout << "A = \n" << A4 << endl;
	cout << "b = " << b4 << endl;
	cout << "x0 = " << x04 << endl;
	cout << "-------------------------\n";

	auto [x4, r4, rn4, i4] = SeidelMethod(A4, b4, x04);

	cout << "Outputs:\n";
	cout << "x = " << x4 << endl;
	cout << "r = " << r4 << endl;
	cout << "rn = " << rn4 << endl;
	cout << "i = " << i4 << endl;
	cout << "-------------------------\n";

	// ------------------------------------------------------------------

	size_t n5 = n4;
	matr A5 = A4;
	vect b5 = b4;
	vect x05 = x04;
	vector<vect> e5(21);
	vect p5(21);
	vect ex5 = GaussMethod(A5, b5);

	for (auto i = 0; i < 21; i++)
		p5[i] = 10 * i;

	cout << "------------------- 5) Error vector -------------------\n";
	cout << "Inputs:\n";
	cout << "x = " << ex5 << " <- (exact solution)" << endl;
	cout << "p = " << p5 << endl;

	for (auto i = 0; i < e5.size(); i++)
	{
		auto [x5, r5, rn5, it5] = SeidelMethod(A5, b5, x05, SM_EPS, p5[i]);
		e5[i] = ex5 - x5;
	}

	cout << "Outputs:\n";
	cout << "e = \n" << e5 << endl;
	cout << "-------------------------\n";

	// ------------------------------------------------------------------

	cout << "------------------- 6) Error vector expansion -------------------\n";
	cout << "Inputs:\n";
	cout << "e = \n" << e5 << endl;
	cout << "-------------------------\n\n\n";

	size_t n6 = n5;
	vector<vect> e6 = e5;
	vector<vect> alpha6(e6.size());
	vector<matr> phi6(e6.size());
	for (auto i = 0; i < e6.size(); i++)
	{
		auto res = basisExpansion(e6[i]);
		alpha6[i] = res.first;
		phi6[i] = res.second;
	}

	cout << "\n\nOutputs:\n";
	cout << "e = \n" << e6 << endl;
	cout << "alpha_p = \n" << alpha6 << endl;
	cout << "-------------------------\n";

	std::ofstream fout("lab1/abs_alpha.csv");
	if (!fout.is_open())
	{
		cout << "Error occuried while opening file!\n";
		return 1;
	}

	for (int i = 0; i < n6 + 1; i++)
		fout << n6 << ";";
	fout << endl;

	for (int j = 0; j < 21; j++)
	{
		fout << p5[j] << ";";
		for (int i = 0; i < n6; i++)
		{

			fout << abs(alpha6[j][i]) << ";";
		}
		fout << endl;
	}
	fout.close();

	return 0;
}