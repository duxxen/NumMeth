#include "../utils/definitions.hpp"
#include "../utils/mathex.h"

typedef function<double(const func&, double, double, long double)> IntegralFunc;

double quad(const func& f, double a, double b, long double eps);
double simp(const func& f, double a, double b, long double eps);

double xi(size_t i, double a, double h);
vect fv(const func& f, const vect& x);

matr createMatr(double p, double q, double a, double b, double h, const IntegralFunc& ifunc = quad);
matr createMatr(const func& p, const func& q, double a, double b, double h, const IntegralFunc& ifunc = quad);
vect createRightVect(const func& f, double a, double b, double h, const IntegralFunc& ifunc = quad);

tuple<vect, double, size_t> SLESolver(const matr& a, const vect& b, long double eps);

int printcsv(const string& path, const matr& mdata, const string& head = "", const string& del = ";");
int printcsv(const string& path, const vect& vx, const vect& vy, const string& head = "", const string& del = ";");

int lab2();
