#include "../utils/definitions.hpp"
#include "../utils/mathex.h"

typedef function<double(const func&, double, double, long double)> IntegralFunc;

double quad(const func& f, double a, double b, long double h);
double simp(const func& f, double a, double b, long double h);

vect createXVect(double a, double b, size_t n);
vect createVectFunc(const func& f, const vect& x);

matr createMatr(double p, double q, double a, double b, double h);
matr createMatr(const func& p, const func& q, double a, double b, double h, const IntegralFunc& ifunc = quad);
vect createRightVect(const func& f, double a, double b, double h, const IntegralFunc& ifunc = quad);

tuple<vect, double, size_t> SLESolver(const matr& a, const vect& b, const vect& x0, long double eps);

int printcsv(const string& path, const matr& mdata, const string& head = "", const string& del = ";");
int printcsv(const string& path, const vect& vx, const vect& vy, const string& head = "", const string& del = ";");

extern string GLOBAL_PATH;
int lab2(const IntegralFunc& ifunc);
int lab2();
