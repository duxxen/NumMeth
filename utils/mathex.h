#pragma once
#include <cmath>
#include <vector>
#include <functional>
#include <cassert>

typedef std::vector<double>					vect;
typedef std::vector<std::vector<double>>	matr;
typedef std::function<double(double)>		func;

const static double M_PI = 3.1415926535897932384626433832795028841971;

int randint(int min = 0, int max = RAND_MAX);
float randfloat(float min = 0, float = 1);

vect randintVect(size_t n, int min = 0, int max = RAND_MAX);
vect randfloatVect(size_t n, float min = 0, float max = 1);

double norm2(const vect& v);
double scalar(const vect& v1, const vect& v2);

vect operator+(const vect& v1, const vect& v2);
vect operator-(const vect& v1, const vect& v2);
vect operator*(double c, const vect& v);
vect operator*(const vect& v, double c);
vect operator/(const vect& v, double c);

matr transposed(const matr& m);

matr operator+(const matr& v1, const matr& v2);
matr operator-(const matr& v1, const matr& v2);
matr operator*(double c, const matr& v);
matr operator*(const matr& v, double c);
matr operator/(const matr& v, double c);

matr operator*(const matr& m1, const matr& m2);
vect operator*(const matr& m, const vect& v);
vect operator*(const vect& v, const matr& m);