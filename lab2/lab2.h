#pragma once
#include "../utils/mathex.h"
#include "../utils/definitions.hpp"

#define Q_EPS 1e-6

typedef function<double(func, double, double, long double)> IntegralFunc;

// returns integral of func(f) at segment[a, b] with accuracy(eps)
double quad(const func& f, double a, double b, long double eps = Q_EPS);

double Simpson(const func& f, double a, double b, long double h);

// returns { matr A, vect b } of thermal conductivity SLE based on { func(p), func(q), func(f), grid(h), range[a, b] } using integral(ifunc) with accuracy(eps)
pair<matr, vect> createSLE(const func& p, const func& q, const func& f, double h, double a, double b, long double eps = Q_EPS, IntegralFunc ifunc = quad);

int lab2();