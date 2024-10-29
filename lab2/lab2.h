#pragma once
#include "../utils/mathex.h"
#include "../utils/definitions.hpp"

#define Q_EPS 1e-6

typedef function<long double(func, long double, long double, long double)> IntegralFunc;

// returns integral of func(f) at segment[a, b] with accuracy(eps)
long double quad(const func& f, long double a, long double b, long double eps = Q_EPS);

long double Simpson(const func& f, long double a, long double b, long double h);

// returns { matr A, vect b } of thermal conductivity SLE based on { func(p), func(q), func(f), grid(h), range[a, b] } using integral(ifunc) with accuracy(eps)
pair<matr, vect> createSLE(const func& p, const func& q, const func& f, long double h, long double a, double b, long double eps = Q_EPS, IntegralFunc ifunc = quad);

// returns { matr A, vect b } of thermal conductivity SLE based on { const(p), const(q), func(f), grid(h), range[a, b] } using integral(ifunc) with accuracy(eps)
pair<matr, vect> createSLE(long double p, long double q, const func& f, long double h, long double a, double b, long double eps = Q_EPS, IntegralFunc ifunc = quad);

vect getFuncVect(const func& u, long double a, long double b, long double h);

int lab2();