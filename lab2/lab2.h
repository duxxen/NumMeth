#pragma once
#include "../utils/mathex.h"
#include "../utils/definitions.hpp"

extern int K, N, C;
extern double a, b;
extern double cp, cq;

double fp(double x);
double fq(double x);
double dfp(double x);

double u(double x);
double du(double x);
double ddu(double x);

double fc(double x);
double ff(double x);

double quad(func f, double a, double b, double q);
double simp(func f, double a, double b, double q);

pair<matr, vect> createCSLE(double h);
pair<matr, vect> createFSLE(double h);

void lab2();