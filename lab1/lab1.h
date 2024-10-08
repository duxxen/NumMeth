#pragma once
#include "../utils/mathex.h"
#include "../utils/definitions.hpp"

// Seidel Method accuracy
#define SM_EPS 1e-6

// converts { matr A, vect b } -> { matr C, vect d } for GaussMethod
pair<matr, vect> toGaussMatrixRecurrent(const matr& A, const vect& b, size_t in);
vect GaussMethod(const matr& A, const vect& b);

// converts { matr A, vect b } -> { matr C, vect d } for SeidelMethod
pair<matr, vect> toSeidelMatrix(const matr& A, const vect& b);

// returns { vect x, vect r, double rn, size_t it } solution of SLE(matr A, vect b) using approx vect(x) with fixed accuracy(eps) or iterations(q)
tuple<vect, vect, double, size_t> SeidelMethod(const matr& A, const vect& b, const vect& x0, long double eps = SM_EPS, int q = -1);

// expands { vect x } -> { vect alpha, matr phi } according to [phi] basis
pair<vect, matr> basisExpansion(const vect& x);

// creates { matr A } thermal conductivity by size(n), grid(h)
matr createThermalMatrix(size_t n, double h);

int lab1();
