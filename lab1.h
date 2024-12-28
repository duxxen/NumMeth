#pragma once
#include "mathex.h"

vect_t SeidelSolve(const matr_t& A, const vect_t& b, double eps = 1e-10);
vect_t SeidelSolve(const matr_t& A, const vect_t& b, const vect_t& x0, double eps = 1e-10);
vect_t SeidelSolve(const matr_t& A, const vect_t& b, size_t maxit, double eps = 1e-10);
vect_t SeidelSolve(const matr_t& A, const vect_t& b, const vect_t& x0, size_t maxit, double eps = 1e-10);

