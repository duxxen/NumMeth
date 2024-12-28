#pragma once
#include "definitions.hpp"

typedef std::vector<double> vect_t;
typedef std::vector<vect_t> matr_t;

vect_t operator+(const vect_t& op1, const vect_t& op2);
vect_t& operator+=(vect_t& rop, const vect_t& op);

vect_t operator-(const vect_t& op1, const vect_t& op2);
vect_t& operator-=(vect_t& rop, const vect_t& op);

vect_t dot(const vect_t& op1, const vect_t& op2);
vect_t dot(const matr_t& mtr, const vect_t& vct);
vect_t& dotadd(vect_t& rvct, const matr_t& mtr, const vect_t& vct);

double norm2(const vect_t& vct);