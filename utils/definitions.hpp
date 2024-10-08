#pragma once
#include <fstream>
#include <cassert>
#include <vector>
#include <functional>
#include <iostream>

using std::cout;
using std::endl;
using std::pair;
using std::tuple;
using std::vector;
using std::function;

template <typename T>
std::ostream& operator <<(std::ostream& out, const std::vector<T>& v)
{
	out << "{ ";
	for (int i = 0; i < v.size(); i++)
		out << v[i] << " ";
	out << "}";

	return out;
}

template <typename T>
std::ostream& operator <<(std::ostream& out, std::vector<std::vector<T>>& m)
{
	for (const auto& v : m)
		out << v << "\n";

	return out;
}