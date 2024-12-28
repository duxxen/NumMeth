#include <tuple>
#include <vector>
#include <cassert>
#include <fstream>
#include <iostream>

template <typename T>
std::ostream& operator <<(std::ostream& out, const std::vector<T>& vct)
{
	out << "[";
	for (const auto& elm : vct)
		out << " " << elm;
	out << " ]";

	return out;
}

template <typename T>
std::ostream& operator <<(std::ostream& out, std::vector<std::vector<T>>& m)
{
	for (const auto& v : m)
		out << v << "\n";

	return out;
}