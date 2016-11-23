#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <array>
#include <iterator>
////////////////////////////////////////////////////////////////////////////////

// STL-friendly fixed size containers

typedef std::array<unsigned char, 2> Vec2b;
typedef std::array<int,           2> Vec2i;
typedef std::array<unsigned,      2> Vec2u;
typedef std::array<float,         2> Vec2f;
typedef std::array<double,        2> Vec2d;

typedef std::array<unsigned char, 3> Vec3b;
typedef std::array<int,           3> Vec3i;
typedef std::array<unsigned,      3> Vec3u;
typedef std::array<float,         3> Vec3f;
typedef std::array<double,        3> Vec3d;

typedef std::array<unsigned char, 4> Vec4b;
typedef std::array<int,           4> Vec4i;
typedef std::array<unsigned,      4> Vec4u;
typedef std::array<float,         4> Vec4f;
typedef std::array<double,        4> Vec4d;

////////////////////////////////////////////////////////////////////////////////

template <typename Scalar, size_t Rows>
inline std::ostream& operator<<(std::ostream &out, std::array<Scalar, Rows> v) {
	out << "{";
	if (!v.empty()) {
		std::copy(v.begin(), v.end() - 1, std::ostream_iterator<Scalar>(out, "; "));
		out << v.back();
	}
	out << "}";
	return out;
}

// -----------------------------------------------------------------------------

template<typename T, size_t N>
std::array<T, N> operator +(std::array<T, N> x, std::array<T, N> y) {
	std::array<T, N> res;
	for (size_t i = 0; i < N; ++i) {
		res[i] = x[i] + y[i];
	}
	return res;
}

template<typename T, size_t N>
std::array<T, N> operator -(std::array<T, N> x, std::array<T, N> y) {
	std::array<T, N> res;
	for (size_t i = 0; i < N; ++i) {
		res[i] = x[i] - y[i];
	}
	return res;
}

template<typename T, size_t N>
std::array<T, N> operator -(std::array<T, N> y) {
	std::array<T, N> res;
	for (size_t i = 0; i < N; ++i) {
		res[i] = - y[i];
	}
	return res;
}

template<typename T, size_t N>
std::array<T, N> operator /(std::array<T, N> x, T s) {
	std::array<T, N> res;
	for (size_t i = 0; i < N; ++i) {
		res[i] = x[i] / s;
	}
	return res;
}

template<typename T, size_t N>
std::array<T, N> operator *(T s, std::array<T, N> x) {
	std::array<T, N> res;
	for (size_t i = 0; i < N; ++i) {
		res[i] = s * x[i];
	}
	return res;
}

// -----------------------------------------------------------------------------

namespace Vec {
	template<typename T, size_t N>
	std::array<T, N> constant(T c) {
		std::array<T, N> x;
		x.fill(c);
		return x;
	}

	template<typename T, size_t N>
	T prod(std::array<T, N> x) {
		T res(1);
		for (size_t i = 0; i < N; ++i) {
			res *= x[i];
		}
		return res;
	}

	template<typename T, size_t N>
	T sqLength(std::array<T, N> x) {
		T res(0);
		for (size_t i = 0; i < N; ++i) {
			res += x[i] * x[i];
		}
		return res;
	}

	template<typename T, size_t N>
	T sqDistance(std::array<T, N> x, std::array<T, N> y) {
		T res(0);
		for (size_t i = 0; i < N; ++i) {
			res += (x[i] - y[i]) * (x[i] - y[i]);
		}
		return res;
	}
}
