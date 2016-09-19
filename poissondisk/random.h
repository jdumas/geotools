#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <random>
#include <chrono>
#include <iostream>
#include <array>
#include "vec.h"
////////////////////////////////////////////////////////////////////////////////

class Random {
public:
	static std::default_random_engine generator;

	static void reset() {
		generator = std::default_random_engine();
	}

	static void initSeed() {
		auto time = std::chrono::high_resolution_clock::now();
		auto count = time.time_since_epoch().count();
		generator.seed(count);
		std::cout << "Random generator initialized with seed: [==[ " << count << " ]==]" << std::endl;
	}

	template<typename T>
	static T uniform_int(T a, T b) {
		std::uniform_int_distribution<T> distribution(a, b);
		return distribution(generator);
	}

	template<typename T>
	static T uniform_real(T a, T b) {
		std::uniform_real_distribution<T> distribution(a, b);
		return distribution(generator);
	}

	template<typename T>
	static T get() {
		return uniform_real<T>(0, 1);
	}

	template<typename T, int n>
	static std::array<T, n> annulus(T r1, T r2) {
		std::array<T, n> p;
		for (int i = 0; i < n; ++i) {
			p[i] = uniform_real<T>(-1, 1);
		}
		T r = uniform_real<T>(r1, r2);
		return r * (p / std::sqrt(Vec::sqLength(p)));
	}
};
