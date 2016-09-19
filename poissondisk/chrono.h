#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>
#include <vector>
#include <algorithm>
////////////////////////////////////////////////////////////////////////////////

// Counter with automatic time display
struct Chrono {
	std::string                                  m_Name;
	std::chrono::steady_clock::time_point        m_TimeStart;
	std::vector<std::pair<std::string, double> > m_Scores;

	// Public time point typedef
	typedef std::chrono::steady_clock::time_point TimePoint;

	// Get current time
	static TimePoint now() { return std::chrono::steady_clock::now(); }

	// Init internal counter
	Chrono(const std::string &s)
		: m_Name(s)
		, m_TimeStart(std::chrono::steady_clock::now())
	{ }

	// Start new subchrono
	void tic(const std::string &subname = "") {
		m_TimeStart = std::chrono::steady_clock::now();
	}

	// Stop last subchrono
	void toc() {
		const double t = Chrono::getElapsedTime(m_TimeStart);
		std::cout << "[" << m_Name << "] " << t << std::endl;
	}

	// Compute time difference
	template<typename clock_t>
	static double getElapsedTime(std::chrono::time_point<clock_t> &t1) {
		using namespace std::chrono;
		auto t2 = clock_t::now();
		auto time_span = duration_cast<duration<double> >(t2 - t1);
		t1 = t2;
		return time_span.count();
	}
};

////////////////////////////////////////////////////////////////////////////////
