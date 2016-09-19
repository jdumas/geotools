////////////////////////////////////////////////////////////////////////////////
#include "poisson_disk.h"
#include "random.h"
#include "vec.h"
// -----------------------------------------------------------------------------
#include <iostream>
#include <array>
#include <cmath>
#include <cassert>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

template<typename real_t, size_t n>
class Grid {

public:
	typedef std::array<real_t, n> vXr;
	typedef std::array<int, n> vXi;

private:
	const real_t m_MinSqDist;
	const real_t m_CellSize;
	const vXr m_Extent;
	const vXi m_Side;
	const vXi m_Coeff;
	std::vector<vXr> m_Content;

public:
	Grid(real_t r, vXr extent)
		: m_MinSqDist(r * r)
		, m_CellSize(r / std::sqrt(n))
		, m_Extent(extent)
		, m_Side(computeSide(r, extent))
		, m_Coeff(computeCoef(m_Side))
		, m_Content(computeProd(m_Side), Vec::constant<real_t, n>(-1))
	{ };

private:
	static vXi computeSide(real_t r, vXr extent) {
		vXi side;
		for (size_t i = 0; i < n; ++i) {
			side[i] = std::floor(1 + extent[i] * std::sqrt(n) / r);
		}
		return side;
	}

	static vXi computeCoef(vXi side) {
		vXi coef;
		coef[0] = 1;
		for (size_t i = 1; i < n; ++i) {
			coef[i] = coef[i - 1] * side[i - 1];
		}
		return coef;
	}

	static int computeProd(vXi side) {
		int m = 1;
		for (int c : side) { m *= c; }
		return m;
	}

private:
	int toGridIndex(vXr p) {
		int u = 0;
		for (size_t i = 0; i < n; ++i) {
			assert(p[i] >= 0 && p[i] <= m_Extent[i]);
			u += m_Coeff[i] * std::floor(p[i] / m_CellSize);
		}
		assert(u >= 0 && u < (int) m_Content.size());
		return u;
	}

	vXi toGridVect(vXr p) {
		vXi u;
		for (size_t i = 0; i < n; ++i) {
			u[i] = std::floor(p[i] / m_CellSize);
		}
		return u;
	}

	bool recurseIsOccupied(vXr p, vXi u, int v, int i) {
		if (i == n) {
			if (m_Content[v][0] < 0) {
				return false; // Cell is empty
			} else {
				return Vec::sqDistance(p, m_Content[v]) <= m_MinSqDist;
			}
		} else {
			const int d = static_cast<int>(std::ceil(std::sqrt(n)));

			for (int j = -d; j <= d; ++j) {
				if (u[i] + j < 0 || u[i] + j >= m_Side[i]) {
					// Don't recurse
				} else if (recurseIsOccupied(p, u, v + j * m_Coeff[i], i + 1)) {
					return true; // Stop recursion
				}
			}
			return false;
		}
	}

public:
	void insertInitPoint(vXr p) {
		int u = toGridIndex(p);
		m_Content[u] = p;
	}

	void insertPoint(vXr p) {
		int u = toGridIndex(p);
		assert(m_Content[u][0] < 0);
		m_Content[u] = p;
	}

	bool isNeighborhoodOccupied(vXr p) {
		return recurseIsOccupied(p, toGridVect(p), toGridIndex(p), 0);
	}
};

////////////////////////////////////////////////////////////////////////////////

template<typename real_t, size_t n>
void PoissonSampling<real_t, n>::box(
	int maxAttempts, std::vector<vXr> &result) const
{
	domain(maxAttempts, PoissonSampling<real_t, n>::Domain(), result);
}

// Implemented after:
// Fast Poisson disk sampling in arbitrary dimensions, R. Bridson, ACM SIGGRAPH 2007 Sketches Program.
template<typename real_t, size_t n>
void PoissonSampling<real_t, n>::domain(
	int maxAttempts,
	const PoissonSampling<real_t, n>::Domain &outputArea,
	std::vector<vXr> &result) const
{
	// Data structures
	const int maxDomainTrials = outputArea.maxTrials();
	std::vector<vXr> active;
	Grid<real_t, n> grid(m_MinDist, m_Extent);

	// Initialization
	if (not result.empty()) {
		// Update containers
		for (vXr p : result) {
			active.push_back(p);
			grid.insertInitPoint(p);
		}
	} else {
		// Start with a random initial point
		vXr firstPoint;
		int j;
		for (j = 0; j < maxDomainTrials; ++j) {
			for (size_t i = 0; i < n; ++i) {
				firstPoint[i] = Random::get<real_t>() * m_Extent[i];
			}
			if (outputArea.contains(firstPoint, m_Extent)) { break; }
		}
		if (j == maxDomainTrials) { return; }

		// Update containers
		result.push_back(firstPoint);
		active.push_back(firstPoint);
		grid.insertPoint(firstPoint);
	}

	// Main loop
	while (not active.empty()) {
		int selectedIndex = Random::uniform_int<int>(0, active.size() - 1);
		const vXr currentPoint = active[selectedIndex];
		int i;

		for (i = 0; i < maxAttempts; ++i) {
			vXr newPoint;
			int j;

			// Try to find a point both in the domain and the annulus (r, 2*r)
			for (j = 0; j < maxDomainTrials; ++j) {
				newPoint = currentPoint + Random::annulus<real_t, n>(m_MinDist, 2 * m_MinDist);
				if (outputArea.contains(newPoint, m_Extent)) { break; }
			}

			if (j == maxDomainTrials) {
				i = maxAttempts;
				break;
			} else if (not grid.isNeighborhoodOccupied(newPoint)) {
				result.push_back(newPoint);
				active.push_back(newPoint);
				grid.insertPoint(newPoint);
			}
		}

		if (i == maxAttempts) {
			// Drop the selected point
			std::swap(active[selectedIndex], active.back());
			active.pop_back();
		}
	}
}

// -----------------------------------------------------------------------------

template<typename real_t, size_t n>
void PoissonSampling<real_t, n>::contour(
	const std::vector<vXr> &poly, std::vector<vXr> &result) const
{
	Grid<real_t, n> grid(m_MinDist, m_Extent);

	for (int i = 0; i < (int) poly.size(); ++i) {
		vXr a = poly[i];
		vXr b = poly[(i + 1) % poly.size()];
		real_t l = distance(a, b);

		real_t minCoef = 0;
		real_t maxCoef = std::min(l, m_MinDist);
		int unsuccessfulTimes = 0;

		while (minCoef < l && unsuccessfulTimes < maxContourTrials) {
			real_t s = Random::uniform_real<real_t>(minCoef, maxCoef);
			vXr p = s/l * a + (1 - s/l) * b;
			if (not grid.isNeighborhoodOccupied(p)) {
				unsuccessfulTimes = 0;
				result.push_back(p);
				grid.insertPoint(p);
				minCoef = s + m_MinDist;
				maxCoef = std::min(l, s + 2 * m_MinDist);
			} else {
				++unsuccessfulTimes;
				if (unsuccessfulTimes >= maxContourTrials) {
					minCoef = maxCoef;
					maxCoef = std::min(l, 2 * maxCoef);
					unsuccessfulTimes = 0;
				}
			}
		}
	}

}

// -----------------------------------------------------------------------------

template<typename real_t, size_t n>
void PoissonSampling<real_t, n>::subset(
	const std::vector<vXr> &soup, std::vector<vXr> &result) const
{
	Grid<real_t, n> grid(m_MinDist, m_Extent);
	std::vector<vXr> active(soup);

	// Initialization
	for (vXr p : result) {
		grid.insertInitPoint(p);
	}

	while (not active.empty()) {
		int selectedIndex = Random::uniform_int<int>(0, active.size() - 1);
		const vXr newPoint = active[selectedIndex];

		if (not grid.isNeighborhoodOccupied(newPoint)) {
			result.push_back(newPoint);
			grid.insertPoint(newPoint);
		}

		// Drop the selected point
		std::swap(active[selectedIndex], active.back());
		active.pop_back();
	}
}

// -----------------------------------------------------------------------------

template<typename real_t, size_t n>
void PoissonSampling<real_t, n>::naive(
	std::vector<vXr> &result) const
{
	Domain defaultDomain;
	const int maxDomainTrials = defaultDomain.maxTrials();
	const real_t minSqDist = m_MinDist * m_MinDist;

	int unsuccessfulTimes = 0;
	while (unsuccessfulTimes < maxDomainTrials) {
		// Generate new random point
		vXr newPoint;
		for (size_t i = 0; i < n; ++i) {
			newPoint[i] = Random::get<real_t>() * m_Extent[i];
		}
		if (not defaultDomain.contains(newPoint, m_Extent)) {
			++unsuccessfulTimes;
			continue;
		}
		// Test against existing neighbors
		bool flag = false;
		for (vXr q : result) {
			if (Vec::sqDistance(newPoint, q) <= minSqDist) {
				flag = true;
				break;
			}
		}
		if (flag == false) {
			unsuccessfulTimes = 0;
			result.push_back(newPoint);
		} else {
			++unsuccessfulTimes;
		}
	}
}
