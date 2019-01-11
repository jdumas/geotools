#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <array>
////////////////////////////////////////////////////////////////////////////////

/*
 * Generate Poisson disk distribution in n-dimensional space.
 *
 * Reference paper:
 *   Bridson, Robert. "Fast Poisson disk sampling in arbitrary dimensions."
 *   SIGGRAPH sketches. 2007.
 */

template<typename real_t, size_t n>
class PoissonSampling {

public:
	typedef std::array<real_t, n> vXr;

public:
	struct Domain {
		virtual ~Domain() = default;

		// The default domain is a n-dimensional box
		virtual bool contains(vXr p, vXr extent) const {
			for (unsigned i = 0; i < n; ++i) {
				if (p[i] < 0 || p[i] > extent[i]) {
					return false;
				}
			}
			return true;
		}

		// Max number of trials when trying to find a point in the domain
		virtual int maxTrials() const {
			return 2000;
		}
	};

private:
	static constexpr int maxContourTrials = 50;

private:
	const real_t m_MinDist;
	const vXr m_Extent;

public:
	PoissonSampling(real_t minDist, vXr extent)
		: m_MinDist(minDist), m_Extent(extent)
	{ }

	// Typically k = 30
	void box(int k, std::vector<vXr> &result) const;
	void domain(int k, const Domain &d, std::vector<vXr> &result) const;
	void contour(const std::vector<vXr> &poly, std::vector<vXr> &result) const;
	void subset(const std::vector<vXr> &soup, std::vector<vXr> &result) const;
	void naive(std::vector<vXr> &result) const;
};
