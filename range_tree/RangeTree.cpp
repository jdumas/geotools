////////////////////////////////////////////////////////////////////////////////
#include "RangeTree.h"
#include "ThreadPool.h"
#include "Common.h"
// -----------------------------------------------------------------------------
#include <algorithm>
#include <numeric>
#include <cassert>
#include <set>
////////////////////////////////////////////////////////////////////////////////

typedef RangeTree::index_t index_t;

// TODO: Use abstract templated type "Callable" to make a generic rangetree that
// can work with a "const double *" as well as an Eigen matrix.

// -----------------------------------------------------------------------------

namespace {

	// Computes the number of bits necessary to write x (assumes x != 0)
	index_t nbits(index_t x) {
#ifdef WIN32
		index_t m = 0;
		while ((1 << m) <= x) { ++m; }
		return m;
#else
		return index_t(__builtin_clz(index_t(1)) - __builtin_clz(x) + 1);
#endif
	}

}

////////////////////////////////////////////////////////////////////////////////
// Forward declarations
////////////////////////////////////////////////////////////////////////////////

template<int MaxDim> class RangeTree2d;
template<int MaxDim> class RangeTree3d;

////////////////////////////////////////////////////////////////////////////////
// Common methods
////////////////////////////////////////////////////////////////////////////////

// Abstract class for internal range tree implementation
class RangeTreeInternal {
public:
	virtual ~RangeTreeInternal() = default;

	// Rebuild the search trees according to a new set of coordinates
	virtual void rebuildIndex (const double *points) = 0;

	// Computes the number of neighbors within a query box
	virtual index_t countPointsInBox (
		const double *points,
		const double *lower,
		const double *upper) const = 0;

	// Retrieve points within a query box (assumes output buffer has been allocated)
	virtual index_t * getPointsInBox (
		const double *points,
		const double *lower,
		const double *upper,
		index_t *neighbors) const = 0;

	// Retrieve points within a query box (assumes output buffer has been allocated)
	virtual void getPointsInBox (
		const double *points,
		const double *lower,
		const double *upper,
		std::vector<index_t> &neighbors) const = 0;

	// Computes the number of neighbors within a query sphere
	virtual index_t countPointsInSphere (
		const double *points,
		const double *lower,
		const double *upper,
		const double *center,
		double sq_dist) const = 0;

	// Retrieve points within a query sphere (assumes output buffer has been allocated)
	virtual index_t * getPointsInSphere (
		const double *points,
		const double *lower,
		const double *upper,
		const double *center,
		double sq_dist,
		index_t *neighbors) const = 0;

	// Retrieve points within a query sphere (assumes output buffer has been allocated)
	virtual void getPointsInSphere (
		const double *points,
		const double *lower,
		const double *upper,
		const double *center,
		double sq_dist,
		std::vector<index_t> &neighbors) const = 0;
};

////////////////////////////////////////////////////////////////////////////////

template<int Dim, int MaxDim>
class RangeTreeLeaves {
	// Allow 2d range-tree to acces to internal members of its child nodes
	friend class RangeTree2d<MaxDim>;
	friend class RangeTree3d<MaxDim>;

protected:
	// Leaves are indices sorted according to coordinate Dim
	std::vector<index_t> m_Leaves;

protected:
	// Default constructor
	RangeTreeLeaves () = default;

	// Base constructor that creates a list of indices of the given size
	RangeTreeLeaves (index_t nb_points)
		: m_Leaves(nb_points)
	{
		std::iota(m_Leaves.begin(), m_Leaves.end(), 0);
	}

	// Base constructor that fills the leaves with the given (sorted) indices
	RangeTreeLeaves (index_t nb_points, const index_t *indices)
		: m_Leaves(indices, indices + nb_points)
	{ }

public:
	// Returns the number of leaves in the tree
	size_t size() const { return m_Leaves.size(); }

protected:
	// Compute the range [first, last) of points within a query box
	std::pair<index_t, index_t> getRangeBox (
		const double *points,
		const double *lower,
		const double *upper) const
	{
		if (m_Leaves.empty()) { return std::make_pair(0, 0); }
		index_t low_left  = 0;
		index_t low_right = index_t(m_Leaves.size()) - 1;
		index_t upp_left  = 0;
		index_t upp_right = index_t(m_Leaves.size()) - 1;
		if (points[MaxDim*m_Leaves[low_right] + Dim] < lower[Dim]) {
			return std::make_pair(0, 0);
		} else if (points[MaxDim*m_Leaves[upp_left] + Dim] > upper[Dim]) {
			return std::make_pair(0, 0);
		} else {
			// Find lower bound
			if (points[MaxDim*m_Leaves[low_left] + Dim] >= lower[Dim]) {
				low_right = low_left;
			} else {
				while (low_right - low_left > 1) {
					index_t m = (low_left + low_right) / 2;
					if (points[MaxDim*m_Leaves[m] + Dim] < lower[Dim]) {
						low_left = m;
					} else {
						low_right = m;
					}
				}
			}

			// Find upper bound
			if (points[MaxDim*m_Leaves[upp_right] + Dim] <= upper[Dim]) {
				upp_left = upp_right;
			} else {
				while (upp_right - upp_left > 1) {
					index_t m = (upp_left + upp_right) / 2;
					if (points[MaxDim*m_Leaves[m] + Dim] <= upper[Dim]) {
						upp_left = m;
					} else {
						upp_right = m;
					}
				}
			}

			return std::make_pair(low_right, upp_left + 1);
		}
	}

};

////////////////////////////////////////////////////////////////////////////////
// Common method for intermediate range-trees
////////////////////////////////////////////////////////////////////////////////

template<int Dim, int MaxDim, typename DerivedTree>
class RangeTreeNodes : public RangeTreeInternal {
protected:
	// Access derived implementation (no virtual methods!)
	inline const DerivedTree * _impl() const {
		return static_cast<const DerivedTree *>(this);
	}

	// Returns the index of the leftmost leaf given an internal node
	index_t leftmostLeaf(index_t node) const {
		index_t level     = nbits(node);
		index_t nb_levels = nbits(index_t(_impl()->m_Leaves.size()));
		index_t mask      = (1u << nb_levels) - 1;
		return (node << (nb_levels - level + 1)) & mask;
	}

	// Returns the index past the rightmost leaf given an internal node
	index_t rightmostLeaf(index_t node) const {
		index_t level     = nbits(node);
		index_t nb_levels = nbits(index_t(_impl()->m_Leaves.size()));
		index_t mask      = (1u << nb_levels) - 1;
		return std::min(index_t(_impl()->m_Leaves.size()),
			((node << (nb_levels - level + 1)) & mask)
				+ (1u << (nb_levels - level + 1)));
	}

	// Test whether two given points lie within a given distance
	static bool distLessThan (const double *p, const double *q, double sq_dist) {
		double d = 0;
		for (int i = 0; i < MaxDim; ++i) {
			d += (p[i] - q[i]) * (p[i] - q[i]);
		}
		return d < sq_dist;
	}

public:
	// Computes the number of neighbors within a query box
	index_t countPointsInBox (
		const double *points,
		const double *lower,
		const double *upper) const override
	{
		index_t accu = 0;
		auto p = _impl()->getRangeBox(points, lower, upper);

		index_t nb_levels = nbits(index_t(_impl()->m_Leaves.size()));
		index_t mask      = (1u << nb_levels) - 1;

		for (index_t leaf = p.first; leaf < p.second;) {
			// Find largest subtree in range for which 'leaf' is the leftmost leaf
			index_t node   = leaf + (1u << nb_levels);
			index_t parent = node >> 1;
			while (parent > 1
				&& leftmostLeaf(parent) >= leaf
				&& rightmostLeaf(parent) <= p.second)
			{
				node = parent;
				parent >>= 1;
			}
			assert(rightmostLeaf(node) > leaf);
			leaf = rightmostLeaf(node);

			if ((node & mask) == node) {
				// Process internal node of the tree
				accu += _impl()->m_Nodes[node].countPointsInBox(
					points, lower, upper);
			} else {
				// Process leaf node of the tree
				node = node & mask;
				const double *x = points + MaxDim*_impl()->m_Leaves[node];
				bool ok = true;
				for (int i = 0; i < Dim; ++i) {
					ok = ok && (x[i] >= lower[i] && x[i] <= upper[i]);
				}
				if (ok) { ++accu; };
			}
		}
		return accu;
	}

	// Retrieve points within a query box (assumes output buffer has been allocated)
	index_t * getPointsInBox (
		const double *points,
		const double *lower,
		const double *upper,
		index_t *neighbors) const override
	{
		auto p = _impl()->getRangeBox(points, lower, upper);

		index_t nb_levels = nbits(index_t(_impl()->m_Leaves.size()));
		index_t mask      = (1u << nb_levels) - 1;

		for (index_t leaf = p.first; leaf < p.second;) {
			// Find largest subtree in range for which 'leaf' is the leftmost leaf
			index_t node   = leaf + (1u << nb_levels);
			index_t parent = node >> 1;
			while (parent > 1
				&& leftmostLeaf(parent) >= leaf
				&& rightmostLeaf(parent) <= p.second)
			{
				node = parent;
				parent >>= 1;
			}
			assert(rightmostLeaf(node) > leaf);
			leaf = rightmostLeaf(node);

			if ((node & mask) == node) {
				// Process internal node of the tree
				neighbors = _impl()->m_Nodes[node].getPointsInBox(
					points, lower, upper, neighbors);
			} else {
				// Process leaf node of the tree
				node = node & mask;
				const double *x = points + MaxDim*_impl()->m_Leaves[node];
				bool ok = true;
				for (int i = 0; i < Dim; ++i) {
					ok = ok && (x[i] >= lower[i] && x[i] <= upper[i]);
				}
				if (ok) {
					neighbors[0] = _impl()->m_Leaves[node];
					++neighbors;
				}
			}
		}
		return neighbors;
	}

	// Retrieve points within a query box (std::vector version)
	void getPointsInBox (
		const double *points,
		const double *lower,
		const double *upper,
		std::vector<index_t> &neighbors) const override
	{
		auto p = _impl()->getRangeBox(points, lower, upper);

		index_t nb_levels = nbits(index_t(_impl()->m_Leaves.size()));
		index_t mask      = (1u << nb_levels) - 1;

		for (index_t leaf = p.first; leaf < p.second;) {
			// Find largest subtree in range for which 'leaf' is the leftmost leaf
			index_t node   = leaf + (1u << nb_levels);
			index_t parent = node >> 1;
			while (parent > 1
				&& leftmostLeaf(parent) >= leaf
				&& rightmostLeaf(parent) <= p.second)
			{
				node = parent;
				parent >>= 1;
			}
			assert(rightmostLeaf(node) > leaf);
			leaf = rightmostLeaf(node);

			if ((node & mask) == node) {
				// Process internal node of the tree
				_impl()->m_Nodes[node].getPointsInBox(
					points, lower, upper, neighbors);
			} else {
				// Process leaf node of the tree
				node = node & mask;
				const double *x = points + MaxDim*_impl()->m_Leaves[node];
				bool ok = true;
				for (int i = 0; i < Dim; ++i) {
					ok = ok && (x[i] >= lower[i] && x[i] <= upper[i]);
				}
				if (ok) {
					neighbors.emplace_back(_impl()->m_Leaves[node]);
				}
			}
		}
	}

	// Computes the number of neighbors within a query sphere
	index_t countPointsInSphere (
		const double *points,
		const double *lower,
		const double *upper,
		const double *center,
		double sq_dist) const override
	{
		index_t accu = 0;
		auto p = _impl()->getRangeBox(points, lower, upper);

		index_t nb_levels = nbits(index_t(_impl()->m_Leaves.size()));
		index_t mask      = (1u << nb_levels) - 1;

		for (index_t leaf = p.first; leaf < p.second;) {
			// Find largest subtree in range for which 'leaf' is the leftmost leaf
			index_t node   = leaf + (1u << nb_levels);
			index_t parent = node >> 1;
			while (parent > 1
				&& leftmostLeaf(parent) >= leaf
				&& rightmostLeaf(parent) <= p.second)
			{
				node = parent;
				parent >>= 1;
			}
			assert(rightmostLeaf(node) > leaf);
			leaf = rightmostLeaf(node);

			if ((node & mask) == node) {
				// Process internal node of the tree
				accu += _impl()->m_Nodes[node].countPointsInSphere(
					points, lower, upper, center, sq_dist);
			} else {
				// Process leaf node of the tree
				node = node & mask;
				const double *x = points + MaxDim*_impl()->m_Leaves[node];
				if (distLessThan(center, x, sq_dist)) {
					++accu;
				}
			}
		}
		return accu;
	}

	// Retrieve points within a query sphere (assumes output buffer has been allocated)
	index_t * getPointsInSphere (
		const double *points,
		const double *lower,
		const double *upper,
		const double *center,
		double sq_dist,
		index_t *neighbors) const override
	{
		auto p = _impl()->getRangeBox(points, lower, upper);

		index_t nb_levels = nbits(index_t(_impl()->m_Leaves.size()));
		index_t mask      = (1u << nb_levels) - 1;

		for (index_t leaf = p.first; leaf < p.second;) {
			// Find largest subtree in range for which 'leaf' is the leftmost leaf
			index_t node   = leaf + (1u << nb_levels);
			index_t parent = node >> 1;
			while (parent > 1
				&& leftmostLeaf(parent) >= leaf
				&& rightmostLeaf(parent) <= p.second)
			{
				node = parent;
				parent >>= 1;
			}
			assert(rightmostLeaf(node) > leaf);
			leaf = rightmostLeaf(node);

			if ((node & mask) == node) {
				// Process internal node of the tree
				neighbors = _impl()->m_Nodes[node].getPointsInSphere(
					points, lower, upper, center, sq_dist, neighbors);
			} else {
				// Process leaf node of the tree
				node = node & mask;
				const double *x = points + MaxDim*_impl()->m_Leaves[node];
				if (distLessThan(center, x, sq_dist)) {
					neighbors[0] = _impl()->m_Leaves[node];
					++neighbors;
				}
			}
		}
		return neighbors;
	}

	// Retrieve points within a query sphere (std::vector version)
	void getPointsInSphere (
		const double *points,
		const double *lower,
		const double *upper,
		const double *center,
		double sq_dist,
		std::vector<index_t> &neighbors) const override
	{
		auto p = _impl()->getRangeBox(points, lower, upper);

		index_t nb_levels = nbits(index_t(_impl()->m_Leaves.size()));
		index_t mask      = (1u << nb_levels) - 1;

		for (index_t leaf = p.first; leaf < p.second;) {
			// Find largest subtree in range for which 'leaf' is the leftmost leaf
			index_t node   = leaf + (1u << nb_levels);
			index_t parent = node >> 1;
			while (parent > 1
				&& leftmostLeaf(parent) >= leaf
				&& rightmostLeaf(parent) <= p.second)
			{
				node = parent;
				parent >>= 1;
			}
			assert(rightmostLeaf(node) > leaf);
			leaf = rightmostLeaf(node);

			if ((node & mask) == node) {
				// Process internal node of the tree
				_impl()->m_Nodes[node].getPointsInSphere(
					points, lower, upper, center, sq_dist, neighbors);
			} else {
				// Process leaf node of the tree
				node = node & mask;
				const double *x = points + MaxDim*_impl()->m_Leaves[node];
				if (distLessThan(center, x, sq_dist)) {
					neighbors.emplace_back(_impl()->m_Leaves[node]);
				}
			}
		}
	}
};

////////////////////////////////////////////////////////////////////////////////
// Specializations for 1d, 2d and 3d range trees
////////////////////////////////////////////////////////////////////////////////

template<int MaxDim>
class RangeTree1d
	: public RangeTreeLeaves<0, MaxDim>
	, public RangeTreeInternal
{
	// Allow base class to access derived member variables and methods
	friend RangeTreeLeaves<0, MaxDim>;

	// Test whether two given points lie within a given distance
	static bool distLessThan (const double *p, const double *q, double sq_dist) {
		double d = 0;
		for (int i = 0; i < MaxDim; ++i) {
			d += (p[i] - q[i]) * (p[i] - q[i]);
		}
		return d < sq_dist;
	}

public:
	// Default empty constructor
	RangeTree1d () = default;

	// Creates a 1d range-tree from a list of 3d points
	RangeTree1d (index_t nb_points, const double *points)
		: RangeTreeLeaves<0, MaxDim>(nb_points)
	{
		std::sort(this->m_Leaves.begin(), this->m_Leaves.end(),
			[points] (index_t i, index_t j) {
				return points[MaxDim*i+0] < points[MaxDim*j+0];
		});
	}

	// Creates a 1d range-tree from a list of 3d points and (sorted) indices
	RangeTree1d (index_t nb_points, const double *points,
		const index_t *sortedByX)
		: RangeTreeLeaves<0, MaxDim>(nb_points, sortedByX)
	{ }

	// Sort the indices according to a new set of coordinates
	void rebuildIndex (const double *points) override {
		std::sort(this->m_Leaves.begin(), this->m_Leaves.end(),
			[points] (index_t i, index_t j) {
				return points[MaxDim*i+0] < points[MaxDim*j+0];
		});
	}

	// Computes the number of neighbors within a query box
	index_t countPointsInBox (
		const double *points,
		const double *lower,
		const double *upper) const override
	{
		auto p = this->getRangeBox(points, lower, upper);
		return p.second - p.first;
	}

	// Retrieve points within a query box (assumes output buffer has been allocated)
	index_t * getPointsInBox (
		const double *points,
		const double *lower,
		const double *upper,
		index_t *neighbors) const override
	{
		auto p = this->getRangeBox(points, lower, upper);
		std::copy(this->m_Leaves.data() + p.first,
			this->m_Leaves.data() + p.second, neighbors);
		return neighbors + p.second - p.first;
	}

	// Retrieve points within a query box (std::vector version)
	void getPointsInBox (
		const double *points,
		const double *lower,
		const double *upper,
		std::vector<index_t> &neighbors) const override
	{
		auto p = this->getRangeBox(points, lower, upper);
		neighbors.insert(neighbors.end(), this->m_Leaves.data() + p.first,
			this->m_Leaves.data() + p.second);
	}

	// Computes the number of neighbors within a query box
	index_t countPointsInSphere(
		const double *points,
		const double *lower,
		const double *upper,
		const double *center,
		double sq_dist) const override
	{
		auto p = this->getRangeBox(points, lower, upper);
		index_t accu = 0;
		for (index_t i = p.first; i < p.second; ++i) {
			const double *x = points + MaxDim*this->m_Leaves[i];
			if (distLessThan(center, x, sq_dist)) {
				++accu;
			}
		}
		return accu;
	}

	// Retrieve points within a query box (assumes output buffer has been allocated)
	index_t * getPointsInSphere (
		const double *points,
		const double *lower,
		const double *upper,
		const double *center,
		double sq_dist,
		index_t *neighbors) const override
	{
		auto p = this->getRangeBox(points, lower, upper);
		for (index_t i = p.first; i < p.second; ++i) {
			const double *x = points + MaxDim*this->m_Leaves[i];
			if (distLessThan(center, x, sq_dist)) {
				neighbors[0] = this->m_Leaves[i];
				++neighbors;
			}
		}
		return neighbors;
	}

	// Retrieve points within a query box (std::vector version)
	void getPointsInSphere (
		const double *points,
		const double *lower,
		const double *upper,
		const double *center,
		double sq_dist,
		std::vector<index_t> &neighbors) const override
	{
		auto p = this->getRangeBox(points, lower, upper);
		for (index_t i = p.first; i < p.second; ++i) {
			const double *x = points + MaxDim*this->m_Leaves[i];
			if (distLessThan(center, x, sq_dist)) {
				neighbors.emplace_back(this->m_Leaves[i]);
			}
		}
	}
};

////////////////////////////////////////////////////////////////////////////////
// 2d range-tree specialization
////////////////////////////////////////////////////////////////////////////////

template<int MaxDim>
class RangeTree2d
	: public RangeTreeLeaves<1, MaxDim>
	, public RangeTreeNodes<1, MaxDim, RangeTree2d<MaxDim> >
{
	// Allow base class to access derived member variables and methods
	friend class RangeTreeNodes<1, MaxDim, RangeTree2d<MaxDim> >;
	friend class RangeTreeLeaves<0, MaxDim>;
	friend class RangeTree3d<MaxDim>;

private:
	// Array of all 1d subtrees
	std::vector<RangeTree1d<MaxDim> > m_Nodes;

public:
	// Default empty constructor
	RangeTree2d () = default;

	// Creates a 2d range-tree from a list of 2d points
	RangeTree2d (index_t nb_points, const double *points)
		: RangeTreeLeaves<1, MaxDim>(nb_points)
		, m_Nodes(1 << nbits(nb_points))
	{
		// Alloc buffers
		std::vector<index_t> sortedByX(this->m_Leaves);
		std::vector<index_t> tempBufferX(this->m_Leaves.size());
		std::vector<unsigned char> predicate(this->m_Leaves.size());

		// Sort indices by X and Y coordinates
		std::sort(sortedByX.begin(), sortedByX.end(),
			[points] (index_t i, index_t j) {
				return points[MaxDim*i+0] < points[MaxDim*j+0];
		});
		std::sort(this->m_Leaves.begin(), this->m_Leaves.end(),
			[points] (index_t i, index_t j) {
				return points[MaxDim*i+1] < points[MaxDim*j+1];
		});

		// Build subtrees
		buildSubtrees<true> (nb_points, points, sortedByX, tempBufferX, predicate);
	}

	// Creates a 2d range-tree from a list of 3d points and indices
	RangeTree2d (index_t nb_points, const double *points,
		const index_t *indicesByX,
		const index_t *indicesByY,
		std::vector<unsigned char> &predicate)
		: RangeTreeLeaves<1, MaxDim>(nb_points, indicesByY)
		, m_Nodes(1 << nbits(nb_points))
	{
		ptx_assert(nb_points != 0);

		// Temporary buffers that will be passed down to subtrees
		// TODO: Use subtrees and move semantics to avoid extra memory allocation
		std::vector<index_t> sortedByX(indicesByX, indicesByX + nb_points);
		std::vector<index_t> tempBufferX(nb_points);

		// Build subtrees
		buildSubtrees<false> (nb_points, points, sortedByX, tempBufferX, predicate);
	}

	// Create a 1d range-tree for each internal node
	template<bool UseThreads>
	void buildSubtrees(index_t nb_points, const double *points,
		std::vector<index_t> &sortedByX,
		std::vector<index_t> &tempBufferX,
		std::vector<unsigned char> &predicate)
	{
		index_t nb_levels = nbits(nb_points);
		index_t mask = (1 << nb_levels) - 1;
		for (index_t level = 0; level < nb_levels; ++level) {

			// Iterate over all internal node of the current level
			auto innerLoop = [&] (index_t i) {
			//for (index_t i = (1u << level); i < (1u << (level + 1)); ++i) {
				index_t idx_start = (i << (nb_levels - level)) & mask;
				index_t idx_end   = idx_start + (1 << (nb_levels - level));
				idx_end = std::min(index_t(this->m_Leaves.size()), idx_end);
				assert(idx_start == this->leftmostLeaf(i));
				assert(idx_end   == this->rightmostLeaf(i));

				if (idx_start < idx_end) {
					// Build 1d subtree
					m_Nodes[i] = RangeTree1d<MaxDim>(idx_end - idx_start, points,
						sortedByX.data() + idx_start);

					// Stable partition of children leaves
					index_t idx_middle = idx_start + (1u << (nb_levels - level - 1));
					idx_middle = std::min(idx_end, idx_middle);

					for (index_t leaf = idx_start; leaf < idx_end; ++leaf) {
						predicate[this->m_Leaves[leaf]] = (leaf < idx_middle);
					}

					index_t idx_left  = idx_start;
					index_t idx_right = idx_middle;
					for (index_t idx = idx_start; idx < idx_end; ++idx) {
						if (predicate[sortedByX[idx]]) {
							assert(idx_left < nb_points);
							tempBufferX[idx_left++] = sortedByX[idx];
						} else {
							assert(idx_right < nb_points);
							tempBufferX[idx_right++] = sortedByX[idx];
						}
					}
				}
			};

			if (UseThreads) {
				ThreadPool::ParallelFor(1u << level, 1u << (level + 1), innerLoop);
			} else {
				ThreadPool::SequentialFor(1u << level, 1u << (level + 1), innerLoop);
			}

			std::swap(sortedByX, tempBufferX);
		}
	}

	// Rebuild the search trees according to a new set of coordinates
	void rebuildIndex (const double *points) override {
		// Sort indices by Y coordinates
		std::sort(this->m_Leaves.begin(), this->m_Leaves.end(),
			[points] (index_t i, index_t j) {
				return points[MaxDim*i+1] < points[MaxDim*j+1];
		});

		if (m_Nodes.size() > 1) {
			// Sort by X coordinates
			m_Nodes[1].m_Leaves = this->m_Leaves;
			std::sort(m_Nodes[1].m_Leaves.begin(), m_Nodes[1].m_Leaves.end(),
				[points] (index_t i, index_t j) {
					return points[MaxDim*i+0] < points[MaxDim*j+0];
			});
		}

		// Used to mark leaves when splitting a node into its 2 subtrees
		std::vector<unsigned char> predicate(this->m_Leaves.size());

		// Build subtrees
		propagateSubtrees<true> (predicate);
	}

	// Create a 1d range-tree for each internal node
	template<bool UseThreads>
	void propagateSubtrees (std::vector<unsigned char> &predicate) {
		index_t nb_levels = nbits(index_t(this->m_Leaves.size()));
		index_t mask = (1 << nb_levels) - 1;
		for (index_t level = 0; level < nb_levels; ++level) {

			// Iterate over all internal node of the current level
			auto innerLoop = [&] (index_t i) {
			//for (index_t i = (1u << level); i < (1u << (level + 1)); ++i) {
				index_t idx_start = (i << (nb_levels - level)) & mask;
				index_t idx_end   = idx_start + (1 << (nb_levels - level));
				idx_end = std::min(index_t(this->m_Leaves.size()), idx_end);
				assert(idx_start == this->leftmostLeaf(i));
				assert(idx_end   == this->rightmostLeaf(i));

				if (idx_start < idx_end) {
					// 1d subtree is supposed to be build for node i
					// m_Node[i].buildSubtrees ();

					// Stable partition of children leaves
					if (level + 1 < nb_levels) {
						index_t idx_middle = idx_start + (1u << (nb_levels - level - 1));
						idx_middle = std::min(idx_end, idx_middle);

						for (index_t leaf = idx_start; leaf < idx_end; ++leaf) {
							predicate[this->m_Leaves[leaf]] = (leaf < idx_middle);
						}

						index_t idx_left  = 0;
						index_t idx_right = 0;
						for (index_t idx = 0; idx < idx_end - idx_start; ++idx) {
							if (predicate[m_Nodes[i].m_Leaves[idx]]) {
								assert(idx_left < m_Nodes[2 * i].m_Leaves.size());
								m_Nodes[2 * i].m_Leaves[idx_left++] = m_Nodes[i].m_Leaves[idx];
							} else {
								assert(idx_right < m_Nodes[2 * i + 1].m_Leaves.size());
								m_Nodes[2 * i + 1].m_Leaves[idx_right++] = m_Nodes[i].m_Leaves[idx];
							}
						}
					}
				}
			};

			if (UseThreads) {
				ThreadPool::ParallelFor(1u << level, 1u << (level + 1), innerLoop);
			} else {
				ThreadPool::SequentialFor(1u << level, 1u << (level + 1), innerLoop);
			}
		}
	}

};

////////////////////////////////////////////////////////////////////////////////
// 3d range-tree specialization
////////////////////////////////////////////////////////////////////////////////

template<int MaxDim>
class RangeTree3d
	: public RangeTreeLeaves<2, 3>
	, public RangeTreeNodes<2, 3, RangeTree3d<MaxDim> >
{
	// Allow base class to access derived member variables and methods
	friend class RangeTreeNodes<2, 3, RangeTree3d<MaxDim> >;

private:
	// Array of all 2d subtrees
	std::vector<RangeTree2d<MaxDim> > m_Nodes;

public:
	// Default empty constructor
	RangeTree3d () = default;

	// Creates a 3d range-tree from a list of 3d points and indices
	RangeTree3d (index_t nb_points, const double *points)
		: RangeTreeLeaves(nb_points)
		, m_Nodes(1u << nbits(nb_points))
	{
		ptx_assert(nb_points != 0);

		// Temporary buffers that will be passed down to subtrees
		std::vector<index_t> sortedByX(m_Leaves);
		std::vector<index_t> sortedByY(m_Leaves);
		std::vector<index_t> tempBufferX(sortedByX.size());
		std::vector<index_t> tempBufferY(sortedByX.size());

		// Used to mark leaves when splitting a node into its 2 subtrees
		std::vector<unsigned char> predicate(m_Leaves.size(), false);

		// Sort all the points only once!
		std::sort(sortedByX.begin(), sortedByX.end(),
			[points] (index_t i, index_t j) {
				return points[MaxDim*i+0] < points[MaxDim*j+0];
		});
		std::sort(sortedByY.begin(), sortedByY.end(),
			[points] (index_t i, index_t j) {
				return points[MaxDim*i+1] < points[MaxDim*j+1];
		});
		std::sort(m_Leaves.begin(), m_Leaves.end(),
			[points] (index_t i, index_t j) {
				return points[MaxDim*i+2] < points[MaxDim*j+2];
		});

		// Build subtrees
		buildSubtrees (nb_points, points, sortedByX, sortedByY,
			tempBufferX, tempBufferY, predicate);
	}

	// Create a 1d range-tree for each internal node
	void buildSubtrees(index_t nb_points, const double *points,
		std::vector<index_t> &sortedByX,
		std::vector<index_t> &sortedByY,
		std::vector<index_t> &tempBufferX,
		std::vector<index_t> &tempBufferY,
		std::vector<unsigned char> &predicate)
	{
		index_t nb_levels = nbits(nb_points);
		index_t mask = (1u << nb_levels) - 1;
		for (index_t level = 0; level < nb_levels; ++level) {

			// Iterate over all internal node of the current level
			auto innerLoop = [&] (index_t i) {
			//ThreadPool::ParallelFor(1u << level, 1u << (level + 1), [&] (index_t i) {
			//for (index_t i = (1u << level); i < (1u << (level + 1)); ++i) {
				index_t idx_start = (i << (nb_levels - level)) & mask;
				index_t idx_end   = idx_start + (1u << (nb_levels - level));
				idx_end = std::min(index_t(m_Leaves.size()), idx_end);
				assert(idx_start == this->leftmostLeaf(i));
				assert(idx_end   == this->rightmostLeaf(i));

				if (idx_start < idx_end) {
					// Build 2d subtree
					m_Nodes[i] = RangeTree2d<MaxDim>(idx_end - idx_start, points,
						sortedByX.data() + idx_start,
						sortedByY.data() + idx_start,
						predicate);

					// Stable partition of children leaves
					index_t idx_middle = idx_start + (1u << (nb_levels - level - 1));
					idx_middle = std::min(idx_end, idx_middle);
					for (index_t leaf = idx_start; leaf < idx_end; ++leaf) {
						predicate[m_Leaves[leaf]] = (leaf < idx_middle);
					}

					index_t idx_left  = idx_start;
					index_t idx_right = idx_middle;
					for (index_t idx = idx_start; idx < idx_end; ++idx) {
						if (predicate[sortedByX[idx]]) {
							assert(idx_left < tempBufferX.size());
							tempBufferX[idx_left++] = sortedByX[idx];
						} else {
							assert(idx_right < tempBufferX.size());
							tempBufferX[idx_right++] = sortedByX[idx];
						}
					}

					idx_left  = idx_start;
					idx_right = idx_middle;
					for (index_t idx = idx_start; idx < idx_end; ++idx) {
						if (predicate[sortedByY[idx]]) {
							assert(idx_left < tempBufferY.size());
							tempBufferY[idx_left++] = sortedByY[idx];
						} else {
							assert(idx_left < tempBufferY.size());
							tempBufferY[idx_right++] = sortedByY[idx];
						}
					}

				}
			};

			ThreadPool::ParallelFor(1u << level, 1u << (level + 1), innerLoop);

			std::swap(tempBufferX, sortedByX);
			std::swap(tempBufferY, sortedByY);
		}
	}

	// Rebuild the search trees according to a new set of coordinates
	void rebuildIndex (const double *points) override {
		// Sort indices by Z coordinate
		std::sort(m_Leaves.begin(), m_Leaves.end(),
			[points] (index_t i, index_t j) {
				return points[MaxDim*i+2] < points[MaxDim*j+2];
		});

		if (m_Nodes.size() > 1) {
			// Sort indices by Y coordinate
			m_Nodes[1].m_Leaves = this->m_Leaves;
			std::sort(m_Nodes[1].m_Leaves.begin(), m_Nodes[1].m_Leaves.end(),
				[points] (index_t i, index_t j) {
					return points[MaxDim*i+1] < points[MaxDim*j+1];
			});

			// Sort indices by X coordinate
			m_Nodes[1].m_Nodes[1].m_Leaves = this->m_Leaves;
			std::sort(m_Nodes[1].m_Nodes[1].m_Leaves.begin(),
				m_Nodes[1].m_Nodes[1].m_Leaves.end(),
				[points] (index_t i, index_t j) {
					return points[MaxDim*i+0] < points[MaxDim*j+0];
			});
		}

		// Used to mark leaves when splitting a node into its 2 subtrees
		std::vector<unsigned char> predicate(m_Leaves.size(), false);

		// Build subtrees
		propagateSubtrees (predicate);
	}

	// Create a 2d range-tree for each internal node
	void propagateSubtrees (std::vector<unsigned char> &predicate) {
		index_t nb_levels = nbits(index_t(this->m_Leaves.size()));
		index_t mask = (1u << nb_levels) - 1;
		for (index_t level = 0; level < nb_levels; ++level) {

			// Iterate over all internal node of the current level
			auto innerLoop = [&] (index_t i) {
			//ThreadPool::ParallelFor(1u << level, 1u << (level + 1), [&] (index_t i) {
			//for (index_t i = (1u << level); i < (1u << (level + 1)); ++i) {
				index_t idx_start = (i << (nb_levels - level)) & mask;
				index_t idx_end   = idx_start + (1u << (nb_levels - level));
				idx_end = std::min(index_t(m_Leaves.size()), idx_end);
				assert(idx_start == this->leftmostLeaf(i));
				assert(idx_end   == this->rightmostLeaf(i));

				if (idx_start < idx_end) {
					// Build 2d subtree
					m_Nodes[i].template propagateSubtrees<false> (predicate);

					if (level + 1 < nb_levels) {
						// Stable partition of children leaves
						index_t idx_middle = idx_start + (1u << (nb_levels - level - 1));
						idx_middle = std::min(idx_end, idx_middle);
						for (index_t leaf = idx_start; leaf < idx_end; ++leaf) {
							predicate[m_Leaves[leaf]] = (leaf < idx_middle);
						}

						index_t idx_left  = 0;
						index_t idx_right = 0;
						for (index_t idx = 0; idx < idx_end - idx_start; ++idx) {
							if (predicate[m_Nodes[i].m_Leaves[idx]]) {
								assert(idx_left < m_Nodes[2 * i].m_Leaves.size());
								m_Nodes[2 * i].m_Leaves[idx_left++] =
									m_Nodes[i].m_Leaves[idx];
							} else {
								assert(idx_right < m_Nodes[2 * i + 1].m_Leaves.size());
								m_Nodes[2 * i + 1].m_Leaves[idx_right++] =
									m_Nodes[i].m_Leaves[idx];
							}
						}

						idx_left  = 0;
						idx_right = 0;
						for (index_t idx = 0; idx < idx_end - idx_start; ++idx) {
							if (predicate[m_Nodes[i].m_Nodes[1].m_Leaves[idx]]) {
								assert(idx_left < m_Nodes[2 * i].m_Nodes[1].m_Leaves.size());
								m_Nodes[2 * i].m_Nodes[1].m_Leaves[idx_left++] =
									m_Nodes[i].m_Nodes[1].m_Leaves[idx];
							} else {
								assert(idx_right < m_Nodes[2 * i + 1].m_Nodes[1].m_Leaves.size());
								m_Nodes[2 * i + 1].m_Nodes[1].m_Leaves[idx_right++] =
									m_Nodes[i].m_Nodes[1].m_Leaves[idx];
							}
						}
					}
				}
			};

			ThreadPool::ParallelFor(1u << level, 1u << (level + 1), innerLoop);
		}
	}
};

////////////////////////////////////////////////////////////////////////////////
// RangeTree member functions implementation
////////////////////////////////////////////////////////////////////////////////

// Constructor
RangeTree::RangeTree (unsigned dim, index_t nb_points, const double *points)
	: m_Dimension(dim)
	, m_NumberOfPoints(nb_points)
	, m_Points(points)
{
	switch (m_Dimension) {
	case 1:
		m_Tree = std::make_shared<RangeTree1d<1> >(nb_points, points);
		break;
	case 2:
		m_Tree = std::make_shared<RangeTree2d<2> >(nb_points, points);
		break;
	case 3:
		m_Tree = std::make_shared<RangeTree3d<3> >(nb_points, points);
		break;
	default:
		throw std::runtime_error("[RangeTree] Invalid Dimension");
	}
}

// Rebuild the search index (if the nb of points is unchanged no reallocation occurs)
void RangeTree::rebuild_index (index_t nb_points, const double *points) {
	if (nb_points != 0 && points != nullptr) {
		// Update pointer to points coordinates
		m_Points = points;
	}
	if (nb_points == 0 || nb_points == m_NumberOfPoints) {
		// Number of points in the set is unchanged, no reallocation
		m_Tree->rebuildIndex(m_Points);
	} else {
		// Number of points in the set has changed, rebuild from scratch
		switch (m_Dimension) {
		case 1:
			m_Tree = std::make_shared<RangeTree1d<1> >(m_NumberOfPoints, m_Points);
			break;
		case 2:
			m_Tree = std::make_shared<RangeTree2d<2> >(m_NumberOfPoints, m_Points);
			break;
		case 3:
			m_Tree = std::make_shared<RangeTree3d<3> >(m_NumberOfPoints, m_Points);
			break;
		default:
			throw std::runtime_error("[RangeTree] Invalid Dimension");
		}
	}
}

// Count points in box, arbitrary box shape
index_t RangeTree::nb_points_in_box (
	const double *query_point, const double *box_dist) const
{
	double lower[3];
	double upper[3];
	for (unsigned i = 0; i < m_Dimension; ++i) {
		lower[i] = query_point[i] - box_dist[i];
		upper[i] = query_point[i] + box_dist[i];
	}
	return m_Tree->countPointsInBox(m_Points, lower, upper);
}

// Count points in box, n-cube box shape
index_t RangeTree::nb_points_in_box (
	const double *query_point, double box_dist) const
{
	double lower[3];
	double upper[3];
	for (unsigned i = 0; i < m_Dimension; ++i) {
		lower[i] = query_point[i] - box_dist;
		upper[i] = query_point[i] + box_dist;
	}
	return m_Tree->countPointsInBox(m_Points, lower, upper);
}

// Retrieve points in box (assumes buffer is allocated), arbitrary box shape
index_t * RangeTree::get_points_in_box (const double *query_point,
	const double *box_dist, index_t * neighbors) const
{
	double lower[3];
	double upper[3];
	for (unsigned i = 0; i < m_Dimension; ++i) {
		lower[i] = query_point[i] - box_dist[i];
		upper[i] = query_point[i] + box_dist[i];
	}
	return m_Tree->getPointsInBox(m_Points, lower, upper, neighbors);
}

// Retrieve points in box (assumes buffer is allocated), n-cube box shape
index_t * RangeTree::get_points_in_box (
	const double *query_point, double box_dist, index_t * neighbors) const
{
	double lower[3];
	double upper[3];
	for (unsigned i = 0; i < m_Dimension; ++i) {
		lower[i] = query_point[i] - box_dist;
		upper[i] = query_point[i] + box_dist;
	}
	return m_Tree->getPointsInBox(m_Points, lower, upper, neighbors);
}

// Retrieve points in box (std::vector version), arbitrary box shape
void RangeTree::get_points_in_box (
	const double *query_point, const double *box_dist,
	std::vector<index_t> & neighbors) const
{
	double lower[3];
	double upper[3];
	for (unsigned i = 0; i < m_Dimension; ++i) {
		lower[i] = query_point[i] - box_dist[i];
		upper[i] = query_point[i] + box_dist[i];
	}
	m_Tree->getPointsInBox(m_Points, lower, upper, neighbors);
}

// Retrieve points in box (std::vector version), n-cube box version
void RangeTree::get_points_in_box (
	const double *query_point, double box_dist,
	std::vector<index_t> & neighbors) const
{
	double lower[3];
	double upper[3];
	for (unsigned i = 0; i < m_Dimension; ++i) {
		lower[i] = query_point[i] - box_dist;
		upper[i] = query_point[i] + box_dist;
	}
	m_Tree->getPointsInBox(m_Points, lower, upper, neighbors);
}

// Count points in sphere
index_t RangeTree::nb_points_in_sphere (
	const double *query_point, double l2_dist) const
{
	double lower[3];
	double upper[3];
	for (unsigned i = 0; i < m_Dimension; ++i) {
		lower[i] = query_point[i] - l2_dist;
		upper[i] = query_point[i] + l2_dist;
	}
	return m_Tree->countPointsInSphere(
		m_Points, lower, upper, query_point, l2_dist*l2_dist);
}

// Retrieve points in sphere (assumes buffer is allocated)
index_t * RangeTree::get_points_in_sphere (
	const double *query_point, double l2_dist, index_t * neighbors) const
{
	double lower[3];
	double upper[3];
	for (unsigned i = 0; i < m_Dimension; ++i) {
		lower[i] = query_point[i] - l2_dist;
		upper[i] = query_point[i] + l2_dist;
	}
	return m_Tree->getPointsInSphere(
		m_Points, lower, upper, query_point, l2_dist*l2_dist, neighbors);
}

// Retrieve points in sphere (std::vector version)
void RangeTree::get_points_in_sphere (
	const double *query_point, double l2_dist,
	std::vector<index_t> & neighbors) const
{
	double lower[3];
	double upper[3];
	for (unsigned i = 0; i < m_Dimension; ++i) {
		lower[i] = query_point[i] - l2_dist;
		upper[i] = query_point[i] + l2_dist;
	}
	m_Tree->getPointsInSphere(
		m_Points, lower, upper, query_point, l2_dist*l2_dist, neighbors);
}

/*
 * TODO:
 * - Generalize to handle non-contiguous data in point coordinates? (e.g. Eigen::Matrices)
 * - If needed, compare with GPU implementation?
 * - write method get_nearest_point (using binary search on box side,
 *   with the first guess given by nearest point in Z), and compare it
 *   with KdTree using L_\infty norm
 */
