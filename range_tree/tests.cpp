////////////////////////////////////////////////////////////////////////////////
#include "TestRangeQueries.h"
#include "RangeTree.h"
#include "Chrono.h"
#include "ThreadPool.h"
#include "Common.h"
// -----------------------------------------------------------------------------
#include <Eigen/Dense>
#include <nanoflann.hpp>
#ifdef USE_GEOGRAM
#include <geogram/points/kd_tree.h>
#endif
#include <random>
////////////////////////////////////////////////////////////////////////////////

namespace Test {
	enum : int {
		NO_FLAG        = 0,
		TEST_GEOGRAM   = 1,
		TEST_NANOFLANN = 2,
		TEST_NAIVE     = 4,
		TEST_OPENMP    = 8,
	};

	void rangeTreeBox    (int n, int m, double dist, int dim = 3, int flags = NO_FLAG);
	void rangeTreeSphere (int n, int m, double dist, int dim = 3, int flags = NO_FLAG);
	[[noreturn]] void kdTree (int n, int m, double dist, int dim = 3, int flags = NO_FLAG);
}

// TODO: Test also performances of box queries in 2D (compare with nanoflann),
// using the fact that in 2D, 2*||(x,y)||_\infty = ||(x+y,x-y)||_1

////////////////////////////////////////////////////////////////////////////////

namespace {

class NaiveRangeSearch {
public:
	// Public index type
	typedef RangeTree::index_t index_t;

private:
	// Dimension of the dataset
	const index_t m_Dimension;

	// Number of points
	const index_t  m_NbPoints;

	// Pointer to the points coordinates
	const double * m_Points;

public:
	// Creates a naive search wrapper from a list of 3d points and indices
	NaiveRangeSearch (
		index_t dimension,
		index_t nb_points,
		const double *points)
		: m_Dimension(dimension)
		, m_NbPoints(nb_points)
		, m_Points(points)
	{ }

	// Box test, arbitrary box shape
	bool in_box(const double *p, const double *c, const double *d) const {
		bool ok = true;
		for (index_t i = 0; i < m_Dimension; ++i) {
			ok = ok && (p[i] >= c[i] - d[i] && p[i] <= c[i] + d[i]);
		}
		return ok;
	}

	// Box test, n-cube box shape
	bool in_box(const double *p, const double *c, double d) const {
		bool ok = true;
		for (index_t i = 0; i < m_Dimension; ++i) {
			ok = ok && (p[i] >= c[i] - d && p[i] <= c[i] + d);
		}
		return ok;
	}

	// Shpere test
	bool in_sphere(const double *p, const double *c, double d) const {
		double l = 0;
		for (index_t i = 0; i < m_Dimension; ++i) {
			l += (p[i] - c[i]) * (p[i] - c[i]);
		}
		return l <= d*d;
	}

	// Computes the number of neighbors within a query box, arbitrary box shape
	index_t nb_points_in_box (const double *query, const double *dist) const {
		index_t accu = 0;
		for (index_t i = 0; i < m_NbPoints; ++i) {
			if (in_box(m_Points + m_Dimension*i, query, dist)) {
				++accu;
			}
		}
		return accu;
	}

	// Computes the number of neighbors within a query box, n-cube box shape
	index_t nb_points_in_box (const double *query, double dist) const {
		index_t accu = 0;
		for (index_t i = 0; i < m_NbPoints; ++i) {
			if (in_box(m_Points + m_Dimension*i, query, dist)) {
				++accu;
			}
		}
		return accu;
	}

	// Retrieve points in box (assumes buffer is allocated), arbitrary box shape
	index_t * get_points_in_box (const double *query, const double *dist,
		index_t * neighbors) const
	{
		for (index_t i = 0; i < m_NbPoints; ++i) {
			if (in_box(m_Points + m_Dimension*i, query, dist)) {
				neighbors[0] = i;
				++neighbors;
			}
		}
		return neighbors;
	}

	// Retrieve points in box (assumes buffer is allocated), n-cube box shape
	index_t * get_points_in_box (const double *query, double dist,
		index_t * neighbors) const
	{
		for (index_t i = 0; i < m_NbPoints; ++i) {
			if (in_box(m_Points + m_Dimension*i, query, dist)) {
				neighbors[0] = i;
				++neighbors;
			}
		}
		return neighbors;
	}

	// Computes the number of neighbors within a query sphere
	index_t nb_points_in_sphere (const double *query, double dist) const {
		index_t accu = 0;
		for (index_t i = 0; i < m_NbPoints; ++i) {
			if (in_sphere(m_Points + m_Dimension*i, query, dist)) {
				++accu;
			}
		}
		return accu;
	}

	// Retrieve points in sphere (assumes buffer is allocated)
	index_t * get_points_in_sphere (const double *query, double dist,
		index_t * neighbors) const {
		for (index_t i = 0; i < m_NbPoints; ++i) {
			if (in_sphere(m_Points + m_Dimension*i, query, dist)) {
				neighbors[0] = i;
				++neighbors;
			}
		}
		return neighbors;
	}

};

}

////////////////////////////////////////////////////////////////////////////////
// Box queries for range-trees
////////////////////////////////////////////////////////////////////////////////

void Test::rangeTreeBox (int n, int m, double box_dist, int dim, int flags) {
	typedef RangeTree::index_t index_t;
	Chrono tm("RangeTree");

	// Generate seeds and query points
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0, 100);
	std::vector<double> pts(size_t(dim*n));
	std::vector<double> queries(size_t(dim*m));
	for (index_t i = 0; i < index_t(dim*n); ++i) {
		pts[i] = distribution(generator);
	}
	for (index_t i = 0; i < index_t(dim*m); ++i) {
		queries[i] = distribution(generator);
	}

	// Build range tree (twice to compare with the noalloc version)
	tm.tic("Building");
	RangeTree rangeTree((unsigned char) dim, (index_t) n, pts.data());
	tm.toc(false);

	tm.tic("Rebuilding");
	rangeTree.rebuild_index((index_t) n, pts.data());
	tm.toc(false);

	// Compare number of neighbors found with the naive search
	if (flags & TEST_NAIVE) {
		NaiveRangeSearch naiveTree((index_t) dim, (index_t) n, pts.data());
		tm.tic("Queries");
		index_t counter = 0;
		for (int i = 0; i < m; ++i) {
			const double *p = queries.data() + dim*i;
			index_t matches = rangeTree.nb_points_in_box(p, box_dist);
			index_t naive = naiveTree.nb_points_in_box(p, box_dist);
			ptx_assert(matches == naive);
			counter += matches;
		}
		std::cout << "Naive + box queries: " << counter
			<< " matches." << std::endl;
	}

	// Perform parallel box queries
	std::vector<std::vector<index_t> > allNeighs((size_t) m);
	tm.tic("Queries");
	if (flags & TEST_OPENMP) {
		// Use openmp parallel for
		#pragma omp parallel for
		for (int i = 0; i < m; ++i) {
			const double *p = queries.data() + index_t(dim)*i;
			rangeTree.get_points_in_box(p, box_dist, allNeighs[i]);
		}
	} else {
		// Use custom thread pool
		ThreadPool::ParallelFor(0u, index_t(m), [&] (index_t i) {
			const double *p = queries.data() + index_t(dim)*i;
			rangeTree.get_points_in_box(p, box_dist, allNeighs[i]);
		});
	}
	tm.toc(false);

	// Count number of neighbors
	index_t nb_neighs = 0;
	for (const auto &t : allNeighs) { nb_neighs += (index_t) t.size(); }
	std::cout << "Number of neighbors: " << nb_neighs << std::endl;

	// Second query round (without allocation)
	for (index_t i = 0; i < index_t(m); ++i) {
		allNeighs[i].clear();
	}
	tm.tic("Queries");
	if (flags & TEST_OPENMP) {
		// Use openmp parallel for
		#pragma omp parallel for
		for (int i = 0; i < m; ++i) {
			const double *p = queries.data() + index_t(dim)*i;
			rangeTree.get_points_in_box(p, box_dist, allNeighs[i]);
		}
	} else {
		// Use custom thread pool
		ThreadPool::ParallelFor(0u, index_t(m), [&] (index_t i) {
			const double *p = queries.data() + index_t(dim)*i;
			rangeTree.get_points_in_box(p, box_dist, allNeighs[i]);
		});
	}
	tm.toc(false);

	// Compare with KdTree from geogram
	#ifdef USE_GEOGRAM
	if (flags & TEST_GEOGRAM) {
		Chrono gm("GeoTree");
		GEO::NearestNeighborSearch_var nnsearch =
			GEO::NearestNeighborSearch::create((GEO::coord_index_t) dim, "BNN");

		gm.tic("Building");
		nnsearch->set_points((index_t) n, pts.data());
		gm.toc(false);

		std::vector<GEO::index_t> nearest(nb_neighs);
		std::vector<double> sq_dist(nb_neighs);
		size_t current = 0;

		gm.tic("Queries");
		for (int i = 0; i < m; ++i) {
			const double *p = queries.data() + dim*i;
			nnsearch->get_nearest_neighbors((index_t) allNeighs[i].size(), p,
				nearest.data() + current,
				sq_dist.data() + current);
			current += allNeighs[i].size();
		}
		gm.toc(false);
	}
	#endif

	// TODO: Compare with KdTree from nanoflann using a custom adaptator
	// (transform L2 distance into L_\infty distance?)
	if (flags & TEST_NANOFLANN) {
		throw std::runtime_error("Not implemented yet.");
	}
}

////////////////////////////////////////////////////////////////////////////////
// Sphere queries for range-trees
////////////////////////////////////////////////////////////////////////////////

#if 0
template <typename DistanceType, typename IndexType = size_t>
class RadiusResultCount
{
public:
	const DistanceType radius;

	size_t m_indices_count;
	std::pair<IndexType,DistanceType> m_worst_index;

	inline RadiusResultCount(DistanceType radius_) : radius(radius_), m_indices_count(0)
	{
		init();
	}

	inline ~RadiusResultCount() { }

	inline void init() { clear(); }
	inline void clear() { m_indices_count = 0; }

	inline size_t size() const { return m_indices_count; }

	inline bool full() const { return true; }

	inline void addPoint(DistanceType dist, IndexType index)
	{
		if (dist<radius) {
			++m_indices_count;
			auto new_index = std::make_pair(index,dist);
			if (m_indices_count == 1 || m_worst_index < new_index) {
				m_worst_index = new_index;
			}
		}
	}

	inline DistanceType worstDist() const { return radius; }

	/** Clears the result set and adjusts the search radius. */
	inline void set_radius_and_clear( const DistanceType r )
	{
		radius = r;
		clear();
	}

	/**
	 * Find the worst result (furtherest neighbor) without copying or sorting
	 * Pre-conditions: size() > 0
	 */
	std::pair<IndexType,DistanceType> worst_item() const
	{
		if (m_indices_count == 0) throw std::runtime_error("Cannot invoke RadiusResultSet::worst_item() on an empty list of results.");
		return m_worst_index;
	}
};
#endif

void Test::rangeTreeSphere (int n, int m, double l2_dist, int dim, int flags) {
	typedef RangeTree::index_t index_t;
	Chrono tm("RangeTree");

	// Generate seeds and query points
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0, 100);
	auto randomGenerator = [&] (int a) { return distribution(generator); };
	Eigen::MatrixXd seeds = Eigen::MatrixXd::NullaryExpr(n, dim, randomGenerator);
	std::vector<double> points(size_t(dim*n));
	std::vector<double> queries(size_t(dim*m));
	// Copy generated coordinates in a row-major buffer
	for (int i = 0; i < n; ++i) {
		for (int k = 0; k < dim; ++k) {
			points[index_t(dim*i+k)] = seeds(i, k);
		}
	}
	for (index_t i = 0; i < index_t(dim*m); ++i) {
		queries[i] = distribution(generator);
	}

	// Build range tree (twice to compare with the noalloc version)
	tm.tic("Building");
	RangeTree rangeTree((unsigned char) dim, (index_t) n, points.data());
	tm.toc(false);

	tm.tic("Rebuilding");
	rangeTree.rebuild_index((index_t) n, points.data());
	tm.toc(false);

	// Compare number of neighbors found with the naive search
	if (flags & TEST_NAIVE) {
		NaiveRangeSearch naiveTree((index_t) dim, (index_t) n, points.data());
		tm.tic("Queries");
		index_t counter = 0;
		for (int i = 0; i < m; ++i) {
			const double *p = queries.data() + index_t(dim)*i;
			index_t matches = rangeTree.nb_points_in_sphere(p, l2_dist);
			index_t naive = naiveTree.nb_points_in_sphere(p, l2_dist);
			ptx_assert(matches == naive);
			counter += matches;
		}
		std::cout << "Naive + sphere queries: " << counter
			<< " matches." << std::endl;
	}

	// Perform parallel box queries
	std::vector<std::vector<index_t> > allNeighs((size_t) m);
	tm.tic("Queries");
	if (flags & TEST_OPENMP) {
		// Use openmp parallel for
		#pragma omp parallel for
		for (int i = 0; i < m; ++i) {
			const double *p = queries.data() + index_t(dim)*i;
			rangeTree.get_points_in_sphere(p, l2_dist, allNeighs[i]);
		}
	} else {
		// Use custom thread pool
		ThreadPool::ParallelFor(0u, index_t(m), [&] (index_t i) {
			const double *p = queries.data() + index_t(dim)*i;
			rangeTree.get_points_in_sphere(p, l2_dist, allNeighs[i]);
		});
	}
	tm.toc(false);

	// Count number of neighbors
	index_t nb_neighs = 0;
	for (const auto &t : allNeighs) { nb_neighs+= (index_t) t.size(); }
	std::cout << "Number of neighbors: " << nb_neighs << std::endl;

	// Second query round (without allocation)
	for (index_t i = 0; i < index_t(m); ++i) {
		allNeighs[i].clear();
	}
	tm.tic("Queries");
	if (flags & TEST_OPENMP) {
		// Use openmp parallel for
		#pragma omp parallel for
		for (int i = 0; i < m; ++i) {
			const double *p = queries.data() + index_t(dim)*i;
			rangeTree.get_points_in_sphere(p, l2_dist, allNeighs[i]);
		}
	} else {
		// Use custom thread pool
		ThreadPool::ParallelFor(0u, index_t(m), [&] (index_t i) {
			const double *p = queries.data() + index_t(dim)*i;
			rangeTree.get_points_in_sphere(p, l2_dist, allNeighs[i]);
		});
	}
	tm.toc(false);

	// Compare with KdTree from geogram
	#ifdef USE_GEOGRAM
	if (flags & TEST_GEOGRAM) {
		Chrono gm("GeoTree");
		GEO::NearestNeighborSearch_var nnsearch =
			GEO::NearestNeighborSearch::create((GEO::coord_index_t) dim, "BNN");

		gm.tic("Building");
		nnsearch->set_points((index_t) n, points.data());
		gm.toc(false);

		std::vector<GEO::index_t> nearest(nb_neighs);
		std::vector<double> sq_dist(nb_neighs);
		std::vector<size_t> offset(m);

		// Exclusive scan to determine write index into global array
		size_t current = 0;
		for (int i = 0; i < m; ++i) {
			offset[i] = current;
			current += allNeighs[i].size();
		}

		gm.tic("Queries");
		ThreadPool::ParallelFor(0u, index_t(m), [&] (index_t i) {
			const double *p = queries.data() + index_t(dim)*i;
			nnsearch->get_nearest_neighbors((index_t) allNeighs[i].size(), p,
				nearest.data() + offset[i],
				sq_dist.data() + offset[i]);
		});
		gm.toc(false);
	}
	#endif

	if (flags & TEST_NANOFLANN) {
		Chrono fm("Nanoflann");

		// Put the seeds in a KdTree using nanoflann
		typedef nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixXd, 3,
			nanoflann::metric_L2> my_kd_tree_t;

		const int leaf_max_size = 10;

		fm.tic("Building");
		my_kd_tree_t mat_index(3, seeds, leaf_max_size);
		mat_index.index->buildIndex();
		fm.toc(false);

		// Search parameters and results
		nanoflann::SearchParams params;
		params.sorted = false;
		typedef my_kd_tree_t::IndexType IndexType;
		std::vector<std::vector<std::pair<IndexType, double> > > ret_matches(m);

		//std::atomic<size_t> nb_all(0);
		tm.tic("Queries");
		ThreadPool::ParallelFor(0u, index_t(m), [&] (index_t i) {
			const double *p = queries.data() + index_t(dim)*i;
			const size_t nb_matches = mat_index.index->radiusSearch(p,
				l2_dist * l2_dist, ret_matches[i], params);
			ptx_assert(nb_matches == allNeighs[i].size());

			//RadiusResultCount<double, my_kd_tree_t::IndexType> resultSet(l2_dist * l2_dist);
			//const size_t nFound = mat_index.index->radiusSearchCustomCallback(
			//	p, resultSet, params);
			//nb_all += nFound;
		});
		tm.toc(false);
		//std::cout << nb_all << std::endl;
	}
}

////////////////////////////////////////////////////////////////////////////////
// Compare KdTree performances
////////////////////////////////////////////////////////////////////////////////

[[noreturn]] void Test::kdTree (int n, int m, double dist, int dim, int flags) {
	// TODO: Write a cleaner test later
	throw std::runtime_error("Not implemented yet");

#if 0
	Chrono tm("KdTree");

	// Generate random distribution
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0, 100);
	auto randomGenerator = [&] (int a) { return distribution(generator); };
	tm.tic("Generate");
	Eigen::MatrixX3d v = Eigen::MatrixX3d::NullaryExpr(n, 3, randomGenerator);
	Eigen::MatrixX3d q = Eigen::MatrixX3d::NullaryExpr(m, 3, randomGenerator);
	tm.toc();

	// Put every point in a KdTree using nanoflann
	typedef nanoflann::KDTreeEigenMatrixAdaptor<
				Eigen::MatrixX3d,
				3,
				nanoflann::metric_L1> my_kd_tree_t;

	tm.tic("Building");
	my_kd_tree_t mat_index(3, v, 16 /* max leaf */ );
	mat_index.index->buildIndex();
	tm.toc();

	// Do a knn search
	const double search_radius = 10.0;
	std::vector<std::pair<long int, double> > ret_matches;

	nanoflann::SearchParams params;
	params.sorted = false;

	tm.tic("Queries");
	size_t totalCount = 0;
	std::vector<GEO::index_t> nnMatches((GEO::index_t) m);
	for (int i = 0; i < m; ++i) {
		double qq[3] = {q(i, 0), q(i, 1), q(i, 2)};
		const size_t nMatches = mat_index.index->radiusSearch(qq,
			search_radius, ret_matches, params);
		/*std::cout << "Query point = " << q.row(i) << std::endl;
		for (size_t k = 0; k < nMatches ; ++k) {
			std::cout << "idx[" << k << "]=" << ret_matches[k].first
				<< " dist[" << k << "]=" << ret_matches[k].second << std::endl;
		}*/
		nnMatches[size_t(i)] = GEO::index_t(nMatches);
		totalCount += nMatches;
	}
	tm.toc();

	std::cout << "radiusSearch(): radius=" << search_radius << " -> " << totalCount << " matches\n";
	std::cout << "\n";

	// Compare with KdTree from geogram
	if (test_geogram) {
		Chrono gm("GeoTree");
		GEO::NearestNeighborSearch_var nnsearch =
			GEO::NearestNeighborSearch::create(3, "BNN");

		gm.tic("Building");
		nnsearch->set_points(GEO::index_t(n), v.data());
		gm.toc();

		gm.tic("Alloc");
		std::vector<GEO::index_t> nearest(totalCount);
		std::vector<double> sq_dist(totalCount);
		gm.toc();
		size_t current = 0;
		gm.tic("Queries");
		for (GEO::index_t i = 0; i < GEO::index_t(m); ++i) {
			nnsearch->get_nearest_neighbors(nnMatches[i], &q(i, 0),
				nearest.data() + current,
				sq_dist.data() + current);
			current += nnMatches[i];
		}
		gm.toc();
	}
#endif
}
