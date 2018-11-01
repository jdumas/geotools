#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>
// -----------------------------------------------------------------------------
#include <vector>
#include <memory>
////////////////////////////////////////////////////////////////////////////////

// Forward declaration of abstract class implementing the internal interface
class RangeTreeInternal;

////////////////////////////////////////////////////////////////////////////////

class RangeTree {

public:
	// Public index type
	typedef unsigned int index_t;

private:
	// Dimension of the dataset (2d or 3d points)
	unsigned m_Dimension;

	// Number of points in the current set
	index_t m_NumberOfPoints;

	// Pointer to the points coordinates
	const double * m_Points;

	// Internal implementation
	std::shared_ptr<RangeTreeInternal> m_Tree;

public:
	// Empty default constructor
	RangeTree () = default;

	// Constructor
	RangeTree (unsigned dim, index_t nb_points, const double *points);

	////////////////////////
	// Index tree methods //
	////////////////////////

	// Return the dimension of the current dataset
	unsigned dimension () const { return m_Dimension; }

	// Rebuild the search index (if the nb of points is unchanged no reallocation occurs)
	void rebuild_index (index_t nb_points = 0, const double *points = nullptr);

	///////////////////////
	// Box query methods //
	///////////////////////

	// Count points in box, arbitrary box shape
	index_t nb_points_in_box (
		const double *query_point, const double *box_dist
	) const;

	// Count points in box, n-cube box shape
	index_t nb_points_in_box (const double *query_point, double box_dist) const;

	// Retrieve points in box (assumes buffer is allocated), arbitrary box shape
	index_t * get_points_in_box (
		const double *query_point, const double *box_dist, index_t * neighbors
	) const;

	// Retrieve points in box (assumes buffer is allocated), n-cube box shape
	index_t * get_points_in_box (
		const double *query_point, double box_dist, index_t * neighbors
	) const;

	// Retrieve points in box (std::vector version), arbitrary box shape
	void get_points_in_box (
		const double *query_point, const double *box_dist,
		std::vector<index_t> & neighbors
	) const;

	// Retrieve points in box (std::vector version), n-cube box version
	void get_points_in_box (
		const double *query_point, double box_dist,
		std::vector<index_t> & neighbors
	) const;

	//////////////////////////
	// Sphere query methods //
	//////////////////////////

	// Count points in sphere
	index_t nb_points_in_sphere (const double *query_point, double l2_dist) const;

	// Retrieve points in sphere (assumes buffer is allocated)
	index_t * get_points_in_sphere (
		const double *query_point, double l2_dist, index_t * neighbors
	) const;

	// Retrieve points in sphere (std::vector version)
	void get_points_in_sphere (
		const double *query_point, double l2_dist,
		std::vector<index_t> & neighbors
	) const;
};
