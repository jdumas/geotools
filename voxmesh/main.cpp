/*
 *  Copyright (c) 2012-2014, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine,
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX
 *     FRANCE
 *
 */

// A mesh file converter using Geogram.


#include "octree.h"
#include <geogram/basic/file_system.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/progress.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/numerics/predicates.h>
#include <algorithm>
#include <array>
#include <iterator>

////////////////////////////////////////////////////////////////////////////////

bool endswith (std::string const &str, std::string const &ending) {
	if (str.length() >= ending.length()) {
		return (0 == str.compare(str.length() - ending.length(), ending.length(), ending));
	} else {
		return false;
	}
}

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

////////////////////////////////////////////////////////////////////////////////

namespace Layout {

	GEO::vec3i index3_from_index(int idx, GEO::vec3i size) {
		return GEO::vec3i(
			idx % size[0],
			(idx / size[0]) % size[1],
			(idx / size[0]) / size[1]
		);
	}

	int index_from_index3(GEO::vec3i vx, GEO::vec3i size) {
		return (vx[2] * size[1] + vx[1]) * size[0] + vx[0];
	}

}

namespace GEO {

	bool filename_has_supported_extension(const std::string &filename) {
		std::vector<std::string> extensions;
		GEO::MeshIOHandlerFactory::list_creators(extensions);
		for (auto &ext : extensions) {
			if (endswith(filename, ext)) {
				return true;
			}
		}
		return false;
	}
}

////////////////////////////////////////////////////////////////////////////////

template<typename T>
class VoxelGrid {
private:
	// Member data
	std::vector<T> m_data;
	GEO::vec3      m_origin;
	double         m_spacing; // voxel size (in mm)
	GEO::vec3i     m_grid_size;

public:
	// Interface
	VoxelGrid(GEO::vec3 origin, GEO::vec3 extent, double voxel_size, int padding);

	GEO::vec3i grid_size() const { return m_grid_size; }
	int num_voxels() const { return m_grid_size[0] * m_grid_size[1]  * m_grid_size[2]; }

	GEO::vec3 origin() const { return m_origin; }
	double spacing() const { return m_spacing; }

	GEO::vec3i index3_from_index(int idx) const { return Layout::index3_from_index(idx, m_grid_size); }
	int index_from_index3(GEO::vec3i vx) const { return Layout::index_from_index3(vx, m_grid_size); }

	GEO::vec3 voxel_center(int x, int y, int z) const;

	const T at(int idx) const { return m_data[idx]; }
	T & at(int idx) { return m_data[idx]; }
	const T * rawbuf() const { return m_data.data(); }
	T * raw_layer(int z) { return m_data.data() + z * m_grid_size[1] * m_grid_size[0]; }
};

// -----------------------------------------------------------------------------

template<typename T>
VoxelGrid<T>::VoxelGrid(GEO::vec3 origin, GEO::vec3 extent, double spacing, int padding)
	: m_origin(origin)
	, m_spacing(spacing)
{
	m_origin -= padding * spacing * GEO::vec3(1, 1, 1);
	m_grid_size[0] = (int) std::ceil(extent[0] / spacing) + 2 * padding;
	m_grid_size[1] = (int) std::ceil(extent[1] / spacing) + 2 * padding;
	m_grid_size[2] = (int) std::ceil(extent[2] / spacing) + 2 * padding;
	GEO::Logger::out("Voxels") << "Grid size: "
		<< m_grid_size[0] << " x " << m_grid_size[1] << " x " << m_grid_size[2] << std::endl;
	m_data.assign(m_grid_size[0] * m_grid_size[1] * m_grid_size[2], T(0));
}

template<typename T>
GEO::vec3 VoxelGrid<T>::voxel_center(int x, int y, int z) const {
	GEO::vec3 pos;
	pos[0] = (x + 0.5) * m_spacing;
	pos[1] = (y + 0.5) * m_spacing;
	pos[2] = (z + 0.5) * m_spacing;
	return pos + m_origin;
}

////////////////////////////////////////////////////////////////////////////////

template<typename T>
class DexelGrid {
private:
	// Member data
	std::vector<std::vector<T> > m_data;
	GEO::vec3      m_origin;
	double         m_spacing; // voxel size (in mm)
	GEO::vec2i     m_grid_size;

public:
	// Interface
	DexelGrid(GEO::vec3 origin, GEO::vec3 extent, double voxel_size, int padding);

	GEO::vec2i grid_size() const { return m_grid_size; }
	int num_dexels() const { return m_grid_size[0] * m_grid_size[1]; }

	GEO::vec3 origin() const { return m_origin; }
	double spacing() const { return m_spacing; }

	GEO::vec2 dexel_center(int x, int y) const;

	const std::vector<T> & at(int x, int y) const { return m_data[x + m_grid_size[0]*y]; }
	std::vector<T> & at(int x, int y) { return m_data[x + m_grid_size[0]*y]; }
};

// -----------------------------------------------------------------------------

template<typename T>
DexelGrid<T>::DexelGrid(GEO::vec3 origin, GEO::vec3 extent, double spacing, int padding)
	: m_origin(origin)
	, m_spacing(spacing)
{
	m_origin -= padding * spacing * GEO::vec3(1, 1, 1);
	m_grid_size[0] = (int) std::ceil(extent[0] / spacing) + 2 * padding;
	m_grid_size[1] = (int) std::ceil(extent[1] / spacing) + 2 * padding;
	GEO::Logger::out("Voxels") << "Grid size: "
		<< m_grid_size[0] << " x " << m_grid_size[1] << std::endl;
	m_data.assign(m_grid_size[0] * m_grid_size[1], std::vector<T>(0));
}

template<typename T>
GEO::vec2 DexelGrid<T>::dexel_center(int x, int y) const {
	GEO::vec2 pos;
	pos[0] = (x + 0.5) * m_spacing;
	pos[1] = (y + 0.5) * m_spacing;
	return pos + GEO::vec2(m_origin[0], m_origin[1]);
}

////////////////////////////////////////////////////////////////////////////////
// NOTE: Function `point_in_triangle_2d` comes from SDFGen by Christopher Batty.
// https://github.com/christopherbatty/SDFGen/blob/master/makelevelset3.cpp
////////////////////////////////////////////////////////////////////////////////

// calculate twice signed area of triangle (0,0)-(x1,y1)-(x2,y2)
// return an SOS-determined sign (-1, +1, or 0 only if it's a truly degenerate triangle)
int orientation(
	double x1, double y1, double x2, double y2, double &twice_signed_area)
{
	twice_signed_area=y1*x2-x1*y2;
	if(twice_signed_area>0) return 1;
	else if(twice_signed_area<0) return -1;
	else if(y2>y1) return 1;
	else if(y2<y1) return -1;
	else if(x1>x2) return 1;
	else if(x1<x2) return -1;
	else return 0; // only true when x1==x2 and y1==y2
}

// -----------------------------------------------------------------------------

// robust test of (x0,y0) in the triangle (x1,y1)-(x2,y2)-(x3,y3)
// if true is returned, the barycentric coordinates are set in a,b,c.
bool point_in_triangle_2d(
	double x0, double y0, double x1, double y1,
	double x2, double y2, double x3, double y3,
	double &a, double &b, double &c)
{
	x1-=x0; x2-=x0; x3-=x0;
	y1-=y0; y2-=y0; y3-=y0;
	int signa=orientation(x2, y2, x3, y3, a);
	if(signa==0) return false;
	int signb=orientation(x3, y3, x1, y1, b);
	if(signb!=signa) return false;
	int signc=orientation(x1, y1, x2, y2, c);
	if(signc!=signa) return false;
	double sum=a+b+c;
	geo_assert(sum!=0); // if the SOS signs match and are nonzero, there's no way all of a, b, and c are zero.
	a/=sum;
	b/=sum;
	c/=sum;
	return true;
}

// -----------------------------------------------------------------------------

// \brief Computes the (approximate) orientation predicate in 2d.
// \details Computes the sign of the (approximate) signed volume of
//  the triangle p0, p1, p2
// \param[in] p0 first vertex of the triangle
// \param[in] p1 second vertex of the triangle
// \param[in] p2 third vertex of the triangle
// \retval POSITIVE if the triangle is oriented positively
// \retval ZERO if the triangle is flat
// \retval NEGATIVE if the triangle is oriented negatively
// \todo check whether orientation is inverted as compared to
//   Shewchuk's version.
inline GEO::Sign orient_2d_inexact(GEO::vec2 p0, GEO::vec2 p1, GEO::vec2 p2) {
	double a11 = p1[0] - p0[0] ;
	double a12 = p1[1] - p0[1] ;

	double a21 = p2[0] - p0[0] ;
	double a22 = p2[1] - p0[1] ;

	double Delta = GEO::det2x2(
		a11, a12,
		a21, a22
	);

	return GEO::geo_sgn(Delta);
}

////////////////////////////////////////////////////////////////////////////////


/**
 * @brief      { Intersect a vertical ray with a triangle }
 *
 * @param[in]  M     { Mesh containing the triangle to intersect }
 * @param[in]  f     { Index of the facet to intersect }
 * @param[in]  q     { Query point (only XY coordinates are used) }
 * @param[out] z     { Intersection }
 *
 * @return     { {-1,0,1} depending on the sign of the intersection. }
 */
template<int X = 0, int Y = 1, int Z = 2>
int intersect_ray_z(const GEO::Mesh &M, GEO::index_t f, const GEO::vec3 &q, double &z) {
	using namespace GEO;

	index_t c = M.facets.corners_begin(f);
	const vec3& p1 = Geom::mesh_vertex(M, M.facet_corners.vertex(c++));
	const vec3& p2 = Geom::mesh_vertex(M, M.facet_corners.vertex(c++));
	const vec3& p3 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));

	double u, v, w;
	if (point_in_triangle_2d(
		q[X], q[Y], p1[X], p1[Y], p2[X], p2[Y], p3[X], p3[Y], u, v, w))
	{
		z = u*p1[Z] + v*p2[Z] + w*p3[Z];
		auto sign = orient_2d_inexact(vec2(p1[X], p1[Y]), vec2(p2[X], p2[Y]), vec2(p3[X], p3[Y]));
		switch (sign) {
		case GEO::POSITIVE: return 1;
		case GEO::NEGATIVE: return -1;
		default: return 0;
		}
	}

	return 0;
}

// -----------------------------------------------------------------------------

template<typename T>
void compute_sign(const GEO::Mesh &M,
	const GEO::MeshFacetsAABB &aabb_tree, VoxelGrid<T> &voxels)
{
	const GEO::vec3i size = voxels.grid_size();

	try {
		GEO::ProgressTask task("Ray marching", 100);

		GEO::vec3 min_corner, max_corner;
		GEO::get_bbox(M, &min_corner[0], &max_corner[0]);

		const GEO::vec3 origin = voxels.origin();
		const double spacing = voxels.spacing();

		GEO::parallel_for([&](int y) {
			if (GEO::Thread::current()->id() == 0) {
				task.progress((int) (100.0 * y / size[1] * GEO::Process::number_of_cores()));
			}
			for (int x = 0; x < size[0]; ++x) {
				GEO::vec3 center = voxels.voxel_center(x, y, 0);

				GEO::Box box;
				box.xyz_min[0] = box.xyz_max[0] = center[0];
				box.xyz_min[1] = box.xyz_max[1] = center[1];
				box.xyz_min[2] = min_corner[2] - spacing;
				box.xyz_max[2] = max_corner[2] + spacing;

				std::vector<std::pair<double, int>> inter;
				auto action = [&M, &inter, &center] (GEO::index_t f) {
					double z;
					if (int s = intersect_ray_z(M, f, center, z)) {
						inter.emplace_back(z, s);
					}
				};
				aabb_tree.compute_bbox_facet_bbox_intersections(box, action);
				std::sort(inter.begin(), inter.end());

				for (int z = 0, s = 0, i = 0; z < size[2]; ++z) {
					GEO::vec3 center = voxels.voxel_center(x, y, z);
					for (; i < inter.size() && inter[i].first < center[2]; ++i) {
						s += inter[i].second;
					}
					const int idx = voxels.index_from_index3(GEO::vec3i(x, y, z));
					voxels.at(idx) = T(s < 0 ? 1 : 0);
				}
			}
		}, 0, size[1]);
	} catch(const GEO::TaskCanceled&) {
		// Do early cleanup
	}
}

// -----------------------------------------------------------------------------

template<typename T>
void compute_sign(const GEO::Mesh &M,
	const GEO::MeshFacetsAABB &aabb_tree, DexelGrid<T> &dexels)
{
	const GEO::vec2i size = dexels.grid_size();

	try {
		GEO::ProgressTask task("Ray marching", 100);

		GEO::vec3 min_corner, max_corner;
		GEO::get_bbox(M, &min_corner[0], &max_corner[0]);

		const GEO::vec3 origin = dexels.origin();
		const double spacing = dexels.spacing();

		GEO::parallel_for([&](int y) {
			if (GEO::Thread::current()->id() == 0) {
				// task.progress((int) (100.0 * y / size[1] * GEO::Process::number_of_cores()));
			}
			for (int x = 0; x < size[0]; ++x) {
				GEO::vec2 center_xy = dexels.dexel_center(x, y);
				GEO::vec3 center(center_xy[0], center_xy[1], 0);

				GEO::Box box;
				box.xyz_min[0] = box.xyz_max[0] = center[0];
				box.xyz_min[1] = box.xyz_max[1] = center[1];
				box.xyz_min[2] = min_corner[2] - spacing;
				box.xyz_max[2] = max_corner[2] + spacing;

				std::vector<std::pair<double, int>> inter;
				auto action = [&M, &inter, &center] (GEO::index_t f) {
					double z;
					if (int s = intersect_ray_z(M, f, center, z)) {
						inter.emplace_back(z, s);
					}
				};
				aabb_tree.compute_bbox_facet_bbox_intersections(box, action);
				std::sort(inter.begin(), inter.end());

				std::vector<double> reduced;
				for (int i = 0, s = 0; i < inter.size(); ++i) {
					const int ds = inter[i].second;
					s += ds;
					if ((s == -1 && ds < 0) || (s == 0 && ds > 0)) {
						reduced.push_back(inter[i].first);
					}
				}

				dexels.at(x, y).resize(reduced.size());
				std::copy_n(reduced.begin(), reduced.size(), dexels.at(x, y).begin());
			}
		}, 0, size[1]);
	} catch(const GEO::TaskCanceled&) {
		// Do early cleanup
	}
}

// -----------------------------------------------------------------------------

void compute_sign(const GEO::Mesh &M, const GEO::MeshFacetsAABB &aabb_tree,
	OctreeGrid &octree, GEO::vec3 origin, double spacing)
{
	Eigen::VectorXf & inside = octree.cellAttributes.create<float>("inside");
	inside.setZero();

	try {
		GEO::ProgressTask task("Ray marching", 100);

		GEO::vec3 min_corner, max_corner;
		GEO::get_bbox(M, &min_corner[0], &max_corner[0]);

		GEO::parallel_for([&](int cellId) {
			auto cell_xyz_min = octree.cellCornerPos(cellId, OctreeGrid::CORNER_X0_Y0_Z0);
			auto extent = octree.cellExtent(cellId);

			GEO::Box box;
			box.xyz_min[0] = origin[0] + spacing * cell_xyz_min[0];
			box.xyz_min[1] = origin[1] + spacing * cell_xyz_min[1];
			box.xyz_max[0] = box.xyz_min[0] + spacing * extent;
			box.xyz_max[1] = box.xyz_min[1] + spacing * extent;
			box.xyz_min[2] = min_corner[2] - spacing;
			box.xyz_max[2] = max_corner[2] + spacing;

			GEO::vec3 center(
				box.xyz_min[0] + 0.5 * spacing * extent,
				box.xyz_min[1] + 0.5 * spacing * extent,
				origin[2] + spacing * cell_xyz_min[2] + 0.5 * spacing * extent
			);

			std::vector<double> inter;
			auto action = [&M, &inter, &center] (GEO::index_t f) {
				double z;
				if (intersect_ray_z(M, f, center, z)) {
					inter.push_back(z);
				}
			};
			aabb_tree.compute_bbox_facet_bbox_intersections(box, action);
			//std::sort(inter.begin(), inter.end());

			int num_before = 0;
			for (double z : inter) {
				if (z < center[2]) { ++num_before; }
			}
			if (num_before % 2 == 1) {
				inside(cellId) = 1.0;
			}
		}, 0, octree.numCells());
	} catch(const GEO::TaskCanceled&) {
		// Do early cleanup
	}
}

////////////////////////////////////////////////////////////////////////////////

typedef unsigned char num_t;

// -----------------------------------------------------------------------------

void paraview_dump(std::string &filename, const VoxelGrid<num_t> &voxels) {
	GEO::vec3i size = voxels.grid_size();

	std::string extension = GEO::FileSystem::extension(filename);
	std::string basename = GEO::FileSystem::base_name(filename, true);
	if (!extension.empty() && extension[0] != '.') {
		extension = "." + extension;
	}
	std::string outname = filename.substr(0, filename.size() - extension.size());
	std::ofstream metafile(outname + ".mhd");
	metafile << "ObjectType = Image\nNDims = 3\n"
		<< "DimSize = " << size[0] << " " << size[1] << " " << size[2] << "\n"
		<< "ElementType = MET_CHAR\nElementDataFile = " + basename + ".raw\n";
	metafile.close();

	std::ofstream rawfile(outname + ".raw", std::ios::binary);
	rawfile.write(reinterpret_cast<const char*>(voxels.rawbuf ()), voxels.num_voxels() * sizeof(num_t));
	rawfile.close();
}

// -----------------------------------------------------------------------------

void triangle_mesh_dump(std::string &filename, const VoxelGrid<num_t> &voxels) {
	using GEO::vec3i;

	vec3i cell_size = voxels.grid_size();
	vec3i node_size = cell_size + vec3i(1, 1, 1);
	int num_cells = voxels.num_voxels();
	int num_nodes = node_size[0] * node_size[1] * node_size[2];

	// Create triangle list from voxel grid
	GEO::vector<GEO::index_t> triangles;
	for (int idx = 0; idx < num_cells; ++idx) {
		vec3i pos = voxels.index3_from_index(idx);

		// Skip empty voxels
		if (voxels.at(idx) == num_t(0)) { continue; }

		// Define corner index
		std::array<GEO::index_t, 8> corners;
		corners[0] = Layout::index_from_index3(pos + vec3i(0, 0, 0), node_size);
		corners[1] = Layout::index_from_index3(pos + vec3i(1, 0, 0), node_size);
		corners[2] = Layout::index_from_index3(pos + vec3i(1, 1, 0), node_size);
		corners[3] = Layout::index_from_index3(pos + vec3i(0, 1, 0), node_size);
		corners[4] = Layout::index_from_index3(pos + vec3i(0, 0, 1), node_size);
		corners[5] = Layout::index_from_index3(pos + vec3i(1, 0, 1), node_size);
		corners[6] = Layout::index_from_index3(pos + vec3i(1, 1, 1), node_size);
		corners[7] = Layout::index_from_index3(pos + vec3i(0, 1, 1), node_size);

		// Subroutine to emit a facet quad
		auto check_facet = [&](int axis, int delta, int v1, int v2, int v3, int v4) {
			// Compute neigh voxel position
			vec3i neigh = pos;
			neigh[axis] += delta;

			// Check whether neigh voxel is empty
			bool neigh_is_empty = false;
			if (neigh[axis] < 0 || neigh[axis] >= cell_size[axis]) {
				neigh_is_empty = true;
			} else {
				int neigh_idx = voxels.index_from_index3(neigh);
				neigh_is_empty = (voxels.at(neigh_idx) == num_t(0));
			}

			// If neigh voxel is empty, emit triangle strips
			if (neigh_is_empty) {
				triangles.insert(triangles.end(), { corners[v1], corners[v2], corners[v3] });
				triangles.insert(triangles.end(), { corners[v3], corners[v2], corners[v4] });
			}
		};

		// Check adjacency and emit facets
		check_facet(0, -1, 0, 4, 3, 7); // left facet
		check_facet(0,  1, 2, 6, 1, 5); // right facet
		check_facet(1, -1, 1, 5, 0, 4); // front facet
		check_facet(1,  1, 3, 7, 2, 6); // back facet
		check_facet(2, -1, 1, 0, 2, 3); // lower facet
		check_facet(2,  1, 4, 5, 7, 6); // upper facet
	}

	// Assign vertex id (and remap triangle list)
	int num_vertices = 0;
	std::vector<int> node_id(num_nodes, -1);
	for (GEO::index_t &c : triangles) {
		if (node_id[c] == -1) {
			node_id[c] = num_vertices++;
		}
		c = node_id[c];
	}

	// Create Geogram mesh
	GEO::Mesh M;
	M.vertices.create_vertices(num_vertices);
	for (int v = 0; v < num_nodes; ++v) {
		if (node_id[v] != -1) {
			vec3i pos = Layout::index3_from_index(v, node_size);
			M.vertices.point(node_id[v]) = GEO::vec3(pos);
		}
	}
	M.facets.assign_triangle_mesh(triangles, true);

	// Connect facets
	M.facets.connect();

	// Rescale to unit box, and set min corner to 0
	// TODO: Add option to normalize, or output original real-world positions
	GEO::vec3 min_corner, max_corner;
	GEO::get_bbox(M, &min_corner[0], &max_corner[0]);
	GEO::vec3 extent = max_corner - min_corner;
	double scaling = std::max(extent[0], std::max(extent[1], extent[2]));
	for (int v = 0; v < M.vertices.nb(); ++v) {
		M.vertices.point(v) = (M.vertices.point(v) - min_corner) / scaling;
	}

	// Save mesh
	GEO::mesh_save(M, filename);
}

// -----------------------------------------------------------------------------

void volume_mesh_dump(std::string &filename, const VoxelGrid<num_t> &voxels) {
	using GEO::vec3i;

	vec3i cell_dims = voxels.grid_size();
	vec3i node_dims = cell_dims + vec3i(1, 1, 1);
	int num_nodes = node_dims[0] * node_dims[1] * node_dims[2];

	auto delta = [](int i) {
		return vec3i((i & 1) ^ ((i >> 1) & 1), (i >> 1) & 1, (i >> 2) & 1);
	};
	auto inv_delta = [](vec3i u) {
		return (u[1] ? 4*u[2] + 3 - u[0] : 4*u[2] + u[0]);
	};
	auto cell_corner_id = [&](int cell_id, int corner_id) {
		auto pos = Layout::index3_from_index(cell_id, cell_dims);
		pos += delta(corner_id);
		return Layout::index_from_index3(pos, node_dims);
	};

	GEO::Mesh mesh;
	mesh.vertices.create_vertices(num_nodes);
	for (int idx = 0; idx < num_nodes; ++idx) {
		vec3i posi = Layout::index3_from_index(idx, node_dims);
		GEO::vec3 pos(posi[0], posi[1], posi[2]);
		mesh.vertices.point(idx) = voxels.origin() + pos * voxels.spacing();
	}

	for (int cell_id = 0; cell_id < voxels.num_voxels(); ++cell_id) {
		if (voxels.at(cell_id) >= 0.5) {
			vec3i diff[8] = {
				vec3i(0,0,0), vec3i(1,0,0), vec3i(0,1,0), vec3i(1,1,0),
				vec3i(0,0,1), vec3i(1,0,1), vec3i(0,1,1), vec3i(1,1,1)
			};
			int v[8];
			for (GEO::index_t lv = 0; lv < 8; ++lv) {
				int corner_id = inv_delta(diff[lv]);
				v[lv] = cell_corner_id(cell_id, corner_id);
			}
			mesh.cells.create_hex(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]);
		}
	}

	mesh.cells.compute_borders();
	mesh.cells.connect();
	mesh.vertices.remove_isolated();

	GEO::mesh_save(mesh, filename);
}

// -----------------------------------------------------------------------------

template<typename T>
void dexel_dump(std::string &filename, const DexelGrid<T> &dexels) {
	using GEO::vec3;

	GEO::Mesh mesh;

	for (int y = 0; y < dexels.grid_size()[1]; ++y) {
		for (int x = 0; x < dexels.grid_size()[0]; ++x) {
			for (int i = 0; 2 * i < dexels.at(x, y).size(); ++i) {
				GEO::vec3 xyz_min, xyz_max;
				xyz_min[0] = dexels.origin()[0] + x * dexels.spacing();
				xyz_min[1] = dexels.origin()[1] + y * dexels.spacing();
				xyz_min[2] = dexels.at(x, y)[2*i+0];
				xyz_max[0] = dexels.origin()[0] + (x+1) * dexels.spacing();
				xyz_max[1] = dexels.origin()[1] + (y+1) * dexels.spacing();
				xyz_max[2] = dexels.at(x, y)[2*i+1];
				vec3 diff[8] = {
					vec3(0,0,0), vec3(1,0,0), vec3(0,1,0), vec3(1,1,0),
					vec3(0,0,1), vec3(1,0,1), vec3(0,1,1), vec3(1,1,1)
				};
				int v = mesh.vertices.nb();
				for (int lv = 0; lv < 8; ++lv) {
					for (int d = 0; d < 3; ++d) {
						diff[lv][d] = xyz_min[d] + diff[lv][d] * (xyz_max[d] - xyz_min[d]);
					}
					diff[lv] += dexels.origin();
					mesh.vertices.create_vertex(&diff[lv][0]);
				}
				mesh.cells.create_hex(v, v+1, v+2, v+3, v+4, v+5, v+6, v+7);
			}
		}
	}

	mesh.cells.compute_borders();
	mesh.cells.connect();
	mesh.vertices.remove_isolated();

	GEO::mesh_save(mesh, filename);
}

////////////////////////////////////////////////////////////////////////////////

void compute_octree(const GEO::Mesh &M, const GEO::MeshFacetsAABB &aabb_tree,
	const std::string &filename, GEO::vec3 min_corner, GEO::vec3 extent,
	double spacing, int padding, bool graded)
{
	GEO::vec3 origin =  min_corner - padding * spacing * GEO::vec3(1, 1, 1);
	Eigen::Vector3i grid_size(
		std::ceil(extent[0] / spacing) + 2 * padding,
		std::ceil(extent[1] / spacing) + 2 * padding,
		std::ceil(extent[2] / spacing) + 2 * padding
	);

	OctreeGrid octree(grid_size);

	// Subdivide cells
	auto should_subdivide = [&](int x, int y, int z, int extent) {
		if (extent == 1) { return false; }
		GEO::Box box;
		box.xyz_min[0] = origin[0] + spacing * x;
		box.xyz_min[1] = origin[1] + spacing * y;
		box.xyz_min[2] = origin[2] + spacing * z;
		box.xyz_max[0] = box.xyz_min[0] + spacing * extent;
		box.xyz_max[1] = box.xyz_min[1] + spacing * extent;
		box.xyz_max[2] = box.xyz_min[2] + spacing * extent;
		bool has_triangles = false;
		auto action = [&has_triangles](int id) { has_triangles = true; };
		aabb_tree.compute_bbox_facet_bbox_intersections(box, action);
		return has_triangles;
	};
	octree.subdivide(should_subdivide, graded);

	// Compute inside/outside info
	compute_sign(M, aabb_tree, octree, origin, spacing);

	// Export
	GEO::Logger::div("Saving");
	GEO::Mesh M_out;
	{
		Eigen::Vector3d o(origin[0], origin[1], origin[2]);
		Eigen::Vector3d s(spacing, spacing, spacing);
		GEO::Logger::out("OctreeGrid") << "Creating volume mesh..." << std::endl;
		octree.createMesh(M_out, o, s);
	}

	GEO::mesh_save(M_out, filename);
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
	// Initialize the Geogram library
	GEO::initialize();

	// Import standard command line arguments, and custom ones
	GEO::CmdLine::import_arg_group("standard");
	GEO::CmdLine::declare_arg("padding", 0, "Number of padded grid cells");
	GEO::CmdLine::declare_arg("resolution", 1.0, "Size of a voxel (in mm)");
	GEO::CmdLine::declare_arg("numvoxels", -1, "Number of voxels along the longest axis");
	GEO::CmdLine::declare_arg("hexmesh", false, "Output a hex mesh of the voxelized data");
	GEO::CmdLine::declare_arg("dexelize", false, "Dexelize the input model");
	GEO::CmdLine::declare_arg("octree", false, "Generate an adaptive octree of the input model");
	GEO::CmdLine::declare_arg("graded", false, "Should the octree be 2:1 graded");

	// Parse command line options and filenames
	std::vector<std::string> filenames;
	if (!GEO::CmdLine::parse(argc, argv, filenames, "input_mesh <output_file>")) {
		return 1;
	}

	int padding = std::max(0, GEO::CmdLine::get_arg_int("padding"));
	double voxel_size = GEO::CmdLine::get_arg_double("resolution");
	int num_voxels = GEO::CmdLine::get_arg_int("numvoxels");
	bool hexmesh = GEO::CmdLine::get_arg_bool("hexmesh");
	bool dexelize = GEO::CmdLine::get_arg_bool("dexelize");
	bool octree = GEO::CmdLine::get_arg_bool("octree");
	bool graded = GEO::CmdLine::get_arg_bool("graded");

	// Default output filename is "output" if unspecified
	if(filenames.size() == 1) {
		filenames.push_back("output");
	}

	// Display input and output filenames
	GEO::Logger::div("Command line");
	GEO::Logger::out("VoxMesh") << "Input file: " << filenames[0] << std::endl;
	GEO::Logger::out("VoxMesh") << "Output file: " << filenames[1] << std::endl;

	// Declare a mesh
	GEO::Mesh M;

	// Load the mesh and display timings
	GEO::Logger::div("Loading");
	{
		GEO::Stopwatch W("Load");
		if(!GEO::mesh_load(filenames[0], M)) {
			return 1;
		}
		geo_assert(M.vertices.dimension() == 3);
	}

	// Initialize voxel grid and AABB tree
	GEO::vec3 min_corner, max_corner;
	GEO::get_bbox(M, &min_corner[0], &max_corner[0]);
	GEO::vec3 extent = max_corner - min_corner;
	if (num_voxels > 0) {
		// Force number of voxels along longest axis
		double max_extent = std::max(extent[0], std::max(extent[1], extent[2]));
		voxel_size = max_extent / num_voxels;
	}
	GEO::MeshFacetsAABB aabb_tree(M);

	// Dexelize the input mesh
	if (dexelize) {
		GEO::Logger::div("Dexelizing");
		DexelGrid<double> dexels(min_corner, extent, voxel_size, padding);
		compute_sign<double>(M, aabb_tree, dexels);

		GEO::Logger::div("Saving");
		dexel_dump(filenames[1], dexels);
		return 0;
	}

	// Compute an octree of the input mesh
	if (octree) {
		GEO::Logger::div("Octree");
		compute_octree(M, aabb_tree, filenames[1], min_corner, extent, voxel_size, padding, graded);
		return 0;
	}

	// Compute inside/outside info
	GEO::Logger::div("Voxelizing");
	VoxelGrid<num_t> voxels(min_corner, extent, voxel_size, padding);
	compute_sign(M, aabb_tree, voxels);

	// Save voxel grid and display timings
	GEO::Logger::div("Saving");
	{
		GEO::Stopwatch W("Save");
		if (hexmesh) {
			volume_mesh_dump(filenames[1], voxels);
		} else if (GEO::filename_has_supported_extension(filenames[1])) {
			triangle_mesh_dump(filenames[1], voxels);
		} else {
			paraview_dump(filenames[1], voxels);
		}
	}

	return 0;
}

