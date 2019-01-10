////////////////////////////////////////////////////////////////////////////////
#include <geogram/basic/common.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/progress.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/process.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_topology.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_fill_holes.h>
#include <geogram/mesh/mesh_preprocessing.h>
#include <geogram/mesh/mesh_degree3_vertices.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/voronoi/RVD.h>
#include <geogram/numerics/predicates.h>
#include <boost/compute/core.hpp>
#include <boost/compute/container.hpp>
#include "mesh_AABB.h"
#include <omp.h>
#include <set>
#include <queue>
#include <limits>
#include <type_traits>
// -----------------------------------------------------------------------------
namespace compute = boost::compute;
////////////////////////////////////////////////////////////////////////////////

namespace GEO {

	/**
	 * \brief Set command line options of the program
	 */
	void set_programm_options() {
		CmdLine::import_arg_group("standard");
		CmdLine::import_arg_group("algo");
		CmdLine::declare_arg("padding", 2, "number of padded grid cells");
		CmdLine::declare_arg("resolution", 1.0, "size of a voxel (in mm)");
		CmdLine::declare_arg("voxelize_only", false, "only compute inside/outside info (no distance)");
		CmdLine::declare_arg("opencl_use_float", false, "kernels use float instead of double (opencl)");
		CmdLine::declare_arg("local_work_size", 8, "horizontal local work group size (opencl)");
		CmdLine::declare_arg("zslab", 5, "number of slices to be computed on the GPU (opencl)");
		CmdLine::declare_arg("use_gpu", true, "use gpu kernels to speedup computation (opencl)");
		CmdLine::declare_arg("check_result", false, "check resulting values against the CPU version (opencl)");
	}

	void get_point_facet_nearest_point(
		const Mesh& M,
		const vec3& p,
		index_t f,
		vec3& nearest_p,
		double& squared_dist
	) {
		geo_debug_assert(M.facets.nb_vertices(f) == 3);
		index_t c = M.facets.corners_begin(f);
		const vec3& p1 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));
		++c;
		const vec3& p2 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));
		++c;
		const vec3& p3 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));
		double lambda1, lambda2, lambda3;  // barycentric coords, not used.
		squared_dist = Geom::point_triangle_squared_distance(
			p, p1, p2, p3, nearest_p, lambda1, lambda2, lambda3
		);
	}

	vec3 get_facet_barycenter(const Mesh& M, index_t f) {
		geo_debug_assert(M.facets.nb_vertices(f) == 3);
		index_t c = M.facets.corners_begin(f);
		const vec3& p1 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));
		++c;
		const vec3& p2 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));
		++c;
		const vec3& p3 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));
		return (p1 + p2 + p3) / 3.0;
	}

}

////////////////////////////////////////////////////////////////////////////////

template<typename real>
class AABBTreeOpenCL {
	compute::context       &m_context;
	compute::command_queue &m_queue;
	compute::kernel         m_dist_kernel;

	compute::vector<real>     cl_vertices;
	compute::vector<unsigned> cl_facets;
	compute::vector<real>     cl_xyz_min;
	compute::vector<real>     cl_xyz_max;
	compute::vector<real>     cl_layer_sq_dist;
	std::vector<real>         m_layer_sq_dist;

	cl_uint m_nb_facets;
	cl_uint3 m_layer_size;

	compute::extents<3> m_local_work_size;
	compute::extents<3> m_global_work_size;

public:
	AABBTreeOpenCL(
		const GEO::Mesh &M, const GEO::MeshFacetsAABB &aabb_tree,
		compute::context &ctx, compute::command_queue &queue);

	void set_layer_size(int nx, int ny, int nz = 1);
	void compute_layer(float *raw_layer, GEO::vec3 origin, real spacing, int oz, int nz);
};

// -----------------------------------------------------------------------------

template<typename real>
AABBTreeOpenCL<real>::AABBTreeOpenCL(
	const GEO::Mesh &M, const GEO::MeshFacetsAABB &aabb_tree,
	compute::context &ctx, compute::command_queue &queue)
	: m_context(ctx)
	, m_queue(queue)
	, cl_vertices(3*M.vertices.nb(), ctx)
	, cl_facets(3*M.facets.nb(), ctx)
	, cl_xyz_min(3*aabb_tree.bboxes_.size(), ctx)
	, cl_xyz_max(3*aabb_tree.bboxes_.size(), ctx)
	, cl_layer_sq_dist(ctx)
	, m_nb_facets(M.facets.nb())
	, m_local_work_size(GEO::CmdLine::get_arg_int("local_work_size"))
{
	using namespace GEO;

	// Host vector for vertices, min_corners, max_corners
	std::vector<real>   vertices(3*M.vertices.nb());
	std::vector<unsigned> facets(3*M.facets.nb());
	std::vector<real> xyz_min(3*aabb_tree.bboxes_.size());
	std::vector<real> xyz_max(3*aabb_tree.bboxes_.size());

	// Fill host vectors
	for (index_t v = 0; v < M.vertices.nb(); ++v) {
		vec3 p = M.vertices.point(v);
		vertices[3*v+0] = (real) p[0];
		vertices[3*v+1] = (real) p[1];
		vertices[3*v+2] = (real) p[2];
	}
	for (index_t f = 0; f < M.facets.nb(); ++f) {
		facets[3*f+0] = M.facets.vertex(f, 0);
		facets[3*f+1] = M.facets.vertex(f, 1);
		facets[3*f+2] = M.facets.vertex(f, 2);
	}
	for (index_t b = 0; b < aabb_tree.bboxes_.size(); ++b) {
		for (int c = 0; c < 3; ++c) {
			xyz_min[3*b+c] = (real) aabb_tree.bboxes_[b].xyz_min[c];
			xyz_max[3*b+c] = (real) aabb_tree.bboxes_[b].xyz_max[c];
		}
	}

	// Copy data to the device
	compute::copy(vertices.begin(), vertices.end(), cl_vertices.begin(), queue);
	compute::copy(facets.begin(), facets.end(), cl_facets.begin(), queue);
	compute::copy(xyz_min.begin(), xyz_min.end(), cl_xyz_min.begin(), queue);
	compute::copy(xyz_max.begin(), xyz_max.end(), cl_xyz_max.begin(), queue);

	// Create program and kernel
	compute::program prog = compute::program::create_with_source_file(
		KERNEL_DIR "search.cl", ctx);

	try {
		// Attempt to compile to program
		if (GEO::CmdLine::get_arg_bool("opencl_use_float")) {
			prog.build("-I " KERNEL_DIR);
		} else {
			prog.build("-D USE_DOUBLE -I " KERNEL_DIR);
		}
	} catch(boost::compute::opencl_error &e){
		// Program failed to compile, print out the build log
		std::cout << prog.build_log() << std::endl;
	}
	m_dist_kernel = prog.create_kernel("squared_dist");

	m_local_work_size[2] = 1;
}

template<typename real>
void AABBTreeOpenCL<real>::set_layer_size(int nx, int ny, int nz) {
	cl_layer_sq_dist.resize(nx * ny * nz);
	m_layer_sq_dist.resize(nx * ny * nz);
	m_layer_size.s[0] = nx;
	m_layer_size.s[1] = ny;
	m_layer_size.s[2] = nz;
	for (int c = 0; c < 3; ++c) {
		m_global_work_size[c] = 0;
		while (m_global_work_size[c] < m_layer_size.s[c]) {
			m_global_work_size[c] += m_local_work_size[c];
		}
	}
	GEO::Logger::out("OpenCL") << "Local work size:  " << m_local_work_size[0]
		<< " x " << m_local_work_size[1] << " x " << m_local_work_size[2] << std::endl;
	GEO::Logger::out("OpenCL") << "Global work size: " << m_global_work_size[0]
		<< " x " << m_global_work_size[1] << " x " << m_global_work_size[2] << std::endl;
}

template<typename real> struct cl_helper;
template<> struct cl_helper<float>  { typedef cl_float3 real3; };
template<> struct cl_helper<double> { typedef cl_double3 real3; };

template<typename real>
void AABBTreeOpenCL<real>::compute_layer(
	float *raw_layer, GEO::vec3 origin, real spacing, int oz, int nz)
{
	typename cl_helper<real>::real3 orig;
	orig.s[0] = (real) origin[0];
	orig.s[1] = (real) origin[1];
	orig.s[2] = (real) origin[2];
	cl_uint3 offset;
	offset.s[0] = 0; offset.s[1] = 0; offset.s[2] = oz;
	// Set kernel arguments and execute
	m_dist_kernel.set_args(cl_vertices, cl_facets, cl_xyz_min, cl_xyz_max,
		cl_layer_sq_dist, m_layer_size, offset, orig, spacing, m_nb_facets);
	m_queue.enqueue_nd_range_kernel(m_dist_kernel, compute::extents<3>(0),
		m_global_work_size, m_local_work_size);
	// Copy data back to host memory
	compute::copy(cl_layer_sq_dist.begin(), cl_layer_sq_dist.end(), m_layer_sq_dist.begin(), m_queue);
	for (int z = 0; z < nz; ++z) {
		for (size_t y = 0; y < m_layer_size.s[1]; ++y) {
			for (size_t x = 0; x < m_layer_size.s[0]; ++x) {
				const size_t idx = x + m_layer_size.s[0] * (y + m_layer_size.s[1] * z);
				raw_layer[idx] = (float) std::sqrt(m_layer_sq_dist[idx]);
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////

typedef std::array<int, 3> Vec3i;

class VoxelGrid {
private:
	// Member data
	std::vector<float> m_data;
	GEO::vec3 m_origin;
	double m_spacing; // voxel size (in mm)
	Vec3i m_grid_size;

public:
	// Interface
	VoxelGrid(GEO::vec3 origin, GEO::vec3 extent, double voxel_size, int padding);

	Vec3i grid_size() const { return m_grid_size; }
	int num_voxels() const { return m_grid_size[0] * m_grid_size[1]  * m_grid_size[2]; }

	GEO::vec3 origin() const { return m_origin; }
	double spacing() const { return m_spacing; }

	Vec3i index3_from_index(int idx) const;
	int index_from_index3(Vec3i vx) const;

	GEO::vec3 voxel_center(int x, int y, int z) const;

	float & at(int idx) { return m_data[idx]; }
	const float * rawbuf() const { return m_data.data(); }
	float * raw_layer(int z) { return m_data.data() + z * m_grid_size[1] * m_grid_size[0]; }
};

////////////////////////////////////////////////////////////////////////////////

VoxelGrid::VoxelGrid(GEO::vec3 origin, GEO::vec3 extent, double spacing, int padding)
	: m_origin(origin)
	, m_spacing(spacing)
{
	m_origin -= padding * spacing * GEO::vec3(1, 1, 1);
	m_grid_size[0] = (int) std::floor(extent[0] / spacing) + 2 * padding;
	m_grid_size[1] = (int) std::floor(extent[1] / spacing) + 2 * padding;
	m_grid_size[2] = (int) std::floor(extent[2] / spacing) + 2 * padding;
	m_data.assign(m_grid_size[0] * m_grid_size[1] * m_grid_size[2], 1.0f);
	GEO::Logger::out("Voxels") << "Grid size: "
		<< m_grid_size[0] << " x " << m_grid_size[1] << " x " << m_grid_size[2] << std::endl;
}

GEO::vec3 VoxelGrid::voxel_center(int x, int y, int z) const {
	GEO::vec3 pos;
	pos[0] = (x + 0.5) * m_spacing;
	pos[1] = (y + 0.5) * m_spacing;
	pos[2] = (z + 0.5) * m_spacing;
	return pos + m_origin;
}

Vec3i VoxelGrid::index3_from_index(int idx) const {
	return {{
		idx % m_grid_size[0],
		(idx / m_grid_size[0]) % m_grid_size[1],
		(idx / m_grid_size[0]) / m_grid_size[1],
	}};
}

int VoxelGrid::index_from_index3(Vec3i vx) const {
	return (vx[2] * m_grid_size[1] + vx[1]) * m_grid_size[0] + vx[0];
}

////////////////////////////////////////////////////////////////////////////////

template<typename real>
void compute_unsigned_distance_field_gpu(const GEO::Mesh &M,
	const GEO::MeshFacetsAABB &aabb_tree, VoxelGrid &voxels)
{
	// Get the default compute device
	compute::device gpu = compute::system::default_device();
	GEO::Logger::out("OpenCL") << "Using device: " << gpu.name() << std::endl;

	// Create a compute context and command queue
	compute::context ctx(gpu);
	compute::command_queue queue(ctx, gpu);

	// Init aabb tree data on the GPU
	AABBTreeOpenCL<real> gpu_tree(M, aabb_tree, ctx, queue);

	const Vec3i size = voxels.grid_size();
	int nz = GEO::CmdLine::get_arg_int("zslab");
	gpu_tree.set_layer_size(size[0], size[1], nz);
	try {
		GEO::ProgressTask task("Sqdist (GPU)", 100);
		for (int z = 0; z < size[2]; z += nz) {
			task.progress((int) (100.0 * z / size[2]));
			int sz = std::min(nz, size[2] - z);
			gpu_tree.compute_layer(voxels.raw_layer(z), voxels.origin(),
				(real) voxels.spacing(), z, sz);
		}
	} catch(const GEO::TaskCanceled&) {
		// Do early cleanup
	}

	if (GEO::CmdLine::get_arg_bool("check_result")) {
		const real eps = std::sqrt(std::numeric_limits<real>::epsilon());
		bool strict = std::is_same<real, double>::value;
		GEO::Logger::out("Epsilon") << "Using epsilon = " << eps << std::endl;
		try {
			for (int idx = 0; idx < voxels.num_voxels(); ++idx) {
				Vec3i vox = voxels.index3_from_index(idx);
				GEO::vec3 pos = voxels.voxel_center(vox[0], vox[1], vox[2]);
				const double sq_dist = aabb_tree.squared_distance(pos);
				if (strict && voxels.at(idx) != (float) std::sqrt(sq_dist)) {
					std::cout << voxels.at(idx) << " " << (float) std::sqrt(sq_dist)
						<< ' ' << voxels.at(idx) - (float) std::sqrt(sq_dist) << std::endl;
				}
				geo_assert(std::abs(voxels.at(idx) - (float) std::sqrt(sq_dist)) < eps);
			}
		} catch(const GEO::TaskCanceled&) {
			// Do early cleanup
		}
	}
}

// -----------------------------------------------------------------------------

// Next power of 2 that is strictly > x
int next_pow2(int x) {
	int n = 1;
	while (n <= x) { n <<= 1; }
	return n;
}

Vec3i morton_to_index3(size_t x) {
	Vec3i pos{0, 0, 0};
	size_t base = 1;
	while (x) {
		pos[0] += base * (x & 1); x >>= 1;
		pos[1] += base * (x & 1); x >>= 1;
		pos[2] += base * (x & 1); x >>= 1;
		base <<= 1;
	}
	return pos;
}

// -----------------------------------------------------------------------------

void compute_unsigned_distance_field_cpu(const GEO::Mesh &M,
	const GEO::MeshFacetsAABB &aabb_tree, VoxelGrid &voxels)
{

	try {
		GEO::ProgressTask task("Sqdist (CPU)", 100);

	#if 0
		#pragma omp parallel for
		for (int idx = 0; idx < voxels.num_voxels(); ++idx) {
			if (omp_get_thread_num() == 0) {
				task.progress((int) (100.0 * idx / voxels.num_voxels() * omp_get_num_threads()));
			}

			Vec3i vox = voxels.index3_from_index(idx);
			GEO::vec3 pos = voxels.voxel_center(vox[0], vox[1], vox[2]);
			double sq_dist = aabb_tree.squared_distance(pos);
			voxels.at(idx) = (float) std::sqrt(sq_dist);
		}
	#elif 1
		GEO::index_t prev_facet = GEO::NO_FACET;
		double sq_dist = std::numeric_limits<double>::max();
		GEO::vec3 nearest_point;
		#pragma omp parallel for private(prev_facet, sq_dist, nearest_point)
		for (int idx = 0; idx < voxels.num_voxels(); ++idx) {
			if (omp_get_thread_num() == 0) {
				task.progress((int) (100.0 * idx / voxels.num_voxels() * omp_get_num_threads()));
			}

			Vec3i vox = voxels.index3_from_index(idx);
			GEO::vec3 query = voxels.voxel_center(vox[0], vox[1], vox[2]);
			if (prev_facet != GEO::NO_FACET) {
				GEO::get_point_facet_nearest_point(M, query, prev_facet, nearest_point, sq_dist);
			}
			aabb_tree.nearest_facet_with_hint(query, prev_facet, nearest_point, sq_dist);
			voxels.at(idx) = (float) std::sqrt(sq_dist);
		}
	#else

		Vec3i num_voxels = voxels.grid_size();
		size_t upper = next_pow2(std::max(num_voxels[0], std::max(num_voxels[1], num_voxels[2])));
		size_t upper3 = upper * upper * upper;

		std::cout << upper << std::endl;

		GEO::index_t prev_facet = GEO::NO_FACET;
		double sq_dist = std::numeric_limits<double>::max();
		GEO::vec3 nearest_point;
		#pragma omp parallel for private(prev_facet, sq_dist, nearest_point)
		for (size_t idx = 0; idx < upper3; ++idx) {
			if (omp_get_thread_num() == 0) {
				task.progress((int) (100.0 * idx / upper3 * omp_get_num_threads()));
			}

			Vec3i vox = morton_to_index3(idx);
			if (vox[0] >= num_voxels[0] || vox[1] >= num_voxels[1] || vox[2] >= num_voxels[2]) {
				continue;
			}
			GEO::vec3 query = voxels.voxel_center(vox[0], vox[1], vox[2]);
			if (prev_facet != GEO::NO_FACET) {
				GEO::get_point_facet_nearest_point(M, query, prev_facet, nearest_point, sq_dist);
			}
			aabb_tree.nearest_facet_with_hint(query, prev_facet, nearest_point, sq_dist);
			voxels.at(voxels.index_from_index3(vox)) = (float) std::sqrt(sq_dist);
		}


	#endif

	} catch(const GEO::TaskCanceled&) {
		// Do early cleanup
	}
}

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
	geo_assert(sum!=0); // if the SOS signs match and are nonkero, there's no way all of a, b, and c are zero.
	a/=sum;
	b/=sum;
	c/=sum;
	return true;
}

// -----------------------------------------------------------------------------

bool intersect_ray_z(const GEO::Mesh &M, GEO::index_t f, const GEO::vec3 &q, double &z) {
	using namespace GEO;

	index_t c = M.facets.corners_begin(f);
	const vec3& p1 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));
	++c;
	const vec3& p2 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));
	++c;
	const vec3& p3 = Geom::mesh_vertex(M, M.facet_corners.vertex(c));

	double u, v, w;
	if (point_in_triangle_2d(
		q[0], q[1], p1[0], p1[1], p2[0], p2[1], p3[0], p3[1], u, v, w))
	{
		z = u*p1[2] + v*p2[2] + w*p3[2];
		return true;
	}

	return false;
}

// -----------------------------------------------------------------------------

void compute_sign(const GEO::Mesh &M,
	const GEO::MeshFacetsAABB &aabb_tree, VoxelGrid &voxels)
{
	const Vec3i size = voxels.grid_size();

	try {
		GEO::ProgressTask task("Ray marching", 100);

		GEO::vec3 min_corner, max_corner;
		GEO::get_bbox(M, &min_corner[0], &max_corner[0]);

		const GEO::vec3 origin = voxels.origin();
		const double spacing = voxels.spacing();

		#pragma omp parallel for
		for (int x = 0; x < size[0]; ++x) {
			if (omp_get_thread_num() == 0) {
				task.progress((int) (100.0 * x / size[0] * omp_get_num_threads()));
			}
			for (int y = 0; y < size[1]; ++y) {
				GEO::vec3 center = voxels.voxel_center(x, y, 0);

				GEO::Box box;
				box.xyz_min[0] = box.xyz_max[0] = center[0];
				box.xyz_min[1] = box.xyz_max[1] = center[1];
				box.xyz_min[2] = min_corner[2];
				box.xyz_max[2] = max_corner[2];

				std::vector<double> inter;
				auto action = [&M, &inter, &center] (GEO::index_t f) {
					double z;
					if (intersect_ray_z(M, f, center, z)) {
						inter.push_back(z);
					}
				};
				aabb_tree.compute_bbox_facet_bbox_intersections(box, action);
				std::sort(inter.begin(), inter.end());

				for (size_t k = 1; k < inter.size(); k += 2) {
					int z1 = int(std::round((inter[k-1] - origin[2])/spacing));
					int z2 = int(std::round((inter[k] - origin[2])/spacing));
					for (int z = z1; z < z2; ++z) {
						if (z >= 0 && z < size[2]) {
							const int idx = voxels.index_from_index3({{x, y, z}});
							voxels.at(idx) *= -1.0f;
						}
					}
				}
			}
		}
	} catch(const GEO::TaskCanceled&) {
		// Do early cleanup
	}
}

////////////////////////////////////////////////////////////////////////////////

void sample_inside(const VoxelGrid &voxels, std::vector<GEO::vec3> &pts, int num_samples) {
	// TODO: Implement
}

////////////////////////////////////////////////////////////////////////////////

void write_samples(std::string &basename, const std::vector<GEO::vec3> &pts) {
	std::ofstream out(basename + ".xyz");
	out.precision(10);
	for (GEO::vec3 x : pts) {
		out << x[0] << ' ' << x[1] << ' ' << x[2] << '\n';
	}
}

// -----------------------------------------------------------------------------

void paraview_dump(std::string &basename, const VoxelGrid &voxels) {
	Vec3i size = voxels.grid_size();

    std::ofstream metafile(basename + ".mhd");
    metafile << "ObjectType = Image\nNDims = 3\n"
    	<< "DimSize = " << size[0] << " " << size[1] << " " << size[2] << "\n"
    	<< "ElementType = MET_FLOAT\nElementDataFile = " + basename + ".raw\n";
    metafile.close();

	std::ofstream rawfile(basename + ".raw", std::ios::binary);
	rawfile.write(reinterpret_cast<const char*>(voxels.rawbuf ()), voxels.num_voxels() * sizeof(float));
	rawfile.close();
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
	using namespace GEO;

	GEO::initialize();

	try {

		Stopwatch W("Total time");

		// Program arguments
		set_programm_options();
		std::vector<std::string> filenames;
		if (!CmdLine::parse(argc, argv, filenames, "meshfile <outputfile>")) {
			return 1;
		}

		int padding = std::max(1, CmdLine::get_arg_int("padding"));
		double voxel_size = CmdLine::get_arg_double("resolution");
		bool voxelize_only = CmdLine::get_arg_bool("voxelize_only");
		//bool sample_inside = CmdLine::get_arg_bool("sample_inside");
		//int num_samples = CmdLine::get_arg_int("num_samples");

		// IO mesh names
		std::string mesh_filename = filenames[0];
		std::string output_basename;
		if (filenames.size() > 1) {
			output_basename = filenames[1];
		} else {
			output_basename = "output";
		}

		Logger::out("I/O") << "Output = " << output_basename << std::endl;

		// Load data
		Logger::div("Loading data");

		Mesh M_in;
		if (!mesh_load(mesh_filename, M_in)) {
			return 1;
		}
		geo_assert(M_in.vertices.dimension() == 3);

		// Initialize voxel grid and AABB tree
		vec3 min_corner, max_corner;
		GEO::get_bbox(M_in, &min_corner[0], &max_corner[0]);
		VoxelGrid voxels(min_corner, max_corner - min_corner, voxel_size, padding);
		GEO::MeshFacetsAABB aabb_tree(M_in);

		// Compute (unsigned) distance field
		if (!voxelize_only) {
			Logger::div("Computing (unsigned) distance field");
			if (GEO::CmdLine::get_arg_bool("use_gpu")) {
				if (GEO::CmdLine::get_arg_bool("opencl_use_float")) {
					compute_unsigned_distance_field_gpu<float>(M_in, aabb_tree, voxels);
				} else {
					compute_unsigned_distance_field_gpu<double>(M_in, aabb_tree, voxels);
				}
			} else {
				compute_unsigned_distance_field_cpu(M_in, aabb_tree, voxels);
			}
		}

		// Compute inside/outside info
		Logger::div("Computing inside/outside info");
		compute_sign(M_in, aabb_tree, voxels);

		// Sample points inside the voxels
		//std::vector<vec3> pts;
		//if (sample_inside) {
			//pts.reserve(num_samples);
		//}

		// Save resulting voxel grid
		Logger::div("Saving result");
		// if (sample_inside) {
		// 	write_samples(output_basename, pts);
		// } else {
			paraview_dump(output_basename, voxels);
		// }

	} catch (const std::exception& e) {
		std::cerr << "Received an exception: " << e.what() << std::endl;
		return 1;
	}

	Logger::out("") << "Everything OK, Returning status 0" << std::endl;
	return 0;
}
