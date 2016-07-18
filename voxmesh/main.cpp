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

#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/progress.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_AABB.h>
#include <algorithm>

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
	GEO::Logger::out("Voxels") << "Grid size: "
		<< m_grid_size[0] << " x " << m_grid_size[1] << " x " << m_grid_size[2] << std::endl;
	m_data.assign(m_grid_size[0] * m_grid_size[1] * m_grid_size[2], 1.0f);
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

////////////////////////////////////////////////////////////////////////////////

bool intersect_ray_z(const GEO::Mesh &M, GEO::index_t f, const GEO::vec3 &q, double &z) {
	using namespace GEO;

	index_t c = M.facets.corners_begin(f);
	const vec3& p1 = Geom::mesh_vertex(M, M.facet_corners.vertex(c++));
	const vec3& p2 = Geom::mesh_vertex(M, M.facet_corners.vertex(c++));
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

		GEO::parallel_for([&](int x) {
			if (GEO::Thread::current()->id() == 0) {
				task.progress((int) (100.0 * x / size[0] * GEO::Process::number_of_cores()));
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
		}, 0, size[0]);

		for (int idx = 0; idx < voxels.num_voxels(); ++idx) {
			// -1 is inside, +1 is outside
			voxels.at(idx) = (voxels.at(idx) < 0 ? 1 : 0);
		}
	} catch(const GEO::TaskCanceled&) {
		// Do early cleanup
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
	// Initialize the Geogram library.
	GEO::initialize();

	// Import standard command line arguments, and custom ones.
	GEO::CmdLine::import_arg_group("standard");
	GEO::CmdLine::declare_arg("padding", 2, "number of padded grid cells");
	GEO::CmdLine::declare_arg("resolution", 1.0, "size of a voxel (in mm)");

	// Parse command line options and filenames.
	std::vector<std::string> filenames;
	if(!GEO::CmdLine::parse(argc, argv, filenames, "in_mesh_file <out_voxel_file>")) {
		return 1;
	}

	int padding = std::max(1, GEO::CmdLine::get_arg_int("padding"));
	double voxel_size = GEO::CmdLine::get_arg_double("resolution");

	// Default output filename is "output" if unspecified.
	if(filenames.size() == 1) {
		filenames.push_back("output");
	}

	// Display input and output filenames.
	GEO::Logger::div("Command line");
	GEO::Logger::out("VoxMesh") << "Input file: " << filenames[0] << std::endl;
	GEO::Logger::out("VoxMesh") << "Output file: " << filenames[1] << std::endl;

	// Declare a mesh.
	GEO::Mesh M;

	// Load the mesh and display timings.
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
	VoxelGrid voxels(min_corner, max_corner - min_corner, voxel_size, padding);
	GEO::MeshFacetsAABB aabb_tree(M);

	// Compute inside/outside info
	GEO::Logger::div("Voxelizing");
	compute_sign(M, aabb_tree, voxels);

	// Save voxel grid and display timings.
	GEO::Logger::div("Saving");
	{
		GEO::Stopwatch W("Save");
		paraview_dump(filenames[1], voxels);
	}

	return 0;
}

