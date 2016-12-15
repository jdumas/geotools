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

#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/progress.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_io.h>

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
	// Initialize the Geogram library
	GEO::initialize();

	// Import standard command line arguments, and custom ones
	GEO::CmdLine::import_arg_group("standard");
	GEO::CmdLine::declare_arg("extent", 1.0, "Maximum physical extent (in mm)");

	// Parse command line options and filenames
	std::vector<std::string> filenames;
	if (!GEO::CmdLine::parse(argc, argv, filenames, "in_mesh_file <out_mesh_file>")) {
		return 1;
	}

	double target_extent = GEO::CmdLine::get_arg_double("extent");

	// Default output filename is "output" if unspecified
	if (filenames.size() == 1) {
		filenames.push_back("output.ply");
	}

	// Display input and output filenames
	GEO::Logger::div("Command line");
	GEO::Logger::out("Normalize") << "Input file: " << filenames[0] << std::endl;
	GEO::Logger::out("Normalize") << "Output file: " << filenames[1] << std::endl;

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

	// Merge identical vertices
	// GEO::mesh_repair(M, GEO::MESH_REPAIR_COLOCATE);

	// Rescale to unit box, and set min corner to 0
	GEO::vec3 min_corner, max_corner;
	GEO::get_bbox(M, &min_corner[0], &max_corner[0]);
	GEO::vec3 extent = max_corner - min_corner;
	double scaling = std::max(extent[0], std::max(extent[1], extent[2]));
	for (int v = 0; v < M.vertices.nb(); ++v) {
		M.vertices.point(v) = target_extent * (M.vertices.point(v) - min_corner) / scaling;
	max_corner = min_corner + 0.5 * extent;

	// -------------------------------------------------------------------------
	// Extract submesh
	// -------------------------------------------------------------------------

	// GEO::Box bbox;
	// for (int i = 0; i < 3; ++i) {
	// 	bbox.xyz_min[i] = min_corner[i];
	// 	bbox.xyz_max[i] = max_corner[i];
	// }
	// GEO::vector<GEO::index_t> to_delete(M.facets.nb(), 0);
	// for (int f = 0; f < M.facets.nb(); ++f) {
	// 	for(GEO::index_t c = M.facets.corners_begin(f);
	// 		c < M.facets.corners_end(f); ++c)
	// 	{
	// 		GEO::vec3 pts = M.vertices.point(M.facet_corners.vertex(c));
	// 		if (! bbox.contains(pts)) {
	// 			to_delete[f] = 1;
	// 			break;
	// 		}
	// 	}
	// }

	// M.facets.delete_elements(to_delete, true);

	// -------------------------------------------------------------------------

	// Save mesh
	GEO::Logger::div("Saving");
	{
		GEO::Stopwatch W("Save");
		GEO::mesh_save(M, filenames[1]);
	}

	return 0;
}

