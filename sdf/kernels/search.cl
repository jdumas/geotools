#include "common.cl"
#include "geom.cl"

////////////////////////////////////////////////////////////////////////////////

real3 get_point(const __global real * corners, index_t idx) {
	return (real3) (corners[3*idx+0], corners[3*idx+1], corners[3*idx+2]);
}

uint get_corner(const __global uint *facets, index_t f, index_t c) {
	return facets[3 * f + c];
}

real get_point_facet_sq_dist(
	const __global real * vertices,
	const __global uint * facets,
	const real3 pos,
	const index_t f)
{
	const real3 p1 = get_point(vertices, get_corner(facets, f, 0));
	const real3 p2 = get_point(vertices, get_corner(facets, f, 1));
	const real3 p3 = get_point(vertices, get_corner(facets, f, 2));
	return point_triangle_squared_distance(pos, p1, p2, p3);
}

////////////////////////////////////////////////////////////////////////////////

real get_nearest_facet_hint(
	const __global real * vertices,
	const __global uint * facets,
	const __global real * min_corners,
	const __global real * max_corners,
	const uint nb_facets,
	const real3 pos)
{
	// Find a good initial value for nearest_f by traversing
	// the boxes and selecting the child such that the center
	// of its bounding box is nearer to the query point.
	// For a large mesh (20M facets) this gains up to 10%
	// performance as compared to picking nearest_f randomly.
	index_t b = 0;
	index_t e = nb_facets;
	index_t n = 1;
	while (e != b + 1) {
		const index_t m = b + (e - b) / 2;
		const index_t childl = 2 * n;
		const index_t childr = 2 * n + 1;
		const real dl = point_box_center_squared_distance(
			pos, get_point(min_corners, childl), get_point(max_corners, childl));
		const real dr = point_box_center_squared_distance(
			pos, get_point(min_corners, childr), get_point(max_corners, childr));
		if (dl < dr) {
			e = m;
			n = childl;
		} else {
			b = m;
			n = childr;
		}
	}
	const index_t nearest_f = b;
	const index_t v = get_corner(facets, nearest_f, 0);
	const real3 nearest_point = get_point(vertices, v);
	return distance2(pos, nearest_point);
}

////////////////////////////////////////////////////////////////////////////////

bool going_up(index_t *node, uint *order, uint *level) {
	while (*node > 1 && *node % 2 == (*order & 1)) {
		*node  = *node / 2;
		*order = (*order >> 1);
		*level = *level - 1;
	}
	if (*node == 1) {
		return true;
	} else {
		if (*order & 1) {
			*node = *node + 1;
		} else {
			*node = *node - 1;
		}
	}
	return false;
}

// -----------------------------------------------------------------------------

uint2 get_interval(index_t b, index_t e, index_t node, uint level) {
	for (uint depth = level; depth > 0; --depth) {
		const index_t m = b + (e - b) / 2;
		if (node & (1u << (depth - 1u))) {
			b = m;
		} else {
			e = m;
		}
	}
	return (uint2) (b, e);
}

// -----------------------------------------------------------------------------

real nearest_facet_iterative(
	const __global real * vertices,
	const __global uint * facets,
	const __global real * min_corners,
	const __global real * max_corners,
	const real3 pos, real sq_dist, const index_t e)
{
	index_t node = 1;
	uint order = 0;
	uint level = 0;
	uint2 s = (uint2) (0, e);
	bool up = false;
	while (true) {
		up = false;
		if(s.s0 + 1 == s.s1) {
			// If node is a leaf: compute point-facet distance
			// and replace current if nearer
			const real cur_sq_dist = get_point_facet_sq_dist(
				vertices, facets, pos, s.s0);
			if (cur_sq_dist < sq_dist) {
				sq_dist = cur_sq_dist;
			}

			// Going up
			if (going_up(&node, &order, &level)) { break; }
			up = true;
		} else {
			// Internal node: check bbox distance with current best, and skip accordingly

			// Go to the next child
			const index_t childl = 2 * node;
			const index_t childr = 2 * node + 1;

			const real dl = point_box_signed_squared_distance(
				pos, get_point(min_corners, childl), get_point(max_corners, childl));
			const real dr = point_box_signed_squared_distance(
				pos, get_point(min_corners, childr), get_point(max_corners, childr));

			bool go_left = false;
			bool go_right = false;
			if (dl < dr) {
				if (dl < sq_dist) {
					go_left = true;
				} else if (dr < sq_dist) {
					go_left = true;
				}
			} else {
				if (dr < sq_dist) {
					go_right = true;
				} else if (dl < sq_dist) {
					go_left = true;
				}
			}

			if (go_left) {
				node = childl;
				order = (order << 1) | (1u);
				++level;
				s.s1 = s.s0 + (s.s1 - s.s0) / 2;
			} else if (go_right) {
				node = childr;
				order = (order << 1);
				++level;
				s.s0 = s.s0 + (s.s1 - s.s0) / 2;
			} else {
				if (going_up(&node, &order, &level)) { break; }
				up = true;
			}
		}

		// When going up the tree, test if we can prune the new branch
		if (up) {
			bool stop = false;
			while (true) {
				real dn = point_box_signed_squared_distance(
					pos, get_point(min_corners, node), get_point(max_corners, node));
				if (dn < sq_dist) {
					break;
				} else {
					if (going_up(&node, &order, &level)) {
						stop = true;
						break;
					}
				}
			}
			if (stop) { break; }
			s = get_interval(0, e, node, level);
		}
	}
	return sq_dist;
}

// -----------------------------------------------------------------------------

real brute_force_test(
	const __global real * vertices,
	const __global uint * facets,
	const __global real * min_corners,
	const __global real * max_corners,
	const real3 pos, const index_t nb_facets)
{
	real sq_dist = get_point_facet_sq_dist(vertices, facets, pos, 0);
	for (index_t f = 0; f < nb_facets; ++f) {
		const real cur_sq_dist = get_point_facet_sq_dist(
			vertices, facets, pos, f);
		if (cur_sq_dist < sq_dist) {
			sq_dist = cur_sq_dist;
		}
	}
	return sq_dist;
}

// -----------------------------------------------------------------------------

real squared_dist_point(
	const __global real * vertices,
	const __global uint * facets,
	const __global real * min_corners,
	const __global real * max_corners,
	const real3 pos, const uint nb_facets)
{
	const real sq_hint = get_nearest_facet_hint(
		vertices, facets, min_corners, max_corners, nb_facets, pos);
	return nearest_facet_iterative(
		vertices, facets, min_corners, max_corners,
		pos, sq_hint, nb_facets);
}

////////////////////////////////////////////////////////////////////////////////

real3 voxel_center(const uint3 coord, const real3 origin, const real spacing) {
	return (to_real3 (coord) + (real3) (0.5)) * spacing + origin;
}

// -----------------------------------------------------------------------------

__kernel void squared_dist(
	const __global real * vertices,
	const __global uint * facets,
	const __global real * min_corners,
	const __global real * max_corners,
	global real * layer_sq_dist,
	const uint3 layer_size,
	const uint3 layer_offset,
	const real3 origin,
	const real spacing,
	const uint nb_facets)
{
	for (size_t z = get_global_id(2); z < layer_size.z; z += get_global_size(2)) {
		for (size_t y = get_global_id(1); y < layer_size.y; y += get_global_size(1)) {
			for (size_t x = get_global_id(0); x < layer_size.x; x += get_global_size(0)) {
				const uint3 layer_coord = (uint3) (x, y, z) + layer_offset;
				const real3 pos = voxel_center(layer_coord, origin, spacing);
				//const real sq_dist = brute_force_test(
				const real sq_dist = squared_dist_point(
					vertices, facets, min_corners, max_corners, pos, nb_facets);
				const size_t idx = x + layer_size.x * (y + layer_size.y * z);
				layer_sq_dist[idx] = sq_dist;
			}
		}
	}
}
