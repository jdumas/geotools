#include "common.cl"

real sqr(const real a) {
	return a * a;
}

real distance2(const real3 a, const real3 b) {
	return dot(b - a, b - a);
}

real length2(const real3 a) {
	return dot(a, a);
}

/**
 * \brief Computes the squared distance between a point and a Box.
 * \param[in] p the point
 * \param[in] B the box
 * \return the squared distance between \p p and \p B
 * \pre p is inside B
 */
real inner_point_box_squared_distance(
	const real3 pos, const real3 xyz_min, const real3 xyz_max)
{
	real result = sqr(pos.x - xyz_min.x);
	result = min(result, sqr(pos.x - xyz_max.x));
	result = min(result, sqr(pos.y - xyz_min.y));
	result = min(result, sqr(pos.y - xyz_max.y));
	result = min(result, sqr(pos.z - xyz_min.z));
	result = min(result, sqr(pos.z - xyz_max.z));
	return result;
}

/**
 * \brief Computes the squared distance between a point and a Box
 *  with negative sign if the point is inside the Box.
 * \param[in] p the point
 * \param[in] B the box
 * \return the signed squared distance between \p p and \p B
 */
real point_box_signed_squared_distance(
	const real3 pos, const real3 xyz_min, const real3 xyz_max)
{
	// TODO: vectorize
	bool inside = true;
	real result = 0;

	// x coord
	if (pos.x < xyz_min.x) {
		inside = false; result += sqr(pos.x - xyz_min.x);
	} else if (pos.x > xyz_max.x) {
		inside = false; result += sqr(pos.x - xyz_max.x);
	}
	// y coord
	if (pos.y < xyz_min.y) {
		inside = false; result += sqr(pos.y - xyz_min.y);
	} else if (pos.y > xyz_max.y) {
		inside = false; result += sqr(pos.y - xyz_max.y);
	}
	// z coord
	if (pos.z < xyz_min.z) {
		inside = false; result += sqr(pos.z - xyz_min.z);
	} else if (pos.z > xyz_max.z) {
		inside = false; result += sqr(pos.z - xyz_max.z);
	}

	if (inside) {
		result = -inner_point_box_squared_distance(pos, xyz_min, xyz_max);
	}
	return result;
}

/**
 * \brief Computes the squared distance between a point and the
 *  center of a box.
 * \param[in] p the point
 * \param[in] B the box
 * \return the squared distance between \p p and the center of \p B
 */
real point_box_center_squared_distance(
	const real3 pos, const real3 xyz_min, const real3 xyz_max)
{
	const real3 diff = pos - (real) (0.5) * (xyz_min + xyz_max);
	return dot(diff, diff);
}

/**
 * \brief Computes the point closest to a given point in a nd triangle
 * \details See
 *  http://www.geometrictools.com/LibMathematics/Distance/Distance.html
 * \param[in] point the query point
 * \param[in] V0 first vertex of the triangle
 * \param[in] V1 second vertex of the triangle
 * \param[in] V2 third vertex of the triangle
 * \param[out] closest_point the point closest to \p point in the
 *  triangle (\p V0, \p V1, \p V2)
 * \param[out] lambda0 barycentric coordinate of the closest point
 *  relative to \p V0
 * \param[out] lambda1 barycentric coordinate of the closest point
 *  relative to \p V1
 * \param[out] lambda2 barycentric coordinate of the closest point
 *  relative to \p V2
 * \tparam VEC the class that represents the points.
 * \return the squared distance between the point and
 *  the triangle (\p V0, \p V1, \p V2)
 */
inline real point_triangle_squared_distance(
	const real3 point, const real3 V0, const real3 V1, const real3 V2)
{
	real3 diff = V0 - point;
	real3 edge0 = V1 - V0;
	real3 edge1 = V2 - V0;
	real a00 = length2(edge0);
	real a01 = dot(edge0, edge1);
	real a11 = length2(edge1);
	real b0 = dot(diff, edge0);
	real b1 = dot(diff, edge1);
	real c = length2(diff);
	real det = fabs(a00 * a11 - a01 * a01);
	real s = a01 * b1 - a11 * b0;
	real t = a01 * b0 - a00 * b1;
	real sqrDistance = 0;

	if(s + t <= det) {
		if(s < 0.0) {
			if(t < 0.0) {   // region 4
				if(b0 < 0.0) {
					t = 0.0;
					if(-b0 >= a00) {
						s = 1.0;
						sqrDistance = a00 + 2.0 * b0 + c;
					} else {
						s = -b0 / a00;
						sqrDistance = b0 * s + c;
					}
				} else {
					s = 0.0;
					if(b1 >= 0.0) {
						t = 0.0;
						sqrDistance = c;
					} else if(-b1 >= a11) {
						t = 1.0;
						sqrDistance = a11 + 2.0 * b1 + c;
					} else {
						t = -b1 / a11;
						sqrDistance = b1 * t + c;
					}
				}
			} else {  // region 3
				s = 0.0;
				if(b1 >= 0.0) {
					t = 0.0;
					sqrDistance = c;
				} else if(-b1 >= a11) {
					t = 1.0;
					sqrDistance = a11 + 2.0 * b1 + c;
				} else {
					t = -b1 / a11;
					sqrDistance = b1 * t + c;
				}
			}
		} else if(t < 0.0) {  // region 5
			t = 0.0;
			if(b0 >= 0.0) {
				s = 0.0;
				sqrDistance = c;
			} else if(-b0 >= a00) {
				s = 1.0;
				sqrDistance = a00 + 2.0 * b0 + c;
			} else {
				s = -b0 / a00;
				sqrDistance = b0 * s + c;
			}
		} else {  // region 0
			// minimum at interior point
			real invDet = (real) (1.0) / det;
			s *= invDet;
			t *= invDet;
			sqrDistance = s * (a00 * s + a01 * t + 2.0 * b0) +
				t * (a01 * s + a11 * t + 2.0 * b1) + c;
		}
	} else {
		real tmp0, tmp1, numer, denom;

		if(s < 0.0) {   // region 2
			tmp0 = a01 + b0;
			tmp1 = a11 + b1;
			if(tmp1 > tmp0) {
				numer = tmp1 - tmp0;
				denom = a00 - 2.0 * a01 + a11;
				if(numer >= denom) {
					s = 1.0;
					t = 0.0;
					sqrDistance = a00 + 2.0 * b0 + c;
				} else {
					s = numer / denom;
					t = 1.0 - s;
					sqrDistance = s * (a00 * s + a01 * t + 2.0 * b0) +
						t * (a01 * s + a11 * t + 2.0 * b1) + c;
				}
			} else {
				s = 0.0;
				if(tmp1 <= 0.0) {
					t = 1.0;
					sqrDistance = a11 + 2.0 * b1 + c;
				}
				else if(b1 >= 0.0) {
					t = 0.0;
					sqrDistance = c;
				} else {
					t = -b1 / a11;
					sqrDistance = b1 * t + c;
				}
			}
		} else if(t < 0.0) {  // region 6
			tmp0 = a01 + b1;
			tmp1 = a00 + b0;
			if(tmp1 > tmp0) {
				numer = tmp1 - tmp0;
				denom = a00 - 2.0 * a01 + a11;
				if(numer >= denom) {
					t = 1.0;
					s = 0.0;
					sqrDistance = a11 + 2.0 * b1 + c;
				} else {
					t = numer / denom;
					s = 1.0 - t;
					sqrDistance = s * (a00 * s + a01 * t + 2.0 * b0) +
						t * (a01 * s + a11 * t + 2.0 * b1) + c;
				}
			} else {
				t = 0.0;
				if(tmp1 <= 0.0) {
					s = 1.0;
					sqrDistance = a00 + 2.0 * b0 + c;
				} else if(b0 >= 0.0) {
					s = 0.0;
					sqrDistance = c;
				} else {
					s = -b0 / a00;
					sqrDistance = b0 * s + c;
				}
			}
		} else { // region 1
			numer = a11 + b1 - a01 - b0;
			if(numer <= 0.0) {
				s = 0.0;
				t = 1.0;
				sqrDistance = a11 + 2.0 * b1 + c;
			} else {
				denom = a00 - 2.0 * a01 + a11;
				if(numer >= denom) {
					s = 1.0;
					t = 0.0;
					sqrDistance = a00 + 2.0 * b0 + c;
				} else {
					s = numer / denom;
					t = 1.0 - s;
					sqrDistance = s * (a00 * s + a01 * t + 2.0 * b0) +
						t * (a01 * s + a11 * t + 2.0 * b1) + c;
				}
			}
		}
	}

	// Account for numerical round-off error.
	if (sqrDistance < 0.0) {
		sqrDistance = 0.0;
	}

	return sqrDistance;
}
