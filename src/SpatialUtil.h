/*

Copyright (c) 2005-2016, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Aboria.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef SPATIAL_UTILS_H_
#define SPATIAL_UTILS_H_

namespace Aboria {
template <unsigned int D> struct bbox;
}

#include "detail/SpatialUtil.h"

namespace Aboria {

///
/// @brief Contains the minimum and maximum extents of a hypercube in @p D
/// dimensional space
///
/// @tparam D the number of spatial dimensions
///
template <unsigned int D> struct bbox {
  typedef Vector<double, D> double_d;

  ///
  /// @brief minimum point in the box (i.e. lower left corner for D=2)
  ///
  double_d bmin;

  ///
  /// @brief maximum point in the box (i.e. upper right corner for D=2)
  ///
  double_d bmax;

  inline CUDA_HOST_DEVICE bbox()
      : bmin(double_d::Constant(detail::get_max<double>())),
        bmax(double_d::Constant(-detail::get_max<double>())) {}

  inline CUDA_HOST_DEVICE bbox(const double_d &p) : bmin(p), bmax(p) {}

  inline CUDA_HOST_DEVICE bbox(const double_d &min, const double_d &max)
      : bmin(min), bmax(max) {}

  ///
  /// @return the bounding box covering both input boxes
  ///
  inline CUDA_HOST_DEVICE bbox operator+(const bbox &arg) {
    bbox bounds;
    for (size_t i = 0; i < D; ++i) {
      bounds.bmin[i] = std::min(bmin[i], arg.bmin[i]);
      bounds.bmax[i] = std::max(bmax[i], arg.bmax[i]);
    }
    return bounds;
  }

  ///
  /// @return true if lhs box is within rhs box
  ///
  inline CUDA_HOST_DEVICE bool operator<(const bbox &arg) {
    bbox bounds;
    bool within = true;
    for (size_t i = 0; i < D; ++i) {
      within |= bmin[i] >= arg.bmin[i];
      within |= bmax[i] < arg.bmax[i];
    }
    return within;
  }

  ///
  /// @return true if lhs box is the same or within rhs box
  ///
  inline CUDA_HOST_DEVICE bool operator<=(const bbox &arg) {
    bbox bounds;
    bool within = true;
    for (size_t i = 0; i < D; ++i) {
      within |= bmin[i] >= arg.bmin[i];
      within |= bmax[i] <= arg.bmax[i];
    }
    return within;
  }

  ///
  /// @return true if box has no volume
  ///
  inline CUDA_HOST_DEVICE bool is_empty() {
    for (size_t i = 0; i < D; ++i) {
      if (bmax[i] < bmin[i])
        return true;
    }
    return false;
  }
};

///
/// @brief print bbox to a stream
///
/// @tparam D the number of spatial dimensions
/// @param out the stream
/// @param b the box to print
///
template <unsigned int D>
std::ostream &operator<<(std::ostream &out, const bbox<D> &b) {
  return out << "bbox(" << b.bmin << "<->" << b.bmax << ")";
}

/// @returns true if the hypersphere defined by \p centre and \p radius
/// intersects with the 1D line defined by the two points \p a and \p b
template <unsigned int D>
bool circle_intersect_line(const Vector<double, D> &centre, const double radius,
                           const Vector<double, D> &a,
                           const Vector<double, D> &b) {

  const Vector<double, D> a_minus_c = a - centre;
  const Vector<double, D> dx = b - a;
  const double dx_norm = dx.norm();
  const double a_minus_c_along_ab = dx.dot(a_minus_c) / dx_norm;

  const double under_sqrt = std::pow(a_minus_c_along_ab, 2) -
                            a_minus_c.squaredNorm() + std::pow(radius, 2);

  bool intersects;
  if (under_sqrt > 0) {
    // two intersections
    const double after_sqrt = std::sqrt(under_sqrt);
    const double d1 = -a_minus_c_along_ab - after_sqrt;
    const double d2 = -a_minus_c_along_ab + after_sqrt;
    intersects = (d1 >= 0 && d1 < dx_norm) || (d2 >= 0 && d2 < dx_norm);
  } else if (under_sqrt == 0) {
    // one intersection
    const double d1 = -a_minus_c_along_ab;
    intersects = d1 >= 0 && d1 < dx_norm;
  } else {
    // no intersections
    intersects = false;
  }

  // check if one of the points is within the sphere
  bool a_in_sphere = (a - centre).squaredNorm() <= std::pow(radius, 2);

  return intersects || a_in_sphere;
}

/// @returns true if the hypersphere defined by \p centre and \p radius
/// intersects with the hypercube defined by \p cube
template <unsigned int D>
bool circle_intersect_cube(const Vector<double, D> &centre, const double radius,
                           const bbox<D> &cube) {

  // loop through vertices of cube and test all the possible edges
  bool an_edge_intersects = false;
  for (int i = 0; i < (1 << D); ++i) {
    // create vertex a
    Vector<double, D> a;
    for (size_t j = 0; j < D; ++j) {
      const bool jth_bit_of_i = i & (1 << j);
      a[j] = jth_bit_of_i ? cube.bmax[j] : cube.bmin[j];
    }
    // number of edges from each vertex is equal to the number of dimensions
    for (size_t j = 0; j < D; ++j) {
      const bool jth_bit_of_i = i & (1 << j);
      // b is equal to a except for the jth bit of i flipped
      // only count if jth bit of i is 0
      if (!jth_bit_of_i) {
        Vector<double, D> b = a;
        b[j] = cube.bmax[j];
        an_edge_intersects |= circle_intersect_line(centre, radius, a, b);
      }
    }
  }

  // check if centre is within cube
  bool centre_in_cube = true;
  for (size_t i = 0; i < D; ++i) {
    centre_in_cube &= centre[i] >= cube.bmin[i] && centre[i] < cube.bmax[i];
  }

  // intersection if either an edge intersects the sphere, or the centre of
  // the sphere is within the hypercube
  return an_edge_intersects || centre_in_cube;
}

} // namespace Aboria

#endif
