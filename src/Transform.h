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

#ifndef TRANSFORM_H_
#define TRANSFORM_H_

#include "LatticeIterator.h"
#include "SpatialUtil.h"
#include "Vector.h"

namespace Aboria {

/// an identity coordinate transform
struct IdentityTransform {
  /// transforms the point @p the the new coordinate system, returning the
  /// result
  template <unsigned int D>
  inline Vector<double, D> operator()(const Vector<double, D> &v) const {
    return v;
  }
  /// transforms the box @p to the new coordinate system. Returns the side
  /// length of the axis-aligned box that bounds the transformed box.
  template <unsigned int D>
  inline Vector<double, D> operator()(const bbox<D> &b) const {
    return b.bmax - b.bmin;
  }
};

/// a linear transform (e.g. a skew coordinate transform)
template <unsigned int D, typename T> class LinearTransform {
  using double_d = Vector<double, D>;
  using int_d = Vector<int, D>;
  using bool_d = Vector<bool, D>;

  /// a user-defined transform that is assumed to be linear
  T m_point_transform;

  /// for any given axis-aligned bounding box, this is the vertex closest to the
  /// maximum eigen-vector (i.e. the vertex closest to the direction of maximum
  /// extension of the transformed box)
  bool_d m_eigen_vertices;

public:
  /// the transform stores the bounding box vertex that is maximally extended
  /// from the centre by the transform, this getter returns this vertex
  const bool_d &get_eigen_vertices() { return m_eigen_vertices; }

  /// takes a user-defined transform @p point_transform that is assumed to be
  /// linear. @p point_transform only has to defined an `operator()(const
  /// Vector<double,D>&)` function
  LinearTransform(const T &point_transform)
      : m_point_transform(point_transform) {

    // find closest vertex to max eigen vector
    bbox<D> unit(double_d::Constant(-1), double_d::Constant(1));
    double_d p;
    double max = 0;
    for (auto i = lattice_iterator<D>(int_d::Constant(0), int_d::Constant(2));
         i != false; ++i) {
      for (size_t j = 0; j < D; ++j) {
        p[j] = (*i)[j] == 1 ? unit.bmax[j] : unit.bmin[j];
      }
      const double pnorm2 = m_point_transform(p).squaredNorm();
      if (pnorm2 > max) {
        m_eigen_vertices = *i;
        max = pnorm2;
      }
    }
  }

  /// transforms the point @p the the new coordinate system, returning the
  /// result
  inline double_d operator()(const double_d &v) const {
    return m_point_transform(v);
  }

  /// transforms the box @p to the new coordinate system. Returns the side
  /// length of the axis-aligned box that bounds the transformed box.
  inline Vector<double, D> operator()(const bbox<D> &b) const {
    // store closest vertices to max eigen vector
    double_d max, min;
    for (size_t i = 0; i < D; ++i) {
      const double centre = 0.5 * (b.bmax[i] + b.bmin[i]);
      if (m_eigen_vertices[i]) {
        max[i] = b.bmax[i] - centre;
        min[i] = b.bmin[i] - centre;
      } else {
        max[i] = b.bmin[i] - centre;
        min[i] = b.bmax[i] - centre;
      }
    }

    // transform them
    max = m_point_transform(max);
    min = m_point_transform(min);

    // work out new bounds
    for (size_t i = 0; i < D; ++i) {
      const double tmp = std::max(max[i], min[i]);
      min[i] = std::min(max[i], min[i]);
      max[i] = tmp;
    }

    return max - min;
  }
};

/// a transform that scales the axis independently
template <unsigned int D> class ScaleTransform {
  using double_d = Vector<double, D>;
  using int_d = Vector<int, D>;
  using bool_d = Vector<bool, D>;

  double_d m_scale;

public:
  ScaleTransform() = default;
  ScaleTransform(const double_d &scale) : m_scale(scale) {}

  /// transforms the point @p the the new coordinate system, returning the
  /// result
  inline double_d operator()(const double_d &v) const { return v * m_scale; }

  /// transforms the box @p to the new coordinate system. Returns the side
  /// length of the axis-aligned box that bounds the transformed box.
  inline Vector<double, D> operator()(const bbox<D> &b) const {
    return (b.bmax - b.bmin) * m_scale;
  }
};

/// create a @ref LinearTransform
template <unsigned int D, typename T>
LinearTransform<D, T> create_linear_transform(const T &transform) {
  return LinearTransform<D, T>(transform);
}

/// create a @ref ScaleTransform
template <unsigned int D>
ScaleTransform<D> create_scale_transform(const Vector<double, D> &scale) {
  return ScaleTransform<D>(scale);
}

} // namespace Aboria

#endif
