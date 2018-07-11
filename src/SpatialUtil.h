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

  inline CUDA_HOST_DEVICE bbox(const double_d &p)
      : bmin(p), bmax(p + std::numeric_limits<double>::epsilon()) {}

  inline CUDA_HOST_DEVICE bbox(const double_d &min, const double_d &max)
      : bmin(min), bmax(max) {}

  inline CUDA_HOST_DEVICE bbox(const double *min, const double *max) {
    for (size_t i = 0; i < D; ++i) {
      bmin[i] = min[i];
      bmax[i] = max[i];
    }
  }

  ///
  /// @return translate the bounding box
  ///
  inline CUDA_HOST_DEVICE bbox &operator+=(const double_d &arg) {
    bmin += arg;
    bmax += arg;
    return *this;
  }

  ///
  /// @return translate the bounding box
  ///
  inline CUDA_HOST_DEVICE bbox operator+(const double_d &arg) const {
    bbox bounds = *this;
    bounds += arg;
    return bounds;
  }

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
      if (bmax[i] < bmin[i] + 3 * std::numeric_limits<double>::epsilon())
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

///
/// @brief Contains the middle and radius of a hypersphere in @p D
/// dimensional space
///
/// @tparam D the number of spatial dimensions
///
/// TODO: NOT IMPLEMENTED
template <unsigned int D> struct bsphere {
  typedef Vector<double, D> double_d;

  ///
  /// @brief middle point in the sphere
  ///
  double_d m_middle;

  ///
  /// @brief radius of the sphere
  ///
  double m_radius;

  inline CUDA_HOST_DEVICE bsphere()
      : m_middle(double_d::Constant(0)), m_radius(-1) {}

  inline CUDA_HOST_DEVICE bsphere(const double_d &p)
      : m_middle(p), m_radius(0) {}

  inline CUDA_HOST_DEVICE bsphere(const double_d &middle, const double &radius)
      : m_middle(middle), m_radius(radius) {}

  ///
  /// @return the sphere covering both input spheres
  ///
  inline CUDA_HOST_DEVICE bsphere operator+(const bsphere &arg) {
    // TODO: NOT IMPLEMENTED
    return arg;
  }

  ///
  /// @return true if sphere is empty
  ///
  inline CUDA_HOST_DEVICE bool is_empty() { return m_radius < 0; }
};

} // namespace Aboria

#endif
