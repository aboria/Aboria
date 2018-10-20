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

#ifndef DISTANCE_H_
#define DISTANCE_H_

#include "Vector.h"
#include <algorithm>
#include <cmath>

namespace Aboria {
namespace detail {

template <int LNormNumber> struct distance_helper {
  CUDA_HOST_DEVICE
  static inline double get_value_to_accumulate(const double arg) {
    switch (LNormNumber) {
    case -1:
      return std::abs(arg);
    case 0:
      return arg != 0;
    case 1:
      return std::abs(arg);
    case 2:
      return std::pow(arg, LNormNumber);
    case 3:
      return std::abs(std::pow(arg, LNormNumber));
    case 4:
      return std::pow(arg, LNormNumber);
    default:
      return std::abs(std::pow(arg, LNormNumber));
    }
  }

  template <unsigned int D, typename VectorType = Vector<double, D>>
  CUDA_HOST_DEVICE static inline VectorType
  get_value_to_accumulate(const VectorType &arg) {
    VectorType ret;
    switch (LNormNumber) {
    case -1:
      for (size_t i = 0; i < D; ++i) {
        ret[i] = std::abs(arg[i]);
      }
    case 0:
      for (size_t i = 0; i < D; ++i) {
        ret[i] = arg[i] != 0;
      }
    case 1:
      for (size_t i = 0; i < D; ++i) {
        ret[i] = std::abs(arg[i]);
      }
    case 2:
      for (size_t i = 0; i < D; ++i) {
        ret[i] = std::pow(arg[i], LNormNumber);
      }
    case 3:
      for (size_t i = 0; i < D; ++i) {
        ret[i] = std::abs(std::pow(arg[i], LNormNumber));
      }
    case 4:
      for (size_t i = 0; i < D; ++i) {
        ret[i] = std::pow(arg[i], LNormNumber);
      }
    default:
      for (size_t i = 0; i < D; ++i) {
        ret[i] = std::abs(std::pow(arg[i], LNormNumber));
      }
    }
    return ret;
  }

  CUDA_HOST_DEVICE
  static inline double do_accumulate(const double accum, const double value) {
    switch (LNormNumber) {
    case -1:
      if (value > accum) {
        return value;
      } else {
        return accum;
      }
    default:
      return accum + value;
    }
  }

  CUDA_HOST_DEVICE
  static inline double accumulate_norm(const double accum, const double arg) {
    return do_accumulate(accum, get_value_to_accumulate(arg));
  }

  CUDA_HOST_DEVICE
  static inline double accumulate_max_norm(const double accum,
                                           const double arg1,
                                           const double arg2) {
    return do_accumulate(accum, std::max(get_value_to_accumulate(arg1),
                                         get_value_to_accumulate(arg2)));
  }

  template <unsigned int D>
  CUDA_HOST_DEVICE static inline double norm2(const Vector<double, D> &vector) {
    double accum = 0;
    for (size_t i = 0; i < D; ++i) {
      accum = accumulate_norm(accum, vector[i]);
    }
    return accum;
  }
};

} // namespace detail
} // namespace Aboria

#endif
