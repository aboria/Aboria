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

#ifndef DETAIL_PRECONDITIONERS_H_
#define DETAIL_PRECONDITIONERS_H_

#include <algorithm>
#include <chrono>
#include <fstream>
#include <unordered_map>

#ifdef HAVE_CAIRO
#include <cairo-svg.h>
#endif

#ifdef HAVE_EIGEN
#include <unsupported/Eigen/SparseExtra>

#include "Operators.h"

namespace Aboria {

namespace detail {
template <typename Function, typename Dest, unsigned int NI, unsigned int NJ,
          typename Blocks, typename Rhs>
void apply_function_to_diagonal_blocks(
    Function &&f, Dest &y, const MatrixReplacement<NI, NJ, Blocks> &mat,
    const Rhs &rhs, std::integral_constant<unsigned int, NI>) {}

template <typename Function, unsigned int NI, unsigned int NJ, typename Blocks>
void apply_function_to_diagonal_blocks(
    Function &&f, const MatrixReplacement<NI, NJ, Blocks> &mat,
    std::integral_constant<unsigned int, NI>) {}

template <typename Function, typename Dest, unsigned int NI, unsigned int NJ,
          typename Blocks, typename Rhs, unsigned int I>
void apply_function_to_diagonal_blocks(
    Function &&f, Dest &x, const MatrixReplacement<NI, NJ, Blocks> &mat,
    const Rhs &b, std::integral_constant<unsigned int, I>) {
  f(x.segment(mat.template start_row<I>(), mat.template size_row<I>()),
    b.segment(mat.template start_col<I>(), mat.template size_col<I>()),
    std::get<I * NJ + I>(mat.m_blocks));
  apply_function_to_diagonal_blocks(
      std::forward<Function>(f), x, mat, b,
      std::integral_constant<unsigned int, I + 1>());
}

template <typename Function, unsigned int NI, unsigned int NJ, typename Blocks,
          unsigned int I>
void apply_function_to_diagonal_blocks(
    Function &&f, const MatrixReplacement<NI, NJ, Blocks> &mat,
    std::integral_constant<unsigned int, I>) {
  f(std::get<I * NJ + I>(mat.m_blocks));
  apply_function_to_diagonal_blocks(
      std::forward<Function>(f), mat,
      std::integral_constant<unsigned int, I + 1>());
}

template <typename Function, typename Dest, unsigned int NI, unsigned int NJ,
          typename Blocks, typename Rhs>
void apply_function_to_diagonal_blocks(
    Function &&function, Dest &x, const MatrixReplacement<NI, NJ, Blocks> &mat,
    const Rhs &b) {
  apply_function_to_diagonal_blocks(std::forward<Function>(function), x, mat, b,
                                    std::integral_constant<unsigned int, 0>());
}

template <typename Function, unsigned int NI, unsigned int NJ, typename Blocks>
void apply_function_to_diagonal_blocks(
    Function &&function, const MatrixReplacement<NI, NJ, Blocks> &mat) {
  apply_function_to_diagonal_blocks(std::forward<Function>(function), mat,
                                    std::integral_constant<unsigned int, 0>());
}

template <typename T> struct storage_vector_type {
  T *m_data;
  size_t m_size;

  storage_vector_type() : m_data(nullptr), m_size(0) {}
  storage_vector_type(const storage_vector_type &other)
      : m_data(nullptr), m_size(0) {
    resize(other.m_size);
    for (size_t i = 0; i < m_size; ++i) {
      m_data[i] = other.m_data[i];
    }
  }

  storage_vector_type(storage_vector_type &&other) = default;

  ~storage_vector_type() {
    if (m_data)
      delete[] m_data;
  }

  storage_vector_type &operator=(const storage_vector_type &other) {
    resize(other.m_size);
    for (size_t i = 0; i < m_size; ++i) {
      m_data[i] = other.m_data[i];
    }
    return *this;
  }
  CUDA_HOST_DEVICE
  T *begin() const { return m_data; }
  CUDA_HOST_DEVICE
  T *end() const { return m_data + m_size; }
  CUDA_HOST_DEVICE
  size_t size() const { return m_size; }
  CUDA_HOST_DEVICE
  void resize(const size_t n) {
    // only supports reductions in size after initial resize
    ASSERT_CUDA(!m_data || n <= m_size);
    if (m_data) {
      delete[] m_data;
    }
    if (n > 0) {
      m_data = new T[n];
    }
    m_size = n;
  }
  CUDA_HOST_DEVICE
  T &operator[](const int i) { return m_data[i]; }
  CUDA_HOST_DEVICE
  const T &operator[](const int i) const { return m_data[i]; }
};
} // namespace detail

} // namespace Aboria

#endif
#endif
