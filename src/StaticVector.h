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

#ifndef STATIC_VECTOR_H_
#define STATIC_VECTOR_H_

#include "CudaInclude.h"
#include <type_traits>

namespace Aboria {

template <class T, std::size_t N> class static_vector {
  // properly aligned uninitialized storage for N T's
  typename std::aligned_storage<sizeof(T), alignof(T)>::type data[N];
  size_t m_size = 0;

public:
  CUDA_HOST_DEVICE
  void push_back(const T &val) {
    ASSERT_CUDA(m_size < N);
    new (data + m_size) T(val);
    ++m_size;
  }

  CUDA_HOST_DEVICE
  void push_back(T &&val) {
    ASSERT_CUDA(m_size < N);
    new (data + m_size) T(val);
    ++m_size;
  }

  // Create an object in aligned storage
  template <typename... Args>
  CUDA_HOST_DEVICE void emplace_back(Args &&... args) {
    ASSERT_CUDA(m_size < N);
    new (data + m_size) T(std::forward<Args>(args)...);
    ++m_size;
  }

  CUDA_HOST_DEVICE
  void resize(size_t n) { m_size = n; }

  CUDA_HOST_DEVICE
  size_t size() const { return m_size; }

  CUDA_HOST_DEVICE
  bool empty() const { return m_size == 0; }

  CUDA_HOST_DEVICE
  void pop_back() {
    --m_size;
    reinterpret_cast<T *>(data + m_size)->~T();
  }

  CUDA_HOST_DEVICE
  T &back() { return *reinterpret_cast<T *>(data + m_size - 1); }

  CUDA_HOST_DEVICE
  const T &back() const {
    return *reinterpret_cast<const T *>(data + m_size - 1);
  }

  // Access an object in aligned storage
  CUDA_HOST_DEVICE
  const T &operator[](std::size_t pos) const {
    return *reinterpret_cast<const T *>(data + pos);
  }

  CUDA_HOST_DEVICE
  T &operator[](std::size_t pos) { return *reinterpret_cast<T *>(data + pos); }

  // Delete objects from aligned storage
  CUDA_HOST_DEVICE
  ~static_vector() {
    for (std::size_t pos = 0; pos < m_size; ++pos) {
      reinterpret_cast<T *>(data + pos)->~T();
    }
  }
};
} // namespace Aboria

#endif
