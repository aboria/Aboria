
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

//
// Acknowledgement: This source was modified from the Thrust example
// bucket_sort2d.cu
//

#ifndef LATTICE_ITERATOR_H_
#define LATTICE_ITERATOR_H_

#include "Vector.h"

namespace Aboria {
template <unsigned int D> class lattice_iterator {
  typedef lattice_iterator<D> iterator;
  typedef Vector<double, D> double_d;
  typedef Vector<int, D> int_d;

  // make a proxy int_d in case you ever
  // want to get a pointer object to the
  // reference (which are both of the
  // same type)
  struct proxy_int_d : public int_d {
    CUDA_HOST_DEVICE
    proxy_int_d() : int_d() {}

    CUDA_HOST_DEVICE
    proxy_int_d(const int_d &arg) : int_d(arg) {}

    CUDA_HOST_DEVICE
    proxy_int_d &operator&() { return *this; }

    CUDA_HOST_DEVICE
    const proxy_int_d &operator&() const { return *this; }

    CUDA_HOST_DEVICE
    const proxy_int_d &operator*() const { return *this; }

    CUDA_HOST_DEVICE
    proxy_int_d &operator*() { return *this; }

    CUDA_HOST_DEVICE
    const proxy_int_d *operator->() const { return this; }

    CUDA_HOST_DEVICE
    proxy_int_d *operator->() { return this; }
  };

  int_d m_min;
  int_d m_max;
  proxy_int_d m_index;
  int_d m_size;
  bool m_valid;

public:
  typedef proxy_int_d pointer;
  typedef std::random_access_iterator_tag iterator_category;
  typedef const proxy_int_d &reference;
  typedef proxy_int_d value_type;
  typedef std::ptrdiff_t difference_type;

  CUDA_HOST_DEVICE
  lattice_iterator() : m_valid(false) {}

  CUDA_HOST_DEVICE
  lattice_iterator(const int_d &min, const int_d &max)
      : m_min(min), m_max(max), m_index(min), m_size(minus(max, min)),
        m_valid(true) {}

  CUDA_HOST_DEVICE
  lattice_iterator(const int_d &min, const int_d &max, const int_d &index)
      : m_min(min), m_max(max), m_index(index), m_size(minus(max, min)),
        m_valid(true) {}

  /*
  CUDA_HOST_DEVICE
  iterator& operator=(const iterator& copy) {
      m_index = copy.m_index;
      m_valid = copy.m_valid;
      ASSERT_CUDA(m_valid?(m_index >= m_min).all()&&(m_index <
  m_max).all():true); return *this;
  }
  */

  CUDA_HOST_DEVICE
  iterator &operator=(const int_d &copy) {
    m_index = copy;
    m_valid = (m_index >= m_min).all() && (m_index < m_max).all();
    return *this;
  }

  CUDA_HOST_DEVICE
  explicit operator size_t() const { return collapse_index_vector(m_index); }

  CUDA_HOST_DEVICE
  const lattice_iterator &get_child_iterator() const { return *this; }

  CUDA_HOST_DEVICE
  reference operator*() const { return dereference(); }

  CUDA_HOST_DEVICE
  reference operator->() const { return dereference(); }

  CUDA_HOST_DEVICE
  iterator &operator++() {
    increment();
    return *this;
  }

  CUDA_HOST_DEVICE
  iterator operator++(int) {
    iterator tmp(*this);
    operator++();
    return tmp;
  }

  CUDA_HOST_DEVICE
  iterator operator+(const int n) {
    iterator tmp(*this);
    tmp.increment(n);
    return tmp;
  }

  CUDA_HOST_DEVICE
  iterator &operator+=(const int n) {
    increment(n);
    return *this;
  }

  CUDA_HOST_DEVICE
  iterator &operator-=(const int n) {
    increment(-n);
    return *this;
  }

  CUDA_HOST_DEVICE
  iterator operator-(const int n) {
    iterator tmp(*this);
    tmp.increment(-n);
    return tmp;
  }

  CUDA_HOST_DEVICE
  size_t operator-(const iterator &start) const {
    int distance;
    if (!m_valid) {
      distance = start.collapse_index_vector(
                     minus(minus(start.m_max, 1), start.m_index)) +
                 1;
    } else if (!start.m_valid) {
      distance = collapse_index_vector(minus(m_index, m_min));
    } else {
      distance = collapse_index_vector(minus(m_index, start.m_index));
    }
    return distance;
  }

  CUDA_HOST_DEVICE
  inline bool operator==(const iterator &rhs) const { return equal(rhs); }

  CUDA_HOST_DEVICE
  inline bool operator==(const bool rhs) const { return equal(rhs); }

  CUDA_HOST_DEVICE
  inline bool operator!=(const iterator &rhs) const { return !operator==(rhs); }

  CUDA_HOST_DEVICE
  inline bool operator!=(const bool rhs) const { return !operator==(rhs); }

private:
  static inline CUDA_HOST_DEVICE int_d minus(const int_d &arg1,
                                             const int_d &arg2) {
    int_d ret;
    for (size_t i = 0; i < D; ++i) {
      ret[i] = arg1[i] - arg2[i];
    }
    return ret;
  }

  static inline CUDA_HOST_DEVICE int_d minus(const int_d &arg1,
                                             const int arg2) {
    int_d ret;
    for (size_t i = 0; i < D; ++i) {
      ret[i] = arg1[i] - arg2;
    }
    return ret;
  }

  CUDA_HOST_DEVICE
  int collapse_index_vector(const int_d &vindex) const {
    int index = 0;
    unsigned int multiplier = 1.0;
    for (int i = D - 1; i >= 0; --i) {
      if (i != D - 1) {
        multiplier *= m_size[i + 1];
      }
      index += multiplier * vindex[i];
    }
    return index;
  }

  CUDA_HOST_DEVICE
  int_d reassemble_index_vector(const int index) const {
    int_d vindex;
    int i = index;
    for (int d = D - 1; d >= 0; --d) {
      double div = (double)i / m_size[d];
      vindex[d] = std::round((div - std::floor(div)) * m_size[d]);
      i = std::floor(div);
    }
    return vindex;
  }

  CUDA_HOST_DEVICE
  bool equal(iterator const &other) const {
    if (!other.m_valid)
      return !m_valid;
    if (!m_valid)
      return !other.m_valid;
    for (size_t i = 0; i < D; ++i) {
      if (m_index[i] != other.m_index[i]) {
        return false;
      }
    }
    return true;
  }

  CUDA_HOST_DEVICE
  bool equal(const bool other) const { return m_valid == other; }

  CUDA_HOST_DEVICE
  reference dereference() const { return m_index; }

  CUDA_HOST_DEVICE
  void increment() {
    for (int i = D - 1; i >= 0; --i) {
      ++m_index[i];
      if (m_index[i] < m_max[i])
        break;
      if (i != 0) {
        m_index[i] = m_min[i];
      } else {
        m_valid = false;
      }
    }
  }

  CUDA_HOST_DEVICE
  void increment(const int n) {
    if (n == 1) {
      increment();
    } else {
      int collapsed_index = collapse_index_vector(m_index);
      m_index = reassemble_index_vector(collapsed_index += n);
    }
  }
};

} // namespace Aboria

#endif
