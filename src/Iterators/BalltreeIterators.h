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

#ifndef BALLTREE_ITERATORS_H_
#define BALLTREE_ITERATORS_H_

#include "Traits.h"

namespace Aboria {

template <typename Query, int LNormNumber> class balltree_query_iterator {
  typedef balltree_query_iterator<Query, LNormNumber> iterator;
  typedef typename Query::child_iterator child_iterator;
  static const unsigned int dimension = Query::dimension;
  typedef Vector<double, dimension> double_d;
  typedef Vector<int, dimension> int_d;
  typedef bbox<dimension> box_type;

public:
  typedef typename child_iterator::value_type value_type;
  typedef typename child_iterator::pointer pointer;
  typedef std::forward_iterator_tag iterator_category;
  typedef typename child_iterator::reference reference;
  typedef std::ptrdiff_t difference_type;

  CUDA_HOST_DEVICE
  balltree_query_iterator() {}

  /// this constructor is used to start the iterator at the head of a bucket
  /// list
  CUDA_HOST_DEVICE
  balltree_query_iterator(const child_iterator &start,
                          const double_d &query_point,
                          const double_d &max_distance,
                          const unsigned tree_depth, const Query *query,
                          const bool ordered = false)
      : m_query_point(query_point), m_inv_max_distance(1.0 / max_distance),
        m_query(query) {
    ASSERT_CUDA(tree_depth <= m_stack_max_size);
    if (start != false) {
      m_stack.push_back(start);
      go_to_next_leaf();
    } else {
#ifndef __CUDA_ARCH__
      LOG(3, "\tballtree_query_iterator (constructor) with query pt = "
                 << m_query_point
                 << "): start is false, no children to search.");
#endif
    }

    if (m_stack.empty()) {
#ifndef __CUDA_ARCH__
      LOG(3, "\tballtree_query_iterator (constructor) with query pt = "
                 << m_query_point
                 << "): search region outside domain or no children to "
                    "search.");
#endif
    } else {
#ifndef __CUDA_ARCH__
      LOG(3, "\tballtree_query_iterator (constructor) with query pt = "
                 << m_query_point
                 << "):  found bbox = " << m_query->get_bounds(m_stack.back()));
#endif
    }
  }

  CUDA_HOST_DEVICE
  balltree_query_iterator(const iterator &copy)
      : m_query_point(copy.m_query_point),
        m_inv_max_distance(copy.m_inv_max_distance), m_query(copy.m_query)

  {
    for (size_t i = 0; i < copy.m_stack.size(); ++i) {
      m_stack.push_back(copy.m_stack[i]);
    }
  }

  CUDA_HOST_DEVICE
  ~balltree_query_iterator() {}

  CUDA_HOST_DEVICE
  child_iterator get_child_iterator() const { return m_stack.back(); }

  CUDA_HOST_DEVICE
  iterator &operator=(const iterator &copy) {
    m_query_point = copy.m_query_point;
    m_inv_max_distance = copy.m_inv_max_distance;

    m_stack.resize(copy.m_stack.size());
    for (size_t i = 0; i < copy.m_stack.size(); ++i) {
      m_stack[i] = copy.m_stack[i];
    }

    m_query = copy.m_query;
    return *this;
  }

  /*
  iterator& operator=(const octtree_depth_first_iterator<Query>& copy) {
      m_stack = copy.m_stack;
      go_to_next_leaf();
      return *this;
  }
  */

  CUDA_HOST_DEVICE
  reference operator*() const { return dereference(); }

  CUDA_HOST_DEVICE
  reference operator->() { return dereference(); }

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
  size_t operator-(iterator start) const {
    size_t count = 0;
    while (start != *this) {
      start++;
      count++;
    }
    return count;
  }

  CUDA_HOST_DEVICE
  inline bool operator==(const iterator &rhs) const { return equal(rhs); }

  CUDA_HOST_DEVICE
  inline bool operator!=(const iterator &rhs) const { return !operator==(rhs); }

  CUDA_HOST_DEVICE
  inline bool operator==(const bool rhs) const { return equal(rhs); }

  CUDA_HOST_DEVICE
  inline bool operator!=(const bool rhs) const { return !operator==(rhs); }

private:
  friend class boost::iterator_core_access;

  CUDA_HOST_DEVICE
  bool child_is_within_query(const child_iterator &node) {
    const box_type &bounds = m_query->get_bounds(node);
    const double_d c = 0.5 * (bounds.bmax + bounds.bmin);
    const double r = (bounds.bmax - c).norm();
    const double_d dx = c - m_query_point;
    const double dx_mod = dx.norm();
    const double_d dx_hat = dx / dx_mod;
    const double closest_radius = dx_mod - r;
    if (closest_radius < 0) {
      return true;
    }
    const double_d dist = dx_hat * closest_radius;

    double accum = 0;
    for (size_t j = 0; j < dimension; j++) {
      accum = detail::distance_helper<LNormNumber>::accumulate_norm(
          accum, dist[j] * m_inv_max_distance[j]);
    }
    return (accum < 1.0);
  }

  CUDA_HOST_DEVICE
  void increment_stack() {
    while (!m_stack.empty()) {
      ++m_stack.back();
#ifndef __CUDA_ARCH__
      LOG(4,
          "\tincrement stack with child " << m_stack.back().get_child_number());
#endif
      if (m_stack.back() == false) {
#ifndef __CUDA_ARCH__
        LOG(4, "\tincrement_stack: pop");
#endif
        m_stack.pop_back();
      } else {
        break;
      }
    }
  }

  CUDA_HOST_DEVICE
  void go_to_next_leaf() {
    bool exit = m_stack.empty();
    while (!exit) {
      child_iterator &node = m_stack.back();
#ifndef __CUDA_ARCH__
      LOG(3, "\tgo_to_next_leaf with child " << node.get_child_number()
                                             << " with bounds "
                                             << m_query->get_bounds(node));
#endif
      if (child_is_within_query(node)) { // could be in this child
#ifndef __CUDA_ARCH__
        LOG(4, "\tthis child is within query, so going to next child");
#endif
        if (m_query->is_leaf_node(*node)) {
          exit = true;
        } else {
#ifndef __CUDA_ARCH__
          LOG(4, "\tdive down");
#endif
          m_stack.push_back(m_query->get_children(node));
        }
      } else { // not in this one, so go to next child, or go up if no more
               // children
#ifndef __CUDA_ARCH__
        LOG(4, "\tthis child is NOT within query, so going to next child");
#endif
        increment_stack();
        exit = m_stack.empty();
      }
    }
#ifndef __CUDA_ARCH__
    if (4 <= ABORIA_LOG_LEVEL) {
      if (m_stack.empty()) {
        LOG(4, "\tgo_to_next_leaf: stack empty, finishing");
      } else {
        LOG(4, "\tgo_to_next_leaf: found leaf, finishing");
      }
    }
#endif
  }

  CUDA_HOST_DEVICE
  void increment() {
#ifndef __CUDA_ARCH__
    LOG(4, "\tincrement (balltree_query_iterator):");
#endif
    increment_stack();
    go_to_next_leaf();

    if (m_stack.empty()) {
#ifndef __CUDA_ARCH__
      LOG(3, "\tend increment (balltree_query_iterator): no more nodes");
#endif
    } else {
#ifndef __CUDA_ARCH__
      LOG(3, "\tend increment (balltree_query_iterator): looking in bbox "
                 << m_query->get_bounds(m_stack.back()));
#endif
    }
  }

  CUDA_HOST_DEVICE
  bool equal(iterator const &other) const {
    if (m_stack.empty() || other.m_stack.empty()) {
      return m_stack.empty() == other.m_stack.empty();
    } else {
      return m_stack.back() == other.m_stack.back();
    }
  }

  CUDA_HOST_DEVICE
  bool equal(const bool other) const { return m_stack.empty() != other; }

  CUDA_HOST_DEVICE
  reference dereference() const { return *m_stack.back(); }

  // unsigned m_stack_size;
  const static unsigned m_stack_max_size = Query::m_max_tree_depth;
  static_vector<child_iterator, m_stack_max_size> m_stack;
  double_d m_query_point;
  double_d m_inv_max_distance;
  const Query *m_query;
};

template <typename Query, int LNormNumber> class balltree_box_query_iterator {
  typedef balltree_box_query_iterator<Query, LNormNumber> iterator;
  typedef typename Query::child_iterator child_iterator;
  static const unsigned int dimension = Query::dimension;
  typedef Vector<double, dimension> double_d;
  typedef Vector<int, dimension> int_d;
  typedef bbox<dimension> box_type;

public:
  typedef child_iterator const value_type;
  typedef child_iterator const *pointer;
  typedef std::forward_iterator_tag iterator_category;
  typedef child_iterator const &reference;
  typedef std::ptrdiff_t difference_type;

  CUDA_HOST_DEVICE
  balltree_box_query_iterator() {}

  /// this constructor is used to start the iterator at the head of a bucket
  /// list
  CUDA_HOST_DEVICE
  balltree_box_query_iterator(const child_iterator &start,
                              const box_type &query_box,
                              const double max_distance,
                              const unsigned tree_depth, const Query *query,
                              const bool ordered = false)
      : m_query_box(query_box), m_max_distance2(std::pow(max_distance, 2)),
        m_query(query) {
    ASSERT_CUDA(tree_depth <= m_stack_max_size);
    if (start != false) {
      m_stack.push_back(start);
      go_to_next_leaf();
    } else {
#ifndef __CUDA_ARCH__
      LOG(3, "\tballtree_box_query_iterator (constructor) with query box = "
                 << m_query_box << "): start is false, no children to search.");
#endif
    }

    if (m_stack.empty()) {
#ifndef __CUDA_ARCH__
      LOG(3, "\tballtree_box_query_iterator (constructor) with query box = "
                 << m_query_box
                 << "): search region outside domain or no children to "
                    "search.");
#endif
    } else {
#ifndef __CUDA_ARCH__
      LOG(3, "\tballtree_box_query_iterator (constructor) with query box = "
                 << m_query_box
                 << "):  found bbox = " << m_query->get_bounds(m_stack.back()));
#endif
    }
  }

  CUDA_HOST_DEVICE
  balltree_box_query_iterator(const iterator &copy)
      : m_query_box(copy.m_query_box), m_max_distance2(copy.m_max_distance2),
        m_query(copy.m_query)

  {
    for (size_t i = 0; i < copy.m_stack.size(); ++i) {
      m_stack.push_back(copy.m_stack[i]);
    }
  }

  CUDA_HOST_DEVICE
  ~balltree_box_query_iterator() {}

  CUDA_HOST_DEVICE
  child_iterator get_child_iterator() const { return m_stack.back(); }

  CUDA_HOST_DEVICE
  iterator &operator=(const iterator &copy) {
    m_query_box = copy.m_query_box;
    m_max_distance2 = copy.m_max_distance2;

    m_stack.resize(copy.m_stack.size());
    for (size_t i = 0; i < copy.m_stack.size(); ++i) {
      m_stack[i] = copy.m_stack[i];
    }

    m_query = copy.m_query;
    return *this;
  }

  /*
  iterator& operator=(const octtree_depth_first_iterator<Query>& copy) {
      m_stack = copy.m_stack;
      go_to_next_leaf();
      return *this;
  }
  */

  CUDA_HOST_DEVICE
  reference operator*() const { return dereference(); }

  CUDA_HOST_DEVICE
  reference operator->() { return dereference(); }

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
  size_t operator-(iterator start) const {
    size_t count = 0;
    while (start != *this) {
      start++;
      count++;
    }
    return count;
  }

  CUDA_HOST_DEVICE
  inline bool operator==(const iterator &rhs) const { return equal(rhs); }

  CUDA_HOST_DEVICE
  inline bool operator!=(const iterator &rhs) const { return !operator==(rhs); }

  CUDA_HOST_DEVICE
  inline bool operator==(const bool rhs) const { return equal(rhs); }

  CUDA_HOST_DEVICE
  inline bool operator!=(const bool rhs) const { return !operator==(rhs); }

private:
  CUDA_HOST_DEVICE
  void increment_stack() {
    while (!m_stack.empty()) {
      ++m_stack.back();
#ifndef __CUDA_ARCH__
      LOG(4,
          "\tincrement stack with child " << m_stack.back().get_child_number());
#endif
      if (m_stack.back() == false) {
#ifndef __CUDA_ARCH__
        LOG(4, "\tincrement_stack: pop");
#endif
        m_stack.pop_back();
      } else {
        break;
      }
    }
  }

  CUDA_HOST_DEVICE
  void go_to_next_leaf() {
    bool exit = m_stack.empty();
    while (!exit) {
      child_iterator &node = m_stack.back();
#ifndef __CUDA_ARCH__
      LOG(3, "\tgo_to_next_leaf with child " << node.get_child_number()
                                             << " with bounds "
                                             << m_query->get_bounds(node));
#endif

      if (detail::spheres_within_distance<LNormNumber>(
              m_query_box, m_query->get_bounds(node),
              m_max_distance2)) { // could be in this child

#ifndef __CUDA_ARCH__
        LOG(4, "\tthis child is within query, so going to next child");
#endif
        if (m_query->is_leaf_node(*node)) {
          exit = true;
        } else {
#ifndef __CUDA_ARCH__
          LOG(4, "\tdive down");
#endif
          m_stack.push_back(m_query->get_children(node));
        }
      } else { // not in this one, so go to next child, or go up if no more
               // children
#ifndef __CUDA_ARCH__
        LOG(4, "\tthis child is NOT within query, so going to next child");
#endif
        increment_stack();
        exit = m_stack.empty();
      }
    }
#ifndef __CUDA_ARCH__
    if (4 <= ABORIA_LOG_LEVEL) {
      if (m_stack.empty()) {
        LOG(4, "\tgo_to_next_leaf: stack empty, finishing");
      } else {
        LOG(4, "\tgo_to_next_leaf: found leaf, finishing");
      }
    }
#endif
  }

  CUDA_HOST_DEVICE
  void increment() {
#ifndef __CUDA_ARCH__
    LOG(4, "\tincrement (balltree_box_query_iterator):");
#endif
    increment_stack();
    go_to_next_leaf();

    if (m_stack.empty()) {
#ifndef __CUDA_ARCH__
      LOG(3, "\tend increment (balltree_box_query_iterator): no more nodes");
#endif
    } else {
#ifndef __CUDA_ARCH__
      LOG(3, "\tend increment (balltree_box_query_iterator): looking in bbox "
                 << m_query->get_bounds(m_stack.back()));
#endif
    }
  }

  CUDA_HOST_DEVICE
  bool equal(iterator const &other) const {
    if (m_stack.empty() || other.m_stack.empty()) {
      return m_stack.empty() == other.m_stack.empty();
    } else {
      return m_stack.back() == other.m_stack.back();
    }
  }

  CUDA_HOST_DEVICE
  bool equal(const bool other) const { return m_stack.empty() != other; }

  CUDA_HOST_DEVICE
  reference dereference() const { return m_stack.back(); }

  // unsigned m_stack_size;
  const static unsigned m_stack_max_size = Query::m_max_tree_depth;
  static_vector<child_iterator, m_stack_max_size> m_stack;
  box_type m_query_box;
  double m_max_distance2;
  const Query *m_query;
};

} // namespace Aboria

#endif
