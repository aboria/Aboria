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

#ifndef SEARCH_H_
#define SEARCH_H_

#include "CudaInclude.h"
#include "Get.h"
#include "NeighbourSearchBase.h"
#include "SpatialUtil.h"
#include "Traits.h"
#include "Vector.h"
#include "detail/Algorithms.h"
#include "detail/Distance.h"

#include "Log.h"
#include <cmath>
#include <iostream>
#include <queue>

namespace Aboria {

template <typename QueryA, typename QueryB, int LNormNumber>
class search_bf_iterator {
  typedef search_bf_iterator<QueryA, QueryB, LNormNumber> iterator;
  typedef typename QueryA::child_iterator child_iterator_a;
  typedef typename QueryB::child_iterator child_iterator_b;
  static const unsigned int dimension = QueryA::dimension;
  static_assert(QueryB::dimension == QueryA::dimension,
                "query dimensions must match");
  typedef Vector<double, dimension> double_d;
  typedef Vector<int, dimension> int_d;
  typedef typename QueryA::traits_type traits_type_a;
  typedef typename QueryB::traits_type traits_type_b;
  typedef
      typename traits_type_a::template tuple<child_iterator_a, child_iterator_b>
          ci_pair;
  typedef typename traits_type_a::template vector<ci_pair> vector_ci_pair;
  typedef typename traits_type_a::template vector<vint3> vector_vint3;
  typedef typename traits_type_a::template vector<int> vector_int;

public:
  typedef vector_ci_pair const value_type;
  typedef vector_ci_pair const *pointer;
  typedef vector_ci_pair const &reference;
  typedef std::forward_iterator_tag iterator_category;
  typedef std::ptrdiff_t difference_type;

  CUDA_HOST_DEVICE
  search_bf_iterator() {}

  /// this constructor is used to start the iterator at the head of a bucket
  /// list
  CUDA_HOST_DEVICE
  search_bf_iterator(const child_iterator_a &start_node_a,
                     const child_iterator_b &start_node_b,
                     const QueryA &query_a, const QueryB &query_b,
                     const double distance, const bool save_internal_leafs)
      : m_distance2(
            detail::distance_helper<LNormNumber>::get_value_to_accumulate(
                distance)),
        m_level_num(1), m_save_internal_leafs(save_internal_leafs),
        m_query_a(&query_a), m_query_b(&query_b)

  {
    if (start_node_a != false && start_node_b != false) {
      if (m_save_internal_leafs) {
        // put all n^2 ci pairs in m_level
        const int size_a = start_node_a.distance_to_end();
        const int size_b = start_node_b.distance_to_end();
        m_level.resize(size_a * size_b);
        detail::tabulate(m_level.begin(), m_level.end(),
                         add_all_buckets{start_node_a, start_node_b, size_b});
      } else {
        // only put in neighbouring buckets
        const int size_a = start_node_a.distance_to_end();
        vector_int vcount(size_a);
        auto count = QueryA::traits_type::make_counting_iterator(0);
        count_neighbouring_buckets fcount{start_node_a, *m_query_a, *m_query_b,
                                          distance};
        detail::transform_exclusive_scan(count, count + size_a, vcount.begin(),
                                         fcount, 0, detail::plus());
        const int n_neighbours = vcount.back() + fcount(size_a);
        m_level.resize(n_neighbours);
        detail::for_each(
            count, count + size_a,
            add_neighbouring_buckets{start_node_a,
                                     iterator_to_raw_pointer(m_level.begin()),
                                     iterator_to_raw_pointer(vcount.begin()),
                                     *m_query_a, *m_query_b, distance});
      }
    } else {
#ifndef __CUDA_ARCH__
      LOG(3, "\tsearch_bf_iterator (constructor): start is false, no "
             "children to search.");
#endif
    }
  }

  CUDA_HOST_DEVICE
  reference leafs() const { return m_leafs; }

  CUDA_HOST_DEVICE
  reference internal_leafs() const { return m_internal_leafs; }

  CUDA_HOST_DEVICE
  reference operator*() const { return dereference(); }

  CUDA_HOST_DEVICE
  pointer operator->() { return &dereference(); }

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

  struct add_all_buckets {
    child_iterator_a first_ci_a;
    child_iterator_b first_ci_b;
    int m_size_b;

    CUDA_HOST_DEVICE
    ci_pair operator()(const int ij) {
      // ij = i*size_b + j;
      const int j = ij % m_size_b;
      const int i = ij / m_size_b;
      auto ci_a = first_ci_a + i;
      auto ci_b = first_ci_b + j;
      return {ci_a, ci_b};
    }
  };

  struct count_neighbouring_buckets {
    child_iterator_a first_ci_a;
    const QueryA m_query_a;
    const QueryB m_query_b;
    const double m_distance;

    CUDA_HOST_DEVICE
    int operator()(const int i) {
      int ret = 0;
      auto ci_a = first_ci_a + i;
      for (auto ci_b = m_query_b.get_buckets_near_bucket(
               m_query_a.get_bounds(ci_a), m_distance);
           ci_b != false; ++ci_b, ++ret)
        ;
      return ret;
    }
  };

  struct add_neighbouring_buckets {
    child_iterator_a first_ci_a;
    ci_pair *m_level;
    int *m_count;
    const QueryA m_query_a;
    const QueryB m_query_b;
    const double m_distance;

    CUDA_HOST_DEVICE void operator()(const int i) {
      int j = 0;
      auto ci_a = first_ci_a + i;
      for (auto ci_b = m_query_b.get_buckets_near_bucket(
               m_query_a.get_bounds(ci_a), m_distance);
           ci_b != false; ++ci_b, ++j) {
        detail::get_impl<0>(m_level[m_count[i] + j]) = *ci_a;
        detail::get_impl<1>(m_level[m_count[i] + j]) = *ci_b;
      }
    }
  };

  struct count_children {
    const QueryA m_query_a;
    const QueryB m_query_b;
    const double m_distance2;

    CUDA_HOST_DEVICE
    auto operator()(const ci_pair &ci2) {
      // op 1: a and b are greater than search, nchild = 0
      //  if not:
      // op 2: a is parent, b is parent, nchild = nca*ncb
      // op 3: a is parent, b is leaf, nchild = nca*1
      // op 4: a is leaf, b is parent, nchild = 1*ncb
      // op 5: a is leaf, b is leaf, nchild = 0
      auto &a = detail::get_impl<0>(ci2);
      auto &b = detail::get_impl<1>(ci2);
      // btw, nchild = 1, really means nchild = 0
      int nchild;
      if (boxes_within_distance(m_query_a.get_bounds(a),
                                m_query_b.get_bounds(b), m_distance2)) {
        const int nchild_a = m_query_a.is_leaf_node(*a)
                                 ? 1
                                 : m_query_a.get_children(a).distance_to_end();
        const int nchild_b = m_query_b.is_leaf_node(*b)
                                 ? 1
                                 : m_query_b.get_children(b).distance_to_end();
        nchild = nchild_a * nchild_b;
      } else {
        nchild = 0;
      }

      return vint3(nchild, nchild == 1, nchild == 0);
    }
  };

  struct copy_children_and_leafs {
    size_t m_end;
    double m_distance2;
    vint3 *m_counts;
    ci_pair *m_leafs;
    ci_pair *m_internal_leafs;
    size_t m_leafs_old_size;
    const QueryA m_query_a;
    const QueryB m_query_b;
    ci_pair *m_next_level;
    ci_pair *m_level;
    const bool m_save_internal_leafs;

    CUDA_HOST_DEVICE
    void operator()(const int i) {
      vint3 my_count;
      if (i == static_cast<int>(m_end)) {
        count_children count{m_query_a, m_query_b, m_distance2};
        my_count = count(m_level[i]);
      } else {
        my_count = m_counts[i + 1] - m_counts[i];
      }

      if (my_count[2]) {
        // pair outside distance
        if (m_save_internal_leafs) {
          const int internal_leafs_index = m_counts[i][2];
          m_leafs[internal_leafs_index] = m_level[i];
        }
      } else if (my_count[1]) {
        // pair are both leafs
        const int leafs_index = m_leafs_old_size + m_counts[i][1];
        m_leafs[leafs_index] = m_level[i];
      } else {
        // pair have children
        auto ci_a = detail::get_impl<0>(m_level[i]);
        auto ci_b = detail::get_impl<1>(m_level[i]);
        {
          int next_level_index = m_counts[i][0];
          for (auto child_a = m_query_a.get_children(ci_a); child_a != false;
               ++child_a) {
            for (auto child_b = m_query_b.get_children(ci_b); child_b != false;
                 ++child_b, ++next_level_index) {
              m_next_level[next_level_index] = ci_pair(child_a, child_b);
            }
          }
        }
      }
    }
  };

private:
  CUDA_HOST_DEVICE static bool boxes_within_distance(const bbox<dimension> &a,
                                                     const bbox<dimension> &b,
                                                     const double distance2) {
    double squaredNorm = 0;
    for (size_t i = 0; i < dimension; ++i) {
      if (a.bmax[i] < b.bmin[i]) {
        detail::distance_helper<LNormNumber>::do_accumulate(
            squaredNorm,
            detail::distance_helper<LNormNumber>::get_value_to_accumulate(
                b.bmin[i] - a.bmax[i]));
      } else if (b.bmax[i] < a.bmin[i]) {
        detail::distance_helper<LNormNumber>::do_accumulate(
            squaredNorm,
            detail::distance_helper<LNormNumber>::get_value_to_accumulate(
                a.bmin[i] - b.bmax[i]));
      }
    }
    return squaredNorm <= distance2;
  }

  CUDA_HOST_DEVICE
  void increment() {
#ifndef __CUDA_ARCH__
    LOG(4, "\tincrement (search_bf_iterator): m_level size = "
               << m_level.size() << " m_leafs size = " << m_leafs.size());
#endif
    // m_level [ci0, ci1, ci2, ...]
    // exclusive scan for # children + # leafs: n_child
    // (std::vector<Vector<int,2>> = [{0,0}, {3,0}, {3,1}, ..., {N-#child
    // cin,NL-#child==0}] resize m_next_level(N) resize m_leafs(NL)

    m_counts.resize(m_level.size());
    count_children fcount_children{*m_query_a, *m_query_b, m_distance2};
    detail::transform_exclusive_scan(m_level.begin(), m_level.end(),
                                     m_counts.begin(), fcount_children,
                                     vint3::Constant(0), detail::plus());

    // resize for new children and leafs
    const vint3 nchildren =
        static_cast<vint3>(m_counts.back()) + fcount_children(m_level.back());
    m_next_level.resize(nchildren[0]);
    const size_t m_leafs_old_size = m_leafs.size();
    m_leafs.resize(m_leafs.size() + nchildren[1]);
    if (m_save_internal_leafs) {
      m_internal_leafs.resize(nchildren[2]);
    }

    // tabulate m_level to copy children to m_next_level, or leafs to
    // m_leafs
    auto count = QueryA::traits_type::make_counting_iterator(0);
    detail::for_each(count, count + m_level.size(),
                     copy_children_and_leafs{
                         m_level.size(), m_distance2,
                         iterator_to_raw_pointer(m_counts.begin()),
                         iterator_to_raw_pointer(m_leafs.begin()),
                         iterator_to_raw_pointer(m_internal_leafs.begin()),
                         m_leafs_old_size, *m_query_a, *m_query_b,
                         iterator_to_raw_pointer(m_next_level.begin()),
                         iterator_to_raw_pointer(m_level.begin()),
                         m_save_internal_leafs});

    // swap level back to m_level and increment level count
    m_level.swap(m_next_level);
    m_level_num++;
#ifndef __CUDA_ARCH__
    LOG(4, "\tend increment (search_bf_iterator):");
#endif
  }

  CUDA_HOST_DEVICE
  bool equal(iterator const &other) const {
    return m_query_a == other.m_query_b && m_level_num == other.m_level_num;
  }

  CUDA_HOST_DEVICE bool equal(const bool other) const {
    return m_level.empty();
  }

  CUDA_HOST_DEVICE
  reference dereference() const { return m_level; }

  double m_distance2;
  vector_ci_pair m_level;
  vector_ci_pair m_next_level;
  vector_vint3 m_counts;
  vector_ci_pair m_leafs;
  vector_ci_pair m_internal_leafs;
  size_t m_level_num;
  const bool m_save_internal_leafs;
  const QueryA *m_query_a;
  const QueryB *m_query_b;
}; // namespace Aboria

/// A const iterator to a set of neighbouring points. This iterator
/// implements a STL forward iterator type
// assume that these iterators, and query functions, are only called from
// device code
template <typename Query, int LNormNumber> class search_iterator {

  typedef typename Query::particle_iterator particle_iterator;
  typedef typename Query::template query_iterator<LNormNumber> query_iterator;
  typedef typename Query::traits_type Traits;
  static const unsigned int dimension = Traits::dimension;

  typedef typename Traits::position position;
  typedef typename Traits::double_d double_d;
  typedef typename Traits::bool_d bool_d;
  typedef typename Traits::int_d int_d;
  typedef typename particle_iterator::value_type p_value_type;
  typedef typename particle_iterator::reference p_reference;
  typedef typename particle_iterator::pointer p_pointer;
  typedef lattice_iterator<dimension> periodic_iterator_type;

  ///
  /// @brief true if the iterator has run out of particles within the
  /// search distance
  ///
  bool m_valid;

  ///
  /// @brief the central point for the distance search
  ///
  double_d m_r;

  ///
  /// @brief a vector to store the distance between m_r and the current
  /// candidate point
  ///
  double_d m_dx;

  ///
  /// @brief pointer to the query object for the spatial data structure
  ///
  const Query *m_query;

  ///
  /// @brief the search distance
  ///
  double m_max_distance;

  ///
  /// @brief the search distance squared
  ///
  double m_max_distance2;

  ///
  /// @brief iterator for searching periodic domains, iterates over the
  /// periodic lattice
  ///
  periodic_iterator_type m_current_periodic;

  ///
  /// @brief the current central search position, taking into account
  /// periodicity
  ///
  double_d m_current_point;

  ///
  /// @brief the current candidate bucket
  ///
  query_iterator m_current_bucket;

  ///
  /// @brief current iterator into candidate bucket
  ///
  particle_iterator m_current_particle;

public:
  typedef p_pointer pointer;
  typedef std::forward_iterator_tag iterator_category;
  typedef p_reference reference;
  typedef p_value_type value_type;
  typedef std::ptrdiff_t difference_type;

  ///
  /// @brief returns iterator for periodic lattice, given the periodicity
  /// of the domain
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  static periodic_iterator_type get_periodic_range(const bool_d is_periodic) {
    int_d start, end;
    for (size_t i = 0; i < dimension; ++i) {
      start[i] = is_periodic[i] ? -1 : 0;
      end[i] = is_periodic[i] ? 2 : 1;
    }
    return periodic_iterator_type(start, end);
  }

  ///
  /// @brief constructs an invalid iterator that can be used as an end()
  /// iterator
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  search_iterator() : m_valid(false) {}

  ///
  /// @brief should generally use this constructor to make a search
  /// iterator. Returns an iterator that will search around the given
  /// point, and iterate through all the particles it finds within the
  /// given maximum distance
  ///
  /// @param query a query object for a spatial data structure
  /// @param r the central point to search around
  /// @param max_distance the maximum distance to search
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  search_iterator(const Query &query, const double_d &r,
                  const double max_distance)
      : m_valid(true), m_r(r), m_query(&query), m_max_distance(max_distance),
        m_max_distance2(
            detail::distance_helper<LNormNumber>::get_value_to_accumulate(
                max_distance)),
        m_current_periodic(get_periodic_range(m_query->get_periodic())),
        m_current_point(
            r + (*m_current_periodic) *
                    (m_query->get_bounds().bmax - m_query->get_bounds().bmin)),
        m_current_bucket(query.template get_buckets_near_point<LNormNumber>(
            m_current_point, max_distance)) {

#if defined(__CUDA_ARCH__)
    CHECK_CUDA((!std::is_same<typename Traits::template vector<double>,
                              std::vector<double>>::value),
               "Cannot use std::vector in device code");

    LOG_CUDA(3, "\tconstructor (search_iterator)");
#else
    LOG(3, "\tconstructor (search_iterator with query pt = "
               << m_r << ", and m_current_point = " << m_current_point << ")");
#endif
    if ((m_valid = get_valid_bucket())) {
      m_current_particle = m_query->get_bucket_particles(*m_current_bucket);
      if ((m_valid = get_valid_candidate())) {
        if (!check_candidate()) {
          increment();
        }
      }
    }
#if defined(__CUDA_ARCH__)
    if (m_valid) {
      LOG_CUDA(3, "\tconstructor (search_iterator) found good candidate");
    } else {
      LOG_CUDA(3, "\tconstructor (search_iterator) didn't find good candidate");
    }
#else
    if (m_valid) {
      LOG_BOLD(3, "\tconstructor (search_iterator with query pt = "
                      << m_r << "): found good canditate at "
                      << get<position>(*m_current_particle));
    } else {
      LOG(3, "\tconstructor (search_iterator with query pt = "
                 << m_r << "): didn't find good candidate");
    }

#endif
  }

  ///
  /// @brief returns the distance $r_b-r_a$ between the current candidate
  /// position $r_b$ an the central search point $r_a$
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  const double_d &dx() const { return m_dx; }

  ///
  /// @brief dereference the iterator, returns a reference to the current
  /// candidate particle
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  reference operator*() const { return dereference(); }

  ///
  /// @brief dereference the iterator, returns a reference to the current
  /// candidate particle
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  reference operator->() { return dereference(); }

  ///
  /// @brief increment the iterator, i.e. move to the next candidate
  /// particle
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  search_iterator &operator++() {
    increment();
    return *this;
  }

  ///
  /// @brief post increment the iterator, i.e. move to the next candidate
  /// particle and return the original iterator
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  search_iterator operator++(int) {
    search_iterator tmp(*this);
    operator++();
    return tmp;
  }

  ///
  /// @brief returns the distance between start and this
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  size_t operator-(search_iterator start) const {
    size_t count = 0;
    while (start != *this) {
      start++;
      count++;
    }
    return count;
  }

  ///
  /// @brief returns how many candidate particle to go until the iterator
  /// becomes invalid
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  size_t distance_to_end() const { return search_iterator() - *this; }

  ///
  /// @brief returns true if @p rhs is the same as this
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  inline bool operator==(const search_iterator &rhs) { return equal(rhs); }

  ///
  /// @brief returns true if @p rhs is not the same as this
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  inline bool operator!=(const search_iterator &rhs) {
    return !operator==(rhs);
  }

  ///
  /// @brief the search iterator can be converted to true if it is
  /// pointing to a valid candidate point, false otherwise
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  inline bool operator==(const bool rhs) const { return equal(rhs); }

  ///
  /// @brief the search iterator can be converted to true if it is
  /// pointing to a valid candidate point, false otherwise
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  inline bool operator!=(const bool rhs) const { return !operator==(rhs); }

private:
  ///
  /// @brief if both iterators are valid, and pointing to the same
  /// particle, then they are equal
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  bool equal(search_iterator const &other) const {
    return m_valid ? m_current_particle == other.m_current_particle
                   : !other.m_valid;
  }

  ///
  /// @brief if this->m_valid is true then the iterator is equal to true
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  bool equal(const bool other) const { return m_valid == other; }

  ///
  /// @brief to be called after incrementing m_current_bucket. If
  /// m_current bucket is no longer valid, then move to the next periodic
  /// lattice. If not more periodic lattices to search through, then the
  /// iterator becomes invalid
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  bool get_valid_bucket() {
    LOG_CUDA(4, "\tget_valid_bucket:");

    while (m_current_bucket == false) {
#ifdef __CUDA_ARCH__
      if (3 <= ABORIA_LOG_LEVEL) {
        printf("\t\tgo_to_next periodic (search_iterator): "
               "m_current_periodic "
               "= (");
        for (int i = 0; i < Traits::dimension; ++i) {
          printf("%d,", (*m_current_periodic)[i]);
        }
        printf("\n,");
      }
#else
      LOG(3, "\tgo_to_next periodic (search_iterator): m_current_periodic = "
                 << *m_current_periodic);
#endif
      ++m_current_periodic;
      if (m_current_periodic == false) {
#ifdef __CUDA_ARCH__
        LOG_CUDA(4, "\tran out of buckets to search (search_iterator):");
#else
        LOG(4, "\tran out of buckets to search (search_iterator):");
#endif
        return false;
      }
      m_current_point =
          m_r + (*m_current_periodic) *
                    (m_query->get_bounds().bmax - m_query->get_bounds().bmin);
      m_current_bucket = m_query->template get_buckets_near_point<LNormNumber>(
          m_current_point, m_max_distance);
    }
    return true;
  }

  ///
  /// @brief to be called after incrementing m_current_particle. If
  /// m_current_particle is invalid, then move to the next bucket
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  bool get_valid_candidate() {
    LOG_CUDA(4, "\tget_valid_candidate:");
    while (m_current_particle == false) {
#ifdef __CUDA_ARCH__
      LOG_CUDA(4, "\tgo_to_next bucket (search_iterator):");
#else
      LOG(4, "\tgo_to_next bucket (search_iterator):");
#endif
      ++m_current_bucket;
      if (!get_valid_bucket()) {
        return false;
      }
      m_current_particle = m_query->get_bucket_particles(*m_current_bucket);
    }
    return true;
  }

  ///
  /// @brief increments m_current_particle and deals with any invalid
  /// iterators
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  bool go_to_next_candidate() {
#ifndef __CUDA_ARCH__
    LOG(4, "\tgo_to_next_candidate (search_iterator):");
#endif
    ++m_current_particle;
    return get_valid_candidate();
  }

  ///
  /// @brief checks that the current particle in m_current_particle is
  /// within the search distance
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  bool check_candidate() {
    LOG_CUDA(4, "\tcheck_candidate:");
    // const double_d& p = get<position>(*m_current_particle) +
    // m_particle_range.get_transpose();
    const double_d &p = get<position>(*m_current_particle);
    // const double_d& transpose = m_particle_range.get_transpose();
    double accum = 0;
    bool outside = false;
    for (size_t i = 0; i < Traits::dimension; i++) {
      m_dx[i] = p[i] - m_current_point[i];
      accum =
          detail::distance_helper<LNormNumber>::accumulate_norm(accum, m_dx[i]);
      if (accum > m_max_distance2) {
        outside = true;
        break;
      }
    }
#ifdef __CUDA_ARCH__
    if (3 <= ABORIA_LOG_LEVEL) {
      printf("\tcheck_candidate: m_r = (");
      for (int i = 0; i < Traits::dimension; ++i) {
        printf("%d,", m_current_point[i]);
      }
      printf(") other r = (");
      for (int i = 0; i < Traits::dimension; ++i) {
        printf("%d,", get<position>(*m_current_particle)[i]);
      }
      printf("). outside = %d", outside);
    }
#else
    LOG(3, "\tcheck_candidate: m_r = " << m_current_point << " other r = "
                                       << get<position>(*m_current_particle)
                                       << ". outside = " << outside);
#endif
    return !outside;
  }

  ///
  /// @brief main increment function, iterates internal iterators until a
  /// candidate particle is found (i.e. one within the search distance)
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  void increment() {
#ifndef __CUDA_ARCH__
    LOG(3, "\tincrement (search_iterator):");
#endif
    bool found_good_candidate = false;
    while (!found_good_candidate && (m_valid = go_to_next_candidate())) {
      found_good_candidate = check_candidate();
#ifndef __CUDA_ARCH__
      LOG(4, "\tfound_good_candidate = " << found_good_candidate);
#endif
    }
#ifndef __CUDA_ARCH__
    LOG(3, "\tend increment (search_iterator): valid = " << m_valid);
#endif
  }

  ///
  /// @brief dereference the iterator, returns a reference to the current
  /// candidate particle
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  reference dereference() const { return *m_current_particle; }
};

template <typename Query> class bucket_pair_iterator {

  typedef typename Query::traits_type Traits;
  static const unsigned int dimension = Query::dimension;

  typedef position_d<dimension> position;
  typedef Vector<double, dimension> double_d;
  typedef Vector<bool, dimension> bool_d;
  typedef Vector<int, dimension> int_d;

  bool m_valid;
  bool m_domain_domain;
  const Query *m_query;
  lattice_iterator<dimension> m_periodic;
  lattice_iterator<dimension> m_i;
  lattice_iterator<dimension> m_j;
  double_d m_position_offset;

public:
  typedef const std::tuple<const int_d &, const int_d &, const double_d &>
      *pointer;
  typedef std::forward_iterator_tag iterator_category;
  typedef const std::tuple<const int_d &, const int_d &, const double_d &>
      reference;
  typedef const std::tuple<const int_d, const int_d, const double_d> value_type;
  typedef std::ptrdiff_t difference_type;

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  static lattice_iterator<dimension> get_periodic_it(const bool_d is_periodic) {
    int_d start, end;
    for (size_t i = 0; i < dimension; ++i) {
      start[i] = is_periodic[i] ? -1 : 0;
      end[i] = is_periodic[i] ? 2 : 1;
    }
    lattice_iterator<dimension> it(start, end);
    it = int_d::Constant(0);
    return it;
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  bucket_pair_iterator() : m_valid(false) {}

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  bucket_pair_iterator(const Query &query)
      : m_valid(true), m_domain_domain(true), m_query(&query),
        m_periodic(get_periodic_it(m_query->get_periodic())),
        m_i(get_regular_buckets(query, *m_periodic)),
        m_position_offset(double_d::Constant(0)) {
#ifndef __CUDA_ARCH__
    LOG(3, "\tcreating bucket_pair_iterator. m_periodic = "
               << *m_periodic << " m_i = " << *m_i << " m_j = " << *m_j);
#endif
    if (m_i + 1 == false) {
#ifndef __CUDA_ARCH__
      LOG(3, "\tbucket_pair_iterator. end of i row");
#endif
      m_domain_domain = false;
      ++m_periodic;
      if (m_periodic == false) {
        m_valid = false;
        return;
      } else {
        m_i = get_regular_buckets(*m_query, *m_periodic);
        m_j = get_neighbouring_buckets(*m_query, *m_i, *m_periodic);
        m_position_offset = (*m_periodic) * (m_query->get_bounds().bmax -
                                             m_query->get_bounds().bmin);
        // domain-domain is always first, so never need m_j = m_i+1
      }
    } else {
      m_j = get_neighbouring_buckets(*m_query, *m_i, *m_periodic);
      m_j = *m_i;
      ++m_j;
    }
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  reference operator*() const { return dereference(); }
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  reference operator->() { return dereference(); }
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  bucket_pair_iterator &operator++() {
    increment();
    return *this;
  }
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  bucket_pair_iterator operator++(int) {
    bucket_pair_iterator tmp(*this);
    operator++();
    return tmp;
  }
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  size_t operator-(bucket_pair_iterator start) const {
    size_t count = 0;
    while (start != *this) {
      start++;
      count++;
    }
    return count;
  }
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  inline bool operator==(const bucket_pair_iterator &rhs) { return equal(rhs); }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  inline bool operator!=(const bucket_pair_iterator &rhs) {
    return !operator==(rhs);
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  inline bool operator==(const bool rhs) const { return equal(rhs); }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  inline bool operator!=(const bool rhs) const { return !operator==(rhs); }

private:
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  bool equal(bucket_pair_iterator const &other) const {
    return m_valid ? other.m_valid && m_i == other.m_i && m_j == other.m_j
                   : !other.m_valid;
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  bool equal(const bool other) const { return m_valid == other; }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  void increment() {
#ifndef __CUDA_ARCH__
    LOG(3, "\tbucket_pair_iterator. increment");
#endif

    m_j++;
    // end of j row
    if (m_j == false) {
#ifndef __CUDA_ARCH__
      LOG(3, "\tbucket_pair_iterator. end of j row");
#endif
      ++m_i;
      if (m_domain_domain ? m_i + 1 == false : m_i == false) {
#ifndef __CUDA_ARCH__
        LOG(3, "\tbucket_pair_iterator. end of i row");
#endif
        m_domain_domain = false;
        ++m_periodic;
        if (m_periodic == false) {
          m_valid = false;
          return;
        } else {
          m_i = get_regular_buckets(*m_query, *m_periodic);
          m_j = get_neighbouring_buckets(*m_query, *m_i, *m_periodic);
          m_position_offset = (*m_periodic) * (m_query->get_bounds().bmax -
                                               m_query->get_bounds().bmin);
          // domain-domain is always first, so never need m_j = m_i+1
        }
      } else {
        m_j = get_neighbouring_buckets(*m_query, *m_i, *m_periodic);
        if (m_domain_domain) {
          m_j = *m_i;
          ++m_j;
          // we already know m_i + 1 is not false
        }
      }
    }
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  lattice_iterator<dimension>
  get_neighbouring_buckets(const Query &query, const int_d &bucket) const {
#ifndef __CUDA_ARCH__
    LOG(4, "\tget_neighbouring_buckets: ");
#endif
    int_d start = bucket - 1;
    int_d end = bucket + 1;

    bool no_buckets = false;
    for (size_t i = 0; i < dimension; i++) {
      if (start[i] < 0) {
        start[i] = 0;
      } else if (start[i] > query.get_end_bucket()[i]) {
        no_buckets = true;
        start[i] = query.get_end_bucket()[i];
      }
      if (end[i] < 0) {
        no_buckets = true;
        end[i] = 0;
      } else if (end[i] > query.get_end_bucket()[i]) {
        end[i] = query.get_end_bucket()[i];
      }
    }
#ifndef __CUDA_ARCH__
    LOG(4, "\tget_neighbouring_buckets: looking in bucket "
               << bucket << ". start = " << start << " end = " << end
               << " no_buckets = " << no_buckets);
#endif
    if (no_buckets) {
      return lattice_iterator<dimension>();
    } else {
      return lattice_iterator<dimension>(start, end + 1);
    }
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  lattice_iterator<dimension>
  get_neighbouring_buckets(const Query &query, const int_d &bucket,
                           const int_d &quadrant) const {
    return get_neighbouring_buckets(
        query,
        // bucket+quadrant*(query.get_end_bucket()+1)
        // TODO: why do I need to cast this???!?!?!?
        (bucket + quadrant * (query.get_end_bucket() + 1))
            .template cast<int>());
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  lattice_iterator<dimension> get_regular_buckets(const Query &query,
                                                  const int_d &quadrant) const {
#ifndef __CUDA_ARCH__
    LOG(4, "\tget_ghost_buckets: " << quadrant);
#endif
    int_d start = int_d::Constant(0);
    int_d end(query.get_end_bucket());

    for (size_t i = 0; i < Traits::dimension; i++) {
      if (!query.get_periodic()[i]) {
        ASSERT_CUDA(quadrant[i] == 0);
      } else {
        ASSERT_CUDA(quadrant[i] <= 1 && quadrant[i] >= -1);
        if (quadrant[i] > 0) {
          start[i] = 0;
          end[i] = 0;
        } else if (quadrant[i] < 0) {
          start[i] = query.get_end_bucket()[i];
          end[i] = query.get_end_bucket()[i];
        }
      }
    }

#ifndef __CUDA_ARCH__
    LOG(4, "\tget_ghost_buckets: looking in quadrant"
               << quadrant << ". start = " << start << " end = " << end);
#endif

    return lattice_iterator<dimension>(start, end + 1);
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  reference dereference() const {
    return reference(*m_i, *m_j, m_position_offset);
  }
};

/*
template <typename query_iterator>
iterator_range<query_iterator>
get_buckets_near_point(const double_d &position, const double
max_distance, detail::cell_list_tag) {
}

template <typename query_iterator>
iterator_range<query_iterator>
get_buckets_near_point(const double_d &position, const double
max_distance, detail::kd_tree_tag) {
}
*/

///
/// @brief returns a @ref search_iterator that iterates over all the
/// particles within a given distance around a point
///
/// @tparam Query the query object type
/// @tparam search_iterator<Query, 1> the search iterator type
/// @param query the query object
/// @param centre the central point of the search
/// @param max_distance the maximum distance to search around @p centre
///
template <int LNormNumber, typename Query,
          typename SearchIterator = search_iterator<Query, LNormNumber>>
CUDA_HOST_DEVICE SearchIterator
distance_search(const Query &query, const typename Query::double_d &centre,
                const double max_distance) {
  return SearchIterator(query, centre, max_distance);
}

///
/// @copydoc distance_search()
///
/// Uses the chebyshev distance
/// <https://en.wikipedia.org/wiki/Chebyshev_distance>
///
///
template <typename Query, typename SearchIterator = search_iterator<Query, -1>>
CUDA_HOST_DEVICE SearchIterator
chebyshev_search(const Query &query, const typename Query::double_d &centre,
                 const double max_distance) {
  return SearchIterator(query, centre, max_distance);
}

///
/// @copydoc distance_search()
///
/// Uses the manhatten distance
/// <https://en.wikipedia.org/wiki/Taxicab_geometry>
///
///
template <typename Query, typename SearchIterator = search_iterator<Query, 1>>
CUDA_HOST_DEVICE SearchIterator
manhatten_search(const Query &query, const typename Query::double_d &centre,
                 const double max_distance) {
  return SearchIterator(query, centre, max_distance);
}

///
/// @copydoc distance_search()
///
/// Uses the euclidean distance
/// <https://en.wikipedia.org/wiki/Euclidean_distance>
///
///
template <typename Query, typename SearchIterator = search_iterator<Query, 2>>
CUDA_HOST_DEVICE SearchIterator
euclidean_search(const Query &query, const typename Query::double_d &centre,
                 const double max_distance) {
  return SearchIterator(query, centre, max_distance);
}

///
/// @brief returns a @ref bucket_pair_iterator that iterates through all
/// the neighbouring buckets (i.e. buckets that are touching) within a
/// domain. Note that this will only work for cell list spatial data
/// structures
///
/// @tparam Query query object type (must be @ref CellListQuery or @ref
/// CellListOrderedQuery)
/// @tparam bucket_pair_iterator<Query> iterator type returned
/// @param query the query object
///
template <typename Query, typename Iterator = bucket_pair_iterator<Query>>
CUDA_HOST_DEVICE Iterator get_neighbouring_buckets(const Query &query) {
  return Iterator(query);
}

///
/// @brief returns a @ref search_iterator that iterates down the levels of the
/// spatial structure tree, finding buckets that are within a given distance
/// from each other.
///
/// @tparam LNormNumber the p-norm number of the distance
/// @tparam QueryA the query object a type
/// @tparam QueryB the query object b type
/// @tparam search_bf_iterator<QueryA,QueryB,LNormNumber> the search
/// iterator type
/// @param query_a the query object a
/// @param query_b the query object b
/// @param max_distance the maximum distance to search around @p centre
///
template <
    int LNormNumber, typename QueryA, typename QueryB,
    typename SearchIterator = search_bf_iterator<QueryA, QueryB, LNormNumber>>
CUDA_HOST_DEVICE SearchIterator distance_pair_search(
    const QueryA &query_a, const QueryB &query_b, const double max_distance) {
  return SearchIterator(query_a.get_children(), query_b.get_children(), query_a,
                        query_b, max_distance, false);
}

///
/// @copydoc distance_pair_search()
///
/// Uses the euclidean distance
/// <https://en.wikipedia.org/wiki/Euclidean_distance>
///
///
template <typename QueryA, typename QueryB,
          typename SearchIterator = search_bf_iterator<QueryA, QueryB, 2>>
CUDA_HOST_DEVICE SearchIterator euclidean_pair_search(
    const QueryA &query_a, const QueryB &query_b, const double max_distance) {
  return SearchIterator(query_a.get_children(), query_b.get_children(), query_a,
                        query_b, max_distance, false);
}

} // namespace Aboria

#endif
