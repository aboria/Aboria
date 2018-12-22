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

#ifndef SUBSAMPLE_TREE_H_
#define SUBSAMPLE_TREE_H_

#include "NeighbourSearchBase.h"
#include "Preconditioners/detail/Preconditioners.h"
#include <boost/iterator/iterator_facade.hpp>
#include <iostream>
#include <set>
#include <vector>

namespace Aboria {

template <typename Query, typename Traits = typename Query::traits_type>
class subsample_iterator {
  using iterator = subsample_iterator<Query, Traits>;
  typedef typename Traits::double_d double_d;
  typedef typename Traits::bool_d bool_d;
  typedef typename Traits::value_type p_value_type;
  typedef typename Traits::raw_reference p_reference;
  typedef typename Traits::raw_pointer p_pointer;

  const Query *m_query;

  const size_t *m_indicies_current;
  const size_t *m_indicies_end;

public:
  typedef Traits traits_type;
  typedef const p_pointer pointer;
  typedef std::random_access_iterator_tag iterator_category;
  typedef const p_reference reference;
  typedef const p_reference value_type;
  typedef std::ptrdiff_t difference_type;

  CUDA_HOST_DEVICE
  subsample_iterator() {}

  CUDA_HOST_DEVICE
  subsample_iterator(const Query *query, const size_t *begin, const size_t *end)
      : m_query(query), m_indicies_current(begin), m_indicies_end(end) {}

  size_t distance_to_end() const { return m_indicies_end - m_indicies_current; }

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
  iterator operator+(int n) {
    iterator tmp(*this);
    tmp.increment(n);
    return tmp;
  }

  CUDA_HOST_DEVICE
  size_t operator-(iterator start) const {
    return m_indicies_current - start.m_indicies_current;
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
  bool equal(iterator const &other) const {
    return m_indicies_current == other.m_indicies_current;
  }

  CUDA_HOST_DEVICE
  bool equal(const bool other) const {
    return (m_indicies_current < m_indicies_end) == other;
  }

  CUDA_HOST_DEVICE
  reference dereference() const {
    return m_query->get_particles_begin()[*m_indicies_current];
  }

  CUDA_HOST_DEVICE
  void increment() { ++m_indicies_current; }

  CUDA_HOST_DEVICE
  void increment(const int n) { m_indicies_current += n; }
};

/// @copydetails NeighbourQueryBase
///
/// @brief This is a query object for the @ref SubsampleTree spatial data
/// structure
///
template <typename Query, typename Traits = typename Query::traits_type>
struct SubsampleTreeQuery {
  const static unsigned int dimension = Query::dimension;
  typedef Traits traits_type;

  using bool_d = typename Traits::bool_d;
  using double_d = typename Traits::double_d;
  using raw_pointer = typename Traits::raw_pointer;
  using position = typename Traits::position;

  template <int LNormNumber>
  using query_iterator = typename Query::template query_iterator<LNormNumber>;

  template <int LNormNumber>
  using bounds_query_iterator =
      typename Query::template bounds_query_iterator<LNormNumber>;

  using all_iterator = typename Query::all_iterator;
  using breadth_first_iterator = typename Query::breadth_first_iterator;
  using child_iterator = typename Query::child_iterator;

  using value_type = typename Query::value_type;
  using reference = typename Query::reference;
  using pointer = typename Query::pointer;
  using particle_iterator = subsample_iterator<Query>;
  using box_type = typename Query::box_type;

  using storage_vector_t = detail::storage_vector_type<size_t>;
  using indicies_t =
      typename traits_type::template vector_type<storage_vector_t>::type;

  const Query &m_query;
  const size_t m_max_points_per_bucket;
  indicies_t m_indicies;

  SubsampleTreeQuery(const Query &query, const size_t max_points)
      : m_query(query), m_max_points_per_bucket(max_points) {

    m_indicies.resize(m_query.number_of_buckets());
    auto tree = make_flat_tree(m_query);
    // go up the tree subsampling points
    for (auto level = tree.rbegin(); level != tree.rend(); ++level) {
      // std::cout << "working on level of size" << level->size() << std::endl;
      detail::for_each(
          level->begin(), level->end(),
          [max = m_max_points_per_bucket, query = m_query,
           indicies = iterator_to_raw_pointer(m_indicies.begin())](
              const child_iterator ci) {
            const size_t index = query.get_bucket_index(*ci);
            const auto p0 = get<position>(query.get_particles_begin());
            // std::cout << "working on bucket " << index << std::endl;
            if (query.is_leaf_node(*ci)) {
              // std::cout << "is leaf: ";
              // get all leaf particles
              auto p = query.get_bucket_particles(*ci);
              const size_t n = p.distance_to_end();
              indicies[index].resize(n);
              int i = 0;
              for (; p != false; ++p) {
                const size_t pindex = &(get<position>(*p)) - p0;
                // std::cout << pindex << " ";
                indicies[index][i++] = pindex;
              }
              // std::cout << std::endl;
            } else {
              // std::cout << "is not leaf: " << std::endl;
              // count children particles
              int nchildren = 0;
              for (auto child_ci = query.get_children(ci); child_ci != false;
                   ++child_ci) {
                const size_t cindex = query.get_bucket_index(*child_ci);
                nchildren += indicies[cindex].size();
              }

              indicies[index].resize(nchildren);

              // get all particles from children
              int ichildren = 0;
              for (auto child_ci = query.get_children(ci); child_ci != false;
                   ++child_ci) {
                const size_t cindex = query.get_bucket_index(*child_ci);
                // get all previously sampled particles
                // for (auto pindex : indicies[cindex]) {
                // std::cout << pindex << " ";
                //}
                detail::copy(indicies[cindex].begin(), indicies[cindex].end(),
                             indicies[index].begin() + ichildren);
                ichildren += indicies[cindex].size();
              }
              // std::cout.flush();
              ASSERT_CUDA(ichildren == nchildren);
              // std::cout << std::endl;
            }

#if defined(__CUDACC__)
            thrust::default_random_engine gen;
#else
            generator_type gen;
#endif
            // advance forward so no random streams intersect
            gen.discard(index * max);

            // sample max_points_per_bucket
            detail::random_unique(indicies[index].begin(),
                                  indicies[index].end(), max, gen);
            indicies[index].resize(std::min(max, indicies[index].size()));
            // std::cout << "chosen: ";
            // for (auto index : indicies[index]) {
            //  std::cout << index << " ";
            //}
            // std::cout << std::endl;
          });
    }
  }

  const box_type &get_bounds() const { return m_query.get_bounds(); }
  const bool_d &get_periodic() const { return m_query.get_periodic(); }

  /*
   * functions for id mapping
   */
  CUDA_HOST_DEVICE
  raw_pointer find(const size_t id) const { return m_query.find(id); }

  /*
   * functions for tree_query_iterator
   */
  bool is_leaf_node(reference bucket) const {
    return m_query.is_leaf_node(bucket);
  }
  bool is_tree() { return m_query.is_tree(); }

  /*
   * end functions for tree_query_iterator
   */

  ///
  /// @copydoc NeighbourQueryBase::get_root() const
  ///
  child_iterator get_root() const { return m_query.get_root(); }

  child_iterator get_children() const { return m_query.get_children(); }

  child_iterator get_children(const child_iterator &ci) const {
    return m_query.get_children(ci);
  }

  ///
  /// @copydoc NeighbourQueryBase::num_children() const
  ///
  size_t num_children() const { return m_query.num_children(); }

  size_t num_children(const child_iterator &ci) const {
    return m_query.num_children(ci);
  }

  const box_type get_bounds(const child_iterator &ci) const {
    return m_query.get_bounds(ci);
  }

  particle_iterator get_bucket_particles(reference bucket) const {
    // now both leafs and internal bucket "have" particles
    const size_t index = get_bucket_index(bucket);
    ASSERT_CUDA(index < m_indicies.size());
    return particle_iterator(&m_query, m_indicies[index].data(),
                             m_indicies[index].data() +
                                 m_indicies[index].size());
  }

  void go_to(const double_d &position, child_iterator &ci) const {
    return m_query.go_to(position, ci);
  }

  child_iterator get_bucket(const double_d &position) const {
    return m_query.get_bucket(position);
  }

  size_t get_parent_index(const child_iterator &ci) const {
    return m_query.get_parent_index(ci);
  }

  const box_type &get_parent_bounds(const child_iterator &ci) const {
    return m_query.get_parent_bounds(ci);
  }

  size_t get_bucket_index(reference bucket) const {
    return m_query.get_bucket_index(bucket);
  }

  size_t number_of_buckets() const { return m_query.number_of_buckets(); }

  template <int LNormNumber>
  query_iterator<LNormNumber>
  get_buckets_near_point(const double_d &position,
                         const double max_distance) const {
    return m_query.get_buckets_near_point<LNormNumber>(position, max_distance);
  }

  ///
  /// @copydoc NeighbourQueryBase::get_buckets_near_bucket()
  ///
  template <int LNormNumber>
  bounds_query_iterator<LNormNumber>
  get_buckets_near_bucket(const box_type &bounds,
                          const double max_distance) const {
    return m_query.get_buckets_near_bucket<LNormNumber>(bounds, max_distance);
  }

  template <int LNormNumber>
  query_iterator<LNormNumber>
  get_buckets_near_point(const double_d &position,
                         const double_d &max_distance) const {
    return m_query.get_buckets_near_point<LNormNumber>(position, max_distance);
  }

  all_iterator get_subtree(const child_iterator &ci) const {
    return m_query.get_subtree(ci);
  }

  all_iterator get_subtree() const { return m_query.get_subtree(); }

  ///
  /// @copydoc NeighbourQueryBase::get_breadth_first(const child_iterator&)
  /// const
  ///
  breadth_first_iterator breadth_first(const child_iterator &ci) const {
    return m_query.breadth_first(ci);
  }

  ///
  /// @copydoc NeighbourQueryBase::get_breadth_first() const
  ///
  breadth_first_iterator breadth_first() const {
    return m_query.breadth_first();
  }

  size_t number_of_particles() const { return m_query.number_of_particles(); }

  raw_pointer get_particles_begin() const {
    return m_query.get_particles_begin();
  }

  unsigned number_of_levels() const { return m_query.number_of_levels(); }
}; // namespace Aboria

template <typename Query>
auto create_subsample_query(const Query &query, const size_t max) {
  return SubsampleTreeQuery<Query>(query, max);
}

} // namespace Aboria

#endif /* KDTREE_H_ */
