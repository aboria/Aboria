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

#ifndef KDTREE_H_
#define KDTREE_H_

#include "Get.h"
#include "Log.h"
#include "NeighbourSearchBase.h"
#include "SpatialUtil.h"
#include "Traits.h"
#include "Vector.h"
#include <boost/iterator/iterator_facade.hpp>
#include <iostream>
#include <set>
#include <vector>

namespace Aboria {

/// \brief KdTree
///
template <typename Traits>
class Kdtree : public neighbour_search_base<Kdtree<Traits>, Traits,
                                            KdtreeQuery<Traits>> {

  typedef typename Traits::double_d double_d;
  typedef typename Traits::position position;
  typedef typename Traits::vector_int vector_int;
  typedef typename Traits::iterator iterator;
  typedef typename Traits::unsigned_int_d unsigned_int_d;
  static const unsigned int dimension = Traits::dimension;

  typedef neighbour_search_base<Kdtree<Traits>, Traits, KdtreeQuery<Traits>>
      base_type;
  friend base_type;

public:
  Kdtree() : base_type() {

    this->m_query.m_root = m_kd_tree.get_root_node();
    this->m_query.m_dummy_root.child1 = this->m_query.m_root;
    this->m_query.m_dummy_root.child2 = this->m_query.m_root;
    this->m_query.m_dummy_root.node_type.sub.divfeat = 0;
    this->m_query.m_dummy_root.node_type.sub.divlow =
        this->m_query.m_bounds.bmin[0];
    this->m_query.m_dummy_root.node_type.sub.divhigh =
        this->m_query.m_bounds.bmin[0];
    this->m_query.m_number_of_buckets = m_kd_tree.size_nodes();
    this->m_query.m_number_of_levels = m_kd_tree.size_levels();
  }

  //~Kdtree() {}

  static constexpr bool ordered() { return true; }

  void print_data_structure() const { print_tree(m_kd_tree.get_root_node()); }

private:
  void set_domain_impl() {
    m_kd_tree.set_leaf_max_size(this->m_n_particles_in_leaf);

    this->m_query.m_bounds.bmin = this->m_bounds.bmin;
    this->m_query.m_bounds.bmax = this->m_bounds.bmax;
    this->m_query.m_periodic = this->m_periodic;
  }

  void update_iterator_impl() {}

  void print_tree(const node_type *nodes) const {}

  void update_positions_impl(iterator update_begin, iterator update_end,
                             const int new_n,
                             const bool call_set_domain = true) {
    ASSERT(update_begin == this->m_particles_begin &&
               update_end == this->m_particles_end,
           "error should be update all");
    const size_t num_points = this->m_alive_indices.size();

    // copy alive indicies to per dimension arrays
    for (size_t i = 0; i < D; ++i) {
      m_sorted_indices[i].resize(n);
      detail::copy(this->m_alive_indices.begin(), this->m_alive_indices.end(),
                   m_sorted_indices[i].begin());
    }

    // sort per dimension arrays by position (i could be compile time?)
    for (size_t i = 0; i < D; ++i) {
      auto get_dim_i = [=](const double_d &p) { return p[i]; };
      detail::sort_by_key(detail::make_transform_iterator(
                              detail::make_permutation_iterator(
                                  get<position>(this->m_particles_begin),
                                  m_sorted_indices[i].begin()),
                              get_dim_i),
                          detail::make_transform_iterator(
                              detail::make_permutation_iterator(
                                  get<position>(this->m_particles_begin),
                                  m_sorted_indices[i].end()),
                              get_dim_i),
                          m_tags.begin());
    }

    build_tree();

    // std::swap(this->m_order,m_kd_tree.get_vind());

    this->m_query.m_root = m_kd_tree.get_root_node();
    this->m_query.m_dummy_root.child1 = this->m_query.m_root;
    this->m_query.m_dummy_root.child2 = this->m_query.m_root;
    this->m_query.m_dummy_root.node_type.sub.divfeat = 0;
    this->m_query.m_dummy_root.node_type.sub.divlow =
        this->m_query.m_bounds.bmin[0];
    this->m_query.m_dummy_root.node_type.sub.divhigh =
        this->m_query.m_bounds.bmin[0];
    this->m_query.m_number_of_buckets = m_kd_tree.size_nodes();
    this->m_query.m_number_of_levels = m_kd_tree.size_levels();

    print_tree(m_kd_tree.get_root_node());
  }

  const KdtreeQuery<Traits> &get_query_impl() const { return m_query; }

  KdtreeQuery<Traits> &get_query_impl() { return m_query; }

private:
  struct calculate_max_span {
    const vector_d *position;
    const int *sorted_indices[D];
    const int *nodes_begin;
    const int *nodes_end;
    double *split;
    int *split_d;
    calculate_max_span(const vector_double_d &_position,
                       const std::array<vector_int, D> &_sorted_indices,
                       const vector_int &_nodes_begin,
                       const vector_int &_nodes_end, vector_double &_split,
                       vector_int &_split_d) {
      position = iterator_to_raw_pointer(_position.begin());
      nodes_begin = iterator_to_raw_pointer(_nodes_begin.begin());
      nodes_end = iterator_to_raw_pointer(_nodes_end).begin());
      split = iterator_to_raw_pointer(_split.begin());
      split_d = iterator_to_raw_pointer(_split_d.begin());
      for (size_t i = 0; i < D; ++i) {
        sorted_indices[i] = make_raw_pointer(_sorted_indices[i]);
      }
    }
    void operator(const int i) {
      double max_span = position[m_sorted_indices[0][m_nodes_end[i]]][0] -
                        position[m_sorted_indices[0][m_nodes_begin[i]]][0];
      int max_d = 0;
      for (size_t d = 1; d < D; ++d) {
        const double span = position[m_sorted_indices[d][nodes_end[i]]][d] -
                            position[m_sorted_indices[d][nodes_begin[i]]][d];
        if (span > max_span) {
          max_d = d;
          max_span = span;
        }
      }
      split[i] = m_sorted_indices[max_d][m_nodes_begin[i]] + max_span / 2;
      split_d[i] = max_d;
    }
  };
  void build_tree() {
    const size_t num_points = this->m_alive_indices.size();
    m_active_nodes.resize(1);
    m_active_nodes[0] =
        0; // >= 0 index of particle, < 0 index of 1st child node
    vector_int node_mask(num_points, 0);
    while (!m_active_nodes.empty()) {
      m_nodes.append(m_active_nodes);

      // get min bounds (TODO: can combine?)
      double_d_vector min_bounds(m_active_nodes.size());
      detail::reduce_by_key(
          node_mask[i].begin(), node_mask[i].end(),
          get<position>(this->m_particles_begin), m_active_nodes.begin(),
          min_bounds.begin(), detail::equal_to<int>(),
          [](const double_d &a, const double_d &b) { return a.min(b); });

      // get max bounds
      double_d_vector max_bounds(m_active_nodes.size());
      detail::reduce_by_key(
          node_mask[i].begin(), node_mask[i].end(),
          get<position>(this->m_particles_begin), m_active_nodes.begin(),
          max_bounds.begin(), detail::equal_to<int>(),
          [](const double_d &a, const double_d &b) { return a.max(b); });

      // calc split (median of largest dimension span)
      vector_double split(min_bounds.size());
      vector_int split_d(min_bounds.size());
      detail::transform(detail::make_zip_iterator(detail::make_tuple(
                            min_bounds.begin(), max_bounds.begin())),
                        detail::make_zip_iterator(detail::make_tuple(
                            min_bounds.end(), max_bounds.end())),
                        detail::make_zip_iterator(
                            detail::make_tuple(split_d.begin(), split.begin())),
                        [](const auto &min_max) {

                        });

      // partition m_sorted_indices by split
      //(Sengupta, S., Harris, M., Zhang, Y., & Owens, J. D. (2007). Scan
      // primitives for GPU computing. Graphics …, 97–106.
      // http://doi.org/10.2312/EGGH/EGGH07/097-106)
      vector_int e(m_sorted_indices[i].size());
      vector_int f(m_sorted_indices[i].size());
      vector_int addr(m_sorted_indices[i].size());
      // [t f t f f t f] #in
      // [0 1 0 1 1 0 1] #e = set 1 in false elts.
      detail::for_each_n(
          m_sorted_indices[i].begin(), m_sorted_indices[i].end(),
          [_position =
               iterator_to_raw_pointer(get<position>(this->m_particles_begin)),
           _node_mask = iterator_to_raw_pointer(node_mask[i].begin()),
           _split_d = iterator_to_raw_pointer(split_d),
           _split = iterator_to_raw_pointer(split)](const int j) {
            const int split_d = _split_d[_node_mask[j]];
            const double split = _split[_node_mask[j]];
            return static_cast<int>(_position[j][split_d] < split);
          });
      // [0 0 1 1 2 3 3] #f = enumerate with false=1
      detail::exclusive_scan_by_key(node_mask[i].begin(), node_mask[i].end(),
                                    e.begin(), f.begin());
      //              4  # add two last elts. in e, f
      detail::for_each_n(m_active_nodes.begin(), m_active_nodes.end(),
                         calculate_max_span(m_sorted_indices, split, split_d));

      //                 # ==total # of falses
      //                 # set as shared variable NF
      //[0 1 2 3 4 5 6]  # each thread knows its id
      //[4 5 5 6 6 6 7]  # t = id - f + NF
      //[4 0 5 1 2 6 3]  # addr = e ? f : t
      //[f f f f t t t]  # out[addr] = in (scatter)

      // remove leaf nodes from active_nodes
      // add remainder to active_nodes
    }
  }

  std::array<vector_int, D> m_sorted_indices;
  vector_int m_active_nodes;
  vector_int m_nodes_begin;
  vector_int m_nodes_end;
  vector_double m_node_split;
  vector_int m_node_d;
  KdtreeQuery<Traits> m_query;
}; // namespace Aboria

template <typename Traits> class KdtreeChildIterator {
  static const unsigned int D = Traits::dimension;
  typedef detail::nanoflann_kd_tree_type<Traits> kd_tree_type;
  typedef Vector<double, D> double_d;
  typedef Vector<int, D> int_d;
  typedef Vector<bool, D> bool_d;
  typedef bbox<D> box_type;

  int m_high;
  const typename kd_tree_type::Node *m_index;
  box_type m_bounds;

public:
  typedef const typename kd_tree_type::Node *pointer;
  typedef std::forward_iterator_tag iterator_category;
  typedef const typename kd_tree_type::Node value_type;
  typedef const typename kd_tree_type::Node &reference;
  typedef std::ptrdiff_t difference_type;

  KdtreeChildIterator() : m_high(2), m_index(nullptr) {}

  KdtreeChildIterator(pointer start, const box_type &bounds)
      : m_high(0), m_index(start), m_bounds(bounds) {
    ASSERT(start != nullptr, "start pointer should not be null");
    ASSERT(start->child1 != nullptr, "start pointer should point to leaf");
  }

  void go_to(const double_d &position) {
    const int i = m_index->node_type.sub.divfeat;
    const double diff = position[i] - 0.5 * (m_index->node_type.sub.divhigh +
                                             m_index->node_type.sub.divhigh);
    ASSERT(position[i] < m_bounds.bmax[i], "position out of bounds");
    ASSERT(position[i] >= m_bounds.bmin[i], "position out of bounds");
    if (diff < 0) {
      m_high = 0;
    } else {
      m_high = 1;
    }
  }

  int get_child_number() const { return m_high; }

  bool is_high() const { return m_high > 0; }

  box_type get_bounds() const {
    box_type ret = m_bounds;
    const int i = m_index->node_type.sub.divfeat;
    if (is_high()) {
      ret.bmin[i] = m_index->node_type.sub.divhigh;
    } else {
      ret.bmax[i] = m_index->node_type.sub.divlow;
    }
    return ret;
  }

  reference operator*() const { return dereference(); }

  reference operator->() const { return dereference(); }

  KdtreeChildIterator &operator++() {
    increment();
    return *this;
  }

  KdtreeChildIterator operator++(int) {
    KdtreeChildIterator tmp(*this);
    operator++();
    return tmp;
  }

  inline bool operator==(const KdtreeChildIterator &rhs) const {
    return equal(rhs);
  }

  inline bool operator!=(const KdtreeChildIterator &rhs) const {
    return !operator==(rhs);
  }

  inline bool operator==(const bool rhs) const { return equal(rhs); }

  inline bool operator!=(const bool rhs) const { return !operator==(rhs); }

private:
  bool equal(KdtreeChildIterator const &other) const {
    return m_index == other.m_index && m_high == other.m_high;
  }

  bool equal(const bool other) const { return (m_high < 2) == other; }

  reference dereference() const {
    if (is_high()) {
      return *m_index->child2;
    } else {
      return *m_index->child1;
    }
  }

  void increment() { ++m_high; }
};

/// @copydetails NeighbourQueryBase
///
/// @brief This is a query object for the @ref Kdtree spatial data
/// structure
///
template <typename Traits> struct KdtreeQuery {
  const static unsigned int dimension = Traits::dimension;
  const static unsigned int m_max_tree_depth = 32 - 2;
  typedef detail::nanoflann_kd_tree_type<Traits> kd_tree_type;
  typedef typename kd_tree_type::Node value_type;
  typedef const value_type &reference;
  typedef const value_type *pointer;

  typedef Traits traits_type;
  typedef typename Traits::raw_pointer raw_pointer;
  typedef typename Traits::double_d double_d;
  typedef typename Traits::bool_d bool_d;
  typedef typename Traits::int_d int_d;
  typedef typename Traits::unsigned_int_d unsigned_int_d;
  template <int LNormNumber>
  using query_iterator = tree_query_iterator<KdtreeQuery, LNormNumber>;
  typedef value_type *root_iterator;
  typedef depth_first_iterator<KdtreeQuery> all_iterator;
  typedef nanoflann_child_iterator<Traits> child_iterator;
  typedef ranges_iterator<Traits> particle_iterator;
  typedef bbox<dimension> box_type;

  bool_d m_periodic;
  bbox<dimension> m_bounds;
  raw_pointer m_particles_begin;
  raw_pointer m_particles_end;
  size_t m_number_of_buckets;
  size_t m_number_of_levels;

  value_type *m_root;
  value_type m_dummy_root;

  size_t *m_id_map_key;
  size_t *m_id_map_value;

  const box_type &get_bounds() const { return m_bounds; }
  const bool_d &get_periodic() const { return m_periodic; }

  /*
   * functions for id mapping
   */
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  raw_pointer find(const size_t id) const {
    const size_t n = number_of_particles();
    size_t *last = m_id_map_key + n;
    size_t *first = detail::lower_bound(m_id_map_key, last, id);
    if ((first != last) && !(id < *first)) {
      return m_particles_begin + m_id_map_value[first - m_id_map_key];
    } else {
      return m_particles_begin + n;
    }
  }

  /*
   * functions for tree_query_iterator
   */
  static bool get_max_levels() { return 5; }
  static bool is_leaf_node(reference bucket) {
    return (bucket.child1 == NULL) && (bucket.child2 == NULL);
  }
  static bool is_tree() { return true; }
  static size_t get_dimension_index(reference bucket) {
    return bucket.node_type.sub.divfeat;
  }
  static double get_cut_low(reference bucket) {
    return bucket.node_type.sub.divlow;
  }
  static double get_cut_high(reference bucket) {
    return bucket.node_type.sub.divhigh;
  }
  static pointer get_child1(pointer bucket) { return bucket->child1; }
  static pointer get_child2(pointer bucket) { return bucket->child2; }
  /*
   * end functions for tree_query_iterator
   */

  child_iterator get_children() const {
    if (m_root == nullptr) { // empty tree, return a false child iterator
      return child_iterator();
    } else if (is_leaf_node(*m_root)) { // leaf root, return a child iterator
                                        // pointing to the leaf
      return ++child_iterator(&m_dummy_root, m_bounds);
    } else {
      return child_iterator(m_root, m_bounds);
    }
  }

  child_iterator get_children(reference bucket, const box_type &bounds) const {
    CHECK(&bucket == m_root, "bucket should be a root bucket");
    return child_iterator(m_root, m_bounds);
  }

  static child_iterator get_children(const child_iterator &ci) {
    if (!is_leaf_node(*ci)) {
      return child_iterator(&(*ci), ci.get_bounds());
    } else {
      return child_iterator();
    }
  }

  static const box_type get_bounds(const child_iterator &ci) {
    return ci.get_bounds();
  }

  template <int LNormNumber>
  static const box_type get_bounds(const query_iterator<LNormNumber> &qi) {
    return qi.get_bounds();
  }

  static const box_type get_bounds(const all_iterator &ai) {
    return ai.get_bounds();
  }

  friend std::ostream &operator<<(std::ostream &os, reference bucket) {
    if (is_leaf_node(bucket)) {
      os << "Leaf node";
    } else {
      os << "Node";
    }
    os << " with bounding box " << get_bucket_bbox(bucket) << std::endl;
    return os;
  }

  particle_iterator get_bucket_particles(reference bucket) const {
#ifndef __CUDA_ARCH__
    LOG(4, "\tget_bucket_particles: looking in bucket with idx = "
               << get_dimension_index(bucket)
               << " start index = " << bucket.node_type.lr.left
               << " end index = " << bucket.node_type.lr.right);
#endif
    if (!is_leaf_node(bucket)) {
      return particle_iterator();
    }

    return particle_iterator(m_particles_begin + bucket.node_type.lr.left,
                             m_particles_begin + bucket.node_type.lr.right);
  }

  /*
  static double_d
  get_bucket_bounds_low(reference bucket) {
      double_d low;
      for (size_t i = 0; i < dimension; ++i) {
          low[i] = bucket.bbox[i].low;
      }
      return low;
  }

  static double_d
  get_bucket_bounds_high(reference bucket) {
      double_d high;
      for (size_t i = 0; i < dimension; ++i) {
          high[i] = bucket.bbox[i].high;
      }
      return high;
  }
  */

  void get_bucket(const double_d &position, pointer &bucket,
                  box_type &bounds) const {
    child_iterator i = get_children();
    i.go_to(position);

    while (!is_leaf_node(*i)) {
      i = get_children(i);
      i.go_to(position);
    }

    bucket = &(*i);
    bounds = i.get_bounds();
  }

  size_t get_bucket_index(reference bucket) const { return bucket.index; }

  size_t number_of_buckets() const { return m_number_of_buckets; }

  template <int LNormNumber>
  query_iterator<LNormNumber>
  get_buckets_near_point(const double_d &position,
                         const double max_distance) const {
#ifndef __CUDA_ARCH__
    LOG(4, "\tget_buckets_near_point: position = "
               << position << " max_distance= " << max_distance);
#endif
    return query_iterator<LNormNumber>(get_children(), position,
                                       double_d::Constant(max_distance),
                                       m_number_of_levels, this);
  }

  template <int LNormNumber>
  query_iterator<LNormNumber>
  get_buckets_near_point(const double_d &position,
                         const double_d &max_distance) const {
#ifndef __CUDA_ARCH__
    LOG(4, "\tget_buckets_near_point: position = "
               << position << " max_distance= " << max_distance);
#endif
    return query_iterator<LNormNumber>(get_children(), position, max_distance,
                                       m_number_of_levels, this);
  }

  iterator_range<root_iterator> get_root_buckets() const {
    return iterator_range<root_iterator>(m_root, m_root + 1);
  }

  all_iterator get_subtree(const child_iterator &ci) const {
    return all_iterator(get_children(ci), m_number_of_levels, this);
  }

  all_iterator get_subtree() const {
    return all_iterator(get_children(), m_number_of_levels, this);
  }

  size_t number_of_particles() const {
    return m_particles_end - m_particles_begin;
  }

  raw_pointer get_particles_begin() const { return m_particles_begin; }

  unsigned number_of_levels() const { return m_number_of_levels; }

  /*
  iterator_range<theta_iterator> get_theta_buckets(const reference bucket)
  const { return iterator_range<theta_iterator>(
  theta_iterator(m_root,bucket), theta_iterator()
              );
  }
  */
};

} // namespace Aboria

#endif /* NANOFLANN_ADAPTOR_H_ */
