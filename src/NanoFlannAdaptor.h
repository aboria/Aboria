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

#ifndef NANOFLANN_ADAPTOR_H_
#define NANOFLANN_ADAPTOR_H_

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

template <typename Traits> class KdtreeNanoflann;

template <typename Traits> struct KdtreeNanoflannQuery;

template <typename Traits> class nanoflann_child_iterator;

} // namespace Aboria

#include "nanoflann/nanoflann.hpp"

namespace Aboria {
namespace detail {

template <typename Traits>
using nanoflann_kd_tree_type = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L_inf_Adaptor<double, KdtreeNanoflann<Traits>>,
    KdtreeNanoflann<Traits>, Traits::dimension, int>;
}

/// \brief Implements neighbourhood searching using a bucket search algorithm,
/// dividing the domain into constant size "buckets".
///
/// This class implements neighbourhood searching using a bucket search
/// algorithm. The domain is first divided up into a regular grid of constant
/// size "buckets", either by using the class constructor to initialise the
/// domain extents, bucket size etc., or by using the reset() member function to
/// reset these parameters.
///
/// After the buckets are created, a set of 3D points can be assigned to their
/// respective buckets using the embed_points() member function. After this,
/// neighbourhood queries around a given point can be performed using
/// find_broadphase_neighbours(), which returns a const iterator to all the
/// points in the same bucket or surrounding buckets of the given point.
///
template <typename Traits>
class KdtreeNanoflann
    : public neighbour_search_base<KdtreeNanoflann<Traits>, Traits,
                                   KdtreeNanoflannQuery<Traits>> {

  typedef typename Traits::double_d double_d;
  typedef typename Traits::position position;
  typedef typename Traits::vector_int vector_int;
  typedef typename Traits::iterator iterator;
  typedef typename Traits::unsigned_int_d unsigned_int_d;
  static const unsigned int dimension = Traits::dimension;

  typedef neighbour_search_base<KdtreeNanoflann<Traits>, Traits,
                                KdtreeNanoflannQuery<Traits>>
      base_type;
  friend base_type;

  typedef detail::nanoflann_kd_tree_type<Traits> kd_tree_type;
  typedef typename kd_tree_type::Node node_type;

public:
  KdtreeNanoflann() : base_type(), m_kd_tree(dimension, *this) {
    // init an empty tree
    std::vector<int> empty;
    m_kd_tree.buildIndex(empty.begin());
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

  //~KdtreeNanoflann() {}

  static constexpr bool ordered() { return true; }

  // Must return the number of data points
  inline size_t kdtree_get_point_count() const {
    return this->m_alive_indices.size();
  }

  // Returns the distance between the vector "p1[0:size-1]"
  // and the data point with index "idx_p2" stored in the class:
  inline double kdtree_distance(const double *p1, const size_t idx_p2,
                                size_t /*size*/) const {
    size_t ret = 0;
    const double_d &p2 = *(get<position>(this->m_particles_begin) + idx_p2);
    for (size_t i = 0; i < dimension; ++i) {
      ret += (p1[i] - p2[i]) * (p1[i] - p2[i]);
    }
    return ret;
  }

  // Returns the dim'th component of the idx'th point in the class:
  // Since this is inlined and the "dim" argument is typically an immediate
  // value, the
  //  "if/else's" are actually solved at compile time.
  inline double kdtree_get_pt(const size_t idx, int dim) const {
    const double_d &p = *(get<position>(this->m_particles_begin) + idx);
    return p[dim];
  }

  // Optional bounding-box computation: return false to default to a standard
  // bbox computation loop.
  //   Return true if the BBOX was already computed by the class and returned in
  //   "bb" so it can be avoided to redo it again. Look at bb.size() to find out
  //   the expected dimensionality (e.g. 2 or 3 for point clouds)
  template <class BBOX> bool kdtree_get_bbox(BBOX &bb) const {
    for (size_t i = 0; i < dimension; ++i) {
      bb[i].low = this->m_bounds.bmin[i];
      bb[i].high = this->m_bounds.bmax[i];
    }
    return true;
  }
  void print_data_structure() const { print_tree(m_kd_tree.get_root_node()); }

private:
  void set_domain_impl() {
    m_kd_tree.set_leaf_max_size(this->m_n_particles_in_leaf);

    this->m_query.m_bounds.bmin = this->m_bounds.bmin;
    this->m_query.m_bounds.bmax = this->m_bounds.bmax;
    this->m_query.m_periodic = this->m_periodic;
  }

  void end_list_of_copies_impl() {}

  void update_iterator_impl() {}

  void print_tree(const node_type *nodes) const {
#ifndef __CUDA_ARCH__
    if (4 <= ABORIA_LOG_LEVEL) {
      std::vector<const node_type *> new_nodes;
      new_nodes.push_back(nodes);
      print_level(new_nodes);
    }
#endif
  }

  void print_level(std::vector<const node_type *> &nodes) const {
#ifndef __CUDA_ARCH__
    if (4 <= ABORIA_LOG_LEVEL) {
      std::vector<const node_type *> new_nodes;
      LOG(4, "printing level with " << nodes.size() << " nodes");
      for (const node_type *ptr : nodes) {
        if (this->m_query.is_leaf_node(*ptr)) {
          const int idx = this->m_query.get_dimension_index(*ptr);
          const int start_index = ptr->node_type.lr.left;
          const int end_index = ptr->node_type.lr.right;
          LOG(4, "\tleaf node with idx = " << idx
                                           << " start index = " << start_index
                                           << " end index = " << end_index);
          LOG(4, "\tparticles in bucket are:");
          for (int i = start_index; i < end_index; ++i) {
            const double_d p = get<position>(this->m_particles_begin)[i];
            LOG(4, "\t\tposition = " << p);
          }
        } else {
          const int idx = this->m_query.get_dimension_index(*ptr);
          const double cut_low = this->m_query.get_cut_low(*ptr);
          const double cut_high = this->m_query.get_cut_high(*ptr);
          LOG(4, "\tNOT leaf node with idx = " << idx << " cut_low= " << cut_low
                                               << " cut_high = " << cut_high);
          const node_type *child1 = this->m_query.get_child1(ptr);
          const node_type *child2 = this->m_query.get_child2(ptr);
          new_nodes.push_back(child1);
          new_nodes.push_back(child2);
        }
      }
      if (new_nodes.size() > 0) {
        print_level(new_nodes);
      } else {
        LOG(4, "finished tree");
      }
    }
#endif
  }

  void update_positions_impl(iterator update_begin, iterator update_end,
                             const int new_n,
                             const bool call_set_domain = true) {
    ASSERT(update_begin == this->m_particles_begin &&
               update_end == this->m_particles_end,
           "error should be update all");

    m_kd_tree.buildIndex(this->m_alive_indices.begin());

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

  /*
  bool add_points_at_end_impl(const size_t dist) {
      return embed_points_impl();
  }

  bool delete_points_impl(const size_t start_index, const size_t n) {
      return embed_points_impl();
  }

  void copy_points_impl(iterator copy_from_iterator, iterator copy_to_iterator)
  { embed_points_impl();
  }
  */

  const KdtreeNanoflannQuery<Traits> &get_query_impl() const { return m_query; }

  KdtreeNanoflannQuery<Traits> &get_query_impl() { return m_query; }

  kd_tree_type m_kd_tree;
  KdtreeNanoflannQuery<Traits> m_query;
};

template <typename Traits> class nanoflann_child_iterator {
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

  nanoflann_child_iterator() : m_high(2), m_index(nullptr) {}

  nanoflann_child_iterator(pointer start, const box_type &bounds)
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

  nanoflann_child_iterator &operator++() {
    increment();
    return *this;
  }

  nanoflann_child_iterator operator++(int) {
    nanoflann_child_iterator tmp(*this);
    operator++();
    return tmp;
  }

  inline bool operator==(const nanoflann_child_iterator &rhs) const {
    return equal(rhs);
  }

  inline bool operator!=(const nanoflann_child_iterator &rhs) const {
    return !operator==(rhs);
  }

  inline bool operator==(const bool rhs) const { return equal(rhs); }

  inline bool operator!=(const bool rhs) const { return !operator==(rhs); }

private:
  bool equal(nanoflann_child_iterator const &other) const {
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
/// @brief This is a query object for the @ref KdtreeNanoflann spatial data
/// structure
///
template <typename Traits> struct KdtreeNanoflannQuery {
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
  template <int LNormNumber, typename Transform = IdentityTransform>
  using query_iterator =
      tree_query_iterator<KdtreeNanoflannQuery, LNormNumber, Transform>;
  typedef value_type *root_iterator;
  typedef depth_first_iterator<KdtreeNanoflannQuery> all_iterator;
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

  ///
  /// @copydoc NeighbourQueryBase::num_children() const
  ///
  size_t num_children() const {
    if (m_root == nullptr) { // empty tree, return a false child iterator
      return 0;
    } else if (is_leaf_node(*m_root)) { // leaf root, return a child iterator
      return 1;
    } else {
      return 2;
    }
  }

  size_t num_children(const child_iterator &ci) const {
    if (is_leaf_node(*ci)) {
      return 0;
    } else {
      return 2;
    }
  }

  static const box_type get_bounds(const child_iterator &ci) {
    return ci.get_bounds();
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

  child_iterator get_bucket(const double_d &position) const {
    child_iterator i = get_children();
    i.go_to(position);

    while (!is_leaf_node(*i)) {
      i = get_children(i);
      i.go_to(position);
    }

    return i;
  }

  size_t get_bucket_index(reference bucket) const { return bucket.index; }

  size_t number_of_buckets() const { return m_number_of_buckets; }

  template <int LNormNumber, typename Transform = IdentityTransform>
  query_iterator<LNormNumber, Transform>
  get_buckets_near_point(const double_d &position, const double max_distance,
                         const Transform &transform = Transform()) const {
#ifndef __CUDA_ARCH__
    LOG(4, "\tget_buckets_near_point: position = "
               << position << " max_distance= " << max_distance);
#endif
    return query_iterator<LNormNumber, Transform>(
        get_children(), position, max_distance, m_number_of_levels, this,
        transform);
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
  iterator_range<theta_iterator> get_theta_buckets(const reference bucket) const
  { return iterator_range<theta_iterator>( theta_iterator(m_root,bucket),
              theta_iterator()
              );
  }
  */
};

} // namespace Aboria

#endif /* NANOFLANN_ADAPTOR_H_ */
