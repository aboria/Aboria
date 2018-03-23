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
  void build_tree() {
    const size_t num_points = this->m_alive_indices.size();

    // setup nodes
    vector_int2     parents_leaf(1, vint2(0,num_points);
    vector_double_d parents_bmin(1, this->m_bounds.bmin);
    vector_double_d parents_bmax(1, this->m_bounds.bmax);

    // setup particles
    vector_int particle_node(num_points, 0);
    vector_int particle_indicies(D*num_points);
    for (size_t i = 0; i < D; ++i) {
    detail::copy(this->m_alive_indices.begin(),this->m_alive_indices.end(),particle_indicies.begin + i*num_points);
    }

    // setup tree
    m_nodes.clear();
    m_nodes_split.clear();
    m_nodes_dim.clear();
    // m_nodes_child  int-> index of first child node (<0 is a leaf, gives index
    // of 1st particle)
    //
    // m_nodes_split double-> (node) location of split (leaf) not used
    //
    // m_nodes_dim   int-> (node) dimension of split (leaf) -index of
    // last+1 particle
    // m_leafs int2 -> indicies of particles
    while (!parents_leaf.empty()) {
      // make a handy transform iterator that gets the particle positions
      /*
      auto positions_begin = detail::make_transform_iterator(
          particle_indicies.begin(),
          [_p = iterator_to_raw_pointer(get<position>(this->m_particles_begin)),
           _alive_indicies =
               iterator_to_raw_pointer(this->m_alive_indices.begin())](
              const int i) { return _p[_alive_indicies[i]]; });

      // get min bounds (TODO: can combine?)
      double_d_vector parents_pmin(num_parents);
      const int num_parents = parents_leaf.size();
      vector_int parents(num_parents);
      detail::reduce_by_key(
          node_mask[i].begin(), node_mask[i].end(), positions_begin,
          parents.begin(), parents_pmin.begin(),
          detail::equal_to<int>(),
          [_p = iterator_to_raw_pointer(get<position>(this->m_particles_begin)),
           _alive_indicies =
               iterator_to_raw_pointer(this->m_alive_indices.begin())](
              const double_d &a, const double_d &b) { return a.min(b); });

      // get max bounds
      double_d_vector parents_pmax(num_parents);
      detail::reduce_by_key(
          node_mask[i].begin(), node_mask[i].end(), positions_begin,
          parents.begin(), parents_pmax.begin(),
          detail::equal_to<int>(),
          [](const double_d &a, const double_d &b) { return a.max(b); });
          */

      // calc split ( median of longest dimension)
      vector_double parents_split_pos(num_parents);
      vector_int    parents_split_dim(num_parents);
      detail::tabulate(
              detail::make_zip_iterator(detail::make_tuple(
          parents_split_pos.begin(), parents_split_dim.begin())),
                detail::make_zip_iterator(detail::make_tuple(
          parents_split_pos.end(), parents_split_dim.end())),
                       [](const int i) {
                         const double_d &minp = get<0>(tpl);
                         const double_d &maxp = get<1>(tpl);
                         const double_d &minb = get<2>(tpl);
                         const double_d &maxb = get<3>(tpl);
                         const int &split_d = get<4>(tpl);
                         const int &split = get<5>(tpl);
                         double max_span = -1;
                         for (size_t i = 0; i < D; ++i) {
                           const double span = maxb[i] - minb[i];
                           if (span > max_span) {
                             max_span = span;
                           }
                         }
                         double max_spread = -1;
                         int cut_d = 0;
                         for (size_t i = 0; i < D; ++i) {
                           const double span = maxb[i] - minb[i];
                           if (span >= (1 - EPS) * max_span) {
                             const double spread = maxp[i] - minp[i];
                             if (spread > max_spread) {
                               max_spread = spread;
                               cut_d = i;
                             }
                           }
                         }
                         const double cut = 0.5 * (bmax[cut_d] + bmin[cut_d]);
                         if (cut < minp[cut_d]) {
                           split = minp[cut_d];
                         } else if (cut > maxp[cut_d]) {
                           split = maxp[cut_d];
                         } else {
                           split = cut;
                         }
                         split_d = cut_d;
                       });

      // partition m_sorted_indices by split
      // TODO: refactor this to a separate algorithm in detail
      //
      //(Sengupta, S., Harris, M., Zhang, Y., & Owens, J. D. (2007). Scan
      // primitives for GPU computing. Graphics …, 97–106.
      // http://doi.org/10.2312/EGGH/EGGH07/097-106)
      vector_int e(node_mask.size());
      vector_int f(node_mask.size());
      vector_int addr(node_mask.size());
      // [t f t f f t f] #in
      // [0 1 0 1 1 0 1] #e = set 1 in false elts.
      detail::transform(
          detail::make_zip_iterator(
              detail::make_tuple(particle_indicies.begin(), node_mask.begin())),
          detail::make_zip_iterator(
              detail::make_tuple(particle_indicies.end(), node_mask.end())),
          e.begin(),
          [_p = iterator_to_raw_pointer(get<position>(this->m_particles_begin)),
           _split = iterator_to_raw_pointer(split),
           _split_d = iterator_to_raw_pointer(split_d)](const auto &i_m) {
            const int &node_index = get<1>(i_m);
            const double &p = _p[get<0>(i_m)][_split_d[node_index]];
            const double &s = _split[node_index];
            return static_cast<int>(p < s);
          });

      // [0 0 1 1 2 3 3] #f = enumerate with false=1
      detail::exclusive_scan_by_key(node_mask[i].begin(), node_mask[i].end(),
                                    e.begin(), f.begin());
      // [4 4 4 4 4 4 4] # add two last elts. in e, f and scan copy to segment
      //                 # ==total # of falses
      vector_int parents_nf(num_parents);
      detail::tabulate(parents_nf.begin, parents_nf.end(),
                       [_e = iterator_to_raw_pointer(e.begin()),
                        _f = iterator_to_raw_pointer(f.begin())](const int i) {
                         const int last_particle_index =
                             _parents_leaf[i][1];
                         return e[last_particle_index] + f[last_particle_index];
                       });

      //[0 1 2 3 4 5 6]  # each thread knows its id
      //[4 5 5 6 6 6 7]  # t = id - f + NF
      //[4 0 5 1 2 6 3]  # addr = e ? f : t
      detail::tabulate(
          addr.begin(), addr.end(),
          [_addr = iterator_to_raw_pointer(addr.begin()),
           _e = iterator_to_raw_pointer(e.begin()),
           _f = iterator_to_raw_pointer(f.begin()),
           _num_false =
               iterator_to_raw_pointer(num_false.begin())](const int i) {
            const int node_index = _node_mask[i];
            _addr[i] =
                _e[i] ? _f : i - _f[i] + _parents_nf[node_index];
          });

      //[f f f f t t t]  # out[addr] = in (scatter)
      vector_int particle_indicies_buffer(particle_indicies.size());
      detail::scatter(particle_indicies.begin(), particle_indicies.end(),
                      addr.begin(), particle_indicies_buffer.begin());
      particle_indicies.swap(particle_indicies_buffer);
      // scatter e to new positions (used to update node mask)
      detail::scatter(e.begin(), e.end(), addr.begin(),
                      particle_indicies_buffer.begin());
      e.swap(particle_indicies_buffer);

      // split the nodes on this level and classify them
      const int num_children = 2*num_parents;
      vector_int2 children_leaf(num_children);
      vector_double_d children_bmin(num_children);
      vector_double_d children_bmax(num_children);
      vector_int children_type(next_nodes_leaf.size());

      detail::tabulate(
          detail::make_zip_iterator(detail::make_tuple(
              children_leaf.begin(), children_bmin.begin(),
              children_bmax.begin(), children_type.begin())),
        detail::make_zip_iterator(detail::make_tuple(
              children_leaf.end(), children_bmin.end(),
              children_bmax.end(), children_type.end())),
          [](const int i) {
            const int active_node_index = i / 2 const int right_node = i % s;
            const int2 leaf = _active_nodes_leaf[active_node_index];
            const int split_d = _active_nodes_split_d[active_node_index];
            const double split = _active_nodes_split[active_node_index];
            const int nf = _active_nodes_n_false[active_node_index];
            const double bmin = _active_nodes_bmin[active_node_index];
            const double bmax = _active_nodes_bmax[active_node_index];
            const int n = leaf[1] - leaf[0];
            if (right_node) {
              leaf[0] = leaf[1] - (n - nf);
              bmin[splid_d] = split;
            } else {
              leaf[1] = leaf[0] + nf;
              bmax[splid_d] = split;
            }
            const int is_node = static_cast<int>(leaf[1] - leaf[0] > threshold);
            return detail::make_tuple(leaf, bmin, bmax, is_node);
          });

      // enumerate nodes and leafs
      vector_int children_num_nodes(num_children);
      detail::exclusive_scan(children_type.begin(), children_type.end(), 
                             children_num_nodes.begin());

      // update mask with children indicies
      detail::transform(detail::make_zip_iterator(
                            detail::make_tuple(node_mask.begin(), e.begin())),
                        node_mask.end(), node_mask.begin(), [](const auto &i) {
                          const int child_index = 2 * get<0>(i) + get<1>(i);
                          return _children_num_nodes[child_index];
                        });

      // remove leaf nodes from node_mask, and particle_indicies
      auto new_end =
          detail::remove_if(
              detail::make_zip_iterator(detail::make_tuple(
                  particle_indicies.begin(), node_mask.begin())),
              detail::make_zip_iterator(detail::make_tuple(
                  particle_indicies.end(), node_mask.end())),
              [](const auto &i) {
                return get<1>(i) < 0;
              })
              .get_iterator_tuple();
      particle_indicies.erase(get<0>(new_end), particle_indicies.end());
      node_mask.erase(get<1>(new_end), node_mask.end());

      // add children to tree
      const int total_children_nodes = *(children_num_nodes.end()-1) 
                                        + *children_type.end();
      const int children_begin = m_nodes_child.size();
      const int children_end = children_begin + num_children;
      m_nodes_child.resize(children_end); 
      detail::tabulate(m_nodes_child.begin()+children_end,m_nodes_child.end(),[](const int i) {
              const int node_index = _num_nodes_and_leafs[i][0];
              const int parent_index = static_cast<int>(std::floor(i/2.0));
              const bool is_node = node_index >= 0;
              const vint2 leaf = _next_nodes_leaf[i];
              const int split_dim = _active_nodes_split_d[parent_index];

                const int child = is_node ? children_end + 2*get<0>(i) : -leaf[0];
              const double split = _active_nodes_split[parent_index];
                const int dim = is_node ? split_dim : -leaf[1];
                return detail::make_tuple(child,split,dim);
              });

      // setup parent nodes (don't copy leafs)
      parent_leaf.resize(total_children_nodes);
      parent_bmin.resize(total_children_nodes);
      parent_bmax.resize(total_children_nodes);
      detail::copy_if(
              detail::make_zip_iterator(detail::make_tuple(
                  children_leaf.begin(), children_bmin.begin(), children_bmax.begin())),
                detail::make_zip_iterator(detail::make_tuple(
                  children_leaf.end(), children_bmin.end(), children_bmax.end())),
                detail::make_zip_iterator(detail::make_tuple(
                  parent_leaf.begin(), parent_bmin.begin(), parent_bmax.begin())),
                children_type.begin(),
              [](const int i) {
                return i == 1;
              });
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
