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

template <typename Traits> struct KdtreeQuery;

/// \brief KdTree
///
template <typename Traits>
class Kdtree : public neighbour_search_base<Kdtree<Traits>, Traits,
                                            KdtreeQuery<Traits>> {

  typedef typename Traits::double_d double_d;
  typedef typename Traits::position position;
  typedef typename Traits::vector_int vector_int;
  typedef typename Traits::vector_int2 vector_int2;
  typedef typename Traits::vector_double vector_double;
  typedef typename Traits::vector_double_d vector_double_d;
  typedef typename Traits::iterator iterator;
  typedef typename Traits::unsigned_int_d unsigned_int_d;
  static const unsigned int dimension = Traits::dimension;

  typedef neighbour_search_base<Kdtree<Traits>, Traits, KdtreeQuery<Traits>>
      base_type;
  friend base_type;

public:
  Kdtree() : base_type() {}

  //~Kdtree() {}

  static constexpr bool ordered() { return true; }

  void print_data_structure() const { print_tree(); }

private:
  void set_domain_impl() {
    this->m_query.m_bounds.bmin = this->m_bounds.bmin;
    this->m_query.m_bounds.bmax = this->m_bounds.bmax;
    this->m_query.m_periodic = this->m_periodic;
  }

  void update_iterator_impl() {}

  void print_tree() {
    int first_index = 0;
    int i = first_index;
    while (i < 0 && i < m_nodes_child.size())
      ++i;
    while (i < m_nodes_child.size()) {
      // scan level to find first node
      int level_size = m_nodes_child[i] - first_index;
      print_level(first_index, level_size);
      first_index = first_index + level_size;
      int i = first_index;
      while (i < 0 && i < m_nodes_child.size())
        ++i;
    }
  }

  void print_level(int index, int level_size) const {
    for (size_t i = 0; i < level_size; ++i) {
      if (m_nodes_child[i] >= 0) {
        std::cout << "|n(" << m_nodes_split_dim[i] << ","
                  << m_nodes_split_pos[i] << ")";
      } else {
        std::cout << "|l(" << -m_nodes_child[i] << "," << -m_nodes_split_dim[i]
                  << ")";
      }
    }
    std::cout << "|" << std::endl;
  }

  void print_particles() {
    for (size_t i = 0; i < m_particle_indicies.size() / dimension; ++i) {
      if (i == 0 || m_particle_node[i] != m_particle_node[i - 1]) {
        std::cout << "(l" << m_particle_node[i] << ": " << std::endl;
      }
      std::cout << this->m_particles_begin[m_particle_indicies[i]] << std::endl;
    }
  }

  void update_positions_impl(iterator update_begin, iterator update_end,
                             const int new_n,
                             const bool call_set_domain = true) {
    ASSERT(update_begin == this->m_particles_begin &&
               update_end == this->m_particles_end,
           "error should be update all");

    const size_t num_points = this->m_alive_indices.size();

    // setup particles
    m_particle_node.resize(dimension * num_points, 0);
    m_particle_indicies.resize(dimension * num_points);
    for (size_t i = 0; i < dimension; ++i) {
      // copy particle indicies that are alive
      detail::copy(this->m_alive_indices.begin(), this->m_alive_indices.end(),
                   m_particle_indicies.begin() + i * num_points);

      // sort indicies by position in dimension i
      detail::sort(
          m_particle_indicies.begin() + i * num_points,
          m_particle_indicies.begin() + (i + 1) * num_points,
          [_p = iterator_to_raw_pointer(get<position>(this->m_particles_begin)),
           _i = i](const int a, const int b) { return _p[a][_i] < _p[b][_i]; });
    }

    // build the tree
    build_tree();

    // copy sorted indicies from 1st dim back to m_alive_indicies
    detail::copy(m_particle_indicies.begin(),
                 m_particle_indicies.begin() + num_points,
                 this->m_alive_indicies.begin());

    this->m_query.m_nodes_child =
        iterator_to_raw_pointer(m_nodes_child.begin());
    this->m_query.m_nodes_split_dim =
        iterator_to_raw_pointer(m_nodes_split_dim.begin());
    this->m_query.m_nodes_split_pos =
        iterator_to_raw_pointer(m_nodes_split_pos.begin());
    this->m_query.m_number_of_buckets = m_nodes_child.size();
    this->m_query.m_number_of_levels = m_number_of_levels;
  }

  const KdtreeQuery<Traits> &get_query_impl() const { return m_query; }

  KdtreeQuery<Traits> &get_query_impl() { return m_query; }

private:
  void build_tree() {
    const size_t num_points = this->m_alive_indices.size();

    // setup nodes
    vector_int2 parents_leaf(1, vint2(0, num_points));
    vector_double_d parents_bmin(1, this->m_bounds.bmin);
    vector_double_d parents_bmax(1, this->m_bounds.bmax);

    vector_int2 children_leaf(1, vint2(0, num_points));
    vector_double_d children_bmin(1, this->m_bounds.bmin);
    vector_double_d children_bmax(1, this->m_bounds.bmax);
    vector_int children_type(1, 1);
    vector_int children_num_nodes(1, 0);

    vector_int particle_indicies2(m_particle_indicies.size());

    // setup tree
    m_nodes_child.clear();
    m_nodes_split_pos.clear();
    m_nodes_split_dim.clear();
    m_number_of_levels = 0;
    // m_nodes_child  int-> index of first child node (<0 is a leaf, gives index
    // of 1st particle)
    //
    // m_nodes_split double-> (node) location of split (leaf) not used
    //
    // m_nodes_dim   int-> (node) dimension of split (leaf) -index of
    // last+1 particle
    // m_leafs int2 -> indicies of particles
    while (!parents_leaf.empty()) {
      // calc split ( median of longest dimension)
      const int num_parents = parents_leaf.size();
      auto parents_it = detail::make_zip_iterator(detail::make_tuple(
          parents_leaf.begin(), parents_bmin.begin(), parents_bmax.begin()));
      vector_double parents_split_pos(num_parents);
      vector_int parents_split_dim(num_parents);
      detail::transform(
          parents_it, parents_it + parents_leaf.size(),
          detail::make_zip_iterator(detail::make_tuple(
              parents_split_pos.begin(), parents_split_dim.begin())),
          [_particle_indicies =
               iterator_to_raw_pointer(m_particle_indicies.begin()),
           _p = iterator_to_raw_pointer(
               get<position>(this->m_particles_begin))](const auto &i) {
            const vint2 leaf = detail::get_impl<0>(i);
            const double_d minb = detail::get_impl<1>(i);
            const double_d maxb = detail::get_impl<2>(i);
            double_d minp;
            double_d maxp;

            double max_span = -1;
            for (size_t i = 0; i < dimension; ++i) {
              const double span = maxb[i] - minb[i];
              if (span > max_span) {
                max_span = span;
              }
            }
            double max_spread = -1;
            int split_d = 0;
            for (size_t i = 0; i < dimension; ++i) {
              const double span = maxb[i] - minb[i];
              if (span >=
                  (1.0 - std::numeric_limits<double>::epsilon()) * max_span) {
                const int min_index =
                    _particle_indicies[i * num_points + leaf[0]];
                const int max_index =
                    _particle_indicies[i * num_points + leaf[1]];
                minp[i] = _p[min_index][i];
                maxp[i] = _p[max_index][i];
                const double spread = maxp[i] - minp[i];
                if (spread > max_spread) {
                  max_spread = spread;
                  split_d = i;
                }
              }
            }
            double split = 0.5 * (maxb[split_d] + minb[split_d]);
            if (split < minp[split_d]) {
              split = minp[split_d];
            } else if (split > maxp[split_d]) {
              split = maxp[split_d];
            }
            return detail::make_tuple(split, split_d);
          });

      // partition particles by split
      // TODO: refactor this to a separate algorithm in detail
      //
      //(Sengupta, S., Harris, M., Zhang, Y., & Owens, J. D. (2007). Scan
      // primitives for GPU computing. Graphics …, 97–106.
      // http://doi.org/10.2312/EGGH/EGGH07/097-106)
      vector_int e(m_particle_node.size());
      vector_int f(m_particle_node.size());
      vector_int addr(m_particle_node.size());
      // [t f t f f t f] #in
      // [0 1 0 1 1 0 1] #e = set 1 in false elts.
      detail::transform(
          detail::make_zip_iterator(detail::make_tuple(
              m_particle_indicies.begin(), m_particle_node.begin())),
          detail::make_zip_iterator(detail::make_tuple(
              m_particle_indicies.end(), m_particle_node.end())),
          e.begin(),
          [_p = iterator_to_raw_pointer(get<position>(this->m_particles_begin)),
           _split = iterator_to_raw_pointer(parents_split_pos.begin()),
           _split_d = iterator_to_raw_pointer(parents_split_dim.begin())](
              const auto &i) {
            const int node_index = get<1>(i);
            if (node_index >= 0) {
              const double &p = _p[get<0>(i)][_split_d[node_index]];
              const double &s = _split[node_index];
              return static_cast<int>(p < s);
            } else {
              return 1; // must be 1 so we dont have to count total false
            }
          });

      // [0 0 1 1 2 3 3] #f = enumerate with false=1
      if (num_parents == 1) {
        for (size_t i = 0; i < dimension; ++i) {
          detail::exclusive_scan_by_key(
              m_particle_node.begin() + i * num_points,
              m_particle_node.begin() + (i + 1) * num_points,
              e.begin() + i * num_points, f.begin() + i * num_points);
        }
      } else {
        detail::exclusive_scan_by_key(m_particle_node.begin(),
                                      m_particle_node.end(), e.begin(),
                                      f.begin());
      }

      // [4 4 4 4 4 4 4] # add two last elts. in e, f and scan copy to
      // segment
      //                 # ==total # of falses
      vector_int parents_nf(num_parents);
      detail::tabulate(
          parents_nf.begin, parents_nf.end(),
          [_e = iterator_to_raw_pointer(e.begin()),
           _f = iterator_to_raw_pointer(f.begin()),
           _parents_leaf =
               iterator_to_raw_pointer(parents_leaf.begin())](const int i) {
            const int last_particle_index = _parents_leaf[i][1];
            return _e[last_particle_index] + _f[last_particle_index];
          });

      //[0 1 2 3 4 5 6]  # each thread knows its id
      //[4 5 5 6 6 6 7]  # t = id - f + NF
      //[4 0 5 1 2 6 3]  # addr = e ? f : t
      detail::tabulate(
          addr.begin(), addr.end(),
          [_addr = iterator_to_raw_pointer(addr.begin()),
           _e = iterator_to_raw_pointer(e.begin()),
           _f = iterator_to_raw_pointer(f.begin()),
           _parents_nf = iterator_to_raw_pointer(parents_nf.begin()),
           _particle_node =
               iterator_to_raw_pointer(m_particle_node.begin())](const int i) {
            const int node_index = _particle_node[i];
            _addr[i] = _e[i] ? _f : i - _f[i] + _parents_nf[node_index];
          });

      //[f f f f t t t]  # out[addr] = in (scatter)
      detail::scatter(m_particle_indicies.begin(), m_particle_indicies.end(),
                      addr.begin(), particle_indicies2.begin());
      m_particle_indicies.swap(particle_indicies2);

      // scatter e to new positions (used to update particle_node)
      detail::scatter(e.begin(), e.end(), addr.begin(),
                      particle_indicies2.begin());
      e.swap(particle_indicies2);

      // split the nodes on this level and classify them
      const int num_children = 2 * num_parents;
      children_leaf.resize(num_children);
      children_bmin.resize(num_children);
      children_bmax.resize(num_children);
      children_type.resize(num_children);
      children_num_nodes.resize(num_children);

      detail::tabulate(
          detail::make_zip_iterator(
              detail::make_tuple(children_leaf.begin(), children_bmin.begin(),
                                 children_bmax.begin(), children_type.begin())),
          detail::make_zip_iterator(
              detail::make_tuple(children_leaf.end(), children_bmin.end(),
                                 children_bmax.end(), children_type.end())),
          [_parents_leaf = iterator_to_raw_pointer(parents_leaf.begin()),
           _parents_split_dim =
               iterator_to_raw_pointer(parents_split_dim.begin()),
           _parents_split_pos =
               iterator_to_raw_pointer(parents_split_pos.begin()),
           _parents_bmin = iterator_to_raw_pointer(parents_bmin.begin()),
           _parents_bmax = iterator_to_raw_pointer(parents_bmax.begin()),
           _threshold = this->m_n_particles_in_leaf](const int i) {
            const int node_index = i / 2;
            const int right_node = i % s;
            const int2 leaf = _parents_leaf[node_index];
            const int split_d = _parents_split_dim[node_index];
            const double split = _parents_split_pos[node_index];
            const int nf = _parents_nf[node_index];
            const double bmin = _parents_bmin[node_index];
            const double bmax = _parents_bmax[node_index];
            const int n = leaf[1] - leaf[0];
            if (right_node) {
              leaf[0] = leaf[1] - (n - nf);
              bmin[splid_d] = split;
            } else {
              leaf[1] = leaf[0] + nf;
              bmax[splid_d] = split;
            }
            const int is_node =
                static_cast<int>(leaf[1] - leaf[0] > _threshold);
            return detail::make_tuple(leaf, bmin, bmax, is_node);
          });

      // enumerate nodes and leafs
      detail::exclusive_scan(children_type.begin(), children_type.end(),
                             children_num_nodes.begin());

      // update mask with children indicies
      detail::transform(
          detail::make_zip_iterator(
              detail::make_tuple(particle_node.begin(), e.begin())),
          detail::make_zip_iterator(
              detail::make_tuple(particle_node.end(), e.end())),
          particle_node.begin(),
          [_children_num_nodes =
               iterator_to_raw_pointer(children_num_nodes.begin()),
           _children_type =
               iterator_to_raw_pointer(children_type.begin())](const auto &i) {
            const int parent_index = get<0>(i);
            if (parent_index >= 0) {
              const int child_index = 2 * parent_index + get<1>(i);
              return _children_type[child_index]
                         ? _children_num_nodes[child_index]
                         : -1;
            } else {
              return parent_index;
            }
          });

      // add children to tree
      ++m_number_of_levels;
      const int total_children_nodes =
          *(children_num_nodes.end() - 1) + *children_type.end();
      const int children_begin = m_nodes_child.size();
      const int children_end = children_begin + num_children;
      m_nodes_child.resize(children_end);
      auto tree_it = detail::make_zip_iterator(
          detail::make_tuple(m_nodes_child.begin(), m_nodes_split_pos.begin(),
                             m_nodes_split_dim.begin()));
      detail::tabulate(
          tree_it + children_begin, tree_it + children_end,
          [_children_num_nodes =
               iterator_to_raw_pointer(children_num_nodes.begin()),
           _children_type = iterator_to_raw_pointer(children_type.begin()),
           _children_leaf = iterator_to_raw_pointer(children_leaf.begin()),
           _parents_split_dim =
               iterator_to_raw_pointer(parents_split_dim.begin()),
           _parents_split_pos = iterator_to_raw_pointer(
               parents_split_pos.begin())](const int i) {
            const int parent_index = i / 2;
            const int next_level_index = _children_num_nodes[i];
            const int is_node = _children_type[i];
            const vint2 leaf = _children_leaf[i];
            const int split_dim = _parents_split_dim[parent_index];

            const int child = is_node ? children_end + 2 * get<0>(i) : -leaf[0];
            const double split = _parents_split_pos[parent_index];
            const int dim = is_node ? split_dim : -leaf[1];
            return detail::make_tuple(child, split, dim);
          });

#ifndef __CUDA_ARCH__
      if (3 <= ABORIA_LOG_LEVEL) {
        print_level(children_begin, num_children);
      }
#endif

      // setup new parent nodes (don't copy leafs)
      parent_leaf.resize(total_children_nodes);
      parent_bmin.resize(total_children_nodes);
      parent_bmax.resize(total_children_nodes);
      auto child_it = detail::make_zip_iterator(detail::make_tuple(
          children_leaf.begin(), children_bmin.begin(), children_bmax.begin()));
      detail::copy_if(child_it, child_it + num_children, children_type.begin(),
                      parent_it);
    }
  }

  vector_int m_nodes_child;
  vector_int m_nodes_split_dim;
  vector_double m_nodes_split_pos;

  vector_int m_particle_indicies;
  vector_int m_particle_node;
  int m_number_of_levels;
  KdtreeQuery<Traits> m_query;
}; // namespace Aboria

template <typename Query> class KdtreeChildIterator {
  typedef bbox<Query::dimension> box_type;

  int m_high;
  const int *m_index;
  box_type m_bounds;

public:
  typedef const int *pointer;
  typedef std::forward_iterator_tag iterator_category;
  typedef const int value_type;
  typedef const int &reference;
  typedef std::ptrdiff_t difference_type;

  KdtreeChildIterator() : m_high(2), m_index(-1) {}

  KdtreeChildIterator(const int start, const box_type &bounds,
                      const Query *query)
      : m_high(0), m_index(start), m_bounds(bounds), m_query(query) {
    ASSERT(start != nullptr, "start pointer should not be null");
  }

  bool is_high() const { return m_high > 0; }

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

  reference dereference() const { return *m_index; }

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
  typedef const int &reference;
  typedef const int *pointer;

  typedef Traits traits_type;
  typedef typename Traits::raw_pointer raw_pointer;
  typedef typename Traits::double_d double_d;
  typedef typename Traits::bool_d bool_d;
  typedef typename Traits::int_d int_d;
  typedef typename Traits::unsigned_int_d unsigned_int_d;
  template <int LNormNumber>
  using query_iterator = tree_query_iterator<KdtreeQuery, LNormNumber>;
  typedef pointer root_iterator;
  typedef depth_first_iterator<KdtreeQuery> all_iterator;
  typedef KdtreeChildIterator<KdtreeQuery> child_iterator;
  typedef ranges_iterator<Traits> particle_iterator;
  typedef bbox<dimension> box_type;

  bool_d m_periodic;
  bbox<dimension> m_bounds;
  raw_pointer m_particles_begin;
  raw_pointer m_particles_end;
  size_t m_number_of_buckets;
  size_t m_number_of_levels;

  int *m_nodes_child;
  int *m_nodes_split_dim;
  double *m_nodes_split_pos;

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
  static bool is_leaf_node(reference bucket) { return bucket < 0; }
  static bool is_tree() { return true; }

  /*
   * end functions for tree_query_iterator
   */

  child_iterator get_children() const {
    return child_iterator(0, m_bounds, this);
  }

  static child_iterator get_children(const child_iterator &ci) {
    if (!is_leaf_node(*ci)) {
      return child_iterator((m_nodes_child + *ci + ci.is_high()),
                            get_bounds(ci));
    } else {
      return child_iterator();
    }
  }

  static const box_type get_bounds(const child_iterator &ci) {
    box_type ret = ci.m_bounds;
    const int i = m_node_split_dim[ci.m_index - m_nodes_begin];
    if (ci.is_high()) {
      ret.bmin[i] = m_node_split_pos[ci.m_index - m_nodes_begin];
    } else {
      ret.bmax[i] = m_node_split_pos[ci.m_index - m_nodes_begin];
    }
    return ret;
  }

  template <int LNormNumber>
  static const box_type get_bounds(const query_iterator<LNormNumber> &qi) {
    return qi.get_bounds();
  }

  static const box_type get_bounds(const all_iterator &ai) {
    return ai.get_bounds();
  }

  particle_iterator get_bucket_particles(reference bucket) const {
#ifndef __CUDA_ARCH__
    LOG(4, "\tget_bucket_particles: looking in bucket with idx = " << bucket);
#endif
    if (!is_leaf_node(bucket)) {
      return particle_iterator();
    }

    return particle_iterator(m_particles_begin - bucket,
                             m_particles_begin -
                                 m_nodes_split_dim[&bucket - m_nodes_child]);
  }

  void go_to(const double_d &position, child_iterator &ci) {
    ASSERT(position[i] < ci.m_bounds.bmax[i], "position out of bounds");
    ASSERT(position[i] >= ci.m_bounds.bmin[i], "position out of bounds");
    const int i = m_node_split_dim[ci.m_index - m_nodes_child];
    const double diff =
        position[i] - m_query->m_node_split_pos[ci.m_index - m_nodes_child];
    if (diff > 0)
      ++ci;
  }

  void get_bucket(const double_d &position, pointer &bucket,
                  box_type &bounds) const {
    child_iterator i = get_children();
    go_to(position, ci);

    while (!is_leaf_node(*i)) {
      i = get_children(i);
      go_to(position, ci);
    }

    bucket = &(*i);
    bounds = get_bounds(i);
  }

  size_t get_bucket_index(reference bucket) const {
    return ci.m_index - m_nodes_child;
  }

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
};

} // namespace Aboria

#endif /* KDTREE_H_ */
