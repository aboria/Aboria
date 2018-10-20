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
// Acknowledgement: This source was modified from the Thrust workshop git
// repository by Jared Hoberock,
// https://github.com/jaredhoberock/thrust-workshop
//

#ifndef OCTTREE_H_
#define OCTTREE_H_

#include "CudaInclude.h"
#include "Particles.h"
#include "SpatialUtil.h"
#include "detail/Algorithms.h"

namespace Aboria {

template <typename Traits> struct HyperOctreeQuery;

///
/// @brief A hyper octree spatial data structure that is paired with a
/// HyperOctreeQuery query type
///
/// An hyper octree is the general n-dimensional analog of an octree (aside: a
/// normal octree is only defined for three dimensions). This data structure is
/// a tree, where each level of the tree splits each bucket along the middle of
/// each dimension. Each new bucket formed by this split becomes a child bucket
/// for that bucket. Once the number of particles within  a box falls below a
/// threshold value (default: 10) then it is a leaf of the tree with no children
///
/// For example, consider a certain level of the tree in two dimensions which
/// has 3 buckets, one a leaf with 9 particles, the other two non-leafs with 44
/// and 11 particles respectivly. The threshold value is set at 10.
/// ~~~
/// -----------------------------------------
/// |                  |                    |
/// |      leaf        |                    |
/// |        9         |        44          |
/// |                  |                    |
/// |                  |                    |
/// |------------------+--------------------|
/// |                  |                    |
/// |                  |                    |
/// |        11        |                    |
/// |                  |                    |
/// |                  |                    |
/// |                  |                    |
/// -----------------------------------------
/// ~~~
/// The next level of the tree might look like this:
/// ~~~
/// -----------------------------------------
/// |                  |   12    |   9      |
/// |                  |         |  leaf    |
/// |                  +---------+----------|
/// |                  |   11    |   12     |
/// |                  |         |          |
/// |--------+---------|---------+----------|
/// |   1    |   8     |                    |
/// | leaf   |  leaf   |                    |
/// |--------+---------|                    |
/// |   0    |   2     |                    |
/// | leaf   |  leaf   |                    |
/// -----------------------------------------
/// ~~~
/// @tparam Traits an instatiation of TraitsCommon
///
template <typename Traits>
class HyperOctree : public neighbour_search_base<HyperOctree<Traits>, Traits,
                                                 HyperOctreeQuery<Traits>> {

  typedef typename Traits::double_d double_d;
  typedef typename Traits::bool_d bool_d;
  typedef typename Traits::position position;
  typedef typename Traits::vector_double_d_const_iterator
      vector_double_d_const_iterator;
  typedef typename Traits::vector_unsigned_int_iterator
      vector_unsigned_int_iterator;
  typedef typename Traits::vector_unsigned_int vector_unsigned_int;
  typedef typename Traits::vector_int vector_int;
  typedef typename Traits::unsigned_int_d unsigned_int_d;
  typedef typename Traits::template vector_type<vint2>::type vector_int2;
  static const unsigned int dimension = Traits::dimension;

  // number of children = 2^d
  static const size_t nchild = (1 << dimension);

  typedef typename Traits::iterator iterator;

  typedef neighbour_search_base<HyperOctree<Traits>, Traits,
                                HyperOctreeQuery<Traits>>
      base_type;

  friend base_type;

public:
  HyperOctree() : base_type(), m_max_level(32 / dimension - 2) {

    // need to init a tree with 1 level (for 0 particles) in case
    // someone does a query on an empty data structure
    build_tree();
    this->m_query.m_number_of_levels = m_number_of_levels;
    this->m_query.m_nodes_begin =
        iterator_to_raw_pointer(this->m_nodes.begin());
    this->m_query.m_leaves_begin =
        iterator_to_raw_pointer(this->m_leaves.begin());
    this->m_query.m_number_of_nodes = m_nodes.size();
  }

  static constexpr bool ordered() { return true; }

  void build_tree();
  struct classify_point;
  struct child_index_to_tag_mask;
  struct classify_node;
  struct write_nodes;
  struct make_leaf;

  void print_data_structure() const {
#ifndef __CUDA_ARCH__
    detail::print_nodes(m_nodes, dimension);
    detail::print_leaves(m_leaves);
#endif
  }

private:
  void set_domain_impl() {
    this->m_query.m_periodic = this->m_periodic;
    this->m_query.m_bounds = this->m_bounds;
  }

  void update_iterator_impl() {}

  void update_positions_impl(iterator update_begin, iterator update_end,
                             const int new_n,
                             const bool call_set_domain = true) {
    ASSERT(update_begin == this->m_particles_begin &&
               update_end == this->m_particles_end,
           "error should be update all");

    const size_t num_points = this->m_alive_indices.size();

    m_tags.resize(num_points);
    if (m_tags.size() > 0) {
      /******************************************
       * 3. Classify points                     *
       ******************************************/
      if (static_cast<size_t>(update_end - update_begin) == num_points) {
        // m_alive_indicies is just a sequential list of indices
        // (i.e. no dead)
        detail::transform(
            get<position>(this->m_particles_begin) + this->m_alive_indices[0],
            get<position>(this->m_particles_begin) + this->m_alive_indices[0] +
                num_points,
            m_tags.begin(), classify_point(this->m_bounds, m_max_level));

      } else {
        // m_alive_indicies contains all alive indicies
        detail::transform(Traits::make_permutation_iterator(
                              get<position>(this->m_particles_begin),
                              this->m_alive_indices.begin()),
                          Traits::make_permutation_iterator(
                              get<position>(this->m_particles_begin),
                              this->m_alive_indices.end()),
                          m_tags.begin(),
                          classify_point(this->m_bounds, m_max_level));
      }

      /******************************************
       * 4. Sort according to classification    *
       ******************************************/
      detail::sort_by_key(m_tags.begin(), m_tags.end(),
                          this->m_alive_indices.begin());
    }

    build_tree();

#ifndef __CUDA_ARCH__
    if (3 <= ABORIA_LOG_LEVEL) {
      detail::print_nodes(m_nodes, dimension);
      detail::print_leaves(m_leaves);
    }
#endif

    this->m_query.m_number_of_levels = m_number_of_levels;
    this->m_query.m_nodes_begin =
        iterator_to_raw_pointer(this->m_nodes.begin());
    this->m_query.m_leaves_begin =
        iterator_to_raw_pointer(this->m_leaves.begin());
    this->m_query.m_number_of_nodes = m_nodes.size();
  }

  /*
  bool add_points_at_end_impl(const size_t dist) {
      const size_t num_points  = this->m_particles_end -
  this->m_particles_begin; auto start_adding_particles =
  this->m_particles_end-dist; m_tags.resize(num_points); auto start_adding_tags
  = this->m_tags.end()-dist;

      //3. Classify new points
      detail::transform(get<position>(start_adding_particles),
              get<position>(this->m_particles_end),
              start_adding_tags,
              classify_point(this->m_bounds, m_max_level));

      // sort and then build tree
      sort_by_tags();
      build_tree();

      this->m_query.m_number_of_levels = m_number_of_levels;
      this->m_query.m_particles_begin =
  iterator_to_raw_pointer(this->m_particles_begin); this->m_query.m_nodes_begin
  = iterator_to_raw_pointer(this->m_nodes.begin());
      this->m_query.m_leaves_begin=
  iterator_to_raw_pointer(this->m_leaves.begin());
      this->m_query.m_number_of_nodes = m_nodes.size();
      this->m_query.m_number_of_particles = this->m_particles_end
                                           -this->m_particles_begin;

      return true;
  }


  bool delete_points_impl(const size_t start_index, const size_t n) {
      const size_t nparticles = this->m_particles_end - this->m_particles_begin;

      m_tags.erase(m_tags.begin()+start_index,m_tags.begin()+start_index+n);
      build_tree();

      this->m_query.m_number_of_levels = m_number_of_levels;
      this->m_query.m_particles_begin =
  iterator_to_raw_pointer(this->m_particles_begin); this->m_query.m_nodes_begin
  = iterator_to_raw_pointer(this->m_nodes.begin());
      this->m_query.m_leaves_begin=
  iterator_to_raw_pointer(this->m_leaves.begin());
      this->m_query.m_number_of_nodes = m_nodes.size();
      this->m_query.m_number_of_particles = this->m_particles_end
                                           -this->m_particles_begin;

      return false;
  }

  void copy_points_impl(iterator copy_from_iterator, iterator copy_to_iterator)
  { auto positions_from = get<position>(copy_from_iterator); auto positions_to =
  get<position>(copy_to_iterator);

      const size_t index_from = copy_from_iterator-this->m_particles_begin;
      const size_t index_to = copy_to_iterator-this->m_particles_begin;

      m_tags[index_to] = m_tags[index_from];
      sort_by_tags();
      build_tree();

      this->m_query.m_number_of_levels = m_number_of_levels;
      this->m_query.m_particles_begin =
  iterator_to_raw_pointer(this->m_particles_begin); this->m_query.m_nodes_begin
  = iterator_to_raw_pointer(this->m_nodes.begin());
      this->m_query.m_leaves_begin=
  iterator_to_raw_pointer(this->m_leaves.begin());
      this->m_query.m_number_of_nodes = m_nodes.size();

  }
  */

  const HyperOctreeQuery<Traits> &get_query_impl() const {
    return this->m_query;
  }

  HyperOctreeQuery<Traits> &get_query_impl() { return m_query; }

  /*
  void sort_by_tags() {
      if (m_tags.size() > 0) {
          //this->m_order.resize(m_tags.size());
          //detail::sequence(this->m_order.begin(), this->m_order.end());
          detail::sort_by_key(m_tags.begin(), m_tags.end(),
                              this->m_alive_indices.begin());
      }
  }
*/

private:
  int max_points;
  int m_max_level;
  unsigned m_number_of_levels;

  vector_int m_tags;
  vector_int m_nodes;
  vector_int2 m_leaves;

  HyperOctreeQuery<Traits> m_query;
};

template <typename Traits> void HyperOctree<Traits>::build_tree() {
  m_nodes.clear();
  m_leaves.clear();
  vector_int active_nodes(1, 0);

  LOG(4, "octree: building tree with max_level = " << m_max_level);

  // Build the tree one level at a time, starting at the root
  for (int level = 1; !active_nodes.empty(); ++level) {

    ASSERT(level <= m_max_level,
           "octtree build_tree: tree has exceeded max levels");

    /******************************************
     * 1. Calculate children                  *
     ******************************************/

    // New children: 2^D quadrants per active node
    vector_int children(nchild * active_nodes.size());

    // For each active node, generate the tag mask for each of its 2^D children
    detail::tabulate(
        children.begin(), children.end(),
        child_index_to_tag_mask(level, m_max_level, active_nodes.data()));

    /******************************************
     * 2. Determine interval for each child   *
     ******************************************/

    // For each child we need interval bounds
    vector_int lower_bounds(children.size());
    vector_int upper_bounds(children.size());

    // Locate lower and upper bounds for points in each quadrant
    detail::lower_bound(m_tags.begin(), m_tags.end(), children.begin(),
                        children.end(), lower_bounds.begin());

    int length = (1 << (m_max_level - level) * dimension) - 1;

    auto plus_length = [=] CUDA_HOST_DEVICE(const int i) { return i + length; };
    detail::upper_bound(
        m_tags.begin(), m_tags.end(),
        Traits::make_transform_iterator(children.begin(), plus_length),
        Traits::make_transform_iterator(children.end(), plus_length),
        upper_bounds.begin());

    /******************************************
     * 3. Mark each child as empty/leaf/node  *
     ******************************************/

    // Mark each child as either empty, a node, or a leaf
    vector_int child_node_kind(children.size(), 0);
    detail::transform(
        Traits::make_zip_iterator(
            Traits::make_tuple(lower_bounds.begin(), upper_bounds.begin())),
        Traits::make_zip_iterator(
            Traits::make_tuple(lower_bounds.end(), upper_bounds.end())),
        child_node_kind.begin(),
        classify_node(this->m_n_particles_in_leaf, level == m_max_level));

    /******************************************
     * 4. Enumerate nodes and leaves          *
     ******************************************/

    // Enumerate the nodes and leaves at this level
    vector_int leaves_on_this_level(child_node_kind.size());
    vector_int nodes_on_this_level(child_node_kind.size());

    // Enumerate nodes at this level
    auto plus = [] CUDA_HOST_DEVICE(const int a, const int b) { return a + b; };
    detail::transform_exclusive_scan(
        child_node_kind.begin(), child_node_kind.end(),
        nodes_on_this_level.begin(), detail::is_a<detail::NODE>(), 0, plus);

    // Enumerate leaves at this level
    detail::transform_exclusive_scan(
        child_node_kind.begin(), child_node_kind.end(),
        leaves_on_this_level.begin(), detail::is_a<detail::LEAF>(), 0, plus);

    int num_nodes_on_this_level =
        nodes_on_this_level.back() +
        (child_node_kind.back() == detail::NODE ? 1 : 0);
    int num_leaves_on_this_level =
        leaves_on_this_level.back() +
        (child_node_kind.back() == detail::LEAF ? 1 : 0);

    /******************************************
     * 5. Add the children to the node list   *
     ******************************************/

    int num_children = child_node_kind.size();

    int children_begin = m_nodes.size();
    m_nodes.resize(m_nodes.size() + num_children);

    detail::transform(Traits::make_zip_iterator(Traits::make_tuple(
                          child_node_kind.begin(), nodes_on_this_level.begin(),
                          leaves_on_this_level.begin())),
                      Traits::make_zip_iterator(Traits::make_tuple(
                          child_node_kind.end(), nodes_on_this_level.end(),
                          leaves_on_this_level.end())),
                      m_nodes.begin() + children_begin,
                      write_nodes(m_nodes.size(), m_leaves.size()));

    /******************************************
     * 6. Add the leaves to the leaf list     *
     ******************************************/

    children_begin = m_leaves.size();

    m_leaves.resize(m_leaves.size() + num_leaves_on_this_level);

    detail::scatter_if(
        Traits::make_transform_iterator(
            Traits::make_zip_iterator(
                Traits::make_tuple(lower_bounds.begin(), upper_bounds.begin())),
            make_leaf()),
        Traits::make_transform_iterator(
            Traits::make_zip_iterator(
                Traits::make_tuple(lower_bounds.end(), upper_bounds.end())),
            make_leaf()),
        leaves_on_this_level.begin(), child_node_kind.begin(),
        m_leaves.begin() + children_begin, detail::is_a<detail::LEAF>());

    /*
    detail::print_active_nodes<dimension>(active_nodes,m_max_level);
    detail::print_children<dimension>(children,m_max_level);
    detail::print_child_bounds(lower_bounds,upper_bounds);
    detail::print_child_node_kind(child_node_kind);
    detail::print_child_enumeration(child_node_kind,nodes_on_this_level,leaves_on_this_level);
    detail::print_nodes(m_nodes,dimension);
    detail::print_leaves(m_leaves);
    */

    /******************************************
     * 7. Set the nodes for the next level    *
     ******************************************/

    // Set active nodes for the next level to be all the childs nodes from this
    // level
    active_nodes.resize(num_nodes_on_this_level);

    detail::copy_if(children.begin(), children.end(), child_node_kind.begin(),
                    active_nodes.begin(), detail::is_a<detail::NODE>());

    m_number_of_levels = level;
  }
}

// Classify a point with respect to the bounding box.
template <typename traits> struct HyperOctree<traits>::classify_point {
  bbox<dimension> box;
  int max_level;

  // Create the classifier
  classify_point(const bbox<dimension> &b, int lvl) : box(b), max_level(lvl) {}

  // Classify a point
  inline CUDA_HOST_DEVICE int operator()(const double_d &p) {
    return detail::point_to_tag(p, box, max_level);
  }
};

template <typename traits> struct HyperOctree<traits>::child_index_to_tag_mask {
  const int level, max_level;

  // mask for lower n bits, where n is the number of dimensions
  const static unsigned mask = nchild - 1;

  typedef typename vector_int::const_pointer ptr_type;
  ptr_type m_nodes;

  child_index_to_tag_mask(int lvl, int max_lvl, ptr_type nodes)
      : level(lvl), max_level(max_lvl), m_nodes(nodes) {}

  inline CUDA_HOST_DEVICE int operator()(int idx) const {
    int tag = m_nodes[idx / nchild];
    int which_child = (idx & mask);
    return detail::child_tag_mask(tag, which_child, level, max_level,
                                  dimension);
  }
};

template <typename traits> struct HyperOctree<traits>::classify_node {
  const int threshold;
  const int last_level;

  classify_node(int threshold, int last_level)
      : threshold(threshold), last_level(last_level) {}

  template <typename tuple_type>
  inline CUDA_HOST_DEVICE int operator()(const tuple_type &t) const {
    int lower_bound = get<0>(t);
    int upper_bound = get<1>(t);
    const int count = upper_bound - lower_bound;
    if (count == 0) {
      return detail::EMPTY;
    } else if (last_level || count <= threshold) {
      return detail::LEAF;
    } else {
      return detail::NODE;
    }
  }
};

template <typename traits> struct HyperOctree<traits>::write_nodes {
  int num_nodes, num_leaves;

  write_nodes(int num_nodes, int num_leaves)
      : num_nodes(num_nodes), num_leaves(num_leaves) {}

  template <typename tuple_type>
  inline CUDA_HOST_DEVICE int operator()(const tuple_type &t) const {
    int node_type = get<0>(t);
    int node_idx = get<1>(t);
    int leaf_idx = get<2>(t);

    if (node_type == detail::EMPTY) {
      return detail::get_empty_id();
    } else if (node_type == detail::LEAF) {
      return detail::get_leaf_id(num_leaves + leaf_idx);
    } else {
      return num_nodes + nchild * node_idx;
    }
  }
};

template <typename traits> struct HyperOctree<traits>::make_leaf {
  typedef vint2 result_type;
  template <typename tuple_type>
  inline CUDA_HOST_DEVICE result_type operator()(const tuple_type &t) const {
    int x = get<0>(t);
    int y = get<1>(t);

    return result_type(x, y);
  }
};

template <unsigned int D> class octree_child_iterator {
  typedef Vector<double, D> double_d;
  typedef Vector<int, D> int_d;
  typedef Vector<bool, D> bool_d;
  typedef bbox<D> box_type;

  int m_high;
  const int *m_index;
  box_type m_bounds;

public:
  typedef const int *pointer;
  typedef std::forward_iterator_tag iterator_category;
  typedef const int *value_type;
  typedef const int &reference;
  typedef std::ptrdiff_t difference_type;

  CUDA_HOST_DEVICE
  octree_child_iterator() : m_high(1 << D), m_index(nullptr) {}

  CUDA_HOST_DEVICE
  octree_child_iterator(const int *start, const box_type &bounds)
      : m_high(0), m_index(start), m_bounds(bounds) {
    ASSERT_CUDA(start != nullptr);
    // Note: big change, octree iterators now visit empty leaves (for fmm)
    //       be careful of follow-on bugs....
    // if (detail::is_empty(*m_index)) {
    //    increment();
    //}
  }

  CUDA_HOST_DEVICE
  void go_to(const double_d &position) {
    int new_high = 0;
    for (size_t i = 0; i < D; ++i) {
      ASSERT_CUDA(position[i] < m_bounds.bmax[i]);
      ASSERT_CUDA(position[i] >= m_bounds.bmin[i]);
      const double mid = 0.5 * (m_bounds.bmax[i] + m_bounds.bmin[i]);

      // Push the bit into the result as we build it
      new_high |= position[i] < mid ? 0 : 1;
      new_high <<= 1;
    }
    // Unshift the last
    new_high >>= 1;

    m_index = m_index - m_high + new_high;
    m_high = new_high;
  }

  CUDA_HOST_DEVICE
  int get_child_number() const { return m_high; }

  CUDA_HOST_DEVICE
  bool is_high(const size_t i) const { return m_high & (1 << (D - 1 - i)); }

  CUDA_HOST_DEVICE
  box_type get_bounds() const {
    box_type ret = m_bounds;
    for (size_t i = 0; i < D; ++i) {
      if (is_high(i)) {
        ret.bmin[i] = 0.5 * (ret.bmax[i] + ret.bmin[i]);
      } else {
        ret.bmax[i] = 0.5 * (ret.bmax[i] + ret.bmin[i]);
      }
    }
    return ret;
  }

  CUDA_HOST_DEVICE
  reference operator*() const { return dereference(); }

  CUDA_HOST_DEVICE
  reference operator->() const { return dereference(); }

  CUDA_HOST_DEVICE
  octree_child_iterator &operator++() {
    increment();
    return *this;
  }

  CUDA_HOST_DEVICE
  octree_child_iterator operator++(int) {
    octree_child_iterator tmp(*this);
    operator++();
    return tmp;
  }

  CUDA_HOST_DEVICE
  inline bool operator==(const octree_child_iterator &rhs) const {
    return equal(rhs);
  }

  CUDA_HOST_DEVICE
  inline bool operator!=(const octree_child_iterator &rhs) const {
    return !operator==(rhs);
  }

  CUDA_HOST_DEVICE
  inline bool operator==(const bool rhs) const { return equal(rhs); }

  CUDA_HOST_DEVICE
  inline bool operator!=(const bool rhs) const { return !operator==(rhs); }

private:
  CUDA_HOST_DEVICE
  bool equal(octree_child_iterator const &other) const {
    return m_index == other.m_index;
  }

  CUDA_HOST_DEVICE
  bool equal(const bool other) const { return (m_high < (1 << D)) == other; }

  CUDA_HOST_DEVICE
  reference dereference() const { return *m_index; }

  CUDA_HOST_DEVICE
  void increment() {
    ++m_index;
    ++m_high;
    // Note: big change, octree iterators now visit empty leaves (for fmm)
    //       be careful of follow-on bugs....
    // do {
    //    ++m_index;
    //    ++m_high;
    //} while (!equal(false) && detail::is_empty(*m_index));
  }
};

/// @copydetails NeighbourQueryBase
///
/// @brief This is a query object for the HyperOctree spatial data structure
///
template <typename Traits> struct HyperOctreeQuery {
  const static unsigned int dimension = Traits::dimension;
  const static unsigned int m_max_tree_depth = 32 / dimension - 2;

  typedef Traits traits_type;
  typedef typename Traits::raw_pointer raw_pointer;
  typedef typename Traits::double_d double_d;
  typedef typename Traits::bool_d bool_d;
  typedef typename Traits::int_d int_d;
  typedef typename Traits::unsigned_int_d unsigned_int_d;
  template <int LNormNumber, typename Transform = IdentityTransform>
  using query_iterator =
      tree_query_iterator<HyperOctreeQuery, LNormNumber, Transform>;
  typedef depth_first_iterator<HyperOctreeQuery> root_iterator;
  typedef depth_first_iterator<HyperOctreeQuery> all_iterator;
  typedef octree_child_iterator<dimension> child_iterator;
  typedef typename child_iterator::reference reference;
  typedef typename child_iterator::pointer pointer;
  typedef typename child_iterator::value_type value_type;
  typedef bbox<dimension> box_type;

  typedef ranges_iterator<Traits> particle_iterator;

  bool_d m_periodic;
  box_type m_bounds;
  raw_pointer m_particles_begin;
  raw_pointer m_particles_end;
  size_t m_number_of_nodes;
  unsigned m_number_of_levels;

  vint2 *m_leaves_begin;
  int *m_nodes_begin;

  size_t *m_id_map_key;
  size_t *m_id_map_value;

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

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  const box_type &get_bounds() const { return m_bounds; }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  const bool_d &get_periodic() const { return m_periodic; }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  static bool is_leaf_node(reference bucket) {
    return detail::is_leaf(bucket) || detail::is_empty(bucket);
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  static bool is_tree() { return true; }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  child_iterator get_children() const {
    return child_iterator(m_nodes_begin, m_bounds);
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  child_iterator get_children(reference bucket, const box_type &bounds) const {
    CHECK(&bucket == m_nodes_begin, "bucket should be a root bucket");
    return child_iterator(m_nodes_begin, m_bounds);
  }

  ///
  /// @copydoc NeighbourQueryBase::num_children() const
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  static size_t num_children() { return 1 << dimension; }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  static size_t num_children(const child_iterator &ci) {
    if (is_leaf_node(*ci)) {
      return 0;
    } else {
      return (1 << dimension);
    }
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  child_iterator get_children(const child_iterator &ci) const {
    if (!is_leaf_node(*ci)) {
      return child_iterator((m_nodes_begin + *ci), ci.get_bounds());
    } else {
      return child_iterator();
    }
  }

  ///
  /// @copydoc NeighbourQueryBase::get_bounds()
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  static box_type get_bounds(const child_iterator &ci) {
    return ci.get_bounds();
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  particle_iterator CUDA_HOST_DEVICE
  get_bucket_particles(reference bucket) const {
    if (!detail::is_leaf(bucket) || detail::is_empty(bucket)) {
      return particle_iterator();
    }

    const int leaf_idx = detail::get_leaf_offset(bucket);
    const vint2 &particle_idxs = m_leaves_begin[leaf_idx];
#ifndef __CUDA_ARCH__
    LOG(4, "\tget_bucket_particles: looking in bucket with start index = "
               << particle_idxs[0] << " end index = " << particle_idxs[1]);
#endif
    return particle_iterator(m_particles_begin + particle_idxs[0],
                             m_particles_begin + particle_idxs[1]);
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  child_iterator get_bucket(const double_d &position) const {
    child_iterator i = get_children();
    i.go_to(position);

    while (!is_leaf_node(*i)) {
      i = get_children(i);
      i.go_to(position);
    }

    return i;
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  size_t get_bucket_index(reference bucket) const {
    return &bucket - m_nodes_begin;
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  size_t number_of_buckets() const { return m_number_of_nodes; }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  template <int LNormNumber, typename Transform = IdentityTransform>
  CUDA_HOST_DEVICE query_iterator<LNormNumber, Transform>
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

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  iterator_range<root_iterator> get_root_buckets() const {
    return iterator_range<root_iterator>(m_nodes_begin, m_nodes_begin + 1);
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  all_iterator get_subtree(const child_iterator &ci) const {
    return all_iterator(get_children(ci), m_number_of_levels, this);
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  all_iterator get_subtree() const {
    return all_iterator(get_children(), m_number_of_levels, this);
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  size_t number_of_particles() const {
    return m_particles_end - m_particles_begin;
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  unsigned number_of_levels() const { return m_number_of_levels; }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  const raw_pointer &get_particles_begin() const { return m_particles_begin; }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  raw_pointer &get_particles_begin() { return m_particles_begin; }
}; // namespace Aboria

} // namespace Aboria
#endif /* OCTTREE_H_ */
