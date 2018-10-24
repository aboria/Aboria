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

#ifndef NEIGHBOUR_SEARCH_BASE_H_
#define NEIGHBOUR_SEARCH_BASE_H_

#include "CudaInclude.h"
#include "Get.h"
#include "Log.h"
#include "SpatialUtil.h"
#include "StaticVector.h"
#include "Traits.h"
#include "Transform.h"
#include "Vector.h"
#include "detail/Algorithms.h"
#include "detail/Distance.h"
#include <stack>

namespace Aboria {

///
/// @brief lightweight object that holds two iterators to the beginning and end
///        of an STL range
///
/// @tparam IteratorType the type of the iterators to be held
///
template <typename IteratorType> struct iterator_range {
  typedef IteratorType iterator;
  IteratorType m_begin;
  IteratorType m_end;

  /*
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  iterator_range() {}

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  iterator_range(const iterator_range &arg) = default;

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  iterator_range(iterator_range &&arg) = default;

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  iterator_range &operator=(const iterator_range &) = default;
  */

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  iterator_range(const IteratorType &begin, const IteratorType &end)
      : m_begin(begin), m_end(end) {}

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  const IteratorType &begin() const { return m_begin; }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  const IteratorType &end() const { return m_end; }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  IteratorType &begin() { return m_begin; }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  IteratorType &end() { return m_end; }
};

///
/// @brief helper function to make an @ref iterator_range from a begin and end
///        pair of iterators
///
/// @tparam IteratorType the type of the iterators
/// @param begin the iterator pointing to the beginning of the range
/// @param end  the iterator pointing to the end of the range
/// @return CUDA_HOST_DEVICE iterator_range<IteratorType>
///
template <typename IteratorType>
CUDA_HOST_DEVICE iterator_range<IteratorType>
make_iterator_range(IteratorType &&begin, IteratorType &&end) {
  return iterator_range<IteratorType>(begin, end);
}

///
/// @brief A base class for the spatial data structure classes.
///
/// It implements generic functionality such as setting up the spatial domain,
/// logging, and the find-by-id map. In particular see @update_positions
///
/// @todo it uses curiously recurring template pattern (CRTP) to implement
/// compile-time polymorphism for historical reasons, but I don't think this is
/// neccessary anymore as all the member functions are quite slow
///
/// @tparam Derived the derived class type
/// @tparam Traits an instantiation of the @TraitsCommon class
/// @tparam QueryType the query type associated with Derived
///
template <typename Derived, typename Traits, typename QueryType>
class neighbour_search_base {
public:
  typedef QueryType query_type;
  typedef typename Traits::double_d double_d;
  typedef typename Traits::bool_d bool_d;
  typedef typename Traits::iterator iterator;
  typedef typename Traits::vector_unsigned_int vector_unsigned_int;
  typedef typename Traits::vector_int vector_int;
  typedef typename Traits::vector_size_t vector_size_t;
  typedef typename Traits::reference reference;
  typedef typename Traits::raw_reference raw_reference;

  const Derived &cast() const { return static_cast<const Derived &>(*this); }
  Derived &cast() { return static_cast<Derived &>(*this); }

  ///
  /// @brief constructs a new spatial data structure
  ///
  /// The spatial domain is set to 1/3 of the maximum and minimum extents
  /// possible using the `double` type. All periodicity is turned off, and the
  /// number of particle per bucket is set to 10
  ///
  neighbour_search_base() : m_id_map(false) {
    LOG_CUDA(2, "neighbour_search_base: constructor, setting default domain");
    const double min = std::numeric_limits<double>::min();
    const double max = std::numeric_limits<double>::max();
    set_domain(double_d::Constant(min / 3.0), double_d::Constant(max / 3.0),
               bool_d::Constant(false), 10, false);
  };

  ///
  /// @brief Returns true if this spatial data structure relies on the order
  ///        of particles in the @ref Particles container. This is overloaded
  ///        by the Derived class
  ///
  /// @return true
  ///
  static constexpr bool ordered() { return true; }

  ///
  /// @brief A function object used to enforce the domain extents on the set
  ///        of particles
  ///
  /// @tparam D the spatial dimension of the particle set
  /// @tparam Reference a raw reference to a particle in the particle set
  ///
  template <unsigned int D, typename Reference> struct enforce_domain_lambda {
    typedef Vector<double, D> double_d;
    typedef Vector<bool, D> bool_d;
    typedef position_d<D> position;
    const double_d low, high;
    const bool_d periodic;

    enforce_domain_lambda(const double_d &low, const double_d &high,
                          const bool_d &periodic)
        : low(low), high(high), periodic(periodic) {}

    ///
    /// @brief updates the position of @p i (periodic domain) or its alive
    ///        flag (non-periodic domain) based on its position relative
    ///        to the spatial domain
    ///
    /// If dimension $j$ is periodic, and $r_j$ is outside the domain,
    /// then $r_j$ is updated to the correct position within the domain.
    /// If dimensino $j$ is non-periodic, and $r_j$ is outside the domain,
    /// then the `alive` variable for @p i is set to `false`
    ///
    ///
    CUDA_HOST_DEVICE
    void operator()(Reference i) const {
      double_d r = Aboria::get<position>(i);
      for (unsigned int d = 0; d < D; ++d) {
        if (!std::isfinite(r[d])) {
#ifdef __CUDA_ARCH__
          LOG_CUDA(2, "removing particle");
#else
          LOG(2, "removing particle with r = " << r);
#endif
          Aboria::get<alive>(i) = uint8_t(false);
        } else if (periodic[d]) {
          while (r[d] < low[d]) {
            r[d] += (high[d] - low[d]);
          }
          while (r[d] >= high[d]) {
            r[d] -= (high[d] - low[d]);
          }
        } else {
          if ((r[d] < low[d]) || (r[d] >= high[d])) {
#ifdef __CUDA_ARCH__
            LOG_CUDA(2, "removing particle");
#else
            LOG(2, "removing particle with r = " << r);
#endif
            Aboria::get<alive>(i) = uint8_t(false);
          }
        }
      }
      Aboria::get<position>(i) = r;
    }
  };

  ///
  /// @brief resets the domain extents, periodicity and number of particles
  ///        within each bucket
  ///
  /// @param min_in the lower extent of the search domain
  /// @param max_in the upper extent of the search domain
  /// @param periodic_in wether or not each dimension is periodic
  /// @param n_particles_in_leaf indicates the average, or maximum number of
  ///        particles in each bucket
  /// @param not_in_constructor used to determine if this function is called
  ///        within the constructor
  ///
  void set_domain(const double_d &min_in, const double_d &max_in,
                  const bool_d &periodic_in,
                  const double n_particles_in_leaf = 10,
                  const bool not_in_constructor = true) {
    LOG(2, "neighbour_search_base: set_domain:");
    m_domain_has_been_set = not_in_constructor;
    m_bounds.bmin = min_in;
    m_bounds.bmax = max_in;
    m_periodic = periodic_in;
    m_n_particles_in_leaf = n_particles_in_leaf;
    if (not_in_constructor) {
      cast().set_domain_impl();
    }
    LOG(2, "\tbounds = " << m_bounds);
    LOG(2, "\tparticles_in_leaf = " << m_n_particles_in_leaf);
    LOG(2, "\tperiodic = " << m_periodic);
  }

  ///
  /// @brief returns an index into the particle set given a particle id
  ///
  /// @param id the id to search for
  /// @return size_t the index into the particle set
  ///
  size_t find_id_map(const size_t id) const {
    const size_t n = m_particles_end - m_particles_begin;
    return detail::lower_bound(m_id_map_key.begin(), m_id_map_key.begin() + n,
                               id) -
           m_id_map_key.begin();
  }

  ///
  /// @brief This function initialises the find-by-id functionality
  ///
  /// Find-by-id works using a key and value vector pair that act as a map
  /// between ids and particle indicies. This pair is sorted by id for
  /// quick(ish) searching, especially in parallel. It is not as good as
  /// `std::map` on a (single-core) CPU, but can be done on a GPU using
  /// `thrust::vector`
  ///
  /// @see find_id_map
  ///
  void init_id_map() {
    m_id_map = true;
    m_id_map_key.clear();
    m_id_map_value.clear();
  }

  ///
  /// @brief print to stdout the find-by-id id-index map
  ///
  void print_id_map() {
    std::cout << "particle ids:\n";
    for (auto i = m_particles_begin; i != m_particles_end; ++i) {
      std::cout << *get<id>(i) << ',';
    }
    std::cout << std::endl;

    std::cout << "alive indices:\n";
    for (auto i = m_alive_indices.begin(); i != m_alive_indices.end(); ++i) {
      std::cout << *i << ',';
    }
    std::cout << std::endl;

    std::cout << "id map (id,index):\n";
    for (size_t i = 0; i < m_id_map_key.size(); ++i) {
      std::cout << "(" << m_id_map_key[i] << "," << m_id_map_value[i] << ")\n";
    }
    std::cout << std::endl;
  }

  ///
  /// @brief This is the base function for updates to any of the spatial
  ///        data structures. It handles the domain enforcement, the deletion
  ///        or addition of particles, and the find-by-id map. It also calls
  ///        the `update_positions_impl` of the Derived class.
  ///
  /// One of the key responsibilities of this function is to set
  /// #m_alive_indices, which is a vector of indices that are still alive after
  /// taking into account the `alive` flags and the position of the particles
  /// within the domain. The Particles::update_positions() function uses this
  /// vector to delete and reorder the particles in the container
  ///
  /// @param begin The `begin` iterator of the particle set
  /// @param end The `end` iterator of the particle set
  /// @param update_begin The `begin` iterator of the range of particles to
  ///        be updated
  /// @param update_end The `end` iterator of the range of particles to be
  /// updated.
  ///                   Note that if any particles are to be deleted, this must
  ///                   be the same as @p end
  /// @param delete_dead_particles If `true` (the default), the function will
  /// ensure that
  ///        the particles with an `false` `alive` variable are removed from the
  ///        particle set. Note that this function does not actually do this
  /// @return true if the particles need to be reordered/deleted
  /// @return false if the particles don't need to be reordered/deleted
  ///
  bool update_positions(iterator begin, iterator end, iterator update_begin,
                        iterator update_end,
                        const bool delete_dead_particles = true) {

    LOG(2, "neighbour_search_base: update_positions: updating "
               << update_end - update_begin << " points");

    const size_t previous_n = m_particles_end - m_particles_begin;
    m_particles_begin = begin;
    m_particles_end = end;
    const size_t dead_and_alive_n = end - begin;

    CHECK(!cast().ordered() || (update_begin == begin && update_end == end),
          "ordered search data structure can only update the entire particle "
          "set");

    const int new_n = dead_and_alive_n - previous_n;

    // make sure update range covers new particles
    CHECK(
        new_n <= 0 || (update_end == end && update_end - update_begin >= new_n),
        "detected " << new_n
                    << " new particles, which are not covered by update range");

    const size_t update_start_index = update_begin - begin;
    const size_t update_end_index = update_end - begin;
    const size_t update_n = update_end_index - update_start_index;
    if (update_n == 0)
      return false;

    // enforce domain
    if (m_domain_has_been_set) {
      detail::for_each(update_begin, update_end,
                       enforce_domain_lambda<Traits::dimension, raw_reference>(
                           get_min(), get_max(), get_periodic()));
    }

    // m_alive_sum will hold a cummulative sum of the living
    // num_dead holds total number of the dead
    // num_alive_new holds total number of the living that are new particles
    int num_dead = 0;
    m_alive_sum.resize(update_n);
    if (delete_dead_particles) {
      detail::exclusive_scan(get<alive>(update_begin), get<alive>(update_end),
                             m_alive_sum.begin(), 0);
      const int num_alive =
          m_alive_sum.back() + static_cast<int>(*get<alive>(update_end - 1));
      num_dead = update_n - num_alive;
      /*
      if (update_n > new_n) {
          const int num_alive_old = m_alive_sum[update_n-new_n+1];
          num_alive_new = num_alive-num_alive_old;
      } else {
          num_alive_new = num_alive;
      }
      */
    } else {
      detail::sequence(m_alive_sum.begin(), m_alive_sum.end());
    }

    CHECK(update_end == end || num_dead == 0,
          "cannot delete dead points if not updating the end of the vector");

    LOG(2, "neighbour_search_base: update_positions: found "
               << num_dead << " dead points, and " << new_n << " new points");

    // m_alive_indices holds particle set indicies that are alive
    m_alive_indices.resize(update_n - num_dead);

#if defined(__CUDACC__)
    typedef typename thrust::detail::iterator_category_to_system<
        typename vector_int::iterator::iterator_category>::type system;
    thrust::counting_iterator<unsigned int, system> count_start(
        update_start_index);
    thrust::counting_iterator<unsigned int, system> count_end(update_end_index);
#else
    auto count_start = Traits::make_counting_iterator(update_start_index);
    auto count_end = Traits::make_counting_iterator(update_end_index);
#endif

    // scatter alive indicies to m_alive_indicies
    detail::scatter_if(count_start, count_end,   // items to scatter
                       m_alive_sum.begin(),      // map
                       get<alive>(update_begin), // scattered if true
                       m_alive_indices.begin());

    if (m_domain_has_been_set) {
      LOG(2, "neighbour_search_base: update_positions_impl:");
      cast().update_positions_impl(update_begin, update_end, new_n);
    }
    if (m_id_map) {
      LOG(2, "neighbour_search_base: update_id_map:");
      // if no new particles, no dead, no reorder, or no init than can assume
      // that previous id map is correct
      if (cast().ordered() || new_n > 0 || num_dead > 0 ||
          m_id_map_key.size() == 0) {
        m_id_map_key.resize(dead_and_alive_n - num_dead);
        m_id_map_value.resize(dead_and_alive_n - num_dead);

        // before update range
        if (update_start_index > 0) {
          detail::sequence(m_id_map_value.begin(),
                           m_id_map_value.begin() + update_start_index);
          detail::copy(get<id>(begin), get<id>(begin) + update_start_index,
                       m_id_map_key.begin());
        }

        // after update range
        ASSERT(update_end_index == dead_and_alive_n,
               "if not updateing last particle then should not get here");

        // update range
        /*
        detail::transform(m_alive_indices.begin(),m_alive_indices.end(),
                     m_id_map_value.begin()+update_start_index,
                     [&](const int index) {
                        const int index_into_update = index -
        update_start_index; const int num_dead_before_index = index_into_update
        - m_alive_sum[index_into_update]; return index - num_dead_before_index;
                     });
                     */
        detail::sequence(m_id_map_value.begin() + update_start_index,
                         m_id_map_value.end(), update_start_index);
        auto raw_id = iterator_to_raw_pointer(get<id>(begin));
        detail::transform(
            m_alive_indices.begin(), m_alive_indices.end(),
            m_id_map_key.begin() + update_start_index,
            [=] CUDA_HOST_DEVICE(const int index) { return raw_id[index]; });
        detail::sort_by_key(m_id_map_key.begin(), m_id_map_key.end(),
                            m_id_map_value.begin());
#ifndef __CUDA_ARCH__
        if (4 <= ABORIA_LOG_LEVEL) {
          print_id_map();
        }
#endif
      }
    }

    query_type &query = cast().get_query_impl();
    query.m_id_map_key = iterator_to_raw_pointer(m_id_map_key.begin());
    query.m_id_map_value = iterator_to_raw_pointer(m_id_map_value.begin());
    query.m_particles_begin = iterator_to_raw_pointer(m_particles_begin);
    query.m_particles_end = iterator_to_raw_pointer(m_particles_end);

    return cast().ordered() || num_dead > 0;
  }

  ///
  /// @brief updates the internal copies held of the `begin` and `end` iterators
  ///         of the particle set
  ///
  /// @param begin the `begin` iterator of the particle set
  /// @param end the `end` iterator of the particle set
  ///
  void update_iterators(iterator begin, iterator end) {
    m_particles_begin = begin;
    m_particles_end = end;
    query_type &query = cast().get_query_impl();
    query.m_particles_begin = iterator_to_raw_pointer(m_particles_begin);
    query.m_particles_end = iterator_to_raw_pointer(m_particles_end);
    LOG(2, "neighbour_search_base: update iterators");
    cast().update_iterator_impl();
  }

  ///
  /// @return the query object
  ///
  const query_type &get_query() const { return cast().get_query_impl(); }

  ///
  /// @return a vector of ints that shows the new order of the particles
  ///        within the particle container
  /// @see update_positions(), #m_alive_indices
  ///
  const vector_int &get_alive_indicies() const { return m_alive_indices; }

  ///
  /// @return true if find-by-id functionality switched on
  ///
  bool get_id_map() const { return m_id_map; }

  ///
  /// @return the minimum extent of the cuboid domain
  ///
  const double_d &get_min() const { return m_bounds.bmin; }

  ///
  /// @return the maximum extent of the cuboid domain
  ///
  const double_d &get_max() const { return m_bounds.bmax; }

  ///
  /// @return the periodicity of the domain
  ///
  const bool_d &get_periodic() const { return m_periodic; }

  ///
  /// @return true if domain has been initialised
  /// @see set_domain()
  ///
  bool domain_has_been_set() const { return m_domain_has_been_set; }

  ///
  /// @return the number of particles in each bucket (average or maximum)
  ///
  double get_max_bucket_size() const { return m_n_particles_in_leaf; }

protected:
  ///
  /// @brief a copy of the `begin` iterator for the particle set
  ///
  iterator m_particles_begin;

  ///
  /// @brief a copy of the `end` iterator for the particle set
  ///
  iterator m_particles_end;

  ///
  /// @brief a vector of ints holding a cummulative sum of the living particles
  /// @see update_positions()
  ///
  vector_int m_alive_sum;

  ///
  /// @brief A vector of particle set indices that are used by
  /// Particles::update_positions() to reorder and delete the particles
  ///
  vector_int m_alive_indices;

  ///
  /// @brief The `key` vector of the id->index map for the find-by-id
  /// functionality
  ///
  vector_size_t m_id_map_key;

  ///
  /// @brief The `value` vector of the id->index map for the find-by-id
  /// functionality
  ///
  vector_size_t m_id_map_value;

  ///
  /// @brief flag set to `true` if find-by-id functionality is turned on
  ///
  bool m_id_map;

  ///
  /// @brief @Vector of bools indicating the periodicity of the domain
  ///
  bool_d m_periodic;

  ///
  /// @brief true if the domain has been set
  /// @see set_domain()
  ///
  bool m_domain_has_been_set;

  ///
  /// @brief the bounds (min->max) of the domain
  ///
  bbox<Traits::dimension> m_bounds;

  ///
  /// @brief holds the set number of particles in each bucket
  ///
  ///
  double m_n_particles_in_leaf;
};

///
/// @brief a lightweight query object that should be used for read-only access
/// to the spatial data structures
///
/// This object is designed to be used within algorithms or gpu kernels, and can
/// be obtained from a Particles type using the Particles::get_query() function
///
/// @code
/// auto query = particles.get_query();
/// std::for_each(particles.begin(), particles.end(), [&] (auto i) {
///   int count = 0;
///   for (auto j = euclidean_search(query, get<position>(i), diameter);
///                   j != false; ++j) {
///     ++count;
///   }
/// });
/// @endcode
///
/// Note that if you wish to run this on a gpu using Trust, you need to copy the
/// query object to the kernel, and mark the lambda function as a device
/// function
///
/// @code
/// auto query = particles.get_query();
/// thrust::for_each(particles.begin(), particles.end(),
///   [=] __host__ __device__ (auto i) {
///     int count = 0;
///     for (auto j = euclidean_search(query, get<position>(i), diameter);
///                   j != false; ++j) {
///       ++count;
///     }
/// });
/// @endcode
///
/// @tparam Traits the @ref TraitsCommon type
template <typename Traits> struct NeighbourQueryBase {

  typedef typename Traits::raw_pointer raw_pointer;
  typedef typename Traits::double_d double_d;
  typedef typename Traits::bool_d bool_d;
  typedef typename Traits::int_d int_d;
  typedef typename Traits::unsigned_int_d unsigned_int_d;
  typedef typename Traits::reference particle_reference;
  typedef typename Traits::const_reference particle_const_reference;
  const static unsigned int dimension = Traits::dimension;
  typedef bbox<dimension> box_type;

  template <int LNormNumber> struct query_iterator {
    ///
    /// A reference type to a bucket
    struct reference;
    ///
    /// A pointer type to a bucket
    struct pointer;
    ///
    /// The value_type of a bucket
    struct value_type;
  };

  ///
  /// @brief An iterator that steps through all the buckets in the tree
  /// (depth-first)
  ///
  struct all_iterator {};
  ///
  /// @brief An iterator that steps through the children of a single bucket in
  /// the tree
  ///
  struct child_iterator {};

  ///
  /// @brief An iterator that steps through the particles within a given
  /// bucket
  ///
  struct particle_iterator {};

  ///
  /// @brief A reference to a bucket in the tree
  ///
  typedef typename query_iterator<2>::reference reference;

  ///
  /// @brief A pointer to a bucket in the tree
  ///
  typedef typename query_iterator<2>::pointer pointer;

  ///
  /// @brief A value_type for a bucket in the tree
  ///
  typedef typename query_iterator<2>::value_type value_type;

  ///
  /// @brief implement find-by-id
  ///
  /// performs a binary search for the id in the map
  ///
  /// @param id the id of the particle to find
  /// @return pointer to the particle found, or a pointer to the end of the
  /// particle set if not found
  ///
  raw_pointer find(const size_t id) const;

  ///
  /// @brief given a reference to a bucket, checks if that bucket is a leaf
  /// (i.e. does not have any children)
  ///
  /// @param bucket the bucket to check
  /// @return true if @p is a bucket
  /// @return false if @p is not a bucket
  ///
  bool is_leaf_node(const value_type &bucket);

  ///
  /// @return true if this data structure is a tree (e.g. kdtree, HyperOctree)
  /// @return false if this data structure is not a tree (e.g. cell list)
  ///
  bool is_tree();

  ///
  /// @brief gets all the children of the root node.
  ///
  /// @return a @ref child_iterator that iterates through all the buckets in the
  /// data structure
  ///
  child_iterator get_children() const;

  ///
  /// @brief returns all the children of the given child_iterator
  ///
  /// @param ci this fuction returns all the children of @p ci
  /// @return same as get_children()
  ///
  child_iterator get_children(const child_iterator &ci) const;

  ///
  /// @brief returns the number of children of the root node
  ///
  size_t num_children() const;

  ///
  /// @brief returns the number of children of a given child iterator @p ci
  ///
  size_t num_children(const child_iterator &ci) const;

  ///
  /// @brief returns the min/max bounds of the given child_iterator @p ci
  ///
  /// @return a @bbox containing the bounds
  ///
  const box_type get_bounds(const child_iterator &ci) const;

  ///
  /// @brief returns entire domain min/max bounds
  ///
  const box_type &get_bounds() const;

  ///
  /// @brief returns the periodicity of the domain
  ///
  const bool_d &get_periodic() const;

  ///
  /// @brief returns a particle_iterator range to all the particles within the
  /// given bucket
  ///
  /// @param bucket a reference to the bucket in question
  ///
  particle_iterator get_bucket_particles(const reference bucket);

  ///
  /// @brief get min/max bounds of given bucket @p bucket
  ///
  bbox<dimension> get_bucket_bbox(const reference bucket) const;

  ///
  /// @brief given a @p position, returns a @p bucket and a min/max @p bounds
  ///
  /// @param position (input)
  /// @param bucket (output)
  /// @param bounds (output)
  ///
  void get_bucket(const double_d &position, pointer &bucket,
                  box_type &bounds) const;

  ///
  /// @brief returns a bucket index/id given a @p bucket reference
  ///
  size_t get_bucket_index(const reference bucket) const;

  ///
  /// @brief returns all the bucket within a distance of a point
  ///
  /// This function can use multiple p-norm distance types by setting @p
  /// LNormNumber, and uses a isotropic distance value given by @p max_distance.
  /// Note that only buckets within the domain are returned, and periodicity is
  /// not taken into account
  ///
  /// @param position the position to search around
  /// @param max_distance the maximum distance of the buckets to be returned
  /// @param transform a transform function object allowing the user to search
  /// in a transformed space (e.g. shew coordinate systems). @see
  /// IdentityTransform for an example
  /// @return an @ref iterator_range containing all the buckets found
  ///
  template <int LNormNumber, typename Transform>
  query_iterator<LNormNumber> get_buckets_near_point(
      const double_d &position, const double max_distance,
      const Transform &transform = IdentityTransform()) const;

  ///
  /// @brief returns all the bucket within an anisotropic distance of a point
  ///
  /// This function can use multiple p-norm distance types by setting @p
  /// LNormNumber, and uses a anisotropic distance value given by @p
  /// max_distance. Note that only buckets within the domain are returned, and
  /// periodicity is not taken into account
  ///
  /// @param position the position to search around
  /// @param max_distance the maximum distance of the buckets to be returned
  /// @return an @ref iterator_range containing all the buckets found
  ///
  template <int LNormNumber = -1>
  query_iterator<LNormNumber>
  get_buckets_near_point(const double_d &position,
                         const double_d &max_distance) const;

  ///
  /// @brief get index of the last bucket in the cell list
  ///
  const int_d &get_end_bucket() const;

  ///
  /// @brief return an @ref all_iterator to the entire tree under
  /// the given child_iterator @p ci
  ///
  /// @param ci the child_iterator to search under
  ///
  all_iterator get_subtree(const child_iterator &ci) const;

  ///
  /// @brief return an @ref all_iterator to the entire tree data structure
  ///
  /// Can use this range to loop through all the buckets in the data structure,
  /// wether it is a tree or not
  ///
  all_iterator get_subtree() const;

  ///
  /// @brief return the total number of buckets in the data structure
  ///
  size_t number_of_buckets() const;

  ///
  /// @brief return the total number of particles in the data structure
  ///
  size_t number_of_particles() const;

  ///
  /// @brief get a pointer to the beginning of the particle container
  ///
  /// can use this and @ref number_of_particles() to loop through all the
  /// particles
  ///
  const raw_pointer &get_particles_begin() const;

  ///
  /// @brief get a pointer to the beginning of the particle container
  ///
  /// can use this and @ref number_of_particles() to loop through all the
  /// particles
  ///
  raw_pointer &get_particles_begin();

  ///
  /// @brief returns the number of levels in the tree
  ///
  unsigned number_of_levels() const;
};

// assume that these iterators, and query functions, are only called from device
// code
template <typename Traits> class ranges_iterator {
  typedef typename Traits::double_d double_d;
  typedef typename Traits::bool_d bool_d;
  typedef typename Traits::value_type p_value_type;
  typedef typename Traits::raw_reference p_reference;
  typedef typename Traits::raw_pointer p_pointer;

public:
  typedef Traits traits_type;
  typedef const p_pointer pointer;
  typedef std::random_access_iterator_tag iterator_category;
  typedef const p_reference reference;
  typedef const p_reference value_type;
  typedef std::ptrdiff_t difference_type;

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  CUDA_HOST_DEVICE
  ranges_iterator() {}

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  CUDA_HOST_DEVICE
  ranges_iterator(const p_pointer &begin, const p_pointer &end)
      : m_current_p(begin), m_end_p(end) {}

  size_t distance_to_end() const { return m_end_p - m_current_p; }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  CUDA_HOST_DEVICE
  reference operator*() const { return dereference(); }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  CUDA_HOST_DEVICE
  reference operator->() const { return dereference(); }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  CUDA_HOST_DEVICE
  ranges_iterator &operator++() {
    increment();
    return *this;
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  CUDA_HOST_DEVICE
  ranges_iterator operator++(int) {
    ranges_iterator tmp(*this);
    operator++();
    return tmp;
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  ranges_iterator operator+(int n) {
    ranges_iterator tmp(*this);
    tmp.increment(n);
    return tmp;
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  CUDA_HOST_DEVICE
  size_t operator-(ranges_iterator start) const {
    return get_by_index<0>(m_current_p) - get_by_index<0>(start.m_current_p);
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  CUDA_HOST_DEVICE
  inline bool operator==(const ranges_iterator &rhs) const {
    return equal(rhs);
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  CUDA_HOST_DEVICE
  inline bool operator!=(const ranges_iterator &rhs) const {
    return !operator==(rhs);
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  CUDA_HOST_DEVICE
  inline bool operator==(const bool rhs) const { return equal(rhs); }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  CUDA_HOST_DEVICE
  inline bool operator!=(const bool rhs) const { return !operator==(rhs); }

private:
  friend class boost::iterator_core_access;

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  CUDA_HOST_DEVICE
  bool equal(ranges_iterator const &other) const {
    return m_current_p == other.m_current_p;
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  CUDA_HOST_DEVICE
  bool equal(const bool other) const {
    return (m_current_p < m_end_p) == other;
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  CUDA_HOST_DEVICE
  reference dereference() const { return *m_current_p; }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  CUDA_HOST_DEVICE
  void increment() { ++m_current_p; }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  CUDA_HOST_DEVICE
  void increment(const int n) { m_current_p += n; }

  p_pointer m_current_p;
  p_pointer m_end_p;
};

/// A const iterator to a set of neighbouring points. This iterator implements
/// a STL forward iterator type
// assume that these iterators, and query functions, are only called from device
// code
template <typename Traits> class linked_list_iterator {
  typedef typename Traits::double_d double_d;
  typedef typename Traits::bool_d bool_d;
  typedef typename Traits::value_type p_value_type;
  typedef typename Traits::raw_reference p_reference;
  typedef typename Traits::raw_pointer p_pointer;

public:
  typedef Traits traits_type;
  typedef const p_pointer pointer;
  typedef std::forward_iterator_tag iterator_category;
  typedef const p_reference reference;
  typedef const p_reference value_type;
  typedef std::ptrdiff_t difference_type;

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  linked_list_iterator() : m_current_index(detail::get_empty_id()) {
#if defined(__CUDA_ARCH__)
    CHECK_CUDA((!std::is_same<typename Traits::template vector<double>,
                              std::vector<double>>::value),
               "Cannot use std::vector in device code");
#endif
  }

  /// this constructor is used to start the iterator at the head of a bucket
  /// list
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  linked_list_iterator(const int index, const p_pointer &particles_begin,
                       int *const linked_list_begin)
      : m_current_index(index), m_particles_begin(particles_begin),
        m_linked_list_begin(linked_list_begin) {}

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  linked_list_iterator(const linked_list_iterator &other)
      : m_current_index(other.m_current_index),
        m_particles_begin(other.m_particles_begin),
        m_linked_list_begin(other.m_linked_list_begin) {}

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  void operator=(const linked_list_iterator &other) {
    m_current_index = other.m_current_index;
    if (get_by_index<0>(m_particles_begin) !=
        get_by_index<0>(other.m_particles_begin)) {
      m_particles_begin = other.m_particles_begin;
    }
    m_linked_list_begin = other.m_linked_list_begin;
  }

  size_t distance_to_end() const {
    size_t count = 0;
    for (auto i = *this; i != false; ++i) {
      ++count;
    }
    return count;
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  reference operator*() const { return dereference(); }
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  reference operator->() { return dereference(); }
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  linked_list_iterator &operator++() {
    increment();
    return *this;
  }
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  linked_list_iterator operator++(int) {
    linked_list_iterator tmp(*this);
    operator++();
    return tmp;
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  linked_list_iterator operator+(int n) {
    linked_list_iterator tmp(*this);
    for (int i = 0; i < n; ++i) {
      tmp.increment();
    }
    return tmp;
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  size_t operator-(linked_list_iterator start) const {
    size_t count = 0;
    while (start != *this) {
      start++;
      count++;
    }
    return count;
  }
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  inline bool operator==(const linked_list_iterator &rhs) const {
    return equal(rhs);
  }
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  inline bool operator!=(const linked_list_iterator &rhs) const {
    return !operator==(rhs);
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  inline bool operator==(const bool rhs) const { return equal(rhs); }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  inline bool operator!=(const bool rhs) const { return !operator==(rhs); }

private:
  friend class boost::iterator_core_access;

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  bool increment() {
#ifndef __CUDA_ARCH__
    LOG(4, "\tincrement (linked_list_iterator):");
#endif
    if (m_current_index != detail::get_empty_id()) {
      m_current_index = m_linked_list_begin[m_current_index];
#ifndef __CUDA_ARCH__
      LOG(4, "\tgoing to new particle m_current_index = " << m_current_index);
#endif
    }

#ifndef __CUDA_ARCH__
    LOG(4, "\tend increment (linked_list_iterator): m_current_index = "
               << m_current_index);
#endif
    if (m_current_index == detail::get_empty_id()) {
      return false;
    } else {
      return true;
    }
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  bool equal(linked_list_iterator const &other) const {
    return m_current_index == other.m_current_index;
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  bool equal(const bool other) const {
    return (m_current_index != detail::get_empty_id()) == other;
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  reference dereference() const {
    return *(m_particles_begin + m_current_index);
  }

  int m_current_index;
  p_pointer m_particles_begin;
  int *m_linked_list_begin;
};

template <typename Traits, typename Iterator> class index_vector_iterator {
  typedef index_vector_iterator<Traits, Iterator> iterator;
  typedef typename Traits::double_d double_d;
  typedef typename Traits::bool_d bool_d;
  typedef typename Traits::value_type p_value_type;
  typedef typename Traits::raw_reference p_reference;
  typedef typename Traits::raw_pointer p_pointer;

public:
  typedef Traits traits_type;
  typedef const p_pointer pointer;
  typedef std::forward_iterator_tag iterator_category;
  typedef const p_reference reference;
  typedef const p_reference value_type;
  typedef std::ptrdiff_t difference_type;

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  index_vector_iterator() {}

  /// this constructor is used to start the iterator at the head of a bucket
  /// list
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  index_vector_iterator(Iterator begin, const p_pointer &particles_begin)
      : m_current_index(begin), m_particles_begin(particles_begin) {}

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  reference operator*() const { return dereference(); }
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  reference operator->() { return dereference(); }
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  iterator &operator++() {
    increment();
    return *this;
  }
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  iterator operator++(int) {
    iterator tmp(*this);
    operator++();
    return tmp;
  }
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  size_t operator-(iterator start) const {
    size_t count = 0;
    while (start != *this) {
      start++;
      count++;
    }
    return count;
  }
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  inline bool operator==(const iterator &rhs) const { return equal(rhs); }
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  inline bool operator!=(const iterator &rhs) const { return !operator==(rhs); }

private:
  friend class boost::iterator_core_access;

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  void increment() {
#ifndef __CUDA_ARCH__
    LOG(4, "\tincrement (index_vector_iterator):");
#endif
    ++m_current_index;
#ifndef __CUDA_ARCH__
    LOG(4, "\tend increment (index_vector_iterator): m_current_index = "
               << m_current_index);
#endif
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  bool equal(iterator const &other) const {
    return m_current_index == other.m_current_index;
  }

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  reference dereference() const {
    return *(m_particles_begin + *m_current_index);
  }

  Iterator m_current_index;
  p_pointer m_particles_begin;
};

template <typename Query> class depth_first_iterator {
  typedef depth_first_iterator<Query> iterator;
  typedef typename Query::child_iterator child_iterator;
  static const unsigned int dimension = Query::dimension;
  typedef Vector<double, dimension> double_d;
  typedef Vector<int, dimension> int_d;
  typedef bbox<dimension> box_type;

public:
  typedef typename child_iterator::value_type value_type;
  typedef child_iterator pointer;
  typedef std::forward_iterator_tag iterator_category;
  typedef typename child_iterator::reference reference;
  typedef std::ptrdiff_t difference_type;

  CUDA_HOST_DEVICE
  depth_first_iterator() {}

  /// this constructor is used to start the iterator at the head of a bucket
  /// list
  CUDA_HOST_DEVICE
  depth_first_iterator(const child_iterator &start_node,
                       const unsigned tree_depth, const Query *query)
      : m_query(query) {
    ASSERT_CUDA(tree_depth <= m_stack_max_size);
    if (start_node != false) {
      m_stack.push_back(start_node);
    } else {
#ifndef __CUDA_ARCH__
      LOG(3, "\tdepth_first_iterator (constructor): start is false, no "
             "children to search.");
#endif
    }
  }

  CUDA_HOST_DEVICE
  depth_first_iterator(const iterator &copy) : m_query(copy.m_query) {
    for (size_t i = 0; i < copy.m_stack.size(); ++i) {
      m_stack.push_back(copy.m_stack[i]);
    }
  }

  ~depth_first_iterator() {}

  CUDA_HOST_DEVICE
  iterator &operator=(const iterator &copy) {
    m_stack.resize(copy.m_stack.size());
    for (size_t i = 0; i < copy.m_stack.size(); ++i) {
      m_stack[i] = copy.m_stack[i];
    }

    m_query = copy.m_query;
    return *this;
  }

  const child_iterator &get_child_iterator() const { return m_stack.back(); }

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
  void increment() {
#ifndef __CUDA_ARCH__
    LOG(4, "\tincrement (depth_first_iterator): depth = "
               << m_stack.size() << " top child number "
               << m_stack.back().get_child_number());
#endif
    if (m_query->is_leaf_node(*m_stack.back())) {
      ++m_stack.back();
      while (!m_stack.empty() && m_stack.back() == false) {
        LOG(4, "\tpop stack (depth_first_iterator):");
        m_stack.pop_back();
      }
    } else {
      LOG(4, "\tpush stack (depth_first_iterator):");
      m_stack.push_back(m_query->get_children(m_stack.back()++));
    }

#ifndef __CUDA_ARCH__
    LOG(4, "\tend increment (depth_first_iterator):");
#endif
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

  const static unsigned m_stack_max_size = Query::m_max_tree_depth;
  static_vector<child_iterator, m_stack_max_size> m_stack;
  const Query *m_query;
};

template <typename Query, int LNormNumber, typename Transform>
class tree_query_iterator {
  typedef tree_query_iterator<Query, LNormNumber, Transform> iterator;
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
  tree_query_iterator() {}

  /// this constructor is used to start the iterator at the head of a bucket
  /// list
  CUDA_HOST_DEVICE
  tree_query_iterator(const child_iterator &start, const double_d &query_point,
                      const double &max_distance, const unsigned tree_depth,
                      const Query *query, const Transform &transform)
      : m_query_point(query_point),
        m_max_distance2(
            detail::distance_helper<LNormNumber>::get_value_to_accumulate(
                max_distance)),
        m_query(query), m_transform(transform) {
    ASSERT_CUDA(tree_depth <= m_stack_max_size);
    if (start != false) {
      m_stack.push_back(start);
      go_to_next_leaf();
    } else {
#ifndef __CUDA_ARCH__
      LOG(3, "\ttree_query_iterator (constructor) with query pt = "
                 << m_query_point
                 << "): start is false, no children to search.");
#endif
    }

    if (m_stack.empty()) {
#ifndef __CUDA_ARCH__
      LOG(3,
          "\ttree_query_iterator (constructor) with query pt = "
              << m_query_point
              << "): search region outside domain or no children to search.");
#endif
    } else {
#ifndef __CUDA_ARCH__
      LOG(3, "\ttree_query_iterator (constructor) with query pt = "
                 << m_query_point
                 << "):  found bbox = " << m_query->get_bounds(m_stack.back()));
#endif
    }
  }

  CUDA_HOST_DEVICE
  tree_query_iterator(const iterator &copy)
      : m_query_point(copy.m_query_point),
        m_max_distance2(copy.m_max_distance2), m_query(copy.m_query),
        m_transform(copy.m_transform)

  {
    for (size_t i = 0; i < copy.m_stack.size(); ++i) {
      m_stack.push_back(copy.m_stack[i]);
    }
  }

  CUDA_HOST_DEVICE
  ~tree_query_iterator() {}

  CUDA_HOST_DEVICE
  child_iterator get_child_iterator() const { return m_stack.back(); }

  CUDA_HOST_DEVICE
  iterator &operator=(const iterator &copy) {
    m_query_point = copy.m_query_point;
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
  friend class boost::iterator_core_access;

  CUDA_HOST_DEVICE
  double get_dist_to_bucket(const box_type &bucket) const {
    const double_d half_bucket_side_length = 0.5 * m_transform(bucket);
    double_d dx =
        m_transform(0.5 * (bucket.bmin + bucket.bmax) - m_query_point);
    for (size_t i = 0; i < dimension; ++i) {
      dx[i] = std::max(std::abs(dx[i]) - half_bucket_side_length[i], 0.0);
    }
    return detail::distance_helper<LNormNumber>::norm2(dx);
  }

  /*
  CUDA_HOST_DEVICE
  double get_dist_to_bucket(const box_type &bucket, std::false_type) const {
    double_d dx = 0.5 * (bucket.bmin + bucket.bmax) - m_query_point;
    for (int i = 0; i < dimension; ++i) {
      dx[i] = std::copysign(
          std::max(std::abs(dx[i]) - 0.5 * (bucket.bmax[i] - bucket.bmin[i]),
                   0.0),
          dx[i]);
    }
    // TODO: this 0.9 is a bit of a fudge factor to make sure we get all
    // relevent buckets, needs improvement
    return 0.9 * detail::distance_helper<LNormNumber>::norm2(m_transform(dx));
  }
  */

  CUDA_HOST_DEVICE
  bool child_is_within_query(const child_iterator &node) {
    const box_type &bounds = m_query->get_bounds(node);
    const double accum = get_dist_to_bucket(bounds);
    // std::cout <<"accum = "<<accum<< std::endl;
    return (accum < m_max_distance2);
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
    LOG(4, "\tincrement (octtree_query_iterator):");
#endif
    increment_stack();
    go_to_next_leaf();

    if (m_stack.empty()) {
#ifndef __CUDA_ARCH__
      LOG(3, "\tend increment (octree_query_iterator): no more nodes");
#endif
    } else {
#ifndef __CUDA_ARCH__
      LOG(3, "\tend increment (octree_query_iterator): looking in bbox "
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
  double m_max_distance2;
  const Query *m_query;
  Transform m_transform;
};

template <typename Query, int LNormNumber, typename Transform>
class lattice_iterator_within_distance {
  typedef lattice_iterator_within_distance<Query, LNormNumber, Transform>
      iterator;
  static const unsigned int dimension = Query::dimension;
  typedef Vector<double, dimension> double_d;
  typedef Vector<int, dimension> int_d;

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

  double_d m_query_point;
  double_d m_half_bucket_length;
  double m_max_distance2;
  int m_quadrant;
  const Query *m_query;
  bool m_valid;
  int_d m_min;
  proxy_int_d m_index;
  Transform m_transform;

public:
  typedef proxy_int_d pointer;
  typedef std::random_access_iterator_tag iterator_category;
  typedef const proxy_int_d &reference;
  typedef proxy_int_d value_type;
  typedef std::ptrdiff_t difference_type;

  CUDA_HOST_DEVICE
  lattice_iterator_within_distance() : m_valid(false) {}

  CUDA_HOST_DEVICE
  lattice_iterator_within_distance(const double_d &query_point,
                                   const double max_distance,
                                   const Query *query,
                                   const Transform &transform)
      : m_query_point(query_point),
        m_max_distance2(
            detail::distance_helper<LNormNumber>::get_value_to_accumulate(
                max_distance)),
        m_quadrant(0), m_query(query), m_valid(true), m_transform(transform) {
    if (outside_domain(query_point)) {
      m_valid = false;
    } else {
      if (std::is_same<Transform, IdentityTransform>::value) {
        m_half_bucket_length =
            0.5 * m_query->m_point_to_bucket_index.m_bucket_side_length;
      } else {
        const auto &bucket_side_length =
            m_query->m_point_to_bucket_index.m_bucket_side_length;
        bbox<dimension> bucket_bounds(-0.5 * bucket_side_length,
                                      0.5 * bucket_side_length);
        m_half_bucket_length = 0.5 * m_transform(bucket_bounds);
      }
      reset_min_and_index();
    }
  }

  CUDA_HOST_DEVICE
  explicit operator size_t() const {
    return m_query->m_point_to_bucket_index.collapse_index_vector(m_index);
  }

  lattice_iterator<dimension> get_child_iterator() const {
    lattice_iterator<dimension> ret = m_query->get_subtree();
    ret = m_index;
    return ret;
  }

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
  size_t operator-(const iterator &start) const {
    int distance = 0;
    iterator tmp = start;
    while (tmp != *this) {
      ++distance;
      ++tmp;
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
  CUDA_HOST_DEVICE
  bool equal(iterator const &other) const {
    if (!other.m_valid)
      return !m_valid;
    if (!m_valid)
      return !other.m_valid;
    for (size_t i = 0; i < dimension; ++i) {
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
  bool ith_quadrant_bit(const int i) const {
    return (1 == ((m_quadrant >> i) & 1));
  }

  double get_min_distance_to_bucket(const int_d &bucket) {
    double_d dx = m_transform(
        m_query->m_point_to_bucket_index.find_bucket_centre(bucket) -
        m_query_point);
    for (size_t i = 0; i < dimension; ++i) {
      dx[i] = std::max(std::abs(dx[i]) - m_half_bucket_length[i], 0.0);
    }
    return detail::distance_helper<LNormNumber>::norm2(dx);
  }

  CUDA_HOST_DEVICE
  void reset_min_and_index() {
    bool no_buckets = true;

    LOG_CUDA(3, "lattice_iterator_within_distance: reset_min_and_index:begin");
    while (m_valid && no_buckets) {
      for (size_t i = 0; i < dimension; ++i) {
        m_min[i] = m_query->m_point_to_bucket_index.get_min_index_by_quadrant(
            m_query_point[i], i, ith_quadrant_bit(i));
      }

      const double accum = get_min_distance_to_bucket(m_min);
      // std::cout <<"accum = "<<accum<< std::endl;

      no_buckets = accum > m_max_distance2;

      // std::cout <<" m_min = "<<m_min<<" m_quadrant = "<<m_quadrant <<
      // std::endl;

      // if good, check that this quadrant is within domain
      if (!no_buckets) {
        for (size_t i = 0; i < dimension; i++) {
          if (ith_quadrant_bit(i)) {
            if (m_min[i] < 0) {
              m_min[i] = 0;
            } else if (m_min[i] > m_query->m_end_bucket[i]) {
              no_buckets = true;
              m_min[i] = m_query->m_end_bucket[i];
            }
          } else {
            if (m_min[i] < 0) {
              no_buckets = true;
              m_min[i] = 0;
            } else if (m_min[i] > m_query->m_end_bucket[i]) {
              m_min[i] = m_query->m_end_bucket[i];
            }
          }
        }
      }
      if (no_buckets) {
        // if no buckets, move onto next quadrant
        ++m_quadrant;
        if (m_quadrant >= (1 << dimension)) {
          m_valid = false;
        }
      } else {
        // we got a non empty quadrent, lets go!
        m_index = m_min;
      }
    }
    // std::cout <<"m_valid = "<<m_valid<<" m_min = "<<m_min<< "m_index =
    // "<<m_index<<" m_quadrant = "<<m_quadrant << std::endl;
    LOG_CUDA(3, "lattice_iterator_within_distance: reset_min_and_index:end");
  }

  CUDA_HOST_DEVICE
  bool outside_domain(const double_d &position) {
    const auto &bounds = m_query->get_bounds();
    double_d dx = m_transform(0.5 * (bounds.bmin + bounds.bmax) - position);
    const double_d half_domain_side_length = 0.5 * m_transform(bounds);
    for (size_t i = 0; i < dimension; ++i) {
      dx[i] = std::max(std::abs(dx[i]) - half_domain_side_length[i], 0.0);
    }
    return detail::distance_helper<LNormNumber>::norm2(dx) > m_max_distance2;
  }

  CUDA_HOST_DEVICE void increment() {
    LOG_CUDA(3, "lattice_iterator_within_distance: increment :begin");
    for (int i = dimension - 1; i >= 0; --i) {

      // increment or decrement index depending on the current
      // quadrant
      bool potential_bucket = true;
      if (ith_quadrant_bit(i)) {
        ++m_index[i];
        potential_bucket = m_index[i] <= m_query->m_end_bucket[i];
      } else {
        --m_index[i];
        potential_bucket = m_index[i] >= 0;
      }

      // std::cout <<"m_min = "<<m_min<< "m_index = "<<m_index<<" m_quadrant =
      // "<<m_quadrant << std::endl;

      // if index is outside domain don't bother calcing
      // distance
      if (potential_bucket) {

        const double accum = get_min_distance_to_bucket(m_index);
        // std::cout <<"accum = "<<accum<< std::endl;

        potential_bucket = accum <= m_max_distance2;
      }

      // if bucket is still good break out of loop
      if (potential_bucket)
        break;

      // must be outside distance or outside domain, so reset index back to min
      m_index[i] = m_min[i];

      // if gone through all dimensions, move to next quadrant
      if (i == 0) {
        ++m_quadrant;
        if (m_quadrant < (1 << dimension)) {
          reset_min_and_index();
        } else {
          // if gone through all quadrants, iterator
          // is now invalid
          m_valid = false;
        }
      }
    }
    LOG_CUDA(3, "lattice_iterator_within_distance: increment :end");
  }
}; // namespace Aboria

} // namespace Aboria

#endif /* BUCKETSEARCH_H_ */
