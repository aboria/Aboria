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

#ifndef BUCKETSEARCH_H_
#define BUCKETSEARCH_H_

#include "CudaInclude.h"
#include "Get.h"
#include "NeighbourSearchBase.h"
#include "SpatialUtil.h"
#include "Traits.h"
#include "Vector.h"
#include "detail/Algorithms.h"

#include "Log.h"
#include <iostream>

namespace Aboria {

template <typename Traits> struct CellListOrdered_params {
  typedef typename Traits::double_d double_d;
  CellListOrdered_params() : side_length(detail::get_max<double>()) {}
  CellListOrdered_params(const double_d &side_length)
      : side_length(side_length) {}
  double_d side_length;
};

template <typename Traits> struct CellListOrderedQuery;

/// @brief A cell list spatial data structure that is paired with a
/// CellListOrderedQuery query type
///
/// This class implements neighbourhood searching using a cell list. That is,
/// the domain is divided up into a regular grid of constant size
/// "buckets", and particles are assigned to their containing bucket.
///
/// The main difference between this class and CellList is that the particle set
/// is reordered according to which cell each particle belongs to. This
/// correlates memory locality with spatial locality, and improves cache access
/// if you are often looping over neighbouring particles. It also mean that
/// particles within a given bucket are sequential in memory.
///
///
template <typename Traits>
class CellListOrdered
    : public neighbour_search_base<CellListOrdered<Traits>, Traits,
                                   CellListOrderedQuery<Traits>> {

  typedef typename Traits::double_d double_d;
  typedef typename Traits::position position;
  typedef typename Traits::vector_double_d_const_iterator
      vector_double_d_const_iterator;
  typedef typename Traits::vector_unsigned_int_iterator
      vector_unsigned_int_iterator;
  typedef typename Traits::vector_unsigned_int vector_unsigned_int;
  typedef typename Traits::unsigned_int_d unsigned_int_d;
  typedef typename Traits::iterator iterator;
  typedef CellListOrdered_params<Traits> params_type;

  typedef neighbour_search_base<CellListOrdered<Traits>, Traits,
                                CellListOrderedQuery<Traits>>
      base_type;

  friend base_type;

public:
  CellListOrdered()
      : base_type(),
        m_size_calculated_with_n(std::numeric_limits<size_t>::max()) {}

  static constexpr bool ordered() { return true; }

  struct delete_points_lambda;

  void print_data_structure() const {
#ifndef __CUDA_ARCH__
    LOG(1, "\tbuckets:");
    for (size_t i = 0; i < m_bucket_begin.size(); ++i) {
      LOG(1, "\ti = " << i << " bucket contents = " << m_bucket_begin[i]
                      << " to " << m_bucket_end[i]);
    }
    LOG(1, "\tend buckets");
    LOG(1, "\tparticles:");
    for (size_t i = 0; i < m_bucket_indices.size(); ++i) {
      LOG(1, "\ti = " << i << " p = "
                      << static_cast<const double_d &>(
                             get<position>(*(this->m_particles_begin + i)))
                      << " bucket = " << m_bucket_indices[i]);
    }
    LOG(1, "\tend particles:");
#endif
  }

private:
  bool set_domain_impl() {
    const size_t n = this->m_alive_indices.size();
    if (n < 0.5 * m_size_calculated_with_n ||
        n > 2 * m_size_calculated_with_n) {
      LOG(2, "CellListOrdered: recalculating bucket size");
      m_size_calculated_with_n = n;
      if (this->m_n_particles_in_leaf > n) {
        m_size = unsigned_int_d::Constant(1);
      } else {
        const double total_volume =
            (this->m_bounds.bmax - this->m_bounds.bmin).prod();
        const double box_volume =
            this->m_n_particles_in_leaf / double(n) * total_volume;
        const double box_side_length =
            std::pow(box_volume, 1.0 / Traits::dimension);
        m_size =
            floor((this->m_bounds.bmax - this->m_bounds.bmin) / box_side_length)
                .template cast<unsigned int>();
        for (size_t i = 0; i < Traits::dimension; ++i) {
          if (m_size[i] == 0) {
            m_size[i] = 1;
          }
        }
      }
      m_bucket_side_length =
          (this->m_bounds.bmax - this->m_bounds.bmin) / m_size;
      m_point_to_bucket_index =
          detail::point_to_bucket_index<Traits::dimension>(
              m_size, m_bucket_side_length, this->m_bounds);

      LOG(2, "\tbucket side length = " << m_bucket_side_length);
      LOG(2, "\tnumber of buckets = " << m_size << " (total=" << m_size.prod()
                                      << ")");

      // setup bucket data structures
      m_bucket_begin.resize(m_size.prod());
      m_bucket_end.resize(m_size.prod());

      this->m_query.m_bucket_begin =
          iterator_to_raw_pointer(m_bucket_begin.begin());
      this->m_query.m_bucket_end =
          iterator_to_raw_pointer(m_bucket_end.begin());
      this->m_query.m_nbuckets = m_bucket_begin.size();

      this->m_query.m_bucket_side_length = m_bucket_side_length;
      this->m_query.m_bounds.bmin = this->m_bounds.bmin;
      this->m_query.m_bounds.bmax = this->m_bounds.bmax;
      this->m_query.m_periodic = this->m_periodic;
      this->m_query.m_end_bucket = m_size - 1;
      this->m_query.m_point_to_bucket_index = m_point_to_bucket_index;
      return true;
    } else {
      return false;
    }
  }

  void update_iterator_impl() {}

  void update_positions_impl(iterator update_begin, iterator update_end,
                             const int new_n,
                             const bool call_set_domain = true) {

    ASSERT(update_begin == this->m_particles_begin &&
               update_end == this->m_particles_end,
           "error should be update all");

    if (call_set_domain) {
      set_domain_impl();
    }

    const size_t n = this->m_alive_indices.size();
    m_bucket_indices.resize(n);
    if (n > 0) {
      // transform the points to their bucket indices
      if (static_cast<size_t>(update_end - update_begin) == n) {
        // m_alive_indicies is just a sequential list of indices
        // (i.e. no dead)
        detail::transform(get<position>(this->m_particles_begin) +
                              this->m_alive_indices[0],
                          get<position>(this->m_particles_begin) +
                              this->m_alive_indices[0] + n,
                          m_bucket_indices.begin(), m_point_to_bucket_index);
      } else {
        // m_alive_indicies contains all alive indicies
        detail::transform(Traits::make_permutation_iterator(
                              get<position>(this->m_particles_begin),
                              this->m_alive_indices.begin()),
                          Traits::make_permutation_iterator(
                              get<position>(this->m_particles_begin),
                              this->m_alive_indices.end()),
                          m_bucket_indices.begin(), m_point_to_bucket_index);
      }

      // sort the points by their bucket index
      detail::sort_by_key(m_bucket_indices.begin(), m_bucket_indices.end(),
                          this->m_alive_indices.begin());
    }

    // find the beginning of each bucket's list of points
    auto search_begin = Traits::make_counting_iterator(0);
    detail::lower_bound(m_bucket_indices.begin(), m_bucket_indices.end(),
                        search_begin, search_begin + m_size.prod(),
                        m_bucket_begin.begin());

    // find the end of each bucket's list of points
    detail::upper_bound(m_bucket_indices.begin(), m_bucket_indices.end(),
                        search_begin, search_begin + m_size.prod(),
                        m_bucket_end.begin());

#ifndef __CUDA_ARCH__
    if (4 <= ABORIA_LOG_LEVEL) {
      LOG(4, "\tbuckets:");
      for (size_t i = 0; i < m_bucket_begin.size(); ++i) {
        LOG(4, "\ti = " << i << " bucket contents = " << m_bucket_begin[i]
                        << " to " << m_bucket_end[i]);
      }
      LOG(4, "\tend buckets");
      LOG(4, "\tparticles:");
      for (size_t i = 0; i < m_bucket_indices.size(); ++i) {
        LOG(4, "\ti = " << i << " p = "
                        << static_cast<const double_d &>(
                               get<position>(*(this->m_particles_begin + i)))
                        << " bucket = " << m_bucket_indices[i]);
      }
      LOG(4, "\tend particles:");
    }
#endif
  }

  const CellListOrderedQuery<Traits> &get_query_impl() const { return m_query; }

  CellListOrderedQuery<Traits> &get_query_impl() { return m_query; }

  // the grid data structure keeps a range per grid bucket:
  // each bucket_begin[i] indexes the first element of bucket i's list of points
  // each bucket_end[i] indexes one past the last element of bucket i's list of
  // points
  vector_unsigned_int m_bucket_begin;
  vector_unsigned_int m_bucket_end;
  vector_unsigned_int m_bucket_indices;
  CellListOrderedQuery<Traits> m_query;

  double_d m_bucket_side_length;
  unsigned_int_d m_size;
  size_t m_size_calculated_with_n;
  detail::point_to_bucket_index<Traits::dimension> m_point_to_bucket_index;
};

/// @copydetails NeighbourQueryBase
///
/// @brief This is a query object for the CellListOrdered spatial data structure
///
template <typename Traits>
struct CellListOrderedQuery : public NeighbourQueryBase<Traits> {

  typedef Traits traits_type;
  typedef typename Traits::raw_pointer raw_pointer;
  typedef typename Traits::double_d double_d;
  typedef typename Traits::bool_d bool_d;
  typedef typename Traits::int_d int_d;
  typedef typename Traits::unsigned_int_d unsigned_int_d;
  const static unsigned int dimension = Traits::dimension;
  template <int LNormNumber, typename Transform = IdentityTransform>
  using query_iterator =
      lattice_iterator_within_distance<CellListOrderedQuery, LNormNumber,
                                       Transform>;
  typedef lattice_iterator<dimension> all_iterator;
  typedef lattice_iterator<dimension> child_iterator;
  typedef typename query_iterator<2>::reference reference;
  typedef typename query_iterator<2>::pointer pointer;
  typedef typename query_iterator<2>::value_type value_type;
  typedef ranges_iterator<Traits> particle_iterator;
  typedef bbox<dimension> box_type;

  ///
  /// @brief pointer to the beginning of the particle set
  ///
  raw_pointer m_particles_begin;

  ///
  /// @brief pointer to the end of the particle set
  ///
  raw_pointer m_particles_end;

  ///
  /// @brief periodicity of the domain
  ///
  bool_d m_periodic;

  ///
  /// @brief dimensions of each bucket
  ///
  double_d m_bucket_side_length;

  ///
  /// @brief index of the last bucket in the cell list
  ///
  int_d m_end_bucket;

  ///
  /// @brief min/max bounds of the domain
  ///
  bbox<dimension> m_bounds;

  ///
  /// @brief function object to transform a point to a bucket index
  ///
  detail::point_to_bucket_index<dimension> m_point_to_bucket_index;

  ///
  /// @brief pointer to the beginning of the buckets
  ///
  unsigned int *m_bucket_begin;

  ///
  /// @brief pointer to the end of the buckets
  ///
  unsigned int *m_bucket_end;

  ///
  /// @brief the number of buckets
  ///
  unsigned int m_nbuckets;

  ///
  /// @brief a pointer to the "key" values of the find-by-id map
  ///
  size_t *m_id_map_key;

  ///
  /// @brief a pointer to the "value" values of the find-by-id map
  ///
  size_t *m_id_map_value;

  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  CellListOrderedQuery() {}

  /*
   * functions for id mapping
   */

  ///
  /// @copydoc NeighbourQueryBase::find()
  ///
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
   * functions for trees
   */

  ///
  /// @copydoc NeighbourQueryBase::is_leaf_node()
  ///
  /// always true for CellListOrdered
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  static bool is_leaf_node(const value_type &bucket) { return true; }

  ///
  /// @copydoc NeighbourQueryBase::is_tree()
  ///
  /// always false for CellListOrdered
  ///
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  static bool is_tree() { return false; }

  ///
  /// @copydoc NeighbourQueryBase::get_children() const
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  child_iterator get_children() const {
    return child_iterator(int_d::Constant(0), m_end_bucket + 1);
  }

  ///
  /// @copydoc NeighbourQueryBase::get_children(const child_iterator&) const
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  child_iterator get_children(const child_iterator &ci) const {
    return child_iterator();
  }

  ///
  /// @copydoc NeighbourQueryBase::num_children(const child_iterator&) const
  ///
  static size_t num_children(const child_iterator &ci) { return 0; }

  ///
  /// @copydoc NeighbourQueryBase::num_children() const
  ///
  size_t num_children() const { return number_of_buckets(); }

  ///
  /// @copydoc NeighbourQueryBase::get_bounds()
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  const box_type get_bounds(const child_iterator &ci) const {
    box_type bounds;
    bounds.bmin = (*ci) * m_bucket_side_length + m_bounds.bmin;
    bounds.bmax = ((*ci) + 1) * m_bucket_side_length + m_bounds.bmin;
    return bounds;
  }

  ///
  /// @copydoc NeighbourQueryBase::get_bounds()
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  const box_type &get_bounds() const { return m_bounds; }

  ///
  /// @copydoc NeighbourQueryBase::get_periodic()
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  const bool_d &get_periodic() const { return m_periodic; }

  ///
  /// @copydoc NeighbourQueryBase::get_bucket_particles()
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  particle_iterator get_bucket_particles(const reference bucket) const {
#ifndef __CUDA_ARCH__
    ASSERT((bucket >= int_d::Constant(0)).all() &&
               (bucket <= m_end_bucket).all(),
           "invalid bucket");
#endif

    const unsigned int bucket_index =
        m_point_to_bucket_index.collapse_index_vector(bucket);
    const unsigned int range_start_index = m_bucket_begin[bucket_index];
    const unsigned int range_end_index = m_bucket_end[bucket_index];

#ifndef __CUDA_ARCH__
    LOG(4, "\tlooking in bucket "
               << bucket << " = " << bucket_index << ". found "
               << range_end_index - range_start_index << " particles");
#endif
    return particle_iterator(m_particles_begin + range_start_index,
                             m_particles_begin + range_end_index);
  }

  ///
  /// @copydoc NeighbourQueryBase::get_bucket_bbox()
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  bbox<dimension> get_bucket_bbox(const reference bucket) const {
    return bbox<dimension>(bucket * m_bucket_side_length + m_bounds.bmin,
                           (bucket + 1) * m_bucket_side_length + m_bounds.bmin);
  }

  ///
  /// @copydoc NeighbourQueryBase::get_bucket()
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  child_iterator get_bucket(const double_d &position) const {
    auto bucket = m_point_to_bucket_index.find_bucket_index_vector(position);
    return child_iterator(bucket, bucket + 1);
  }

  ///
  /// @copydoc NeighbourQueryBase::get_bucket_index()
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  size_t get_bucket_index(const reference bucket) const {
    return m_point_to_bucket_index.collapse_index_vector(bucket);
  }

  ///
  /// @copydoc NeighbourQueryBase::get_buckets_near_point()
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  template <int LNormNumber, typename Transform = IdentityTransform>
  CUDA_HOST_DEVICE query_iterator<LNormNumber, Transform>
  get_buckets_near_point(const double_d &position, const double max_distance,
                         const Transform &transform = Transform()) const {
#ifndef __CUDA_ARCH__
    LOG(4, "\tget_buckets_near_point: position = "
               << position << " max_distance = " << max_distance);
#endif
    return query_iterator<LNormNumber, Transform>(position, max_distance, this,
                                                  transform);
  }

  ///
  /// @copydoc NeighbourQueryBase::get_end_bucket()
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  const int_d &get_end_bucket() const { return m_end_bucket; }

  ///
  /// @copydoc NeighbourQueryBase::get_subtree(const child_iterator&) const
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  all_iterator get_subtree(const child_iterator &ci) const {
    return all_iterator();
  }

  ///
  /// @copydoc NeighbourQueryBase::get_subtree() const
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  all_iterator get_subtree() const {
    return all_iterator(int_d::Constant(0), m_end_bucket + 1);
  }

  ///
  /// @copydoc NeighbourQueryBase::number_of_buckets()
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  size_t number_of_buckets() const { return (m_end_bucket + 1).prod(); }

  ///
  /// @copydoc NeighbourQueryBase::number_of_particles()
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  size_t number_of_particles() const {
    return (m_particles_end - m_particles_begin);
  }

  ///
  /// @copydoc NeighbourQueryBase::get_particles_begin() const
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  const raw_pointer &get_particles_begin() const { return m_particles_begin; }

  ///
  /// @copydoc NeighbourQueryBase::get_particles_begin()
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  raw_pointer &get_particles_begin() { return m_particles_begin; }

  ///
  /// @copydoc NeighbourQueryBase::number_of_levels()
  ///
  /// always 2 for CellListOrdered
  ///
  unsigned number_of_levels() const { return 2; }
};

} // namespace Aboria

#endif /* BUCKETSEARCH_H_ */
