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

#ifndef BUCKETSEARCH_SERIAL_H_
#define BUCKETSEARCH_SERIAL_H_

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

template <typename Traits> struct CellListQuery;

/// @brief A cell list spatial data structure that is paired with a
/// CellListQuery query type
///
/// This class implements neighbourhood searching using a bucket search, or
/// cell list algorithm. The domain is first divided up into a regular grid of
/// constant size "buckets".
///
/// After the buckets are created, a set of 3D points can be assigned to their
/// respective buckets. After this,
/// neighbourhood queries around a given point can be performed be looking
/// at all the points in the same bucket or surrounding buckets of the given
/// point.
///
template <typename Traits>
class CellList : public neighbour_search_base<CellList<Traits>, Traits,
                                              CellListQuery<Traits>> {

  typedef typename Traits::double_d double_d;
  typedef typename Traits::position position;
  typedef typename Traits::vector_int vector_int;
  typedef typename Traits::iterator iterator;
  typedef typename Traits::unsigned_int_d unsigned_int_d;

  typedef neighbour_search_base<CellList<Traits>, Traits, CellListQuery<Traits>>
      base_type;

  friend base_type;

public:
  ///
  /// @brief constructs the default base structure
  ///
  CellList()
      : base_type(),
        m_size_calculated_with_n(std::numeric_limits<size_t>::max()),
        m_serial(detail::concurrent_processes<Traits>() == 1) {}

  ///
  /// @brief This structure is not ordered. That is, the order of the particles
  ///        is not important to the search algorithm. Instead, the indicies of
  ///        the particles are kept in a linked list
  ///
  static constexpr bool ordered() { return false; }

  struct delete_points_in_bucket_lambda;
  struct insert_points_lambda_sequential_serial;
  struct insert_points_lambda_non_sequential_serial;
  struct insert_points_lambda_sequential;
  struct insert_points_lambda_non_sequential;
  struct copy_points_in_bucket_lambda;

  ///
  /// @brief Print the data structure to stdout
  ///
  void print_data_structure() const {
#ifndef __CUDA_ARCH__
    LOG(1, "\tbuckets:");
    for (size_t i = 0; i < m_buckets.size(); ++i) {
      if (m_buckets[i] != detail::get_empty_id()) {
        LOG(1, "\ti = " << i << " bucket contents = " << m_buckets[i]);
      }
    }
    LOG(1, "\tend buckets");
    LOG(1, "\tlinked list:");
    for (size_t i = 0; i < m_linked_list.size(); ++i) {
      if (m_serial) {
        LOG(1, "\ti = " << i << " p = "
                        << static_cast<const double_d &>(
                               get<position>(*(this->m_particles_begin + i)))
                        << " contents = " << m_linked_list[i]
                        << ". reverse = " << m_linked_list_reverse[i]);
      } else {
        LOG(1, "\ti = " << i << " p = "
                        << static_cast<const double_d &>(
                               get<position>(*(this->m_particles_begin + i)))
                        << " contents = " << m_linked_list[i]);
      }
    }
    LOG(1, "\tend linked list:");
#endif
  }

private:
  ///
  /// @brief set the domain, called from base class
  /// @return true if this function wishes to trigger a reconstruction of the
  /// data structure (i.e. resize the buckets and re-insert all the particles to
  /// the linked list) set
  ///
  bool set_domain_impl() {
    const size_t n = this->m_particles_end - this->m_particles_begin;
    if (n < 0.5 * m_size_calculated_with_n ||
        n > 2 * m_size_calculated_with_n) {
      m_size_calculated_with_n = n;
      LOG(2, "CellList: recalculating bucket size");
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

      m_buckets.assign(m_size.prod(), detail::get_empty_id());

      // TODO: should always be true?
      m_use_dirty_cells = true;

      this->m_query.m_buckets_begin =
          iterator_to_raw_pointer(m_buckets.begin());

      this->m_query.m_bucket_side_length = this->m_bucket_side_length;
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

  ///
  /// @brief check that the data structure is internally consistent
  ///
  void check_data_structure() {
    int num_particles = 0;
    for (size_t i = 0; i < m_buckets.size(); i++) {
      int j = m_buckets[i];
      int old_j = detail::get_empty_id();
      while (j != detail::get_empty_id()) {
        ASSERT(m_linked_list_reverse[j] == old_j,
               "m_linked_list_reverse not right: m_linked_list_reverse[j] = "
                   << m_linked_list_reverse[j] << ", j = " << j
                   << ", old_j = " << old_j);
        ASSERT(m_dirty_buckets[j] == i,
               "m_dirty_buckets not right: m_dirty_buckets[j] = "
                   << m_dirty_buckets[j] << ", i = " << i << ", j = " << j);
        ASSERT(j < m_linked_list.size(), "j index too large");
        ASSERT(j >= 0, "j index less than zero");
        num_particles++;
        old_j = j;
        j = m_linked_list[old_j];
      }
    }
    ASSERT(num_particles == m_linked_list.size(),
           "m_linked_list size inconsistent");
    ASSERT(num_particles == m_linked_list_reverse.size(),
           "m_linked_list_reverse size inconsistent");
    ASSERT(num_particles == m_dirty_buckets.size(),
           "m_dirty_buckets size inconsistent");
    for (size_t i = 0; i < m_linked_list.size(); ++i) {
      ASSERT(i < m_linked_list.size(), "i index too large");
      ASSERT(m_dirty_buckets[i] < m_buckets.size(),
             "m_dirty_buckets not right");
      ASSERT(m_dirty_buckets[i] >= 0, "m_dirty_buckets not right");
      ASSERT(m_linked_list[i] < m_linked_list.size() ||
                 m_linked_list[i] == detail::get_empty_id(),
             "m_linked_list not right: m_linked_list[i] = "
                 << m_linked_list[i]
                 << ", m_linked_list.size() = " << m_linked_list.size());
      ASSERT(m_linked_list[i] >= 0 ||
                 m_linked_list[i] == detail::get_empty_id(),
             "m_linked_list not right");
      ASSERT(m_linked_list_reverse[i] < m_linked_list_reverse.size() ||
                 m_linked_list_reverse[i] == detail::get_empty_id(),
             "m_linked_list_reverse not right");
      ASSERT(m_linked_list_reverse[i] >= 0 ||
                 m_linked_list_reverse[i] == detail::get_empty_id(),
             "m_linked_list_reverse not right");
    }
  }

  ///
  /// @brief called by base class, does nothing for this data structure
  ///
  void update_iterator_impl() {
    // check_data_structure();
  }

  ///
  /// @brief called by base class, deals with an updated particle range.

  /// Needs to handle:
  ///     dead particles (only if @p update_end == @p m_particles_end),
  ///     particles with different positions,
  ///     and new particles
  ///
  /// @param update_begin iterator to the beginning of the update range
  /// @param update_end iterator to the end of the update range
  /// @param num_new_particles_added number of new particles to be added
  /// @param call_set_domain set to true (the default) if set_domain_impl()
  /// needs to be called
  ///
  void update_positions_impl(iterator update_begin, iterator update_end,
                             const int num_new_particles_added,
                             const bool call_set_domain = true) {
    // if call_set_domain == false then set_domain_impl() has already
    // been called, and returned true
    const bool reset_domain = call_set_domain ? set_domain_impl() : true;
    const size_t n_update = update_end - update_begin;
    const size_t n_alive = this->m_alive_indices.size();
    const size_t n_dead_in_update = n_update - n_alive;
    const size_t n_all = this->m_particles_end - this->m_particles_begin;
    const size_t n = n_all - n_dead_in_update;
    const int update_begin_index = update_begin - this->m_particles_begin;

    LOG(2, "BucketSearchSerial: update_positions, n_update = "
               << n_update << " n_alive = " << n_alive << " n = " << n);

    if (n_update == n_all || reset_domain) {
      // updating everthing so clear out entire ds
      if (!reset_domain) {
        if (m_dirty_buckets.size() <
            static_cast<float>(m_buckets.size()) /
                detail::concurrent_processes<Traits>()) {
          for (int i : m_dirty_buckets) {
            m_buckets[i] = detail::get_empty_id();
          }
        } else {
          m_buckets.assign(m_buckets.size(), detail::get_empty_id());
        }
      }
      m_linked_list.assign(n, detail::get_empty_id());
      if (m_serial) {
        m_linked_list_reverse.assign(n, detail::get_empty_id());
      }
      m_dirty_buckets.resize(n);
    } else {
      // only updating some so only clear out indices in update range
      const int start_index_deleted = update_begin - this->m_particles_begin;
      const int end_index_deleted =
          update_end - this->m_particles_begin - num_new_particles_added;
      if (m_serial) {
        for (int toi = start_index_deleted; toi < end_index_deleted; ++toi) {
          const int toi_back = m_linked_list_reverse[toi];
          const int toi_forward = m_linked_list[toi];
          ASSERT(
              toi_back == detail::get_empty_id() ||
                  (toi_back < static_cast<int>(m_linked_list_reverse.size()) &&
                   toi_back >= 0),
              "invalid index of " << toi_back
                                  << ". Note linked list reverse size is "
                                  << m_linked_list_reverse.size());
          ASSERT(toi_forward == detail::get_empty_id() ||
                     (toi_forward <
                          static_cast<int>(m_linked_list_reverse.size()) &&
                      toi_forward >= 0),
                 "invalid index of " << toi_back
                                     << ". Note linked list reverse size is "
                                     << m_linked_list_reverse.size());
          if (toi_back != detail::get_empty_id()) {
            m_linked_list[toi_back] = toi_forward;
          } else {
            m_buckets[m_dirty_buckets[toi]] = toi_forward;
          }
          if (toi_forward != detail::get_empty_id()) {
            m_linked_list_reverse[toi_forward] = toi_back;
          }
        }
      } else {
        m_deleted_buckets.resize(end_index_deleted - start_index_deleted);
        detail::copy(m_dirty_buckets.begin() + start_index_deleted,
                     m_dirty_buckets.begin() + end_index_deleted,
                     m_deleted_buckets.begin());
        detail::sort(m_deleted_buckets.begin(), m_deleted_buckets.end());
        detail::for_each(
            m_deleted_buckets.begin(),
            detail::unique(m_deleted_buckets.begin(), m_deleted_buckets.end()),
            delete_points_in_bucket_lambda(
                start_index_deleted, end_index_deleted,
                iterator_to_raw_pointer(m_linked_list.begin()),
                iterator_to_raw_pointer(m_buckets.begin())));
      }
      // increase capacity of linked list vectors
      m_linked_list.resize(n, detail::get_empty_id());
      if (m_serial) {
        m_linked_list_reverse.resize(n, detail::get_empty_id());
      }
      m_dirty_buckets.resize(n);
    }

    if (reset_domain) {
      // resetting domain so need to insert all particles
      insert_points(0, update_begin_index);
    }
    // then insert points that are still alive within update range
    if (n_dead_in_update == 0) {
      insert_points(update_begin_index, n_update);
    } else {
      insert_points(this->m_alive_indices.begin(), this->m_alive_indices.end(),
                    update_begin_index);
    }

#ifndef __CUDA_ARCH__
    if (4 <= ABORIA_LOG_LEVEL) {
      LOG(4, "\tbuckets:");
      for (size_t i = 0; i < m_buckets.size(); ++i) {
        if (m_buckets[i] != detail::get_empty_id()) {
          LOG(4, "\ti = " << i << " bucket contents = " << m_buckets[i]);
        }
      }
      LOG(4, "\tend buckets");
      LOG(4, "\tlinked list:");
      for (size_t i = 0; i < m_linked_list.size(); ++i) {
        if (m_serial) {
          LOG(4, "\ti = " << i << " p = "
                          << static_cast<const double_d &>(
                                 get<position>(*(this->m_particles_begin + i)))
                          << " contents = " << m_linked_list[i]
                          << ". reverse = " << m_linked_list_reverse[i]);
        } else {
          LOG(4, "\ti = " << i << " p = "
                          << static_cast<const double_d &>(
                                 get<position>(*(this->m_particles_begin + i)))
                          << " contents = " << m_linked_list[i]);
        }
      }
      LOG(4, "\tend linked list:");
    }
#endif

    // check_data_structure();

    this->m_query.m_linked_list_begin =
        iterator_to_raw_pointer(this->m_linked_list.begin());
  }

  ///
  /// @brief insert non-consecutive points into data structure
  ///
  /// @param start_adding iterator to the start of the list of indices to be
  /// inserted
  /// @param stop_adding iterator to the end of the list of indices to be
  /// inserted
  /// @param start the index of the particle set to start adding (the start of
  /// the update range)
  ///
  void insert_points(typename vector_int::iterator start_adding,
                     typename vector_int::iterator stop_adding,
                     const int start) {
    const int n = stop_adding - start_adding;
#if defined(__CUDACC__)
    typedef typename thrust::detail::iterator_category_to_system<
        typename vector_int::iterator::iterator_category>::type system;
    thrust::counting_iterator<int, system> count(0);
#else
    auto count = Traits::make_counting_iterator(0);
#endif

    if (m_serial) { // running in serial
      detail::for_each(
          count, count + n,
          insert_points_lambda_non_sequential_serial(
              iterator_to_raw_pointer(get<position>(this->m_particles_begin)),
              iterator_to_raw_pointer(start_adding), m_point_to_bucket_index,
              iterator_to_raw_pointer(m_buckets.begin()),
              iterator_to_raw_pointer(m_dirty_buckets.begin()),
              iterator_to_raw_pointer(m_linked_list.begin()),
              iterator_to_raw_pointer(m_linked_list_reverse.begin()), start));
    } else { // running in parallel
      detail::for_each(
          count, count + n,
          insert_points_lambda_non_sequential(
              iterator_to_raw_pointer(get<position>(this->m_particles_begin)),
              iterator_to_raw_pointer(start_adding), m_point_to_bucket_index,
              iterator_to_raw_pointer(m_buckets.begin()),
              iterator_to_raw_pointer(m_dirty_buckets.begin()),
              iterator_to_raw_pointer(m_linked_list.begin()), start));
    }
  }

  void insert_points(const int start, const int n) {
#if defined(__CUDACC__)
    typedef typename thrust::detail::iterator_category_to_system<
        typename vector_int::iterator::iterator_category>::type system;
    thrust::counting_iterator<int, system> count(0);
#else
    auto count = Traits::make_counting_iterator(0);
#endif

    if (m_serial) { // running in serial
      detail::for_each(
          count, count + n,
          insert_points_lambda_sequential_serial(
              iterator_to_raw_pointer(get<position>(this->m_particles_begin)),
              m_point_to_bucket_index,
              iterator_to_raw_pointer(m_buckets.begin()),
              iterator_to_raw_pointer(m_dirty_buckets.begin()),
              iterator_to_raw_pointer(m_linked_list.begin()),
              iterator_to_raw_pointer(m_linked_list_reverse.begin()), start));
    } else { // running in parallel

      detail::for_each(
          count, count + n,
          insert_points_lambda_sequential(
              iterator_to_raw_pointer(get<position>(this->m_particles_begin)),
              m_point_to_bucket_index,
              iterator_to_raw_pointer(m_buckets.begin()),
              iterator_to_raw_pointer(m_dirty_buckets.begin()),
              iterator_to_raw_pointer(m_linked_list.begin()), start));
    }
  }

  ///
  /// @brief returns the query object via the base class
  ///
  const CellListQuery<Traits> &get_query_impl() const { return m_query; }

  ///
  /// @brief returns the query object via the base class
  ///
  CellListQuery<Traits> &get_query_impl() { return m_query; }

  ///
  /// @brief vector of the beginning nodes for each bucket linked list
  ///
  vector_int m_buckets;

  ///
  /// @brief vector of the linked list nodes for each particle
  ///
  vector_int m_linked_list;

  ///
  /// @brief vector of the backwards linked list nodes for each particle
  ///
  vector_int m_linked_list_reverse;

  ///
  /// @brief vector of which buckets have particles in them
  ///
  vector_int m_dirty_buckets;

  ///
  /// @brief internal storage for detecting which buckets have particles to be
  /// deleted
  ///
  vector_int m_deleted_buckets;

  ///
  /// @brief internal storage for detecting which buckets have particles to be
  /// copied
  ///
  vector_int m_copied_buckets;

  ///
  /// @brief the query object
  ///
  CellListQuery<Traits> m_query;

  ///
  /// @brief use m_dirty_buckets for efficiency
  ///
  bool m_use_dirty_cells;

  ///
  /// @brief last resizing of the buckets occurred with n particles
  ///
  size_t m_size_calculated_with_n;

  ///
  /// @brief running with multiple processes, or on gpu
  ///
  bool m_serial;

  ///
  /// @brief how many buckets in each dimension
  ///
  unsigned_int_d m_size;

  ///
  /// @brief length of buckets in each dimension
  ///
  double_d m_bucket_side_length;

  ///
  /// @brief struct to convert a point to a bucket index
  ///
  detail::point_to_bucket_index<Traits::dimension> m_point_to_bucket_index;
};

///
/// @brief non-threadsafe function object to insert a non-consecutive list of
/// particles into
///  the data structure
///
/// @tparam Traits The @ref TraitsCommon type
///
template <typename Traits>
struct CellList<Traits>::insert_points_lambda_non_sequential_serial {
  typedef typename Traits::double_d double_d;
  typedef typename detail::point_to_bucket_index<Traits::dimension> ptobl_type;
  double_d *m_positions;
  int *m_alive_indices;
  ptobl_type m_point_to_bucket_index;
  int *m_buckets;
  int *m_dirty_buckets;
  int *m_linked_list;
  int *m_linked_list_reverse;
  int start;

  ///
  /// @brief copy all the neccessary info into the function object
  ///
  insert_points_lambda_non_sequential_serial(
      double_d *m_positions, int *m_alive_indices,
      const ptobl_type &m_point_to_bucket_index, int *m_buckets,
      int *m_dirty_buckets, int *m_linked_list, int *m_linked_list_reverse,
      int start)
      : m_positions(m_positions), m_alive_indices(m_alive_indices),
        m_point_to_bucket_index(m_point_to_bucket_index), m_buckets(m_buckets),
        m_dirty_buckets(m_dirty_buckets), m_linked_list(m_linked_list),
        m_linked_list_reverse(m_linked_list_reverse), start(start) {}

  ///
  /// @brief insert a particle at index @p i into the data structure
  ///
  /// It is assumed that this is run in serial (not thread-safe)
  CUDA_HOST_DEVICE
  void operator()(const int i) {
    // use actual index to insert into ds
    const int new_index = i + start;
    // use m_alive_index to get position
    const double_d &r = m_positions[m_alive_indices[i]];
    const unsigned int bucketi = m_point_to_bucket_index.find_bucket_index(r);
    const int bucket_entry = m_buckets[bucketi];

    // Insert into own cell
    m_buckets[bucketi] = new_index;
    m_dirty_buckets[new_index] = bucketi;
    m_linked_list[new_index] = bucket_entry;
    m_linked_list_reverse[new_index] = detail::get_empty_id();
    if (bucket_entry != detail::get_empty_id())
      m_linked_list_reverse[bucket_entry] = new_index;
  }
};

///
/// @brief non-threadsafe insert a consecutive vector of particles into the data
/// structure
///
/// @tparam Traits the @ref TraitsCommon type
///
template <typename Traits>
struct CellList<Traits>::insert_points_lambda_sequential_serial {
  typedef typename Traits::double_d double_d;
  typedef typename detail::point_to_bucket_index<Traits::dimension> ptobl_type;
  double_d *m_positions;
  ptobl_type m_point_to_bucket_index;
  int *m_buckets;
  int *m_dirty_buckets;
  int *m_linked_list;
  int *m_linked_list_reverse;
  int start;

  ///
  /// @brief copy in all the required info
  ///
  insert_points_lambda_sequential_serial(
      double_d *m_positions, const ptobl_type &m_point_to_bucket_index,
      int *m_buckets, int *m_dirty_buckets, int *m_linked_list,
      int *m_linked_list_reverse, int start)
      : m_positions(m_positions),
        m_point_to_bucket_index(m_point_to_bucket_index), m_buckets(m_buckets),
        m_dirty_buckets(m_dirty_buckets), m_linked_list(m_linked_list),
        m_linked_list_reverse(m_linked_list_reverse), start(start) {}

  ///
  /// @brief insert a particle at index @p i into the data structure
  ///
  /// It is assumed that this is run in serial (not thread-safe)
  CUDA_HOST_DEVICE
  void operator()(const int i) {
    // use actual index to insert into ds
    const int new_index = i + start;
    const double_d &r = m_positions[new_index];
    const unsigned int bucketi = m_point_to_bucket_index.find_bucket_index(r);
    ASSERT_CUDA(bucketi < m_point_to_bucket_index.m_bucket_index.m_size.prod());
    // std::cout << "inserting particle in index "<<new_index<<" at "<<r << "
    // into bucket "<<bucketi<<std::endl;
    const int bucket_entry = m_buckets[bucketi];

    // Insert into own cell
    m_buckets[bucketi] = new_index;
    m_dirty_buckets[new_index] = bucketi;
    m_linked_list[new_index] = bucket_entry;
    m_linked_list_reverse[new_index] = detail::get_empty_id();
    if (bucket_entry != detail::get_empty_id())
      m_linked_list_reverse[bucket_entry] = new_index;
  }
};

///
/// @brief thread-safe function object to insert a list of non-consecutive
/// particles into the data structure
///
/// @tparam Traits the @ref TraitsCommon type
///
template <typename Traits>
struct CellList<Traits>::insert_points_lambda_non_sequential {
  typedef typename Traits::double_d double_d;
  typedef typename detail::point_to_bucket_index<Traits::dimension> ptobl_type;
  double_d *m_positions;
  int *m_alive_indices;
  ptobl_type m_point_to_bucket_index;
  int *m_buckets;
  int *m_dirty_buckets;
  int *m_linked_list;
  int start;

  ///
  /// @brief copy all the required info
  ///
  insert_points_lambda_non_sequential(double_d *m_positions,
                                      int *m_alive_indices,
                                      const ptobl_type &m_point_to_bucket_index,
                                      int *m_buckets, int *m_dirty_buckets,
                                      int *m_linked_list, int start)
      : m_positions(m_positions), m_alive_indices(m_alive_indices),
        m_point_to_bucket_index(m_point_to_bucket_index), m_buckets(m_buckets),
        m_dirty_buckets(m_dirty_buckets), m_linked_list(m_linked_list),
        start(start) {}

  ///
  /// @brief insert a particle with index @p i into the data structure
  ///
  /// implements a lock-free linked list using atomic CAS
  ///
  CUDA_HOST_DEVICE
  void operator()(const int i) {
    // if (!m_alive[i]) return;

    const int new_index = i + start;
    const unsigned int bucketi = m_point_to_bucket_index.find_bucket_index(
        m_positions[m_alive_indices[i]]);

// printf("insert_points_lambda: i = %d, bucketi = %d",i,bucketi);

// try inserting at head of the list
#if defined(__CUDA_ARCH__)
    int next;
    do {
      next = m_buckets[bucketi];
    } while (atomicCAS(m_buckets + bucketi, next, new_index) != new_index);
#else
    int next = m_buckets[bucketi];
    while (!__atomic_compare_exchange_n(m_buckets + bucketi, &next, new_index,
                                        true, __ATOMIC_RELAXED,
                                        __ATOMIC_RELAXED))
      ;
#endif

    // successful
    m_dirty_buckets[new_index] = bucketi;
    m_linked_list[new_index] = next;

    // m_linked_list_reverse[i] = detail::get_empty_id();
    // if (next != detail::get_empty_id()) m_linked_list_reverse[next] = i;
  }
};

///
/// @brief thread-safe function object to insert a list of consecutive
/// particles into the data structure
///
/// @tparam Traits the @ref TraitsCommon type
///
template <typename Traits>
struct CellList<Traits>::insert_points_lambda_sequential {
  typedef typename Traits::double_d double_d;
  typedef typename detail::point_to_bucket_index<Traits::dimension> ptobl_type;
  double_d *m_positions;
  ptobl_type m_point_to_bucket_index;
  int *m_buckets;
  int *m_dirty_buckets;
  int *m_linked_list;
  int start;

  ///
  /// @brief copy all the info
  ///
  insert_points_lambda_sequential(double_d *m_positions,
                                  const ptobl_type &m_point_to_bucket_index,
                                  int *m_buckets, int *m_dirty_buckets,
                                  int *m_linked_list, int start)
      : m_positions(m_positions),
        m_point_to_bucket_index(m_point_to_bucket_index), m_buckets(m_buckets),
        m_dirty_buckets(m_dirty_buckets), m_linked_list(m_linked_list),
        start(start) {}

  ///
  /// @brief insert the particle at index @p i
  ///
  /// implements a lock-free linked list using atomic CAS
  ///
  CUDA_HOST_DEVICE
  void operator()(const int i) {
    // if (!m_alive[i]) return;

    const int new_index = i + start;
    const unsigned int bucketi =
        m_point_to_bucket_index.find_bucket_index(m_positions[new_index]);

// printf("insert_points_lambda: i = %d, bucketi = %d",i,bucketi);

// try inserting at head of the list
#if defined(__CUDA_ARCH__)
    int next;
    do {
      next = m_buckets[bucketi];
    } while (atomicCAS(m_buckets + bucketi, next, new_index) != new_index);
#else
    int next = m_buckets[bucketi];
    while (!__atomic_compare_exchange_n(m_buckets + bucketi, &next, new_index,
                                        true, __ATOMIC_RELAXED,
                                        __ATOMIC_RELAXED))
      ;
#endif

    // successful
    m_dirty_buckets[new_index] = bucketi;
    m_linked_list[new_index] = next;

    // m_linked_list_reverse[i] = detail::get_empty_id();
    // if (next != detail::get_empty_id()) m_linked_list_reverse[next] = i;
  }
};

///
/// @brief function object to delete a range of particles within the buckets
///
/// @tparam Traits the @ref TraitsCommon type
///
template <typename Traits>
struct CellList<Traits>::delete_points_in_bucket_lambda {
  int *m_linked_list;
  int *m_buckets;
  int start_index_deleted;
  int end_index_deleted;

  delete_points_in_bucket_lambda(int start_index_deleted, int end_index_deleted,
                                 int *m_linked_list, int *m_buckets)
      : m_linked_list(m_linked_list), m_buckets(m_buckets),
        start_index_deleted(start_index_deleted),
        end_index_deleted(end_index_deleted) {}

  ///
  /// @brief goes through all the particles in bucket @p celli and deletes
  ///  particles within the range given by m_start_index_deleted < i <
  ///  m_end_index_deleted
  ///
  CUDA_HOST_DEVICE
  void operator()(const int celli) {
    // go through linked list
    int i_minus_1 = detail::get_empty_id();
    int i = m_buckets[celli];
    while (i != detail::get_empty_id()) {
      // unlink each contiguous deleted range
      if (i >= start_index_deleted && i < end_index_deleted) {
        // get first forward index not in range
        do {
          i = m_linked_list[i];
        } while (i != detail::get_empty_id() && i >= start_index_deleted &&
                 i < end_index_deleted);

        // update links
        if (i_minus_1 != detail::get_empty_id()) {
          m_linked_list[i_minus_1] = i;
        } else {
          m_buckets[celli] = i;
        }
        if (i == detail::get_empty_id()) {
          break;
        }
      }
      i_minus_1 = i;
      i = m_linked_list[i];
    }
  }
};

///
/// @brief function object to copy a range of particles to another index range
///
/// @tparam Traits the @ref TraitsCommon type
///
template <typename Traits>
struct CellList<Traits>::copy_points_in_bucket_lambda {
  int *m_linked_list;
  int *m_buckets;
  size_t start_index_deleted;
  size_t start_index_copied;
  size_t end_index_copied;

  copy_points_in_bucket_lambda(size_t start_index_deleted,
                               size_t start_index_copied,
                               size_t end_index_copied, int *m_linked_list,
                               int *m_buckets)
      : start_index_deleted(start_index_deleted),
        start_index_copied(start_index_copied),
        end_index_copied(end_index_copied), m_linked_list(m_linked_list),
        m_buckets(m_buckets) {}

  ///
  /// @brief goes through bucket @p celli and updates the particle indices to
  /// point to the new range
  ///
  CUDA_HOST_DEVICE
  void operator()(const int celli) {
    // go through linked list
    int i_minus_1 = detail::get_empty_id();
    int i = m_buckets[celli];
    while (i != detail::get_empty_id()) {
      // update each copied index
      if (i >= start_index_copied && i < end_index_copied) {
        int i_plus_1 = m_linked_list[i];
        i = i - start_index_copied + start_index_deleted;
        if (i_minus_1 != detail::get_empty_id()) {
          m_linked_list[i_minus_1] = i;
        } else {
          m_buckets[celli] = i;
        }
        m_linked_list[i] = i_plus_1;
      }
      i_minus_1 = i;
      i = m_linked_list[i];
    }
  }
};

/// @copydetails NeighbourQueryBase
///
/// @brief This is a query object for the CellList spatial data structure
///
template <typename Traits>
struct CellListQuery : public NeighbourQueryBase<Traits> {

  typedef Traits traits_type;
  typedef typename Traits::raw_pointer raw_pointer;
  typedef typename Traits::double_d double_d;
  typedef typename Traits::bool_d bool_d;
  typedef typename Traits::int_d int_d;
  typedef typename Traits::unsigned_int_d unsigned_int_d;
  typedef typename Traits::reference particle_reference;
  typedef typename Traits::const_reference particle_const_reference;
  const static unsigned int dimension = Traits::dimension;
  template <int LNormNumber, typename Transform = IdentityTransform>
  using query_iterator =
      lattice_iterator_within_distance<CellListQuery, LNormNumber, Transform>;

  typedef lattice_iterator<dimension> all_iterator;
  typedef lattice_iterator<dimension> child_iterator;
  typedef typename query_iterator<2>::reference reference;
  typedef typename query_iterator<2>::pointer pointer;
  typedef typename query_iterator<2>::value_type value_type;
  typedef linked_list_iterator<Traits> particle_iterator;
  typedef bbox<dimension> box_type;

  ///
  /// @brief periodicity of domain
  ///
  bool_d m_periodic;

  ///
  /// @brief bucket length in each dimension
  ///
  double_d m_bucket_side_length;

  ///
  /// @brief index of last bucket
  ///
  int_d m_end_bucket;

  ///
  /// @brief domain min/max bounds
  ///
  bbox<dimension> m_bounds;

  ///
  /// @brief function object to calculate a bucket index from a position
  ///
  detail::point_to_bucket_index<dimension> m_point_to_bucket_index;

  ///
  /// @brief pointer to the beginning of the particle set
  ///
  raw_pointer m_particles_begin;

  ///
  /// @brief pointer to the end of the particle set
  ///
  raw_pointer m_particles_end;

  ///
  /// @brief pointer to the beginning of the buckets
  ///
  int *m_buckets_begin;

  ///
  /// @brief pointer to the beginning of the linked list
  ///
  int *m_linked_list_begin;

  ///
  /// @brief pointer to the find-by-id map key
  ///
  size_t *m_id_map_key;

  ///
  /// @brief pointer to the find-by-id map value
  ///
  size_t *m_id_map_value;

  ///
  /// @brief constructor checks that we are not using std::vector and cuda
  /// at the same time
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  CellListQuery() {
#if defined(__CUDA_ARCH__)
    CHECK_CUDA((!std::is_same<typename Traits::template vector<double>,
                              std::vector<double>>::value),
               "Cannot use std::vector in device code");
#endif
  }

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
  /// CellList implements a "flat" tree, so all buckets are
  /// leafs
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  static bool is_leaf_node(const value_type &bucket) { return true; }

  ///
  /// @copydoc NeighbourQueryBase::is_tree()
  ///
  /// CellList is not a proper tree structure, so will return false always
  ///
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  static bool is_tree() { return false; }

  ///
  /// @copydoc NeighbourQueryBase::get_children()
  ///
  /// for CellList this is all the buckets (flat tree with one root
  /// node, all the rest leafs)
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

#ifndef __CUDA_ARCH__
    LOG(4, "\tget_bucket_particles: looking in bucket " << bucket << " = "
                                                        << bucket_index);
#endif
    return particle_iterator(m_buckets_begin[bucket_index], m_particles_begin,
                             m_linked_list_begin);
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
  /// @see lattice_iterator_within_distance
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
  /// @copydoc NeighbourQueryBase::get_particles_begin()
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
  unsigned number_of_levels() const { return 2; }
};

} // namespace Aboria

#endif /* BUCKETSEARCH_H_ */
