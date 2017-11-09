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

#include <boost/iterator/iterator_facade.hpp>
#include "Traits.h"
#include "Vector.h"
#include "Log.h"
#include "Get.h"
#include "NeighbourSearchBase.h"
#include "detail/SpatialUtil.h"
#include <vector>
#include <iostream>
#include <set>


namespace Aboria {

template <typename Traits>
struct bucket_search_serial_params {
    typedef typename Traits::double_d double_d;
    bucket_search_serial_params(): 
        side_length(detail::get_max<double>()) {}
    bucket_search_serial_params(const double_d& side_length):
        side_length(side_length) {}
    double_d side_length;
};

template <typename Traits>
class bucket_search_serial_query; 

/// \brief Implements neighbourhood searching using a bucket search algorithm, dividing
/// the domain into constant size "buckets".
///
/// This class implements neighbourhood searching using a bucket search algorithm. The 
/// domain is first divided up into a regular grid of constant size "buckets", either by
/// using the class constructor to initialise the domain extents, bucket size etc., or by 
/// using the reset() member function to reset these parameters.
///
/// After the buckets are created, a set of 3D points can be assigned to their respective buckets
/// using the embed_points() member function. After this, neighbourhood queries around
/// a given point can be performed using find_broadphase_neighbours(), which returns a const 
/// iterator to all the points in the same bucket or surrounding buckets of the given point.
///
template <typename Traits>
class bucket_search_serial: 
    public neighbour_search_base<bucket_search_serial<Traits>,
                                 Traits,
                                 bucket_search_serial_query<Traits>> {

    typedef typename Traits::double_d double_d;
    typedef typename Traits::position position;
    typedef typename Traits::vector_int vector_int;
    typedef typename Traits::iterator iterator;
    typedef typename Traits::unsigned_int_d unsigned_int_d;

    typedef bucket_search_serial_params<Traits> params_type;

    typedef neighbour_search_base<bucket_search_serial<Traits>,
                                 Traits,
                                 bucket_search_serial_query<Traits>> base_type;

    friend base_type;

public:
    bucket_search_serial():
        m_size_calculated_with_n(std::numeric_limits<size_t>::max()),
        m_serial(detail::concurrent_processes<Traits>() == 1),
        base_type() {}

    static constexpr bool ordered() {
        return false;
    }

    struct delete_points_in_bucket_lambda;
    struct insert_points_lambda_sequential_serial;
    struct insert_points_lambda_non_sequential_serial;
    struct insert_points_lambda_sequential;
    struct insert_points_lambda_non_sequential;
    struct copy_points_in_bucket_lambda;

    void print_data_structure() const {
        #ifndef __CUDA_ARCH__
            LOG(1,"\tbuckets:");
            for (int i = 0; i<m_buckets.size(); ++i) {
                if (m_buckets[i] != detail::get_empty_id()) {
                    LOG(1,"\ti = "<<i<<" bucket contents = "<<m_buckets[i]);
                }
            }
            LOG(1,"\tend buckets");
            LOG(1,"\tlinked list:");
            for (int i = 0; i<m_linked_list.size(); ++i) {
                if (m_serial) {
                    LOG(1,"\ti = "<<i<<" p = "<<
                        static_cast<const double_d&>(get<position>(*(this->m_particles_begin+i)))<<
                        " contents = "<<m_linked_list[i]<<". reverse = "<<m_linked_list_reverse[i]);
                } else {
                     LOG(1,"\ti = "<<i<<" p = "<<
                        static_cast<const double_d&>(get<position>(*(this->m_particles_begin+i)))<<
                        " contents = "<<m_linked_list[i]);

                }
            }
            LOG(1,"\tend linked list:");
#endif
    }



private:
    bool set_domain_impl() {
        const size_t n = this->m_particles_end - this->m_particles_begin;
        if (n < 0.5*m_size_calculated_with_n || n > 2*m_size_calculated_with_n) {
            m_size_calculated_with_n = n;
            LOG(2,"bucket_search_serial: recalculating bucket size");
            if (this->m_n_particles_in_leaf > n) {
                m_size = unsigned_int_d(1);
            } else {
                const double total_volume = (this->m_bounds.bmax-this->m_bounds.bmin).prod();
                const double box_volume = this->m_n_particles_in_leaf/double(n)*total_volume;
                const double box_side_length = std::pow(box_volume,1.0/Traits::dimension);
                m_size = 
                    floor((this->m_bounds.bmax-this->m_bounds.bmin)/box_side_length)
                    .template cast<unsigned int>();
                for (int i=0; i<Traits::dimension; ++i) {
                    if (m_size[i] == 0) {
                        m_size[i] = 1;
                    }
                }
            }
            m_bucket_side_length = (this->m_bounds.bmax-this->m_bounds.bmin)/m_size;
            m_point_to_bucket_index = 
                detail::point_to_bucket_index<Traits::dimension>(m_size,m_bucket_side_length,this->m_bounds);

            LOG(2,"\tbucket side length = "<<m_bucket_side_length);
            LOG(2,"\tnumber of buckets = "<<m_size<<" (total="<<m_size.prod()<<")");

            m_buckets.assign(m_size.prod(), detail::get_empty_id());

            if (!m_serial) { 
                m_buckets_begin.resize(m_size.prod());
                m_buckets_end.resize(m_size.prod());
            }

            //TODO: should always be true?
            m_use_dirty_cells = true;

            this->m_query.m_buckets_begin = iterator_to_raw_pointer(m_buckets.begin());

            this->m_query.m_bucket_side_length = this->m_bucket_side_length;
            this->m_query.m_bounds.bmin = this->m_bounds.bmin;
            this->m_query.m_bounds.bmax = this->m_bounds.bmax;
            this->m_query.m_periodic = this->m_periodic;
            this->m_query.m_end_bucket = m_size-1;
            this->m_query.m_point_to_bucket_index = m_point_to_bucket_index;
            return true;
        } else {
            return false;
        }
    }

    void check_data_structure() {
        int num_particles = 0;
        for (int i=0; i<m_buckets.size(); i++) {
            int j = m_buckets[i];
            int old_j = detail::get_empty_id();
            while (j != detail::get_empty_id()) {
                ASSERT(m_linked_list_reverse[j] == old_j, "m_linked_list_reverse not right: m_linked_list_reverse[j] = "<<m_linked_list_reverse[j]<<", j = "<<j<<", old_j = "<<old_j); 
                ASSERT(m_dirty_buckets[j] == i, "m_dirty_buckets not right: m_dirty_buckets[j] = "<<m_dirty_buckets[j]<<", i = "<<i<<", j = "<<j); 
                ASSERT(j < m_linked_list.size(), "j index too large");
                ASSERT(j >= 0, "j index less than zero");
                num_particles++;
                old_j = j;
                j = m_linked_list[old_j];
            }
        }
        ASSERT(num_particles == m_linked_list.size(), "m_linked_list size inconsistent");
        ASSERT(num_particles == m_linked_list_reverse.size(), "m_linked_list_reverse size inconsistent");
        ASSERT(num_particles == m_dirty_buckets.size(), "m_dirty_buckets size inconsistent");
        for (int i=0; i<m_linked_list.size(); ++i) {
            ASSERT(i < m_linked_list.size(), "i index too large");
            ASSERT(i >= 0, "i index less than zero");
            ASSERT(m_dirty_buckets[i] < m_buckets.size(), "m_dirty_buckets not right");
            ASSERT(m_dirty_buckets[i] >= 0, "m_dirty_buckets not right");
            ASSERT(m_linked_list[i] < m_linked_list.size() || m_linked_list[i] == detail::get_empty_id(), "m_linked_list not right: m_linked_list[i] = "<<m_linked_list[i]<<", m_linked_list.size() = "<<m_linked_list.size());
            ASSERT(m_linked_list[i] >= 0 || m_linked_list[i] == detail::get_empty_id(), "m_linked_list not right");
            ASSERT(m_linked_list_reverse[i] < m_linked_list_reverse.size() || m_linked_list_reverse[i] == detail::get_empty_id(), "m_linked_list_reverse not right");
            ASSERT(m_linked_list_reverse[i] >= 0 || m_linked_list_reverse[i] == detail::get_empty_id(), "m_linked_list_reverse not right");
        }
    }

    /*
    void end_list_of_copies_impl() {
        int n = this->m_particles_end-this->m_particles_begin;
        ASSERT(n <= m_linked_list.size(),"particle size should not be greater than linked list size");
        // any resizing should be due to "copy_points"
        if (n < m_linked_list.size()) {
            LOG(3,"bucket_search_serial: resizing arrays in update_iterator")
            const bool resize_buckets = set_domain_impl();
            if (resize_buckets) {
                // buckets all changed, so start from scratch
                LOG(3,"bucket_search_serial: re_embed points")
                embed_points_impl(false);
            } else {
                m_linked_list.resize(n);
                if (m_serial) m_linked_list_reverse.resize(n);
                m_dirty_buckets.resize(n);
                this->m_query.m_linked_list_begin = iterator_to_raw_pointer(
                        this->m_linked_list.begin());
            }
        }
        #ifndef __CUDA_ARCH__
        if (4 <= ABORIA_LOG_LEVEL) { 
            LOG(4,"\tbuckets:");
            for (int i = 0; i<m_buckets.size(); ++i) {
                if (m_buckets[i] != detail::get_empty_id()) {
                    LOG(4,"\ti = "<<i<<" bucket contents = "<<m_buckets[i]);
                }
            }
            LOG(4,"\tend buckets");
            LOG(4,"\tlinked list:");
            for (int i = 0; i<m_linked_list.size(); ++i) {
                if (m_serial) {
                    LOG(4,"\ti = "<<i<<" p = "<<
                        static_cast<const double_d&>(get<position>(*(this->m_particles_begin+i)))<<
                        " contents = "<<m_linked_list[i]<<". reverse = "<<m_linked_list_reverse[i]);
                } else {
                     LOG(4,"\ti = "<<i<<" p = "<<
                        static_cast<const double_d&>(get<position>(*(this->m_particles_begin+i)))<<
                        " contents = "<<m_linked_list[i]);

                }
            }
            LOG(4,"\tend linked list:");
        }
        #endif


    }
    */

    void update_iterator_impl() {
        //check_data_structure();
    }


    // need to handle dead particles (do not insert into ds)
    void update_positions_impl(iterator update_begin, iterator update_end,
                               const int num_new_particles_added,
                               const bool call_set_domain=true) {
        // if call_set_domain == false then set_domain_impl() has already
        // been called, and returned true
        const bool reset_domain = call_set_domain ? set_domain_impl() : true;
        const size_t n_update = update_end-update_begin;
        const size_t n_alive = this->m_alive_indices.size();
        const size_t n_dead_in_update = n_update-n_alive;
        const size_t n_all = this->m_particles_end-this->m_particles_begin;
        const size_t n = n_all-n_dead_in_update;
        const int update_begin_index = update_begin-this->m_particles_begin;
        const int update_end_index = update_end-this->m_particles_begin;

        LOG(2,"BucketSearchSerial: update_positions, n_update = "<<n_update<<" n_alive = "<<n_alive<<" n = "<<n);

        if (n_update == n_all || reset_domain) {
            // updating everthing so clear out entire ds
            if (!reset_domain) {
                if (m_dirty_buckets.size() < static_cast<float>(m_buckets.size())
                        / detail::concurrent_processes<Traits>()) {
                    for (int i: m_dirty_buckets) {
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
            const int start_index_deleted = update_begin-this->m_particles_begin;
            const int end_index_deleted = update_end-this->m_particles_begin-num_new_particles_added;
            if (m_serial) {
                for (int toi = start_index_deleted; toi < end_index_deleted; ++toi) {
                    const int toi_back = m_linked_list_reverse[toi];
                    const int toi_forward = m_linked_list[toi];
                    ASSERT(toi_back == detail::get_empty_id() ||
                            (toi_back < m_linked_list_reverse.size() && toi_back >= 0),
                            "invalid index of "<<toi_back <<". Note linked list reverse size is "
                            << m_linked_list_reverse.size());
                    ASSERT(toi_forward == detail::get_empty_id() ||
                            (toi_forward < m_linked_list_reverse.size() && toi_forward >= 0),
                            "invalid index of "<<toi_back <<". Note linked list reverse size is "
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
                m_deleted_buckets.resize(end_index_deleted-start_index_deleted);
                detail::copy(m_dirty_buckets.begin()+start_index_deleted,
                        m_dirty_buckets.begin()+end_index_deleted,
                        m_deleted_buckets.begin());
                detail::sort(m_deleted_buckets.begin(),m_deleted_buckets.end());
                detail::for_each(m_deleted_buckets.begin(),
                        detail::unique(m_deleted_buckets.begin(),
                            m_deleted_buckets.end()),
                        delete_points_in_bucket_lambda(start_index_deleted,
                                                       end_index_deleted,
                            iterator_to_raw_pointer(m_linked_list.begin()),
                            iterator_to_raw_pointer(m_buckets.begin())
                            ));
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
            insert_points(0,update_begin_index);
        }
        // then insert points that are still alive within update range
        if (n_dead_in_update == 0) {
            insert_points(update_begin_index,n_update);
        } else {
            insert_points(this->m_alive_indices.begin(),
                      this->m_alive_indices.end(),
                      update_begin_index); 
        }
        
#ifndef __CUDA_ARCH__
        if (4 <= ABORIA_LOG_LEVEL) { 
            LOG(4,"\tbuckets:");
            for (int i = 0; i<m_buckets.size(); ++i) {
                if (m_buckets[i] != detail::get_empty_id()) {
                    LOG(4,"\ti = "<<i<<" bucket contents = "<<m_buckets[i]);
                }
            }
            LOG(4,"\tend buckets");
            LOG(4,"\tlinked list:");
            for (int i = 0; i<m_linked_list.size(); ++i) {
                if (m_serial) {
                    LOG(4,"\ti = "<<i<<" p = "<<
                        static_cast<const double_d&>(get<position>(*(this->m_particles_begin+i)))<<
                        " contents = "<<m_linked_list[i]<<". reverse = "<<m_linked_list_reverse[i]);
                } else {
                     LOG(4,"\ti = "<<i<<" p = "<<
                        static_cast<const double_d&>(get<position>(*(this->m_particles_begin+i)))<<
                        " contents = "<<m_linked_list[i]);

                }
            }
            LOG(4,"\tend linked list:");
        }
#endif

        //check_data_structure();

        this->m_query.m_linked_list_begin = iterator_to_raw_pointer(this->m_linked_list.begin());
    }

    void insert_points(typename vector_int::iterator start_adding, 
                       typename vector_int::iterator stop_adding,
                       const int start) {
        const int n = stop_adding-start_adding;
#if defined(__CUDACC__)
        typedef typename thrust::detail::iterator_category_to_system<
            typename vector_int::iterator::iterator_category
            >::type system;
        detail::counting_iterator<int,system> count(0);
#else
        detail::counting_iterator<int> count(0);
#endif

        if (m_serial) { // running in serial
            detail::for_each(count,count+n,
                    insert_points_lambda_non_sequential_serial(
                        iterator_to_raw_pointer(
                            get<position>(this->m_particles_begin)),
                        iterator_to_raw_pointer(start_adding),
                        m_point_to_bucket_index,
                        iterator_to_raw_pointer(m_buckets.begin()),
                        iterator_to_raw_pointer(m_dirty_buckets.begin()),
                        iterator_to_raw_pointer(m_linked_list.begin()),
                        iterator_to_raw_pointer(m_linked_list_reverse.begin()),
                        start));
        } else { // running in parallel
            detail::for_each(count,count+n,
                    insert_points_lambda_non_sequential(
                        iterator_to_raw_pointer(
                            get<position>(this->m_particles_begin)),
                        iterator_to_raw_pointer(start_adding),
                        m_point_to_bucket_index,
                        iterator_to_raw_pointer(m_buckets.begin()),
                        iterator_to_raw_pointer(m_dirty_buckets.begin()),
                        iterator_to_raw_pointer(m_linked_list.begin()),
                        start));
        }
    }


    void insert_points(const int start,const int n) {
#if defined(__CUDACC__)
        typedef typename thrust::detail::iterator_category_to_system<
            typename vector_int::iterator::iterator_category
            >::type system;
        detail::counting_iterator<int,system> count(0);
#else
        detail::counting_iterator<int> count(0);
#endif

        if (m_serial) { // running in serial
            auto raw_positions = iterator_to_raw_pointer(
                                    get<position>(this->m_particles_begin));
            detail::for_each(count,count+n,
                    insert_points_lambda_sequential_serial(
                        iterator_to_raw_pointer(
                            get<position>(this->m_particles_begin)),
                        m_point_to_bucket_index,
                        iterator_to_raw_pointer(m_buckets.begin()),
                        iterator_to_raw_pointer(m_dirty_buckets.begin()),
                        iterator_to_raw_pointer(m_linked_list.begin()),
                        iterator_to_raw_pointer(m_linked_list_reverse.begin()),
                        start));
        } else { // running in parallel

            detail::for_each(count,count+n,
                    insert_points_lambda_sequential(
                        iterator_to_raw_pointer(
                            get<position>(this->m_particles_begin)),
                        m_point_to_bucket_index,
                        iterator_to_raw_pointer(m_buckets.begin()),
                        iterator_to_raw_pointer(m_dirty_buckets.begin()),
                        iterator_to_raw_pointer(m_linked_list.begin()),
                        start));
        }
    }

    const bucket_search_serial_query<Traits>& get_query_impl() const {
        return m_query;
    }

    bucket_search_serial_query<Traits>& get_query_impl() {
        return m_query;
    }

    vector_int m_buckets;
    vector_int m_buckets_begin;
    vector_int m_buckets_end;
    vector_int m_linked_list;
    vector_int m_linked_list_reverse;
    vector_int m_dirty_buckets;
    vector_int m_deleted_buckets;
    vector_int m_copied_buckets;
    bucket_search_serial_query<Traits> m_query;
    bool m_use_dirty_cells;

    size_t m_size_calculated_with_n;
    bool m_serial;
    unsigned_int_d m_size;
    double_d m_bucket_side_length;
    detail::point_to_bucket_index<Traits::dimension> m_point_to_bucket_index;

};

template <typename Traits>
struct bucket_search_serial<Traits>::insert_points_lambda_non_sequential_serial {
    typedef typename Traits::double_d double_d;
    typedef typename detail::point_to_bucket_index<Traits::dimension> ptobl_type;
    double_d* m_positions;
    int* m_alive_indices;
    int* m_buckets;
    int* m_dirty_buckets;
    int* m_linked_list;
    int* m_linked_list_reverse;
    int start;
    ptobl_type m_point_to_bucket_index;


    insert_points_lambda_non_sequential_serial(double_d* m_positions,
                         int* m_alive_indices,
                         const ptobl_type& m_point_to_bucket_index, 
                         int* m_buckets,
                         int* m_dirty_buckets,
                         int* m_linked_list,
                         int* m_linked_list_reverse,
                         int start
                         ):
        m_positions(m_positions),
        m_alive_indices(m_alive_indices),
        m_point_to_bucket_index(m_point_to_bucket_index),
        m_buckets(m_buckets),
        m_dirty_buckets(m_dirty_buckets),
        m_linked_list(m_linked_list),
        m_linked_list_reverse(m_linked_list_reverse),
        start(start)
    {}

    // implements a lock-free linked list using atomic cas
    CUDA_HOST_DEVICE
    void operator()(const int i) {
        // use actual index to insert into ds
        const int new_index = i+start;
        // use m_alive_index to get position
        const double_d& r = m_positions[m_alive_indices[i]];
        const unsigned int bucketi = m_point_to_bucket_index.find_bucket_index(r);
        const int bucket_entry = m_buckets[bucketi];

        // Insert into own cell
        m_buckets[bucketi] = new_index;
        m_dirty_buckets[new_index] = bucketi;
        m_linked_list[new_index] = bucket_entry;
        m_linked_list_reverse[new_index] = detail::get_empty_id();
        if (bucket_entry != detail::get_empty_id()) m_linked_list_reverse[bucket_entry] = new_index;
    }
};

template <typename Traits>
struct bucket_search_serial<Traits>::insert_points_lambda_sequential_serial {
    typedef typename Traits::double_d double_d;
    typedef typename detail::point_to_bucket_index<Traits::dimension> ptobl_type;
    double_d* m_positions;
    int* m_buckets;
    int* m_dirty_buckets;
    int* m_linked_list;
    int* m_linked_list_reverse;
    int start;
    ptobl_type m_point_to_bucket_index;


    insert_points_lambda_sequential_serial(double_d* m_positions,
                         const ptobl_type& m_point_to_bucket_index, 
                         int* m_buckets,
                         int* m_dirty_buckets,
                         int* m_linked_list,
                         int* m_linked_list_reverse,
                         int start
                         ):
        m_positions(m_positions),
        m_point_to_bucket_index(m_point_to_bucket_index),
        m_buckets(m_buckets),
        m_dirty_buckets(m_dirty_buckets),
        m_linked_list(m_linked_list),
        m_linked_list_reverse(m_linked_list_reverse),
        start(start)
    {}

    // implements a lock-free linked list using atomic cas
    CUDA_HOST_DEVICE
    void operator()(const int i) {
        // use actual index to insert into ds
        const int new_index = i+start;
        const double_d& r = m_positions[new_index];
        const unsigned int bucketi = m_point_to_bucket_index.find_bucket_index(r);
        //std::cout << "inserting particle in index "<<new_index<<" at "<<r << " into bucket "<<bucketi<<std::endl;
        const int bucket_entry = m_buckets[bucketi];

        // Insert into own cell
        m_buckets[bucketi] = new_index;
        m_dirty_buckets[new_index] = bucketi;
        m_linked_list[new_index] = bucket_entry;
        m_linked_list_reverse[new_index] = detail::get_empty_id();
        if (bucket_entry != detail::get_empty_id()) m_linked_list_reverse[bucket_entry] = new_index;

    }
};




template <typename Traits>
struct bucket_search_serial<Traits>::insert_points_lambda_non_sequential {
    typedef typename Traits::double_d double_d;
    typedef typename detail::point_to_bucket_index<Traits::dimension> ptobl_type;
    double_d* m_positions;
    int* m_alive_indices;
    int* m_buckets;
    int* m_dirty_buckets;
    int* m_linked_list;
    int start;
    ptobl_type m_point_to_bucket_index;


    insert_points_lambda_non_sequential(double_d* m_positions,
                         int* m_alive_indices,
                         const ptobl_type& m_point_to_bucket_index, 
                         int* m_buckets,
                         int* m_dirty_buckets,
                         int* m_linked_list,
                         int start
                         ):
        m_positions(m_positions),
        m_alive_indices(m_alive_indices),
        m_point_to_bucket_index(m_point_to_bucket_index),
        m_buckets(m_buckets),
        m_dirty_buckets(m_dirty_buckets),
        m_linked_list(m_linked_list),
        start(start)
    {}

    // implements a lock-free linked list using atomic cas
    CUDA_HOST_DEVICE
    void operator()(const int i) {
        //if (!m_alive[i]) return;

        const int new_index = i+start;
        const unsigned int bucketi = m_point_to_bucket_index.find_bucket_index(m_positions[m_alive_indices[i]]);

        //printf("insert_points_lambda: i = %d, bucketi = %d",i,bucketi);

        //try inserting at head of the list
        #if defined(__CUDA_ARCH__)
        int next;
        do {
            next = m_buckets[bucketi];
        } while (atomicCAS(m_buckets + bucketi, next, new_index) != new_index);
        #else
        int next = m_buckets[bucketi];
        while (!__atomic_compare_exchange_n(m_buckets + bucketi, &next, new_index, 
                                    true, __ATOMIC_RELAXED, __ATOMIC_RELAXED));
        #endif

        //successful
        m_dirty_buckets[new_index] = bucketi;
        m_linked_list[new_index] = next;

        //m_linked_list_reverse[i] = detail::get_empty_id();
        //if (next != detail::get_empty_id()) m_linked_list_reverse[next] = i;
    }
};

template <typename Traits>
struct bucket_search_serial<Traits>::insert_points_lambda_sequential {
    typedef typename Traits::double_d double_d;
    typedef typename detail::point_to_bucket_index<Traits::dimension> ptobl_type;
    double_d* m_positions;
    int* m_buckets;
    int* m_dirty_buckets;
    int* m_linked_list;
    int start;
    ptobl_type m_point_to_bucket_index;


    insert_points_lambda_sequential(double_d* m_positions,
                         const ptobl_type& m_point_to_bucket_index, 
                         int* m_buckets,
                         int* m_dirty_buckets,
                         int* m_linked_list,
                         int start
                         ):
        m_positions(m_positions),
        m_point_to_bucket_index(m_point_to_bucket_index),
        m_buckets(m_buckets),
        m_dirty_buckets(m_dirty_buckets),
        m_linked_list(m_linked_list),
        start(start)
    {}

    // implements a lock-free linked list using atomic cas
    CUDA_HOST_DEVICE
    void operator()(const int i) {
        //if (!m_alive[i]) return;

        const int new_index = i+start;
        const unsigned int bucketi = m_point_to_bucket_index.find_bucket_index(m_positions[new_index]);

        //printf("insert_points_lambda: i = %d, bucketi = %d",i,bucketi);

        //try inserting at head of the list
        #if defined(__CUDA_ARCH__)
        int next;
        do {
            next = m_buckets[bucketi];
        } while (atomicCAS(m_buckets + bucketi, next, new_index) != new_index);
        #else
        int next = m_buckets[bucketi];
        while (!__atomic_compare_exchange_n(m_buckets + bucketi, &next, new_index, 
                                    true, __ATOMIC_RELAXED, __ATOMIC_RELAXED));
        #endif

        //successful
        m_dirty_buckets[new_index] = bucketi;
        m_linked_list[new_index] = next;

        //m_linked_list_reverse[i] = detail::get_empty_id();
        //if (next != detail::get_empty_id()) m_linked_list_reverse[next] = i;
    }
};




template <typename Traits>
struct bucket_search_serial<Traits>::delete_points_in_bucket_lambda {
    int *m_linked_list;
    int *m_buckets;
    size_t start_index_deleted;
    size_t end_index_deleted;

    delete_points_in_bucket_lambda(size_t start_index_deleted,
                         size_t end_index_deleted,
                         int* m_linked_list,
                         int* m_buckets):
        m_linked_list(m_linked_list),
        m_buckets(m_buckets),
        start_index_deleted(start_index_deleted),
        end_index_deleted(end_index_deleted)
    {}

    CUDA_HOST_DEVICE
    void operator()(const int celli) {
        // go through linked list
        int i_minus_1 = detail::get_empty_id();
        int i = m_buckets[celli];
        while (i != detail::get_empty_id()) {
            // unlink each contiguous deleted range
            if (i >= start_index_deleted && i < end_index_deleted) {
                //get first forward index not in range
                do {
                    i = m_linked_list[i];
                } while (i != detail::get_empty_id() && 
                         i >= start_index_deleted 
                         && i < end_index_deleted);

                //update links
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

template <typename Traits>
struct bucket_search_serial<Traits>::copy_points_in_bucket_lambda {
    int* m_linked_list;
    int* m_buckets;
    size_t start_index_deleted;
    size_t start_index_copied;
    size_t end_index_copied;

    copy_points_in_bucket_lambda(size_t start_index_deleted,
                       size_t start_index_copied,
                       size_t end_index_copied,
                       int* m_linked_list,
                       int* m_buckets):
        start_index_deleted(start_index_deleted),
        start_index_copied(start_index_copied),
        end_index_copied(end_index_copied),
        m_linked_list(m_linked_list),
        m_buckets(m_buckets)
    {}

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

/*
template <typename Traits>
struct bucket_search_serial<Traits>::copy_points_lambda {
    int* m_linked_list;
    int* m_buckets;
    raw_pointer m_particles_begin;
    int m_n;

    copy_points_lambda(size_t start_index_deleted,
                       size_t start_index_copied,
                       size_t end_index_copied,
                       int* m_linked_list,
                       int* m_buckets):
        start_index_deleted(start_index_deleted),
        start_index_copied(start_index_copied),
        end_index_copied(end_index_copied),
        m_linked_list(m_linked_list),
        m_buckets(m_buckets)
    {}

    CUDA_HOST_DEVICE
    void copy_points_per_particle(const size_t to, const size_t from) {
    }

    CUDA_HOST_DEVICE
    void copy_points_per_bucket(const size_t to, const size_t from) {
        // go through from linked list
        int i_minus_1 = detail::get_empty_id();
        int i = m_buckets[m_dirty_buckets[from]];
        while (i != detail::get_empty_id()) {
            // update each copied index
            if (i == from) {
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

        // go through to linked list
        i_minus_1 = detail::get_empty_id();
        i = m_buckets[m_dirty_buckets[to]];
        while (i != detail::get_empty_id()) {
            // delete to index
            if (i == to) {
                int i_plus_1 = m_linked_list[i];
                m_linked_list[i_minus_1] = i_plus_1;
            }
            i_minus_1 = i;
            i = m_linked_list[i];
        }
    }
        

    CUDA_HOST_DEVICE
    void operator()(reference i) {
        if (!get<alive>(i)) {
            const size_t index = &get<position>(i) 
                - get<position>(m_particles_begin);

            // do search for next alive index
            size_t delete_index = index;
            do {
                delete_index = m_n-m_delete_indicies[delete_index];
                if (delete_index <= index) return;
            } while (!get<alive>(m_particles_begin)[delete_index]);

            // might not need to update search
            if (m_update_search) {
                if (m_serial) {
                    copy_points_per_particle(index,delete_index);
                } else {
                    copy_points_per_bucket(index,delete_index);
                }
            }
                
            // copy particle info
            i = *(m_particles_begin+delete_index);
        }

        
    }
};
*/


// assume that query functions, are only called from device code
// TODO: most of this code shared with bucket_search_parallel_query, need to combine them
template <typename Traits>
struct bucket_search_serial_query {

    typedef Traits traits_type;
    typedef typename Traits::raw_pointer raw_pointer;
    typedef typename Traits::double_d double_d;
    typedef typename Traits::bool_d bool_d;
    typedef typename Traits::int_d int_d;
    typedef typename Traits::unsigned_int_d unsigned_int_d;
    typedef typename Traits::reference particle_reference;
    typedef typename Traits::const_reference particle_const_reference;
    const static unsigned int dimension = Traits::dimension;
    template <int LNormNumber>
    using query_iterator = lattice_iterator_within_distance<bucket_search_serial_query,LNormNumber>;

    typedef lattice_iterator<dimension> all_iterator;
    typedef lattice_iterator<dimension> child_iterator;
    typedef typename query_iterator<2>::reference reference;
    typedef typename query_iterator<2>::pointer pointer;
    typedef typename query_iterator<2>::value_type value_type;
    typedef linked_list_iterator<Traits> particle_iterator;
    typedef detail::bbox<dimension> box_type;

    bool_d m_periodic;
    double_d m_bucket_side_length; 
    int_d m_end_bucket;
    detail::bbox<dimension> m_bounds;
    detail::point_to_bucket_index<dimension> m_point_to_bucket_index;

    raw_pointer m_particles_begin;
    raw_pointer m_particles_end;
    int *m_buckets_begin;
    int *m_linked_list_begin;
    int *m_id_map_key;
    int *m_id_map_value;

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    bucket_search_serial_query():
        m_periodic(),
        m_particles_begin(),
        m_buckets_begin()
    {
    #if defined(__CUDA_ARCH__)
        CHECK_CUDA((!std::is_same<typename Traits::template vector<double>,
                                  std::vector<double>>::value),
                   "Cannot use std::vector in device code");
    #endif
    }

    /*
     * functions for id mapping
     */
    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    raw_pointer find(const size_t id) const {
        const size_t n = number_of_particles();
        int *last = m_id_map_key+n;
        int *first = detail::lower_bound(m_id_map_key,last,id);
        if ((first != last) && !(id < *first)) {
            return m_particles_begin + m_id_map_value[first-m_id_map_key];
        } else {
            return m_particles_begin + n;
        }
    }
    

    /*
     * functions for updating search ds
     */
    /*
    CUDA_HOST_DEVICE
    void copy_points(raw_pointer copy_from, raw_pointer copy_to) {
        LOG(4,"neighbour_search_base: copy_points: fromi = "<<copy_from-m_particles_begin<<" toi = "<<copy_to-m_particles_begin);
        ASSERT((copy_to-m_particles_begin>=0) && 
                (m_particles_end-copy_to>0),"invalid copy to pointer");
        ASSERT((copy_from-m_particles_begin>=0) && 
                (m_particles_end-copy_from>0),"invalid copy from pointer");
        
        if (m_domain_has_been_set) {
            cast().copy_points_impl(copy_from_iterator,copy_to_iterator);
        }

        if (m_id_map) {
            copy_id_map(*get<id>(copy_from_iterator), 
                        *get<id>(copy_to_iterator), 
                        copy_to_iterator-m_particles_begin);
        }
    }
    */

    /*
     * functions for trees
     */
    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    static bool is_leaf_node(const value_type& bucket) {
        return true;
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    static bool is_tree() {
        return false;
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    child_iterator get_children() const {
        return child_iterator(int_d(0),m_end_bucket+1);
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    child_iterator get_children(const child_iterator& ci) const {
        return child_iterator();
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    const box_type get_bounds(const child_iterator& ci) const {
        box_type bounds;
        bounds.bmin = (*ci)*m_bucket_side_length + m_bounds.bmin;
        bounds.bmax = ((*ci)+1)*m_bucket_side_length + m_bounds.bmin;
        return bounds;
    }
    
    // dodgy hack cause nullptr cannot be converted to pointer
    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    static const pointer get_child1(const pointer& bucket) {
        CHECK(false,"this should not be called")
	    return pointer(-1);
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    static const pointer get_child2(const pointer& bucket) {
        CHECK(false,"this should not be called")
	    return pointer(-1);
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    const box_type& get_bounds() const { return m_bounds; }
    
    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    const bool_d& get_periodic() const { return m_periodic; }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    iterator_range<particle_iterator> 
    CUDA_HOST_DEVICE
    get_bucket_particles(const reference bucket) const {
#ifndef __CUDA_ARCH__
        ASSERT((bucket>=int_d(0)).all() && (bucket <= m_end_bucket).all(), "invalid bucket");
#endif
        
        const unsigned int bucket_index = m_point_to_bucket_index.collapse_index_vector(bucket);

#ifndef __CUDA_ARCH__
        LOG(4,"\tget_bucket_particles: looking in bucket "<<bucket<<" = "<<bucket_index);
#endif
        return iterator_range<particle_iterator>(
                particle_iterator(m_buckets_begin[bucket_index],
                    m_particles_begin,
                    m_linked_list_begin),
                particle_iterator());
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    detail::bbox<dimension> get_bucket_bbox(const reference bucket) const {
        return detail::bbox<dimension>(
                bucket*m_bucket_side_length + m_bounds.bmin,
                (bucket+1)*m_bucket_side_length + m_bounds.bmin
                );
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    box_type get_root_bucket_bounds(reference bucket) const {
        box_type bounds;
        bounds.bmin = bucket*m_bucket_side_length + m_bounds.bmin;
        bounds.bmax = (bucket+1)*m_bucket_side_length + m_bounds.bmin;
        return bounds;
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    void get_bucket(const double_d &position, pointer& bucket, box_type& bounds) const {
        bucket = m_point_to_bucket_index.find_bucket_index_vector(position);
        bounds.bmin = bucket*m_bucket_side_length + m_bounds.bmin;
        bounds.bmax = (bucket+1)*m_bucket_side_length + m_bounds.bmin;
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    size_t get_bucket_index(const reference bucket) const {
        return m_point_to_bucket_index.collapse_index_vector(bucket);
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    template <int LNormNumber=-1>
    CUDA_HOST_DEVICE
    iterator_range<query_iterator<LNormNumber>> 
    get_buckets_near_point(const double_d &position, const double max_distance) const {
#ifndef __CUDA_ARCH__
        LOG(4,"\tget_buckets_near_point: position = "<<position<<" max_distance = "<<max_distance);
#endif

        return iterator_range<query_iterator<LNormNumber>>(
                            query_iterator<LNormNumber>(position,double_d(max_distance),this),
                            query_iterator<LNormNumber>());
    }

    
     

    ABORIA_HOST_DEVICE_IGNORE_WARN
    template <int LNormNumber=-1>
    CUDA_HOST_DEVICE
    iterator_range<query_iterator<LNormNumber>> 
    get_buckets_near_point(const double_d &position, const double_d &max_distance) const {
#ifndef __CUDA_ARCH__
        LOG(4,"\tget_buckets_near_point: position = "<<position<<" max_distance = "<<max_distance);
#endif
        return iterator_range<query_iterator<LNormNumber>>(
                            query_iterator<LNormNumber>(position,max_distance,this),
                            query_iterator<LNormNumber>());
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    const int_d& get_end_bucket() const {
        return m_end_bucket;
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    iterator_range<all_iterator> get_subtree(const child_iterator& ci) const {
        return iterator_range<all_iterator>(
                all_iterator(),
                all_iterator());
    }
    
    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    iterator_range<all_iterator> get_subtree() const {
        return iterator_range<all_iterator>(
                all_iterator(int_d(0),m_end_bucket+1),
                all_iterator()
                );
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    size_t number_of_buckets() const {
        return (m_end_bucket+1).prod();
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    size_t number_of_particles() const {
        return (m_particles_end-m_particles_begin);
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    const raw_pointer& get_particles_begin() const {
        return m_particles_begin;
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    raw_pointer& get_particles_begin() {
        return m_particles_begin;
    }

    /*
    CUDA_HOST_DEVICE
    iterator_range<theta_iterator> get_theta_buckets(const reference bucket) const {
        
        int_d start = bucket-int_d(2);
        int_d end = bucket+int_d(2);

        bool no_buckets = false;
        for (int i=0; i<Traits::dimension; i++) {
            if (start[i] < 0) {
                start[i] = 0;
            } else if (start[i] > m_end_bucket[i]) {
                no_buckets = true;
                start[i] = m_end_bucket[i];
            }
            if (end[i] < 0) {
                no_buckets = true;
                end[i] = 0;
            } else if (end[i] > m_end_bucket[i]) {
                end[i] = m_end_bucket[i];
            }
        }
        if (no_buckets) {
            return iterator_range<theta_iterator>(
                theta_iterator(end,end,endend),
                ++theta_iterator(end,end,end)
                );

        } else {
            return iterator_range<theta_iterator>(
                theta_iterator(start,end,start),
                ++theta_iterator(start,end,end)
                );
        }
    }
    */

    

};

}

#endif /* BUCKETSEARCH_H_ */
