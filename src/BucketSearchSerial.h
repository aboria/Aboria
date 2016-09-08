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
                                 bucket_search_serial_params<Traits>,
                                 linked_list_iterator<Traits>,
                                 bucket_search_serial_query<Traits>> {

    typedef typename Traits::double_d double_d;
    typedef typename Traits::position position;
    typedef typename Traits::vector_double_d_const_iterator vector_double_d_const_iterator;
    typedef typename Traits::vector_unsigned_int_iterator vector_unsigned_int_iterator;
    typedef typename Traits::vector_unsigned_int vector_unsigned_int;
    typedef typename Traits::unsigned_int_d unsigned_int_d;
    typedef bucket_search_serial_params<Traits> params_type;

    friend neighbour_search_base<bucket_search_serial<Traits>,
                                 Traits,
                                 bucket_search_serial_params<Traits>,
                                 linked_list_iterator<Traits>,
                                 bucket_search_serial_query<Traits>>;


    void set_domain_impl(const params_type params) {
        m_bucket_side_length = params.side_length; 
        m_size = 
            floor((this->m_bounds.bmax-this->m_bounds.bmin)/m_bucket_side_length)
            .template cast<unsigned int>();
        m_bucket_side_length = (this->m_bounds.bmax-this->m_bounds.bmin)/m_size;
        m_point_to_bucket_index = 
            detail::point_to_bucket_index<Traits::dimension>(m_size,m_bucket_side_length,this->m_bounds);
 
	    LOG(2,"\tbucket_side_length = "<<m_bucket_side_length);
	    LOG(2,"\tnumber of buckets = "<<m_size<<" (total="<<m_size.prod()<<")");

        m_buckets.assign(m_size.prod(), CELL_EMPTY);
        use_dirty_cells = false;

        this->m_query.m_buckets = iterator_to_raw_pointer(this->m_buckets);
        this->m_query.m_nbuckets = m_bucket_begin.size();

        this->m_query.m_bounds.bmin = this->m_bounds.bmin;
        this->m_query.m_bounds.bmax = this->m_bounds.bmax;
        this->m_query.m_periodic = this->m_periodic;
        this->m_query.m_size = m_size;
        this->m_query.m_bucket_side_length = m_bucket_side_length;
        this->m_query.m_point_to_bucket_index = m_point_to_bucket_index;
    }

    void embed_points_impl() {
        const size_t n = this->m_particles_end - this->m_particles_begin;

        /*
         * clear head of linked lists (m_buckets)
         */
        if (use_dirty_cells) {
            if (m_dirty_buckets.size()<m_buckets.size()) {
                for (int i: m_dirty_buckets) {
                    m_buckets[i] = CELL_EMPTY;
                    for (int j: ghosting_indices_pb[i]) {
                        m_buckets[j] = CELL_EMPTY;
                    }
                }
            } else {
                m_buckets.assign(m_buckets.size(), CELL_EMPTY);
            }
        }
        use_dirty_cells = true;

        m_linked_list.assign(n, CELL_EMPTY);
        m_linked_list_reverse.assign(n, CELL_EMPTY);
        m_m_dirty_buckets.assign(n,CELL_EMPTY);
        int i = 0;
        for (size_t i; i<n; ++i) {
            const int celli = find_cell_index(get<position>(this->m_particles_begin[i]));
            const int cell_entry = m_buckets[celli];

            // Insert into own cell
            m_buckets[celli] = i;
            m_dirty_buckets[i] = celli;
            m_linked_list[i] = cell_entry;
            m_linked_list_reverse[i] = CELL_EMPTY;
            if (cell_entry != CELL_EMPTY) m_linked_list_reverse[cell_entry] = i;
        }

        this->m_query.m_particles_begin = iterator_to_raw_pointer(this->m_particles_begin);
        this->m_query.m_linked_list = iterator_to_raw_pointer(this->m_linked_list);
    }


    void add_points_at_end_impl(const size_t dist) {
        const size_t n = this->m_particles_end - this->m_particles_begin;
        const size_t start_adding = n-dist;
        ASSERT(m_linked_list.size() = start_adding);
        ASSERT(m_linked_list_reverse.size() = start_adding);
        ASSERT(m_dirty_buckets.size() = start_adding);
        m_linked_list.resize(n,CELL_EMPTY);
        m_linked_list_reverse.resize(n,CELL_EMPTY);
        m_dirty_buckets.resize(n,CELL_EMPTY);

        for (size_t i = start_adding; i<n; ++i) {
            const int celli = find_cell_index(get<position>(this->m_particles_begin[i]));
            const int cell_entry = m_buckets[celli];

            // Insert into own cell
            m_buckets[celli] = i;
            m_dirty_buckets[i] = celli;
            m_linked_list[i] = cell_entry;
            m_linked_list_reverse[i] = CELL_EMPTY;
            if (cell_entry != CELL_EMPTY) m_linked_list_reverse[cell_entry] = i;
        }

        this->m_query.m_particles_begin = iterator_to_raw_pointer(this->m_particles_begin);
        this->m_query.m_linked_list = iterator_to_raw_pointer(this->m_linked_list);
    }

    void delete_points_at_end_impl(const size_t dist) {
        const size_t n = this->m_particles_end - this->m_particles_begin;
        const size_t start_delete = n-dist;
        ASSERT(m_linked_list.size()-n = dist);
        ASSERT(m_linked_list_reverse.size()-n = dist);
        ASSERT(m_dirty_buckets.size()-n = dist);
        const size_t oldn = m_linked_list.size();
        for (size_t i = n; i<oldn; ++i) {
            untrack_point(i); 
        }
        m_linked_list.resize(n);
        m_linked_list_reverse.resize(n);
        m_dirty_buckets.resize(n);

        this->m_query.m_particles_begin = iterator_to_raw_pointer(this->m_particles_begin);
        this->m_query.m_linked_list = iterator_to_raw_pointer(this->m_linked_list);
    }

    void update_point(const_iterator update_iterator) {
        const unsigned int i = std::distance(begin_iterator,update_iterator);
        const bool particle_based = true;

        const int forwardi = m_linked_list[i];
        const int backwardsi = m_linked_list_reverse[i];

        if (forwardi != CELL_EMPTY) m_linked_list_reverse[forwardi] = backwardsi;
        if (backwardsi != CELL_EMPTY) {
            m_linked_list[backwardsi] = forwardi;
        } else {
            const int celli = m_dirty_buckets[i];
            ASSERT(m_buckets[celli]==i,"inconsistant m_buckets data structures!");
            m_buckets[celli] = forwardi;
        }

        const int celli = find_cell_index(return_vect3d(*update_iterator));
        const int cell_entry = m_buckets[celli];

        // Insert into own cell
        m_buckets[celli] = i;
        m_dirty_buckets[i] = celli;
        m_linked_list[i] = cell_entry;
        m_linked_list_reverse[i] = CELL_EMPTY;
        if (cell_entry != CELL_EMPTY) m_linked_list_reverse[cell_entry] = i;

    }

    void untrack_point(const size_t i) {
        ASSERT((i>=0) && (i<m_linked_list.size()),"invalid untrack index");

        const int forwardi = m_linked_list[i];
        const int backwardsi = m_linked_list_reverse[i];

        if (forwardi != CELL_EMPTY) m_linked_list_reverse[forwardi] = backwardsi;
        if (backwardsi != CELL_EMPTY) {
            m_linked_list[backwardsi] = forwardi;
        } else {
            const int celli = m_dirty_buckets[i];
            ASSERT(m_buckets[celli]==i,"inconsistant m_buckets data structures!");
            m_buckets[celli] = forwardi;
        }
    }


    void copy_points_impl(const_iterator copy_from_iterator, const_iterator copy_to_iterator) {
        const size_t toi = std::distance(this->m_particles_begin,copy_to_iterator);
        const size_t fromi = std::distance(this->m_particles_begin,copy_from_iterator);
                const int forwardi = m_linked_list[fromi];

        m_linked_list[toi] = forwardi;
        m_linked_list_reverse[toi] = fromi;
        m_linked_list[fromi] = toi;
        m_linked_list_reverse[forwardi] = toi;
        m_dirty_buckets[toi] = m_dirty_buckets[fromi];
    }

    const bucket_search_parallel_query<Traits>& get_query_impl() const {
        return m_query;
    }

    /*
    const bucket_search_parallel_query<Traits>& get_query() const {
        return m_query;
    }
    */


template<typename T, typename F>
void BucketSearch<T,F>::embed_points(const T _begin_iterator, const T _end_iterator) {
	
}

template<typename T, typename F>
void BucketSearch<T,F>::add_point(const T point_to_add_iterator) {
	
}

template<typename T, typename F>
void BucketSearch<T,F>::

template<typename T, typename F>
void BucketSearch<T,F>::reset(const double_d& _low, const double_d& _high, double _max_interaction_radius, const Vect3b& _periodic) {
	
}
template<typename T, typename F>
typename BucketSearch<T,F>::const_iterator BucketSearch<T,F>::find_broadphase_neighbours(const double_d& r,const int my_index, const bool self) const {
	return const_iterator(this,correct_position_for_periodicity(r),my_index,self);
}
template<typename T, typename F>
typename BucketSearch<T,F>::const_iterator BucketSearch<T,F>::end() const {
	return const_iterator();
}
template<typename T, typename F>
double_d BucketSearch<T,F>::correct_position_for_periodicity(const double_d& source_r, const double_d& to_correct_r) const {
	double_d corrected_r = to_correct_r - source_r;
	for (int i = 0; i < NDIM; ++i) {
		if (!periodic[i]) continue;
		if (corrected_r[i] > domain_size[i]/2.0) corrected_r[i] -= domain_size[i];
		else if (corrected_r[i] < -domain_size[i]/2.0) corrected_r[i] += domain_size[i];
	}
	return corrected_r + source_r;
}

template<typename T, typename F>
double_d BucketSearch<T,F>::correct_position_for_periodicity(const double_d& to_correct_r) const {
	double_d corrected_r = to_correct_r;
	for (int i = 0; i < NDIM; ++i) {
		if (!periodic[i]) continue;
		while (corrected_r[i] >= high[i]) corrected_r[i] -= domain_size[i];
		while (corrected_r[i] < low[i]) corrected_r[i] += domain_size[i];
	}
	return corrected_r;
}



}

#endif /* BUCKETSEARCH_H_ */
