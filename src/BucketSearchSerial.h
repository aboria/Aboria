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
    bucket_search_serial():m_size_calculated_with_n(-1),base_type() {}
    static constexpr bool cheap_copy_and_delete_at_end() {
        return true;
    }

private:
    void set_domain_impl() {
        const size_t n = this->m_particles_end - this->m_particles_begin;
        if (n < 0.5*m_size_calculated_with_n || n > 2*m_size_calculated_with_n) {
            m_size_calculated_with_n = n;
            LOG(2,"bucket_search_serial: recalculating bucket size");
            if (this->m_n_particles_in_leaf > n) {
                m_size = unsigned_int_d(1);
            } else {
                const double total_volume = (this->m_bounds.bmax-this->m_bounds.bmin).prod();
                const double box_volume = double(this->m_n_particles_in_leaf)/double(n)*total_volume;
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
            //TODO: should always be true?
            m_use_dirty_cells = true;

            this->m_query.m_buckets_begin = iterator_to_raw_pointer(m_buckets.begin());

            this->m_query.m_bucket_side_length = this->m_bucket_side_length;
            this->m_query.m_bounds.bmin = this->m_bounds.bmin;
            this->m_query.m_bounds.bmax = this->m_bounds.bmax;
            this->m_query.m_periodic = this->m_periodic;
            this->m_query.m_end_bucket = m_size-1;
            this->m_query.m_point_to_bucket_index = m_point_to_bucket_index;
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



    void update_iterator_impl() {
        this->m_query.m_particles_begin = iterator_to_raw_pointer(this->m_particles_begin);
        //check_data_structure();
    }


    void embed_points_impl() {
        set_domain_impl(); 
        const size_t n = this->m_particles_end - this->m_particles_begin;

        /*
         * clear head of linked lists (m_buckets)
         */
        if (m_use_dirty_cells) {
            //TODO: wont ever be true??
            if (m_dirty_buckets.size()<m_buckets.size()) {
                for (int i: m_dirty_buckets) {
                    m_buckets[i] = detail::get_empty_id();
                }
            } else {
                m_buckets.assign(m_buckets.size(), detail::get_empty_id());
            }
        }
        m_use_dirty_cells = true;

        m_linked_list.assign(n, detail::get_empty_id());
        m_linked_list_reverse.assign(n, detail::get_empty_id());
        m_dirty_buckets.assign(n,detail::get_empty_id());
        for (size_t i=0; i<n; ++i) {
            const double_d& r = get<position>(this->m_particles_begin)[i];
            const unsigned int bucketi = m_point_to_bucket_index.find_bucket_index(r);
            ASSERT(bucketi < m_buckets.size() && bucketi >= 0, "bucket index out of range");
            const int bucket_entry = m_buckets[bucketi];

            // Insert into own bucket
            m_buckets[bucketi] = i;
            m_dirty_buckets[i] = bucketi;
            m_linked_list[i] = bucket_entry;
            m_linked_list_reverse[i] = detail::get_empty_id();
            if (bucket_entry != detail::get_empty_id()) m_linked_list_reverse[bucket_entry] = i;
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
                LOG(4,"\ti = "<<i<<" p = "<<get<position>(*(this->m_particles_begin+i))<<" contents = "<<m_linked_list[i]<<". reverse = "<<m_linked_list_reverse[i]);
            }
            LOG(4,"\tend linked list:");
        }
#endif

        //check_data_structure();

        this->m_query.m_particles_begin = iterator_to_raw_pointer(this->m_particles_begin);
        this->m_query.m_linked_list_begin = iterator_to_raw_pointer(this->m_linked_list.begin());
    }


    void add_points_at_end_impl(const size_t dist) {
        set_domain_impl(); 
        const size_t n = this->m_particles_end - this->m_particles_begin;
        const size_t start_adding = n-dist;
        ASSERT(m_linked_list.size() == start_adding, "m_linked_list not consistent with dist");
        ASSERT(m_linked_list_reverse.size() == start_adding, "m_linked_list_reverse not consistent with dist");
        ASSERT(m_dirty_buckets.size() == start_adding, "m_dirty_buckets not consistent with dist");
        m_linked_list.resize(n,detail::get_empty_id());
        m_linked_list_reverse.resize(n,detail::get_empty_id());
        m_dirty_buckets.resize(n,detail::get_empty_id());

        for (size_t i = start_adding; i<n; ++i) {
            const double_d& r = get<position>(this->m_particles_begin)[i];
            const unsigned int bucketi = m_point_to_bucket_index.find_bucket_index(r);
            ASSERT(bucketi < m_buckets.size() && bucketi >= 0, "bucket index out of range");
            const int bucket_entry = m_buckets[bucketi];

            // Insert into own cell
            m_buckets[bucketi] = i;
            m_dirty_buckets[i] = bucketi;
            m_linked_list[i] = bucket_entry;
            m_linked_list_reverse[i] = detail::get_empty_id();
            if (bucket_entry != detail::get_empty_id()) m_linked_list_reverse[bucket_entry] = i;
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
            for (int i = 0; i<m_linked_list.size(); ++i) {
                LOG(4,"\ti = "<<i<<" p = "<<get<position>(*(this->m_particles_begin+i))<<" contents = "<<m_linked_list[i]<<". reverse = "<<m_linked_list_reverse[i]);
            }
        }
#endif

        //check_data_structure();

        this->m_query.m_particles_begin = iterator_to_raw_pointer(this->m_particles_begin);
        this->m_query.m_linked_list_begin = iterator_to_raw_pointer(this->m_linked_list.begin());
    }

    void delete_points_at_end_impl(const size_t dist) {
        set_domain_impl(); 
        const size_t n = this->m_particles_end - this->m_particles_begin;
        ASSERT(m_linked_list.size()-n == dist, "m_linked_list not consistent with dist");
        ASSERT(m_linked_list_reverse.size()-n == dist, "m_linked_list_reverse not consistent with dist");
        ASSERT(m_dirty_buckets.size()-n == dist, "m_dirty_buckets not consistent with dist");
        const size_t oldn = m_linked_list.size();
        for (size_t i = n; i<oldn; ++i) {
            if (m_dirty_buckets[i] == detail::get_empty_id()) continue;
            const int celli = m_dirty_buckets[i];

            //get first backwards index < n
            int backwardsi = i;
            while (backwardsi >= n) {
                backwardsi = m_linked_list_reverse[backwardsi];
                if (backwardsi == detail::get_empty_id()) break;
            }

            //get first forward index < n
            int forwardi = i;
            while (forwardi >= n) {
                forwardi = m_linked_list[forwardi];
                if (forwardi == detail::get_empty_id()) break;
            }

            if (forwardi != detail::get_empty_id()) {
                m_linked_list_reverse[forwardi] = backwardsi;
            }
            if (backwardsi != detail::get_empty_id()) {
                m_linked_list[backwardsi] = forwardi;
            } else if (forwardi != detail::get_empty_id()) {
                m_buckets[celli] = forwardi;
            } else {
                m_buckets[celli] = detail::get_empty_id();
            }
        }
        m_linked_list.resize(n);
        m_linked_list_reverse.resize(n);
        m_dirty_buckets.resize(n);


        //check_data_structure();

        this->m_query.m_particles_begin = iterator_to_raw_pointer(this->m_particles_begin);
        this->m_query.m_linked_list_begin = iterator_to_raw_pointer(this->m_linked_list.begin());
    }

    void update_point(iterator update_iterator) {
        const size_t i = std::distance(this->m_particles_begin,update_iterator);
        const bool particle_based = true;

        const int forwardi = m_linked_list[i];
        const int backwardsi = m_linked_list_reverse[i];

        if (forwardi != detail::get_empty_id()) {
            m_linked_list_reverse[forwardi] = backwardsi;
        }
        if (backwardsi != detail::get_empty_id()) {
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
        m_linked_list_reverse[i] = detail::get_empty_id();
        if (cell_entry != detail::get_empty_id()) m_linked_list_reverse[cell_entry] = i;

    }

    void untrack_point(const size_t i) {
        ASSERT((i>=0) && (i<m_linked_list.size()),"invalid untrack index");

        const int forwardi = m_linked_list[i];
        const int backwardsi = m_linked_list_reverse[i];

        ASSERT(forwardi == detail::get_empty_id() ||
                ((forwardi>=0) && (forwardi<m_linked_list_reverse.size())),"invalid untrack index (forwards)");
        ASSERT(backwardsi == detail::get_empty_id() ||
                ((backwardsi>=0) && (backwardsi<m_linked_list.size())),"invalid untrack index (backwards)");

        if (forwardi != detail::get_empty_id()) {
            m_linked_list_reverse[forwardi] = backwardsi;
        }
        if (backwardsi != detail::get_empty_id()) {
            m_linked_list[backwardsi] = forwardi;
        } else {
            const int celli = m_dirty_buckets[i];
            ASSERT(m_buckets[celli]==i,"inconsistant m_buckets data structures!");
            m_buckets[celli] = forwardi;
        }
    }


    void copy_points_impl(iterator copy_from_iterator, iterator copy_to_iterator) {
        const size_t toi = std::distance(this->m_particles_begin,copy_to_iterator);
        const size_t fromi = std::distance(this->m_particles_begin,copy_from_iterator);
        ASSERT(toi != fromi,"toi and fromi are the same");

        // unlink old toi pointers
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

        // setup fromi <-> toi 
        const int forwardi = m_linked_list[fromi];
        const int bucketi = m_dirty_buckets[fromi];
        ASSERT(toi < m_linked_list.size() && toi >= 0,"invalid index");
        ASSERT(bucketi < m_buckets.size() && bucketi >= 0,"invalid index");
        ASSERT(fromi < m_linked_list_reverse.size() && fromi >= 0,"invalid index");
        ASSERT(forwardi == detail::get_empty_id() ||
            (forwardi < m_linked_list_reverse.size() && forwardi>= 0),
            "invalid index of "<<forwardi<<". Note linked list reverse size is "
            << m_linked_list_reverse.size());

        m_linked_list[toi] = forwardi;
        m_linked_list_reverse[toi] = fromi;
        m_linked_list[fromi] = toi;
        if (forwardi != detail::get_empty_id()) { //check this
            m_linked_list_reverse[forwardi] = toi; 
        }
        m_dirty_buckets[toi] = bucketi;

        //check_data_structure();
    }

    const bucket_search_serial_query<Traits>& get_query_impl() const {
        return m_query;
    }

    vector_int m_buckets;
    vector_int m_linked_list;
    vector_int m_linked_list_reverse;
    vector_int m_dirty_buckets;
    bucket_search_serial_query<Traits> m_query;
    bool m_use_dirty_cells;

    size_t m_size_calculated_with_n;
    unsigned_int_d m_size;
    double_d m_bucket_side_length;
    detail::point_to_bucket_index<Traits::dimension> m_point_to_bucket_index;

};

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
    const static unsigned int dimension = Traits::dimension;
    typedef lattice_iterator<dimension> query_iterator;
    typedef lattice_iterator<dimension> all_iterator;
    typedef lattice_iterator<dimension> child_iterator;
    typedef typename query_iterator::reference reference;
    typedef typename query_iterator::pointer pointer;
    typedef typename query_iterator::value_type value_type;
    typedef linked_list_iterator<Traits> particle_iterator;
    typedef detail::bbox<dimension> box_type;

    bool_d m_periodic;
    double_d m_bucket_side_length; 
    int_d m_end_bucket;
    detail::bbox<dimension> m_bounds;
    detail::point_to_bucket_index<dimension> m_point_to_bucket_index;

    raw_pointer m_particles_begin;
    int *m_buckets_begin;
    int *m_linked_list_begin;

    inline
    CUDA_HOST_DEVICE
    bucket_search_serial_query():
        m_periodic(),
        m_particles_begin(),
        m_buckets_begin()
    {}

    /*
     * functions for trees
     */
    static bool is_leaf_node(const value_type& bucket) {
        return true;
    }

    static bool is_tree() {
        return false;
    }

    child_iterator get_children() const {
        return child_iterator(int_d(0),m_end_bucket+1);
    }

    child_iterator get_children(const child_iterator& ci) const {
        return child_iterator();
    }

    const box_type get_bounds(const child_iterator& ci) const {
        box_type bounds;
        bounds.bmin = (*ci)*m_bucket_side_length + m_bounds.bmin;
        bounds.bmax = ((*ci)+1)*m_bucket_side_length + m_bounds.bmin;
        return bounds;
    }
    
    // dodgy hack cause nullptr cannot be converted to pointer
    static const pointer get_child1(const pointer& bucket) {
        CHECK(false,"this should not be called")
	    return pointer(-1);
    }
    static const pointer get_child2(const pointer& bucket) {
        CHECK(false,"this should not be called")
	    return pointer(-1);
    }

    const box_type& get_bounds() const { return m_bounds; }
    const bool_d& get_periodic() const { return m_periodic; }

    CUDA_HOST_DEVICE
    iterator_range<particle_iterator> 
    get_bucket_particles(const reference bucket) const {
        ASSERT((bucket>=int_d(0)).all() && (bucket <= m_end_bucket).all(), "invalid bucket");
        
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

    CUDA_HOST_DEVICE
    detail::bbox<dimension> get_bucket_bbox(const reference bucket) const {
        return detail::bbox<dimension>(
                bucket*m_bucket_side_length + m_bounds.bmin,
                (bucket+1)*m_bucket_side_length + m_bounds.bmin
                );
    }

    CUDA_HOST_DEVICE
    box_type get_root_bucket_bounds(reference bucket) const {
        box_type bounds;
        bounds.bmin = bucket*m_bucket_side_length + m_bounds.bmin;
        bounds.bmax = (bucket+1)*m_bucket_side_length + m_bounds.bmin;
        return bounds;
    }

    CUDA_HOST_DEVICE
    void get_bucket(const double_d &position, pointer& bucket, box_type& bounds) const {
        bucket = m_point_to_bucket_index.find_bucket_index_vector(position);
        bounds.bmin = bucket*m_bucket_side_length + m_bounds.bmin;
        bounds.bmax = (bucket+1)*m_bucket_side_length + m_bounds.bmin;
    }

    CUDA_HOST_DEVICE
    size_t get_bucket_index(const reference bucket) const {
        return m_point_to_bucket_index.collapse_index_vector(bucket);
    }

    template <int LNormNumber=-1>
    CUDA_HOST_DEVICE
    iterator_range<query_iterator> 
    get_buckets_near_point(const double_d &position, const double max_distance) const {
        return get_buckets_near_point(position,double_d(max_distance));
    }
     

    template <int LNormNumber=-1>
    CUDA_HOST_DEVICE
    iterator_range<query_iterator> 
    get_buckets_near_point(const double_d &position, const double_d &max_distance) const {
#ifndef __CUDA_ARCH__
        LOG(4,"\tget_buckets_near_point: position = "<<position<<" max_distance = "<<max_distance);
#endif
 
        value_type bucket = m_point_to_bucket_index.find_bucket_index_vector(position);
        int_d start = m_point_to_bucket_index.find_bucket_index_vector(position-max_distance);
        int_d end = m_point_to_bucket_index.find_bucket_index_vector(position+max_distance);

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
#ifndef __CUDA_ARCH__
        LOG(4,"\tget_buckets_near_point: looking in bucket "<<bucket<<". start = "<<start<<" end = "<<end<<" no_buckets = "<<no_buckets);
#endif
        if (no_buckets) {
            return iterator_range<query_iterator>(
                    query_iterator()
                    ,query_iterator()
                    );
        } else {
            return iterator_range<query_iterator>(
                    query_iterator(start,end+1)
                    ,query_iterator()
                    );
        }
    }


    iterator_range<all_iterator> get_subtree(const child_iterator& ci) const {
        return iterator_range<all_iterator>(
                all_iterator(),
                all_iterator());
    }
    
    iterator_range<all_iterator> get_subtree() const {
        return iterator_range<all_iterator>(
                all_iterator(int_d(0),m_end_bucket+1),
                all_iterator()
                );
    }

    size_t number_of_buckets() const {
        return (m_end_bucket+1).prod();
    }

    raw_pointer get_particles_begin() const {
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
