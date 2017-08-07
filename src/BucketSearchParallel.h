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
// Acknowledgement: This source was modified from the Thrust example bucket_sort2d.cu
//


#ifndef BUCKETSEARCH_H_
#define BUCKETSEARCH_H_

#include "detail/Algorithms.h"
#include "detail/SpatialUtil.h"
#include "NeighbourSearchBase.h"
#include "Traits.h"
#include "CudaInclude.h"
#include "Vector.h"
#include "Get.h"

#include <iostream>
#include "Log.h"

namespace Aboria {

   
template <typename Traits>
struct bucket_search_parallel_params {
    typedef typename Traits::double_d double_d;
    bucket_search_parallel_params(): 
        side_length(detail::get_max<double>()) {}
    bucket_search_parallel_params(const double_d& side_length):
        side_length(side_length) {}
    double_d side_length;
};

template <typename Traits>
class bucket_search_parallel_query; 

template <typename Traits>
class bucket_search_parallel: 
    public neighbour_search_base<bucket_search_parallel<Traits>,
                                 Traits,
                                 bucket_search_parallel_query<Traits>> {

    typedef typename Traits::double_d double_d;
    typedef typename Traits::position position;
    typedef typename Traits::vector_double_d_const_iterator vector_double_d_const_iterator;
    typedef typename Traits::vector_unsigned_int_iterator vector_unsigned_int_iterator;
    typedef typename Traits::vector_unsigned_int vector_unsigned_int;
    typedef typename Traits::unsigned_int_d unsigned_int_d;
    typedef typename Traits::iterator iterator;
    typedef bucket_search_parallel_params<Traits> params_type;

    typedef neighbour_search_base<bucket_search_parallel<Traits>,
                                 Traits,
                                 bucket_search_parallel_query<Traits>> base_type;

    friend base_type;


public:
    bucket_search_parallel():m_size_calculated_with_n(std::numeric_limits<size_t>::max()),base_type() {}

    static constexpr bool ordered() {
        return true;
    }

    struct delete_points_lambda;

private:


    bool set_domain_impl() {
        const size_t n = this->m_particles_end - this->m_particles_begin;
        if (n < 0.5*m_size_calculated_with_n || n > 2*m_size_calculated_with_n) {
            LOG(2,"bucket_search_serial: recalculating bucket size");
            m_size_calculated_with_n = n;
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

            // setup bucket data structures
            m_bucket_begin.resize(m_size.prod());
            m_bucket_end.resize(m_size.prod());

            this->m_query.m_bucket_begin = iterator_to_raw_pointer(m_bucket_begin.begin());
            this->m_query.m_bucket_end = iterator_to_raw_pointer(m_bucket_end.begin());
            this->m_query.m_nbuckets = m_bucket_begin.size();

            this->m_query.m_bucket_side_length = m_bucket_side_length;
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

    void update_iterator_impl() {
        this->m_query.m_particles_begin = iterator_to_raw_pointer(this->m_particles_begin);
        this->m_query.m_particles_end = iterator_to_raw_pointer(this->m_particles_end);
    }

    void embed_points_impl(const bool call_set_domain=true) {
        if (call_set_domain) {
            set_domain_impl();
        }
        const size_t n = this->m_particles_end - this->m_particles_begin;
        m_bucket_indices.resize(n);
        if (n > 0) {
            build_bucket_indices(
                    get<position>(this->m_particles_begin),
                    get<position>(this->m_particles_end),m_bucket_indices.begin());
            sort_by_bucket_index();
        }
        build_buckets();

        this->m_query.m_particles_begin = iterator_to_raw_pointer(this->m_particles_begin);
        this->m_query.m_particles_end = iterator_to_raw_pointer(this->m_particles_end);
    }

    iterator_range<vector_unsigned_int::const_iterator>
    get_order_impl(iterator begin, iterator end) {
        set_domain_impl();
        const size_t n = begin - end;
        m_bucket_indices.resize(n);
        if (n > 0) {
            build_bucket_indices(
                    get<position>(begin),
                    get<position>(end),
                    m_bucket_indices.begin());
            sort_by_bucket_index();
        }
        return iterator_range<vector_unsigned_int::const_iterator>(
                                                m_indices.begin(),
                                                m_indices.end());
    }


    void add_points_at_end_impl(const size_t dist) {
        const bool embed_all = set_domain_impl();
        auto start_adding = embed_all?this->m_particles_begin:
                                      (this->m_particles_end-dist);
        const size_t total = m_bucket_indices.size() + dist;
        auto positions_start_adding = embed_all?get<position>(this->m_particles_begin):
                                            get<position>(this->m_particles_end)- dist;
        m_bucket_indices.resize(total);
        auto bucket_indices_start_adding = embed_all?m_bucket_indices.begin():
                                                      m_bucket_indices.end() - dist;
        build_bucket_indices(positions_start_adding,
                             get<position>(this->m_particles_end),
                             bucket_indices_start_adding);
        sort_by_bucket_index();
        build_buckets();

#ifndef __CUDA_ARCH__
        if (4 <= ABORIA_LOG_LEVEL) { 
            for (int i = 0; i<m_bucket_indices.size(); ++i) {
                LOG(4,"\tp = "<<get<position>(*(this->m_particles_begin+i))<<" index = "<<m_bucket_indices[i]<<". m_bucket_begin[index] = "<<m_bucket_begin[m_bucket_indices[i]]<<". m_bucket_end[index] = "<<m_bucket_end[m_bucket_indices[i]]);
            }
        }
#endif

        this->m_query.m_particles_begin = iterator_to_raw_pointer(this->m_particles_begin);
        this->m_query.m_particles_end = iterator_to_raw_pointer(this->m_particles_end);
    }

    bool delete_points_impl(const size_t i, const size_t n) {
        const bool resize_buckets = set_domain_impl();
        if (resize_buckets) {
            // buckets all changed, so start from scratch
            embed_points_impl(false);
            return true;
        } else {
            // only redo buckets that changed
            const size_t start_bucket = m_bucket_indices[i];
            const size_t end_bucket = m_bucket_indices[i+n]+1;

            m_bucket_indices.erase(m_bucket_indices.begin()+i,
                                   m_bucket_indices.begin()+i+n);
            
            //set begins in deleted range to i
            detail::fill(m_buckets_begin.begin()+start_bucket,
                         m_buckets_begin.begin()+end_bucket,
                         i);
            
            //set ends in deleted range to i
            detail::fill(m_buckets_end.begin()+start_bucket,
                         m_buckets_end.begin()+end_bucket,
                         i);

            //minus n from begins after deleted range
            detail::transform(m_buckets_begin.begin()+end_bucket,
                              m_buckets_begin.end(),
                              m_buckets_begin.begin()+end_bucket,
                              detail::_1 - n);

            //minus n from ends after deleted range
            detail::transform(m_buckets_end.begin()+end_bucket,
                              m_buckets_end.end(),
                              m_buckets_end.begin()+end_bucket,
                              detail::_1 - n);

            return false;
        }
        

    }

    void copy_points_impl(iterator copy_from_iterator, iterator copy_to_iterator) {
        ERROR("data structure depends on particle ordering, cannot copy");
    }


    const bucket_search_parallel_query<Traits>& get_query_impl() const {
        return m_query;
    }

    /*
    const bucket_search_parallel_query<Traits>& get_query() const {
        return m_query;
    }
    */
    
    void build_buckets() {
        // find the beginning of each bucket's list of points
        detail::counting_iterator<unsigned int> search_begin(0);
        detail::lower_bound(m_bucket_indices.begin(),
                m_bucket_indices.end(),
                search_begin,
                search_begin + m_size.prod(),
                m_bucket_begin.begin());

        // find the end of each bucket's list of points
        detail::upper_bound(m_bucket_indices.begin(),
                m_bucket_indices.end(),
                search_begin,
                search_begin + m_size.prod(),
                m_bucket_end.begin());
    }

    void build_bucket_indices(
        vector_double_d_const_iterator positions_begin,
        vector_double_d_const_iterator positions_end,
        vector_unsigned_int_iterator bucket_indices_begin) {
        // transform the points to their bucket indices
        detail::transform(positions_begin,
                positions_end,
                bucket_indices_begin,
                m_point_to_bucket_index);
    }

    void sort_by_bucket_index() {
        // sort the points by their bucket index
        if (m_bucket_indices.size() > 0) {
            m_indices.resize(m_bucket_indices.size());
            detail::sequence(m_indices.begin(), m_indices.end());
            detail::sort_by_key(m_bucket_indices.begin(),
                                m_bucket_indices.end(),
                                m_indices.begin());
        }
    }

 

    // the grid data structure keeps a range per grid bucket:
    // each bucket_begin[i] indexes the first element of bucket i's list of points
    // each bucket_end[i] indexes one past the last element of bucket i's list of points
    vector_unsigned_int m_bucket_begin;
    vector_unsigned_int m_bucket_end;
    vector_unsigned_int m_bucket_indices;
    vector_unsigned_int m_indices;
    bucket_search_parallel_query<Traits> m_query;

    double_d m_bucket_side_length;
    unsigned_int_d m_size;
    size_t m_size_calculated_with_n;
    detail::point_to_bucket_index<Traits::dimension> m_point_to_bucket_index;
};


// assume that query functions, are only called from device code
template <typename Traits>
struct bucket_search_parallel_query {

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
    typedef ranges_iterator<Traits> particle_iterator;
    typedef detail::bbox<dimension> box_type;

    raw_pointer m_particles_begin;
    raw_pointer m_particles_end;

    bool_d m_periodic;
    double_d m_bucket_side_length; 
    int_d m_end_bucket;
    detail::bbox<dimension> m_bounds;
    detail::point_to_bucket_index<dimension> m_point_to_bucket_index;

    unsigned int *m_bucket_begin;
    unsigned int *m_bucket_end;
    unsigned int m_nbuckets;

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    bucket_search_parallel_query():
        m_periodic(),
        m_particles_begin(),
        m_bucket_begin()
    {}

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

    //const double_d& get_min_bucket_size() const { return m_bucket_side_length; }
    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    const box_type& get_bounds() const { return m_bounds; }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    const bool_d& get_periodic() const { return m_periodic; }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    iterator_range<particle_iterator> get_bucket_particles(const reference bucket) const {
#ifndef __CUDA_ARCH__
        ASSERT((bucket>=int_d(0)).all() && (bucket <= m_end_bucket).all(), "invalid bucket");
#endif

        const unsigned int bucket_index = m_point_to_bucket_index.collapse_index_vector(bucket);
        const unsigned int range_start_index = m_bucket_begin[bucket_index]; 
        const unsigned int range_end_index = m_bucket_end[bucket_index]; 

#ifndef __CUDA_ARCH__
        LOG(4,"\tlooking in bucket "<<bucket<<" = "<<bucket_index<<". found "<<range_end_index-range_start_index<<" particles");
#endif
        return iterator_range<particle_iterator>(
                particle_iterator(m_particles_begin + range_start_index),
                particle_iterator(m_particles_begin + range_end_index));
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
    iterator_range<query_iterator> 
    get_buckets_near_point(const double_d &position, const double max_distance) const {
        return get_buckets_near_point(position,double_d(max_distance));
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    template <int LNormNumber=-1>
    CUDA_HOST_DEVICE
    iterator_range<query_iterator> 
    get_buckets_near_point(const double_d &position, const double_d& max_distance) const {
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
        LOG(4,"\tget_buckets_near_point: looking in bucket "<<bucket<<". start = "<<start<<" end = "<<end);
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
    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    bool get_children_buckets(const bucket_reference &bucket, std::array<value_type,2>& children) {
        return false;
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    iterator_range<query_iterator> get_root_buckets() const {
        return iterator_range<query_iterator>(
                query_iterator(int_d(0),m_end_bucket,int_d(0)),
                ++query_iterator(int_d(0),m_end_bucket,m_end_bucket)
                );
    }
    */
};

   


}


#endif /* BUCKETSEARCH_H_ */
