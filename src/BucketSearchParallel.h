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
                                 bucket_search_parallel_params<Traits>,
                                 ranges_iterator<Traits>,
                                 bucket_search_parallel_query<Traits>> {

    typedef typename Traits::double_d double_d;
    typedef typename Traits::position position;
    typedef typename Traits::vector_double_d_const_iterator vector_double_d_const_iterator;
    typedef typename Traits::vector_unsigned_int_iterator vector_unsigned_int_iterator;
    typedef typename Traits::vector_unsigned_int vector_unsigned_int;
    typedef typename Traits::unsigned_int_d unsigned_int_d;
    typedef typename Traits::iterator iterator;
    typedef bucket_search_parallel_params<Traits> params_type;

    friend neighbour_search_base<bucket_search_parallel<Traits>,
                                 Traits,
                                 bucket_search_parallel_params<Traits>,
                                 ranges_iterator<Traits>,
                                 bucket_search_parallel_query<Traits>>;

public:
    static constexpr bool unordered() {
        return false;
    }

private:

    void set_domain_impl(const params_type params) {
        m_bucket_side_length = params.side_length; 
        m_size = 
            floor((this->m_bounds.bmax-this->m_bounds.bmin)/m_bucket_side_length)
            .template cast<unsigned int>();
        for (int i=0; i<Traits::dimension; ++i) {
            if (m_size[i] == 0) {
                m_size[i] = this->m_bounds.bmax[i]-this->m_bounds.bmin[i];
            }
        }
        m_bucket_side_length = (this->m_bounds.bmax-this->m_bounds.bmin)/m_size;
        m_point_to_bucket_index = 
            detail::point_to_bucket_index<Traits::dimension>(m_size,m_bucket_side_length,this->m_bounds);
 
	    LOG(2,"\tbucket_side_length = "<<m_bucket_side_length);
	    LOG(2,"\tnumber of buckets = "<<m_size<<" (total="<<m_size.prod()<<")");

        // setup bucket data structures
        m_bucket_begin.resize(m_size.prod());
        m_bucket_end.resize(m_size.prod());

        this->m_query.m_bucket_begin = iterator_to_raw_pointer(m_bucket_begin.begin());
        this->m_query.m_bucket_end = iterator_to_raw_pointer(m_bucket_end.begin());
        this->m_query.m_nbuckets = m_bucket_begin.size();

        this->m_query.m_bounds.bmin = this->m_bounds.bmin;
        this->m_query.m_bounds.bmax = this->m_bounds.bmax;
        this->m_query.m_periodic = this->m_periodic;
        this->m_query.m_size = m_size;
        this->m_query.m_bucket_side_length = m_bucket_side_length;
        this->m_query.m_point_to_bucket_index = m_point_to_bucket_index;
    }

    void update_iterator_impl() {
        this->m_query.m_particles_begin = iterator_to_raw_pointer(this->m_particles_begin);
        this->m_query.m_particles_end = iterator_to_raw_pointer(this->m_particles_end);
    }

    void embed_points_impl() {
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


    void add_points_at_end_impl(const size_t dist) {
        auto start_adding = this->m_particles_end-dist;
        const size_t total = m_bucket_indices.size() + dist;
        auto positions_end = get<position>(this->m_particles_end);
        auto positions_start_adding = positions_end - dist;
        m_bucket_indices.resize(total);
        auto bucket_indices_start_adding = m_bucket_indices.end() - dist;

        build_bucket_indices(positions_start_adding,positions_end,bucket_indices_start_adding);
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

    void delete_points_at_end_impl(const size_t dist) {
        const size_t n = this->m_particles_end - this->m_particles_begin;
        ASSERT(m_bucket_indices.size()-n == dist,"m_bucket_indices size not consistent with dist argument");
        const size_t oldn = m_bucket_indices.size();
        m_bucket_indices.resize(n);

        build_buckets();

        this->m_query.m_particles_begin = iterator_to_raw_pointer(this->m_particles_begin);
        this->m_query.m_particles_end = iterator_to_raw_pointer(this->m_particles_end);
    }

    void copy_points_impl(iterator copy_from_iterator, iterator copy_to_iterator) {
        auto positions_from = get<position>(copy_from_iterator);
        auto positions_to = get<position>(copy_to_iterator);

        const size_t toi = std::distance(this->m_particles_begin,copy_to_iterator);
        const size_t fromi = std::distance(this->m_particles_begin,copy_from_iterator);

        build_bucket_indices(positions_from,positions_from+1,
               m_bucket_indices.begin() + fromi);
        build_bucket_indices(positions_to,positions_to+1,
               m_bucket_indices.begin() + toi);
        sort_by_bucket_index();
        build_buckets();
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
            detail::sort_by_key(m_bucket_indices.begin(),
                m_bucket_indices.end(),
                this->m_particles_begin);
        }
    }

 

    // the grid data structure keeps a range per grid bucket:
    // each bucket_begin[i] indexes the first element of bucket i's list of points
    // each bucket_end[i] indexes one past the last element of bucket i's list of points
    vector_unsigned_int m_bucket_begin;
    vector_unsigned_int m_bucket_end;
    vector_unsigned_int m_bucket_indices;
    bucket_search_parallel_query<Traits> m_query;

    double_d m_bucket_side_length; 
    unsigned_int_d m_size;
    detail::point_to_bucket_index<Traits::dimension> m_point_to_bucket_index;
};


// assume that query functions, are only called from device code
template <typename Traits>
struct bucket_search_parallel_query {

    typedef typename Traits::raw_pointer raw_pointer;
    typedef typename Traits::double_d double_d;
    typedef typename Traits::bool_d bool_d;
    typedef typename Traits::int_d int_d;
    typedef typename Traits::unsigned_int_d unsigned_int_d;
    typedef typename Traits::reference reference;
    typedef typename Traits::position position;

    raw_pointer m_particles_begin;
    raw_pointer m_particles_end;

    bool_d m_periodic;
    double_d m_bucket_side_length; 
    unsigned_int_d m_size;
    detail::bbox<Traits::dimension> m_bounds;
    detail::point_to_bucket_index<Traits::dimension> m_point_to_bucket_index;


    unsigned int *m_bucket_begin;
    unsigned int *m_bucket_end;
    unsigned int m_nbuckets;

    inline
    CUDA_HOST_DEVICE
    bucket_search_parallel_query():
        m_periodic(),
        m_particles_begin(),
        m_bucket_begin()
    {}

    CUDA_HOST_DEVICE
    iterator_range<ranges_iterator<Traits>> get_neighbours(const double_d& position) const {
        return iterator_range<ranges_iterator<Traits>>(find_broadphase_neighbours(position,1,false),
                                              ranges_iterator<Traits>(m_particles_end));
    }


    CUDA_HOST_DEVICE
    iterator_range<linked_list_iterator<Traits>> get_bucket_particles(const bucket_iterator::reference bucket) const {
        linked_list_iterator<Traits> search_iterator(m_particles_begin,
                            m_bucket_side_length,
                            m_linked_list_begin,
                            m_buckets_begin,get_bucket_bbox(bucket).centre());
        const int bucket_index = m_point_to_bucket_index.collapse_index_vector(*bucket);
        search_iterator.add_bucket(bucket_index,transpose);
        return iterator_range<linked_list_iterator<Traits>>(
                search_iterator,
                linked_list_iterator<Traits>());
    }

    CUDA_HOST_DEVICE
    detail::bbox<dimension>& get_bucket_bbox(const bucket_iterator::reference bucket) const {
        return detail::bbox<dimension>((*bucket)*m_bucket_side_length + m_bounds.bmin);
    }

    CUDA_HOST_DEVICE
    iterator_range<bucket_iterator> get_all_buckets() const {
        return iterator_range<bucket_iterator>(
                bucket_iterator(0),
                bucket_iterator(m_size.prod()));
    }


    CUDA_HOST_DEVICE
    ranges_iterator<Traits> find_broadphase_neighbours(
            const double_d& r, 
            const double bucket_radius, 
            const bool self) const {
        
        ASSERT((r >= m_bounds.bmin).all() && (r < m_bounds.bmax).all(), "Error, search position "<<r<<" is outside neighbourhood search bounds " << m_bounds);
        const unsigned_int_d my_bucket = m_point_to_bucket_index.find_bucket_index_vector(r);

#ifndef __CUDA_ARCH__
        LOG(3,"BucketSearch: find_broadphase_neighbours: around r = "<<r<<". my_index = "<<my_index<<" self = "<<self);
        LOG(3,"\tbounds = "<<m_bounds);
	    LOG(3,"\tperiodic = "<<m_periodic);
	    LOG(3,"\tbucket_side_length = "<<m_bucket_side_length);
	    LOG(3,"\tnumber of buckets = "<<m_size<<" (total="<<m_size.prod()<<")");
#endif

        ranges_iterator<Traits> search_iterator(m_particles_end,m_bucket_side_length,r);
        lattice_iterator<Traits> bucket_iterator(my_bucket+int_d(-1),
                                                 my_bucket+ind_d(1),
                                                 my_bucket+int_d(-1));


        while ((*bucket_iterator)[0] <= my_bucket[0]+1) { 
            int_d& other_bucket = *bucket_iterator; 

            // handle end cases
            double_d transpose(0);
            bool outside = false;
            for (int i=0; i<Traits::dimension; i++) {
                if (other_bucket[i] < 0) {
                    if (m_periodic[i]) {
                        other_bucket[i] = m_size[i]-1;
                        transpose[i] = -(m_bounds.bmax-m_bounds.bmin)[i];
                    } else {
                        outside = true;
                        break;
                    }
                }
                if (other_bucket[i] == m_size[i]) {
                    if (m_periodic[i]) {
                        other_bucket[i] = 0;
                        transpose[i] = (m_bounds.bmax-m_bounds.bmin)[i];
                    } else {
                        outside = true;
                        break;
                    }
                }
            }

            if (!outside) {

                const unsigned int other_bucket_index = m_point_to_bucket_index.collapse_index_vector(other_bucket);
                const unsigned int range_start_index = m_bucket_begin[other_bucket_index]; 
                const unsigned int range_end_index = m_bucket_end[other_bucket_index]; 

#ifndef __CUDA_ARCH__
                LOG(4,"\tlooking in bucket "<<other_bucket<<" = "<<other_bucket_index<<". found "<<range_end_index-range_start_index<<" particles");
#endif

                if (range_end_index-range_start_index > 0) {

                    //std::cout << "adding range for my_bucket = "<<my_bucket<<" other_bucket = "<<other_bucket<<" range_start_index = "<<range_start_index<<" range_end_index = "<<range_end_index<< " transpose = "<<transpose<<std::endl;
                    search_iterator.add_range(
                            m_particles_begin + range_start_index,
                            m_particles_begin + range_end_index,
                            transpose);
                }

            }

            // go to next candidate bucket
            bucket_iterator++;
            
        }
        
        return search_iterator;
    }


};

   


}


#endif /* BUCKETSEARCH_H_ */
