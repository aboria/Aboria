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


#ifndef SEARCH_H_
#define SEARCH_H_

#include "detail/Algorithms.h"
#include "detail/SpatialUtil.h"
#include "detail/Distance.h"
#include "NeighbourSearchBase.h"
#include "Traits.h"
#include "CudaInclude.h"
#include "Vector.h"
#include "Get.h"

#include <iostream>
#include <queue>
#include <cmath>
#include "Log.h"

namespace Aboria {

/// A const iterator to a set of neighbouring points. This iterator implements
/// a STL forward iterator type
// assume that these iterators, and query functions, are only called from device code
template <typename Query, int LNormNumber>
class search_iterator {

    typedef typename Query::particle_iterator particle_iterator;
    typedef typename Query::query_iterator query_iterator;
    typedef typename Query::traits_type Traits;
    static const unsigned int dimension = Traits::dimension;

    typedef typename Traits::position position;
    typedef typename Traits::double_d double_d;
    typedef typename Traits::bool_d bool_d;
    typedef typename particle_iterator::value_type p_value_type;
    typedef typename particle_iterator::reference p_reference;
    typedef typename particle_iterator::pointer p_pointer;

    bool m_valid;
    double_d m_r;
    double_d m_dx;
    const Query *m_query;
    double m_max_distance;
    double m_max_distance2;
    iterator_range<query_iterator> m_bucket_range;
    query_iterator m_current_bucket;
    iterator_range_with_transpose<particle_iterator> m_particle_range;
    particle_iterator m_current_particle;

public:
    typedef const tuple_ns::tuple<p_reference,const double_d&>* pointer;
	typedef std::forward_iterator_tag iterator_category;
    typedef const tuple_ns::tuple<p_reference,const double_d&> reference;
    typedef const tuple_ns::tuple<p_reference,const double_d&> value_type;
	typedef std::ptrdiff_t difference_type;

    CUDA_HOST_DEVICE
    search_iterator():
        m_valid(false)
    {}

    CUDA_HOST_DEVICE
    search_iterator(const Query &query,
                    const double_d &r,
                    const double max_distance):
        m_valid(true),
        m_r(r),
        m_query(&query),
        m_max_distance(max_distance),
        m_max_distance2(detail::distance_helper<LNormNumber>::get_value_to_accumulate(max_distance)),
        m_bucket_range(query.get_buckets_near_point(r,max_distance)),
        m_current_bucket(m_bucket_range.begin()),
        m_particle_range(query.get_bucket_particles(*m_current_bucket)),
        m_current_particle(m_particle_range.begin())
    {
        get_valid_candidate();
        if (m_valid && !check_candidate()) {
            increment();
        }
        LOG(4,"\tconstructor (search_iterator): r = "<<m_r<<"m_current_particle position= "<<get<position>(*m_current_particle));
    }
    
    CUDA_HOST_DEVICE
    reference operator *() const {
        return dereference();
    }
    CUDA_HOST_DEVICE
    reference operator ->() {
        return dereference();
    }
    CUDA_HOST_DEVICE
    search_iterator& operator++() {
        increment();
        return *this;
    }
    CUDA_HOST_DEVICE
    search_iterator operator++(int) {
        search_iterator tmp(*this);
        operator++();
        return tmp;
    }
    CUDA_HOST_DEVICE
    size_t operator-(search_iterator start) const {
        size_t count = 0;
        while (start != *this) {
            start++;
            count++;
        }
        return count;
    }
    CUDA_HOST_DEVICE
    inline bool operator==(const search_iterator& rhs) {
        return equal(rhs);
    }
    CUDA_HOST_DEVICE
    inline bool operator!=(const search_iterator& rhs){
        return !operator==(rhs);
    }

 private:
    friend class boost::iterator_core_access;


    CUDA_HOST_DEVICE
    bool equal(search_iterator const& other) const {
        return m_valid ? 
                    m_current_particle == other.m_current_particle
                    : 
                    !other.m_valid;
    }


    /*
    bool get_valid_bucket() {
        for (; m_current_bucket != m_bucket_range.end(); ++m_current_bucket) {
            detail::bbox<dimension> bbox = m_query->get_bucket_bbox(*m_current_bucket);
            double accum = 0;
            bool intersect = true;
            for (int i = 0; i < dimension; ++i) {
                const double low = m_r[i]-m_max_distance;
                const double high = m_r[i]+m_max_distance;
                if ( (low < bbox.bmin[i] && high < bbox.bmin[i]) 
                   ||(low > bbox.bmax[i] && high > bbox.bmax[i])) { 
                    intersect = false;
                    break;
                }
            }
            if (intersect) {
                return true;
            }
        }
        // must have exhausted buckets
        return false;
    }
    */

    CUDA_HOST_DEVICE
    void get_valid_candidate() {
        while (m_current_particle == m_particle_range.end()) {
            ++m_current_bucket;
            if (m_current_bucket == m_bucket_range.end()) {
#ifndef __CUDA_ARCH__
                LOG(4,"\tran out of buckets to search (search_iterator): m_current_bucket = "<<*m_current_bucket); 
#endif
                m_valid = false;
                break; 
            }
#ifndef __CUDA_ARCH__
            LOG(4,"\tgo_to_next bucket (search_iterator): new bucket = "<<*m_current_bucket); 
#endif
            m_particle_range = m_query->get_bucket_particles(*m_current_bucket);
            m_current_particle = m_particle_range.begin();
        }
    }

    CUDA_HOST_DEVICE
    void go_to_next_candidate() {
#ifndef __CUDA_ARCH__
        LOG(4,"\tgo_to_next_candidate (search_iterator):"); 
#endif
        ++m_current_particle;
        get_valid_candidate();
    }

    

    CUDA_HOST_DEVICE
    bool check_candidate() {
        //const double_d& p = get<position>(*m_current_particle) + m_particle_range.get_transpose();
        const double_d& p = get<position>(*m_current_particle); 
        const double_d& transpose = m_particle_range.get_transpose();
        double accum = 0;
        bool outside = false;
        for (int i=0; i < Traits::dimension; i++) {
            m_dx[i] = p[i] + transpose[i] - m_r[i];
            accum = detail::distance_helper<LNormNumber>::accumulate_norm(accum, m_dx[i]);
            if (accum > m_max_distance2) {
                outside = true;
                break;
            } 
        }
#ifndef __CUDA_ARCH__
        LOG(4,"\tcheck_candidate: m_r = "<<m_r<<" other r = "<<get<position>(*m_current_particle)<<" trans = "<<m_particle_range.get_transpose()<<". outside = "<<outside); 
#endif
        return !outside;
    }

    CUDA_HOST_DEVICE
    void increment() {
#ifndef __CUDA_ARCH__
        LOG(4,"\tincrement (search_iterator):"); 
#endif
        bool found_good_candidate = false;
        while (!found_good_candidate && m_valid) {
            go_to_next_candidate();
            if (m_valid) {
                found_good_candidate = check_candidate();
            }
#ifndef __CUDA_ARCH__
            LOG(4,"\tfound_good_candidate = "<<found_good_candidate); 
#endif
            
        }
#ifndef __CUDA_ARCH__
        LOG(4,"\tend increment (search_iterator): invalid = " << m_valid); 
#endif
    }


    CUDA_HOST_DEVICE
    reference dereference() const
    { return reference(*m_current_particle,m_dx); }

    
};


/*
template <typename query_iterator>
iterator_range<query_iterator> 
get_buckets_near_point(const double_d &position, const double max_distance, detail::cell_list_tag) {
}

template <typename query_iterator>
iterator_range<query_iterator> 
get_buckets_near_point(const double_d &position, const double max_distance, detail::kd_tree_tag) {
}
*/


template<typename Query,
         typename SearchIterator = search_iterator<Query,-1>>
iterator_range<SearchIterator> 
chebyshev_search(const Query& query, 
           const typename Query::double_d& centre,
           const double max_distance) {
    return iterator_range<SearchIterator>(
                 SearchIterator(query,centre,max_distance)
                ,SearchIterator()
            );
}

template<typename Query,
         typename SearchIterator = search_iterator<Query,1>>
iterator_range<SearchIterator> 
manhatten_search(const Query& query, 
           const typename Query::double_d& centre,
           const double max_distance) {
    return iterator_range<SearchIterator>(
                 SearchIterator(query,centre,max_distance)
                ,SearchIterator()
            );
}

template<typename Query,
         typename SearchIterator = search_iterator<Query,2>>
iterator_range<SearchIterator> 
euclidean_search(const Query& query, 
           const typename Query::double_d& centre,
           const double max_distance) {
    return iterator_range<SearchIterator>(
                 SearchIterator(query,centre,max_distance)
                ,SearchIterator()
            );
}

template<int LNormNumber,
         typename Query,
         typename SearchIterator = search_iterator<Query,LNormNumber>>
iterator_range<SearchIterator> 
distance_search(const Query& query, 
           const typename Query::double_d& centre,
           const double max_distance) {
    return iterator_range<SearchIterator>(
                 SearchIterator(query,centre,max_distance)
                ,SearchIterator()
            );
}


}

#endif
