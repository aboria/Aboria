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
#include "NeighbourSearchBase.h"
#include "Traits.h"
#include "CudaInclude.h"
#include "Vector.h"
#include "Get.h"

#include <iostream>
#include <queue>
#include "Log.h"

namespace Aboria {

/// A const iterator to a set of neighbouring points. This iterator implements
/// a STL forward iterator type
// assume that these iterators, and query functions, are only called from device code
template <typename Query>
class box_search_iterator {

    typedef typename Query::particle_iterator particle_iterator;
    typedef typename Query::bucket_iterator bucket_iterator;
    typedef typename Query::traits_type Traits;

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
    iterator_range<bucket_iterator> m_bucket_range;
    bucket_iterator m_current_bucket;
    iterator_range_with_transpose<particle_iterator> m_particle_range;
    particle_iterator m_current_particle;

public:
    typedef const tuple_ns::tuple<p_reference,const double_d&>* pointer;
	typedef std::forward_iterator_tag iterator_category;
    typedef const tuple_ns::tuple<p_reference,const double_d&> reference;
    typedef const tuple_ns::tuple<p_reference,const double_d&> value_type;
	typedef std::ptrdiff_t difference_type;

    CUDA_HOST_DEVICE
    box_search_iterator():
        m_valid(false)
    {}

    CUDA_HOST_DEVICE
    box_search_iterator(const Query &query,const double_d &r):
        m_valid(true),
        m_r(r),
        m_query(&query),
        m_bucket_range(query.get_near_buckets(query.get_bucket(r))),
        m_current_bucket(m_bucket_range.begin()),
        m_particle_range(query.get_bucket_particles(*m_current_bucket)),
        m_current_particle(m_particle_range.begin())
    {
        get_valid_candidate();
        if (!check_candidate()) {
            increment();
        }
        LOG(4,"\tconstructor (box_search_iterator): r = "<<m_r<<" m_current_bucket = "<<*m_current_bucket<<"m_current_particle position= "<<get<position>(*m_current_particle));
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
    box_search_iterator& operator++() {
        increment();
        return *this;
    }
    CUDA_HOST_DEVICE
    box_search_iterator operator++(int) {
        box_search_iterator tmp(*this);
        operator++();
        return tmp;
    }
    CUDA_HOST_DEVICE
    size_t operator-(box_search_iterator start) const {
        size_t count = 0;
        while (start != *this) {
            start++;
            count++;
        }
        return count;
    }
    CUDA_HOST_DEVICE
    inline bool operator==(const box_search_iterator& rhs) {
        return equal(rhs);
    }
    CUDA_HOST_DEVICE
    inline bool operator!=(const box_search_iterator& rhs){
        return !operator==(rhs);
    }

 private:
    friend class boost::iterator_core_access;


    CUDA_HOST_DEVICE
    bool equal(box_search_iterator const& other) const {
        return m_valid ? 
                    m_current_particle == other.m_current_particle
                    : 
                    !other.m_valid;
    }

    CUDA_HOST_DEVICE
    bool get_valid_candidate() {
        while (m_current_particle == m_particle_range.end()) {
            ++m_current_bucket;
            if (m_current_bucket == m_bucket_range.end()) {
                m_valid = false;
                break; 
            }
            m_particle_range = m_query->get_bucket_particles(*m_current_bucket);
            m_current_particle = m_particle_range.begin();
        }
        return m_valid;
    }

    CUDA_HOST_DEVICE
    bool go_to_next_candidate() {
#ifndef __CUDA_ARCH__
        LOG(4,"\tgo_to_next_candidate (box_search_iterator):"); 
#endif
        ++m_current_particle;
        return get_valid_candidate();
    }

    CUDA_HOST_DEVICE
    bool check_candidate() {
        //const double_d& p = get<position>(*m_current_particle) + m_particle_range.get_transpose();
        const double_d& p = get<position>(*m_current_particle); 
        const double_d& transpose = m_particle_range.get_transpose();
        bool outside = false;
        for (int i=0; i < Traits::dimension; i++) {
            m_dx[i] = p[i] + transpose[i] - m_r[i];
            if (std::abs(m_dx[i]) > m_query->get_min_bucket_size()[i]) {
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
        LOG(4,"\tincrement (box_search_iterator):"); 
#endif
        bool found_good_candidate = false;
        while (!found_good_candidate && go_to_next_candidate()) {
            found_good_candidate = check_candidate();
#ifndef __CUDA_ARCH__
            LOG(4,"\tfound_good_candidate = "<<found_good_candidate); 
#endif
        }
#ifndef __CUDA_ARCH__
        LOG(4,"\tend increment (box_search_iterator): invalid = " << m_valid); 
#endif
    }


    CUDA_HOST_DEVICE
    reference dereference() const
    { return reference(*m_current_particle,m_dx); }

    
};



template<typename Query,
         typename SearchIterator = box_search_iterator<Query>>
iterator_range<SearchIterator> 
box_search(const Query& query, 
           const typename Query::double_d& box_centre) {
    return iterator_range<SearchIterator>(
                 SearchIterator(query,box_centre)
                ,SearchIterator()
            );
}


}

#endif
