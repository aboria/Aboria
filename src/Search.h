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
template <typename Iterator>
class box_search_iterator {
    typedef typename Iterator::traits_type Traits;
    typedef typename Traits::position position;
    typedef typename Traits::double_d double_d;
    typedef typename Traits::bool_d bool_d;
    typedef typename Iterator::value_type p_value_type;
    typedef typename Iterator::reference p_reference;
    typedef typename Iterator::pointer p_pointer;

public:
    typedef const tuple_ns::tuple<p_reference,const double_d&>* pointer;
	typedef std::forward_iterator_tag iterator_category;
    typedef const tuple_ns::tuple<p_reference,const double_d&> reference;
    typedef const tuple_ns::tuple<p_reference,const double_d&> value_type;
	typedef std::ptrdiff_t difference_type;

    CUDA_HOST_DEVICE
    box_search_iterator()
    {}

    CUDA_HOST_DEVICE
    box_search_iterator(const double_d &r, const double_d &box_side_length):
        m_r(r),
        m_box_side_length(box_side_length)
    {
    }
   

    CUDA_HOST_DEVICE
    void add_range(iterator_range<Iterator> &&range) {
        if (range.begin()!=range.end()) {
            m_buckets_to_search.push(std::move(range));
            if (m_buckets_to_search.size()==1) {
                // make sure we have a good candidate ready to go
                // if none in range then increment makes sure we 
                // go back to being invalid
                m_current_p = range.begin();
                if (!check_candidate()) {
                    increment();
                }
            }
        }
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
    bool is_invalid() const {
        return m_buckets_to_search.empty();
    }

    CUDA_HOST_DEVICE
    bool equal(box_search_iterator const& other) const {
        return is_invalid() ? 
                    other.is_invalid() 
                    : 
                    m_current_p == other.m_current_p;
    }

    CUDA_HOST_DEVICE
    bool go_to_next_candidate() {
        if (is_invalid()) return false;
#ifndef __CUDA_ARCH__
        LOG(4,"\tgo_to_next_candidate (box_search_iterator):"); 
#endif

        ++m_current_p;
        if (m_current_p == m_buckets_to_search.front().end()) {
            m_buckets_to_search.pop();

#ifndef __CUDA_ARCH__
            LOG(4,"\tend of range, moving to next range "); 
#endif
            if (!m_buckets_to_search.empty()) {
                ASSERT(m_buckets_to_search.front().begin() != m_buckets_to_search.front().end(),"error, empty range found");
                m_current_p = m_buckets_to_search.front().begin();
            } else {

#ifndef __CUDA_ARCH__
                LOG(4,"\tfinished ranges"); 
#endif
                return false;
            }
        }
        return true;
    }

    CUDA_HOST_DEVICE
    bool check_candidate() {
        const double_d& p = get<position>(*m_current_p) + m_current_p.get_transpose();
        m_dx = p - m_r;

        bool outside = false;
        for (int i=0; i < Traits::dimension; i++) {
            if (std::abs(m_dx[i]) > m_box_side_length[i]) {
                outside = true;
                break;
            } 
        }
#ifndef __CUDA_ARCH__
        LOG(4,"\tcheck_candidate: m_r = "<<m_r<<" other r = "<<p<<" trans = "<<m_current_p.get_transpose()<<". m_box_side_length = "<<m_box_side_length<<". outside = "<<outside); 
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
        LOG(4,"\tend increment (box_search_iterator)"); 
#endif
    }


    CUDA_HOST_DEVICE
    reference dereference() const
    { return reference(*m_current_p,m_dx); }


    double_d m_r;
    double_d m_dx;
    double_d m_box_side_length;


    const static unsigned int max_nbuckets = detail::ipow(3,Traits::dimension); 
    std::queue<iterator_range<Iterator>> m_buckets_to_search;
    Iterator m_current_p;

};

template<typename Query>
iterator_range<box_search_iterator<typename Query::particle_iterator>> 
box_search(const Query query, 
           const typename Query::double_d box_centre, 
           const typename Query::double_d box_sides) {

    ASSERT(box_sides <= query.get_min_bucket_size(),"box query with greater than neighbour search min box size not currently supported");
    typedef box_search_iterator<typename Query::particle_iterator> search_iterator;

    iterator_range<search_iterator> search_range(
            search_iterator(box_centre,box_sides)
            ,search_iterator()
            );
    for (const auto &i: query.get_near_buckets(query.get_bucket(box_centre))) {
        search_range.begin().add_range(query.get_bucket_particles(i));
    }
    return search_range;
}

template<typename Query>
iterator_range<box_search_iterator<typename Query::particle_iterator>> 
box_search(const Query query, 
           const typename Query::double_d box_centre) {

    typedef box_search_iterator<typename Query::particle_iterator> search_iterator;

    iterator_range<search_iterator> search_range(
            search_iterator(box_centre,query.get_min_bucket_size())
            ,search_iterator()
            );
    for (const auto &i: query.get_near_buckets(query.get_bucket(box_centre))) {
        search_range.begin().add_range(query.get_bucket_particles(i));
    }
    return search_range;
}


}

#endif
