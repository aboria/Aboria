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


#ifndef NEIGHBOUR_SEARCH_BASE_H_
#define NEIGHBOUR_SEARCH_BASE_H_

#include "detail/Algorithms.h"
#include "detail/SpatialUtil.h"
#include "Traits.h"
#include "CudaInclude.h"
#include "Vector.h"
#include "Get.h"
#include "Log.h"

namespace Aboria {


template <typename IteratorType>
struct iterator_range {
    IteratorType m_begin;
    IteratorType m_end;
    CUDA_HOST_DEVICE
    iterator_range(IteratorType&& begin, IteratorType&& end):
        m_begin(begin),m_end(end) 
    {}
    CUDA_HOST_DEVICE
    const IteratorType &begin() const { return m_begin; }
    CUDA_HOST_DEVICE
    const IteratorType &end() const { return m_end; }
};

template <typename IteratorType>
iterator_range<IteratorType> make_iterator_range(IteratorType&& begin, IteratorType&& end) {
    return iterator_range<IteratorType>(begin,end);
}

template <typename Derived, typename Traits, typename Params, typename ConstIterator, typename QueryType>
class neighbour_search_base {
public:

    typedef ConstIterator query_iterator;
    typedef QueryType query_type;
    typedef Params params_type;
    typedef typename Traits::position position;
    typedef typename Traits::double_d double_d;
    typedef typename Traits::bool_d bool_d;
    typedef typename Traits::iterator iterator;



    const Derived& cast() const { return static_cast<const Derived&>(*this); }
    Derived& cast() { return static_cast<Derived&>(*this); }

    neighbour_search_base() {
        LOG_CUDA(2,"neighbour_search_base: constructor, setting default domain");
        const double min = std::numeric_limits<double>::min();
        const double max = std::numeric_limits<double>::max();
        set_domain(double_d(min/3.0),double_d(max/3.0),bool_d(false),double_d(max/3.0-min/3.0),false); 
    };

    static constexpr bool unordered() {
        return true;
    }

    /// resets the domain extents, periodicity and bucket size
    /// \param low the lower extent of the search domain
    /// \param high the upper extent of the search domain
    /// \param _max_interaction_radius the side length of each bucket
    /// \param periodic a boolean vector indicating wether each dimension
    void set_domain(const double_d &min_in, const double_d &max_in, const bool_d& periodic_in, const Params& params, const bool not_in_constructor=true) {
        LOG_CUDA(2,"neighbour_search_base: set_domain:");
        m_bounds.bmin = min_in;
        m_bounds.bmax = max_in;
        m_periodic = periodic_in;
        if (not_in_constructor) {
            LOG(2,"\tbounds = "<<m_bounds);
            LOG(2,"\tperiodic = "<<m_periodic);
            cast().set_domain_impl(params);
        }
    }


    /// embed a set of points into the buckets, assigning each 3D point into the bucket
    /// that contains that point. Any points already assigned to the buckets are 
    /// removed. 
    /// \param begin_iterator an iterator to the beginning of the set of points
    /// \param end_iterator an iterator to the end of the set of points
    /// \see embed_points_incremental() 
    void embed_points(iterator begin, iterator end) {
        m_particles_begin = begin;
        m_particles_end = end;

        CHECK(!m_bounds.is_empty(), "trying to embed particles into an empty domain. use the function `set_domain` to setup the spatial domain first.");

        const size_t n = m_particles_end - m_particles_begin;
	    LOG(2,"neighbour_search_base: embed_points: embedding "<<n<<" points");

        cast().embed_points_impl();

    }

    void add_points_at_end(const iterator &begin, 
                           const iterator &start_adding, 
                           const iterator &end) {
        ASSERT(start_adding-begin == m_particles_end-m_particles_begin, "prior number of particles embedded into domain is not consistent with distance between begin and start_adding");

        m_particles_begin = begin;
        m_particles_end = end;

        ASSERT(!m_bounds.is_empty(), "trying to embed particles into an empty domain. use the function `set_domain` to setup the spatial domain first.");

        const size_t dist = end - start_adding;
        if (dist > 0) {
            LOG(2,"neighbour_search_base: add_points_at_end: embedding "<<dist<<" new points. Total number = "<<end-begin);
            cast().add_points_at_end_impl(dist);
        }
    }

    void copy_points(iterator copy_from_iterator, iterator copy_to_iterator) {
            ASSERT((copy_to_iterator-m_particles_begin>=0) && 
                    (m_particles_end-copy_to_iterator>0),"invalid copy to iterator");
            ASSERT((copy_from_iterator-m_particles_begin>=0) && 
                    (m_particles_end-copy_from_iterator>0),"invalid copy from iterator");
            if (copy_to_iterator==copy_from_iterator) return;
            cast().copy_points_impl(copy_from_iterator,copy_from_iterator);
    }


    void delete_points_at_end(const iterator &begin, 
                           const iterator &end) {
        ASSERT(end-begin <= m_particles_end-m_particles_begin, "prior number of particles embedded into domain is not consistent with distance between begin and start_adding");

        const size_t oldn = m_particles_end-m_particles_begin;
        m_particles_begin = begin;
        m_particles_end = end;

        ASSERT(!m_bounds.is_empty(), "trying to embed particles into an empty domain. use the function `set_domain` to setup the spatial domain first.");

        const size_t dist = oldn-(end-begin);
        if (dist > 0) {
            LOG(2,"neighbour_search_base: delete_points_at_end: deleting "<<dist<<" points. Total number = "<<end-begin);
            cast().add_points_at_end_impl(dist);
        }
    }
 

    const query_type& get_query() const {
        return cast().get_query_impl();
    }

    const double_d& get_min() const { return m_bounds.bmin; }
    const double_d& get_max() const { return m_bounds.bmax; }
    const bool_d& get_periodic() const { return m_periodic; }

protected:
    iterator m_particles_begin;
    iterator m_particles_end;
    bool_d m_periodic;
    detail::bbox<Traits::dimension> m_bounds;
};


template <typename Traits>
class ranges_iterator {
    typedef typename Traits::position position;
    typedef typename Traits::double_d double_d;
    typedef typename Traits::bool_d bool_d;
    typedef typename Traits::value_type p_value_type;
    typedef typename Traits::reference p_reference;
    typedef typename Traits::raw_pointer raw_pointer;

public:
    typedef const tuple_ns::tuple<const p_value_type&,const double_d&>* pointer;
	typedef std::forward_iterator_tag iterator_category;
    typedef const tuple_ns::tuple<p_reference,const double_d&> reference;
    typedef const tuple_ns::tuple<p_reference,const double_d&> value_type;
	typedef std::ptrdiff_t difference_type;

    CUDA_HOST_DEVICE
    ranges_iterator(const raw_pointer end):
        m_end(end),
        m_nbuckets(0),
        m_node(end)
    {}

    CUDA_HOST_DEVICE
    ranges_iterator(const raw_pointer end, const double_d box_side_length, const double_d &r):
        m_end(end),
        m_box_side_length(box_side_length),
        m_nbuckets(0),
        m_r(r),
        m_node(end)
    {}

    CUDA_HOST_DEVICE
    void add_range(raw_pointer begin, raw_pointer end, const double_d &transpose) {
#ifndef __CUDA_ARCH__
        LOG(4,"\tranges_iterator::add_range. Adding "<<end-begin<<" particles with transpose = "<<transpose<<". Number of bucket ranges already here  = "<<m_nbuckets);
#endif
        m_begins[m_nbuckets] = begin;
        m_ends[m_nbuckets] = end;
        m_transpose[m_nbuckets] = transpose;
        m_nbuckets++;

        if (m_node == m_end) {
            m_current_index = m_nbuckets-1;
            m_node = m_begins[m_current_index];
            if (!check_candidate()) {
                increment(); 
            }
        }
    }

    CUDA_HOST_DEVICE
    bool equal(ranges_iterator const& other) const {
        return m_node == other.m_node;
    }

    CUDA_HOST_DEVICE
    reference dereference() const { 
        return reference(*m_node,m_dx); 
    }

    CUDA_HOST_DEVICE
    bool go_to_next_candidate() {
        m_node++;
        if (m_node == m_ends[m_current_index]) {
            m_current_index++;
            //std::cout << "moving on to next index i = "<<m_current_index<<" with range "<<m_begins[m_current_index]-m_bucket_sort->m_particles_begin<<" to "<<m_ends[m_current_index]-m_bucket_sort->m_particles_begin<<std::endl;
            if (m_current_index < m_nbuckets) {
                m_node = m_begins[m_current_index];
                //std::cout << "particle index = "<<m_node-m_bucket_sort->m_particles_begin<<std::endl;
            } else {
                m_node = m_end;
                return false;
            }
        }
        return true;
    }

    CUDA_HOST_DEVICE
    bool check_candidate() {
        //std::cout << "check my_r = "<<m_r<<" r = "<<get<position>(*m_node)<<" trans = "<<m_transpose[m_current_index]<<" index = "<<m_current_index<<std::endl;
        const double_d p = get<position>(*m_node) + m_transpose[m_current_index];
        m_dx = p - m_r;

        bool outside = false;
        for (int i=0; i < Traits::dimension; i++) {
            if (std::abs(m_dx[i]) > m_box_side_length[i]) {
                outside = true;
                break;
            } 
        }

        return !outside;
    }

    CUDA_HOST_DEVICE
    void increment() {
        bool found_good_candidate = false;
        while (!found_good_candidate && go_to_next_candidate()) {
            found_good_candidate = check_candidate();
        }
    }

    CUDA_HOST_DEVICE
    reference operator *() const {
        return dereference();
    }

    CUDA_HOST_DEVICE
    reference operator ->() const {
        return dereference();
    }

    CUDA_HOST_DEVICE
    ranges_iterator& operator++() {
        increment();
        return *this;
    }

    CUDA_HOST_DEVICE
    ranges_iterator operator++(int) {
        ranges_iterator tmp(*this);
        operator++();
        return tmp;
    }

    CUDA_HOST_DEVICE
    size_t operator-(ranges_iterator start) const {
        size_t count = 0;
        while (start != *this) {
            ++start; ++count;
        }
        return count;
    }

    CUDA_HOST_DEVICE
    inline bool operator==(const ranges_iterator& rhs) {
        return equal(rhs);
    }

    CUDA_HOST_DEVICE
    inline bool operator!=(const ranges_iterator& rhs){
        return !operator==(rhs);
    }

private:
    friend class boost::iterator_core_access;

    raw_pointer m_end;
    double_d m_box_side_length;
    
    double_d m_r;
    double_d m_dx;
    double_d m_search_side;
    raw_pointer m_node;
    raw_pointer m_node_end;
    
    const static unsigned int max_nbuckets = detail::ipow(3,Traits::dimension); 
    unsigned int m_nbuckets; 
    raw_pointer m_begins[max_nbuckets];
    raw_pointer m_ends[max_nbuckets];
    double_d m_transpose[max_nbuckets];
    int m_current_index = -1;
};

/// A const iterator to a set of neighbouring points. This iterator implements
/// a STL forward iterator type
template <typename Traits>
class linked_list_iterator {
    typedef typename Traits::position position;
    typedef typename Traits::double_d double_d;
    typedef typename Traits::bool_d bool_d;
    typedef typename Traits::value_type p_value_type;
    typedef typename Traits::reference p_reference;
    typedef typename Traits::raw_pointer raw_pointer;

public:
    typedef const tuple_ns::tuple<const p_value_type&,const double_d&>* pointer;
	typedef std::forward_iterator_tag iterator_category;
    typedef const tuple_ns::tuple<p_reference,const double_d&> reference;
    typedef const tuple_ns::tuple<p_reference,const double_d&> value_type;
	typedef std::ptrdiff_t difference_type;

    linked_list_iterator(): 
        m_cell_empty(detail::get_empty_id()) {
        m_node = &m_cell_empty;
    }

   
    /// this constructor is used to start the iterator at the head of a bucket 
    /// list
    /// \param bucket_sort a pointer to the parent BucketSearch class
    /// \param centre the neighbourhood query point
    linked_list_iterator(
            const raw_pointer begin,
            const double_d box_side_length, 
            size_t* linked_list_begin,
            size_t* buckets_begin,
            const double_d& r):
        m_begin(begin),
        m_box_side_length(box_side_length),
        m_linked_list_begin(linked_list_begin),
        m_buckets_begin(buckets_begin),
        m_r(r),
        m_bucket_i(0),
        m_nbuckets(0),
        m_cell_empty(detail::get_empty_id()) {

        m_node = &m_cell_empty;
    }

    void add_bucket(const size_t bucket_index, const double_d& transpose) {
        m_buckets_to_search[m_nbuckets] = bucket_index;
        m_transpose[m_nbuckets] = transpose;
        ++m_nbuckets;

        if (*m_node == detail::get_empty_id()) {
            increment();
        }
    }

     bool go_to_next_candidate() {
        if (*m_node != detail::get_empty_id()) {
            m_node = m_linked_list_begin + *m_node;
            //std::cout << "going to new particle *mnode = "<<*m_node<<std::endl;
        }
        while ((*m_node == detail::get_empty_id()) && (m_bucket_i < m_nbuckets)) {
            //std::cout << "going to new_cell with offset = "<<*surrounding_cell_offset_i<<std::endl;
            m_node = m_buckets_begin + m_buckets_to_search[m_bucket_i];
            ++m_bucket_i;
        }
        if (m_bucket_i == m_nbuckets) {
            return false;
        } else {
            return true;
        }
    }

    
    reference operator *() const {
        return dereference();
    }
    reference operator ->() {
        return dereference();
    }
    linked_list_iterator& operator++() {
        increment();
        return *this;
    }
    linked_list_iterator operator++(int) {
        linked_list_iterator tmp(*this);
        operator++();
        return tmp;
    }
    size_t operator-(linked_list_iterator start) const {
        size_t count = 0;
        while (start != *this) {
            start++;
            count++;
        }
        return count;
    }
    inline bool operator==(const linked_list_iterator& rhs) {
        return equal(rhs);
    }
    inline bool operator!=(const linked_list_iterator& rhs){
        return !operator==(rhs);
    }

 private:
    friend class boost::iterator_core_access;

    bool equal(linked_list_iterator const& other) const {
        //std::cout <<" testing equal *m_node = "<<*m_node<<" other.m_node = "<<*(other.m_node)<<std::endl;
        return *m_node == *(other.m_node);
    }


    bool check_candidate() {
        const double_d& p = get<position>(m_begin)[*m_node] + m_transpose[m_bucket_i];
        m_dx = p - m_r;

        bool outside = false;
        for (int i=0; i < Traits::dimension; i++) {
            if (std::abs(m_dx[i]) > m_box_side_length[i]) {
                outside = true;
                break;
            } 
        }
        return !outside;

    }

    void increment() {
        bool found_good_candidate = false;
        while (!found_good_candidate && go_to_next_candidate()) {
            found_good_candidate = check_candidate();
        }
    }


    reference dereference() const
    { return reference(*(m_begin + *m_node),m_dx); }


    size_t* m_node;
    double_d m_r;
    double_d m_dx;
    double_d m_box_side_length;

    const static unsigned int max_nbuckets = detail::ipow(3,Traits::dimension); 
    double_d m_transpose[max_nbuckets];
    size_t m_buckets_to_search[max_nbuckets];
    size_t m_nbuckets;
    size_t m_empty_bucket;
    size_t m_bucket_i;
    size_t m_cell_empty;

    raw_pointer m_begin;
    size_t* m_buckets_begin;
    size_t* m_linked_list_begin;

};



}


#endif /* BUCKETSEARCH_H_ */
