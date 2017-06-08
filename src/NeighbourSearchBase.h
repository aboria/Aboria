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
#include "detail/Distance.h"
#include "Traits.h"
#include "CudaInclude.h"
#include "Vector.h"
#include "Get.h"
#include "Log.h"
#include <stack>

namespace Aboria {


template <typename IteratorType>
struct iterator_range_with_transpose {
    typedef IteratorType iterator;
    typedef typename iterator::traits_type traits_type;
    typedef typename traits_type::double_d double_d;
    IteratorType m_begin;
    IteratorType m_end;
    CUDA_HOST_DEVICE
    iterator_range_with_transpose()
    {}
    /*
    CUDA_HOST_DEVICE
    iterator_range_with_transpose(IteratorType&& begin, IteratorType&& end, const double_d& transpose):
        m_begin(std::move(begin)),m_end(std::move(end)),m_transpose(transpose) 
    {}
    */
    CUDA_HOST_DEVICE
    iterator_range_with_transpose(const IteratorType& begin, const IteratorType& end, const double_d &transpose):
        m_begin(begin),m_end(end),m_transpose(transpose) 
    {}
    CUDA_HOST_DEVICE
    iterator_range_with_transpose(const IteratorType& begin, const IteratorType& end):
        m_begin(begin),m_end(end),m_transpose(0) 
    {}
    CUDA_HOST_DEVICE
    const IteratorType &begin() const { return m_begin; }
    CUDA_HOST_DEVICE
    const IteratorType &end() const { return m_end; }
    CUDA_HOST_DEVICE
    IteratorType &begin() { return m_begin; }
    CUDA_HOST_DEVICE
    IteratorType &end() { return m_end; }

    CUDA_HOST_DEVICE
    const double_d& get_transpose() { return m_transpose; }
    double_d m_transpose;

};

template <typename IteratorType>
struct iterator_range{
    typedef IteratorType iterator;
    IteratorType m_begin;
    IteratorType m_end;
    CUDA_HOST_DEVICE
    iterator_range()
    {}
    CUDA_HOST_DEVICE
    iterator_range(const IteratorType& begin, const IteratorType& end):
        m_begin(begin),m_end(end)
    {}
    CUDA_HOST_DEVICE
    const IteratorType &begin() const { return m_begin; }
    CUDA_HOST_DEVICE
    const IteratorType &end() const { return m_end; }
    CUDA_HOST_DEVICE
    IteratorType &begin() { return m_begin; }
    CUDA_HOST_DEVICE
    IteratorType &end() { return m_end; }
};

template <typename IteratorType>
iterator_range<IteratorType> make_iterator_range(IteratorType&& begin, IteratorType&& end) {
    return iterator_range<IteratorType>(begin,end);
}

template <typename Derived, typename Traits, typename QueryType>
class neighbour_search_base {
public:

    typedef QueryType query_type;
    typedef typename Traits::double_d double_d;
    typedef typename Traits::bool_d bool_d;
    typedef typename Traits::iterator iterator;



    const Derived& cast() const { return static_cast<const Derived&>(*this); }
    Derived& cast() { return static_cast<Derived&>(*this); }

    neighbour_search_base() {
        LOG_CUDA(2,"neighbour_search_base: constructor, setting default domain");
        const double min = std::numeric_limits<double>::min();
        const double max = std::numeric_limits<double>::max();
        set_domain(double_d(min/3.0),double_d(max/3.0),bool_d(false),10,false); 
    };

    static constexpr bool cheap_copy_and_delete_at_end() {
        return true;
    }

    /// resets the domain extents, periodicity and bucket size
    /// \param low the lower extent of the search domain
    /// \param high the upper extent of the search domain
    /// \param _max_interaction_radius the side length of each bucket
    /// \param periodic a boolean vector indicating wether each dimension
    void set_domain(const double_d &min_in, const double_d &max_in, const bool_d& periodic_in, const unsigned int n_particles_in_leaf=10, const bool not_in_constructor=true) {
        LOG(2,"neighbour_search_base: set_domain:");
        m_bounds.bmin = min_in;
        m_bounds.bmax = max_in;
        m_periodic = periodic_in;
        m_n_particles_in_leaf = n_particles_in_leaf; 
        if (not_in_constructor) {
            cast().set_domain_impl();
        }
        LOG(2,"\tbounds = "<<m_bounds);
	    LOG(2,"\tparticles_in_leaf = "<<m_n_particles_in_leaf);
        LOG(2,"\tperiodic = "<<m_periodic);
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

    void update_iterators(iterator begin, iterator end) {
        m_particles_begin = begin;
        m_particles_end = end;

	    LOG(2,"neighbour_search_base: update iterators");
        cast().update_iterator_impl();
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
            LOG(2,"neighbour_search_base: add_points_at_end: done.");

        }
        ASSERT(m_particles_begin==begin, "did not update m_particles_begin correctly");
        ASSERT(m_particles_end==end, "did not update m_particles_end correctly");
    }

    void copy_points(iterator copy_from_iterator, iterator copy_to_iterator) {
        LOG(4,"neighbour_search_base: copy_points: fromi = "<<copy_from_iterator-m_particles_begin<<" toi = "<<copy_to_iterator-m_particles_begin);
        ASSERT((copy_to_iterator-m_particles_begin>=0) && 
                (m_particles_end-copy_to_iterator>0),"invalid copy to iterator");
        ASSERT((copy_from_iterator-m_particles_begin>=0) && 
                (m_particles_end-copy_from_iterator>0),"invalid copy from iterator");
        if (copy_to_iterator==copy_from_iterator) return;
        cast().copy_points_impl(copy_from_iterator,copy_to_iterator);
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
            cast().delete_points_at_end_impl(dist);
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
    unsigned int m_n_particles_in_leaf; 
};



// assume that these iterators, and query functions, are only called from device code
template <typename Traits>
class ranges_iterator {
    typedef typename Traits::double_d double_d;
    typedef typename Traits::bool_d bool_d;
    typedef typename Traits::value_type p_value_type;
    typedef typename Traits::raw_reference p_reference;
    typedef typename Traits::raw_pointer p_pointer;

public:
    typedef Traits traits_type;
    typedef const p_pointer pointer;
	typedef std::random_access_iterator_tag iterator_category;
    typedef const p_reference reference;
    typedef const p_reference value_type;
	typedef std::ptrdiff_t difference_type;

    CUDA_HOST_DEVICE
    ranges_iterator() {}

    CUDA_HOST_DEVICE
    ranges_iterator(const p_pointer& begin):
        m_current_p(begin)
    {}

        
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
        return tuple_ns::get<0>(m_current_p.get_tuple()) 
                - tuple_ns::get<0>(start.m_current_p.get_tuple());
    }

    CUDA_HOST_DEVICE
    inline bool operator==(const ranges_iterator& rhs) const {
        return equal(rhs);
    }

    CUDA_HOST_DEVICE
    inline bool operator!=(const ranges_iterator& rhs) const {
        return !operator==(rhs);
    }

private:
    friend class boost::iterator_core_access;

    CUDA_HOST_DEVICE
    bool equal(ranges_iterator const& other) const {
        return m_current_p == other.m_current_p;
    }

    CUDA_HOST_DEVICE
    reference dereference() const { 
        return *m_current_p; 
    }

    CUDA_HOST_DEVICE
    void increment() {
        ++m_current_p;
    }

    p_pointer m_current_p;
};

/// A const iterator to a set of neighbouring points. This iterator implements
/// a STL forward iterator type
// assume that these iterators, and query functions, are only called from device code
template <typename Traits>
class linked_list_iterator {
    typedef typename Traits::double_d double_d;
    typedef typename Traits::bool_d bool_d;
    typedef typename Traits::value_type p_value_type;
    typedef typename Traits::raw_reference p_reference;
    typedef typename Traits::raw_pointer p_pointer;

public:
    typedef Traits traits_type;
    typedef const p_pointer pointer;
	typedef std::forward_iterator_tag iterator_category;
    typedef const p_reference reference;
    typedef const p_reference value_type;
	typedef std::ptrdiff_t difference_type;

    CUDA_HOST_DEVICE
    linked_list_iterator(): 
        m_current_index(detail::get_empty_id()) {
    }

   
    /// this constructor is used to start the iterator at the head of a bucket 
    /// list
    CUDA_HOST_DEVICE
    linked_list_iterator(
            const int index,
            const p_pointer& particles_begin,
            int* const linked_list_begin):
        m_current_index(index),
        m_particles_begin(particles_begin),
        m_linked_list_begin(linked_list_begin)
    {}

    CUDA_HOST_DEVICE
    linked_list_iterator(const linked_list_iterator& other):
        m_current_index(other.m_current_index),
        m_particles_begin(other.m_particles_begin),
        m_linked_list_begin(other.m_linked_list_begin)
    {}

    CUDA_HOST_DEVICE
    void operator=(const linked_list_iterator& other) {
        m_current_index = other.m_current_index;
        if (tuple_ns::get<0>(m_particles_begin.get_tuple()) != 
            tuple_ns::get<0>(other.m_particles_begin.get_tuple())) {
            m_particles_begin = other.m_particles_begin;
        }
        m_linked_list_begin = other.m_linked_list_begin;
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
    linked_list_iterator& operator++() {
        increment();
        return *this;
    }
    CUDA_HOST_DEVICE
    linked_list_iterator operator++(int) {
        linked_list_iterator tmp(*this);
        operator++();
        return tmp;
    }
    CUDA_HOST_DEVICE
    size_t operator-(linked_list_iterator start) const {
        size_t count = 0;
        while (start != *this) {
            start++;
            count++;
        }
        return count;
    }
    CUDA_HOST_DEVICE
    inline bool operator==(const linked_list_iterator& rhs) const {
        return equal(rhs);
    }
    CUDA_HOST_DEVICE
    inline bool operator!=(const linked_list_iterator& rhs) const {
        return !operator==(rhs);
    }

 private:
    friend class boost::iterator_core_access;

    CUDA_HOST_DEVICE
    bool increment() {
#ifndef __CUDA_ARCH__
        LOG(4,"\tincrement (linked_list_iterator):"); 
#endif
        if (m_current_index != detail::get_empty_id()) {
            m_current_index = m_linked_list_begin[m_current_index];
#ifndef __CUDA_ARCH__
            LOG(4,"\tgoing to new particle m_current_index = "<<m_current_index);
#endif
        } 
        
        if (m_current_index == detail::get_empty_id()) {
            return false;
        } else {
            return true;
        }
#ifndef __CUDA_ARCH__
        LOG(4,"\tend increment (linked_list_iterator): m_current_index = "<<m_current_index); 
#endif
    }

    CUDA_HOST_DEVICE
    bool equal(linked_list_iterator const& other) const {
        return m_current_index == other.m_current_index;
    }


    CUDA_HOST_DEVICE
    reference dereference() const
    { return *(m_particles_begin + m_current_index); }


    int m_current_index;
    p_pointer m_particles_begin;
    int* m_linked_list_begin;

};


template <typename Traits, typename Iterator>
class index_vector_iterator {
    typedef index_vector_iterator<Traits,Iterator> iterator;
    typedef typename Traits::double_d double_d;
    typedef typename Traits::bool_d bool_d;
    typedef typename Traits::value_type p_value_type;
    typedef typename Traits::raw_reference p_reference;
    typedef typename Traits::raw_pointer p_pointer;

public:
    typedef Traits traits_type;
    typedef const p_pointer pointer;
	typedef std::forward_iterator_tag iterator_category;
    typedef const p_reference reference;
    typedef const p_reference value_type;
	typedef std::ptrdiff_t difference_type;

    CUDA_HOST_DEVICE
    index_vector_iterator() 
    {}
       
    /// this constructor is used to start the iterator at the head of a bucket 
    /// list
    CUDA_HOST_DEVICE
    index_vector_iterator(
            Iterator begin,
            const p_pointer& particles_begin):
        m_current_index(begin),
        m_particles_begin(particles_begin)
    {}


    CUDA_HOST_DEVICE
    reference operator *() const {
        return dereference();
    }
    CUDA_HOST_DEVICE
    reference operator ->() {
        return dereference();
    }
    CUDA_HOST_DEVICE
    iterator& operator++() {
        increment();
        return *this;
    }
    CUDA_HOST_DEVICE
    iterator operator++(int) {
        iterator tmp(*this);
        operator++();
        return tmp;
    }
    CUDA_HOST_DEVICE
    size_t operator-(iterator start) const {
        size_t count = 0;
        while (start != *this) {
            start++;
            count++;
        }
        return count;
    }
    CUDA_HOST_DEVICE
    inline bool operator==(const iterator& rhs) const {
        return equal(rhs);
    }
    CUDA_HOST_DEVICE
    inline bool operator!=(const iterator& rhs) const {
        return !operator==(rhs);
    }

 private:
    friend class boost::iterator_core_access;

    CUDA_HOST_DEVICE
    void increment() {
#ifndef __CUDA_ARCH__
        LOG(4,"\tincrement (index_vector_iterator):"); 
#endif
        ++m_current_index;
#ifndef __CUDA_ARCH__
        LOG(4,"\tend increment (index_vector_iterator): m_current_index = "<<m_current_index); 
#endif
    }

    CUDA_HOST_DEVICE
    bool equal(iterator const& other) const {
        return m_current_index == other.m_current_index;
    }


    CUDA_HOST_DEVICE
    reference dereference() const
    { return *(m_particles_begin + *m_current_index); }


    Iterator m_current_index;
    p_pointer m_particles_begin;
};

template <typename Query>
class depth_first_iterator {
    typedef depth_first_iterator<Query> iterator;
    typedef typename Query::child_iterator child_iterator;
    static const unsigned int dimension = Query::dimension;
    typedef Vector<double,dimension> double_d;
    typedef Vector<int,dimension> int_d;
    typedef detail::bbox<dimension> box_type;

public:
    typedef typename child_iterator::value_type value_type;
    typedef child_iterator pointer;
	typedef std::forward_iterator_tag iterator_category;
    typedef typename child_iterator::reference reference;
	typedef std::ptrdiff_t difference_type;

    CUDA_HOST_DEVICE
    depth_first_iterator()
    {}
       
    /// this constructor is used to start the iterator at the head of a bucket 
    /// list
    CUDA_HOST_DEVICE
    depth_first_iterator(const child_iterator& start_node,
                        const Query *query
                  ):
        m_query(query)
    {
        if (start_node != false) {
            m_stack.push(start_node); 
        } else {
#ifndef __CUDA_ARCH__
            LOG(3,"\tdepth_first_iterator (constructor): start is false, no children to search.");
#endif
        }
    }


    box_type get_bounds() const {
        return m_stack.top().get_bounds();
    }

    const child_iterator& get_child_iterator() const {
        return m_stack.top();
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
    iterator& operator++() {
        increment();
        return *this;
    }
    CUDA_HOST_DEVICE
    iterator operator++(int) {
        iterator tmp(*this);
        operator++();
        return tmp;
    }
    CUDA_HOST_DEVICE
    size_t operator-(iterator start) const {
        size_t count = 0;
        while (start != *this) {
            start++;
            count++;
        }
        return count;
    }
    CUDA_HOST_DEVICE
    inline bool operator==(const iterator& rhs) const {
        return equal(rhs);
    }
    CUDA_HOST_DEVICE
    inline bool operator!=(const iterator& rhs) const {
        return !operator==(rhs);
    }

 private:
    friend class boost::iterator_core_access;

    CUDA_HOST_DEVICE
    void increment() {
#ifndef __CUDA_ARCH__
        LOG(4,"\tincrement (depth_first_iterator): depth = "<<m_stack.size()<<" top child number "<<m_stack.top().get_child_number()); 
#endif
        if (m_query->is_leaf_node(*m_stack.top())) {
            ++m_stack.top();
            while (!m_stack.empty() && m_stack.top() == false) {
                LOG(4,"\tpop stack (depth_first_iterator):"); 
                m_stack.pop();
            }
        } else {
            LOG(4,"\tpush stack (depth_first_iterator):"); 
            m_stack.push(m_query->get_children(m_stack.top()++));
        }

#ifndef __CUDA_ARCH__
        LOG(4,"\tend increment (depth_first_iterator):"); 
#endif
    }

    CUDA_HOST_DEVICE
    bool equal(iterator const& other) const {
        if (m_stack.empty() || other.m_stack.empty()) {
            return m_stack.empty() == other.m_stack.empty();
        } else {
            return m_stack.top() == other.m_stack.top();
        }
    }

    CUDA_HOST_DEVICE
    reference dereference() const
    { return *m_stack.top(); }


    std::stack<child_iterator> m_stack;
    const Query *m_query;
};

template <typename Query, int LNormNumber>
class tree_query_iterator {
    typedef tree_query_iterator<Query,LNormNumber> iterator;
    typedef typename Query::child_iterator child_iterator;
    static const unsigned int dimension = Query::dimension;
    typedef Vector<double,dimension> double_d;
    typedef Vector<int,dimension> int_d;
    typedef detail::bbox<dimension> box_type;

public:
    typedef typename child_iterator::value_type value_type;
    typedef typename child_iterator::pointer pointer;
	typedef std::forward_iterator_tag iterator_category;
    typedef typename child_iterator::reference reference;
	typedef std::ptrdiff_t difference_type;

    CUDA_HOST_DEVICE
    tree_query_iterator()
    {}
       
    /// this constructor is used to start the iterator at the head of a bucket 
    /// list
    CUDA_HOST_DEVICE
    tree_query_iterator(const child_iterator& start,
                  const double_d& query_point,
                  const double_d& max_distance,
                  const Query *query,
                  const bool ordered=false
                  ):
        m_query_point(query_point),
        m_inv_max_distance(1.0/max_distance),
        m_query(query)
    {
        if (start != false) {
            m_stack.push(start);
            go_to_next_leaf();
        } else {
#ifndef __CUDA_ARCH__
            LOG(3,"\ttree_query_iterator (constructor) with query pt = "<<m_query_point<<"): start is false, no children to search.");
#endif
        }

        if (m_stack.empty()) {
#ifndef __CUDA_ARCH__
            LOG(3,"\ttree_query_iterator (constructor) with query pt = "<<m_query_point<<"): search region outside domain or no children to search.");
#endif
        } else {
#ifndef __CUDA_ARCH__
            LOG(3,"\tocttree_query_iterator (constructor) with query pt = "<<m_query_point<<"):  found bbox = "<<m_query->get_bounds(m_stack.top()));
#endif
        }
    }

    tree_query_iterator(const iterator& copy):
        m_query_point(copy.m_query_point),
        m_inv_max_distance(copy.m_inv_max_distance),
        m_query(copy.m_query)
    {
        m_stack = copy.m_stack;
        //std::copy(copy.m_stack.begin(),copy.m_stack.end(),m_stack.begin()); 
    }

    box_type get_bounds() const {
        return m_stack.top().get_bounds();
    }

    iterator& operator=(const iterator& copy) {
        m_query_point = copy.m_query_point;
        m_inv_max_distance = copy.m_inv_max_distance;
        m_stack = copy.m_stack;
        m_query = copy.m_query;
        return *this;
        //std::copy(copy.m_stack.begin(),copy.m_stack.end(),m_stack.begin()); 
    }

    /*
    iterator& operator=(const octtree_depth_first_iterator<Query>& copy) {
        m_stack = copy.m_stack;
        go_to_next_leaf();
        return *this;
    }
    */


    CUDA_HOST_DEVICE
    reference operator *() const {
        return dereference();
    }
    CUDA_HOST_DEVICE
    reference operator ->() {
        return dereference();
    }
    CUDA_HOST_DEVICE
    iterator& operator++() {
        increment();
        return *this;
    }
    CUDA_HOST_DEVICE
    iterator operator++(int) {
        iterator tmp(*this);
        operator++();
        return tmp;
    }
    CUDA_HOST_DEVICE
    size_t operator-(iterator start) const {
        size_t count = 0;
        while (start != *this) {
            start++;
            count++;
        }
        return count;
    }
    CUDA_HOST_DEVICE
    inline bool operator==(const iterator& rhs) const {
        return equal(rhs);
    }
    CUDA_HOST_DEVICE
    inline bool operator!=(const iterator& rhs) const {
        return !operator==(rhs);
    }

 private:
    friend class boost::iterator_core_access;

    bool child_is_within_query(const child_iterator& node) {
        const box_type& bounds = m_query->get_bounds(node);
        double accum = 0;
        for (int j = 0; j < dimension; j++) {
            const bool less_than_bmin = m_query_point[j]<bounds.bmin[j];
            const bool more_than_bmax = m_query_point[j]>bounds.bmax[j];

            // dist 0 if between min/max, or distance to min/max if not
            const double dist = (less_than_bmin^more_than_bmax)
                                *(less_than_bmin?
                                        (bounds.bmin[j]-m_query_point[j]):
                                        (m_query_point[j]-bounds.bmax[j]));

            accum = detail::distance_helper<LNormNumber>::
                        accumulate_norm(accum,dist*m_inv_max_distance[j]); 
        }
        return (accum < 1.0);
    }

    void increment_stack() {
        while (!m_stack.empty()) {
            ++m_stack.top();
            LOG(4,"\tincrement stack with child "<<m_stack.top().get_child_number());
            if (m_stack.top() == false) {
                LOG(4,"\tincrement_stack: pop");
                m_stack.pop();
            } else {
                break;
            }
        } 
    }


    void go_to_next_leaf() {
        bool exit = m_stack.empty();
        while (!exit) {
            child_iterator& node = m_stack.top();
            LOG(3,"\tgo_to_next_leaf with child "<<node.get_child_number()<<" with bounds "<<node.get_bounds());
            if (child_is_within_query(node)) { // could be in this child
                LOG(4,"\tthis child is within query, so going to next child");
                if (m_query->is_leaf_node(*node)) {
                    exit = true;
                } else {
                    LOG(4,"\tdive down");
                    m_stack.push(m_query->get_children(node));
                }
            } else { // not in this one, so go to next child, or go up if no more children
                LOG(4,"\tthis child is NOT within query, so going to next child");
                increment_stack();
                exit = m_stack.empty();
            }
        }
#ifndef __CUDA_ARCH__
        if (4 <= ABORIA_LOG_LEVEL) { 
            if (m_stack.empty()) {
                LOG(4,"\tgo_to_next_leaf: stack empty, finishing");
            } else {
                LOG(4,"\tgo_to_next_leaf: found leaf, finishing");
            }
        }
#endif

    }

    void increment() {
#ifndef __CUDA_ARCH__
        LOG(4,"\tincrement (octtree_query_iterator):"); 
#endif
        increment_stack();
        go_to_next_leaf();

        if (m_stack.empty()) {
#ifndef __CUDA_ARCH__
            LOG(3,"\tend increment (octree_query_iterator): no more nodes"); 
#endif
        } else {
#ifndef __CUDA_ARCH__
            LOG(3,"\tend increment (octree_query_iterator): looking in bbox "<<m_query->get_bounds(m_stack.top())); 
#endif
        }
    }
    

    CUDA_HOST_DEVICE
    bool equal(iterator const& other) const {
        if (m_stack.empty() || other.m_stack.empty()) {
            return m_stack.empty() == other.m_stack.empty();
        } else {
            return m_stack.top() == other.m_stack.top();
        }
    }


    CUDA_HOST_DEVICE
    reference dereference() const
    { return *m_stack.top(); }


    std::stack<child_iterator> m_stack;
    double_d m_query_point;
    double_d m_inv_max_distance;
    const Query *m_query;
};




/*
template <typename Query, int LNormNumber>
class tree_query_iterator {
    typedef tree_query_iterator<Query,LNormNumber> iterator;
    static const unsigned int dimension = Query::dimension;
    typedef Vector<double,dimension> double_d;
    typedef Vector<int,dimension> int_d;

public:
    typedef typename Query::value_type const value_type;
    typedef const value_type* pointer;
	typedef std::forward_iterator_tag iterator_category;
    typedef const value_type& reference;
	typedef std::ptrdiff_t difference_type;

    CUDA_HOST_DEVICE
    tree_query_iterator():
        m_node(nullptr)
    {}
       
    /// this constructor is used to start the iterator at the head of a bucket 
    /// list
    CUDA_HOST_DEVICE
    tree_query_iterator(const value_type* start_node,
                  const double_d& query_point,
                  const double max_distance,
                  const Query *query
                  ):

        m_query_point(query_point),
        m_max_distance2(
                detail::distance_helper<LNormNumber>
                        ::get_value_to_accumulate(max_distance)),
        m_dists(0),
        m_query(query),
        m_node(nullptr)
    {
        //m_stack.reserve(m_query->get_max_levels());
        if (start_node == nullptr) {
                LOG(4,"\ttree_query_iterator (constructor) empty tree, returning default iterator");
        } else {
            double accum = 0;
            for (int i = 0; i < dimension; ++i) {
                const double val = m_query_point[i];
                if (val < m_query->get_bounds_low()[i]) {
                    m_dists[i] = val - m_query->get_bounds_low()[i];
                } else if (m_query_point[i] > m_query->get_bounds_high()[i]) {
                    m_dists[i] = val - m_query->get_bounds_high()[i];
                }
                accum = detail::distance_helper<LNormNumber>::accumulate_norm(accum,m_dists[i]); 
            }
            if (accum <= m_max_distance2) {
                LOG(4,"\ttree_query_iterator (constructor) with query pt = "<<m_query_point<<"): searching root node");
                m_node = start_node;
                go_to_next_leaf();
            } else {
                LOG(4,"\ttree_query_iterator (constructor) with query pt = "<<m_query_point<<"): search region outside domain");
            }
        }
    }


    tree_query_iterator(const iterator& copy):
        m_query_point(copy.m_query_point),
        m_max_distance2(copy.m_max_distance2),
        m_node(copy.m_node),
        m_dists(copy.m_dists)
    {
        m_stack = copy.m_stack;
        //std::copy(copy.m_stack.begin(),copy.m_stack.end(),m_stack.begin()); 
    }

    iterator& operator=(const iterator& copy) {
        m_query_point = copy.m_query_point;
        m_max_distance2 = copy.m_max_distance2;
        m_node=copy.m_node;
        m_dists=copy.m_dists;
        m_stack = copy.m_stack;
        return *this;
        //std::copy(copy.m_stack.begin(),copy.m_stack.end(),m_stack.begin()); 
    }

    iterator& operator=(const tree_depth_first_iterator<Query>& copy) {
        m_node=copy.m_node;
#ifndef NDEBUG
        const double_d low = copy.m_query->get_bounds_low(*m_node);
        const double_d high = copy.m_query->get_bounds_high(*m_node);
        ASSERT((low <= m_query_point).all() && (high > m_query_point).all(),"query point not in depth_first_iterator")
#endif
        std::copy(copy.m_stack.begin(),copy.m_stack.end(),m_stack.begin()); 
        return *this;
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
    iterator& operator++() {
        increment();
        return *this;
    }
    CUDA_HOST_DEVICE
    iterator operator++(int) {
        iterator tmp(*this);
        operator++();
        return tmp;
    }
    CUDA_HOST_DEVICE
    size_t operator-(iterator start) const {
        size_t count = 0;
        while (start != *this) {
            start++;
            count++;
        }
        return count;
    }
    CUDA_HOST_DEVICE
    inline bool operator==(const iterator& rhs) const {
        return equal(rhs);
    }
    CUDA_HOST_DEVICE
    inline bool operator!=(const iterator& rhs) const {
        return !operator==(rhs);
    }

 private:
    friend class boost::iterator_core_access;

    void go_to_next_leaf() {
        while(!m_query->is_leaf_node(*m_node)) {
            ASSERT(m_query->get_child1(m_node) != NULL,"no child1");
            ASSERT(m_query->get_child2(m_node) != NULL,"no child2");
            // Which child branch should be taken first?
            const size_t idx = m_query->get_dimension_index(*m_node);
            const double val = m_query_point[idx];
            const double diff_cut_high = val - m_query->get_cut_high(*m_node);
            const double diff_cut_low = val - m_query->get_cut_low(*m_node);

            pointer bestChild;
            pointer otherChild;
            double cut_dist,bound_dist;


            LOG(4,"\ttree_query_iterator (go_to_next_leaf) with query pt = "<<m_query_point<<"): idx = "<<idx<<" cut_high = "<<m_query->get_cut_high(*m_node)<<" cut_low = "<<m_query->get_cut_low(*m_node));
            
            if ((diff_cut_low+diff_cut_high)<0) {
                LOG(4,"\ttree_query_iterator (go_to_next_leaf) low child is best");
                bestChild = m_query->get_child1(m_node);
                otherChild = m_query->get_child2(m_node);
                cut_dist = diff_cut_high;
                //bound_dist = val - m_query->get_bounds_high(*m_node);
            } else {
                LOG(4,"\ttree_query_iterator (go_to_next_leaf) high child is best");
                bestChild = m_query->get_child2(m_node);
                otherChild = m_query->get_child1(m_node);
                cut_dist = diff_cut_low;
                //bound_dist = val - m_query->get_bounds_low(*m_node);
            }
            
            // if other child possible save it to stack for later
            double save_dist = m_dists[idx];
            m_dists[idx] = cut_dist;

            // calculate norm of m_dists 
            double accum = 0;
            for (int i = 0; i < dimension; ++i) {
                accum = detail::distance_helper<LNormNumber>::accumulate_norm(accum,m_dists[i]); 
            }
            if (accum <= m_max_distance2) {
                LOG(4,"\ttree_query_iterator (go_to_next_leaf) push other child for later");
                m_stack.push(std::make_pair(otherChild,m_dists));
            }

            // restore m_dists and move to bestChild
            m_dists[idx] = save_dist;
            m_node = bestChild;
        }
        LOG(4,"\ttree_query_iterator (go_to_next_leaf) found a leaf node");
    }
    
    void pop_new_child_from_stack() {
        std::tie(m_node,m_dists) = m_stack.top();
        m_stack.pop();
    }

    CUDA_HOST_DEVICE
    void increment() {
#ifndef __CUDA_ARCH__
        LOG(4,"\tincrement (tree_iterator):"); 
#endif
        if (m_stack.empty()) {
            m_node = nullptr;
        } else {
            pop_new_child_from_stack();
            go_to_next_leaf();
        }

#ifndef __CUDA_ARCH__
        LOG(4,"\tend increment (tree_iterator): m_node = "<<m_node); 
#endif
    }

    CUDA_HOST_DEVICE
    bool equal(iterator const& other) const {
        return m_node == other.m_node;
    }


    CUDA_HOST_DEVICE
    reference dereference() const
    { return *m_node; }


    std::stack<std::pair<pointer,double_d>> m_stack;
    double_d m_query_point;
    double m_max_distance2;
    const value_type* m_node;
    double_d m_dists;
    const Query *m_query;
};
*/



/*
template <typename Query>
class tree_theta_iterator {
    typedef tree_theta_iterator<Query> iterator;
    static const unsigned int dimension = Query::dimension;
    typedef Vector<double,dimension> double_d;
    typedef Vector<int,dimension> int_d;
    typedef Query::reference node_reference;
    typedef Query::pointer node_pointer;
    typedef Query::value_type node_value_type;

public:
    typedef const node_pointer pointer;
	typedef std::forward_iterator_tag iterator_category;
    typedef const node_reference reference;
    typedef const node_value_type value_type;
	typedef std::ptrdiff_t difference_type;

    CUDA_HOST_DEVICE
    tree_theta_iterator():
        m_node(nullptr)
    {}

    tree_theta_iterator(const node_pointer node,
                        const Query *query,
                        const std::stack<const node_pointer>& stack):
        m_node(node),
        m_query(query),
        m_stack(stack),
        m_centre(0.5*(m_query->get_bounds_low(m_node) 
                      + m_query->get_bounds_high(m_node))),
        m_r2((m_query->get_bounds_high(m_node) - m_centre).squaredNorm()),
        m_r(std::sqrt(m_r2))
        m_theta2(std::pow(0.5,2))
    {}
       
    tree_theta_iterator(const tree_depth_first_iterator<Query> copy):
        tree_theta_iterator(copy.m_node,copy.m_query,copy.m_stack)
    {}


    tree_theta_iterator(const iterator& copy):
        m_query_point(copy.m_query_point),
        m_max_distance2(copy.m_max_distance2),
        m_node(copy.m_node),
        m_dists(copy.m_dists)
        m_stack(copy.m_stack)
    {}

    iterator& operator=(const iterator& copy) {
        m_query_point = copy.m_query_point;
        m_max_distance2 = copy.m_max_distance2;
        m_node=copy.m_node;
        m_dists=copy.m_dists;
        m_stack = copy.m_stack;
        return *this;
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
    iterator& operator++() {
        increment();
        return *this;
    }
    CUDA_HOST_DEVICE
    iterator operator++(int) {
        iterator tmp(*this);
        operator++();
        return tmp;
    }
    CUDA_HOST_DEVICE
    size_t operator-(iterator start) const {
        size_t count = 0;
        while (start != *this) {
            start++;
            count++;
        }
        return count;
    }
    CUDA_HOST_DEVICE
    inline bool operator==(const iterator& rhs) const {
        return equal(rhs);
    }
    CUDA_HOST_DEVICE
    inline bool operator!=(const iterator& rhs) const {
        return !operator==(rhs);
    }

 private:

    bool theta_condition(const node_reference& node) {
        double_d other_size = m_query->get_bounds_low(node)-m_query->get_bounds_low(node);
        double d = 0.5*(other_high + other_low - m_low - m_high).norm(); 
        double other_r2 = 0.25*(other_high-other_low).squaredNorm();
        if (other_r2 < m_r2) {
            const double other_r = std::sqrt(other_r2);
            return m_r2 > m_theta2*std::pow(d-other_r,2);
        } else {
            return other_r2 < m_theta2*std::pow(d-m_r,2);
        }
    }

    void go_to_next_node() {
        m_theta_condition = true;
        while(!m_query->is_leaf_node(*m_node) && m_theta_condition) {
            ASSERT(m_query->get_child1(m_node) != NULL,"no child1");
            ASSERT(m_query->get_child2(m_node) != NULL,"no child2");
            node_pointer child1 = m_query->get_child1(m_node);
            node_pointer child2 = m_query->get_child1(m_node);
            if (theta_condition(child1)) {
                m_stack.push(child2);
                m_node = child1;
            } else if (theta_condition(child2)) {
                m_node = child2;
            } else {
                pop_new_child_from_stack();
            }

            bool child1_theta = theta_condition(child1);
            bool child2_theta = theta_condition(child2);
            if (child1_theta && child2_theta) {
                m_node = child1;
                //keep going
            } else if (child1_theta) {
                m_stack.push(child1);
                m_node = child2;
                m_theta_condition = false;
                //return
            } else if (child2_theta) {
                m_stack.push(child2);
                m_node = child1;
                m_theta_condition = false;
                //return
            } else {
                CHECK(false,"should not get here!");
            }

        }
        LOG(4,"\ttree_interaction_iterator (go_to_next_leaf) found a candidate node. m_theta_condition = "<<m_theta_condition);
    }
    
    // returns true if the new node satisfies the theta condition 
    void pop_new_child_from_stack() {
        m_node = m_stack.top();
        m_stack.pop();
    }

    CUDA_HOST_DEVICE
    void increment() {
#ifndef __CUDA_ARCH__
        LOG(4,"\tincrement (tree_interaction_iterator):"); 
#endif
        if (m_stack.empty()) {
            m_node = nullptr;
        } else {
            do {
                pop_new_child_from_stack();
            } while (!theta_condition(m_node));
            go_to_next_leaf();
        }

#ifndef __CUDA_ARCH__
        LOG(4,"\tend increment (tree_interaction_iterator): m_node = "<<m_node); 
#endif
    }

    CUDA_HOST_DEVICE
    bool equal(iterator const& other) const {
        return m_node == other.m_node;
    }


    CUDA_HOST_DEVICE
    reference dereference() const
    { return reference(*m_node,m_theta); }


    std::stack<pointer> m_stack;
    double_d m_low;
    double_d m_high;
    double m_r2;
    double m_r;
    bool m_theta_condition;
    double m_theta2;
    const value_type* m_node;
    const Query *m_query;
};
*/


/*
// assume that these iterators, and query functions, can be called from device code
template <unsigned int D,unsigned int LNormNumber>
class lattice_query_iterator {
    typedef lattice_query_iterator<D,LNormNumber> iterator;
    typedef Vector<double,D> double_d;
    typedef Vector<int,D> int_d;

    int_d m_min;
    int_d m_max;
    int_d m_index;
    double_d m_index_point;
    double_d m_query_point;
    int_d m_query_index;
    double m_max_distance2;
    double_d m_box_size;
    bool m_check_distance;
public:
    typedef const int_d* pointer;
	typedef std::forward_iterator_tag iterator_category;
    typedef const int_d& reference;
    typedef const int_d value_type;
	typedef std::ptrdiff_t difference_type;

    CUDA_HOST_DEVICE
    lattice_query_iterator(
                     const int_d& min,
                     const int_d& max,
                     const int_d& start,
                     const double_d& start_point,
                     const double_d& query_point,
                     const int_d& query_index,
                     const double max_distance,
                     const double_d& box_size):
        m_min(min),
        m_max(max),
        m_index(start),
        m_index_point(start_point),
        m_query_point(query_point),
        m_query_index(query_index),
        m_max_distance2(
                detail::distance_helper<LNormNumber>
                        ::get_value_to_accumulate(max_distance)),
        m_box_size(box_size),
        m_check_distance(true)
    {
        if (distance_helper<LNormNumber>::norm(m_index_point)<m_max_distance2) {
            increment();
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
    iterator& operator++() {
        increment();
        return *this;
    }

    CUDA_HOST_DEVICE
    iterator operator++(int) {
        iterator tmp(*this);
        operator++();
        return tmp;
    }


    CUDA_HOST_DEVICE
    size_t operator-(iterator start) const {
        size_t count = 0;
        while (start != *this) {
            start++;
            count++;
        }
        return count;
    }

    CUDA_HOST_DEVICE
    inline bool operator==(const iterator& rhs) const {
        return equal(rhs);
    }

    CUDA_HOST_DEVICE
    inline bool operator!=(const iterator& rhs) const {
        return !operator==(rhs);
    }

private:
    friend class boost::iterator_core_access;

    CUDA_HOST_DEVICE
    bool equal(iterator const& other) const {
        return (m_index == other.m_index).all();
    }

    CUDA_HOST_DEVICE
    reference dereference() const { 
        return m_index; 
    }

    CUDA_HOST_DEVICE
    void increment() {
        if (m_check_distance) {
            double distance;
            do  {
                for (int i=0; i<D; i++) {
                    ++m_index[i];
                    if (m_index[i] == m_query_index[i]) {
                        m_check_distance = false;
                    } else {
                        m_index_dist[i] += m_box_size[i];
                    }
                    if (m_index[i] <= m_max[i]) break;
                    if (i != D-1) {
                        m_index[i] = m_min[i];
                    } 
                }
                if (m_index[D-1] > m_max[D-1]) {
                    distance = 0;
                } else {
                    distance = distance_helper<LNormNumber>::norm(m_index_dist);
                }
            } while (distance > m_max_distance2);
        } else {
            for (int i=0; i<D; i++) {
                if (m_index[i] == m_query_index[i]) {
                    m_check_distance = true;
                }
                ++m_index[i];
                if (m_index[i] <= m_max[i]) break;
                if (i != D-1) {
                    m_index[i] = m_min[i];
                } 
            }
        }
    }
};

*/

template <unsigned int D>
class lattice_iterator {
    typedef lattice_iterator<D> iterator;
    typedef Vector<double,D> double_d;
    typedef Vector<int,D> int_d;

    // make a proxy int_d in case you ever
    // want to get a pointer object to the
    // reference (which are both of the
    // same type)
    struct proxy_int_d: public int_d {
        proxy_int_d():
            int_d() 
        {}

        proxy_int_d(const int_d& arg):
            int_d(arg) 
        {}

        proxy_int_d& operator&() {
            return *this;
        }
        
        const proxy_int_d& operator&() const {
            return *this;
        }

        const proxy_int_d& operator*() const {
            return *this;
        }

        proxy_int_d& operator*() {
            return *this;
        }

        const proxy_int_d* operator->() const {
            return this;
        }

       proxy_int_d* operator->() {
            return this;
        }
    };

   
    int_d m_min;
    int_d m_max;
    proxy_int_d m_index;
    int_d m_size;
    bool m_valid;
public:
    typedef proxy_int_d pointer;
	typedef std::random_access_iterator_tag iterator_category;
    typedef const proxy_int_d& reference;
    typedef proxy_int_d value_type;
	typedef std::ptrdiff_t difference_type;

    lattice_iterator():
        m_valid(false)
    {}

    lattice_iterator(const int_d &min, 
                     const int_d &max):
        m_min(min),
        m_max(max),
        m_index(min),
        m_size(minus(max,min)),
        m_valid(true)
    {}

    explicit operator size_t() const {
        return collapse_index_vector(m_index);
    }

    const lattice_iterator& get_child_iterator() const {
        return *this;
    }

    reference operator *() const {
        return dereference();
    }

    reference operator ->() const {
        return dereference();
    }

    iterator& operator++() {
        increment();
        return *this;
    }

    iterator operator++(int) {
        iterator tmp(*this);
        operator++();
        return tmp;
    }

    iterator operator+(const int n) {
        iterator tmp(*this);
        tmp.increment(n);
        return tmp;
    }

    iterator& operator+=(const int n) {
        increment(n);
        return *this;
    }

    iterator& operator-=(const int n) {
        increment(-n);
        return *this;
    }

    iterator operator-(const int n) {
        iterator tmp(*this);
        tmp.increment(-n);
        return tmp;
    }

    size_t operator-(const iterator& start) const {
        int distance;
        if (!m_valid) {
            const int_d test = minus(minus(start.m_max,1),start.m_index);
            distance = start.collapse_index_vector(minus(minus(start.m_max,1),start.m_index))+1;
        } else if (!start.m_valid) {
            distance = collapse_index_vector(minus(m_index,m_min));
        } else {
            distance = collapse_index_vector(minus(m_index,start.m_index));
        }
        assert(distance > 0);
        return distance;
    }

    inline bool operator==(const iterator& rhs) const {
        return equal(rhs);
    }

    inline bool operator==(const bool rhs) const {
        return equal(rhs);
    }

    inline bool operator!=(const iterator& rhs) const {
        return !operator==(rhs);
    }

    inline bool operator!=(const bool rhs) const {
        return !operator==(rhs);
    }

private:

    static inline 
    int_d minus(const int_d& arg1, const int_d& arg2) {
        int_d ret;
        for (size_t i = 0; i < D; ++i) {
            ret[i] = arg1[i]-arg2[i];
        }
        return ret;
    }

    static inline 
    int_d minus(const int_d& arg1, const int arg2) {
        int_d ret;
        for (size_t i = 0; i < D; ++i) {
            ret[i] = arg1[i]-arg2;
        }
        return ret;
    }

    int collapse_index_vector(const int_d &vindex) const {
        int index = 0;
        unsigned int multiplier = 1.0;
        for (int i = D-1; i>=0; --i) {
            if (i != D-1) {
                multiplier *= m_size[i+1];
            }
            index += multiplier*vindex[i];
        }
        return index;
    }

    int_d reassemble_index_vector(const int index) const {
        int_d vindex;
        int i = index;
        for (int d = D-1; d>=0; --d) {
            double div = (double)i / m_size[d];
            vindex[d] = std::round((div-std::floor(div)) * m_size[d]);
            i = std::floor(div);
        }
        return vindex;
    }

    bool equal(iterator const& other) const {
        if (!other.m_valid) return !m_valid;
        if (!m_valid) return !other.m_valid;
        for (size_t i = 0; i < D; ++i) {
           if (m_index[i] != other.m_index[i]) {
               return false;
           }
        }
        return true;
    }

    bool equal(const bool other) const {
        return m_valid==other;
    }

    reference dereference() const { 
        return m_index; 
    }

    void increment() {
        for (int i=D-1; i>=0; --i) {
            ++m_index[i];
            if (m_index[i] < m_max[i]) break;
            if (i != 0) {
                m_index[i] = m_min[i];
            } else {
                m_valid = false;
            }
        }
    }

    void increment(const int n) {
        int collapsed_index = collapse_index_vector(m_index);
        m_index = reassemble_index_vector(collapsed_index += n);
    }
};




}


#endif /* BUCKETSEARCH_H_ */
