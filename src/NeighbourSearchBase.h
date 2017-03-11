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
    iterator_range(IteratorType&& begin, IteratorType&& end):
        m_begin(std::move(begin)),m_end(std::move(end)) 
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

    static constexpr bool unordered() {
        return true;
    }

    /// resets the domain extents, periodicity and bucket size
    /// \param low the lower extent of the search domain
    /// \param high the upper extent of the search domain
    /// \param _max_interaction_radius the side length of each bucket
    /// \param periodic a boolean vector indicating wether each dimension
    void set_domain(const double_d &min_in, const double_d &max_in, const bool_d& periodic_in, const unsigned int n_particles_in_leaf=10, const bool not_in_constructor=true) {
        LOG_CUDA(2,"neighbour_search_base: set_domain:");
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
	typedef std::forward_iterator_tag iterator_category;
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
        size_t count = 0;
        while (start != *this) {
            ++start; ++count;
        }
        return count;
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
    bool increment() {
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

template <unsigned int D, typename Query, int LNormNumber>
class tree_query_iterator {
    typedef tree_query_iterator<D,Query,LNormNumber> iterator;
    static const unsigned int dimension = D;
    typedef Vector<double,D> double_d;
    typedef Vector<int,D> int_d;

public:
    typedef const typename Query::pointer pointer;
	typedef std::forward_iterator_tag iterator_category;
    typedef const typename Query::reference reference;
    typedef const typename Query::value_type value_type;
	typedef std::ptrdiff_t difference_type;

    CUDA_HOST_DEVICE
    tree_query_iterator():
        m_node(nullptr)
    {}
       
    /// this constructor is used to start the iterator at the head of a bucket 
    /// list
    CUDA_HOST_DEVICE
    tree_query_iterator(pointer start_node,
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
        m_node(start_node)
    {
        m_stack.reserve(m_query->get_max_levels());
        go_to_next_leaf();
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

    /*
    iterator& operator=(iterator&& copy):
        m_stack=std::move(copy.m_stack),
        m_query_point(copy.m_query_point),
        m_max_distance2(copy.m_max_distance2),
        m_node(copy.m_node),
        m_dists(copy.m_dists)
    {}
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

    void go_to_next_leaf() {
        while(!m_query->is_leaf_node(*m_node)) {
            /* Which child branch should be taken first? */
            const size_t idx = m_query->get_dimension_index(*m_node);
            double val = m_query_point[idx];
            double diff1 = val - m_query->get_cut_low(*m_node);
            double diff2 = val - m_query->get_cut_high(*m_node);

            pointer bestChild;
            pointer otherChild;
            double cut_dist;
            if ((diff1+diff2)<0) {
                bestChild = m_query->get_child1(m_node);
                otherChild = m_query->get_child2(m_node);
                cut_dist = std::abs(diff2);
            } else {
                bestChild = m_query->get_child2(m_node);
                otherChild = m_query->get_child1(m_node);
                cut_dist = std::abs(diff1);
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
                m_stack.push(std::make_pair(otherChild,m_dists));
            }

            // restore m_dists and move to bestChild
            m_dists[idx] = save_dist;
            m_node = bestChild;
        }
    }
    
    void pop_new_child_from_stack() {
        std::tie(m_node,m_dists) = m_stack.top();
        m_stack.pop();
    }

    CUDA_HOST_DEVICE
    bool increment() {
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


// assume that these iterators, and query functions, can be called from device code
template <unsigned int D>
class lattice_iterator {
    typedef Vector<double,D> double_d;
    typedef Vector<int,D> int_d;

    int_d m_min;
    int_d m_max;
    int_d m_index;
    detail::bucket_index<D> m_bucket_index;
public:
    typedef const int_d* pointer;
	typedef std::random_access_iterator_tag iterator_category;
    typedef const int_d& reference;
    typedef const int_d value_type;
	typedef std::ptrdiff_t difference_type;

    lattice_iterator()
    {}

    CUDA_HOST_DEVICE
    lattice_iterator(const int_d &min, 
                     const int_d &max, 
                     const int_d &start):
        m_min(min),
        m_max(max),
        m_index(start),
        m_bucket_index(max-min)
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
    lattice_iterator& operator++() {
        increment();
        return *this;
    }

    CUDA_HOST_DEVICE
    lattice_iterator operator++(int) {
        lattice_iterator tmp(*this);
        operator++();
        return tmp;
    }

    CUDA_HOST_DEVICE
    lattice_iterator operator+(const int n) {
        lattice_iterator tmp(*this);
        tmp.increment(n);
        return tmp;
    }

    CUDA_HOST_DEVICE
    lattice_iterator& operator+=(const int n) {
        increment(n);
        return *this;
    }

    CUDA_HOST_DEVICE
    lattice_iterator& operator-=(const int n) {
        increment(-n);
        return *this;
    }

    CUDA_HOST_DEVICE
    lattice_iterator operator-(const int n) {
        lattice_iterator tmp(*this);
        tmp.increment(-n);
        return tmp;
    }

    CUDA_HOST_DEVICE
    size_t operator-(lattice_iterator start) const {
        const int distance = m_bucket_index.collapse_index_vector(m_index-start.m_index);
        ASSERT(distance > 0, "start iterator not before this iterator!");
        return distance;
    }

    CUDA_HOST_DEVICE
    inline bool operator==(const lattice_iterator& rhs) const {
        return equal(rhs);
    }

    CUDA_HOST_DEVICE
    inline bool operator!=(const lattice_iterator& rhs) const {
        return !operator==(rhs);
    }

private:
    friend class boost::iterator_core_access;

    CUDA_HOST_DEVICE
    bool equal(lattice_iterator const& other) const {
        return (m_index == other.m_index).all();
    }

    CUDA_HOST_DEVICE
    reference dereference() const { 
        return m_index; 
    }

    CUDA_HOST_DEVICE
    void increment() {
        for (int i=0; i<D; i++) {
            ++m_index[i];
            if (m_index[i] <= m_max[i]) break;
            if (i != D-1) {
                m_index[i] = m_min[i];
            } 
        }
    }

    CUDA_HOST_DEVICE
    void increment(const int n) {
        int collapsed_index = m_bucket_index.collapse_index_vector(m_index);
        m_index = m_bucket_index.reassemble_index_vector(collapsed_index += n);
    }
};

/*

template <typename Traits>
class lattice_iterator_with_hole {
    typedef typename Traits::double_d double_d;
    typedef typename Traits::int_d int_d;

    int_d m_min;
    int_d m_max;
    int_d m_hole_min;
    int_d m_hole_max;
    int_d m_index;
    detail::bucket_index<Traits::dimension> m_bucket_index;
public:
    typedef const int_d* pointer;
	typedef std::forward_iterator_tag iterator_category;
    typedef const int_d& reference;
    typedef const int_d value_type;
	typedef std::ptrdiff_t difference_type;

    CUDA_HOST_DEVICE
    lattice_iterator_with_hole(const int_d &min, 
                     const int_d &max, 
                     const int_d &hmin, 
                     const int_d &hmax, 
                     const int_d &start):
        m_min(min),
        m_max(max),
        m_hole_min(hmin),
        m_hole_max(hmax),
        m_index(start),
        m_bucket_index(max-min)
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
    lattice_iterator_with_hole& operator++() {
        increment();
        return *this;
    }

    CUDA_HOST_DEVICE
    lattice_iterator_with_hole operator++(int) {
        lattice_iterator_with_hole tmp(*this);
        operator++();
        return tmp;
    }

    
    CUDA_HOST_DEVICE
    size_t operator-(lattice_iterator_with_hole start) const {
        size_t count = 0;
        while (start != *this) {
            start++;
            count++;
        }
        return count;
    }

    CUDA_HOST_DEVICE
    inline bool operator==(const lattice_iterator_with_hole& rhs) const {
        return equal(rhs);
    }


    CUDA_HOST_DEVICE
    inline bool operator!=(const lattice_iterator_with_hole& rhs) const {
        return !operator==(rhs);
    }

private:
    friend class boost::iterator_core_access;

    CUDA_HOST_DEVICE
    bool equal(lattice_iterator_with_hole const& other) const {
        return (m_index == other.m_index).all();
    }

    CUDA_HOST_DEVICE
    reference dereference() const { 
        return m_index; 
    }

    CUDA_HOST_DEVICE
    void increment() {
        for (int i=0; i<Traits::dimension; i++) {
            ++m_index[i];
            if (m_index[i] == m_hole_min[i]) {
                m_index[i] = m_hole_max[i]+1;
            }
            if (m_index[i] <= m_max[i]) break;
            if (i != Traits::dimension-1) {
                m_index[i] = m_min[i];
            }
        }
    }

   
};
*/






}


#endif /* BUCKETSEARCH_H_ */
