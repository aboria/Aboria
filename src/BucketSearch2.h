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

#include "Traits.h"
#include "CudaInclude.h"
#include "Vector.h"
#include "SpatialUtil.h"
#include "Get.h"

namespace Aboria {

template <typename traits>
class BucketSearch {
    UNPACK_TRAITS(traits)

public:

    /// A const iterator to a set of neighbouring points. This iterator implements
    /// a STL forward iterator type
	class const_iterator;

    BucketSearch() {};

    void embed_points(particles_iterator begin, particles_iterator end) {
        m_particles_begin = begin;
        m_particles_end = end;
        m_positions_begin = get<position>(m_particles_begin);
        m_positions_end = get<position>(m_particles_end);

        CHECK(!m_bounds.is_empty(), "trying to embed particles into an empty domain. use the function `set_domain` to setup the spatial domain first.");

        m_bucket_indices.resize(m_particles_begin-m_particles_end);

        build_bucket_indices(m_positions_begin,m_positions_end,m_bucket_indices.begin());
        sort_by_bucket_index();
        build_buckets();
    }


    void add_points_at_end(const particles_iterator &begin, const particles_iterator &start_adding, const particles_iterator &end);

    /// return a const forward iterator to all the points in the neighbourhood of \p r. If 
    /// this function is being used to find all the point pairs within the same point container, then
    /// a naive looping through and using find_broadphase_neighbours() will find each pair twice. 
    /// This can be avoided by setting self=true and supplying the index of each point with my_index
    const_iterator find_broadphase_neighbours(const double_d& r, 
                                              const int my_index, 
                                              const bool self) const;

    const_iterator end() { return const_iterator(this); }


    void set_domain(const double_d &min_in, const double_d &max_in, const bool_d& periodic_in, const double_d& side_length) {
        LOG(2,"BucketSearch: set_domain:");
        m_bounds.bmin = min_in;
        m_bounds.bmax = max_in;
        m_periodic = periodic_in;
        m_bucket_side_length = side_length;
        m_size = ((m_bounds.bmax-m_bounds.bmin)/m_bucket_side_length).template cast<unsigned int>();
        m_bucket_side_length = (m_bounds.bmax-m_bounds.bmin)/m_size;
	    LOG(2,"\tbounds = "<<m_bounds);
	    LOG(2,"\tperiodic = "<<m_periodic);
	    LOG(2,"\tbucket_side_length = "<<m_bucket_side_length);

        // setup bucket data structures
        m_bucket_begin.resize(m_size.prod());
        m_bucket_end.resize(m_size.prod());
    }


    const double_d& get_min() const { return m_bounds.bmin; }
    const double_d& get_max() const { return m_bounds.bmax; }
    const double_d& get_side_length() const { return m_bucket_side_length; }
    const bool_d& get_periodic() const { return m_periodic; }


private:
    void build_buckets();
    void build_bucket_indices(
        vector_double_d_const_iterator positions_begin,
        vector_double_d_const_iterator positions_end,
        vector_unsigned_int_iterator bucket_indices_begin);
    void sort_by_bucket_index();
 

	inline unsigned int collapse_index_vector(const unsigned_int_d &vindex) const {
        unsigned int index = vindex[0];
        unsigned int multiplier = 1.0;
        for (int i=1; i<dimension; i++) {
            multiplier *= m_size[i];
		    ASSERT((vindex[i] > 0) && (vindex[i] < m_size[i]), "index is outside of dimension "<<i<<": "<<vindex);
            index += multiplier*vindex[i];
        }
        return index;
    }


    // hash a point in the unit square to the index of
    // the grid bucket that contains it
	inline unsigned int find_bucket_index(const double_d &r) const {
        // find the raster indices of p's bucket
        unsigned int index = (r[0]-m_bounds.bmin[0])/m_bucket_side_length[0];
        unsigned int multiplier = 1.0;
        for (int i=1; i<dimension; i++) {
            multiplier *= m_size[i];
            const unsigned int raster_d = (r[i]-m_bounds.bmin[i])/m_bucket_side_length[i];
		    ASSERT((raster_d > 0) && (raster_d < m_size[i]), "position is outside of dimension "<<i<<": "<<r);
            index += multiplier*raster_d;
        }
        return index;
    }
     
    struct point_to_bucket_index {
        const BucketSearch& bs; 
        CUDA_HOST_DEVICE
        point_to_bucket_index(const BucketSearch& bs):bs(bs) {}

        CUDA_HOST_DEVICE
        unsigned int operator()(const double_d& v) const {
            return bs.find_bucket_index(v);
        }
    };

    particles_iterator m_particles_begin;
    particles_iterator m_particles_end;
    vector_double_d_const_iterator m_positions_begin;
    vector_double_d_const_iterator m_positions_end;
    bool_d m_periodic;
    double_d m_bucket_side_length; 
    unsigned_int_d m_size;
    bbox<dimension> m_bounds;

    // the grid data structure keeps a range per grid bucket:
    // each bucket_begin[i] indexes the first element of bucket i's list of points
    // each bucket_end[i] indexes one past the last element of bucket i's list of points
    vector_unsigned_int m_bucket_begin;
    vector_unsigned_int m_bucket_end;
    vector_unsigned_int m_bucket_indices;

    unsigned_int_d m_surrounding_buckets_offsets;
};


template <typename traits>
void BucketSearch<traits>::build_bucket_indices(
        vector_double_d_const_iterator positions_begin,
        vector_double_d_const_iterator positions_end,
        vector_unsigned_int_iterator bucket_indices_begin
        ) {
    // transform the points to their bucket indices
    transform(positions_begin,
            positions_end,
            bucket_indices_begin,
            point_to_bucket_index(*this));
}

template <typename traits>
void BucketSearch<traits>::sort_by_bucket_index() {
    // sort the points by their bucket index
    traits::sort_by_key(m_bucket_indices.begin(),
            m_bucket_indices.end(),
            m_particles_begin);
}


template <typename traits>
void BucketSearch<traits>::build_buckets() {

    // find the beginning of each bucket's list of points
    typename traits::template counting_iterator<unsigned int> search_begin(0);
    traits::lower_bound(m_bucket_indices.begin(),
            m_bucket_indices.end(),
            search_begin,
            search_begin + m_size.prod(),
            m_bucket_begin.begin());

    // find the end of each bucket's list of points
    traits::upper_bound(m_bucket_indices.begin(),
            m_bucket_indices.end(),
            search_begin,
            search_begin + m_size.prod(),
            m_bucket_end.begin());
}


template <typename traits>
void BucketSearch<traits>::add_points_at_end(const particles_iterator &begin, const particles_iterator &start_adding, const particles_iterator &end) {
    m_particles_begin = begin;
    m_particles_end = end;
    m_positions_begin = get<position>(m_particles_begin);
    m_positions_end = get<position>(m_particles_end);

    CHECK(!m_bounds.is_empty(), "trying to embed particles into an empty domain. use the function `set_domain` to setup the spatial domain first.");

    const size_t dist = start_adding-end;
    vector_double_d_const_iterator positions_start_adding = m_positions_end - dist;
    vector_unsigned_int_iterator bucket_indices_start_adding = m_bucket_indices.end() - dist;
    build_bucket_indices(positions_start_adding,m_positions_end,bucket_indices_start_adding);
    sort_by_bucket_index();
    build_buckets();
}


template <typename traits>
typename BucketSearch<traits>::const_iterator
BucketSearch<traits>::find_broadphase_neighbours(
        const double_d& r, 
        const int my_index, 
        const bool self) const {
    
    const unsigned int my_bucket = find_bucket_index(r);

    const_iterator search_iterator(this,r);
    int_d bucket_offset(-1);
    constexpr unsigned int last_d = dimension-1;
    bool still_going = true;
    while (still_going) {
        unsigned_int_d other_bucket = my_bucket + bucket_offset; 

        // handle end cases
        double_d transpose(0);
        for (int i=0; i<dimension; i++) {
            if (other_bucket[i] == std::numeric_limits<unsigned int>::max()) {
                if (m_periodic[i]) {
                    other_bucket[i] = m_size[i]-1;
                    transpose[i] = -m_size[i];
                } else {
                    break;
                }
            }
            if (other_bucket[i] == m_size[i]) {
                if (m_periodic[i]) {
                    other_bucket[i] = 0;
                    transpose[i] = m_size[i];
                } else {
                    break;
                }
            }
        }

        const unsigned int other_bucket_index = collapse_index_vector(other_bucket);
        search_iterator.add_range(
                m_particles_begin + m_bucket_begin[other_bucket_index],
                m_particles_begin + m_bucket_end[other_bucket_index],
                transpose);

        // go to next candidate bucket
        for (int i=0; i<dimension; i++) {
            bucket_offset[i]++;
            if (bucket_offset[i] <= 1) break;
            if (i == last_d) still_going = false;
            bucket_offset[i] = 0;
        }
    }
    
    return search_iterator;
}



template <typename traits>
class BucketSearch<traits>::const_iterator {
public:
    typedef const std::tuple<const particles_value_type&,const double_d&>* pointer;
	typedef std::forward_iterator_tag iterator_category;
    typedef const std::tuple<particles_reference_type,const double_d&> reference;
    typedef const std::tuple<particles_reference_type,const double_d&> value_type;
	typedef std::ptrdiff_t difference_type;

    const_iterator(const BucketSearch<traits>* bucket_sort):
        m_bucket_sort(bucket_sort),
        m_node(bucket_sort->m_particles_end) {}

    const_iterator(const BucketSearch<traits>* bucket_sort, const double_d &r):
        m_bucket_sort(bucket_sort),
        m_r(r) {}

    void add_range(particles_iterator begin, particles_iterator end, const double_d &transpose) {
        m_begins.push_back(begin);
        m_ends.push_back(end);
        m_transpose.push_back(transpose);
        if (m_begins.size() == 1) {
            m_current_index = 0;
            m_node = m_begins[m_current_index];
        }
    }

    bool equal(const_iterator const& other) const {
        return m_node == other.m_node;
    }

    reference dereference() const { 
        return reference(*m_node,m_dx); 
    }

    bool go_to_next_candidate() {
        m_node++;
        if (m_node == m_ends[m_current_index]) {
            m_current_index++;
            if (m_current_index < m_begins.size()) {
                m_node == m_begins[m_current_index];
            } else {
                m_node = m_bucket_sort->m_particles_end;
                return false;
            }
        }
        return true;
    }

    void increment() {
        bool found_good_candidate = false;
        while (!found_good_candidate && go_to_next_candidate()) {

            const double_d p = get<position>(*m_node) + m_transpose[m_current_index];
            m_dx = p - m_r;

            bool outside = false;
            for (int i=0; i < dimension; i++) {
                if (std::abs(m_dx[i]) < m_bucket_sort->m_bucket_side_length[i]) {
                    outside = true;
                    break;
                } 
            }

            found_good_candidate = !outside;
        }
    }

    reference operator *() {
        return dereference();
    }
    reference operator ->() {
        return dereference();
    }
    const_iterator& operator++() {
        increment();
        return *this;
    }
    const_iterator operator++(int) {
        const_iterator tmp(*this);
        operator++();
        return tmp;
    }

    size_t operator-(const_iterator start) const {
        size_t count = 0;
        while (start != *this) {
            start++;
            count++;
        }
        return count;
    }

    inline bool operator==(const const_iterator& rhs) {
        return equal(rhs);
    }

    inline bool operator!=(const const_iterator& rhs){
        return !operator==(rhs);
    }

private:
    friend class boost::iterator_core_access;

    
    const BucketSearch* m_bucket_sort;
    double_d m_r;
    double_d m_dx;
    particles_iterator m_node;
    std::vector<particles_iterator> m_begins;
    std::vector<particles_iterator> m_ends;
    std::vector<double_d> m_transpose;
    unsigned int m_current_index = -1;
};

}


#endif /* BUCKETSEARCH_H_ */
