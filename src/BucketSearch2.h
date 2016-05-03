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

#include "CudaInclude.h"
#include "SpatialUtil.h"
#include "Particles.h"



template <typename traits>
class BucketSearch {
    const static unsigned int m_dimension = traits::dimension;
    typedef typename traits::vector_double_d vector_double_d;
    typedef typename traits::vector_int_d vector_int_d;
    typedef typename traits::vector_unsigned_int_d vector_unsigned_int_d;
    typedef typename traits::vector_bool_d vector_bool_d;
    typedef typename Vector<double,m_dimension> double_d;
    typedef typename Vector<int,m_dimension> int_d;
    typedef typename Vector<unsigned int,m_dimension> unsigned_int_d;
    typedef typename Vector<bool,m_dimension> bool_d;
    typedef typename traits::vector_int vector_int;
    typedef typename traits::vector_unsigned_int vector_unsigned_int;
    typedef typename traits::iterator particles_iterator;
    typedef typename traits::const_iterator const_particles_iterator;
    typedef typename traits::value_type particles_value_type;

public:

    /// A const iterator to a set of neighbouring points. This iterator implements
    /// a STL forward iterator type
	class const_iterator;

    BucketSearch() {};

    void embed_points(particles_iterator &begin, particles_iterator& end) {
        m_particles_begin = begin;
        m_particles_end = end;
        m_positions_begin = get<0>(m_particles_begin);
        m_positions_end = get<0>(m_particles_end);

        CHECK(!bbox.is_empty(), "trying to embed particles into an empty domain. use the function `set_domain` to setup the spatial domain first.");

        m_bucket_begin.resize(m_size.prod());
        m_bucket_end.resize(m_size.prod());
        m_bucket_indices.resize(distance(m_particles_begin,m_particles_end));

        build_buckets();
    }

    /// return a const forward iterator to all the points in the neighbourhood of \p r. If 
    /// this function is being used to find all the point pairs within the same point container, then
    /// a naive looping through and using find_broadphase_neighbours() will find each pair twice. 
    /// This can be avoided by setting self=true and supplying the index of each point with my_index
    const_iterator find_broadphase_neighbours(const Vect3d& r, 
                                              const int my_index, 
                                              const bool self) const;

    /// return an end() iterator to compare against the result of find_broadphase_neighbours in order
    /// to determine the end of the neighbour list
	const_iterator end() const;

    void set_domain(double_d &min_in, double_d &max_in, bool_d&periodic_in, double_d& side_length) {
        LOG(2,"BucketSearch: set_domain:");
        m_bounds.bmin = min_in;
        m_bounds.bmax = max_in;
        m_periodic = periodic_in;
        m_bucket_side_length = side_length;
        m_size =  ((m_bounds.bmax-m_bounds.bmin)/m_bucket_side_length).cast<int>();
        m_bucket_side_length = (m_bounds.bmax-m_bounds.bmin)/m_size;
	    LOG(2,"\tbounds = "<<m_bounds);
	    LOG(2,"\tperiodic = "<<m_periodic);
	    LOG(2,"\tbucket_side_length = "<<m_bucket_side_length);

        embed_points(m_particles_begin,m_particles_end);
    }

    void get_domain(double_d &min_out, double_d &max_out, bool_d&periodic_out, double_d& side_length_out) {
        min_out = m_bounds.bmin;
        max_out = m_bounds.bmax;
        periodic_out = m_periodic;
        side_length_out = m_bucket_side_length;
    }

private:
    void build_buckets();

	inline unsigned int collapse_index_vector(const unsigned_int_d &vindex) const {
        unsigned int index = vindex[0];
        unsigned int multiplier = 1.0;
        for (int i=1; i<m_dimension; i++) {
            multiplier *= m_size[i];
		    ASSERT((vindex[i] > 0) && (vindex[i] < m_size[i]), "index is outside of dimension "<<i<<": "<<vindex);
            index += multiplier*vindex[i];
        }
        return index;
    }

        // find the raster indices of p's bucket
        unsigned int index = (v[0]-box.min[0])/length[0];
        unsigned int multiplier = 1.0;
        for (int i=1; i<m_dimension; i++) {
            multiplier *= m_size[i];
            const unsigned_int raster_d = (v[i]-box.min[i])/m_bucket_side_length[i];
		    ASSERT((raster_d > 0) && (raster_d < m_size[i]), "position is outside of dimension "<<i<<": "<<r);
            index += multiplier*raster_d;
        }
        return index;
    }
 

    // hash a point in the unit square to the index of
    // the grid bucket that contains it
	inline unsigned int find_bucket_index(const Vect3d &r) const {
        // find the raster indices of p's bucket
        unsigned int index = (v[0]-box.min[0])/length[0];
        unsigned int multiplier = 1.0;
        for (int i=1; i<m_dimension; i++) {
            multiplier *= size[i];
            const unsigned_int raster_d = (v[i]-box.min[i])/length[i];
		    ASSERT((raster_d > 0) && (raster_d < m_size[i]), "position is outside of dimension "<<i<<": "<<r);
            index += multiplier*raster_d;
        }
        return index;
    }
     
    template<unsigned int D>
    struct point_to_bucket_index {
        CUDA_HOST_DEVICE
        point_to_bucket_index() {}

        CUDA_HOST_DEVICE
        unsigned int operator()(const double_d& v) const {
            return find_bucket_index(v);
        }
    };

    iterator& m_particles_begin;
    iterator& m_particles_end;
    vector_double_d::iterator& m_positions_begin;
    vector_double_d::iterator& m_positions_end;
    bool_d m_periodic;
    double_d m_bucket_side_length; 
    unsigned_int_d m_size;
    bbox m_bounds;

    // the grid data structure keeps a range per grid bucket:
    // each bucket_begin[i] indexes the first element of bucket i's list of points
    // each bucket_end[i] indexes one past the last element of bucket i's list of points
    vector_unsigned_int m_bucket_begin;
    vector_unsigned_int m_bucket_end;
    vector_unsigned_int m_bucket_indices;

    unsigned_int_d m_surrounding_buckets_offsets;
};


template <typename traits>
BucketSearch<traits>::build_buckets() {

    // transform the points to their bucket indices
    transform(m_positions_begin,
            m_positions_end,
            m_bucket_indices.begin(),
            point_to_bucket_index());

    // sort the points by their bucket index
    sort_by_key(m_bucket_indices.begin(),
            m_bucket_indices.end(),
            m_positions_begin);

    // find the beginning of each bucket's list of points
    counting_iterator<unsigned int> search_begin(0);
    lower_bound(m_bucket_indices.begin(),
            m_bucket_indices.end(),
            search_begin,
            search_begin + w*h,
            m_bucket_begin.begin());

    // find the end of each bucket's list of points
    upper_bound(m_bucket_indices.begin(),
            m_bucket_indices.end(),
            search_begin,
            search_begin + w*h,
            m_bucket_end.begin());
}


template <typename traits>
std::pair<BucketSearch<traits>::const_particles_iterator,
          BucketSearch<traits>::const_particles_iterator> 
BucketSearch<traits>::find_broadphase_neighbours(
        const Vect3d& r, 
        const int my_index, 
        const bool self) const {
    
    const unsigned int my_bucket = point_to_bucket_index(r);
    const_particles_iterator search_iterator(this,r);
    int_d bucket_offset(-1);
    constexpr unsigned int last_d = m_dimension-1;
    bool still_going = true;
    while (still_going) {
        unsigned_int_d other_bucket = my_bucket + bucket_offset; 

        // handle end cases
        for (int i=0; i<m_dimension; i++) {
            if (other_bucket[i] == std::numeric_limits<unsigned int>::max()) {
                if (m_periodic[i]) {
                    other_bucket[i] = m_size[i]-1;
                } else {
                    break;
                }
            }
            if (other_bucket[i] == m_size[i]) {
                if (m_periodic[i]) {
                    other_bucket[i] = 0;
                } else {
                    break;
                }
            }
        }

        const unsigned_int other_bucket_index = collapse_index_vector(other_bucket);
        search_iterator.add_range(m_bucket_begin[other_bucket_index],m_bucket_end[other_bucket_index]);

        // go to next candidate bucket
        for (int i=0; i<m_dimension; i++) {
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
    typedef const std::tuple<const particles_value_type&,const double_d&> value_type;
    typedef const std::tuple<const particles_value_type&,const double_d&> reference;
	typedef std::ptrdiff_t difference_type;

    const_iterator(const BucketSearch<traits>* bucket_sort, double_d &r):
        m_bucket_sort(bucket_sort),
        m_r(r) {}

    void add_range(particles_iterator &begin, particles_iterator &end) {
        m_begins.push_back(begin);
        m_ends.push_back(end);
        if (m_begins.size() == 1) {
            m_current_index = 0;
            m_node = m_begins[m_current_index];
        }
    }

    bool equal(const_iterator const& other) const {
        return m_node == other.m_node;
    }

    reference dereference() const { 
        return std::tie(bucket_sort->begin_iterator[*m_node],dx); 
    }

    void go_to_next_candidate() {
        m_node++;
        if (m_node == m_ends[m_current_index]) {
            m_current_index++;
            if (m_current_index < m_begins.size()) {
                m_node == m_begins[m_current_index];
            } else {
                m_node = m_bucket_sort->m_particles_end;
            }
        }
    }

    void increment() {
        bool found_good_candidate = false;
        while (!(found_good_candidate || (m_node == m_bucket_sort->m_particles_end))) {
            go_to_next_candidate();

            const double_d p = get<position>(m_node);
            m_dx = centre-bucket_sort->correct_position_for_periodicity(r, p);

            bool outside = false;
            for (int i=0; i < m_dimension; i++) {
                if (std::abs(dx[i]) < bucket_sort->m_bucket_side_length[i]) {
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
    unsigned int m_current_index = -1;
};


#endif /* BUCKETSEARCH_H_ */
