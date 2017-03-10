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


#ifndef NANOFLANN_ADAPTOR_H_
#define NANOFLANN_ADAPTOR_H_

#include <boost/iterator/iterator_facade.hpp>
#include "Traits.h"
#include "Vector.h"
#include "Log.h"
#include "Get.h"
#include "NeighbourSearchBase.h"
#include "detail/SpatialUtil.h"
#include <vector>
#include <iostream>
#include <set>

namespace Aboria {

template <typename Traits>
class nanoflann_adaptor; 

template <typename Traits>
class nanoflann_adaptor_query; 

}

#include "nanoflann/nanoflann.hpp"

namespace Aboria {

namespace detail {

template< typename order_iterator, typename value_iterator >
void reorder_destructive( order_iterator order_begin, order_iterator order_end, value_iterator v )  {
    typedef typename std::iterator_traits< value_iterator >::value_type value_t;
    typedef typename std::iterator_traits< order_iterator >::value_type index_t;
    typedef typename std::iterator_traits< order_iterator >::difference_type diff_t;

    diff_t remaining = order_end - 1 - order_begin;
    for ( index_t s = index_t(); remaining > 0; ++ s ) {
        index_t d = order_begin[s];
        if ( d == (diff_t) -1 ) continue;
        -- remaining;
        value_t temp = v[s];
        for ( index_t d2; d != s; d = d2 ) {
            std::swap( temp, v[d] );
            std::swap( order_begin[d], d2 = (diff_t) -1 );
            -- remaining;
        }
        v[s] = temp;
    }
}


template <typename Traits>
using nanoflann_kd_tree_type = 
        nanoflann::KDTreeSingleIndexAdaptor<
            nanoflann::L_inf_Adaptor<double, nanoflann_adaptor<Traits> > ,
            nanoflann_adaptor<Traits>,
            Traits::dimension 
        >;
}

/// \brief Implements neighbourhood searching using a bucket search algorithm, dividing
/// the domain into constant size "buckets".
///
/// This class implements neighbourhood searching using a bucket search algorithm. The 
/// domain is first divided up into a regular grid of constant size "buckets", either by
/// using the class constructor to initialise the domain extents, bucket size etc., or by 
/// using the reset() member function to reset these parameters.
///
/// After the buckets are created, a set of 3D points can be assigned to their respective buckets
/// using the embed_points() member function. After this, neighbourhood queries around
/// a given point can be performed using find_broadphase_neighbours(), which returns a const 
/// iterator to all the points in the same bucket or surrounding buckets of the given point.
///
template <typename Traits>
class nanoflann_adaptor: 
    public neighbour_search_base<nanoflann_adaptor<Traits>,
                                 Traits,
                                 nanoflann_adaptor_query<Traits>> {

    typedef typename Traits::double_d double_d;
    typedef typename Traits::position position;
    typedef typename Traits::vector_int vector_int;
    typedef typename Traits::iterator iterator;
    typedef typename Traits::unsigned_int_d unsigned_int_d;
    static const unsigned int dimension = Traits::dimension;

    typedef neighbour_search_base<nanoflann_adaptor<Traits>,
                                 Traits,
                                 nanoflann_adaptor_query<Traits>> base_type;
    friend base_type;

    typedef detail::nanoflann_kd_tree_type<Traits> kd_tree_type;


public:

    nanoflann_adaptor():
        base_type(), 
        m_kd_tree(nullptr) 
    {}

    ~nanoflann_adaptor() {
        delete m_kd_tree;
    }


    static constexpr bool cheap_copy_and_delete_at_end() {
        return false;
    }

    // Must return the number of data points
	inline size_t kdtree_get_point_count() const { 
        return std::distance(this->m_particles_begin,this->m_particles_end);
    }

	// Returns the distance between the vector "p1[0:size-1]" 
    // and the data point with index "idx_p2" stored in the class:
	inline double kdtree_distance(const double *p1, const size_t idx_p2,size_t /*size*/) const
	{
        size_t ret = 0;
        const double_d& p2 = *(get<position>(this->m_particles_begin)+idx_p2);
        for (int i = 0; i < dimension; ++i) {
           ret += (p1[i]-p2[i])*(p1[i]-p2[i]); 
        }
		return ret;
	}

	// Returns the dim'th component of the idx'th point in the class:
	// Since this is inlined and the "dim" argument is typically an immediate value, the
	//  "if/else's" are actually solved at compile time.
	inline double kdtree_get_pt(const size_t idx, int dim) const
	{
        const double_d& p = *(get<position>(this->m_particles_begin)+idx);
        return p[dim];
	}

	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX& bb) const { 
        for (int i = 0; i < dimension; ++i) {
            bb[i].low = this->m_bounds.bmin[i];
        }
	    return true;
    }

private:
    void set_domain_impl() {
        CHECK(!this->m_periodic.any(),"kd-tree does not work (yet) with periodic boundaries");

        this->m_query.m_bounds.bmin = this->m_bounds.bmin;
        this->m_query.m_bounds.bmax = this->m_bounds.bmax;
        this->m_query.m_periodic = this->m_periodic;
    }


    void update_iterator_impl() {
        this->m_query.m_particles_begin = iterator_to_raw_pointer(this->m_particles_begin);
    }


    void embed_points_impl() {
        if (m_kd_tree != nullptr) {
            delete m_kd_tree;
        } else {
            m_kd_tree = new kd_tree_type(
                                dimension, 
                                *this, 
                                nanoflann::KDTreeSingleIndexAdaptorParams(
                                    this->m_n_particles_in_leaf
                                    ) 
                            );
        }
	    m_kd_tree->buildIndex();

        detail::reorder_destructive(
                m_kd_tree->get_vind().begin(), 
                m_kd_tree->get_vind().end(), 
                this->m_particles_begin);

        this->m_query.m_root = m_kd_tree->get_root();
        this->m_query.m_particles_begin = iterator_to_raw_pointer(this->m_particles_begin);
    }


    void add_points_at_end_impl(const size_t dist) {
        embed_points_impl();
    }

    void delete_points_at_end_impl(const size_t dist) {
        embed_points_impl();
    }

    void copy_points_impl(iterator copy_from_iterator, iterator copy_to_iterator) {
        embed_points_impl();
    }

    const nanoflann_adaptor_query<Traits>& get_query_impl() const {
        return m_query;
    }


    kd_tree_type* m_kd_tree;
    nanoflann_adaptor_query<Traits> m_query;
};



// this is NOT going to work from device code because we are adapting
// a host code only library
template <typename Traits>
struct nanoflann_adaptor_query {
    const static unsigned int dimension = Traits::dimension;
    typedef detail::nanoflann_kd_tree_type<Traits> kd_tree_type;
    typedef typename kd_tree_type::Node value_type;
    typedef value_type& reference;
    typedef value_type* pointer;

    typedef Traits traits_type;
    typedef typename Traits::raw_pointer raw_pointer;
    typedef typename Traits::double_d double_d;
    typedef typename Traits::bool_d bool_d;
    typedef typename Traits::int_d int_d;
    typedef typename Traits::unsigned_int_d unsigned_int_d;
    typedef tree_query_iterator<dimension,nanoflann_adaptor_query,-1> query_iterator;
    typedef ranges_iterator<Traits> particle_iterator;

    bool_d m_periodic;
    detail::bbox<dimension> m_bounds;
    raw_pointer m_particles_begin;

    pointer m_root;

    /*
     * functions for tree_query_iterator
     */
    bool is_leaf_node(reference bucket) {
        return (bucket->child1 == NULL) && (bucket->child2 == NULL);
    }
    size_t get_dimension_index(reference bucket) {
        return bucket->node_type.sub.divfeat;
    }
    double get_cut_low(reference bucket) {
        return bucket->node_type.sub.divlow;
    }
    double get_cut_high(reference bucket) {
        return bucket->node_type.sub.divhigh;
    }
    pointer get_child1(pointer bucket) {
	    return bucket->child1;
    }
    pointer get_child2(pointer bucket) {
	    return bucket->child2;
    }
    /*
     * end functions for tree_query_iterator
     */
           

    iterator_range_with_transpose<particle_iterator> 
    get_bucket_particles(const reference bucket) const {
        ASSERT(!m_periodic.any(), "ERROR: kdtree doesnt work with periodic (yet)");
        double_d transpose(0); 

#ifndef __CUDA_ARCH__
        LOG(4,"\tget_bucket_particles: looking in bucket "<<bucket);
#endif        
        
        return iterator_range_with_transpose<particle_iterator>(
                        particle_iterator(m_particles_begin + bucket.node_type.lr.left),
                        particle_iterator(m_particles_begin + bucket.node_type.lr.right),
                        transpose);
    }

    detail::bbox<dimension> 
    get_bucket_bbox(const reference bucket) const {
        return detail::bbox<dimension>(bucket->bbox.low,bucket->bbox.high);
    }

    iterator_range<query_iterator> 
    get_buckets_near_point(const double_d &position, const double max_distance) const {
#ifndef __CUDA_ARCH__
        LOG(4,"\tget_buckets_near_point: position = "<<position<<" max_distance= "<<max_distance);
#endif
        return iterator_range<query_iterator>(
                query_iterator(position,max_distance),
                query_iterator()
                );
    }

    /*
    bool get_children_buckets(const reference &bucket, std::array<value_type,2>& children) {
		if ((bucket->child1 == NULL)&&(bucket->child2 == NULL)) {
            return false;
        } else {
            children[0] = bucket.child1;
            children[1] = bucket.child2;
            return true;
        }
    }

    iterator_range<query_iterator> get_root_buckets() const {
        m_query_nodes.clear();
        m_query_nodes.push_back(m_kd_tree.get_root_node());
        return iterator_range<query_iterator>(
                m_query_nodes.begin(),
                m_query_nodes.end()
                );
    }
    */

};

}

#endif /* NANOFLANN_ADAPTOR_H_ */
