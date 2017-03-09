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
class nanoflann_adaptor_query; 

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
                                 nanoflann_adaptor_params<Traits>,
                                 nanoflann_adaptor_query<Traits>> {

    typedef typename Traits::double_d double_d;
    typedef typename Traits::position position;
    typedef typename Traits::vector_int vector_int;
    typedef typename Traits::iterator iterator;
    typedef typename Traits::unsigned_int_d unsigned_int_d;

    typedef neighbour_search_base<nanoflann_adaptor<Traits>,
                                 Traits,
                                 nanoflann_adaptor_params<Traits>,
                                 nanoflann_adaptor_query<Traits>> base_type;
    friend base_type;

    typdef nanoflann_adaptor<Traits>
    typedef nanoflann::KDTreeSingleIndexAdaptor<
		L_inf_Adaptor<double, nanoflann_adaptor<Traits> > ,
		nanoflann_adaptor<Traits>,
		dimension 
		> kd_tree_type;


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
        return std::distance(m_particles_begin,m_particles_end);
    }

	// Returns the distance between the vector "p1[0:size-1]" 
    // and the data point with index "idx_p2" stored in the class:
	inline double kdtree_distance(const double *p1, const size_t idx_p2,size_t /*size*/) const
	{
        size_t ret = 0;
        const double_d& p2 = *(get<position>(m_particles_begin)+idx_p2);
        for (int i = 0; i < dimension; ++i) {
           ret += (p1[i]-p2[i])*(p1[i]-p2[i]); 
        }
		return ret;
	}

	// Returns the dim'th component of the idx'th point in the class:
	// Since this is inlined and the "dim" argument is typically an immediate value, the
	//  "if/else's" are actually solved at compile time.
	inline T kdtree_get_pt(const size_t idx, int dim) const
	{
        const double_d& p = *(get<position>(m_particles_begin)+idx);
        return p[dim];
	}

	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX& bb) const { 
        for (int i = 0; i < dimension; ++i) {
            bb[i].low = m_bounds.bmin[i];
        }
	    return true;
    }

private:
    void set_domain_impl() {
        CHECK(!m_periodic.any(),"kd-tree does not work (yet) with periodic boundaries");

        this->m_query.m_buckets_begin = iterator_to_raw_pointer(m_buckets.begin());
        this->m_query.m_bucket_side_length = this->m_bucket_side_length;
        this->m_query.m_bounds.bmin = this->m_bounds.bmin;
        this->m_query.m_bounds.bmax = this->m_bounds.bmax;
        this->m_query.m_periodic = this->m_periodic;
        this->m_query.m_end_bucket = m_size-1;
        this->m_query.m_point_to_bucket_index = m_point_to_bucket_index;
    }


    void update_iterator_impl() {
        this->m_query.m_particles_begin = iterator_to_raw_pointer(this->m_particles_begin);
    }


    void embed_points_impl() {
        if (m_kd_tree != nullptr) {
            delete kd_tree;
        } else {
            kd_tree = new kd_tree_type(dimension, 
                                        *this, 
                                        nanoflann::KDTreeSingleIndexAdaptorParams(m_n_particles_in_leaf) 
                                        );
        }
	    kd_tree->buildIndex();

        this->m_query.m_particles_begin = iterator_to_raw_pointer(this->m_particles_begin);
        this->m_query.m_linked_list_begin = iterator_to_raw_pointer(this->m_linked_list.begin());
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

    typedef Traits traits_type;
    typedef typename Traits::raw_pointer raw_pointer;
    typedef typename Traits::double_d double_d;
    typedef typename Traits::bool_d bool_d;
    typedef typename Traits::int_d int_d;
    typedef typename Traits::unsigned_int_d unsigned_int_d;
    typedef typename Traits::reference reference;
    typedef typename Traits::position position;
    const static unsigned int dimension = Traits::dimension;
    typedef std::vector<NodePtr> node_list_type;
    typedef node_list_type::iterator bucket_iterator;
    typedef typename bucket_iterator::reference bucket_reference;
    typedef typename bucket_iterator::value_type bucket_value_type;
    typedef index_vector_iterator<Traits> particle_iterator;

    bool_d m_periodic;
    detail::bbox<dimension> m_bounds;

    raw_pointer m_particles_begin;
    node_list_type m_query_nodes;
    NodePtr m_root_node;
    double m_max_distance;
    double_d m_ref_point;

    inline
    nanoflann_adaptor_query():
        m_periodic(),
        m_particles_begin(),
        m_buckets_begin()
    {}

    inline void addLeaf(const NodePtr node) {
#ifndef NDEBUG
        double dist = 0;
        for (int i = 0; i < dimension; ++i) {
            dist += std::min(
                     distance.accum_dist(m_ref_point[i], node->bbox.low, dimension)
                    ,distance.accum_dist(m_ref_point[i], node->bbox.high, dimension)
                    );
        }
        ASSERT(dist < m_max_distance, "trying to add a leaf not within max distance, something wrong with kdtree");
#endif
        m_query_nodes.push_back(node);
    }

    inline DistanceType worstDist() const { return m_max_distance; }

    iterator_range_with_transpose<particle_iterator> get_bucket_particles(const bucket_reference &bucket) const {
        ASSERT(!m_periodic.any(), "ERROR: kdtree doesnt work with periodic (yet)");
        double_d transpose(0); 

#ifndef __CUDA_ARCH__
        LOG(4,"\tget_bucket_particles: looking in bucket "<<bucket);
#endif        
    }

    detail::bbox<dimension> get_bucket_bbox(const bucket_reference &bucket) const {
        return detail::bbox<dimension>(bucket.bbox.low,bucket.bbox.high);
    }

    iterator_range<bucket_iterator> get_buckets_near_point(const double_d &position, const double max_distance) const {
#ifndef __CUDA_ARCH__
        LOG(4,"\tget_buckets_near_point: position = "<<position<<" max_radius = "<<);
#endif
 
        m_query_nodes.clear();

        // this fills m_query_nodes
        kd_tree.findNeighborLeafs(*this, position.data(), nanoflann::SearchParams());

        return iterator_range<bucket_iterator>(
                m_query_nodes.begin(),
                m_query_nodes.end()
                );
    }

    bool get_children_buckets(const bucket_reference &bucket, std::array<bucket_value_type,2>& children) {
		if ((bucket->child1 == NULL)&&(bucket->child2 == NULL)) {
            return false;
        } else {
            children[0] = bucket.child1;
            children[1] = bucket.child2;
            return true;
        }
    }

    iterator_range<bucket_iterator> get_root_buckets() const {
        return iterator_range<bucket_iterator>(
                m_root_node,
                ++m_root_node,
                );
    }

};

}

#endif /* NANOFLANN_ADAPTOR_H_ */
