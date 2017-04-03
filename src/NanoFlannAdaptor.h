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

// reorders v such that v_new[i] == v[order[i]]
template< typename order_iterator, typename value_iterator >
void reorder_destructive( order_iterator order_begin, order_iterator order_end, value_iterator v )  {
    typedef typename std::iterator_traits< value_iterator >::value_type value_t;
    typedef typename std::iterator_traits< value_iterator >::reference reference;
    typedef typename std::iterator_traits< order_iterator >::value_type index_t;
    typedef typename std::iterator_traits< order_iterator >::difference_type diff_t;

    diff_t size = order_end - order_begin;
    
    size_t i, j, k;
    value_t temp;
    for(i = 0; i < size; i++){
        if(i != order_begin[i]){
            temp = v[i];
            k = i;
            while(i != (j = order_begin[k])){
                // every move places a value in it's final location
                // NOTE: need the static cast or assignment has no effect (TODO: why?)
                static_cast<reference>(v[k])  = v[j];
                order_begin[k] = k;
                k = j;
            }
            v[k] = temp;
            order_begin[k] = k;
        }
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
    typedef typename kd_tree_type::Node node_type;


public:

    nanoflann_adaptor():
        base_type(), 
        m_kd_tree(dimension,*this) 
    {}

    //~nanoflann_adaptor() {}


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
            bb[i].high = this->m_bounds.bmax[i];
        }
	    return true;
    }

private:
    void set_domain_impl() {
        m_kd_tree.set_leaf_max_size(this->m_n_particles_in_leaf);

        this->m_query.m_bounds.bmin = this->m_bounds.bmin;
        this->m_query.m_bounds.bmax = this->m_bounds.bmax;
        this->m_query.m_periodic = this->m_periodic;
    }


    void update_iterator_impl() {
        this->m_query.m_particles_begin = iterator_to_raw_pointer(this->m_particles_begin);
    }

    void print_tree(const node_type* nodes) {
#ifndef __CUDA_ARCH__
        if (4 <= ABORIA_LOG_LEVEL) { 
            std::vector<const node_type*> new_nodes;
            new_nodes.push_back(nodes);
            print_level(new_nodes);
        }
#endif
    }

    void print_level(std::vector<const node_type*> &nodes) {
#ifndef __CUDA_ARCH__
        if (4 <= ABORIA_LOG_LEVEL) { 
            std::vector<const node_type*> new_nodes;
            LOG(4,"printing level with "<<nodes.size()<<" nodes");
            for (const node_type* ptr: nodes) {
                if (this->m_query.is_leaf_node(*ptr)) {
                    const int idx = this->m_query.get_dimension_index(*ptr);
                    const double_d bbox_low = this->m_query.get_bucket_bounds_low(*ptr);
                    const double_d bbox_high = this->m_query.get_bucket_bounds_high(*ptr);
                    const int start_index = ptr->node_type.lr.left;
                    const int end_index = ptr->node_type.lr.right;
                    LOG(4,"\tleaf node with idx = "<<idx<<" box low =  "<<bbox_low<<" and high = "<<bbox_high<<" start index = "<<start_index<<" end index = "<<end_index);
                    LOG(4,"\tparticles in bucket are:");
                    for (int i = start_index; i < end_index; ++i) {
                        const double_d p = get<position>(this->m_particles_begin)[i];
                        LOG(4,"\t\tposition = "<<p);
                    }
                } else {
                    const int idx = this->m_query.get_dimension_index(*ptr);
                    const double_d bbox_low = this->m_query.get_bucket_bounds_low(*ptr);
                    const double_d bbox_high = this->m_query.get_bucket_bounds_high(*ptr);
                    const double cut_low = this->m_query.get_cut_low(*ptr);
                    const double cut_high = this->m_query.get_cut_high(*ptr);
                    LOG(4,"\tNOT leaf node with idx = "<<idx<<" cut_low= "<<cut_low<<" cut_high = "<<cut_high);
                    const node_type* child1 = this->m_query.get_child1(ptr);
                    const node_type* child2 = this->m_query.get_child2(ptr);
                    new_nodes.push_back(child1);
                    new_nodes.push_back(child2);
                }
            }
            if (new_nodes.size() > 0) {
                print_level(new_nodes);
            } else {
                LOG(4,"finished tree");
            }
        }
#endif
    }


    void embed_points_impl() {
	    m_kd_tree.buildIndex();

        detail::reorder_destructive(
                m_kd_tree.get_vind().begin(), 
                m_kd_tree.get_vind().end(), 
                this->m_particles_begin);


        this->m_query.m_root = m_kd_tree.get_root_node();
        this->m_query.m_particles_begin = iterator_to_raw_pointer(this->m_particles_begin);
        this->m_query.m_number_of_buckets = m_kd_tree.size_nodes();

        print_tree(m_kd_tree.get_root_node());
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


    kd_tree_type m_kd_tree;
    nanoflann_adaptor_query<Traits> m_query;
};


// this is NOT going to work from device code because we are adapting
// a host code only library
template <typename Traits>
struct nanoflann_adaptor_query {
    const static unsigned int dimension = Traits::dimension;
    typedef detail::nanoflann_kd_tree_type<Traits> kd_tree_type;
    typedef typename kd_tree_type::Node value_type;

    typedef Traits traits_type;
    typedef typename Traits::raw_pointer raw_pointer;
    typedef typename Traits::double_d double_d;
    typedef typename Traits::bool_d bool_d;
    typedef typename Traits::int_d int_d;
    typedef typename Traits::unsigned_int_d unsigned_int_d;
    typedef tree_query_iterator<dimension,nanoflann_adaptor_query,-1> query_iterator;
    typedef value_type* root_iterator;
    typedef tree_depth_first_iterator<dimension,nanoflann_adaptor_query> all_iterator;
    typedef ranges_iterator<Traits> particle_iterator;
    typedef typename query_iterator::reference reference;
    typedef typename query_iterator::value_type value_type;
    typedef linked_list_iterator<Traits> particle_iterator;



    bool_d m_periodic;
    detail::bbox<dimension> m_bounds;
    raw_pointer m_particles_begin;
    size_t m_number_of_buckets;

    value_type* m_root;

    const double_d& get_bounds_low() const { return m_bounds.bmin; }
    const double_d& get_bounds_high() const { return m_bounds.bmax; }
    const bool_d& get_periodic() const { return m_periodic; }

    /*
     * functions for tree_query_iterator
     */
    static bool get_max_levels() {
        return 5;
    }
    static bool is_leaf_node(const value_type& bucket) {
        return (bucket.child1 == NULL) && (bucket.child2 == NULL);
    }
    static size_t get_dimension_index(const value_type& bucket) {
        return bucket.node_type.sub.divfeat;
    }
    static double get_cut_low(const value_type& bucket) {
        return bucket.node_type.sub.divlow;
    }
    static double get_cut_high(const value_type& bucket) {
        return bucket.node_type.sub.divhigh;
    }
    static const value_type* get_child1(const value_type* bucket) {
	    return bucket->child1;
    }
    static const value_type* get_child2(const value_type* bucket) {
	    return bucket->child2;
    }
    /*
     * end functions for tree_query_iterator
     */

    friend std::ostream& operator<<(std::ostream& os, const value_type& bucket) {
        if (is_leaf_node(bucket)) {
            os << "Leaf node";
        } else {
            os << "Node";
        } 
        os <<" with bounding box " << get_bucket_bbox(bucket) << std::endl;
        return os;
    }   
           

    iterator_range<particle_iterator> 
    get_bucket_particles(const value_type& bucket) const {
#ifndef __CUDA_ARCH__
        LOG(4,"\tget_bucket_particles: looking in bucket with bounding box low =  "<<get_bucket_bounds_low(bucket)<<" and high = "<<get_bucket_bounds_high(bucket)<<" idx = "<<get_dimension_index(bucket)<<" start index = "<<bucket.node_type.lr.left<<" end index = "<<bucket.node_type.lr.right);
#endif        
        
        return iterator_range<particle_iterator>(
                        particle_iterator(m_particles_begin + bucket.node_type.lr.left),
                        particle_iterator(m_particles_begin + bucket.node_type.lr.right));
    }

    static double_d
    get_bucket_bounds_low(const value_type& bucket) {
        double_d low;
        for (int i = 0; i < dimension; ++i) {
            low[i] = bucket.bbox[i].low;
        }
        return low;
    }

    static double_d
    get_bucket_bounds_high(const value_type& bucket) {
        double_d high;
        for (int i = 0; i < dimension; ++i) {
            high[i] = bucket.bbox[i].high;
        }
        return high;
    }

     CUDA_HOST_DEVICE
    value_type& get_bucket(const double_d &position) const {
        value_type* node = m_root;
        while(!is_leaf_node(*node)) {
            ASSERT(get_child1(m_node) != nullptr,"no child1");
            ASSERT(get_child2(m_node) != nullptr,"no child2");
            const size_t idx = get_dimension_index(*m_node);
            const double diff_cut_high = position[idx] - get_cut_high(*m_node);
            const double diff_cut_low = position[idx]- get_cut_low(*m_node);

            if ((diff_cut_low+diff_cut_high)<0) {
                node = get_child1(node);
            } else {
                node = get_child2(node);
            }
        }
        return *node;
    }

    CUDA_HOST_DEVICE
    size_t get_bucket_index(const value_type& bucket) const {
        return bucket.index;
    }

    size_t number_of_buckets() const {
        return m_number_of_buckets;
    }

    template <int LNormNumber=-1>
    iterator_range<query_iterator> 
    get_buckets_near_point(const double_d &position, const double max_distance) const {
#ifndef __CUDA_ARCH__
        LOG(4,"\tget_buckets_near_point: position = "<<position<<" max_distance= "<<max_distance);
#endif
        return iterator_range<query_iterator>(
                query_iterator(m_root,position,max_distance,this),
                query_iterator()
                );
    }

    

    iterator_range<root_iterator> get_root_buckets() const {
        return iterator_range<root_iterator>(&m_root, &m_root+1);
    }

    iterator_range<all_iterator> get_all_buckets() const {
        return iterator_range<all_iterator>(all_iterator(m_root),all_iterator());
    }


    /*
    CUDA_HOST_DEVICE
    iterator_range<theta_iterator> get_theta_buckets(const reference bucket) const {
        return iterator_range<theta_iterator>(
                theta_iterator(m_root,bucket),
                theta_iterator()
                );
    }
    */

};






}

#endif /* NANOFLANN_ADAPTOR_H_ */
