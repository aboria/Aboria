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
// Acknowledgement: This source was modified from the Thrust workshop git repository by Jared Hoberock, https://github.com/jaredhoberock/thrust-workshop
//


#ifndef OCTTREE_H_
#define OCTTREE_H_

#include "CudaInclude.h"
#include "detail/SpatialUtil.h"
#include "Particles.h"

namespace Aboria {

template <typename Traits>
class octtree_query; 

template <typename Traits>
class octtree: 
    public neighbour_search_base<octtree<Traits>,
                                 Traits,
                                 octtree_query<Traits>> {

    typedef typename Traits::double_d double_d;
    typedef typename Traits::bool_d bool_d;
    typedef typename Traits::position position;
    typedef typename Traits::vector_double_d_const_iterator vector_double_d_const_iterator;
    typedef typename Traits::vector_unsigned_int_iterator vector_unsigned_int_iterator;
    typedef typename Traits::vector_unsigned_int vector_unsigned_int;
    typedef typename Traits::vector_int vector_int;
    typedef typename Traits::unsigned_int_d unsigned_int_d;
    typedef typename Traits::template vector_type<int2>::type vector_int2;
    static const unsigned int dimension = Traits::dimension;

    // number of children = 2^d
    static const size_t nchild = (1  << dimension);

    typedef typename Traits::iterator iterator;

    typedef neighbour_search_base<octtree<Traits>,
                                 Traits,
                                 octtree_query<Traits>> base_type;

    friend base_type;


public:
    octtree():base_type() {}
    static constexpr bool cheap_copy_and_delete_at_end() {
        return false;
    }

private:


    void set_domain_impl() {
        const size_t n = this->m_particles_end - this->m_particles_begin;
    }
    void update_iterator_impl() {
        this->m_query.m_particles_begin = iterator_to_raw_pointer(this->m_particles_begin);
        this->m_query.m_particles_end = iterator_to_raw_pointer(this->m_particles_end);
    }

    void embed_points_impl() {
        const size_t num_points  = this->m_particles_end - this->m_particles_begin;

        m_tags.resize(num_points);

        /******************************************
         * 3. Classify points                     *
         ******************************************/

        detail::transform(get<position>(this->m_particles_begin), 
                get<position>(this->m_particles_end), 
                m_tags.begin(), 
                classify_point(m_bounds, m_max_level));

                // Now that we have the geometric information, we can sort the
        // points accordingly.
        sort_by_tags();
        
        build_tree();

    }

    void add_points_at_end_impl(const size_t dist) {
        const size_t num_points  = this->m_particles_end - this->m_particles_begin;
        auto start_adding_particles = this->m_particles_end-dist;
        m_tags.resize(num_points);
        auto start_adding_tags = this->m_tags.end()-dist;

        /******************************************
         * 3. Classify new points                 *
         ******************************************/
        detail::transform(get<position>(start_adding_particles), 
                get<position>(this->m_particles_end), 
                start_adding_tags, 
                classify_point(m_bounds, m_max_level));

        // sort and then build tree
        sort_by_tags();
        build_tree();
    }


    void delete_points_at_end_impl(const size_t dist) {
        const size_t n = this->m_particles_end - this->m_particles_begin;

        m_tags.resize(n);
        build_tree();
    }
     
    void copy_points_impl(iterator copy_from_iterator, iterator copy_to_iterator) {
        auto positions_from = get<position>(copy_from_iterator);
        auto positions_to = get<position>(copy_to_iterator);
    }
    const octtree_query<Traits>& get_query_impl() const {
        return this->m_query;
    }

    void sort_by_tags() {
        /******************************************
         * 4. Sort according to classification    *
         ******************************************/
        if (m_tags.size() > 0) {
            m_indices.resize(m_tags.size());
            detail::sequence(m_indices.begin(), m_indices.end());
            detail::sort_by_key(m_tags.begin(), m_tags.end(), m_indices.begin());
            detail::reorder_destructive(m_indices.begin(), 
                                        m_indices.end(), 
                                        this->m_particles_begin);
        }
    }
private:
    void build_tree();

    struct classify_point;
    struct child_index_to_tag_mask;
    struct classify_node;
    struct write_nodes;
    struct make_leaf;

    int max_points;
    int m_max_level;

    vector_int m_tags;
    vector_int m_indices;
    vector_int m_nodes;
    vector_int2 m_leaves;
    detail::bbox<dimension> m_bounds;
};



template <typename traits>
void octtree<traits>::build_tree() {
    m_nodes.clear();
    m_leaves.clear();
    vector_int active_nodes(1,0);

    // Build the tree one level at a time, starting at the root
    for (int level = 1 ; !active_nodes.empty() && level <= m_max_level ; ++level)
    {
        /******************************************
         * 1. Calculate children                  *
         ******************************************/

        // New children: 2^D quadrants per active node
        vector_int children(nchild*active_nodes.size());

        // For each active node, generate the tag mask for each of its 2^D children
        detail::tabulate(children.begin(), children.end(),
                child_index_to_tag_mask(level, m_max_level, active_nodes.data()));

        /******************************************
         * 2. Determine interval for each child   *
         ******************************************/

        // For each child we need interval bounds
        vector_int lower_bounds(children.size());
        vector_int upper_bounds(children.size());

        // Locate lower and upper bounds for points in each quadrant
        detail::lower_bound(m_tags.begin(),
                m_tags.end(),
                children.begin(),
                children.end(),
                lower_bounds.begin());

        int length = (1 << (m_max_level - level) * 2) - 1;

        detail::upper_bound(m_tags.begin(),
                m_tags.end(),
                detail::make_transform_iterator(children.begin(), detail::_1 + length),
                detail::make_transform_iterator(children.end(), detail::_1 + length),
                upper_bounds.begin());

        /******************************************
         * 3. Mark each child as empty/leaf/node  *
         ******************************************/

        // Mark each child as either empty, a node, or a leaf
        vector_int child_node_kind(children.size(), 0);
        detail::transform(lower_bounds.begin(), lower_bounds.end(),
                upper_bounds.begin(),
                child_node_kind.begin(),
                classify_node(this->m_n_particles_in_leaf, level == m_max_level));


        /******************************************
         * 4. Enumerate nodes and leaves          *
         ******************************************/

        // Enumerate the nodes and leaves at this level
        vector_int leaves_on_this_level(child_node_kind.size());
        vector_int nodes_on_this_level(child_node_kind.size());

        // Enumerate nodes at this level
        detail::transform_exclusive_scan(child_node_kind.begin(), 
                child_node_kind.end(), 
                nodes_on_this_level.begin(), 
                detail::is_a<detail::NODE>(), 
                0, 
                std::plus<int>());

        // Enumerate leaves at this level
        detail::transform_exclusive_scan(child_node_kind.begin(), 
                child_node_kind.end(), 
                leaves_on_this_level.begin(), 
                detail::is_a<detail::LEAF>(), 
                0, 
                std::plus<int>());

        int num_nodes_on_this_level = nodes_on_this_level.back() + (child_node_kind.back() == detail::NODE ? 1 : 0);
        int num_leaves_on_this_level = leaves_on_this_level.back() + (child_node_kind.back() == detail::LEAF ? 1 : 0);


        /******************************************
         * 5. Add the children to the node list   *
         ******************************************/

        int num_children = child_node_kind.size();

        int children_begin = m_nodes.size();
        m_nodes.resize(m_nodes.size() + num_children);

        detail::transform(
                detail::make_zip_iterator(
                    detail::make_tuple(child_node_kind.begin(), nodes_on_this_level.begin(), leaves_on_this_level.begin())),
                detail::make_zip_iterator(
                    detail::make_tuple(child_node_kind.end(), nodes_on_this_level.end(), leaves_on_this_level.end())),
                m_nodes.begin() + children_begin,
                write_nodes(m_nodes.size(), m_leaves.size()));


        /******************************************
         * 6. Add the leaves to the leaf list     *
         ******************************************/

        children_begin = m_leaves.size();

        m_leaves.resize(m_leaves.size() + num_leaves_on_this_level);

        detail::scatter_if(detail::make_transform_iterator(
                    detail::make_zip_iterator(
                        detail::make_tuple(lower_bounds.begin(), upper_bounds.begin())),
                    make_leaf()),
                detail::make_transform_iterator(
                    detail::make_zip_iterator(
                        detail::make_tuple(lower_bounds.end(), upper_bounds.end())),
                    make_leaf()),
                leaves_on_this_level.begin(),
                child_node_kind.begin(),
                m_leaves.begin() + children_begin,
                detail::is_a<detail::LEAF>());


        /******************************************
         * 7. Set the nodes for the next level    *
         ******************************************/

        // Set active nodes for the next level to be all the childs nodes from this level
        active_nodes.resize(num_nodes_on_this_level);

        detail::copy_if(children.begin(),
                children.end(),
                child_node_kind.begin(),
                active_nodes.begin(),
                detail::is_a<detail::NODE>());

    }
}

// Classify a point with respect to the bounding box.
template <typename traits>
struct octtree<traits>::classify_point {
    detail::bbox<dimension> box;
    int max_level;

    // Create the classifier
    classify_point(const detail::bbox<dimension> &b, int lvl) : box(b), max_level(lvl) {}

    // Classify a point
    inline CUDA_HOST_DEVICE
    int operator()(const double_d &p) { return point_to_tag(p, box, max_level); }
};


template <typename traits>
struct octtree<traits>::child_index_to_tag_mask {
    const int level, max_level;

    // mask for lower n bits, where n is the number of dimensions
    const static unsigned mask = nchild - 1;

    typedef typename vector_int::const_pointer ptr_type;
    ptr_type m_nodes;

    child_index_to_tag_mask(int lvl, int max_lvl, ptr_type nodes) : level(lvl), max_level(max_lvl), m_nodes(nodes) {}

    inline CUDA_HOST_DEVICE
    int operator()(int idx) const
    {
        int tag = m_nodes[idx/nchild];
        int which_child = (idx&mask);
        return detail::child_tag_mask(tag, which_child, level, max_level);
    }
};


template <typename traits>
struct octtree<traits>::classify_node
{
    const int threshold;
    const int last_level;

    classify_node(int threshold, int last_level) : threshold(threshold), last_level(last_level) {}

    inline CUDA_HOST_DEVICE
    int operator()(int lower_bound, int upper_bound) const
    {
        const int count = upper_bound - lower_bound;
        if (count == 0)
        {
            return detail::EMPTY;
        }
        else if (last_level || count < threshold)
        {
            return detail::LEAF;
        }
        else
        {
            return detail::NODE;
        }
    }
};

template <typename traits>
struct octtree<traits>::write_nodes {
    int num_nodes, num_leaves;

    write_nodes(int num_nodes, int num_leaves) : 
        num_nodes(num_nodes), num_leaves(num_leaves) 
    {}

    template <typename tuple_type>
    inline CUDA_HOST_DEVICE
    int operator()(const tuple_type &t) const
    {
        int node_type = get<0>(t);
        int node_idx  = get<1>(t);
        int leaf_idx  = get<2>(t);

        if (node_type == detail::EMPTY)
        {
            return detail::get_empty_id();
        }
        else if (node_type == detail::LEAF)
        {
            return detail::get_leaf_id(num_leaves + leaf_idx);
        }
        else
        {
            return num_nodes + nchild * node_idx;
        }
    }
};

template <typename traits>
struct octtree<traits>::make_leaf {
    typedef int2 result_type;
    template <typename tuple_type>
        inline CUDA_HOST_DEVICE
        result_type operator()(const tuple_type &t) const
        {
            int x = get<0>(t);
            int y = get<1>(t);

            return result_type(x, y);
        }
};

template <typename Traits>
struct octtree_query {
    const static unsigned int dimension = Traits::dimension;
    typedef detail::nanoflann_kd_tree_type<Traits> kd_tree_type;
    typedef typename kd_tree_type::Node value_type;
    typedef const value_type& reference;
    typedef const value_type* pointer;

    typedef Traits traits_type;
    typedef typename Traits::raw_pointer raw_pointer;
    typedef typename Traits::double_d double_d;
    typedef typename Traits::bool_d bool_d;
    typedef typename Traits::int_d int_d;
    typedef typename Traits::unsigned_int_d unsigned_int_d;
    typedef tree_query_iterator<nanoflann_adaptor_query,-1> query_iterator;
    typedef value_type* root_iterator;
    typedef tree_depth_first_iterator<nanoflann_adaptor_query> all_iterator;
    typedef ranges_iterator<Traits> particle_iterator;

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
    static bool is_leaf_node(reference bucket) {
        return (bucket.child1 == NULL) && (bucket.child2 == NULL);
    }
    static bool is_tree() {
        return true;
    }
    static size_t get_dimension_index(reference bucket) {
        return bucket.node_type.sub.divfeat;
    }
    static double get_cut_low(reference bucket) {
        return bucket.node_type.sub.divlow;
    }
    static double get_cut_high(reference bucket) {
        return bucket.node_type.sub.divhigh;
    }
    static pointer get_child1(pointer bucket) {
	    return bucket->child1;
    }
    static pointer get_child2(pointer bucket) {
	    return bucket->child2;
    }
    /*
     * end functions for tree_query_iterator
     */

    friend std::ostream& operator<<(std::ostream& os, reference bucket) {
        if (is_leaf_node(bucket)) {
            os << "Leaf node";
        } else {
            os << "Node";
        } 
        os <<" with bounding box " << get_bucket_bbox(bucket) << std::endl;
        return os;
    }   
           

    iterator_range<particle_iterator> 
    get_bucket_particles(reference bucket) const {
#ifndef __CUDA_ARCH__
        LOG(4,"\tget_bucket_particles: looking in bucket with bounding box low =  "<<get_bucket_bounds_low(bucket)<<" and high = "<<get_bucket_bounds_high(bucket)<<" idx = "<<get_dimension_index(bucket)<<" start index = "<<bucket.node_type.lr.left<<" end index = "<<bucket.node_type.lr.right);
#endif        
        
        return iterator_range<particle_iterator>(
                        particle_iterator(m_particles_begin + bucket.node_type.lr.left),
                        particle_iterator(m_particles_begin + bucket.node_type.lr.right));
    }

    static double_d
    get_bucket_bounds_low(reference bucket) {
        double_d low;
        for (int i = 0; i < dimension; ++i) {
            low[i] = bucket.bbox[i].low;
        }
        return low;
    }

    static double_d
    get_bucket_bounds_high(reference bucket) {
        double_d high;
        for (int i = 0; i < dimension; ++i) {
            high[i] = bucket.bbox[i].high;
        }
        return high;
    }

     CUDA_HOST_DEVICE
    reference get_bucket(const double_d &position) const {
        pointer node = m_root;
        while(!is_leaf_node(*node)) {
            ASSERT(get_child1(node) != nullptr,"no child1");
            ASSERT(get_child2(node) != nullptr,"no child2");
            const size_t idx = get_dimension_index(*node);
            const double diff_cut_high = position[idx] - get_cut_high(*node);
            const double diff_cut_low = position[idx]- get_cut_low(*node);

            if ((diff_cut_low+diff_cut_high)<0) {
                node = get_child1(node);
            } else {
                node = get_child2(node);
            }
        }
        return *node;
    }

    CUDA_HOST_DEVICE
    size_t get_bucket_index(reference bucket) const {
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
                query_iterator(m_root,position,double_d(max_distance),this),
                query_iterator()
                );
    }

    template <int LNormNumber=-1>
    iterator_range<query_iterator> 
    get_buckets_near_point(const double_d &position, const double_d &max_distance) const {
#ifndef __CUDA_ARCH__
        LOG(4,"\tget_buckets_near_point: position = "<<position<<" max_distance= "<<max_distance);
#endif
        return iterator_range<query_iterator>(
                query_iterator(m_root,position,max_distance,this),
                query_iterator()
                );
    }

    iterator_range<root_iterator> get_root_buckets() const {
        return iterator_range<root_iterator>(m_root, m_root+1);
    }

    iterator_range<all_iterator> get_subtree(reference bucket) const {
        return iterator_range<all_iterator>(all_iterator(&bucket,this),all_iterator());
    }

    raw_pointer get_particles_begin() const {
        return m_particles_begin;
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




// assume that query functions, are only called from device code
template <typename Traits>
struct octtree_query {

    typedef Traits traits_type;
    typedef typename Traits::raw_pointer raw_pointer;
    typedef typename Traits::double_d double_d;
    typedef typename Traits::bool_d bool_d;
    typedef typename Traits::int_d int_d;
    typedef typename Traits::unsigned_int_d unsigned_int_d;
    const static unsigned int dimension = Traits::dimension;
    typedef lattice_iterator<dimension> query_iterator;
    typedef lattice_iterator<dimension> root_iterator;
    typedef lattice_iterator<dimension> all_iterator;
    typedef typename query_iterator::reference reference;
    typedef typename query_iterator::pointer pointer;
    typedef typename query_iterator::value_type value_type;
    typedef ranges_iterator<Traits> particle_iterator;

    raw_pointer m_particles_begin;
    raw_pointer m_particles_end;

    bool_d m_periodic;
    double_d m_bucket_side_length; 
    int_d m_end_bucket;
    detail::bbox<dimension> m_bounds;
    detail::point_to_bucket_index<dimension> m_point_to_bucket_index;

    unsigned int *m_bucket_begin;
    unsigned int *m_bucket_end;
    unsigned int m_nbuckets;

    inline
    CUDA_HOST_DEVICE
    bucket_search_parallel_query():
        m_periodic(),
        m_particles_begin(),
        m_bucket_begin()
    {}

    /*
     * functions for trees
     */
    static bool is_leaf_node(const value_type& bucket) {
        return true;
    }

    static bool is_tree() {
        return false;
    }

    // dodgy hack cause nullptr cannot be converted to pointer
    static const pointer get_child1(const pointer& bucket) {
        CHECK(false,"this should not be called")
	    return pointer(-1);
    }
    static const pointer get_child2(const pointer& bucket) {
        CHECK(false,"this should not be called")
	    return pointer(-1);
    }

    //const double_d& get_min_bucket_size() const { return m_bucket_side_length; }
    const double_d& get_bounds_low() const { return m_bounds.bmin; }
    const double_d& get_bounds_high() const { return m_bounds.bmax; }
    const bool_d& get_periodic() const { return m_periodic; }

    CUDA_HOST_DEVICE
    iterator_range<particle_iterator> get_bucket_particles(const reference bucket) const {
        ASSERT((bucket>=int_d(0)).all() && (bucket <= m_end_bucket).all(), "invalid bucket");

        const unsigned int bucket_index = m_point_to_bucket_index.collapse_index_vector(bucket);
        const unsigned int range_start_index = m_bucket_begin[bucket_index]; 
        const unsigned int range_end_index = m_bucket_end[bucket_index]; 

#ifndef __CUDA_ARCH__
        LOG(4,"\tlooking in bucket "<<bucket<<" = "<<bucket_index<<". found "<<range_end_index-range_start_index<<" particles");
#endif
        return iterator_range<particle_iterator>(
                particle_iterator(m_particles_begin + range_start_index),
                particle_iterator(m_particles_begin + range_end_index));
    }

    CUDA_HOST_DEVICE
    detail::bbox<dimension> get_bucket_bbox(const reference bucket) const {
        return detail::bbox<dimension>(
                bucket*m_bucket_side_length + m_bounds.bmin,
                (bucket+1)*m_bucket_side_length + m_bounds.bmin
                );
    }

    double_d get_bucket_bounds_low(const reference bucket) const {
        return bucket*m_bucket_side_length + m_bounds.bmin;
    }

    double_d get_bucket_bounds_high(const reference bucket) const {
        return (bucket+1)*m_bucket_side_length + m_bounds.bmin;
    }

    CUDA_HOST_DEVICE
    value_type get_bucket(const double_d &position) const {
        return m_point_to_bucket_index.find_bucket_index_vector(position);
    }

    CUDA_HOST_DEVICE
    size_t get_bucket_index(const reference bucket) const {
        return m_point_to_bucket_index.collapse_index_vector(bucket);
    }

    template <int LNormNumber=-1>
    CUDA_HOST_DEVICE
    iterator_range<query_iterator> 
    get_buckets_near_point(const double_d &position, const double max_distance) const {
        return get_buckets_near_point(position,double_d(max_distance));
    }

    template <int LNormNumber=-1>
    CUDA_HOST_DEVICE
    iterator_range<query_iterator> 
    get_buckets_near_point(const double_d &position, const double_d& max_distance) const {
#ifndef __CUDA_ARCH__
        LOG(4,"\tget_buckets_near_point: position = "<<position<<" max_distance = "<<max_distance);
#endif
 
        value_type bucket = m_point_to_bucket_index.find_bucket_index_vector(position);
        int_d start = m_point_to_bucket_index.find_bucket_index_vector(position-max_distance);
        int_d end = m_point_to_bucket_index.find_bucket_index_vector(position+max_distance);

        bool no_buckets = false;
        for (int i=0; i<Traits::dimension; i++) {
            if (start[i] < 0) {
                start[i] = 0;
            } else if (start[i] > m_end_bucket[i]) {
                no_buckets = true;
                start[i] = m_end_bucket[i];
            }
            if (end[i] < 0) {
                no_buckets = true;
                end[i] = 0;
            } else if (end[i] > m_end_bucket[i]) {
                end[i] = m_end_bucket[i];
            }
        }
#ifndef __CUDA_ARCH__
        LOG(4,"\tget_buckets_near_point: looking in bucket "<<bucket<<". start = "<<start<<" end = "<<end);
#endif
        if (no_buckets) {
            return iterator_range<query_iterator>(
                    query_iterator(end,end,end)
                    ,query_iterator(end,end,end)
                    );
        } else {
            return iterator_range<query_iterator>(
                    query_iterator(start,end,start)
                    ,++query_iterator(start,end,end)
                    );
        }
    }

    CUDA_HOST_DEVICE
    iterator_range<root_iterator> get_root_buckets() const {
        return iterator_range<query_iterator>(
                root_iterator(int_d(0),m_end_bucket,int_d(0)),
                ++root_iterator(int_d(0),m_end_bucket,m_end_bucket)
                );
    }

    iterator_range<all_iterator> get_subtree(reference bucket) const {
        return iterator_range<all_iterator>(
                all_iterator(bucket,bucket,bucket),
                ++all_iterator(bucket,bucket,bucket));
    }


    size_t number_of_buckets() const {
        return (m_end_bucket+1).prod();
    }

    raw_pointer get_particles_begin() const {
        return m_particles_begin;
    }


    /*
    CUDA_HOST_DEVICE
    bool get_children_buckets(const bucket_reference &bucket, std::array<value_type,2>& children) {
        return false;
    }

    CUDA_HOST_DEVICE
    iterator_range<query_iterator> get_root_buckets() const {
        return iterator_range<query_iterator>(
                query_iterator(int_d(0),m_end_bucket,int_d(0)),
                ++query_iterator(int_d(0),m_end_bucket,m_end_bucket)
                );
    }
    */
};

   


}





}
#endif /* OCTTREE_H_ */
