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
        set_domain_impl();
        const size_t n = this->m_particles_end - this->m_particles_begin;
    }

    void add_points_at_end_impl(const size_t dist) {
        set_domain_impl();
        auto start_adding = this->m_particles_end-dist;
    }


    void delete_points_at_end_impl(const size_t dist) {
        set_domain_impl();
        const size_t n = this->m_particles_end - this->m_particles_begin;
    }
     
    void copy_points_impl(iterator copy_from_iterator, iterator copy_to_iterator) {
        auto positions_from = get<position>(copy_from_iterator);
        auto positions_to = get<position>(copy_to_iterator);
    }
    const octtree_query<Traits>& get_query_impl() const {
        return m_query;
    }

     
private:
    void build_tree();

    struct classify_point;
    struct child_index_to_tag_mask;
    struct classify_node;
    struct write_nodes;
    struct make_leaf;


    int max_points;
    int max_level;
    int threshold;
    bool_d periodic;

    vector_int m_tags;
    vector_int m_nodes;
    vector_int m_indices;
    detail::bbox<dimension> m_bounds;
    vector_int2 m_leaves;
};


template <typename traits>
void octtree<traits>::embed_points(vector_double_d& points) {

    const int num_points = points.size();

    m_nodes.clear();
    m_leaves.clear();
    m_tags.resize(num_points);
    m_indices.resize(num_points);

    /******************************************
     * 2. Compute bounding box TOOK THIS OUT                *
     ******************************************/

    //detail::bbox<dimension> points_bounds = detail::reduce(points.begin(), points.end(), detail::bbox<dimension>(), std::plus<detail::bbox<dimension> >());
    //assert(points_bounds < m_bounds);

    /******************************************
     * 3. Classify points                     *
     ******************************************/

    detail::transform(points.begin(), 
            points.end(), 
            m_tags.begin(), 
            classify_point(m_bounds, max_level));

    /******************************************
     * 4. Sort according to classification    *
     ******************************************/


    // Now that we have the geometric information, we can sort the
    // points accordingly.
    detail::sequence(m_indices.begin(), m_indices.end());
    detail::sort_by_key(m_tags.begin(), m_tags.end(), m_indices.begin());

    build_tree();
}



template <typename traits>
void octTree<traits>::build_tree() {
  vector_int active_nodes(1,0);

  // Build the tree one level at a time, starting at the root
  for (int level = 1 ; !active_nodes.empty() && level <= max_level ; ++level)
  {
    /******************************************
     * 1. Calculate children                  *
     ******************************************/

    // New children: 2^D quadrants per active node
    constexpr size_t nchild = detail::ipow(2,dimension);
    vector_int children(nchild*active_nodes.size());

    // For each active node, generate the tag mask for each of its 2^D children
    detail::tabulate(children.begin(), children.end(),
                   child_index_to_tag_mask(level, max_level, active_nodes.data()));

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
  
    int length = (1 << (max_level - level) * 2) - 1;

    detail::upper_bound(m_tags.begin(),
                m_tags.end(),
                make_transform_iterator(children.begin(), detail::_1 + length),
                make_transform_iterator(children.end(), detail::_1 + length),
                upper_bounds.begin());

    /******************************************
     * 3. Mark each child as empty/leaf/node  *
     ******************************************/

    // Mark each child as either empty, a node, or a leaf
    vector_int child_node_kind(children.size(), 0);
    detail::transform(lower_bounds.begin(), lower_bounds.end(),
                    upper_bounds.begin(),
                    child_node_kind.begin(),
                    classify_node(threshold, level == max_level));


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

    // number of tree children (fixed from number of dimensions)
    const static int nchild = detail::ipow(2,dimension);

    // mask for lower n bits, where n is the number of dimensions
    const static unsigned mask = (1  << dimension) - 1;

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
    int threshold;
    int last_level;

    classify_node(int threshold, int last_level) : threshold(threshold), last_level(last_level) {}

    inline CUDA_HOST_DEVICE
        int operator()(int lower_bound, int upper_bound) const
        {
            int count = upper_bound - lower_bound;
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
                return num_nodes + 4 * node_idx;
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



}
#endif /* OCTTREE_H_ */
