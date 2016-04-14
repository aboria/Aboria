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
#include "OctTreeUtil.h"
#include "Particles.h"

namespace Aboria {


// Utility functions to encode leaves and children in single int
// are defined in util.h:
//   bool is_empty(int id);
//   bool is_node(int id);
//   bool is_leaf(int id);
//   int get_empty_id();
//   int get_leaf_id(int offset);
//   int get_leaf_offset(int id);
//   int child_tag_mask(int tag, int which_child, int level, int max_level);


template <typename traits>
class OctTree {
    typedef typename traits::template vector_type<Vect3d>::type vector_Vect3d;
    typedef typename traits::template vector_type<int>::type vector_int;
    typedef typename traits::template vector_type<Vect2i>::type vector_int2;

public:
    OctTree() {};

    void embed_points(vector_Vect3d& points);

    //void get_neighbours(const Vect3d &centre_point, const double radius, std::vector<range> &neighbours);

    //template<typename targets_traits, typename F>
    //void evaluate_kernel_fmm(targets_traits::iterator targets_begin, targets_traits::iterator targets_end, F &functor);

    void set_domain(Vect3d &min_in, Vect3d &max_in, Vect3b &periodic_in) {
        min = min_in;
        max = max_in;
        periodic = periodic_in;
    }
    void get_domain(Vect3d &min_out, Vect3d &max_out, Vect3b &periodic_out) {
        min_out = min;
        max_out = max;
        periodic_out = periodic;
    }
    void set_max_points(int arg) { max_points = arg; }
    void get_max_points(int arg) { return max_points; }
    void set_max_level(int arg) { max_level = arg; }
    void get_max_level(int arg) { return max_level; }
    void set_threshold(int arg) { threshold = arg; }
    void get_threshold(int arg) { return threshold; }

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
    Vect3d min,max;
    Vect3b periodic;

    vector_int tags;
    vector_int nodes;
    vector_int indices;
    bbox bounds;
    vector_int2 &leaves;
};


template <typename traits>
void OctTree<traits>::embed_points(vector_Vect3d& points) {

    const int num_points = points.size();

    nodes.clear();
    leaves.clear();
    tags.resize(num_points);
    indices.resize(num_points);

    /******************************************
     * 2. Compute bounding box                *
     ******************************************/

    bounds = reduce(points.begin(), points.end(), bbox(), plus<bbox>());

    /******************************************
     * 3. Classify points                     *
     ******************************************/

    transform(points.begin(), 
            points.end(), 
            tags.begin(), 
            classify_point(bounds, max_level));

    /******************************************
     * 4. Sort according to classification    *
     ******************************************/


    // Now that we have the geometric information, we can sort the
    // points accordingly.
    sequence(indices.begin(), indices.end());
    sort_by_key(tags.begin(), tags.end(), indices.begin());
}



template <typename traits>
void OctTree<traits>::build_tree() {
  vector_int active_nodes(1,0);

  // Build the tree one level at a time, starting at the root
  for (int level = 1 ; !active_nodes.empty() && level <= max_level ; ++level)
  {
    /******************************************
     * 1. Calculate children                  *
     ******************************************/

    // New children: 4 quadrants per active node = 4 children
    vector_int children(4*active_nodes.size());

    // For each active node, generate the tag mask for each of its 4 children
    tabulate(children.begin(), children.end(),
                   child_index_to_tag_mask(level, max_level, active_nodes.data()));

    /******************************************
     * 2. Determine interval for each child   *
     ******************************************/

    // For each child we need interval bounds
    vector_int lower_bounds(children.size());
    vector_int upper_bounds(children.size());

    // Locate lower and upper bounds for points in each quadrant
    lower_bound(tags.begin(),
                      tags.end(),
                      children.begin(),
                      children.end(),
                      lower_bounds.begin());
  
    int length = (1 << (max_level - level) * 2) - 1;

    upper_bound(tags.begin(),
                      tags.end(),
                      make_transform_iterator(children.begin(), _1 + length),
                      make_transform_iterator(children.end(), _1 + length),
                      upper_bounds.begin());


    /******************************************
     * 3. Mark each child as empty/leaf/node  *
     ******************************************/

    // Mark each child as either empty, a node, or a leaf
    vector_int child_node_kind(children.size(), 0);
    transform(lower_bounds.begin(), lower_bounds.end(),
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
    transform_exclusive_scan(child_node_kind.begin(), 
                                   child_node_kind.end(), 
                                   nodes_on_this_level.begin(), 
                                   is_a<NODE>(), 
                                   0, 
                                   plus<int>());
  
    // Enumerate leaves at this level
    transform_exclusive_scan(child_node_kind.begin(), 
                                   child_node_kind.end(), 
                                   leaves_on_this_level.begin(), 
                                   is_a<LEAF>(), 
                                   0, 
                                   plus<int>());

    int num_nodes_on_this_level = nodes_on_this_level.back() + (child_node_kind.back() == NODE ? 1 : 0);
    int num_leaves_on_this_level = leaves_on_this_level.back() + (child_node_kind.back() == LEAF ? 1 : 0);


    /******************************************
     * 5. Add the children to the node list   *
     ******************************************/

    int num_children = child_node_kind.size();

    int children_begin = nodes.size();
    nodes.resize(nodes.size() + num_children);
      
    transform(
            make_zip_iterator(
                make_tuple(child_node_kind.begin(), nodes_on_this_level.begin(), leaves_on_this_level.begin())),
             make_zip_iterator(
                make_tuple(child_node_kind.end(), nodes_on_this_level.end(), leaves_on_this_level.end())),
                nodes.begin() + children_begin,
                write_nodes(nodes.size(), leaves.size()));


    /******************************************
     * 6. Add the leaves to the leaf list     *
     ******************************************/

    children_begin = leaves.size();

    leaves.resize(leaves.size() + num_leaves_on_this_level);

    scatter_if(make_transform_iterator(
                         make_zip_iterator(
                             make_tuple(lower_bounds.begin(), upper_bounds.begin())),
                         make_leaf()),
                     make_transform_iterator(
                         make_zip_iterator(
                             make_tuple(lower_bounds.end(), upper_bounds.end())),
                         make_leaf()),
                     leaves_on_this_level.begin(),
                     child_node_kind.begin(),
                     leaves.begin() + children_begin,
                     is_a<LEAF>());


    /******************************************
     * 7. Set the nodes for the next level    *
     ******************************************/

    // Set active nodes for the next level to be all the childs nodes from this level
    active_nodes.resize(num_nodes_on_this_level);
  
    copy_if(children.begin(),
                  children.end(),
                  child_node_kind.begin(),
                  active_nodes.begin(),
                  is_a<NODE>());

  }
}

// Classify a point with respect to the bounding box.
template <typename traits>
struct OctTree<traits>::classify_point {
    bbox box;
    int max_level;

    // Create the classifier
    classify_point(const bbox &b, int lvl) : box(b), max_level(lvl) {}

    // Classify a point
    inline CUDA_HOST_DEVICE
        int operator()(const Vect3d &p) { return point_to_tag(p, box, max_level); }
};


template <typename traits>
struct OctTree<traits>::child_index_to_tag_mask {
    int level, max_level;
    typedef typename vector_int::const_pointer ptr_type;
    ptr_type nodes;

    child_index_to_tag_mask(int lvl, int max_lvl, ptr_type nodes) : level(lvl), max_level(max_lvl), nodes(nodes) {}

    inline CUDA_HOST_DEVICE
        int operator()(int idx) const
        {
            int tag = nodes[idx/4];
            int which_child = (idx&3);
            return child_tag_mask(tag, which_child, level, max_level);
        }
};


template <typename traits>
struct OctTree<traits>::classify_node
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
                return EMPTY;
            }
            else if (last_level || count < threshold)
            {
                return LEAF;
            }
            else
            {
                return NODE;
            }
        }
};

template <typename traits>
struct OctTree<traits>::write_nodes {
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

            if (node_type == EMPTY)
            {
                return get_empty_id();
            }
            else if (node_type == LEAF)
            {
                return get_leaf_id(num_leaves + leaf_idx);
            }
            else
            {
                return num_nodes + 4 * node_idx;
            }
        }
};

template <typename traits>
struct OctTree<traits>::make_leaf {
    typedef Vect2i result_type;
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
