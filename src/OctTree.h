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

// Operator which merges two bounding boxes.
struct merge_bboxes
{
  inline CUDA_HOST_DEVICE
  bbox operator()(const bbox &b0, const bbox &b1) const
  {
    bbox bounds;
    bounds.xmin = min(b0.xmin, b1.xmin);
    bounds.xmax = max(b0.xmax, b1.xmax);
    bounds.ymin = min(b0.ymin, b1.ymin);
    bounds.ymax = max(b0.ymax, b1.ymax);
    return bounds;
  }
};

bbox compute_bounding_box(const thrust::device_vector<float2> &points)
{
  return thrust::reduce(points.begin(), points.end(), bbox(), merge_bboxes());
}


// Classify a point with respect to the bounding box.
struct classify_point
{
  bbox box;
  int max_level;

  // Create the classifier
  classify_point(const bbox &b, int lvl) : box(b), max_level(lvl) {}

  // Classify a point
  inline CUDA_HOST_DEVICE
  int operator()(const float2 &p) { return point_to_tag(p, box, max_level); }
};

void compute_tags(const thrust::device_vector<float2> &points, const bbox &bounds, int max_level, thrust::device_vector<int> &tags)
{
  thrust::transform(points.begin(), 
                    points.end(), 
                    tags.begin(), 
                    classify_point(bounds, max_level));
}


void sort_points_by_tag(thrust::device_vector<int> &tags, thrust::device_vector<int> &indices)
{
  thrust::sequence(indices.begin(), indices.end());
  thrust::sort_by_key(tags.begin(), tags.end(), indices.begin());
}


struct child_index_to_tag_mask
{
  int level, max_level;
  thrust::device_ptr<const int> nodes;
  
  child_index_to_tag_mask(int lvl, int max_lvl, thrust::device_ptr<const int> nodes) : level(lvl), max_level(max_lvl), nodes(nodes) {}
  
  inline CUDA_HOST_DEVICE
  int operator()(int idx) const
  {
    int tag = nodes[idx/4];
    int which_child = (idx&3);
    return child_tag_mask(tag, which_child, level, max_level);
  }
};


void compute_child_tag_masks(const thrust::device_vector<int> &active_nodes,
                             int level,
                             int max_level,
                             thrust::device_vector<int> &children)
{
  // For each active node, generate the tag mask for each of its 4 children
  thrust::tabulate(children.begin(), children.end(),
                   child_index_to_tag_mask(level, max_level, active_nodes.data()));
}


void find_child_bounds(const thrust::device_vector<int> &tags,
                       const thrust::device_vector<int> &children,
                       int level,
                       int max_level,
                       thrust::device_vector<int> &lower_bounds,
                       thrust::device_vector<int> &upper_bounds)
{
  // Locate lower and upper bounds for points in each quadrant
  thrust::lower_bound(tags.begin(),
                      tags.end(),
                      children.begin(),
                      children.end(),
                      lower_bounds.begin());
  
  int length = (1 << (max_level - level) * 2) - 1;

  using namespace thrust::placeholders;

  thrust::upper_bound(tags.begin(),
                      tags.end(),
                      thrust::make_transform_iterator(children.begin(), _1 + length),
                      thrust::make_transform_iterator(children.end(), _1 + length),
                      upper_bounds.begin());
}


struct classify_node
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


void classify_children(const thrust::device_vector<int> &lower_bounds,
                       const thrust::device_vector<int> &upper_bounds,
                       int level,
                       int max_level,
                       int threshold,
                       thrust::device_vector<int> &child_node_kind)
{
  thrust::transform(lower_bounds.begin(), lower_bounds.end(),
                    upper_bounds.begin(),
                    child_node_kind.begin(),
                    classify_node(threshold, level == max_level));
}


std::pair<int,int> enumerate_nodes_and_leaves(const thrust::device_vector<int> &child_node_kind,
                                              thrust::device_vector<int> &nodes_on_this_level,
                                              thrust::device_vector<int> &leaves_on_this_level)
{
  // Enumerate nodes at this level
  thrust::transform_exclusive_scan(child_node_kind.begin(), 
                                   child_node_kind.end(), 
                                   nodes_on_this_level.begin(), 
                                   is_a<NODE>(), 
                                   0, 
                                   thrust::plus<int>());
  
  // Enumerate leaves at this level
  thrust::transform_exclusive_scan(child_node_kind.begin(), 
                                   child_node_kind.end(), 
                                   leaves_on_this_level.begin(), 
                                   is_a<LEAF>(), 
                                   0, 
                                   thrust::plus<int>());

  std::pair<int,int> num_nodes_and_leaves_on_this_level;

  num_nodes_and_leaves_on_this_level.first = nodes_on_this_level.back() + (child_node_kind.back() == NODE ? 1 : 0);
  num_nodes_and_leaves_on_this_level.second = leaves_on_this_level.back() + (child_node_kind.back() == LEAF ? 1 : 0);

  return num_nodes_and_leaves_on_this_level;
}


struct write_nodes
{
  int num_nodes, num_leaves;

  write_nodes(int num_nodes, int num_leaves) : 
    num_nodes(num_nodes), num_leaves(num_leaves) 
  {}

  template <typename tuple_type>
  inline CUDA_HOST_DEVICE
  int operator()(const tuple_type &t) const
  {
    int node_type = thrust::get<0>(t);
    int node_idx  = thrust::get<1>(t);
    int leaf_idx  = thrust::get<2>(t);

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


void create_child_nodes(const thrust::device_vector<int> &child_node_kind,
                        const thrust::device_vector<int> &nodes_on_this_level,
                        const thrust::device_vector<int> &leaves_on_this_level,
                        int num_leaves,
                        thrust::device_vector<int> &nodes)
{
  int num_children = child_node_kind.size();

  int children_begin = nodes.size();
  nodes.resize(nodes.size() + num_children);
  
  thrust::transform(thrust::make_zip_iterator(
                        thrust::make_tuple(
                            child_node_kind.begin(), nodes_on_this_level.begin(), leaves_on_this_level.begin())),
                    thrust::make_zip_iterator(
                        thrust::make_tuple(
                            child_node_kind.end(), nodes_on_this_level.end(), leaves_on_this_level.end())),
                    nodes.begin() + children_begin,
                    write_nodes(nodes.size(), num_leaves));
}


struct make_leaf
{
  typedef int2 result_type;
  template <typename tuple_type>
  inline CUDA_HOST_DEVICE
  int2 operator()(const tuple_type &t) const
  {
    int x = thrust::get<0>(t);
    int y = thrust::get<1>(t);

    return make_int2(x, y);
  }
};


void create_leaves(const thrust::device_vector<int> &child_node_kind,
                   const thrust::device_vector<int> &leaves_on_this_level,
                   const thrust::device_vector<int> &lower_bounds,
                   const thrust::device_vector<int> &upper_bounds,
                   int num_leaves_on_this_level,
                   thrust::device_vector<int2> &leaves)
{
  int children_begin = leaves.size();

  leaves.resize(leaves.size() + num_leaves_on_this_level);

  thrust::scatter_if(thrust::make_transform_iterator(
                         thrust::make_zip_iterator(
                             thrust::make_tuple(lower_bounds.begin(), upper_bounds.begin())),
                         make_leaf()),
                     thrust::make_transform_iterator(
                         thrust::make_zip_iterator(
                             thrust::make_tuple(lower_bounds.end(), upper_bounds.end())),
                         make_leaf()),
                     leaves_on_this_level.begin(),
                     child_node_kind.begin(),
                     leaves.begin() + children_begin,
                     is_a<LEAF>());
}


void activate_nodes_for_next_level(const thrust::device_vector<int> &children,
                                   const thrust::device_vector<int> &child_node_kind,
                                   int num_nodes_on_this_level,
                                   thrust::device_vector<int> &active_nodes)
{
  // Set active nodes for the next level to be all the childs nodes from this level
  active_nodes.resize(num_nodes_on_this_level);
  
  thrust::copy_if(children.begin(),
                  children.end(),
                  child_node_kind.begin(),
                  active_nodes.begin(),
                  is_a<NODE>());
}

template <typename traits>
class OctTree {
    typedef traits::vector_Vect3d vector_Vect3d;
    typedef traits::vector_int vector_int;
    typedef traits::iterator iterator;
    typedef traits::range range;

public:
    OctTree() {};

    void embed_points(iterator begin, iterator end);

    void get_neighbours(const Vect3d &centre_point, const double radius, std::vector<range> &neighbours);

    template<typename targets_traits, typename F>
    void evaluate_kernel_fmm(targets_traits::iterator targets_begin, targets_traits::iterator targets_end, F &functor);

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
    void set_max_levels(int arg) { max_levels = arg; }
    void get_max_levels(int arg) { return max_levels; }

private:
    void build_tree(const vector_int &tags, vector_int &nodes, vector_int2 &leaves);

    iterator begin;
    iterator end;
    int max_points;
    int max_levels;
    Vect3d min,max;
    Vect3b periodic;
}



void OctTree::build_tree(const vector_int &tags, vector_int &nodes, vector_int2 &leaves) {
  thrust::device_vector<int> active_nodes(1,0);

  // Build the tree one level at a time, starting at the root
  for (int level = 1 ; !active_nodes.empty() && level <= max_level ; ++level)
  {
    /******************************************
     * 1. Calculate children                  *
     ******************************************/

    // New children: 4 quadrants per active node = 4 children
    thrust::device_vector<int> children(4*active_nodes.size());
    compute_child_tag_masks(active_nodes, level, max_level, children);

    /******************************************
     * 2. Determine interval for each child   *
     ******************************************/

    // For each child we need interval bounds
    thrust::device_vector<int> lower_bounds(children.size());
    thrust::device_vector<int> upper_bounds(children.size());
    find_child_bounds(tags, children, level, max_level, lower_bounds, upper_bounds);

    /******************************************
     * 3. Mark each child as empty/leaf/node  *
     ******************************************/

    // Mark each child as either empty, a node, or a leaf
    thrust::device_vector<int> child_node_kind(children.size(), 0);
    classify_children(lower_bounds, upper_bounds, level, max_level, threshold, child_node_kind);

    /******************************************
     * 4. Enumerate nodes and leaves          *
     ******************************************/

    // Enumerate the nodes and leaves at this level
    thrust::device_vector<int> leaves_on_this_level(child_node_kind.size());
    thrust::device_vector<int> nodes_on_this_level(child_node_kind.size());

    // Enumerate nodes and leaves at this level
    std::pair<int,int> num_nodes_and_leaves_on_this_level =
      enumerate_nodes_and_leaves(child_node_kind, nodes_on_this_level, leaves_on_this_level);

    /******************************************
     * 5. Add the children to the node list   *
     ******************************************/

    create_child_nodes(child_node_kind, nodes_on_this_level, leaves_on_this_level, leaves.size(), nodes);

    /******************************************
     * 6. Add the leaves to the leaf list     *
     ******************************************/

    create_leaves(child_node_kind, leaves_on_this_level, lower_bounds, upper_bounds, num_nodes_and_leaves_on_this_level.second, leaves);

    /******************************************
     * 7. Set the nodes for the next level    *
     ******************************************/

    activate_nodes_for_next_level(children, child_node_kind, num_nodes_and_leaves_on_this_level.first, active_nodes);
  }
}

void run_experiment(thrust::device_vector<float2> *points,
                    thrust::device_vector<int> *nodes,
                    thrust::device_vector<int2> *leaves,
                    const int threshold,
                    const int max_level)
{
  const int num_points = points->size();
  /******************************************
   * 1. Generate points                     *
   ******************************************/

  generate_random_points(*points);

  /******************************************
   * 2. Compute bounding box                *
   ******************************************/

  bbox bounds = compute_bounding_box(*points);

  /******************************************
   * 3. Classify points                     *
   ******************************************/

  thrust::device_vector<int> tags(num_points);
  
  compute_tags(*points, bounds, max_level, tags);

  /******************************************
   * 4. Sort according to classification    *
   ******************************************/

  thrust::device_vector<int> indices(num_points);

  // Now that we have the geometric information, we can sort the
  // points accordingly.
  sort_points_by_tag(tags, indices);

  /******************************************
   * 5. Build the tree                      *
   ******************************************/

  build_tree(tags, bounds, max_level, threshold, *nodes, *leaves);
}

int main()
{
  const int num_points = 4*1024*1024;
  const int threshold = 32; // A node with fewer than threshold points is a leaf.
  const int max_level = 10;

  thrust::device_vector<float2> points(num_points);
  thrust::device_vector<int> nodes;
  thrust::device_vector<int2> leaves;

  std::cout << "Warming up...\n" << std::endl;

  // validate and warm up the JIT
  run_experiment(&points, &nodes, &leaves, threshold, max_level);

  std::cout << "Timing...\n" << std::endl;

  int num_trials = 25;
  double mean_msecs = time_invocation_cuda(num_trials, run_experiment, &points, &nodes, &leaves, threshold, max_level);
  double mean_secs = mean_msecs / 1000;
  double millions_of_points = double(num_points) / 1000000;

  std::cout << millions_of_points / mean_secs << " millions of points generated and treeified per second." << std::endl;
}


}
#endif /* OCTTREE_H_ */
