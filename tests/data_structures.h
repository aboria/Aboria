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

#ifndef SPATIAL_DATA_STRUCTURES_H_
#define SPATIAL_DATA_STRUCTURES_H_

#include <chrono>
#include <cxxtest/TestSuite.h>
typedef std::chrono::system_clock Clock;

#include "Aboria.h"

using namespace Aboria;

class SpatialDataStructuresTest : public CxxTest::TestSuite {
public:
  void test_documentation(void) {
    //[spatial_data_structures
    /*`
    [section Using the spatial data structures]

    The neighbour-searching functionality within Aboria uses the underlying
    spatial data structure (i.e. cell list, octtree, kdtree) to perform its
    task. You can directly access this underlying data structure as well, in
    order for you to write your own spatial algorithms.

    Aboria uses iterators to provide a generic way to interact with the spatial
    data structure, so that the code that you write is independent of the
    particlular data structure that you use. An iterator is a generic C++ object
    that "iterates" through a 1D container. For example, the classic iterator is
    the loop index `i` in the code below.
    */

    size_t N = 10;
    std::vector<double> v(N);
    for (size_t i = 0; i < N; ++i) {
      v[i] = i;
    }

    /*`
    The STL library abstracts away the particular type of `i`, and defines a set
    of iterators for each container. For example, the `std::vector` class has
    its own iterator, which you can use as follows.
    */

    size_t index = 0;
    for (std::vector<double>::iterator i = v.begin(), i != v.end(), ++i) {
      *i = index++;
    }

    /*`
    Or in more compact notation
    */

    index = 0;
    for (auto i = v.begin(), i != v.end(), ++i) {
      *i = index++;
    }

    /*`
    Before we can start using the data structure iterators in Aboria, we need a
    particle set. Lets create a random set of points in 2D:
    */

    N = 100;
    typedef Particles<std::tuple<>, 2> Particles_t;
    typedef typename Particles_t::position position;
    Particles_t particles(N);

    std::uniform_real_distribution<double> uniform(0, 1);
    for (size_t i = 0; i < N; ++i) {
      auto &gen = get<generator>(particles)[i];
      get<position>(particles)[i] = vdouble2(uniform(gen), uniform(gen));
    }

    particles.init_neighbour_search(vdouble2(0, 0), vdouble2(1, 1),
                                    vdouble2(false, false));

    /*`
    In order to start interacting with the spatial data structures, we need to
    get its query object from the particle set. A query object is a lightweight
    object that has all the information neccessary to access either the particle
    set itself or the underlying spatial data structures. It was designed
    because the [classref Aboria::Particles] container and the neighbour search
    classes (e.g. [classref Aboria::NeighbourSearchBase]) were unsuitable for
    copying to a gpu in order to perform calculations there, so a simpler class,
    the query class, was created with the required functionality.
    */

    typedef Particles_t::query_type Query_t;
    Query_t query = particles.get_query();

    /*`
    The base class for all the query objects is [classref
    Aboria::NeighbourQueryBase], and all the query classes for the individual
    data structures are derived from this. Now that we have `query`, we can
    create a [classref Aboria::NeighbourQueryBase::child_iterator]. This is the
    lowest level data structure iterator, and allows you to iterate through a
    set of child nodes attached to a single parent node within a tree structure.
    Note that all the spatial data structures in Aboria are considered trees.
    For the octtree and kdtree data structures, this description is obvious, but
    the cell list is also treated as a tree, in this case a tree with one root
    node having N children, where N is the total number of buckets in the cell
    list.

    You can create a child_iterator by using the [funcref
    Aboria::NeighbourQueryBase::get_children()] function. This creates a
    child_iterator that loops through the children of the root node of the tree.
    We will use this iterator to loop through all these children and print out
    the spatial bounds of each node.

    */

    for (auto i = query.get_children(); i != false; ++i) {
      std::cout << query.get_bounds(i) << std::endl;
    }

    /*`
    Above we use the [funcref Aboria::NeighbourQueryBase::get_bounds(const
    child_iterator& ci) const] function to get the bounds of the child iterator.
    This returns a [classref Aboria::detail::bbox] class that contains the
    minimum and maximum spatial extents of the node pointed to by `i`.

    Note that here we are using the default spatial data structure, a cell list
    provided by [classref Aboria::CellList], so the "tree" here will only have 2
    levels, and the loop above will loop through the second (i.e. non-root)
    level. We can also create a proper tree structure using the hyper oct-tree
    data structure given by [classref Aboria::octtree], like so:
    */

    typedef Particles<std::tuple<>, 2, std::vector, octtree> ParticlesOcttree_t;
    ParticlesOcttree_t particles_octtree(N);

    for (size_t i = 0; i < N; ++i) {
      auto &gen = get<generator>(particles_octtree)[i];
      get<position>(particles_octtree)[i] =
          vdouble2(uniform(gen), uniform(gen));
    }

    particles_octtree.init_neighbour_search(vdouble2(0, 0), vdouble2(1, 1),
                                            vdouble2(false, false));

    /*`
    Now `particles_octtree` contains a full oct-tree, dividing the spatial
    domain into a hierarchical set of boxes that make up our tree data
    structure. The simplest iteration we might want to do on the tree is a
    depth-first iteration, which is easiest achieved by iteration. The [classref
    Aboria::NeighbourQueryBase::get_children(const child_iterator& ci) const]
    function can be used to get the children of a `child_iterator`, and using a
    C++ lambda function to provide the recursion we can implement a depth-first
    iteration like so
    */

    auto depth_first = [&](const auto &parent) {
      std::cout << query.get_bounds(parent) << std::endl;
      for (auto i = query.get_children(parent); i != false; ++i) {
        depth_first(i);
      }
    };

    for (auto i = query.get_children(); i != false; ++i) {
      std::cout << query.get_bounds() << std::endl;
      depth_first(i);
    }

    /*`

    This construction might be a bit clumsy to use in practice however, so
    Aboria provides a special depth-first iterator to allow you to write
    a loop equivalent to the recursive depth-first code given above.

    */

    for (auto i = query.get_subtree(); i != false; ++i) {
      std::cout << query.get_bounds() << std::endl;
    }

    /*`
    The [funcref Aboria::NeighbourQueryBase::get_subtree()] function returns a
    [classref Aboria::NeighbourQueryBase::all_iterator] that performs a
    depth-first iteration over the tree. Note that you can also pass in a
    child_iterator to [funcref Aboria::NeighbourQueryBase::get_subtree(const
    child_iterator& ci) const] to iterate over the sub-tree below a particular
    node of the tree.

    You might also want to distinguish between leaf nodes (nodes with no
    children) and non-leaf nodes. You can do this with the [funcref
    Aboria::NeighbourQueryBase:is_leaf()] function, which can be used like so
    */

    for (auto i = query.get_subtree(); i != false; ++i) {
      if (query.is_leaf_node(*i)) {
        std::cout << "leaf node with bounds = " << query.get_bounds()
                  << std::endl;
      } else {
        std::cout << "non-leaf node with bounds = " << query.get_bounds()
                  << std::endl;
      }
    }

    /*`
    Leaf nodes in the tree are the only nodes that contain particles. You can
    loop through all the particles in a given leaf node using the [funcref
    Aboria::NeighbourQueryBase::get_bucket_particles()] function. Note that this
    function returns an [classref Aboria::iterator_range], a lightweight object
    containing  `begin()` and `end()` functions that return iterators over a
    given range. These can be used in C++ range-based loops.
    */

    for (auto i = query.get_subtree(); i != false; ++i) {
      if (query.is_leaf_node(*i)) {
        std::cout << "leaf node with bounds = " << query.get_bounds()
                  << std::endl;
        for (auto &p : get_bucket_particles(*i)) {
          std::cout << "\t has particle with position" << get<position>(p);
        }
      } else {
        std::cout << "non-leaf node with bounds = " << query.get_bounds()
                  << std::endl;
      }
    }

    /*`

    Neighbour searching is a key functionality area for Aboria, so
    */

    /*`
    [endsect]
    */
    //]
  }
};

#endif /* SPATIAL_DATA_STRUCTURES_H_ */
