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

#ifdef HAVE_CAIRO
#include <cairo-svg.h>
#endif

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
    for (std::vector<double>::iterator i = v.begin(); i != v.end(); ++i) {
      *i = index++;
    }

    /*`
    Or in more compact notation
    */

    index = 0;
    for (auto i = v.begin(); i != v.end(); ++i) {
      *i = index++;
    }

    /*`
    The iterators in Aboria are similar to STL iterators in that they can step
    through a given range of objects using the `++` operator, with some slight
    differences that we will describe below.


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
    classes (e.g. [classref Aboria::NeighbourQueryBase]) were unsuitable for
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

    [note All the spatial data structures in Aboria are considered trees.
    For the HyperOctree and kdtree data structures, this description is obvious,
    but the cell list is also treated as a tree, in this case a tree with one
    root node having N children, where N is the total number of buckets in the
    cell list.]

    [note Aboria tends to use the terms nodes and buckets fairly interchangably.
    ]

    You can create a child_iterator by using the [memberref
    Aboria::NeighbourQueryBase::get_children] function. This creates a
    child_iterator that loops through the children of the root node of the tree.
    We will use this iterator to loop through all these children and print out
    the spatial bounds of each node. In order to determine when we have reached
    the end of the children, we can compare the iterator to `false`. This
    pattern is widely used in Aboria, rather than specific `end` iterators as
    used in the STL.

    */

    for (auto i = query.get_children(); i != false; ++i) {
      std::cout << query.get_bounds(i) << std::endl;
    }

    /*`
    Above we use the [memberref Aboria::NeighbourQueryBase::get_bounds] function
    to get the bounds of the child iterator. This returns a [classref
    Aboria::bbox] class that contains the minimum and
    maximum spatial extents of the node pointed to by `i`.

    Note that here we are using the default spatial data structure, a cell list
    provided by [classref Aboria::CellList], so the "tree" here will only have 2
    levels, and the loop above will loop through the second (i.e. non-root)
    level. We can also create a proper tree structure using the hyper oct-tree
    data structure given by [classref Aboria::HyperOctree], like so:
    */

    typedef Particles<std::tuple<>, 2, std::vector, HyperOctree>
        ParticlesOcttree_t;
    ParticlesOcttree_t particles_octtree(N);

    for (size_t i = 0; i < N; ++i) {
      auto &gen = get<generator>(particles_octtree)[i];
      get<position>(particles_octtree)[i] =
          vdouble2(uniform(gen), uniform(gen));
    }

    particles_octtree.init_neighbour_search(vdouble2(0, 0), vdouble2(1, 1),
                                            vdouble2(false, false));

    auto query_octtree = particles_octtree.get_query();

    /*`
    Now `particles_octtree` contains a full oct-tree, dividing the spatial
    domain into a hierarchical set of boxes that make up our tree data
    structure. The simplest iteration we might want to do on the tree is a
    depth-first iteration, which is easiest achieved by recursion. The
    [memberref Aboria::NeighbourQueryBase::get_children] function can be used to
    get the children of a [classref Aboria::NeighbourQueryBase::child_iterator],
    and using a C++ lambda function to provide the recursion we can implement a
    depth-first iteration like so
    */

    std::cout << "recursive depth-first" << std::endl;
    for (auto i = query_octtree.get_children(); i != false; ++i) {
      std::function<void(decltype(i) &)> depth_first;
      depth_first = [&](const auto &parent) {
        std::cout << query_octtree.get_bounds(parent) << std::endl;
        for (auto i = query_octtree.get_children(parent); i != false; ++i) {
          depth_first(i);
        }
      };
      depth_first(i);
    }

    /*`
    This construction might be a bit clumsy to use in practice however, so
    Aboria provides a special depth-first iterator
    [classref Aboria::NeighbourQueryBase::all_iterator] to allow you to write a
    loop equivalent to the recursive depth-first code given above.

    The [memberref Aboria::NeighbourQueryBase::get_subtree] function returns a
    [classref Aboria::NeighbourQueryBase::all_iterator] that performs a
    depth-first iteration over the tree. Note that you can also pass in a
    child_iterator to [memberref Aboria::NeighbourQueryBase::get_subtree] to
    iterate over the sub-tree below a particular node of the tree.
    */

    std::cout << "subtree depth-first" << std::endl;
    for (auto i = query_octtree.get_subtree(); i != false; ++i) {
      std::cout << query_octtree.get_bounds(i.get_child_iterator())
                << std::endl;
    }

    /*`
    You might also want to distinguish between leaf nodes (nodes with no
    children) and non-leaf nodes. You can do this with the
    [memberref Aboria::NeighbourQueryBase::is_leaf_node] function, which takes a
    reference to a node (rather than an iterator), and can be used like so
    */

    std::cout << "subtree depth-first showing leaf nodes" << std::endl;
    for (auto i = query_octtree.get_subtree(); i != false; ++i) {
      auto ci = i.get_child_iterator();
      if (query_octtree.is_leaf_node(*ci)) {
        std::cout << "leaf node with bounds = " << query_octtree.get_bounds(ci)
                  << std::endl;
      } else {
        std::cout << "non-leaf node with bounds = "
                  << query_octtree.get_bounds(ci) << std::endl;
      }
    }

    /*`
    Leaf nodes in the tree are the only nodes that contain particles. You can
    loop through all the particles in a given leaf node using the [memberref
    Aboria::NeighbourQueryBase::get_bucket_particles] function, which returns
    an iterator. Note for non-leaf nodes, the [memberref
    Aboria::NeighbourQueryBase::get_bucket_particles] will return an iterator
    that is immediatelly false, so this loop is safe even for non-leaf nodes.
    */

    std::cout << "subtree depth-first showing leaf nodes and particles"
              << std::endl;
    for (auto i = query_octtree.get_subtree(); i != false; ++i) {
      auto ci = i.get_child_iterator();
      if (query_octtree.is_leaf_node(*ci)) {
        std::cout << "leaf node with bounds = " << query_octtree.get_bounds(ci)
                  << std::endl;
        for (auto j = query_octtree.get_bucket_particles(*ci); j != false;
             ++j) {
          std::cout << "\t has particle with position" << get<position>(*j)
                    << std::endl;
        }
      } else {
        std::cout << "non-leaf node with bounds = "
                  << query_octtree.get_bounds(ci) << std::endl;
      }
    }

    /*`
    Aboria also provides functions to query leaf nodes, or buckets, within a
    certain distance of a point, and these are used internally for the neighbour
    search functionality discussed in earlier sections. You can use the
    [memberref Aboria::NeighbourQueryBase::get_buckets_near_point] function,
    which returns a [classref Aboria::NeighbourQueryBase::query_iterator] of all
    the buckets with a given distance of a point. This function also takes a
    template argument `P`, which refers to the p-norm distance that it uses
    (i.e. P=2 is the standard euclidean distance).

    [caution The distance search provided by [memberref
    Aboria::NeighbourQueryBase::get_buckets_near_point] does not respect the
    periodicity of the domain, so if you did a search near a lhs edge of a
    periodic domain, it would not pick up buckets on the neighbouring periodic
    rhs edge.]


    */

    const int P = 2;
    const vdouble2 search_point = vdouble2(0.5, 0.5);
    const double search_radius = 0.1;
    std::cout << "searching within " << search_point << " of point "
              << search_point << std::endl;
    for (auto i = query_octtree.get_buckets_near_point<P>(search_point,
                                                          search_radius);
         i != false; ++i) {
      auto ci = i.get_child_iterator();
      std::cout << "\t found bucket at " << query_octtree.get_bounds(ci)
                << std::endl;
      for (auto j = query_octtree.get_bucket_particles(*ci); j != false; ++j) {
        std::cout << "\t\t found particle at " << get<position>(*j)
                  << std::endl;
      }
    }

    /*`
    [endsect]
    */
    //]
  }

  template <template <typename, typename> class Vector,
            template <typename> class SearchMethod>
  void helper_data_structure() {

    const size_t N = 20;
    typedef Particles<std::tuple<>, 2, Vector, SearchMethod> Particles_t;
    typedef typename Particles_t::position position;
    Particles_t particles(N);

    std::uniform_real_distribution<double> uniform(0, 1);
    for (size_t i = 0; i < N; ++i) {
      auto &gen = get<generator>(particles)[i];
      get<position>(particles)[i] = vdouble2(uniform(gen), uniform(gen));
    }

    particles.init_neighbour_search(vdouble2(0, 0), vdouble2(1, 1),
                                    vdouble2(false, false), 5);

    auto query = particles.get_query();

    int child_count_recurse = 0;
    for (auto i = query.get_children(); i != false; ++i) {
      std::function<int(decltype(i) &)> depth_first;
      depth_first = [&](const auto &parent) {
        int child_count = 1;
        auto parent_bounds = query.get_bounds(parent);
        for (auto i = query.get_children(parent); i != false; ++i) {
          child_count += depth_first(i);
          auto bounds = query.get_bounds(i);
          for (size_t i = 0; i < 2; ++i) {
            TS_ASSERT_LESS_THAN_EQUALS(bounds.bmax[i], parent_bounds.bmax[i]);
            TS_ASSERT_LESS_THAN_EQUALS(parent_bounds.bmin[i], bounds.bmin[i]);
          }
        }
        if (query.is_leaf_node(*parent)) {
          TS_ASSERT_EQUALS(child_count, 1);
          for (auto j = query.get_bucket_particles(*i); j != false; ++j) {
            auto p = get<position>(*j);
            for (size_t i = 0; i < 2; ++i) {
              TS_ASSERT_LESS_THAN_EQUALS(p[i], parent_bounds.bmax[i]);
              TS_ASSERT_LESS_THAN_EQUALS(parent_bounds.bmin[i], p[i]);
            }
          }
        }
        return child_count;
      };
      child_count_recurse += depth_first(i);
    }

    int child_count_subtree = 0;
    for (auto i = query.get_subtree(); i != false; ++i) {
      child_count_subtree++;
      auto bounds = query.get_bounds(i.get_child_iterator());
      int num_particles = 0;
      for (auto j = query.get_bucket_particles(*i); j != false; ++j) {
        num_particles++;
        auto p = get<position>(*j);
        for (size_t i = 0; i < 2; ++i) {
          TS_ASSERT_LESS_THAN_EQUALS(p[i], bounds.bmax[i]);
          TS_ASSERT_LESS_THAN_EQUALS(bounds.bmin[i], p[i]);
        }
      }
      if (!query.is_leaf_node(*i)) {
        TS_ASSERT_EQUALS(num_particles, 0);
      }
    }

    TS_ASSERT_EQUALS(child_count_recurse, child_count_subtree);
  }

  template <template <typename> class SearchMethod> void draw_data_structure() {
    using Particles_t = Particles<std::tuple<>, 2, std::vector, SearchMethod>;
    Particles_t particles(500);
    using position = typename Particles_t::position;
    std::normal_distribution<double> normal(0.5, 0.2);
    for (size_t i = 0; i < particles.size(); ++i) {
      auto &gen = get<generator>(particles)[i];
      get<position>(particles)[i] = vdouble2(normal(gen), normal(gen));
    }
#ifdef HAVE_CAIRO
    std::string search_name(typeid(typename Particles_t::search_type).name());
    std::string extension(".svg");
    const int image_size = 512;
    cairo_surface_t *surface = cairo_svg_surface_create(
        (search_name + extension).c_str(), image_size, image_size);
    cairo_svg_surface_restrict_to_version(surface, CAIRO_SVG_VERSION_1_2);
    cairo_t *cr = cairo_create(surface);

    cairo_scale(cr, image_size, image_size);
    particles.init_neighbour_search(vdouble2::Constant(0),
                                    vdouble2::Constant(1),
                                    vdouble2::Constant(false), 5);

    auto query = particles.get_query();

    cairo_set_source_rgba(cr, 0.5, 0, 0, 0.5);
    vdouble2 lw(0.005, 0.005);
    // cairo_device_to_user_distance(cr, &lw[0], &lw[1]);
    cairo_set_line_width(cr, lw[0]);
    for (auto i = query.get_subtree(); i != false; ++i) {
      auto ci = i.get_child_iterator();
      if (query.is_leaf_node(*ci)) {
        // draw its outline
        auto bounds = query.get_bounds(ci);
        cairo_move_to(cr, bounds.bmin[0], bounds.bmin[1]);
        cairo_line_to(cr, bounds.bmax[0], bounds.bmin[1]);
        cairo_line_to(cr, bounds.bmax[0], bounds.bmax[1]);
        cairo_line_to(cr, bounds.bmin[0], bounds.bmax[1]);
        cairo_close_path(cr);
        cairo_stroke(cr);
      }
    }

    cairo_set_source_rgba(cr, 0, 0, 0, 1.0);
    const double PI = boost::math::constants::pi<double>();
    for (size_t i = 0; i < particles.size(); ++i) {
      const auto &pos = get<position>(particles)[i];
      cairo_arc(cr, pos[0], pos[1], lw[0], 0, 2 * PI);
      cairo_fill(cr);
    }

    const vdouble2 search_point = vdouble2(0.75, 0.57);
    const double search_radius = 0.2;

    // draw search region and point
    cairo_set_source_rgba(cr, 0, 0, 0.5, 0.5);
    cairo_arc(cr, search_point[0], search_point[1], search_radius, 0, 2 * PI);
    cairo_stroke(cr);
    cairo_arc(cr, search_point[0], search_point[1], lw[0], 0, 2 * PI);
    cairo_fill(cr);

    cairo_set_source_rgba(cr, 0, 0, 0.5, 0.2);
    for (auto i = query.template get_buckets_near_point<2>(search_point,
                                                           search_radius);
         i != false; ++i) {
      auto ci = i.get_child_iterator();
      auto bounds = query.get_bounds(ci);
      // colour in search buckets
      cairo_move_to(cr, bounds.bmin[0], bounds.bmin[1]);
      cairo_line_to(cr, bounds.bmax[0], bounds.bmin[1]);
      cairo_line_to(cr, bounds.bmax[0], bounds.bmax[1]);
      cairo_line_to(cr, bounds.bmin[0], bounds.bmax[1]);
      cairo_close_path(cr);
      cairo_fill(cr);
      for (auto j = query.get_bucket_particles(*ci); j != false; ++j) {
        const auto &pos = get<position>(*j);
        cairo_arc(cr, pos[0], pos[1], lw[0], 0, 2 * PI);
        cairo_fill(cr);
      }
    }

    cairo_destroy(cr);
    cairo_surface_destroy(surface);
#endif // HAVE_CAIRO
  }

  void test_visualise_data_structures() {
    draw_data_structure<CellList>();
    draw_data_structure<CellListOrdered>();
    draw_data_structure<Kdtree>();
    draw_data_structure<KdtreeNanoflann>();
    draw_data_structure<HyperOctree>();
  }

  void test_CellList() {
    std::cout << "CellList" << std::endl;
    helper_data_structure<std::vector, CellList>();
  }
  void test_CellListOrdered() {
    std::cout << "CellListOrdered" << std::endl;
    helper_data_structure<std::vector, CellListOrdered>();
  }
  void test_HyperOctree() {
    std::cout << "Octtree" << std::endl;
    helper_data_structure<std::vector, HyperOctree>();
  }
  void test_kdtree() {
    std::cout << "kd tree" << std::endl;
    helper_data_structure<std::vector, Kdtree>();
  }
};

#endif /* SPATIAL_DATA_STRUCTURES_H_ */
