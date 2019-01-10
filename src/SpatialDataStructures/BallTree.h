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

#ifndef BALL_TREE_H_
#define BALL_TREE_H_

#include "Get.h"
#include "Iterators/BalltreeIterators.h"
#include "Log.h"
#include "NeighbourSearchBase.h"
#include "SpatialUtil.h"
#include "Traits.h"
#include "Vector.h"
#include "miniball/cpp/main/Seb.h"
#include <boost/iterator/iterator_facade.hpp>
#include <iostream>
#include <set>
#include <vector>

namespace Aboria {

template <typename Traits> struct BalltreeQuery;

/// \brief KdTree
///
template <typename Traits>
class Balltree : public neighbour_search_base<Balltree<Traits>, Traits,
                                              BalltreeQuery<Traits>> {

  typedef typename Traits::double_d double_d;
  typedef typename Traits::position position;
  typedef typename Traits::vector_int vector_int;
  typedef typename Traits::vector_int2 vector_int2;
  typedef typename Traits::vector_bool2 vector_bool2;

  typedef typename Traits::template vector_type<vint3>::type vector_int3;

  typedef typename Traits::vector_double vector_double;
  typedef typename Traits::vector_double_d vector_double_d;
  typedef typename Traits::iterator iterator;
  typedef typename Traits::unsigned_int_d unsigned_int_d;
  static const unsigned int dimension = Traits::dimension;

  typedef neighbour_search_base<Balltree<Traits>, Traits, BalltreeQuery<Traits>>
      base_type;
  friend base_type;

public:
  Balltree() : base_type(), m_number_of_levels(0) {

    this->m_query.m_nodes_first_child =
        iterator_to_raw_pointer(m_nodes_first_child.begin());
    this->m_query.m_nodes_centre =
        iterator_to_raw_pointer(m_nodes_centre.begin());
    this->m_query.m_nodes_radius =
        iterator_to_raw_pointer(m_nodes_radius.begin());
    this->m_query.m_nodes_particles =
        iterator_to_raw_pointer(m_nodes_particles.begin());
    this->m_query.m_number_of_buckets = m_nodes_first_child.size();
    this->m_query.m_number_of_levels = m_number_of_levels;
  }

  Balltree(const Balltree &other)
      : base_type(other), m_nodes_first_child(other.m_nodes_first_child),
        m_nodes_particles(other.m_nodes_particles),
        m_nodes_centre(other.m_nodes_centre),
        m_nodes_radius(other.m_nodes_radius),
        m_particle_indicies(other.m_particle_indicies),
        m_number_of_levels(other.m_number_of_levels), m_query(other.m_query) {

    this->m_query.m_nodes_first_child =
        iterator_to_raw_pointer(m_nodes_first_child.begin());
    this->m_query.m_nodes_centre =
        iterator_to_raw_pointer(m_nodes_centre.begin());
    this->m_query.m_nodes_radius =
        iterator_to_raw_pointer(m_nodes_radius.begin());
    this->m_query.m_nodes_particles =
        iterator_to_raw_pointer(m_nodes_particles.begin());
  }
  //~Balltree() {}

  static constexpr bool ordered() { return true; }

  void print_data_structure() const { print_tree(); }

private:
  void set_domain_impl() {
    this->m_query.m_bounds.bmin = this->m_bounds.bmin;
    this->m_query.m_bounds.bmax = this->m_bounds.bmax;
    this->m_query.m_periodic = this->m_periodic;
  }

  void update_iterator_impl() {}

  void print_tree() const {
    int first_index = 0;
    size_t i = first_index;
    while (i < m_nodes_first_child.size() && m_nodes_first_child[i] < 0)
      ++i;
    while (i < m_nodes_first_child.size()) {
      // scan level to find first node
      int level_size = m_nodes_first_child[i] - first_index;
      print_level(first_index, level_size);
      first_index += level_size;
      i = first_index;
      while (i < m_nodes_first_child.size() && m_nodes_first_child[i] < 0)
        ++i;
    }
    print_level(first_index, m_nodes_first_child.size() - first_index);
  }

  void print_level(int index, int level_size) const {
    for (int i = index; i < index + level_size; ++i) {
      if (m_nodes_first_child[i] >= 0) {
        std::cout << "|n(" << m_nodes_first_child[i] << ","
                  << m_nodes_particles[i][0] << "," << m_nodes_particles[i][1]
                  << ")";
      } else {
        std::cout << "|l(" << m_nodes_particles[i][0] << ","
                  << m_nodes_particles[i][1] << ")";
      }
    }
    std::cout << "|" << std::endl;
  }

  void update_positions_impl(iterator update_begin, iterator update_end,
                             const int new_n,
                             const bool call_set_domain = true) {
    ASSERT(update_begin == this->m_particles_begin &&
               update_end == this->m_particles_end,
           "error should be update all");

    const size_t num_points = this->m_alive_indices.size();

    // setup particles
    LOG(3, "update_positions_impl(balltree): setup particles");
    m_particle_indicies.resize(num_points);
    // copy particle indicies that are alive
    detail::copy(this->m_alive_indices.begin(), this->m_alive_indices.end(),
                 m_particle_indicies.begin());

    /*
       for (size_t i = 0; i < dimension; ++i) {
       std::cout << "dimension " << i << std::endl;
       for (size_t j = i * num_points; j < (i + 1) * num_points; ++j) {
       std::cout << "particle_indicies[" << j
       << "] = " << m_particle_indicies[j] << " ("
       << get<position>(
       this->m_particles_begin)[m_particle_indicies[j]]
       << ")" << std::endl;
       }
       }
       */

    // build the tree
    LOG(3, "update_positions_impl(balltree): build tree");
    build_tree();

    // copy sorted indicies from 1st dim back to m_alive_indicies
    LOG(3, "update_positions_impl(balltree): finished build tree");
    detail::copy(m_particle_indicies.begin(),
                 m_particle_indicies.begin() + num_points,
                 this->m_alive_indices.begin());

    this->m_query.m_nodes_first_child =
        iterator_to_raw_pointer(m_nodes_first_child.begin());
    this->m_query.m_nodes_centre =
        iterator_to_raw_pointer(m_nodes_centre.begin());
    this->m_query.m_nodes_radius =
        iterator_to_raw_pointer(m_nodes_radius.begin());
    this->m_query.m_nodes_particles =
        iterator_to_raw_pointer(m_nodes_particles.begin());
    this->m_query.m_number_of_buckets = m_nodes_first_child.size();
    this->m_query.m_number_of_levels = m_number_of_levels;
  }

  void update_alive_impl(iterator update_begin, iterator update_end) {
    LOG(3, "update_alive_impl(balltree): updating build tree");

    const int update_start_index = update_begin - this->m_particles_begin;

    // update particle indicies for top node
    for (int c = 0; c < 2; ++c) {
      const int new_index = m_nodes_particles[0][c] >= update_start_index
                                ? this->m_alive_sum[m_nodes_particles[0][c] -
                                                    update_start_index] +
                                      update_start_index
                                : m_nodes_particles[0][c];
      m_nodes_particles[0][c] = new_index;
    }

    // update particle indicies for all the rest
    // also update if node is leaf or not, and remove empty children
    for (auto i = m_query.breadth_first(); i != false; ++i) {
      detail::for_each(i->begin(), i->end(), [&](const auto &ci) {
        const int index = m_query.get_child_index(*ci);
        for (int c = 0; c < 2; ++c) {
          const int new_index =
              m_nodes_particles[index][c] >= update_start_index
                  ? this->m_alive_sum[m_nodes_particles[index][c] -
                                      update_start_index] +
                        update_start_index
                  : m_nodes_particles[index][c];
          m_nodes_particles[index][c] = new_index;
          if (!m_query.is_leaf_node(*ci)) {
            const int num_particles =
                m_nodes_particles[index][1] - m_nodes_particles[index][0];
            if (num_particles <= this->m_n_particles_in_leaf) {
              // convert to leaf if new number of particles is below threshold
              m_nodes_first_child[index] = -1;
            } else {
              // if any children empty erase them by copying the *other* child
              // to this node
              for (int c = 0; c < 2; ++c) {
                const int child_index = m_nodes_first_child[index] + c;
                if (m_nodes_particles[child_index][1] ==
                    m_nodes_particles[child_index][0]) {
                  const int other_child_index =
                      m_nodes_first_child[index] + (1 - c);
                  m_nodes_first_child[index] =
                      m_nodes_first_child[other_child_index];
                  m_nodes_centre[index] = m_nodes_centre[other_child_index];
                  m_nodes_radius[index] = m_nodes_radius[other_child_index];
                }
              }
            }
          }
        }
      });

      // filter out unchanged sub-trees
      i.filter([&](const auto &ci) {
        const int index = m_query.get_child_index(*ci);
        return m_query.is_leaf_node(*ci) ||
               m_nodes_particles[index][1] < update_start_index;
      });
    }
#ifndef __CUDA_ARCH__
    if (3 <= ABORIA_LOG_LEVEL) {
      print_tree();
    }
#endif
  }

  const BalltreeQuery<Traits> &get_query_impl() const { return m_query; }

  BalltreeQuery<Traits> &get_query_impl() { return m_query; }

private:
  struct accessor {
    const double_d *_p;
    const int *_particle_indicies;
    const int start;
    const int end;

    const double_d &operator[](const int i) const {
      return _p[_particle_indicies[i + start]];
    }

    const size_t size() const { return end - start; }
  };
  struct sort_by_split {
    double_d *_p;
    int *_particle_indicies;

    void get_centre_and_radius(const int p_start, const int p_end,
                               double_d &centre, double &radius) const {

      if (p_end - p_start <= 1) {
        centre = _p[_particle_indicies[p_start]];
        radius = 0.0;
      } else {
        Seb::Smallest_enclosing_ball<double, double_d, accessor> miniball(
            dimension, accessor{_p, _particle_indicies, p_start, p_end});

        auto ci = miniball.center_begin();
        for (int i = 0; i < dimension; ++i, ++ci) {
          centre[i] = *ci;
        }
        radius = miniball.radius();

        /*
           const double_d &a = _p[_particle_indicies[p_start]];
           const double_d &b = _p[_particle_indicies[p_end - 1]];
           centre = 0.5 * (a + b);
           radius = (a - centre).squaredNorm();
           for (int i = p_start; i < p_end; ++i) {
           const double_d &p = _p[_particle_indicies[i]];
           const double dx2 = (p - centre).squaredNorm();
           if (dx2 > radius) {
           radius = dx2;
           }
           }
           radius = std::sqrt(radius) + std::numeric_limits<double>::epsilon();
           */
      }
      // ensure that every point is in bounds
      double radius2 = std::pow(radius, 2);
      for (int i = p_start; i < p_end; ++i) {
        const double_d &p = _p[_particle_indicies[i]];
        const double r2 = (p - centre).squaredNorm();
        if (r2 > radius2) {
          radius = std::sqrt(r2);
          radius2 = std::pow(radius, 2);
          /*
          std::cout << "EROROROROR " << p << " " << centre << " "
                    << radius - (p - centre).norm() << " "
                    << std::numeric_limits<double>::epsilon() << std::endl;
                    */
        }
      }
    }

    // note: centre and radius define an exact circle enclosing the points
    // (i.e. will pass through the points that are maximally separated)
    int partition(const int p_start, const int p_end, const double_d &centre,
                  const double &radius) const {
      // if only one particle do nothing
      if (p_end - p_start <= 1) {
        return p_start;
      }
      vint2 outer;
      double max_dx2 = 0;
      for (int i = p_start; i < p_end; ++i) {
        for (int j = p_start; j < p_end; ++j) {
          const double_d &pi = _p[_particle_indicies[i]];
          const double_d &pj = _p[_particle_indicies[j]];
          const double dx2 = (pi - pj).squaredNorm();
          if (dx2 > max_dx2) {
            outer[0] = i;
            outer[1] = j;
            max_dx2 = dx2;
          }
        }
      }
      /*
      vint2 outer2;
      // why does this need to be 0.7???
      const double r2 = std::pow(0.7 * radius, 2);
      double max_dx2_2 = 0;
      for (int i = p_start; i < p_end; ++i) {
        const double_d &pi = _p[_particle_indicies[i]];
        if ((pi - centre).squaredNorm() < r2)
          continue;
        for (int j = p_start; j < p_end; ++j) {
          const double_d &pj = _p[_particle_indicies[j]];
          const double dx2 = (pi - pj).squaredNorm();
          if (dx2 > max_dx2_2) {
            outer2[0] = i;
            outer2[1] = j;
            max_dx2_2 = dx2;
          }
        }
      }
      if ((outer != outer2).any() &&
          (outer[0] != outer2[1] || outer[1] != outer2[0])) {
        std::cout << "EROREROREOREREORR: " << outer << " " << outer2 << " "
                  << max_dx2 << " " << max_dx2_2 << " " << radius << " "
                  << (_p[_particle_indicies[outer[0]]] - centre).norm() - radius
                  << " "
                  << (_p[_particle_indicies[outer[1]]] - centre).norm() - radius
                  << std::endl;
      }
      */
      // point projection of each particle onto line defined by AB:
      // dot(AP,AB) / dot(AB,AB)
      const double_d &a = _p[_particle_indicies[outer[0]]];
      const double_d &b = _p[_particle_indicies[outer[1]]];
      /*
         const double_d ab = b - a;
         const double dot_ab_ab = ab.squaredNorm();
         for (int i = p_start; i < p_end; ++i) {
         const double_d &pi = _p[_particle_indicies[i]];
         _partition[i] = (pi - a).squaredNorm() <= (pi - b).squaredNorm();
         }
         */

      // partition particles by closest point
      auto split_it = detail::partition(
          _particle_indicies + p_start, _particle_indicies + p_end,
          [&](const int i) {
            const double_d &pi = _p[i];
            return (pi - a).squaredNorm() <= (pi - b).squaredNorm();
          });
      return split_it - _particle_indicies;
    }
  };
  void build_tree() {
    const size_t num_points = this->m_alive_indices.size();

    // setup nodes
    vector_int3 parents_leaf(1, vint3(0, 0, num_points));
    vector_int parents_split(1);

    // do intial sorting by split direction
    auto sort_function = sort_by_split{
        iterator_to_raw_pointer(get<position>(this->m_particles_begin)),
        iterator_to_raw_pointer(m_particle_indicies.begin())};

    // setup tree with a single node
    m_nodes_first_child.resize(1);
    m_nodes_particles.resize(1);
    m_nodes_particles[0] = vint2(0, num_points);
    m_nodes_centre.resize(1);
    m_nodes_radius.resize(1);
    sort_function.get_centre_and_radius(m_nodes_particles[0][0],
                                        m_nodes_particles[0][1],
                                        m_nodes_centre[0], m_nodes_radius[0]);

    m_number_of_levels = 1;

    if (num_points > this->m_n_particles_in_leaf) {
      m_nodes_first_child[0] = 1;
      parents_leaf[0][1] =
          sort_function.partition(parents_leaf[0][0], parents_leaf[0][2],
                                  m_nodes_centre[0], m_nodes_radius[0]);
    } else {
      m_nodes_first_child[0] = -1;
      parents_leaf.clear();
    }

    while (!parents_leaf.empty()) {
      m_number_of_levels++;
      const int num_parents = parents_leaf.size();
      LOG(3, "build_tree(balltree): building level " << m_number_of_levels);

      // categorise (leaf/non-leaf) parent's children
      vector_bool2 parents_child_is_non_leaf(num_parents);
      detail::transform(
          parents_leaf.begin(), parents_leaf.end(),
          parents_child_is_non_leaf.begin(),
          [_threshold = this->m_n_particles_in_leaf](const vint3 &leaf) {
            const int nparent = leaf[2] - leaf[0];
            const int nchild1 = leaf[1] - leaf[0];
            const int nchild2 = nparent - nchild1;
            return vbool2(nchild1 > _threshold, nchild2 > _threshold);
          });

      /*
         std::cout << "parents_child_is_non_leaf: ";
         for (int i = 0; i < num_parents; ++i) {
         std::cout << parents_child_is_non_leaf[i] << " ";
         }
         std::cout << std::endl;
         */

      // find indicies of first children for next iteration
      vector_int parents_first_child_index(num_parents);
      detail::transform_exclusive_scan(
          parents_child_is_non_leaf.begin(), parents_child_is_non_leaf.end(),
          parents_first_child_index.begin(),
          [](const vbool2 &i) { return i[0] + i[1]; }, 0, detail::plus());

      /*
         std::cout << "parents_first_child_index: ";
         for (int i = 0; i < num_parents; ++i) {
         std::cout << parents_first_child_index[i] << " ";
         }
         std::cout << std::endl;
         */

      // form children in main data structure
      const int num_children = 2 * parents_leaf.size();
      const int level_start_index = m_nodes_first_child.size();
      m_nodes_first_child.resize(level_start_index + num_children);
      m_nodes_particles.resize(level_start_index + num_children);
      m_nodes_centre.resize(level_start_index + num_children);
      m_nodes_radius.resize(level_start_index + num_children);

      auto new_level_it =
          Traits::make_zip_iterator(Traits::make_tuple(
              m_nodes_first_child.begin(), m_nodes_particles.begin(),
              m_nodes_centre.begin(), m_nodes_radius.begin())) +
          level_start_index;

      detail::tabulate(
          new_level_it, new_level_it + num_children,
          [_next_level_start_index = level_start_index + num_children,
           _parents_leaf = iterator_to_raw_pointer(parents_leaf.begin()),
           _parents_first_child =
               iterator_to_raw_pointer(parents_first_child_index.begin()),
           _parents_child_is_non_leaf =
               iterator_to_raw_pointer(parents_child_is_non_leaf.begin()),
           _sort_function = sort_function](const int i) {
            const int parent_index = i / 2;
            const vint3 &leaf = _parents_leaf[parent_index];
            const vbool2 &is_non_leaf =
                _parents_child_is_non_leaf[parent_index];

            int first_child;
            vint2 particles;

            if (i % 2) {
              // second child
              particles = vint2(leaf[1], leaf[2]);
              if (is_non_leaf[1]) {
                // has a child
                first_child =
                    _next_level_start_index +
                    2 * (_parents_first_child[parent_index] + is_non_leaf[0]);

              } else {
                // no children
                first_child = -1;
              }
            } else {
              // first child
              particles = vint2(leaf[0], leaf[1]);
              if (is_non_leaf[0]) {
                // has a child
                first_child = _next_level_start_index +
                              2 * _parents_first_child[parent_index];

              } else {
                // no children
                first_child = -1;
              }
            }

            double_d centre;
            double radius;
            _sort_function.get_centre_and_radius(particles[0], particles[1],
                                                 centre, radius);

            return Traits::make_tuple(first_child, particles, centre, radius);
          });

      /*
         std::cout << "m_nodes_first_child: ";
         for (int i = level_start_index; i < level_start_index + num_children;
         ++i) {
         std::cout << m_nodes_first_child[i] << " ";
         }
         std::cout << std::endl;
         std::cout << "m_nodes_particles: ";
         for (int i = level_start_index; i < level_start_index + num_children;
         ++i) {
         std::cout << m_nodes_particles[i] << " ";
         }
         std::cout << std::endl;
         */

      // setup parents for next iteration (i.e. remove all leaf nodes)

      const int new_num_parents = parents_first_child_index.back() +
                                  parents_child_is_non_leaf.back()[0] +
                                  parents_child_is_non_leaf.back()[1];

      vector_int3 parents_leaf2(new_num_parents);

#if defined(__CUDACC__)
      typedef typename thrust::detail::iterator_category_to_system<
          typename Traits::vector_int::iterator::iterator_category>::type
          system;
      thrust::counting_iterator<int, system> count(0);
#else
      auto count = Traits::make_counting_iterator(0);
#endif
      detail::for_each(
          count, count + num_parents,
          [_parents_leaf = iterator_to_raw_pointer(parents_leaf.begin()),
           _nodes_centre = iterator_to_raw_pointer(m_nodes_centre.begin() +
                                                   level_start_index),
           _nodes_radius = iterator_to_raw_pointer(m_nodes_radius.begin() +
                                                   level_start_index),
           _parents_leaf2 = iterator_to_raw_pointer(parents_leaf2.begin()),
           _parents_child_is_non_leaf =
               iterator_to_raw_pointer(parents_child_is_non_leaf.begin()),
           _parents_first_child_index =
               iterator_to_raw_pointer(parents_first_child_index.begin()),
           _sort_function = sort_function](const int i) {
            int child_index = _parents_first_child_index[i];
            const vint3 &leaf = _parents_leaf[i];
            for (int c = 0; c < 2; ++c) {
              if (_parents_child_is_non_leaf[i][c]) {
                const double_d &centre = _nodes_centre[2 * i + c];
                const double &radius = _nodes_radius[2 * i + c];
                const int split = _sort_function.partition(leaf[c], leaf[c + 1],
                                                           centre, radius);
                _parents_leaf2[child_index++] =
                    vint3(leaf[c], split, leaf[c + 1]);
              }
            }
          });

      parents_leaf.swap(parents_leaf2);
    }
#ifndef __CUDA_ARCH__
    if (3 <= ABORIA_LOG_LEVEL) {
      print_tree();
    }
#endif
  }

  vector_int m_nodes_first_child;
  vector_int2 m_nodes_particles;
  vector_double_d m_nodes_centre;
  vector_double m_nodes_radius;

  vector_int m_particle_indicies;
  int m_number_of_levels;
  BalltreeQuery<Traits> m_query;
}; // namespace Aboria

template <typename Query> class BalltreeChildIterator {
  typedef bbox<Query::dimension> box_type;
  template <typename Traits> friend struct BalltreeQuery;

public:
  struct value_type {
    int high;
    const int *parent;
  };
  value_type m_data;

  typedef const value_type *pointer;
  typedef std::forward_iterator_tag iterator_category;
  typedef const value_type &reference;
  typedef std::ptrdiff_t difference_type;

  BalltreeChildIterator() : m_data{2, nullptr} {}

  BalltreeChildIterator(const int *parent)
      : m_data{parent != nullptr ? 0 : 2, parent} {}

  bool is_high() const { return m_data.high > 0; }

  int get_child_number() const { return m_data.high; }

  int distance_to_end() const { return 2 - m_data.high; }

  reference operator*() const { return dereference(); }

  reference operator->() const { return dereference(); }

  BalltreeChildIterator &operator++() {
    increment();
    return *this;
  }

  BalltreeChildIterator operator++(int) {
    BalltreeChildIterator tmp(*this);
    operator++();
    return tmp;
  }

  CUDA_HOST_DEVICE
  BalltreeChildIterator operator+(const int n) {
    BalltreeChildIterator tmp(*this);
    tmp.increment(n);
    return tmp;
  }

  inline bool operator==(const BalltreeChildIterator &rhs) const {
    return equal(rhs);
  }

  inline bool operator!=(const BalltreeChildIterator &rhs) const {
    return !operator==(rhs);
  }

  inline bool operator==(const bool rhs) const { return equal(rhs); }

  inline bool operator!=(const bool rhs) const { return !operator==(rhs); }

private:
  bool equal(BalltreeChildIterator const &other) const {
    return m_data.parent == other.m_data.parent &&
           m_data.high == other.m_data.high;
  }

  bool equal(const bool other) const { return (m_data.high < 2) == other; }

  reference dereference() const { return m_data; }

  void increment() { ++m_data.high; }

  void increment(const int n) { m_data.high += n; }
};

/// @copydetails NeighbourQueryBase
///
/// @brief This is a query object for the @ref Balltree spatial data
/// structure
///
template <typename Traits> struct BalltreeQuery {
  friend class Balltree<Traits>;
  const static unsigned int dimension = Traits::dimension;
  const static unsigned int m_max_tree_depth = 32 - 2;

  typedef Traits traits_type;
  typedef typename Traits::raw_pointer raw_pointer;
  typedef typename Traits::double_d double_d;
  typedef typename Traits::bool_d bool_d;
  typedef typename Traits::int_d int_d;
  typedef typename Traits::unsigned_int_d unsigned_int_d;
  typedef typename Traits::position position;

  template <int LNormNumber>
  using query_iterator = balltree_query_iterator<BalltreeQuery, LNormNumber>;

  template <int LNormNumber>
  using bounds_query_iterator =
      balltree_box_query_iterator<BalltreeQuery, LNormNumber>;

  typedef depth_first_iterator<BalltreeQuery> all_iterator;
  typedef BreadthFirstSearchIterator<BalltreeQuery> breadth_first_iterator;
  typedef BalltreeChildIterator<BalltreeQuery> child_iterator;

  typedef typename child_iterator::value_type value_type;
  typedef typename child_iterator::reference reference;
  typedef typename child_iterator::pointer pointer;
  typedef ranges_iterator<Traits> particle_iterator;
  typedef bbox<dimension> box_type;

  bool_d m_periodic;
  bbox<dimension> m_bounds;
  raw_pointer m_particles_begin;
  raw_pointer m_particles_end;

  /// number of entries in m_nodes_first_child, note that this is one more
  /// than the total number of buckets when a user iterates through them all
  size_t m_number_of_buckets;
  size_t m_number_of_levels;

  int *m_nodes_first_child;
  vint2 *m_nodes_particles;
  double_d *m_nodes_centre;
  double *m_nodes_radius;
  int m_dummy_root{-1};

  size_t *m_id_map_key;
  size_t *m_id_map_value;

  const box_type &get_bounds() const { return m_bounds; }
  const bool_d &get_periodic() const { return m_periodic; }

private:
  inline int get_parent_index(reference b) const {
    return b.parent - m_nodes_first_child;
  }

protected:
  inline int get_child_index(reference b) const { return *b.parent + b.high; }

public:
  /*
   * functions for id mapping
   */
  ABORIA_HOST_DEVICE_IGNORE_WARN
  CUDA_HOST_DEVICE
  raw_pointer find(const size_t id) const {
    const size_t n = number_of_particles();
    size_t *last = m_id_map_key + n;
    size_t *first = detail::lower_bound(m_id_map_key, last, id);
    if ((first != last) && !(id < *first)) {
      return m_particles_begin + m_id_map_value[first - m_id_map_key];
    } else {
      return m_particles_begin + n;
    }
  }

  /*
   * functions for tree_query_iterator
   */
  bool is_leaf_node(reference bucket) const {
    return m_nodes_first_child[get_child_index(bucket)] < 0;
  }
  static bool is_tree() { return true; }

  /*
   * end functions for tree_query_iterator
   */

  ///
  /// @copydoc NeighbourQueryBase::get_root() const
  ///
  child_iterator get_root() const { return ++child_iterator(&m_dummy_root); }

  child_iterator get_children() const {
    return child_iterator(m_nodes_first_child);
  }

  child_iterator get_children(const child_iterator &ci) const {
    if (!is_leaf_node(*ci)) {
      return child_iterator(m_nodes_first_child + get_child_index(*ci));
    } else {
      return child_iterator();
    }
  }

  ///
  /// @copydoc NeighbourQueryBase::num_children() const
  ///
  size_t num_children() const { return m_nodes_first_child != nullptr ? 2 : 0; }

  size_t num_children(const child_iterator &ci) const {
    if (is_leaf_node(*ci)) {
      return 0;
    } else {
      return 2;
    }
  }

  ///
  /// @copydoc NeighbourQueryBase::get_bounds(const child_iterator &ci) const
  ///
  /// actual bounds are a hyper-sphere with outer points given by
  /// m_nodes_particles. Since a bounding box in Aboria is assumed to be an
  /// axis-aligned hyper-rectangle, return smallest one of these covering
  /// the hyper-sphere
  ///
  const box_type get_bounds(const child_iterator &ci) const {
    const int cindex = get_child_index(*ci);
    const double_d &c = m_nodes_centre[cindex];
    const double &r = m_nodes_radius[cindex];
    return box_type(c - r, c + r + std::numeric_limits<double>::epsilon());
  }

  particle_iterator get_bucket_particles(reference bucket) const {
#ifndef __CUDA_ARCH__
    LOG(4, "\tget_bucket_particles: looking in bucket with idx = "
               << *bucket.parent + bucket.high);
#endif
    if (!is_leaf_node(bucket)) {
      return particle_iterator();
    }

    const int cindex = get_child_index(bucket);

    return particle_iterator(m_particles_begin + m_nodes_particles[cindex][0],
                             m_particles_begin + m_nodes_particles[cindex][1]);
  }

  void go_to(const double_d &position, child_iterator &ci) const {
    if ((*ci).parent == m_dummy_root) {
      // if root node do nothing
      return;
    }

    // find line between two children's centres
    const double_d &a = m_nodes_centre[get_child_index(*ci)];
    const double_d &b = m_nodes_centre[get_child_index(*(ci + 1))];
    const double_d ab = b - a;

    // project position onto this line
    const double position_projection =
        (position - a).dot(ab) / ab.squaredNorm();

    // if position is on the side of the other child, advance iterator
    if (position_projection > 0.5)
      ++ci;
  }

  child_iterator get_bucket(const double_d &position) const {
    child_iterator i = get_children();
    go_to(position, i);

    while (!is_leaf_node(*i)) {
      i = get_children(i);
      go_to(position, i);
    }

    return i;
  }

  size_t get_parent_index(const child_iterator &ci) const {
    return get_parent_index(*ci);
  }

  const box_type &get_parent_bounds(const child_iterator &ci) const {
    return ci.m_data.bounds;
  }

  size_t get_bucket_index(reference bucket) const {
    // minus one since we disregard the root index
    return get_child_index(bucket) - 1;
  }

  size_t number_of_buckets() const {
    // minus one since we disregard the root index
    return m_number_of_buckets - 1;
  }

  template <int LNormNumber>
  query_iterator<LNormNumber>
  get_buckets_near_point(const double_d &position,
                         const double max_distance) const {
#ifndef __CUDA_ARCH__
    LOG(4, "\tget_buckets_near_point: position = "
               << position << " max_distance= " << max_distance);
#endif
    return query_iterator<LNormNumber>(get_children(), position,
                                       double_d::Constant(max_distance),
                                       m_number_of_levels, this);
  }

  ///
  /// @copydoc NeighbourQueryBase::get_buckets_near_bucket()
  ///
  template <int LNormNumber>
  bounds_query_iterator<LNormNumber>
  get_buckets_near_bucket(const box_type &bounds,
                          const double max_distance) const {
#ifndef __CUDA_ARCH__
    LOG(4, "\tget_buckets_near_bucket: bounds = "
               << bounds << " max_distance = " << max_distance);
#endif

    return bounds_query_iterator<LNormNumber>(
        get_children(), bounds, max_distance, m_number_of_levels, this);
  }

  template <int LNormNumber>
  query_iterator<LNormNumber>
  get_buckets_near_point(const double_d &position,
                         const double_d &max_distance) const {
#ifndef __CUDA_ARCH__
    LOG(4, "\tget_buckets_near_point: position = "
               << position << " max_distance= " << max_distance);
#endif
    return query_iterator<LNormNumber>(get_children(), position, max_distance,
                                       m_number_of_levels, this);
  }

  all_iterator get_subtree(const child_iterator &ci) const {
    return all_iterator(get_children(ci), m_number_of_levels, this);
  }

  all_iterator get_subtree() const {
    return all_iterator(get_children(), m_number_of_levels, this);
  }

  ///
  /// @copydoc NeighbourQueryBase::get_breadth_first(const child_iterator&)
  /// const
  ///
  breadth_first_iterator breadth_first(const child_iterator &ci) const {
    return breadth_first_iterator(ci, this);
  }

  ///
  /// @copydoc NeighbourQueryBase::get_breadth_first() const
  ///
  breadth_first_iterator breadth_first() const {
    return breadth_first_iterator(get_children(), this);
  }

  size_t number_of_particles() const {
    return m_particles_end - m_particles_begin;
  }

  raw_pointer get_particles_begin() const { return m_particles_begin; }

  unsigned number_of_levels() const { return m_number_of_levels; }
};

} // namespace Aboria

#endif /* KDTREE_H_ */
