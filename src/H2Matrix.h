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

#ifndef H2_MATRIX_H_
#define H2_MATRIX_H_

#include "DualBreadthFirstSearch.h"
#include "detail/FastMultipoleMethod.h"
#include <Eigen/SparseCore>

namespace Aboria {

template <typename Expansions, typename RowParticles, typename ColParticles>
class H2Matrix {
  using row_query_t = typename RowParticles::query_type;
  using col_query_t = typename ColParticles::query_type;
  using traits_t = row_query_t::traits_type;
  using row_reference = row_query_t::reference;
  using col_reference = col_query_t::reference;
  using row_pointer = row_query_t::pointer;
  using col_pointer = col_query_t::pointer;

  static const unsigned int dimension = row_query_t::dimension;
  using position = position_d<dimension>;
  using row_child_iterator = row_query_t::child_iterator;
  using col_child_iterator = col_query_t::child_iterator;

  using dual_bf_search_t = DualBreadthFirstSearch<row_query_t, col_query_t>;
  using level_t = dual_bf_search_t::value_type;
  using pair_t = level_t::value_type;
  using tree_t = traits_type::template vector_type<level_t>::type;

  using p_vector_t = Expansions::p_vector_type;
  using m_vector_t = Expansions::m_vector_type;
  using l_vector_t = Expansions::l_vector_type;
  using l2p_matrix_t = Expansions::l2p_matrix_type;
  using p2m_matrix_t = Expansions::p2m_matrix_type;
  using p2p_matrix_t = Expansions::p2p_matrix_type;
  using l2l_matrix_t = Expansions::l2l_matrix_type;
  using m2m_matrix_t = Expansions::m2m_matrix_type;
  using m2l_matrix_t = Expansions::m2l_matrix_type;

  struct dual_node_t {
    pair_t m_ci_pair;
    int m_index_of_children;
    enum { P2P, M2L } tag;
    union m_op {
      p2p_matrix_t p2p;
      m2l_matrix_t m2l;
    };

    node_t(const pair_t &ci_pair) : m_ci_pair(ci_pair) {}
  };

  struct row_node_t {
    row_child_iterator m_ci;
    int m_index_of_children;
    l_vector_t m_l;
    l2l_matrix_t l2l;
    p2l_matrix_t p2l;

    row_node_t(const row_child_iterator &ci) : m_ci(ci) {}
  };

  struct col_node_t {
    col_child_iterator m_ci;
    int m_index_of_children;
    m_vector_t m_m;
    l2l_matrix_t m2m;
    p2l_matrix_t m2p;

    col_node_t(const col_child_iterator &ci) : m_ci(ci) {}
  };

  using dual_level_t =
      typename traits_type::template vector_type<dual_node_t>::type;
  using row_level_t =
      typename traits_type::template vector_type<row_node_t>::type;
  using col_level_t =
      typename traits_type::template vector_type<col_node_t>::type;

  using dual_tree_t =
      typename traits_type::template vector_type<dual_level_t>::type;
  using row_tree_t =
      typename traits_type::template vector_type<row_level_t>::type;
  using col_tree_t =
      typename traits_type::template vector_type<col_level_t>::type;

  typedef typename Query::particle_iterator particle_iterator;
  typedef typename particle_iterator::reference particle_reference;
  typedef bbox<dimension> box_type;

  // vectors used to cache values
  mutable row_tree_t m_row_tree;
  mutable col_tree_t m_col_tree;
  dual_node_t m_dual_tree;

  Expansions m_expansions;

  const row_query_t *m_row_query;
  const col_query_t *m_col_query;

public:
  template <typename RowParticles>
  H2Matrix(const RowParticles &row_particles, const ColParticles &col_particles,
           const Expansions &expansions)
      : m_row_query(&row_particles.get_query()),
        m_col_query(&col_particles.get_query()), m_expansions(expansions) {

    // generate h2 matrix
    LOG(2, "H2Matrix: creating matrix using "
               << row_particles.size() << " row particles and "
               << col_particles.size() << " column particles");

    const bool row_equals_col = static_cast<const void *>(&row_particles) ==
                                static_cast<const void *>(&col_particles);

    // construct row_tree
    {
      LOG(2, "\tconstructing row tree...");
      bf_search_t bf_search(m_row_query);
      m_row_tree.clear();
      while (!bf_search) {
        ++bf_search;
        m_row_tree.emplace_back(bf_search->size());
        detail::tabulate(
            m_row_tree->begin(), m_row_tree->end(), [](const int i) {
              row_node_t node;
              node.m_ci = (*bf_search)[i];
              node.m_index_of_children = bf_search->get_num_children()[i];
              if (m_row_query->is_leaf(node.m_ci)) {
                node.p2m.resize();
                m_expansions.P2M_matrix(node.p2m, *m_row_query, node.m_ci);
              } else {
                for (child) {
                  node.m2m m_expansions.M2M_matrix(chidl.node.p2m, *m_row_query,
                                                   child.node.m_ci, node.m_ci);
                }
              }
              return row_node_t(, bf_search->get_num_children[i], )
            });
      }
    }

    // construct col_tree
    {
      LOG(2, "\tconstructing col tree...");
      bf_search_t bf_search(m_col_query);
      m_col_tree.clear();
      while (!bf_search) {
        ++bf_search;
        m_col_tree.emplace_back(bf_search->size());
        detail::copy(bf_search->begin(), bf_search->end(),
                     m_col_tree.back().begin());
      }
    }

    // construct dual_tree
    {
      LOG(2, "\tconstructing dual tree...");
      dual_bf_search_t dual_bf_search(m_row_query, m_col_query);
      m_dual_tree.clear();
      while (!dual_bf_search) {
        ++dual_bf_search;
        m_dual_tree.emplace_back(dual_bf_search->size());
        detail::tabulate(m_dual_tree.back().begin(), m_dual_tree.back().end(),
                         [](const int i) {
                           auto &pair = (*dual_bf_search)[i];
                           const node_t new_node = pair;
                           if (new_node.is_leaf()) {
                             // "turn-off" node
                             pair = level_t::value_type();
                           }
                           return new_node;
                         });
      }
      ASSERT(m_dual_tree.size() >= m_row_tree.size());
      ASSERT(m_dual_tree.size() >= m_col_tree.size());
    }

    LOG(2, "\tdone");
  }

  H2Matrix(const H2Matrix &matrix) = default;

  H2Matrix(H2Matrix &&matrix) = default;

  ~H2Matrix() = default;

  // target_vector += A*source_vector
  template <typename VectorTypeTarget, typename VectorTypeSource>
  void matrix_vector_multiply(VectorTypeTarget &target_vector,
                              const VectorTypeSource &source_vector) const {

    // for all leaf nodes setup source vector
    for (auto &bucket : m_query->get_subtree()) {
      if (m_query->is_leaf_node(bucket)) { // leaf node
        const size_t index = m_query->get_bucket_index(bucket);
        for (size_t i = 0; i < m_col_indices[index].size(); ++i) {
          m_source_vector[index][i] = source_vector[m_col_indices[index][i]];
        }
      }
    }

    // upward sweep of the column tree
    for (auto col_level = m_col_tree.rbegin(); col_level != m_col_tree.rend();
         ++col_level) {
      mvm_upward_sweep(col_level, source_vector);
    }

    // downward sweep of row tree calculating interactions.
    for (auto row_level = m_row_tree.begin(),
              auto dual_level = m_dual_tree.begin();
         row_level != m_row_tree.rend(); ++row_level) {
      mvm_downward_sweep(row_level, target_vector);
      mvm_pair_interactions(dual_level, target_vector);
    }
  }

  /// Convert H2 Matrix to an extended sparse matrix
  ///
  /// This creates an Eigen sparse matrix A that can be applied to
  /// by the vector $[x W g]$ to produce $[y 0 0]$, i.e. $A*[x W g]' = [y 0
  /// 0]$, where $x$ and $y$ are the input/output vector, i.e. $y =
  /// H2Matrix*x$, and $W$ and $g$ are the multipole and local coefficients
  /// respectively.
  ///
  /// A is given by the block matrix representation:
  ///
  /// Sushnikova, D., & Oseledets, I. V. (2014). Preconditioners for
  /// hierarchical matrices based on their extended sparse form. Retrieved
  /// from http://arxiv.org/abs/1412.1253
  ///
  ///             |P2P  0     0      L2P| |x  |   |y|
  /// A*[x W g] = |0   M2L    L2L-I  -I |*|W  | = |0|
  ///             |0   M2M-I  0       0 | |g_n|   |0|
  ///             |P2M  -I    0       0 | |g_l|   |0|
  ///
  /// where $g_l$ is the $g$ vector for every leaf of the tree, and $g_n$ is
  /// the remaining nodes. Note M2M = L2L'.
  ///
  /// TODO: this will only work for rows == columns
  sparse_matrix_type gen_extended_matrix() const {

    const size_t size_x = m_col_particles->size();
    ASSERT(m_vector_type::RowsAtCompileTime != Eigen::Dynamic,
           "not using compile time m size");
    const size_t vector_size = m_vector_type::RowsAtCompileTime;

    const size_t size_W = m_W.size() * vector_size;
    const size_t size_g = m_g.size() * vector_size;
    const size_t n = size_x + size_W + size_g;

    std::vector<int> reserve(n, 0);

    // sizes for first x columns
    auto reserve0 = std::begin(reserve);
    for (size_t i = 0; i < m_source_vector.size(); ++i) {
      const size_t n = m_source_vector[i].size();
      // p2m
      std::transform(reserve0, reserve0 + n, reserve0,
                     [](const int count) { return count + vector_size; });
      // p2p
      for (size_t j = 0; j < m_strong_connectivity[i].size(); ++j) {
        size_t index = m_query->get_bucket_index(*m_strong_connectivity[i][j]);
        const size_t target_size = m_target_vector[index].size();
        std::transform(
            reserve0, reserve0 + n, reserve0,
            [target_size](const int count) { return count + target_size; });
      }
      reserve0 += m_source_vector[i].size();
    }

    // sizes for W columns
    // m2l
    for (size_t i = 0; i < m_source_vector.size(); ++i) {
      for (size_t j = 0; j < m_weak_connectivity[i].size(); ++j) {
        std::transform(reserve0 + i * vector_size,
                       reserve0 + (i + 1) * vector_size,
                       reserve0 + i * vector_size,
                       [](const int count) { return count + vector_size; });
      }
    }
    // m2m-I or -I
    std::transform(reserve0, reserve0 + size_W, reserve0,
                   [](const int count) { return count + 1; });

    auto subtree_range = m_query->get_subtree();
    for (auto bucket_it = subtree_range.begin();
         bucket_it != subtree_range.end(); ++bucket_it) {
      const size_t i = m_query->get_bucket_index(*bucket_it);
      // row is the index of the parent
      const size_t row_index = size_x + size_W + i * vector_size;
      // m2m
      if (!m_query->is_leaf_node(*bucket_it)) { // non leaf node
        for (child_iterator cj =
                 m_query->get_children(bucket_it.get_child_iterator());
             cj != false; ++cj) {
          const size_t j = m_query->get_bucket_index(*cj);
          std::transform(reserve0 + j * vector_size,
                         reserve0 + (j + 1) * vector_size,
                         reserve0 + j * vector_size,
                         [](const int count) { return count + vector_size; });
        }
      }
    }

    // sizes for g columns
    reserve0 += size_W;

    // looping over columns
    for (auto bucket_it = subtree_range.begin();
         bucket_it != subtree_range.end(); ++bucket_it) {
      const size_t i = m_query->get_bucket_index(*bucket_it);
      // col is the index of the parent
      const size_t col_index = size_x + size_W + i * vector_size;

      // L2P & -I
      const size_t target_size = m_target_vector[i].size();
      std::transform(reserve0 + i * vector_size,
                     reserve0 + (i + 1) * vector_size,
                     reserve0 + i * vector_size,
                     [&](const int count) { return count + target_size + 1; });

      // L2L
      if (!m_query->is_leaf_node(*bucket_it)) { // non leaf node
        for (child_iterator cj =
                 m_query->get_children(bucket_it.get_child_iterator());
             cj != false; ++cj) {
          std::transform(reserve0 + i * vector_size,
                         reserve0 + (i + 1) * vector_size,
                         reserve0 + i * vector_size,
                         [&](const int count) { return count + vector_size; });
        }
      }
    }

    LOG(2, "\tcreating " << n << "x" << n << " extended sparse matrix");
    LOG(3, "\tnote: vector_size = " << vector_size);
    LOG(3, "\tnote: size_W is = " << size_W);
    LOG(3, "\tnote: size_g is = " << size_g);
    for (size_t i = 0; i < n; ++i) {
      LOG(4, "\tfor column " << i << ", reserving " << reserve[i] << " rows");
    }

    // create matrix and reserve space
    sparse_matrix_type A(n, n);
    A.reserve(reserve);

    // fill in P2P
    size_t row_index = 0;
    // loop over rows, filling in columns as we go
    for (size_t i = 0; i < m_target_vector.size(); ++i) {
      for (size_t j = 0; j < m_strong_connectivity[i].size(); ++j) {
        size_t source_index =
            m_query->get_bucket_index(*m_strong_connectivity[i][j]);
        // loop over n particles (the n rows in this bucket)
        for (size_t p = 0; p < m_target_vector[i].size(); ++p) {
          // p2p - loop over number of source particles (the columns)
          for (size_t sp = 0; sp < m_source_vector[source_index].size(); ++sp) {
            A.insert(row_index + p, m_ext_indicies[source_index] + sp) =
                m_p2p_matrices[i][j](p, sp);
          }
        }
      }
      row_index += m_target_vector[i].size();
    }

    // fill in P2M
    // loop over cols this time
    for (size_t i = 0; i < m_source_vector.size(); ++i) {
      // loop over n particles (the n cols in this bucket)
      for (size_t p = 0; p < m_source_vector[i].size(); ++p) {
        // p2m - loop over size of multipoles (the rows)
        const size_t row_index = size_x + size_W + i * vector_size;
        for (int m = 0; m < vector_size; ++m) {
          A.insert(row_index + m, m_ext_indicies[i] + p) =
              m_p2m_matrices[i](m, p);
        }
      }
    }

    // fill in W columns
    // m2l
    // loop over rows
    for (size_t i = 0; i < m_source_vector.size(); ++i) {
      for (size_t j = 0; j < m_weak_connectivity[i].size(); ++j) {
        const size_t row_index = size_x + i * vector_size;
        const size_t index =
            m_query->get_bucket_index(*m_weak_connectivity[i][j]);
        const size_t col_index = size_x + index * vector_size;
        for (int im = 0; im < vector_size; ++im) {
          for (int jm = 0; jm < vector_size; ++jm) {
            A.insert(row_index + im, col_index + jm) =
                m_m2l_matrices[i][j](im, jm);
          }
        }
      }
    }

    // m2m - I or -I
    // looping over rows (parents)
    for (auto bucket_it = subtree_range.begin();
         bucket_it != subtree_range.end(); ++bucket_it) {
      const size_t i = m_query->get_bucket_index(*bucket_it);
      // row is the index of the parent
      const size_t row_index = size_x + size_W + i * vector_size;

      // -I
      const size_t col_index = size_x + i * vector_size;
      for (int im = 0; im < vector_size; ++im) {
        A.insert(row_index + im, col_index + im) = -1;
      }

      // m2m
      if (!m_query->is_leaf_node(*bucket_it)) { // non leaf node
        for (child_iterator cj =
                 m_query->get_children(bucket_it.get_child_iterator());
             cj != false; ++cj) {
          const size_t j = m_query->get_bucket_index(*cj);
          size_t col_index = size_x + j * vector_size;
          // std::cout << "inserting at ("<<row_index<<","<<col_index<<")" <<
          // std::endl;
          for (int im = 0; im < vector_size; ++im) {
            for (int jm = 0; jm < vector_size; ++jm) {
              A.insert(row_index + im, col_index + jm) =
                  m_l2l_matrices[j](jm, im);
            }
          }
        }
      }
    }

    // fill in g columns
    // looping over rows
    // L2P
    row_index = 0;
    for (size_t i = 0; i < m_target_vector.size(); ++i) {
      if (m_target_vector[i].size() > 0) { // leaf node
        size_t col_index = size_x + size_W + i * vector_size;
        for (size_t p = 0; p < m_target_vector[i].size(); ++p) {
          for (int m = 0; m < vector_size; ++m) {
            A.insert(row_index + p, col_index + m) = m_l2p_matrices[i](p, m);
          }
        }
        row_index += m_target_vector[i].size();
      }
    }

    // looping over columns
    for (auto bucket_it = subtree_range.begin();
         bucket_it != subtree_range.end(); ++bucket_it) {
      const size_t i = m_query->get_bucket_index(*bucket_it);
      // col is the index of the parent
      const size_t col_index = size_x + size_W + i * vector_size;

      // -I
      const size_t row_index = size_x + i * vector_size;
      for (int im = 0; im < vector_size; ++im) {
        A.insert(row_index + im, col_index + im) = -1;
      }

      // L2L
      if (!m_query->is_leaf_node(*bucket_it)) { // non leaf node
        for (child_iterator cj =
                 m_query->get_children(bucket_it.get_child_iterator());
             cj != false; ++cj) {
          const size_t j = m_query->get_bucket_index(*cj);
          size_t row_index = size_x + j * vector_size;
          // std::cout << "inserting at ("<<row_index<<","<<col_index<<")" <<
          // std::endl;
          for (int im = 0; im < vector_size; ++im) {
            for (int jm = 0; jm < vector_size; ++jm) {
              A.insert(row_index + im, col_index + jm) =
                  m_l2l_matrices[j](im, jm);
            }
          }
        }
      }
    }

    A.makeCompressed();

#ifndef NDEBUG
    for (size_t i = 0; i < n; ++i) {
      size_t count = 0;
      for (sparse_matrix_type::InnerIterator it(A, i); it; ++it) {
        count++;
        ASSERT(!std::isnan(it.value()),
               "entry (" << it.row() << "," << it.col() << ") is nan");
      }
      ASSERT(count == reserve[i], "column " << i << " reserved size "
                                            << reserve[i] << " and final count "
                                            << count << " do not agree");
      ASSERT(count > 0, "column " << i << " has zero entries");
    }
    /*
    row_index = 0;
    size_t col_index = 0;
    for (size_t i = 0; i < m_target_vector.size(); ++i) {
        for (size_t pi = 0; pi < m_target_vector[i].size(); ++pi,++row_index)
    { double sum = 0; size_t count = 0; col_index = 0; for (size_t j = 0; j <
    m_source_vector.size(); ++j) { for (size_t pj = 0; pj <
    m_source_vector[j].size();
    ++pj,++col_index) { for (sparse_matrix_type::InnerIterator
    it(A,col_index); it ;++it) { if (it.row() == row_index) {
                            //std::cout << "("<<it.row()<<","<<it.col()<<") =
    "<<it.value() << std::endl; sum += it.value()*m_source_vector[j][pj];
                            count++;
                        }
                    }
                }
            }
            for (size_t j = 0; j < m_W.size(); ++j) {
                for (int mj = 0; mj < vector_size; ++mj,++col_index) {
                    for (sparse_matrix_type::InnerIterator it(A,col_index); it
    ;++it) { if (it.row() == row_index) {
                            //std::cout << "("<<it.row()<<","<<it.col()<<") =
    "<<it.value() << std::endl; sum += it.value()*m_W[j][mj]; count++;
                        }
                    }

                }
            }
            for (size_t j = 0; j < m_g.size(); ++j) {
                for (int mj = 0; mj < vector_size; ++mj,++col_index) {
                    for (sparse_matrix_type::InnerIterator it(A,col_index); it
    ;++it) { if (it.row() == row_index) {
                            //std::cout << "("<<it.row()<<","<<it.col()<<") =
    "<<it.value() << std::endl; sum += it.value()*m_g[j][mj]; count++;
                        }
                    }

                }
            }
            std::cout << "X: row "<<row_index<<" sum = "<<sum<<" should be
    "<<m_target_vector[i][pi] << std::endl; std::cout << "count is
    "<<count<<std::endl; ASSERT(count > 0,"count == 0");
        }
    }
    row_index = size_x;
    for (size_t i = 0; i < m_W.size(); ++i) {
        for (int mi = 0; mi < vector_size; ++mi,++row_index) {
            double sum = 0;
            size_t count = 0;
            col_index = 0;
            for (size_t j = 0; j < m_source_vector.size(); ++j) {
                for (size_t pj = 0; pj < m_source_vector[j].size();
    ++pj,++col_index) { for (sparse_matrix_type::InnerIterator
    it(A,col_index); it ;++it) { if (it.row() == row_index) {
                            //std::cout << "("<<it.row()<<","<<it.col()<<") =
    "<<it.value() << std::endl; sum += it.value()*m_source_vector[j][pj];
                            count++;
                        }
                    }
                }
            }
            for (size_t j = 0; j < m_W.size(); ++j) {
                for (int mj = 0; mj < vector_size; ++mj,++col_index) {
                    for (sparse_matrix_type::InnerIterator it(A,col_index); it
    ;++it) { if (it.row() == row_index) {
                            //std::cout << "("<<it.row()<<","<<it.col()<<") =
    "<<it.value() << std::endl; sum += it.value()*m_W[j][mj]; count++;
                        }
                    }

                }
            }
            for (size_t j = 0; j < m_g.size(); ++j) {
                for (int mj = 0; mj < vector_size; ++mj,++col_index) {
                    for (sparse_matrix_type::InnerIterator it(A,col_index); it
    ;++it) { if (it.row() == row_index) {
                            //std::cout << "("<<it.row()<<","<<it.col()<<") =
    "<<it.value() << std::endl; sum += it.value()*m_g[j][mj]; count++;
                        }
                    }

                }
            }
            std::cout << "W: row "<<row_index<<" sum = "<<sum<<" should be
    "<<0<< std::endl; std::cout << "count is "<<count<<std::endl; ASSERT(count
    > 0,"count == 0");
        }
    }
    row_index = size_x+size_W;
    for (size_t i = 0; i < m_g.size(); ++i) {
        for (int mi = 0; mi < vector_size; ++mi,++row_index) {
            double sum = 0;
            size_t count = 0;
            col_index = 0;
            for (size_t j = 0; j < m_source_vector.size(); ++j) {
                for (size_t pj = 0; pj < m_source_vector[j].size();
    ++pj,++col_index) { for (sparse_matrix_type::InnerIterator
    it(A,col_index); it ;++it) { if (it.row() == row_index) {
                            //std::cout << "("<<it.row()<<","<<it.col()<<") =
    "<<it.value() << std::endl; sum += it.value()*m_source_vector[j][pj];
                            count++;
                        }
                    }
                }
            }
            for (size_t j = 0; j < m_W.size(); ++j) {
                for (int mj = 0; mj < vector_size; ++mj,++col_index) {
                    for (sparse_matrix_type::InnerIterator it(A,col_index); it
    ;++it) { if (it.row() == row_index) {
                            //std::cout << "("<<it.row()<<","<<it.col()<<") =
    "<<it.value() << std::endl; sum += it.value()*m_W[j][mj]; count++;
                        }
                    }

                }
            }
            for (size_t j = 0; j < m_g.size(); ++j) {
                for (int mj = 0; mj < vector_size; ++mj,++col_index) {
                    for (sparse_matrix_type::InnerIterator it(A,col_index); it
    ;++it) { if (it.row() == row_index) {
                            //std::cout << "("<<it.row()<<","<<it.col()<<") =
    "<<it.value() << std::endl; sum += it.value()*m_g[j][mj]; count++;
                        }
                    }

                }
            }
            std::cout << "g: row "<<row_index<<" sum = "<<sum<<" should be
    "<<0
    << std::endl; std::cout << "count is "<<count<<std::endl; ASSERT(count >
    0,"count == 0");
        }
    }
    */

#endif

    LOG(2, "\tdone");
    return A;
  }

  /// Convert H2 Matrix to an stripped down extended sparse matrix
  ///
  /// This creates an Eigen sparse matrix A that can be used as a
  /// preconditioner for the true extended sparse matrix
  ///
  /// A is given by the block matrix representation:
  ///
  /// Sushnikova, D., & Oseledets, I. V. (2014). Preconditioners for
  /// hierarchical matrices based on their extended sparse form. Retrieved
  /// from http://arxiv.org/abs/1412.1253
  ///
  ///             |I    0     0 | |x  |
  /// A*[x W g] = |0   M2L    0 |*|W  |
  ///             |0    0     I | |g  |
  ///
  /// TODO: this will only work for rows == columns
  sparse_matrix_type gen_stripped_extended_matrix() const {

    const size_t size_x = m_col_particles->size();
    ASSERT(m_vector_type::RowsAtCompileTime != Eigen::Dynamic,
           "not using compile time m size");
    const size_t vector_size = m_vector_type::RowsAtCompileTime;

    const size_t size_W = m_W.size() * vector_size;
    const size_t size_g = m_g.size() * vector_size;
    const size_t n = size_x + size_W + size_g;

    std::vector<int> reserve(n, 0);

    // sizes for first x columns
    auto reserve0 = std::begin(reserve);
    std::transform(reserve0, reserve0 + size_x, reserve0,
                   [](const int count) { return count + 1; });

    reserve0 += size_x;

    // sizes for W columns
    // m2l
    for (size_t i = 0; i < m_source_vector.size(); ++i) {
      if (m_weak_connectivity[i].size() == 0) {
        std::transform(reserve0 + i * vector_size,
                       reserve0 + (i + 1) * vector_size,
                       reserve0 + i * vector_size,
                       [](const int count) { return count + 1; });
      }

      for (size_t j = 0; j < m_weak_connectivity[i].size(); ++j) {
        std::transform(reserve0 + i * vector_size,
                       reserve0 + (i + 1) * vector_size,
                       reserve0 + i * vector_size,
                       [](const int count) { return count + vector_size; });
      }
    }

    // sizes for g columns
    reserve0 += size_W;

    std::transform(reserve0, reserve0 + size_g, reserve0,
                   [](const int count) { return count + 1; });

    LOG(2, "\tcreating " << n << "x" << n
                         << " approximate extended sparse matrix");
    LOG(3, "\tnote: vector_size = " << vector_size);
    LOG(3, "\tnote: size_W is = " << size_W);
    LOG(3, "\tnote: size_g is = " << size_g);
    for (size_t i = 0; i < n; ++i) {
      LOG(4, "\tfor column " << i << ", reserving " << reserve[i] << " rows");
    }

    // create matrix and reserve space
    sparse_matrix_type A(n, n);
    A.reserve(reserve);

    // fill in x columns
    // loop over rows
    size_t row_index = 0;
    for (size_t i = 0; i < m_target_vector.size(); ++i) {
      for (size_t p = 0; p < m_target_vector[i].size(); ++p) {
        A.insert(row_index + p, row_index + p) = 1;
      }
      row_index += m_target_vector[i].size();
    }

    // fill in W columns
    // m2l
    // loop over rows
    for (size_t i = 0; i < m_source_vector.size(); ++i) {
      const size_t row_index = size_x + i * vector_size;

      if (m_weak_connectivity[i].size() == 0) {
        for (int m = 0; m < vector_size; ++m) {
          A.insert(row_index + m, row_index + m) = 1;
        }
      }

      for (size_t j = 0; j < m_weak_connectivity[i].size(); ++j) {
        const size_t index =
            m_query->get_bucket_index(*m_weak_connectivity[i][j]);
        const size_t col_index = size_x + index * vector_size;
        for (int im = 0; im < vector_size; ++im) {
          for (int jm = 0; jm < vector_size; ++jm) {
            if (m_m2l_matrices[i][j](im, jm) < 1e-10) {
              std::cout << "tooo small!!!!! " << m_m2l_matrices[i][j](im, jm)
                        << std::endl;
            }
            A.insert(row_index + im, col_index + jm) =
                m_m2l_matrices[i][j](im, jm);
          }
        }
      }
    }

    // fill in g columns
    // looping over rows
    for (size_t i = 0; i < m_target_vector.size(); ++i) {
      const size_t row_index = size_x + size_W + i * vector_size;
      for (int m = 0; m < vector_size; ++m) {
        A.insert(row_index + m, row_index + m) = 1;
      }
    }

    A.makeCompressed();

#ifndef NDEBUG
    for (size_t i = 0; i < n; ++i) {
      size_t count = 0;
      for (sparse_matrix_type::InnerIterator it(A, i); it; ++it) {
        count++;
        ASSERT(!std::isnan(it.value()),
               "entry (" << it.row() << "," << it.col() << ") is nan");
      }
      ASSERT(count == reserve[i], "column " << i << " reserved size "
                                            << reserve[i] << " and final count "
                                            << count << " do not agree");
      ASSERT(count > 0, "column " << i << " has zero entries");
    }

#endif

    LOG(2, "\tdone");
    return A;
  }

  // this function returns a vector of indicies map that correspond to a
  // mapping between the extended column vector e_x to the standard column
  // particle vector x
  //
  // i.e. x = e_x[map]
  index_vector_type gen_column_map() const {
    index_vector_type map(m_col_particles->size());
    for (size_t i = 0; i < m_col_indices.size(); ++i) {
      for (size_t p = 0; p < m_col_indices[i].size(); ++p) {
        map[m_col_indices[i][p]] = m_ext_indicies[i] + p;
      }
    }
    return map;
  }

  // this function returns a vector of indicies map that correspond to a
  // mapping between the extended row vector e_b to the standard row particle
  // vector b
  //
  // i.e. b = e_b[map]
  index_vector_type gen_row_map() const {
    index_vector_type map(m_row_size);
    size_t index = 0;
    for (size_t i = 0; i < m_row_indices.size(); ++i) {
      for (size_t p = 0; p < m_row_indices[i].size(); ++p, index++) {
        map[m_row_indices[i][p]] = index;
      }
    }
    return map;
  }

  /*
  size_t get_column_extended_vector_size() const {
      const size_t size_x = m_col_particles.size();
      const size_t vector_size = m_vector_type::RowsAtCompileTime;
      const size_t size_W = m_W.size()*vector_size;
      const size_t size_g = m_g.size()*vector_size;
      const size_t n = size_x + size_W + size_g;
      return n;
  }
  size_t get_row_extended_vector_size() const {
      const size_t size_x = m_row_size;
      const size_t vector_size = m_vector_type::RowsAtCompileTime;
      const size_t size_W = m_W.size()*vector_size;
      const size_t size_g = m_g.size()*vector_size;
      const size_t n = size_x + size_W + size_g;
      return n;
  }
  */

  template <typename VectorTypeSource>
  column_vector_type
  gen_extended_vector(const VectorTypeSource &source_vector) const {
    const size_t size_x = source_vector.size();
    ASSERT(size_x == m_col_particles->size(),
           "source vector not same size as column particles");
    const size_t vector_size = m_vector_type::RowsAtCompileTime;
    const size_t size_W = m_W.size() * vector_size;
    const size_t size_g = m_g.size() * vector_size;
    const size_t n = size_x + size_W + size_g;

    column_vector_type extended_vector(n);

    // x
    for (auto &bucket : m_query->get_subtree()) {
      if (m_query->is_leaf_node(bucket)) { // leaf node
        const size_t index = m_query->get_bucket_index(bucket);
        for (size_t i = 0; i < m_col_indices[index].size(); ++i) {
          extended_vector[m_ext_indicies[index] + i] =
              source_vector[m_col_indices[index][i]];
        }
      }
    }

    // zero remainder
    for (int i = size_x; i < n; ++i) {
      extended_vector[i] = 0;
    }

    return extended_vector;
  }

  column_vector_type get_internal_state() const {
    const size_t size_x = m_ext_indicies.back() + m_source_vector.back().size();
    ASSERT(size_x == m_col_particles->size(),
           "source vector not same size as column particles");
    const size_t vector_size = m_vector_type::RowsAtCompileTime;
    const size_t size_W = m_W.size() * vector_size;
    const size_t size_g = m_g.size() * vector_size;
    const size_t n = size_x + size_W + size_g;

    column_vector_type extended_vector(n);
    LOG(2, "get_internal_state:");

    // x
    LOG(4, "x");
    for (size_t i = 0; i < m_source_vector.size(); ++i) {
      for (size_t p = 0; p < m_source_vector[i].size(); ++p) {
        extended_vector[m_ext_indicies[i] + p] = m_source_vector[i][p];
        LOG(4, "row " << m_ext_indicies[i] + p << " value "
                      << extended_vector[m_ext_indicies[i] + p]);
      }
    }

    // W
    LOG(4, "W");
    for (size_t i = 0; i < m_W.size(); ++i) {
      for (int m = 0; m < vector_size; ++m) {
        extended_vector[size_x + i * vector_size + m] = m_W[i][m];
        LOG(4, "row " << size_x + i * vector_size + m << " value "
                      << extended_vector[size_x + i * vector_size + m]);
      }
    }

    // g
    LOG(4, "g");
    for (size_t i = 0; i < m_g.size(); ++i) {
      for (int m = 0; m < vector_size; ++m) {
        extended_vector[size_x + size_W + i * vector_size + m] = m_g[i][m];
        LOG(4,
            "row " << size_x + size_W + i * vector_size + m << " value "
                   << extended_vector[size_x + size_W + i * vector_size + m]);
      }
    }

    return extended_vector;
  }

  column_vector_type
  filter_extended_vector(const column_vector_type &extended_vector) const {
    const size_t size_x = m_col_particles->size();
    const size_t vector_size = m_vector_type::RowsAtCompileTime;
    const size_t size_W = m_W.size() * vector_size;
    const size_t size_g = m_g.size() * vector_size;
    const size_t n = size_x + size_W + size_g;

    column_vector_type filtered_vector(size_x);

    // x
    for (auto &bucket : m_query->get_subtree()) {
      if (m_query->is_leaf_node(bucket)) { // leaf node
        const size_t index = m_query->get_bucket_index(bucket);
        for (size_t i = 0; i < m_row_indices[index].size(); ++i) {
          filtered_vector[m_row_indices[index][i]] =
              extended_vector[m_ext_indicies[index] + i];
        }
      }
    }

    return filtered_vector;
  }

private:
  template <typename RowParticles>
  void generate_matrices(
      const child_iterator_vector_type &parents_strong_connections,
      const box_type &box_parent, const child_iterator &ci,
      const RowParticles &row_particles, const ColParticles &col_particles) {

    const box_type &target_box = m_query->get_bounds(ci);
    size_t target_index = m_query->get_bucket_index(*ci);
    LOG(3, "generate_matrices with bucket " << target_box);

    detail::theta_condition<dimension> theta(target_box.bmin, target_box.bmax);

    // add strongly connected buckets to current connectivity list
    if (parents_strong_connections.empty()) {
      for (child_iterator cj = m_query->get_children(); cj != false; ++cj) {
        const box_type &source_box = m_query->get_bounds(cj);
        if (theta.check(source_box.bmin, source_box.bmax)) {
          // add strongly connected buckets to current connectivity list
          m_strong_connectivity[target_index].push_back(cj);
        } else {
          // from weakly connected buckets,
          // add connectivity and generate m2l matricies
          size_t source_index = m_query->get_bucket_index(*cj);
          m_m2l_matrices[target_index].emplace_back();
          m_expansions.M2L_matrix(*(m_m2l_matrices[target_index].end() - 1),
                                  target_box, source_box);
          m_weak_connectivity[target_index].push_back(cj);
        }
      }
    } else {

      // transfer matrix with parent if not at start
      m_expansions.L2L_matrix(m_l2l_matrices[target_index], target_box,
                              box_parent);

      for (const child_iterator &source : parents_strong_connections) {
        if (m_query->is_leaf_node(*source)) {
          m_strong_connectivity[target_index].push_back(source);
        } else {
          for (child_iterator cj = m_query->get_children(source); cj != false;
               ++cj) {
            const box_type &source_box = m_query->get_bounds(cj);
            if (theta.check(source_box.bmin, source_box.bmax)) {
              m_strong_connectivity[target_index].push_back(cj);
            } else {
              size_t source_index = m_query->get_bucket_index(*cj);
              m_m2l_matrices[target_index].emplace_back();
              m_expansions.M2L_matrix(*(m_m2l_matrices[target_index].end() - 1),
                                      target_box, source_box);
              m_weak_connectivity[target_index].push_back(cj);
            }
          }
        }
      }
    }
    if (!m_query->is_leaf_node(*ci)) { // leaf node
      for (child_iterator cj = m_query->get_children(ci); cj != false; ++cj) {
        generate_matrices(m_strong_connectivity[target_index], target_box, cj,
                          row_particles, col_particles);
      }
    } else {
      m_expansions.P2M_matrix(m_p2m_matrices[target_index], target_box,
                              m_col_indices[target_index], col_particles);
      m_expansions.L2P_matrix(m_l2p_matrices[target_index], target_box,
                              m_row_indices[target_index], row_particles);

      child_iterator_vector_type strong_copy =
          m_strong_connectivity[target_index];
      m_strong_connectivity[target_index].clear();
      for (child_iterator &source : strong_copy) {
        if (m_query->is_leaf_node(*source)) {
          m_p2p_matrices[target_index].emplace_back();
          m_strong_connectivity[target_index].push_back(source);
          size_t source_index = m_query->get_bucket_index(*source);
          m_expansions.P2P_matrix(*(m_p2p_matrices[target_index].end() - 1),
                                  m_row_indices[target_index],
                                  m_col_indices[source_index], row_particles,
                                  col_particles);
        } else {
          auto range = m_query->get_subtree(source);
          for (all_iterator i = range.begin(); i != range.end(); ++i) {
            if (m_query->is_leaf_node(*i)) {
              m_strong_connectivity[target_index].push_back(
                  i.get_child_iterator());
              m_p2p_matrices[target_index].emplace_back();
              size_t index = m_query->get_bucket_index(*i);
              m_expansions.P2P_matrix(*(m_p2p_matrices[target_index].end() - 1),
                                      m_row_indices[target_index],
                                      m_col_indices[index], row_particles,
                                      col_particles);
            }
          }
        }
      }
    }
  }

  template <typename RowParticles>
  void generate_row_matrices(const child_iterator &ci,
                             const RowParticles &row_particles,
                             const ColParticles &col_particles) {

    const box_type &target_box = m_query->get_bounds(ci);
    size_t target_index = m_query->get_bucket_index(*ci);
    LOG(3, "generate_row_matrices with bucket " << target_box);

    if (!m_query->is_leaf_node(*ci)) { // leaf node
      for (child_iterator cj = m_query->get_children(ci); cj != false; ++cj) {
        generate_row_matrices(cj, row_particles, col_particles);
      }
    } else {
      m_expansions.L2P_matrix(m_l2p_matrices[target_index], target_box,
                              m_row_indices[target_index], row_particles);

      for (child_iterator &source : m_strong_connectivity[target_index]) {
        ASSERT(m_query->is_leaf_node(*source), "should be leaf node");
        m_p2p_matrices[target_index].emplace_back();
        size_t source_index = m_query->get_bucket_index(*source);
        m_expansions.P2P_matrix(*(m_p2p_matrices[target_index].end() - 1),
                                m_row_indices[target_index],
                                m_col_indices[source_index], row_particles,
                                col_particles);
      }
    }
  }

  void mvm_upward_sweep(typename col_tree_t::reverse_iterator ci) const {

    const size_t my_index = m_query->get_bucket_index(*ci);
    const box_type &target_box = m_query->get_bounds(ci);
    LOG(3, "calculate_dive_P2M_and_M2M with bucket " << target_box);
    m_vector_type &W = m_W[my_index];
    if (m_query->is_leaf_node(*ci)) { // leaf node
      W = m_p2m_matrices[my_index] * m_source_vector[my_index];
    } else {
      // do M2M
      W = m_vector_type::Zero();
      for (child_iterator cj = m_query->get_children(ci); cj != false; ++cj) {
        const size_t child_index = m_query->get_bucket_index(*cj);
        m_vector_type &child_W = mvm_upward_sweep(cj);
        W += m_l2l_matrices[child_index].transpose() * child_W;
      }
    }
    return W;
  }

  void mvm_downward_sweep(const m_vector_type &g_parent,
                          const child_iterator &ci) const {
    const box_type &target_box = m_query->get_bounds(ci);
    LOG(3, "calculate_dive_M2L_and_L2L with bucket " << target_box);
    size_t target_index = m_query->get_bucket_index(*ci);
    m_vector_type &g = m_g[target_index];

    // L2L
    g = m_l2l_matrices[target_index] * g_parent;

    // M2L (weakly connected buckets)
    for (size_t i = 0; i < m_weak_connectivity[target_index].size(); ++i) {
      const child_iterator &source_ci = m_weak_connectivity[target_index][i];
      size_t source_index = m_query->get_bucket_index(*source_ci);
      g += m_m2l_matrices[target_index][i] * m_W[source_index];
    }

    if (!m_query->is_leaf_node(*ci)) { // dive down to next level
      for (child_iterator cj = m_query->get_children(ci); cj != false; ++cj) {
        mvm_downward_sweep(g, cj);
      }
    } else {
      m_target_vector[target_index] = m_l2p_matrices[target_index] * g;

      // direct evaluation (strongly connected buckets)
      for (size_t i = 0; i < m_strong_connectivity[target_index].size(); ++i) {
        const child_iterator &source_ci =
            m_strong_connectivity[target_index][i];
        size_t source_index = m_query->get_bucket_index(*source_ci);
        m_target_vector[target_index] +=
            m_p2p_matrices[target_index][i] * m_source_vector[source_index];
      }
    }
  }
};

/*
template <typename Expansions, typename RowParticlesType, typename
ColParticlesType> H2Matrix<Expansions,ColParticlesType> make_h2_matrix(const
RowParticlesType& row_particles, const ColParticlesType& col_particles, const
Expansions& expansions) { return
H2Matrix<Expansions,ColParticlesType>(row_particles,col_particles,expansions);
}
*/

} // namespace Aboria

#endif
