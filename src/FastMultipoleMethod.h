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

#ifndef FAST_MULTIPOLE_METHOD_H_
#define FAST_MULTIPOLE_METHOD_H_

#include "detail/FastMultipoleMethod.h"
#ifdef HAVE_EIGEN
#include "detail/Kernels.h"
#endif

namespace Aboria {

template <typename Expansions, typename Kernel, typename RowParticles,
          typename ColParticles>
class FastMultipoleMethod {
protected:
  typedef typename RowParticles::query_type row_query_type;
  typedef typename ColParticles::query_type col_query_type;
  typedef typename row_query_type::reference row_reference;
  typedef typename row_query_type::pointer row_pointer;
  typedef typename row_query_type::child_iterator row_child_iterator;
  typedef typename row_query_type::particle_iterator row_particle_iterator;
  typedef typename row_particle_iterator::reference row_particle_reference;
  typedef typename col_query_type::reference col_reference;
  typedef typename col_query_type::pointer col_pointer;
  typedef typename col_query_type::child_iterator col_child_iterator;
  typedef typename col_query_type::particle_iterator col_particle_iterator;
  typedef typename col_particle_iterator::reference col_particle_reference;

  typedef typename Expansions::l_expansion_type l_expansion_type;
  typedef typename Expansions::m_expansion_type m_expansion_type;
  typedef typename ColParticles::traits_type traits_type;
  typedef typename traits_type::template vector_type<l_expansion_type>::type
      l_storage_type;
  typedef typename traits_type::template vector_type<m_expansion_type>::type
      m_storage_type;
  typedef typename traits_type::template vector_type<
      col_child_iterator
      // typename std::remove_const<child_iterator>::type
      >::type col_child_iterator_vector_type;
  typedef typename traits_type::template vector_type<
      col_child_iterator_vector_type>::type connectivity_type;
  typedef typename traits_type::double_d double_d;
  typedef typename traits_type::position position;
  static const unsigned int dimension = traits_type::dimension;
  typedef bbox<dimension> box_type;

  mutable m_storage_type m_W;
  mutable l_storage_type m_g;
  mutable connectivity_type m_connectivity;

  const ColParticles *m_col_particles;
  const RowParticles *m_row_particles;
  const col_query_type *m_col_query;
  const row_query_type *m_row_query;
  Expansions m_expansions;
  Kernel m_kernel;
  const int m_num_tasks;

public:
  FastMultipoleMethod(const RowParticles &row_particles,
                      const ColParticles &col_particles,
                      const Expansions &expansions, const Kernel &kernel)
      : m_col_particles(&col_particles), m_row_particles(&row_particles),
        m_col_query(&col_particles.get_query()),
        m_row_query(&row_particles.get_query()), m_expansions(expansions),
        m_kernel(kernel),
#ifdef HAVE_OPENMP
        m_num_tasks(omp_get_max_threads())
#else
        m_num_tasks(1)
#endif
  {
  }

  // target_vector += A*source_vector
  template <typename VectorTypeTarget, typename VectorTypeSource>
  void matrix_vector_multiply(VectorTypeTarget &target_vector,
                              const VectorTypeSource &source_vector) const {
    CHECK(target_vector.size() == source_vector.size(),
          "source and target vector not same length")
    m_W.resize(m_col_query->number_of_buckets());
    m_g.resize(m_row_query->number_of_buckets());
    m_connectivity.resize(m_row_query->number_of_buckets());

// upward sweep of tree
//
/*
#pragma omp parallel default(none)                                             \
    shared(source_vector, target_vector, m_num_tasks, m_W, m_g, m_col_query,   \
           m_row_query, m_expansions, m_kernel, m_connectivity)
           */
#pragma omp parallel shared(source_vector, target_vector)
    {
#pragma omp single
      {
        const int nchild_col = m_col_query->num_children();
        for (col_child_iterator ci =
                 m_col_particles->get_query().get_children();
             ci != false; ++ci) {
          /*
      #pragma omp task default(none) firstprivate(ci) \ shared(source_vector,
      m_num_tasks, m_W, m_col_query, m_expansions)
          */
#pragma omp task default(shared) firstprivate(ci) shared(source_vector)
          calculate_dive_P2M_and_M2M(ci, source_vector,
                                     m_num_tasks - nchild_col);
        }

#pragma omp taskwait

        // downward sweep of tree.
        //
        const int nchild_row = m_row_query->num_children();
        for (row_child_iterator ci = m_row_query->get_children(); ci != false;
             ++ci) {
          /*
#pragma omp task default(none) firstprivate(ci)                                \
shared(target_vector, source_vector, m_num_tasks, m_W, m_g, m_col_query,   \
      m_row_query, m_expansions, m_kernel, m_connectivity)
      */
#pragma omp task default(shared) firstprivate(ci)                              \
    shared(target_vector, source_vector)
          {
            col_child_iterator_vector_type dummy;
            l_expansion_type g{};
            calculate_dive_M2L_and_L2L(target_vector, dummy, g, box_type(), ci,
                                       source_vector, m_num_tasks - nchild_row);
          }
        }
#pragma omp taskwait
      }
    }
  }

private:
  template <typename VectorType>
  m_expansion_type &calculate_dive_P2M_and_M2M(const col_child_iterator &ci,
                                               const VectorType &source_vector,
                                               const int num_tasks) const {
    const size_t my_index = m_col_query->get_bucket_index(*ci);
    const box_type &my_box = m_col_query->get_bounds(ci);
    LOG(3, "calculate_dive_P2M_and_M2M with bucket " << my_box);
    m_expansion_type &W = m_W[my_index];
    typedef detail::VectorTraits<typename m_expansion_type::value_type>
        vector_traits;
    std::fill(std::begin(W), std::end(W), vector_traits::Zero());
    if (m_col_query->is_leaf_node(*ci)) { // leaf node
      detail::calculate_P2M(W, my_box, m_col_query->get_bucket_particles(*ci),
                            source_vector, m_col_query->get_particles_begin(),
                            m_expansions);
    } else {
      if (num_tasks > 0) {
        const int nchildren = m_col_query->num_children(ci);
        for (col_child_iterator cj = m_col_query->get_children(ci); cj != false;
             ++cj) {
          /*
#pragma omp task default(none) firstprivate(cj)                                \
shared(source_vector, W, my_box, m_W, m_col_query, m_expansions)
*/
#pragma omp task default(shared) firstprivate(cj)                              \
    shared(source_vector, W, my_box)
          {
            m_expansion_type &child_W = calculate_dive_P2M_and_M2M(
                cj, source_vector, num_tasks - nchildren);

            m_expansion_type localW{};
            const box_type &child_box = m_col_query->get_bounds(cj);
            m_expansions.M2M(localW, my_box, child_box, child_W);
            for (size_t i = 0; i < localW.size(); ++i) {
              detail::VectorTraits<typename m_expansion_type::value_type>::
                  AtomicIncrement(W[i], localW[i]);
            }
          }
        }
#pragma omp taskwait
      } else {
        for (col_child_iterator cj = m_col_query->get_children(ci); cj != false;
             ++cj) {
          m_expansion_type &child_W =
              calculate_dive_P2M_and_M2M(cj, source_vector, num_tasks);

          const box_type &child_box = m_col_query->get_bounds(cj);
          m_expansions.M2M(W, my_box, child_box, child_W);
        }
      }
    }
    return W;
  }

  /*
    template <typename VectorType>
    void calculate_dive_P2M_and_M2M(const tree_t &tree,
                                    const VectorType &source_vector) const {
      m_multipoles.resize(m_col_query.number_of_buckets());
      for (auto level = tree.rbegin(); level != --tree.rend(); ++level) {
        detail::for_each(level.begin(), level.end(), [](child_iterator ci) {
          const auto parent_box = m_col_query->get_parent_bounds(ci);
          const size_t parent_index = m_col_query->get_parent_index(ci);
          auto &parent_multipole = m_multipoles[parent_index];
          const auto box = m_col_query->get_bounds(ci);
          LOG(3, "calculate_dive_P2M_and_M2M with bucket " << box);
          const size_t index = m_col_query->get_bucket_index(*ci);
          auto &multipole = m_multipoles[index];
          if (m_col_query->is_leaf_node(*ci)) { // leaf node
            detail::calculate_P2M(
                multipole, box, m_col_query->get_bucket_particles(*ci),
                source_vector, m_col_query->get_particles_begin(),
  m_expansions);
          }
          m_expansions.M2M(parent_multipole, parent_box, box, multipole);
          }
      });
    }
  }

  template <typename VectorTypeTarget, typename VectorTypeSource>
  void calculate_dive_M2L_and_L2L(const tree_t &source_tree,
                                  const tree_t &target_tree,
                                  VectorTypeTarget &target_vector,
                                  const VectorTypeSource &source_vector) const {
    auto target_level = target_tree.begin();
    auto source_level = source_tree.begin();
    int2_vector_t current_level(1, vint4(0, 0, 0, 0));
    int2_vector_t next_level;

    for (; target_level != target_tree.end(); ++target_level, ++source_level) {
      // execute L2L and L2P on target_tree
      detail::for_each(
          ++target_level.begin(), target_level.end(), [](child_iterator ci) {
            const auto parent_box = m_col_query->get_parent_bounds(ci);
            const size_t parent_index = m_col_query->get_parent_index(ci);
            auto &parent_multipole = m_multipoles[parent_index];
            const auto box = m_col_query->get_bounds(ci);
            LOG(3, "calculate_L2L with bucket " << box);
            const size_t index = m_col_query->get_bucket_index(*ci);
            auto &multipole = m_multipoles[index];
            m_expansions.L2L(g, box, box_parent, g_parent);
            if (m_col_query->is_leaf_node(*ci)) { // leaf node
              detail::calculate_L2P(target_vector, local, box,
                                    m_row_query->get_bucket_particles(*ci),
                                    m_row_query->get_particles_begin(),
                                    m_expansions);
            }
          });

      // determine number of children ( P2P (leaf+leaf) M2L(node+node+theta) = 0
      // children, 1 leaf  = nc, 0 leaf =  nc*nc)

      num_children.resize(current_level.size());
      auto count_children = [](const vint4 &ij) {
        int num_children = 0;
        const auto &ci = target_tree[ij[0]][ij[1]];
        const auto &cj = source_tree[ij[2]][ij[3]];
        const bool ci_is_leaf = m_row_query->is_leaf_node(*ci);
        const bool cj_is_leaf = m_col_query->is_leaf_node(*cj);
        if (ci_is_leaf && cj_is_leaf) {
          return 0;
        } else if (detail::theta_condition < dimension()) {
          return 0;
        } else if (ci_is_leaf) {
          return m_col_query->number_of_children(cj);
        } else if (cj_is_leaf) {
          return m_row_query->number_of_children(ci);
        } else {
          return m_row_query->number_of_children(ci) *
                 m_col_query->number_of_children(cj);
        }
      };
      detail::transform(current_level.begin(), current_level.end(),
                        num_children.begin(), count_children);

      // partion pairs by number of children = 0
      auto ppoint = detail::partition(current_level.begin(),
  current_level.end(), num_children.begin(),
                                      [](const int nc) { return nc > 0; });

      // execute P2P and M2L ops (nc = 0)
      detail::for_each(ppoint, current_level.end(), []());

      // enumerate the number of children
      detail::exclusive_scan(num_children.begin(),num_children.end());

      // create next level

      // count number of children

      // create new level
      next_level.resize(num_children);
      detail::for_each(
          traits_t::make_zip_iterator(
              traits_t::make_tuple(current_level.begin(),
  num_children.begin())), traits_t::make_zip_iterator(
              traits_t::make_tuple(current_level.end(), num_children.end())),
          next_level.end(), [](auto i) {
        auto ij = i.template get<0>();
        auto ci = target_level[ij[0]];
        auto cj = source_level[ij[1]];
        int next_index = i.template get<1>();
        for (int ci_index = target_next_index[ij[0]]; ci != false;
             ++ci, ++ci_index) {

          size_t target_box = m_row_query->get_bounds(ci);
          detail::theta_condition<dimension> theta(target_box.bmin,
                                                   target_box.bmax);
          for (int cj_index = source_next_index[ij[1]]; cj != false;
               ++cj, ++cj_index) {

            size_t source_box = m_col_query->get_bounds(cj);
            if (theta.check(source_box.bmin, source_box.bmax)) {
              if (is_leaf(ci, cj)) {
                // do P2P or M2L
                m_expansions.M2L(g, target_box, source_box, m_W[source_index]);
              } else {
                _next_level[next_index++] = vint2(ci_index, cj_index);
              }
            }
          }
          return num_children;
        });

        // swap to current level
        current_level.swap(next_level);

        // expand parents
        const box_type &source_box = m_row_query->get_bounds(ci);
        LOG(3, "calculate_dive_M2L_and_L2L with bucket " << target_box);
        size_t target_index = m_row_query->get_bucket_index(*ci);
        l_expansion_type &g = m_g[target_index];
        typedef detail::VectorTraits<typename l_expansion_type::value_type>
            vector_traits;
        std::fill(std::begin(g), std::end(g), vector_traits::Zero());
        typename connectivity_type::reference connected_buckets =
            m_connectivity[target_index];
        connected_buckets.clear();

        if (connected_buckets_parent.empty()) {
          for (col_child_iterator cj = m_col_query->get_children(); cj != false;
               ++cj) {
            const box_type &source_box = m_col_query->get_bounds(cj);
            if (theta.check(source_box.bmin, source_box.bmax)) {
              connected_buckets.push_back(cj);
            } else {
              size_t source_index = m_col_query->get_bucket_index(*cj);
              m_expansions.M2L(g, target_box, source_box, m_W[source_index]);
            }
          }
        } else {
          // expansion from parent
          m_expansions.L2L(g, target_box, box_parent, g_parent);

          // expansions from weakly connected buckets on this level
          // and store strongly connected buckets to connectivity list
          for (const col_child_iterator &source : connected_buckets_parent) {
            if (m_col_query->is_leaf_node(*source)) {
              connected_buckets.push_back(source);
            } else {
              for (col_child_iterator cj = m_col_query->get_children(source);
                   cj != false; ++cj) {
                const box_type &source_box = m_col_query->get_bounds(cj);
                if (theta.check(source_box.bmin, source_box.bmax)) {
                  connected_buckets.push_back(cj);
                } else {
                  size_t source_index = m_col_query->get_bucket_index(*cj);
                  m_expansions.M2L(g, target_box, source_box,
  m_W[source_index]);
                }
              }
            }
          }
        }
        if (!m_row_query->is_leaf_node(*ci)) { // leaf node
          for (row_child_iterator cj = m_row_query->get_children(ci); cj !=
  false;
               ++cj) {
            calculate_dive_M2L_and_L2L(target_vector, connected_buckets, g,
                                       target_box, cj, source_vector);
          }
        } else if (target_vector.size() > 0) {
          detail::calculate_L2P(target_vector, g, target_box,
                                m_row_query->get_bucket_particles(*ci),
                                m_row_query->get_particles_begin(),
  m_expansions);

          for (col_child_iterator &cj : connected_buckets) {
            if (m_col_query->is_leaf_node(*cj)) {
              LOG(3, "calculate_P2P: target = " << target_box << " source = "
                                                << m_col_query->get_bounds(cj));
              detail::calculate_P2P(target_vector, source_vector,
                                    m_row_query->get_bucket_particles(*ci),
                                    m_col_query->get_bucket_particles(*cj),
                                    m_row_query->get_particles_begin(),
                                    m_col_query->get_particles_begin(),
  m_kernel); } else { for (auto j = m_col_query->get_subtree(cj); j != false;
  ++j) { if (m_col_query->is_leaf_node(*j)) {
                  detail::calculate_P2P(target_vector, source_vector,
                                        m_row_query->get_bucket_particles(*ci),
                                        m_col_query->get_bucket_particles(*j),
                                        m_row_query->get_particles_begin(),
                                        m_col_query->get_particles_begin(),
                                        m_kernel);
                }
              }
            }
          }
        }
        */

  template <typename VectorTypeTarget, typename VectorTypeSource>
  void calculate_dive_M2L_and_L2L(
      VectorTypeTarget &target_vector,
      const col_child_iterator_vector_type &connected_buckets_parent,
      const l_expansion_type &g_parent, const box_type &box_parent,
      const row_child_iterator &ci, const VectorTypeSource &source_vector,
      const int num_tasks) const {
    const box_type &target_box = m_row_query->get_bounds(ci);
    LOG(3, "calculate_dive_M2L_and_L2L with bucket " << target_box);
    size_t target_index = m_row_query->get_bucket_index(*ci);
    l_expansion_type &g = m_g[target_index];
    typedef detail::VectorTraits<typename l_expansion_type::value_type>
        vector_traits;
    std::fill(std::begin(g), std::end(g), vector_traits::Zero());
    typename connectivity_type::reference connected_buckets =
        m_connectivity[target_index];
    connected_buckets.clear();

    detail::theta_condition<dimension> theta(target_box.bmin, target_box.bmax);

    if (connected_buckets_parent.empty()) {
      for (col_child_iterator cj = m_col_query->get_children(); cj != false;
           ++cj) {
        const box_type &source_box = m_col_query->get_bounds(cj);
        if (theta.check(source_box.bmin, source_box.bmax)) {
          connected_buckets.push_back(cj);
        } else {
          size_t source_index = m_col_query->get_bucket_index(*cj);
          m_expansions.M2L(g, target_box, source_box, m_W[source_index]);
        }
      }
    } else {
      // expansion from parent
      m_expansions.L2L(g, target_box, box_parent, g_parent);

      // expansions from weakly connected buckets on this level
      // and store strongly connected buckets to connectivity list
      for (const col_child_iterator &source : connected_buckets_parent) {
        if (m_col_query->is_leaf_node(*source)) {
          connected_buckets.push_back(source);
        } else {
          for (col_child_iterator cj = m_col_query->get_children(source);
               cj != false; ++cj) {
            const box_type &source_box = m_col_query->get_bounds(cj);
            if (theta.check(source_box.bmin, source_box.bmax)) {
              connected_buckets.push_back(cj);
            } else {
              size_t source_index = m_col_query->get_bucket_index(*cj);
              m_expansions.M2L(g, target_box, source_box, m_W[source_index]);
            }
          }
        }
      }
    }
    if (!m_row_query->is_leaf_node(*ci)) { // leaf node
      if (num_tasks > 0) {
        const int nchildren = m_row_query->num_children(ci);
        for (row_child_iterator cj = m_row_query->get_children(ci); cj != false;
             ++cj) {
          /*
#pragma omp task default(none) firstprivate(cj) shared(                        \
target_vector, connected_buckets, g, target_box, source_vector, m_W, m_g,  \
m_col_query, m_row_query, m_expansions, m_kernel, m_connectivity)
*/
#pragma omp task default(shared) firstprivate(cj)                              \
    shared(target_vector, connected_buckets, g, target_box, source_vector)
          calculate_dive_M2L_and_L2L(target_vector, connected_buckets, g,
                                     target_box, cj, source_vector,
                                     num_tasks - nchildren);
        }
#pragma omp taskwait
      } else {
        for (row_child_iterator cj = m_row_query->get_children(ci); cj != false;
             ++cj) {
          calculate_dive_M2L_and_L2L(target_vector, connected_buckets, g,
                                     target_box, cj, source_vector, num_tasks);
        }
      }
    } else if (target_vector.size() > 0) {
      detail::calculate_L2P(target_vector, g, target_box,
                            m_row_query->get_bucket_particles(*ci),
                            m_row_query->get_particles_begin(), m_expansions);

      for (col_child_iterator &cj : connected_buckets) {
        if (m_col_query->is_leaf_node(*cj)) {
          LOG(3, "calculate_P2P: target = " << target_box << " source = "
                                            << m_col_query->get_bounds(cj));
          detail::calculate_P2P(target_vector, source_vector,
                                m_row_query->get_bucket_particles(*ci),
                                m_col_query->get_bucket_particles(*cj),
                                m_row_query->get_particles_begin(),
                                m_col_query->get_particles_begin(), m_kernel);
        } else {
          for (auto j = m_col_query->get_subtree(cj); j != false; ++j) {
            if (m_col_query->is_leaf_node(*j)) {
              detail::calculate_P2P(target_vector, source_vector,
                                    m_row_query->get_bucket_particles(*ci),
                                    m_col_query->get_bucket_particles(*j),
                                    m_row_query->get_particles_begin(),
                                    m_col_query->get_particles_begin(),
                                    m_kernel);
            }
          }
        }
      }
    }
  }
};

#ifdef HAVE_EIGEN
template <unsigned int D, unsigned int N, typename Function,
          typename KernelHelper = detail::position_kernel_helper<D, Function>,
          typename Block = typename KernelHelper::Block>
detail::BlackBoxExpansions<D, N, Function, Block::RowsAtCompileTime,
                           Block::ColsAtCompileTime>
make_black_box_expansion(const Function &function) {
  return detail::BlackBoxExpansions<D, N, Function, Block::RowsAtCompileTime,
                                    Block::ColsAtCompileTime>(function);
}
#else
template <unsigned int D, unsigned int N, typename Function>
detail::BlackBoxExpansions<D, N, Function, 1, 1>
make_black_box_expansion(const Function &function) {
  return detail::BlackBoxExpansions<D, N, Function, 1, 1>(function);
}
#endif

template <typename Expansions, typename Kernel, typename RowParticles,
          typename ColParticles>
FastMultipoleMethod<Expansions, Kernel, RowParticles, ColParticles>
make_fmm(const RowParticles &row_particles, const ColParticles &col_particles,
         const Expansions &expansions, const Kernel &kernel) {
  return FastMultipoleMethod<Expansions, Kernel, RowParticles, ColParticles>(
      row_particles, col_particles, expansions, kernel);
}

/*
template <typename Expansions, typename Kernel, typename ColParticles,
typename VectorType>
FastMultipoleMethodWithSource<Expansions,Kernel,ColParticles>
make_fmm_with_source(const ColParticles &col_particles,
                     const Expansions& expansions,
                     const Kernel& kernel,
                     const VectorType& source_vector) {
    return FastMultipoleMethodWithSource<Expansions,Kernel,ColParticles>
                            (col_particles,expansions,kernel,source_vector);
}
*/

} // namespace Aboria

#endif
