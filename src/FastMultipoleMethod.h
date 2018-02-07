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
class FastMultipoleMethodBase {
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
  typedef detail::bbox<dimension> box_type;

  mutable m_storage_type m_W;
  mutable l_storage_type m_g;
  mutable connectivity_type m_connectivity;

  const ColParticles *m_col_particles;
  const RowParticles *m_row_particles;
  const col_query_type *m_col_query;
  const row_query_type *m_row_query;
  Expansions m_expansions;
  Kernel m_kernel;

  FastMultipoleMethodBase(const RowParticles &row_particles,
                          const ColParticles &col_particles,
                          const Expansions &expansions, const Kernel &kernel)
      : m_col_particles(&col_particles), m_row_particles(&row_particles),
        m_col_query(&col_particles.get_query()),
        m_row_query(&row_particles.get_query()), m_expansions(expansions),
        m_kernel(kernel) {}

  template <typename VectorType>
  m_expansion_type &
  calculate_dive_P2M_and_M2M(const col_child_iterator &ci,
                             const VectorType &source_vector) const {
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
      for (col_child_iterator cj = m_col_query->get_children(ci); cj != false;
           ++cj) {
        m_expansion_type &child_W =
            calculate_dive_P2M_and_M2M(cj, source_vector);
        const box_type &child_box = m_col_query->get_bounds(cj);
        m_expansions.M2M(W, my_box, child_box, child_W);
      }
    }
    return W;
  }

  /*
    template <typename RowParticles>
    void generate_levels(
        const col_child_iterator_vector_type &parents_strong_connections,
        const box_type &box_parent, const child_iterator &ci_parent,
        const row_child_iterator &ci, const RowParticles &row_particles,
        const ColParticles &col_particles, const size_t level, bool active) {

      const box_type &target_box = m_query->get_bounds(ci);
      size_t target_index = m_query->get_bucket_index(*ci);
      LOG(3, "generate_levels with bucket " << target_box);

      // add current iterator to level
      if (level >= m_levels.size()) {
        m_levels.resize(level + 1);
      }

      // setup theta condition test
      detail::theta_condition<dimension> theta(target_box.bmin,
    target_box.bmax);

      // add strongly connected buckets to current connectivity list
      if (parents_strong_connections.empty()) {
        for (col_child_iterator cj = m_col_query->get_children(); cj != false;
             ++cj) {
          const box_type &source_box = m_query->get_bounds(cj);
          if (theta.check(source_box.bmin, source_box.bmax)) {
            // add strongly connected buckets to current connectivity list
            m_strong_connectivity[target_index].push_back(cj);
          } else {
            // add weakly connected buckets to current connectivity list
            m_weak_connectivity[target_index].push_back(cj);
            active = true;
          }
        }
      } else {
        if (active) {
          m_parent_connectivity[target_index] = ci_parent;
        }
        for (const col_child_iterator &source : parents_strong_connections) {
          if (m_col_query->is_leaf_node(*source)) {
            m_strong_connectivity[target_index].push_back(source);
          } else {
            for (col_child_iterator cj = m_col_query->get_children(source);
                 cj != false; ++cj) {
              const box_type &source_box = m_col_query->get_bounds(cj);
              if (theta.check(source_box.bmin, source_box.bmax)) {
                m_strong_connectivity[target_index].push_back(cj);
              } else {
                m_weak_connectivity[target_index].push_back(cj);
                active = true;
              }
            }
          }
        }
      }
      if (!m_row_query->is_leaf_node(*ci)) { // non leaf node
        for (child_iterator cj = m_query->get_children(ci); cj != false; ++cj) {
          generate_levels(m_strong_connectivity[target_index], target_box, ci,
    cj, row_particles, col_particles, level + 1, active);
        }
      } else {
        ASSERT(active == true,
               "have a non-active leaf node, don't know what to do in this
    case"); child_iterator_vector_type strong_copy =
            m_strong_connectivity[target_index];
        m_strong_connectivity[target_index].clear();
        for (child_iterator &source : strong_copy) {
          if (m_query->is_leaf_node(*source)) {
            m_strong_connectivity[target_index].push_back(source);
          } else {
            auto range = m_query->get_subtree(source);
            for (all_iterator i = range.begin(); i != range.end(); ++i) {
              if (m_query->is_leaf_node(*i)) {
                m_strong_connectivity[target_index].push_back(
                    i.get_child_iterator());
              }
            }
          }
        }
      }

      if (active) {
        m_levels[level].push_back(ci);
      }
    }
    */

  template <typename VectorTypeTarget, typename VectorTypeSource>
  void calculate_dive_M2L_and_L2L(
      VectorTypeTarget &target_vector,
      const col_child_iterator_vector_type &connected_buckets_parent,
      const l_expansion_type &g_parent, const box_type &box_parent,
      const row_child_iterator &ci,
      const VectorTypeSource &source_vector) const {
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
      for (row_child_iterator cj = m_row_query->get_children(ci); cj != false;
           ++cj) {
        calculate_dive_M2L_and_L2L(target_vector, connected_buckets, g,
                                   target_box, cj, source_vector);
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
          for (col_reference ref_j : m_col_query->get_subtree(cj)) {
            if (m_col_query->is_leaf_node(ref_j)) {
              detail::calculate_P2P(target_vector, source_vector,
                                    m_row_query->get_bucket_particles(*ci),
                                    m_col_query->get_bucket_particles(ref_j),
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

template <typename Expansions, typename Kernel, typename RowParticles,
          typename ColParticles>
class FastMultipoleMethod
    : public FastMultipoleMethodBase<Expansions, Kernel, RowParticles,
                                     ColParticles> {
  typedef FastMultipoleMethodBase<Expansions, Kernel, RowParticles,
                                  ColParticles>
      base_type;
  typedef typename base_type::row_child_iterator row_child_iterator;
  typedef typename base_type::col_child_iterator col_child_iterator;
  typedef typename base_type::col_child_iterator_vector_type
      col_child_iterator_vector_type;
  typedef typename base_type::l_expansion_type l_expansion_type;
  typedef typename base_type::m_expansion_type m_expansion_type;
  typedef typename base_type::box_type box_type;
  typedef typename base_type::double_d double_d;
  typedef typename base_type::position position;
  static const unsigned int dimension = base_type::dimension;

public:
  FastMultipoleMethod(const RowParticles &row_particles,
                      const ColParticles &col_particles,
                      const Expansions &expansions, const Kernel &kernel)
      : base_type(row_particles, col_particles, expansions, kernel) {}

  // target_vector += A*source_vector
  template <typename VectorTypeTarget, typename VectorTypeSource>
  void matrix_vector_multiply(VectorTypeTarget &target_vector,
                              const VectorTypeSource &source_vector) const {
    CHECK(target_vector.size() == source_vector.size(),
          "source and target vector not same length")
    this->m_W.resize(this->m_col_query->number_of_buckets());
    this->m_g.resize(this->m_row_query->number_of_buckets());
    this->m_connectivity.resize(this->m_row_query->number_of_buckets());

    // upward sweep of tree
    //
    for (col_child_iterator ci =
             this->m_col_particles->get_query().get_children();
         ci != false; ++ci) {
      this->calculate_dive_P2M_and_M2M(ci, source_vector);
    }

    // downward sweep of tree.
    //
    for (row_child_iterator ci = this->m_row_query->get_children(); ci != false;
         ++ci) {
      col_child_iterator_vector_type dummy;
      l_expansion_type g = {};
      this->calculate_dive_M2L_and_L2L(target_vector, dummy, g, box_type(), ci,
                                       source_vector);
    }
    /*
    else {
        for (child_iterator ci = this->m_query->get_children(); ci != false;
    ++ci) { child_iterator_vector_type dummy; std::vector<double> dummy2;
            expansion_type g = {};
            this->calculate_dive_M2L_and_L2L(dummy2,dummy,g,box_type(),ci,source_vector);
        }

        for (size_t i = 0; i < row_particles.size(); ++i) {
            const double_d& p = get<position>(row_particles)[i];
            pointer bucket;
            box_type box;
            this->m_query->get_bucket(p,bucket,box);
            LOG(3,"evaluating expansion at point "<<p<<" with box "<<box);
            const size_t index = this->m_query->get_bucket_index(*bucket);

            double sum = Expansions::L2P(p,box,this->m_g[index]);
            for (child_iterator& ci: this->m_connectivity[index]) {
                if (this->m_query->is_leaf_node(*ci)) {
                    sum += detail::calculate_P2P_position(p
                        ,this->m_query->get_bucket_particles(*ci)
                        ,this->m_expansions,source_vector,this->m_query->get_particles_begin());
                } else {
                    for (reference subtree_reference:
    this->m_query->get_subtree(ci)) { if
    (this->m_query->is_leaf_node(subtree_reference)) { sum +=
    detail::calculate_P2P_position(p
                                    ,this->m_query->get_bucket_particles(subtree_reference)
                                    ,this->m_expansions,source_vector,this->m_query->get_particles_begin());
                        }
                    }
                }
            }
            target_vector[i] += sum;
        }
    }
        */
  }
};

/*
template <typename Expansions, typename Kernel, typename RowParticles,
typename ColParticles> class FastMultipoleMethodWithSource: public
FastMultipoleMethodBase<Expansions,Kernel,RowParticles,ColParticles> { typedef
FastMultipoleMethodBase<Expansions,Kernel,RowParticles,ColParticles>
base_type; typedef typename base_type::row_child_iterator row_child_iterator;
    typedef typename base_type::col_child_iterator col_child_iterator;
    typedef typename base_type::col_child_iterator_vector_type
col_child_iterator_vector_type; typedef typename base_type::expansion_type
expansion_type; typedef typename base_type::box_type box_type; typedef
typename base_type::double_d double_d; typedef typename base_type::position
position; static const unsigned int dimension = base_type::dimension; public:
    template <typename VectorType>
    FastMultipoleMethodWithSource(const RowParticles& col_particles,
                        const ColParticles& col_particles,
                        const Expansions& expansions,
                        const Kernel& kernel,
                        const VectorType& source_vector):
        base_type(row_particles,col_particles,expansions,kernel)
    {
        this->m_W.resize(this->m_col_query->number_of_buckets());
        this->m_g.resize(this->m_row_query->number_of_buckets());
        this->m_connectivity.resize(this->m_row_query->number_of_buckets());

        // upward sweep of tree
        //
        for (col_child_iterator ci = this->m_col_query->get_children(); ci !=
false; ++ci) { this->calculate_dive_P2M_and_M2M(ci,source_vector);
        }

        // downward sweep of tree.
        //
        for (row_child_iterator ci = this->m_row_query->get_children(); ci !=
false; ++ci) { col_child_iterator_vector_type dummy; VectorType dummy2;
            expansion_type g = {};
            this->calculate_dive_M2L_and_L2L(dummy2,dummy,g,box_type(),ci,source_vector);
        }
    }


    // evaluate expansions for given point
    template <typename ParticleRef, typename VectorType>
    double evaluate_at_particle(const ParticleRef& p, const VectorType&
source_vector) const { pointer bucket; box_type box;
        this->m_query->get_bucket(p,bucket,box);
        LOG(3,"evaluating expansion at point "<<p<<" with box "<<box);
        const size_t index = this->m_query->get_bucket_index(*bucket);

        double sum = Expansions::L2P(p,box,this->m_g[index]);
        for (child_iterator& ci: this->m_connectivity[index]) {
            if (this->m_query->is_leaf_node(*ci)) {
                sum += detail::calculate_P2P_position(p
                    ,this->m_query->get_bucket_particles(*ci)
                    ,this->m_kernel,source_vector,this->m_query->get_particles_begin());
            } else {
                for (reference subtree_reference:
this->m_query->get_subtree(ci)) { if
(this->m_query->is_leaf_node(subtree_reference)) { sum +=
detail::calculate_P2P_position(p
                                ,this->m_query->get_bucket_particles(subtree_reference)
                                ,this->m_kernel,source_vector,this->m_query->get_particles_begin());
                    }
                }
            }
        }
        return sum;
    }

};
*/

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
