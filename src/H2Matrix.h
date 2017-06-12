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

#include "detail/FastMultipoleMethod.h"

namespace Aboria {

template <typename Expansions, typename Query>
class H2Matrix {
    typedef typename Query::traits_type traits_type;
    typedef typename Query::reference reference;
    typedef typename Query::pointer pointer;
    static const unsigned int dimension = Query::dimension;
    typedef position_d<dimension> position;
    typedef typename Query::child_iterator child_iterator;
    typedef typename Query::all_iterator all_iterator;
    typedef typename traits_type::template vector_type<
        child_iterator
        >::type child_iterator_vector_type;
    typedef typename traits_type::template vector_type<
        child_iterator_vector_type
        >::type connectivity_type;
    
    typedef typename Expansions::p_vector_type p_vector_type;
    typedef typename Expansions::m_vector_type m_vector_type;
    typedef typename Expansions::l2p_matrix_type l2p_matrix_type;
    typedef typename Expansions::p2m_matrix_type p2m_matrix_type;
    typedef typename Expansions::p2p_matrix_type p2p_matrix_type;
    typedef typename Expansions::l2l_matrix_type l2l_matrix_type;
    typedef typename Expansions::m2l_matrix_type m2l_matrix_type;
    typedef typename traits_type::template vector_type<p_vector_type>::type p_vectors_type;
    typedef typename traits_type::template vector_type<m_vector_type>::type m_vectors_type;
    typedef typename traits_type::template vector_type<l2p_matrix_type>::type l2p_matrices_type;
    typedef typename traits_type::template vector_type<p2m_matrix_type>::type p2m_matrices_type;
    typedef typename traits_type::template vector_type<p2p_matrix_type>::type p2p_matrices_type;
    typedef typename traits_type::template vector_type<p2p_matrices_type>::type vector_of_p2p_matrices_type;
    typedef typename traits_type::template vector_type<l2l_matrix_type>::type l2l_matrices_type;
    typedef typename traits_type::template vector_type<m2l_matrix_type>::type m2l_matrices_type;
    typedef typename traits_type::template vector_type<m2l_matrices_type>::type vector_of_m2l_matrices_type;
    typedef typename traits_type::template vector_type<size_t>::type indices_type;
    typedef typename traits_type::template vector_type<indices_type>::type vector_of_indices_type;

    typedef typename Query::particle_iterator particle_iterator;
    typedef typename particle_iterator::reference particle_reference;
    typedef detail::bbox<dimension> box_type;
    m_vectors_type m_W;
    m_vectors_type m_g;
    p_vectors_type m_source_vector;
    p_vectors_type m_target_vector;

    l2p_matrices_type m_l2p_matrices;
    p2m_matrices_type m_p2m_matrices;
    l2l_matrices_type m_l2l_matrices;
    vector_of_p2p_matrices_type m_p2p_matrices;
    vector_of_m2l_matrices_type m_m2l_matrices;
    vector_of_indices_type m_row_indices;
    vector_of_indices_type m_col_indices;

    connectivity_type m_strong_connectivity; 
    connectivity_type m_weak_connectivity; 

    Expansions m_expansions;

    const Query* m_query;

public:

    template <typename RowParticles, typename ColParticles>
    H2Matrix(const RowParticles &row_particles, const ColParticles &col_particles, const Expansions& expansions):
        m_query(&col_particles.get_query()),m_expansions(expansions)
    {
        //generate h2 matrix 
        const size_t n = m_query->number_of_buckets();
        const bool row_equals_col = &row_particles == &col_particles;
        m_W.resize(n);
        m_g.resize(n);
        m_l2l_matrices.resize(n);
        m_m2l_matrices.resize(n);
        m_row_indices.resize(n);
        m_col_indices.resize(n);
        m_source_vector.resize(n);
        m_target_vector.resize(n);
        m_l2p_matrices.resize(n);
        m_p2m_matrices.resize(n);
        m_p2p_matrices.resize(n);
        m_strong_connectivity.resize(n);
        m_weak_connectivity.resize(n);

        // setup row and column indices
        if (!row_equals_col) {
            for (int i = 0; i < row_particles.size(); ++i) {
                //get target index
                pointer bucket;
                box_type box;
                m_query->get_bucket(get<position>(row_particles)[i],bucket,box);
                const size_t index = m_query->get_bucket_index(*bucket); 
                m_row_indices[index].push_back(i);
            }
        }
        for (auto& bucket: m_query->get_subtree()) {
            if (m_query->is_leaf_node(bucket)) { // leaf node
                const size_t index = m_query->get_bucket_index(bucket); 
                //get target_index
                for (auto& p: m_query->get_bucket_particles(bucket)) {
                    const size_t i = &get<position>(p)
                                   - &get<position>(col_particles)[0];
                    m_col_indices[index].push_back(i);
                    if (row_equals_col) {
                        m_row_indices[index].push_back(i);
                    }
                }

                m_source_vector[index].resize(m_col_indices[index].size());
                m_target_vector[index].resize(m_row_indices[index].size());
            }
        }

        // downward sweep of tree to generate matrices
        for (child_iterator ci = m_query->get_children(); ci != false; ++ci) {
            const box_type& target_box = m_query->get_bounds(ci);
            generate_matrices(child_iterator_vector_type(),box_type(),ci,row_particles,col_particles);
        }
    }


    // target_vector += A*source_vector
    template <typename VectorType>
    void matrix_vector_multiply(VectorType& target_vector, const VectorType& source_vector) {

        // for all leaf nodes setup source vector
        for (auto& bucket: m_query->get_subtree()) {
            if (m_query->is_leaf_node(bucket)) { // leaf node
                const size_t index = m_query->get_bucket_index(bucket); 
                for (int i = 0; i < m_col_indices[index].size(); ++i) {
                    m_source_vector[index][i] = source_vector[m_col_indices[index][i]];
                }
            }
        }

        // upward sweep of tree
        for (child_iterator ci = m_query->get_children(); ci != false; ++ci) {
            mvm_upward_sweep(ci);
        }

        // downward sweep of tree. 
        for (child_iterator ci = m_query->get_children(); ci != false; ++ci) {
            m_vector_type g = m_vector_type::Zero();
            mvm_downward_sweep(g,ci);
        }

        // for all leaf nodes copy to target vector
        for (auto& bucket: m_query->get_subtree()) {
            if (m_query->is_leaf_node(bucket)) { // leaf node
                const size_t index = m_query->get_bucket_index(bucket); 
                for (int i = 0; i < m_row_indices[index].size(); ++i) {
                    target_vector[m_row_indices[index][i]] += m_target_vector[index][i];
                }
            }
        }
    }

    // evaluate expansions for given point
    template <typename VectorType>
    double matrix_vector_multiply_eval_point(const Vector<double,dimension>& p, const VectorType& source_vector) {
        pointer bucket;
        box_type box;
        m_query->get_bucket(p,bucket,box);
        LOG(3,"evaluating expansion at point "<<p<<" with box "<<box);
        const size_t index = m_query->get_bucket_index(*bucket); 

        double sum = Expansions::L2P(p,box,m_g[index]);
        // direct evaluation (strongly connected buckets)
        for (int i = 0; i < m_strong_connectivity[index].size(); ++i) {
            child_iterator& source_ci = m_strong_connectivity[index][i];
            size_t source_index = m_query->get_bucket_index(*source_ci);
            sum += detail::calculate_P2P_position(p
                    ,m_query->get_bucket_particles(*source_ci)
                    ,m_expansions,source_vector,m_query->get_particles_begin());
        }
        return sum;
    }

private:
    template <typename RowParticles, typename ColParticles>
    void generate_matrices(
            const child_iterator_vector_type& parents_strong_connections,
            const box_type& box_parent, 
            const child_iterator& ci,
            const RowParticles &row_particles,
            const ColParticles &col_particles
            ) {

        const box_type& target_box = m_query->get_bounds(ci);
        size_t target_index = m_query->get_bucket_index(*ci);
        LOG(3,"generate_matrices with bucket "<<target_box);
        
        detail::theta_condition<dimension> theta(target_box.bmin,target_box.bmax);

        // add strongly connected buckets to current connectivity list
        if (parents_strong_connections.empty()) {
            for (child_iterator cj = m_query->get_children(); cj != false; ++cj) {
                const box_type& source_box = m_query->get_bounds(cj);
                if (theta.check(source_box.bmin,source_box.bmax)) {
                    // add strongly connected buckets to current connectivity list
                    m_strong_connectivity[target_index].push_back(cj);
                } else {
                    // from weakly connected buckets, 
                    // add connectivity and generate m2l matricies
                    size_t source_index = m_query->get_bucket_index(*cj);
                    m_m2l_matrices[target_index].push_back(m2l_matrix_type());
                    m_expansions.M2L_matrix(
                            *(m_m2l_matrices[target_index].end()-1)
                            ,target_box,source_box);
                    m_weak_connectivity[target_index].push_back(cj);

                }
            }
        } else {

            // transfer matrix with parent if not at start
            m_expansions.L2L_matrix(m_l2l_matrices[target_index],target_box,box_parent);

            for (const child_iterator& source: parents_strong_connections) {
                if (m_query->is_leaf_node(*source)) {
                    m_strong_connectivity[target_index].push_back(source);
                } else {
                    for (child_iterator cj = m_query->get_children(source); cj != false; ++cj) {
                        const box_type& source_box = m_query->get_bounds(cj);
                        if (theta.check(source_box.bmin,source_box.bmax)) {
                            m_strong_connectivity[target_index].push_back(cj);
                        } else {
                            size_t source_index = m_query->get_bucket_index(*cj);
                            m_m2l_matrices[target_index].push_back(m2l_matrix_type());
                            m_expansions.M2L_matrix(
                                    *(m_m2l_matrices[target_index].end()-1)
                                    ,target_box,source_box);
                            m_weak_connectivity[target_index].push_back(cj);

                        }
                    }
                }
            }
        }
        if (!m_query->is_leaf_node(*ci)) { // leaf node
            for (child_iterator cj = m_query->get_children(ci); cj != false; ++cj) {
                generate_matrices(m_strong_connectivity[target_index]
                                    ,target_box,cj,row_particles,col_particles);
            }
        } else {
            m_expansions.P2M_matrix(m_p2m_matrices[target_index], 
                    target_box,
                    m_col_indices[target_index],
                    col_particles);
            m_expansions.L2P_matrix(m_l2p_matrices[target_index], 
                    target_box,
                    m_row_indices[target_index],
                    row_particles);

            
            child_iterator_vector_type strong_copy = m_strong_connectivity[target_index];
            m_strong_connectivity[target_index].clear();
            for (child_iterator& source: strong_copy) {
                if (m_query->is_leaf_node(*source)) {
                    m_p2p_matrices[target_index].push_back(p2p_matrix_type());
                    m_strong_connectivity[target_index].push_back(source);
                    size_t source_index = m_query->get_bucket_index(*source);
                    m_expansions.P2P_matrix(
                            *(m_p2p_matrices[target_index].end()-1),
                            m_row_indices[target_index],m_col_indices[source_index],
                            row_particles,col_particles);
                } else {
                    auto range = m_query->get_subtree(source);
                    for (all_iterator i = range.begin(); i!=range.end(); ++i) {
                        if (m_query->is_leaf_node(*i)) {
                            m_strong_connectivity[target_index].push_back(i.get_child_iterator());
                            m_p2p_matrices[target_index].push_back(p2p_matrix_type());
                            size_t index = m_query->get_bucket_index(*i);
                            m_expansions.P2P_matrix(
                                *(m_p2p_matrices[target_index].end()-1),
                                m_row_indices[target_index],m_col_indices[index],
                                row_particles,col_particles);
                        }
                    }
                }
            }
        }
    }

    m_vector_type& mvm_upward_sweep(const child_iterator& ci) {
        const size_t my_index = m_query->get_bucket_index(*ci);
        const box_type& target_box = m_query->get_bounds(ci);
        LOG(3,"calculate_dive_P2M_and_M2M with bucket "<<target_box);
        m_vector_type& W = m_W[my_index];
        if (m_query->is_leaf_node(*ci)) { // leaf node
            W = m_p2m_matrices[my_index]*m_source_vector[my_index];
        } else { 
            // do M2M 
            W = m_vector_type::Zero();
            for (child_iterator cj = m_query->get_children(ci); cj != false; ++cj) {
                const size_t child_index = m_query->get_bucket_index(*cj);
                m_vector_type& child_W = mvm_upward_sweep(cj);
                W += m_l2l_matrices[child_index].transpose()*child_W;
            }
        }
        return W;
    }

    void mvm_downward_sweep(
            const m_vector_type& g_parent, 
            const child_iterator& ci) {
        const box_type& target_box = m_query->get_bounds(ci);
        LOG(3,"calculate_dive_M2L_and_L2L with bucket "<<target_box);
        size_t target_index = m_query->get_bucket_index(*ci);
        m_vector_type& g = m_g[target_index];

        // L2L
        g = m_l2l_matrices[target_index]*g_parent;

        // M2L (weakly connected buckets)
        for (int i = 0; i < m_weak_connectivity[target_index].size(); ++i) {
            child_iterator& source_ci = m_weak_connectivity[target_index][i];
            size_t source_index = m_query->get_bucket_index(*source_ci);
            g += m_m2l_matrices[target_index][i]*m_W[source_index];
        }

        if (!m_query->is_leaf_node(*ci)) { // dive down to next level
            for (child_iterator cj = m_query->get_children(ci); cj != false; ++cj) {
                mvm_downward_sweep(g,cj);
            }
        } else {
            m_target_vector[target_index] = m_l2p_matrices[target_index]*g;

            // direct evaluation (strongly connected buckets)
            for (int i = 0; i < m_strong_connectivity[target_index].size(); ++i) {
                child_iterator& source_ci = m_strong_connectivity[target_index][i];
                size_t source_index = m_query->get_bucket_index(*source_ci);
                m_target_vector[target_index] += 
                        m_p2p_matrices[target_index][i]*m_source_vector[source_index];
            }
        }

    }


};


template <typename Expansions, typename RowParticlesType, typename ColParticlesType,
         typename NeighbourQuery = typename ColParticlesType::query_type>
H2Matrix<Expansions,NeighbourQuery>
make_h2_matrix(const RowParticlesType& row_particles, const ColParticlesType& col_particles, const Expansions& expansions) {
    return H2Matrix<Expansions,NeighbourQuery>(row_particles,col_particles,expansions);
}

}

#endif

