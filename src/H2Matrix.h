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

template <typename Expansions, typename ColParticles>
class H2Matrix {
    typedef typename ColParticles::query_type query_type;
    typedef typename query_type::traits_type traits_type;
    typedef typename query_type::reference reference;
    typedef typename query_type::pointer pointer;
    typedef typename query_type::child_iterator child_iterator;
    typedef typename Expansions::expansion_type expansion_type;
    typedef typename Expansions::matrix_type matrix_type;
    typedef typename traits_type::template vector_type<expansion_type>::type storage_type;
    typedef typename traits_type::template vector_type<
        child_iterator
        >::type child_iterator_vector_type;
    typedef typename traits_type::template vector_type<
        child_iterator_vector_type
        >::type connectivity_type;
    typedef typename traits_type::template vector_type<matrix_type>::type vector_of_matrix_type;
    typedef typename traits_type::template vector_type<
       vector_of_matrix_type 
        >::type vector_of_vector_of_matrix_type;
    typedef typename Expansions::dynamic_vector_type dynamic_vector_type;
    typedef typename Expansions::l2p_matrix_type l2p_matrix_type;
    typedef typename traits_type::template vector_type<dynamic_vector_type>::type vector_of_vector_type;
    typedef typename traits_type::template vector_type<l2p_matrix_type>::type l2p_matrices_type;
    typedef typename traits_type::template vector_type<size_t>::type vector_of_size_t;
    typedef typename traits_type::template vector_type<vector_of_size_t>::type vector_of_vector_of_size_t;

    typedef typename query_type::particle_iterator particle_iterator;
    typedef typename particle_iterator::reference particle_reference;
    static const unsigned int dimension = traits_type::dimension;
    typedef detail::bbox<dimension> box_type;
    storage_type m_W;
    storage_type m_g;
    vector_of_vector_type m_source_vector;
    vector_of_vector_type m_target_vector;

    vector_of_dyn_matrix_type m_l2p_matrices;
    vector_of_dyn_matrix_type m_p2m_matrices;
    vector_of_vector_of_dyn_matrix_type m_k_direct_matrices;
    vector_of_matrix_type m_l2l_matrices;
    vector_of_vector_of_matrix_type m_m2l_matrices;
    vector_of_vector_of_size_t m_row_indices;
    vector_of_vector_of_size_t m_col_indices;

    connectivity_type m_strong_connectivity; 
    connectivity_type m_weak_connectivity; 

    Expansions m_expansions;

public:

    template <typename RowParticles>
    H2Matrix(const ColParticles &col_particles, const RowParticles &row_particles,const Expansions& expansions):
        m_col_particles(col_particles),m_expansions(expansions)
    {
        //generate h2 matrix 
        query_type& col_query = m_col_particles.get_query();
        const size_t n = query->number_of_buckets();
        const bool row_equals_col = &row_particles == &col_particles;
        m_W.resize(n);
        m_g.resize(n);
        m_l2l_matrices.resize(n);
        m_m2l_matricies.resize(n);
        m_row_indices.resize(n);
        m_col_indices.resize(n);
        m_l2p_matricies.resize(n);
        m_p2m_matricies.resize(n);
        m_strong_connectivity.resize(n);
        m_weak_connectivity.resize(n);

        // setup row and column indices
        if (!row_equals_col) {
            for (int i = 0; i < row_particles.size(); ++i) {
                //get target index
                pointer bucket;
                box_type box;
                m_query->get_bucket(p,bucket,box);
                const size_t index = m_query->get_bucket_index(*bucket); 
                m_row_indices[index].push_back(i);
            }
        }
        for (auto& bucket: m_query->get_subtree()) {
            if (m_query->is_leaf_node(*ci)) { // leaf node
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
            }
        }

        // downward sweep of tree to generate matrices
        for (child_iterator ci = m_query->get_children(); ci != false; ++ci) {
            const box_type& target_box = m_query->get_bounds(ci);
            LOG(3,"downward sweep with target bucket "<<target_box);
            size_t target_index = m_query->get_bucket_index(*ci);

            detail::theta_condition<dimension> theta(target_box.bmin,target_box.bmax);
            for (child_iterator cj = m_query->get_children(); cj != false; ++cj) {
                const box_type& source_box = m_query->get_bounds(cj);
                if (theta.check(source_box.bmin,source_box.bmax)) {
                    // add strongly connected buckets to current connectivity list
                    m_strong_connectivity[target_index].push_back(cj);
                } else {
                    // from weakly connected buckets, 
                    // add connectivity and generate m2l matricies
                    size_t source_index = m_query->get_bucket_index(*cj);
                    m_m2l_matrices[target_index].push_back(matrix_type());
                    m_expansions.M2L_matrix(
                            *(m_m2l_matrices[target_index].end()-1)
                            ,target_box,source_box);
                    m_weak_connectivity[target_index].push_back(cj);

                }
            }
            // if a leaf node then do dive
            if (!m_query->is_leaf_node(*ci)) { // leaf node
                for (child_iterator cj = m_query->get_children(ci); cj != false; ++cj) {
                    generate_matrices(target_index,target_box,cj
                        ,row_particles);
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
            }
        }
    }

    template <typename RowParticles>
    void generate_matrices(
            const size_t& parent_index,
            const box_type& box_parent, 
            const child_iterator& ci,
            const RowParticles &row_particles
            ) {

        const box_type& my_box = m_query->get_bounds(ci);
        size_t target_index = m_query->get_bucket_index(*ci);
        LOG(3,"generate_matrices with bucket "<<my_box);
        
        detail::theta_condition<dimension> theta(my_box.bmin,my_box.bmax);

        // transfer matrix with parent if not at start
        m_expansions.L2L_matrix(m_l2l_matrices[target_index],box_parent,my_box);

        // add strongly connected buckets to current connectivity list
        for (child_iterator& source: m_strong_connectivity[parent_index]) {
            if (m_query->is_leaf_node(*source)) {
                m_strong_connectivity[target_index].push_back(source);
            } else {
                for (child_iterator cj = m_query->get_children(source); cj != false; ++cj) {
                    const box_type& source_box = m_query->get_bounds(cj);
                    if (theta.check(source_box.bmin,source_box.bmax)) {
                        m_strong_connectivity[target_index].push_back(cj);
                    } else {
                        size_t source_index = m_query->get_bucket_index(*cj);
                        m_m2l_matrices[target_index].push_back(matrix_type());
                        m_expansions.M2L_matrix(
                                *(m_m2l_matrices[target_index].end()-1)
                                ,target_box,source_box);
                        m_weak_connectivity[target_index].push_back(cj);

                    }
                }
            }
        }
        if (!m_query->is_leaf_node(*ci)) { // leaf node
            for (child_iterator cj = m_query->get_children(ci); cj != false; ++cj) {
                generate_matrices(target_index,my_box,cj,row_particles);
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

            for (child_iterator& source: m_strong_connectivity[target_index]) {
                if (m_query->is_leaf_node(*source)) {
                    m_k_direct_matrices[target_index].push_back(dyn_matrix_type());
                    size_t source_index = m_query->get_bucket_index(*source);
                    m_expansions.K_direct_matrix(
                            *(m_k_direct_matrices[target_index].end()-1),
                            m_row_indices[target_index],m_col_indices[source_index]
                            row_particles,col_particles);
                } else {
                    for (reference bucket: m_query->get_subtree(source)) {
                        if (m_query->is_leaf_node(bucket)) {
                            size_t source_index = m_query->get_bucket_index(bucket);
                            m_expansions.K_direct_matrix(
                                *(m_k_direct_matrices[target_index].end()-1),
                                m_row_indices[target_index],m_col_indices[source_index]
                                row_particles,col_particles);
                        }
                    }
                }
            }
        }
    }

    expansion_type& calculate_dive_P2M_and_M2M(const child_iterator& ci) {
        const size_t my_index = m_query->get_bucket_index(*ci);
        const box_type& my_box = m_query->get_bounds(ci);
        LOG(3,"calculate_dive_P2M_and_M2M with bucket "<<my_box);
        expansion_type& W = m_W[my_index];
        if (m_query->is_leaf_node(*ci)) { // leaf node
            W = m_p2m_matrices[my_index]*m_source_vector[my_index];
        } else { 
            // do M2M 
            for (child_iterator cj = m_query->get_children(ci); cj != false; ++cj) {
                const size_t child_index = m_query->get_bucket_index(*cj);
                expansion_type& child_W = calculate_dive_P2M_and_M2M(cj,source_vector);
                W += m_transfer_matrices[child_index].transpose()*child_W;
            }
        }
        return W;
    }

    void calculate_dive_M2L_and_L2L(
            const typename connectivity_type::reference connected_buckets_parent,
            const expansion_type& g_parent, 
            const child_iterator& ci) {
        const box_type& target_box = m_query->get_bounds(ci);
        LOG(3,"calculate_dive_M2L_and_L2L with bucket "<<target_box);
        size_t target_index = m_query->get_bucket_index(*ci);
        expansion_type& g = m_g[target_index];
#ifndef NDEBUG
        for (int i = 0; i < g.size(); ++i) {
            ASSERT(std::abs(g[i]) < std::numeric_limits<double>::epsilon(),"g not zeroed");
        }
#endif
        // L2L
        g = m_l2l_matrices[target_index]*g_parent;

        // M2L
        for (int i = 0; i < m_weak_connectivity[target_box].size(); ++i) {
            size_t source_index = m_weak_connectivity[target_box][i];
            g += m_m2l_matrices[target_index][i]*m_W[source_index];
        }

        if (!m_query->is_leaf_node(*ci)) { // leaf node
            for (child_iterator cj = m_query->get_children(ci); cj != false; ++cj) {
                calculate_dive_M2L_and_L2L(connected_buckets,g,cj);
            }
        } else {
            m_target_vector[target_index] = m_l2p_matrices[target_index]*g;

            for (int i = 0; i < m_weak_connectivity[target_box].size(); ++i) {
                //***********************************
                size_t source_index = m_weak_connectivity[target_box][i];
                if (m_query->is_leaf_node(*cj)) {
                    size_t source_index = m_query->get_bucket_index(*cj);
                    m_target_vector[target_box] += m_k_direct_matrix[target_index][source_index]*m_source_vector[source_index];
                    sum += detail::calculate_K_direct(p
                        ,m_query->get_bucket_particles(*cj)
                        ,m_expansions,source_vector,m_query->get_particles_begin());
                } else {
                    for (reference subtree_reference: m_query->get_subtree(ci)) {
                        if (m_query->is_leaf_node(subtree_reference)) {
                            size_t source_index = m_query->get_bucket_index(*cj);
                            m_target_vector[target_box] += m_k_direct_matrix[source_index]*m_source_vector[source_index];
                            sum += detail::calculate_K_direct(p
                                    ,m_query->get_bucket_particles(subtree_reference)
                                    ,m_expansions,source_vector,m_query->get_particles_begin());
                        }
                    }
                }
            }
            //M2P???
        }

    }


    template <typename ColParticles, typename VectorType>
    void matrix_vector_multiply(const ColParticles& row_particles, 
            VectorType& target_vector, const VectorType& source_vector) {
        // upward sweep of tree
        //
        for (child_iterator ci = m_query->get_children(); ci != false; ++ci) {
            calculate_dive_P2M_and_M2M(ci,source_vector);
        }

        // downward sweep of tree. We treat the first layer specially since don't
        // have to do a L2L expansion from parents
        //
        // TODO: can I use a fft to speed this up, since it is on a regular grid?:
        // "A Matrix Version of the Fast Multipole Method"
        //
        for (child_iterator ci = m_query->get_children(); ci != false; ++ci) {
            const box_type& target_box = m_query->get_bounds(ci);
            LOG(3,"calculate_M2L with target bucket "<<target_box);
            size_t target_index = m_query->get_bucket_index(*ci);
            expansion_type& g = m_g[target_index];
#ifndef NDEBUG
            for (int i = 0; i < g.size(); ++i) {
                ASSERT(std::abs(g[i]) < std::numeric_limits<double>::epsilon(),"g not zeroed");
            }
#endif
            typename connectivity_type::reference 
                connected_buckets = m_connectivity[target_index];

            // do M2L transfers
            for (child_iterator& source: connected_buckets) {
                size_t source_index = m_query->get_bucket_index(*source);
                g = m_transfer_matrices[target_index][source_index]*m_W[source_index];
            }
         
            // if a leaf node then do dive
            if (!m_query->is_leaf_node(*ci)) { // leaf node
                for (child_iterator cj = m_query->get_children(ci); cj != false; ++cj) {
                    calculate_dive_M2L_and_L2L(connected_buckets,g,target_box,cj);
                }
            } else if (row_equals_col) {
            //M2P???
            }
        }

        if (!row_equals_col) {
            for (particle: row_particles) {
                target_vector[i] = matrix_vector_multiply_eval_point(get<position>(particle),source_vector);
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
        for (child_iterator& ci: m_connectivity[index]) { 
            if (m_query->is_leaf_node(*ci)) {
                sum += detail::calculate_K_direct(p
                    ,m_query->get_bucket_particles(*ci)
                    ,m_expansions,source_vector,m_query->get_particles_begin());
            } else {
                for (reference subtree_reference: m_query->get_subtree(ci)) {
                    if (m_query->is_leaf_node(subtree_reference)) {
                        sum += detail::calculate_K_direct(p
                                ,m_query->get_bucket_particles(subtree_reference)
                                ,m_expansions,source_vector,m_query->get_particles_begin());
                    }
                }
            }
        }
        return sum;
    }

};


template <unsigned int D, unsigned int N, typename Function> 
detail::BlackBoxExpansions<D,N,Function> make_black_box_expansion(const Function& function) {
    return detail::BlackBoxExpansions<D,N,Function>(function);
}


template <typename Expansions, typename NeighbourQuery>
FastMultipoleMethod<Expansions,NeighbourQuery>
make_fmm_query(const NeighbourQuery &query, const Expansions& expansions) {
    return FastMultipoleMethod<Expansions,NeighbourQuery>(query,expansions);
}

}

#endif

