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

   
template <typename Expansions, typename RowParticles, typename ColParticles>
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
    typedef typename traits_type::template vector_type<expansion_type>::type matrix_storage_type;
    typedef typename traits_type::template vector_type<
        matrix_storage_type
        >::type vector_matrix_storage_type;

    typedef typename query_type::particle_iterator particle_iterator;
    typedef typename particle_iterator::reference particle_reference;
    static const unsigned int dimension = traits_type::dimension;
    typedef detail::bbox<dimension> box_type;
    storage_type m_W;
    storage_type m_g;
    matrix_storage_type m_transfer_matrices;
    vector_matrix_storage_type m_m2l_matrices;
    connectivity_type m_connectivity; 
    const RowParticles& m_row_particles; 
    const ColParticles& m_col_particles; 
    Expansions m_expansions;

public:

    H2Matrix(const RowParticles &row_particles, const ColParticles &col_particles, const Expansions& expansions):
        m_row_particles(&row_particles),m_col_particles(col_particles),m_expansions(expansions)
    {
        //generate h2 matrix 
        query_type& row_query = m_row_particles.get_query();
        const size_t n = query->number_of_buckets();
        m_W.resize(n);
        m_g.resize(n);
        m_transfer_matrices.resize(n);
        m_m2l_matricies.resize(n);
        m_connectivity.resize(n);

        // downward sweep of tree to generate matrices
        for (child_iterator ci = m_query->get_children(); ci != false; ++ci) {
            const box_type& target_box = m_query->get_bounds(ci);
            LOG(3,"downward sweep with target bucket "<<target_box);
            size_t target_index = m_query->get_bucket_index(*ci);

            typename connectivity_type::reference 
                connected_buckets = m_connectivity[target_index];
            typename vector_matrix_storage_type::reference
                connected_matrices = m_m2l_matrices[target_index];
            typename vector_matrix_storage_type::reference
                transfer_matrix = m_transfer_matrices[target_index];
            ASSERT(connected_buckets.size() == 0,"root bucket not empty");
            ASSERT(connected_matrices.size() == 0,"root bucket not empty");

            detail::theta_condition<dimension> theta(target_box.bmin,target_box.bmax);

            // add strongly connected buckets to current connectivity list
            for (child_iterator cj = m_query->get_children(); cj != false; ++cj) {
                const box_type& source_box = m_query->get_bounds(cj);
                if (theta.check(source_box.bmin,source_box.bmax)) {
                    connected_buckets.push_back(cj);
                }
            }
            // if a leaf node then do dive
            if (!m_query->is_leaf_node(*ci)) { // leaf node
                for (child_iterator cj = m_query->get_children(ci); cj != false; ++cj) {
                    generate_matrices(connected_buckets,target_box,cj);
                }
            }

            connected_buckets.empty();
            // from weakly connected buckets, add connectivity and generate m2l matricies
            for (child_iterator cj = m_query->get_children(); cj != false; ++cj) {
                const box_type& source_box = m_query->get_bounds(cj);
                if (!theta.check(source_box.bmin,source_box.bmax)) {
                    size_t source_index = m_query->get_bucket_index(*cj);
                    connected_matricies.push_back(matrix_type());
                    m_expansions.M2L_matrix(*(connected_matrices.end()-1),target_box,source_box);
                    connected_buckets.push_back(cj);
                }
            }
        }
    }

    void generate_matrices(
            const typename connectivity_type::reference connected_buckets_parent,
            const box_type& box_parent, 
            const child_iterator& ci) {

        const box_type& my_box = m_query->get_bounds(ci);
        LOG(3,"generate_matrices with bucket "<<my_box);
        size_t target_index = m_query->get_bucket_index(*ci);
        typename connectivity_type::reference 
            connected_buckets = m_connectivity[target_index];
        typename vector_matrix_storage_type::reference
            connected_matrices = m_m2l_matrices[target_index];
        typename vector_matrix_storage_type::reference
            transfer_matrix = m_transfer_matrices[target_index];
        ASSERT(connected_buckets.size() == 0,"con bucket not empty");
        ASSERT(connected_matrices.size() == 0,"root bucket not empty");

        detail::theta_condition<dimension> theta(my_box.bmin,my_box.bmax);

        // transfer matrix with parent
        m_expansions.M2M_matrix(transfer_matrix,box_parent,my_box);

        // add strongly connected buckets to current connectivity list
        for (child_iterator& source: connected_buckets_parent) {
            if (m_query->is_leaf_node(*source)) {
                connected_buckets.push_back(source);
            } else {
                for (child_iterator cj = m_query->get_children(source); cj != false; ++cj) {
                    const box_type& source_box = m_query->get_bounds(cj);
                    if (theta.check(source_box.bmin,source_box.bmax)) {
                        connected_buckets.push_back(cj);
                    }
                }
            }
        }
        if (!m_query->is_leaf_node(*ci)) { // leaf node
            for (child_iterator cj = m_query->get_children(ci); cj != false; ++cj) {
                generate_matrices(connected_buckets,my_box,cj);
            }
        }

        connected_buckets.empty();
        // from weakly connected buckets, add connectivity and generate m2l matricies
        for (child_iterator& source: connected_buckets_parent) {
            if (!m_query->is_leaf_node(*source)) {
                for (child_iterator cj = m_query->get_children(source); cj != false; ++cj) {
                    const box_type& source_box = m_query->get_bounds(cj);
                    if (!theta.check(source_box.bmin,source_box.bmax)) {
                        size_t source_index = m_query->get_bucket_index(*cj);
                        connected_matricies.push_back(matrix_type());
                        m_expansions.M2L_matrix(*(connected_matrices.end()-1),target_box,source_box);
                        connected_buckets.push_back(cj);
                    }
                }
            }
        }
    }

    template <typename VectorType>
    expansion_type& calculate_dive_P2M_and_M2M(const child_iterator& ci, 
                                               const VectorType& source_vector) {
        const size_t my_index = m_query->get_bucket_index(*ci);
        const box_type& my_box = m_query->get_bounds(ci);
        LOG(3,"calculate_dive_P2M_and_M2M with bucket "<<my_box);
        expansion_type& W = m_W[my_index];
        if (m_query->is_leaf_node(*ci)) { // leaf node
            detail::calculate_P2M(W, my_box, 
                    m_query->get_bucket_particles(*ci),
                    source_vector,m_query->get_particles_begin(),m_expansions);
        } else { 
            // do M2M 
            for (child_iterator cj = m_query->get_children(ci); cj != false; ++cj) {
                const size_t child_index = m_query->get_bucket_index(*cj);
                expansion_type& child_W = calculate_dive_P2M_and_M2M(cj,source_vector);
                const box_type& child_box = m_query->get_bounds(cj);
                W += m_transfer_matrices[child_index]*child_W;
            }
        }
        return W;
    }

    void calculate_dive_M2L_and_L2L(
            const typename connectivity_type::reference connected_buckets_parent,
            const expansion_type& g_parent, 
            const box_type& box_parent, 
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
        typename connectivity_type::reference 
            connected_buckets = m_connectivity[target_index];
        ASSERT(connected_buckets.size() == 0,"con bucket not empty");
        //connected_buckets.reserve(connected_buckets_parent.size());

        // L2L
        g = m_transfer_matrices[target_box].transpose()*g_parent;

        // M2L
        for (child_iterator& source: connected_buckets) {
            size_t source_index = m_query->get_bucket_index(*cj);
            g += m_m2l_matrices[target_index][source_index]*m_W[source_index];
        }

        if (!m_query->is_leaf_node(*ci)) { // leaf node
            for (child_iterator cj = m_query->get_children(ci); cj != false; ++cj) {
                calculate_dive_M2L_and_L2L(connected_buckets,g,target_box,cj);
            }
        } else if (row_equals_col) {
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
                const box_type& source_box = m_query->get_bounds(source);
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

