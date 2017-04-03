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

namespace Aboria {

template <typename Expansions, typename NeighbourQuery, 
          typename SourceVectorType>
struct calculate_P2M_and_M2M {
    typedef typename NeighbourQuery::traits_type traits_type;
    typedef typename traits_type::raw_pointer SourceParticleIterator;
    typedef typename Expansions::expansion_type expansion_type;
    typedef typename NeighbourQuery::reference reference;
    typedef typename NeighbourQuery::particle_iterator particle_iterator;
    typedef typename std::vector<expansion_type> StorageVectorType;
    static const unsigned int dimension = traits_type::dimension;
    typedef detail::bbox<dimension> box_type;


    const NeighbourQuery &m_query;
    const SourceVectorType &m_source_vector;
    StorageVectorType &m_W;
        
    calculate_P2M_and_M2M(const NeighbourQuery& query, 
                          const SourceVectorType& source_vector,
                          StorageVectorType& W):
        m_query(query),
        m_source_vector(source_vector),
        m_W(W)
    {}

    expansion_type& operator()(const reference bucket) {
        box_type my_box(m_query.get_bucket_bounds_low(bucket),
                            m_query.get_bucket_bounds_high(bucket));
        size_t my_index = m_query.get_bucket_index(bucket);
        expansion_type& W = m_W[my_index];

        if (m_query.is_leaf_node(bucket)) { // leaf node
            detail::calculate_P2M<Expansions>(W, my_box, 
                          m_query.get_bucket_particles(bucket),
                          m_source_vector,m_query.get_particles_begin());
        } else { // assume binary tree for now
            const reference child1 = *m_query.get_child1(&bucket);
            const reference child2 = *m_query.get_child2(&bucket);
            box_type child1_box(m_query.get_bucket_bounds_low(child1),
                                m_query.get_bucket_bounds_high(child1));
            box_type child2_box(m_query.get_bucket_bounds_low(child2),
                                m_query.get_bucket_bounds_high(child2));

            expansion_type& child1_W = this->operator()(child1);
            expansion_type& child2_W = this->operator()(child2);
            Expansions::M2M(W,my_box,child1_box,child1_W);
            Expansions::M2M(W,my_box,child1_box,child2_W);
        }
        return W;
    }
};

template <typename Expansions, typename NeighbourQuery, 
          typename SourceVectorType, typename Function>
struct calculate_M2L_and_L2L {
    typedef typename NeighbourQuery::reference reference;
    typedef typename NeighbourQuery::value_type* pointer;
    typedef typename NeighbourQuery::traits_type traits_type;
    typedef typename Expansions::expansion_type expansion_type;
    typedef typename NeighbourQuery::particle_iterator particle_iterator;
    typedef std::vector<expansion_type> StorageVectorType;
    static const unsigned int dimension = traits_type::dimension;
    typedef detail::bbox<dimension> box_type;
    typedef iterator_range<typename NeighbourQuery::root_iterator> root_iterator_range;

    const NeighbourQuery &m_query;
    StorageVectorType &m_g;
    StorageVectorType &m_W;
    Function &m_K;
    std::vector<std::vector<pointer>>& m_connectivity;

    calculate_M2L_and_L2L(const NeighbourQuery& query, 
                StorageVectorType &W,
                  StorageVectorType& g,
                  Function& K,
                  std::vector<std::vector<pointer>>& connectivity):
        m_query(query),
        m_W(W),
        m_g(g),
        m_K(K),
        m_connectivity(connectivity)
    {}

    void calculate_dive(const std::vector<pointer>& connected_buckets_parent,
                        const expansion_type& g_parent, 
                        const box_type& box_parent, 
                        const reference bucket) {
        box_type target_box(m_query.get_bucket_bounds_low(bucket),
                        m_query.get_bucket_bounds_high(bucket));
        size_t target_index = m_query.get_bucket_index(bucket);
        expansion_type& g = m_g[target_index];
        std::vector<pointer>& connected_buckets = m_connectivity[target_index];
        //connected_buckets.reserve(connected_buckets_parent.size());

        // expansion from parent
        Expansions::L2L(g,target_box,box_parent,g_parent);

        detail::theta_condition<dimension> theta(target_box.bmin,target_box.bmax);
        // expansions from weakly connected buckets on this level
        // and store strongly connected buckets to connectivity list
        for (const pointer source_pointer: connected_buckets_parent) {
            const reference source_bucket = *source_pointer;
            box_type source_box(m_query.get_bucket_bounds_low(source_bucket),
                                m_query.get_bucket_bounds_high(source_bucket));
            size_t source_index = m_query.get_bucket_index(source_bucket);
            if (theta.check(source_box.bmin,source_box.bmax)) {
                connected_buckets.push_back(&source_bucket);
            } else {
                Expansions::M2L(g,target_box,source_box,m_W[source_index],m_K);
            }
        }

        if (!m_query.is_leaf_node(bucket)) { // leaf node
            calculate_dive(connected_buckets,g,target_box,*m_query.get_child1(&bucket));
            calculate_dive(connected_buckets,g,target_box,*m_query.get_child2(&bucket));
        }
    }

    // only called for root nodes, uses the "calculate_L2L_dive" function
    // to recurse down the tree
    void operator()(const reference bucket) {
        // do a single-level FMM on the root nodes (no heirarchy above this)
        // TODO: can I use a fft to speed this up:
        // "A Matrix Version of the Fast Multipole Method"
        size_t target_index = m_query.get_bucket_index(bucket);
        box_type target_box(m_query.get_bucket_bounds_low(bucket),
                        m_query.get_bucket_bounds_high(bucket));
        detail::theta_condition<dimension> theta(target_box.bmin,target_box.bmax);
        root_iterator_range root_buckets = m_query.get_root_buckets();
        std::vector<pointer>& connected_buckets = m_connectivity[target_index];
        for (const reference source_bucket: root_buckets) {
            box_type source_box(m_query.get_bucket_bounds_low(source_bucket),
                                m_query.get_bucket_bounds_high(source_bucket));
            size_t source_index = m_query.get_bucket_index(source_bucket);
            if (theta.check(source_box.bmin,source_box.bmax)) {
                connected_buckets.push_back(&source_bucket);
            } else {
                Expansions::M2L(m_g[target_index],target_box,source_box,m_W[source_index],m_K);
            }
        }

        // now dive into the tree and do a proper FMM
        if (!m_query.is_leaf_node(bucket)) { 
            expansion_type& g = m_g[target_index];
            calculate_dive(connected_buckets,g,target_box,*m_query.get_child1(&bucket));
            calculate_dive(connected_buckets,g,target_box,*m_query.get_child2(&bucket));
        }
    }

};

   
template <typename Expansions, typename Function, typename NeighbourQuery>
class FastMultipoleMethod {
    typedef typename NeighbourQuery::traits_type traits_type;
    typedef typename NeighbourQuery::reference reference;
    typedef typename NeighbourQuery::value_type* pointer;
    typedef typename Expansions::expansion_type expansion_type;
    typedef typename traits_type::template vector_type<expansion_type>::type storage_type;
    typedef typename traits_type::template vector_type<pointer>::type bucket_pointer_vector_type;
    typedef typename traits_type::template vector_type<bucket_pointer_vector_type>::type connectivity_type;
    typedef iterator_range<typename NeighbourQuery::root_iterator> root_iterator_range;
    typedef typename NeighbourQuery::particle_iterator particle_iterator;
    typedef typename particle_iterator::reference particle_reference;
    static const unsigned int dimension = traits_type::dimension;
    typedef detail::bbox<dimension> box_type;
    storage_type m_W;
    storage_type m_g;
    connectivity_type m_connectivity; 
    const NeighbourQuery *m_query;
    Function m_K;

public:

    FastMultipoleMethod(const NeighbourQuery &query, const Function& K):
        m_query(&query),m_K(K)
    {}

    template <typename VectorType>
    void calculate_expansions(const VectorType& source_vector) {
        root_iterator_range root_buckets = m_query->get_root_buckets();

        m_W.resize(m_query->number_of_buckets());
        m_g.resize(m_query->number_of_buckets());

        // upward sweep of tree
        // calculate P2M and M2M expansions for source buckets (recursive up
        // the tree so use function object)
        std::for_each(root_buckets.begin(),
                      root_buckets.end(),
                      calculate_P2M_and_M2M<Expansions,NeighbourQuery,
                                            VectorType>(
                                                *m_query,
                                                source_vector,
                                                m_W));


        // downward sweep of tree
        // calculate L2L translations for source buckets (recursive 
        // down the tree so use function object)
        std::for_each(root_buckets.begin(),
                      root_buckets.end(),
                      calculate_M2L_and_L2L<Expansions,NeighbourQuery,
                                            VectorType,Function>(
                                                *m_query,
                                                m_W,
                                                m_g,
                                                m_K,
                                                m_connectivity));

    }


    // evaluate expansions for given point
    template <typename VectorType>
    double evaluate_expansion(const Vector<double,dimension>& p, const VectorType& source_vector) {
        const reference bucket = m_query->get_bucket(p);
        const size_t index = m_query->get_bucket_index(bucket); 
        box_type box(m_query->get_bucket_bounds_low(bucket),
                     m_query->get_bucket_bounds_high(bucket));

        double sum = Expansions::L2P(p,box,m_g[index]);
        for (pointer source_pointer: m_connectivity[index]) { 
            sum += detail::calculate_K_direct(p
                    ,m_query->get_bucket_particles(*source_pointer)
                    ,m_K,source_vector,m_query->get_particles_begin());
        }
        return sum;
    }

};

}

#endif

