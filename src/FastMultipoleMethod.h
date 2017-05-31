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
    typedef typename NeighbourQuery::child_iterator child_iterator;
    typedef typename NeighbourQuery::particle_iterator particle_iterator;
    typedef typename std::vector<expansion_type> StorageVectorType;
    static const unsigned int dimension = traits_type::dimension;
    typedef detail::bbox<dimension> box_type;


    const NeighbourQuery &m_query;
    const SourceVectorType &m_source_vector;
    StorageVectorType &m_W;
    Expansions &m_expansions;
        
    calculate_P2M_and_M2M(const NeighbourQuery& query, 
                          Expansions& expansions,
                          const SourceVectorType& source_vector,
                          StorageVectorType& W):
        m_query(query),
        m_source_vector(source_vector),
        m_W(W),
        m_expansions(expansions)
    {}

    expansion_type& calculate_dive(const child_iterator& it) {
        const size_t my_index = m_query.get_bucket_index(*it);
        const box_type& my_box = it.get_bounds();
        expansion_type& W = m_W[my_index];
        if (m_query.is_leaf_node(*it)) { // leaf node
            detail::calculate_P2M(W, my_box, 
                          m_query.get_bucket_particles(*it),
                          m_source_vector,m_query.get_particles_begin(),m_expansions);
        } else { 
            for (child_iterator ci = m_query->get_children(it); 
                    ci != false; ++ci) {
                expansion_type& child_W = calculate_dive(ci);
                const box_type& child_box = ci.get_bounds();
                m_expansions.M2M(W,my_box,child_box,child_W);
            }
        }
        return W;
    }
    expansion_type& operator()(reference bucket) {
        const box_type& my_box = m_query.get_root_bucket_bounds(bucket);
        size_t my_index = m_query.get_bucket_index(bucket);
        expansion_type& W = m_W[my_index];
#ifndef NDEBUG
        for (int i = 0; i < W.size(); ++i) {
            ASSERT(std::abs(W[i]) < std::numeric_limits<double>::epsilon(),"W not zeroed");
        }
#endif
        if (m_query.is_leaf_node(bucket)) { // leaf node
            detail::calculate_P2M(W, my_box, 
                          m_query.get_bucket_particles(bucket),
                          m_source_vector,m_query.get_particles_begin(),m_expansions);
        } else { 
            for (child_iterator ci = m_query->get_children(bucket,my_box); 
                    ci != false; ++ci) {
                expansion_type& child_W = calculate_dive(ci);
                const box_type& child_box = ci.get_bounds();
                m_expansions.M2M(W,my_box,child_box,child_W);
            }
        }
        return W;
    }
};

template <typename Expansions, typename NeighbourQuery, 
          typename SourceVectorType, typename ConnectivityType>
struct calculate_M2L_and_L2L {
    typedef typename NeighbourQuery::reference reference;
    typedef typename NeighbourQuery::pointer pointer;
    typedef typename NeighbourQuery::traits_type traits_type;
    typedef typename Expansions::expansion_type expansion_type;
    typedef typename NeighbourQuery::particle_iterator particle_iterator;
    typedef typename NeighbourQuery::child_iterator child_iterator;
    typedef std::vector<expansion_type> StorageVectorType;
    static const unsigned int dimension = traits_type::dimension;
    typedef detail::bbox<dimension> box_type;
    typedef iterator_range<typename NeighbourQuery::root_iterator> root_iterator_range;
    typedef std::queue<typename std::remove_const<pointer>::type> queue_type;

    const NeighbourQuery &m_query;
    StorageVectorType &m_g;
    const StorageVectorType &m_W;
    Expansions &m_expansions;
    ConnectivityType& m_connectivity;

    calculate_M2L_and_L2L(const NeighbourQuery& query, 
                  Expansions& expansions,
                  const StorageVectorType& W,
                  StorageVectorType& g,
                  ConnectivityType& connectivity):
        m_query(query),
        m_W(W),
        m_g(g),
        m_expansions(expansions),
        m_connectivity(connectivity)
    {}

    void calculate_dive(const typename ConnectivityType::reference connected_buckets_parent,
                        const expansion_type& g_parent, 
                        const box_type& box_parent, 
                        const child_iterator& it) {
        const box_type& target_box = it.get_bounds();
        size_t target_index = m_query.get_bucket_index(*it);
        expansion_type& g = m_g[target_index];
#ifndef NDEBUG
        for (int i = 0; i < g.size(); ++i) {
            ASSERT(std::abs(g[i]) < std::numeric_limits<double>::epsilon(),"g not zeroed");
        }
#endif
        typename ConnectivityType::reference connected_buckets = m_connectivity[target_index];
        ASSERT(connected_buckets.size() == 0,"con bucket not empty");
        //connected_buckets.reserve(connected_buckets_parent.size());

        // expansion from parent
        m_expansions.L2L(g,target_box,box_parent,g_parent);

        detail::theta_condition<dimension> theta(target_box.bmin,target_box.bmax);
        // expansions from weakly connected buckets on this level
        // and store strongly connected buckets to connectivity list
        for (child_iterator& source: connected_buckets_parent) {
            if (m_query.is_leaf_node(*source)) {
                connected_buckets.push_back(source);
            } else {
                for (child_iterator ci = m_query->get_children(source); 
                    ci != false; ++ci) {
                    const box_type& child_box = ci.get_bounds();
                    if (theta.check(child_box.bmin,child_box.bmax)) {
                        connected_buckets.push_back(ci);
                    } else {
                        //std::cout << "bucket from "<<child1_box.bmin<<" to "<<child1_box.bmax<<"is not connected to target box from "<<target_box.bmin<<" to "<<target_box.bmax<<std::endl;
                        size_t child1_index = m_query.get_bucket_index(*ci);
                        m_expansions.M2L(g,target_box,child_box,m_W[child1_index]);
                    }
                }
            }
        }
        if (!m_query.is_leaf_node(*it)) { // leaf node
            for (child_iterator& ci = m_query->get_children(it); 
                    ci != false; ++ci) {
                calculate_dive(connected_buckets,g,target_box,ci);
            }
        }
    }

    // only called for root nodes, uses the "calculate_L2L_dive" function
    // to recurse down the tree
    void operator()(reference bucket) {
        // do a single-level FMM on the root nodes (no heirarchy above this)
        // TODO: can I use a fft to speed this up, since it is on a regular grid?:
        // "A Matrix Version of the Fast Multipole Method"
        size_t target_index = m_query.get_bucket_index(bucket);
        const box_type& target_box = m_query.get_root_bucket_bounds(bucket);
        detail::theta_condition<dimension> theta(target_box.bmin,target_box.bmax);
        expansion_type& g = m_g[target_index];
        root_iterator_range root_buckets = m_query.get_root_buckets();
        typename ConnectivityType::reference connected_buckets = m_connectivity[target_index];
        for (reference source_bucket: root_buckets) {
            const box_type& source_box = m_query.get_root_bucket_bounds(source_bucket);
            size_t source_index = m_query.get_bucket_index(source_bucket);
            if (theta.check(source_box.bmin,source_box.bmax)) {
                connected_buckets.push_back(&source_bucket);
            } else {
                m_expansions.M2L(g,target_box,source_box,m_W[source_index]);
            }
        }

        // now dive into the tree and do a proper FMM
        if (!m_query.is_leaf_node(bucket)) { 
            for (child_iterator ci = m_query->get_children(bucket,target_box); 
                    ci != false; ++ci) {
                calculate_dive(connected_buckets,g,target_box,ci);
            }
        }
    }

};

   
template <typename Expansions, typename NeighbourQuery>
class FastMultipoleMethod {
    typedef typename NeighbourQuery::traits_type traits_type;
    typedef typename NeighbourQuery::reference reference;
    typedef typename NeighbourQuery::pointer pointer;
    typedef typename NeighbourQuery::child_iterator child_iterator;
    typedef typename Expansions::expansion_type expansion_type;
    typedef typename traits_type::template vector_type<expansion_type>::type storage_type;
    typedef typename traits_type::template vector_type<
        typename std::remove_reference<child_iterator>::type
        >::type child_iterator_vector_type;
    typedef typename traits_type::template vector_type<child_iterator_vector_type>::type connectivity_type;
    typedef iterator_range<typename NeighbourQuery::root_iterator> root_iterator_range;
    typedef typename NeighbourQuery::particle_iterator particle_iterator;
    typedef typename particle_iterator::reference particle_reference;
    static const unsigned int dimension = traits_type::dimension;
    typedef detail::bbox<dimension> box_type;
    storage_type m_W;
    storage_type m_g;
    connectivity_type m_connectivity; 
    const NeighbourQuery *m_query;
    Expansions m_expansions;

public:

    FastMultipoleMethod(const NeighbourQuery &query, const Expansions& expansions):
        m_query(&query),m_expansions(expansions)
    {}

    template <typename VectorType>
    void calculate_expansions(const VectorType& source_vector) {
        root_iterator_range root_buckets = m_query->get_root_buckets();

        const size_t n = m_query->number_of_buckets();
        m_W.resize(n);
        m_g.resize(n);
        m_connectivity.resize(n);

        // upward sweep of tree
        // calculate P2M and M2M expansions for source buckets (recursive up
        // the tree so use function object)
        std::for_each(root_buckets.begin(),
                      root_buckets.end(),
                      calculate_P2M_and_M2M<Expansions,NeighbourQuery,
                                            VectorType>(
                                                *m_query,
                                                m_expansions,
                                                source_vector,
                                                m_W
                                                ));


        // downward sweep of tree
        // calculate L2L translations for source buckets (recursive 
        // down the tree so use function object)
        std::for_each(root_buckets.begin(),
                      root_buckets.end(),
                      calculate_M2L_and_L2L<Expansions,NeighbourQuery,
                                            VectorType,connectivity_type>(
                                                *m_query,
                                                m_expansions,
                                                m_W,
                                                m_g,
                                                m_connectivity));

    }


    // evaluate expansions for given point
    template <typename VectorType>
    double evaluate_expansion(const Vector<double,dimension>& p, const VectorType& source_vector) {
        pointer bucket;
        box_type box;
        m_query->get_bucket(p,bucket,box);
        const size_t index = m_query->get_bucket_index(*bucket); 

        /*
        std::cout <<"(";
        for (int i = 0; i < m_g[index].size(); ++i) {
            std::cout <<m_g[index][i]<< ",";
        }
        std::cout <<")";
        */
            
        double sum = Expansions::L2P(p,box,m_g[index]);
        for (child_iterator& source_pointer: m_connectivity[index]) { 
            if (m_query->is_leaf_node(*source_pointer)) {
                sum += detail::calculate_K_direct(p
                    ,m_query->get_bucket_particles(*source_pointer)
                    ,m_expansions,source_vector,m_query->get_particles_begin());
            } else {
                for (reference subtree_reference: m_query->get_subtree(source_pointer)) {
                    if (m_query->is_leaf_node(*source_pointer)) {
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

