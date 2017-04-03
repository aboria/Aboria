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

#include "detail/Algorithms.h"
#include "detail/SpatialUtil.h"
#include "NeighbourSearchBase.h"
#include "Traits.h"
#include "CudaInclude.h"
#include "Vector.h"
#include "Get.h"

#include <iostream>
#include <unordered_map>
#include "Log.h"

namespace Aboria {

namespace detail {

    template <typename InOutIterator, unsigned int D, unsigned int N, typename PositionIterator, typename InputIterator> 
    struct BlackBoxExpansions {
        typedef detail::bbox<D> box_type;
        const size_t ncheb = std::pow(N,D); 
        typedef std::array<double,ncheb> vector_type;
        typedef Vector<double,D> double_d;


        void P2M(vector_type accum, 
                 box_type& box, 
                 const double_d& position,
                 const double& source ) {

            detail::ChebyshevRnSingle<D,N> cheb_rn(position,box);
            lattice_iterator<dimension> mj(m_start,m_end,m_start);
            for (int j=0; j<ncheb; ++j,++mj) {
                accum_begin[j] += cheb_rn(*mj)*source;
            }
        }

        void M2M(vector_type& accum, 
                 const box_type& target_box, 
                 const box_type& source_box, 
                 const vector_type& source) {

            int_d start(0);
            int_d end(N-1);
            std::array<double_d,ncheb> source_positions;

            lattice_iterator<dimension> mi(m_start,m_end,m_start);
            for (int i=0; i<ncheb; ++i,++mi) {
                const double_d pi_unit_box = detail::chebyshev_node_nd(*mi);
                source_positions[i] = 0.5*(pi_unit_box+1)*(source_box.bmax-source_box.bmin) 
                                        + source_box.bmin;
            }
             
            cheb_rn.calculate_Sn_with_bbox(source_positions.begin(),target_box,ncheb,N);

            lattice_iterator<D> mi(m_start,m_end,m_start);
            for (int i=0; i<ncheb; ++i,++mi) {
                for (int j=0; j<ncheb; ++j) {
                    //TODO: calculate cheb funs in here with ChebSnSingle?
                    accum[i] += cheb_rn(*mi,j)*source[j];
                }
            }
        }


        template <typename Function>
        void M2L(vector_type& accum, 
                 const box_type& target_box, 
                 const box_type& source_box, 
                 const vector_type& source,
                 Function& K) {

            int_d start(0);
            int_d end(N-1);

            lattice_iterator<dimension> mi(m_start,m_end,m_start);
            for (int i=0; i<ncheb; ++i,++mj) {
                const double_d pi_unit_box = detail::chebyshev_node_nd(*mi);
                const double_d pi = 0.5*(pi_unit_box+1)*(target_box.bmax-target_box.bmin) 
                                                                    + target_box.bmin;

                lattice_iterator<dimension> mj(m_start,m_end,m_start);
                for (int j=0; j<ncheb; ++j,++mj) {
                    const double_d pj_unit_box = detail::chebyshev_node_nd(*mj);
                    const double_d pj = 0.5*(pj_unit_box+1)*(source_box.bmax-source_box.bmin) 
                                                                    + source_box.bmin;
                    accum[i] += K(pj-pi,pi,pj)*source[j];
                }
            }


        }

        void L2L(vector_type& accum, 
                 const box_type& target_box, 
                 const box_type& source_box, 
                 const vector_type& source) {
            M2M(accum,target_box,source_box,source);
        }

        double L2P(const double_d& p,
                   const box_type& box, 
                   const vector_type& source) {
            detail::ChebyshevRnSingle cheb_rn(p,box);
            lattice_iterator<dimension> mj(m_start,m_end,m_start);
            double sum = 0;
            for (int j=0; j<ncheb; ++j,++mj) {
                sum += cheb_rn(*mj)*source[j];
            }
            return sum;
        }

    };


    template <typename Traits, 
              typename SourceVectorType, 
              typename Expansions,
                    typename VectorType=Expansions::vector_type,
                    typename SourceParticleIterator=Traits::raw_pointer, 
                    unsigned int D=Traits::dimension>
    void calculate_P2M(VectorType& sum, 
                        const detail::bbox<D>& box, 
                        const iterator_range<ranges_iterator<Traits>>& range, 
                        const SourceVectorType& source_vector,
                        const SourceParticleIterator& source_particles_begin) {
        typedef Traits::position position;
        const size_t N = std::distance(range.begin(),range.end());
        const double_d* pbegin = &get<position>(*range.begin());
        const size_t index = pbegin - &get<position>(source_particles_begin)[0];
        for (int i = 0; i < N; ++i) {
            const Vector<double,D>& pi = pbegin[i]; 
            Expansions::P2M(sum,box,pi,source_vector[index+i]);  
        }
    }

    // assume serial processing of particles, this could be more efficient for ranges iterators
    template <typename Iterator, 
              typename SourceVectorType, 
              typename Expansions,
                 typename VectorType=Expansions::vector_type, 
                 typename Traits=Iterator::traits_type,
                 typename SourceParticleIterator=Traits::raw_pointer, 
                 unsigned int D=Traits::dimension>
                 typename = typename
        std::enable_if<!std::is_same<Iterator,ranges_iterator<traits_type>>>
    void calculate_P2M(VectorType& sum, 
                        detail::bbox<D>& my_box, 
                        iterator_range<Iterator>& range, 
                        const SourceVectorType& source_vector,
                        const SourceParticleIterator& source_particles_begin) {

        for (particle_reference i: range) {
            const Vector<double,D>& pi = get<position>(i); 
            const size_t index = &pi- &get<position>(source_particles_begin)[0];
            Expansions::P2M(sum,box,pi,source_vector[index]);  
        }

    }


    template <typename Traits, 
              typename Function, 
              typename SourceVectorType,
                    typename SourceParticleIterator=Traits::raw_pointer, 
                    unsigned int D=Traits::dimension>
    double calculate_K_direct(const Vector<double,D>& position,
                            iterator_range<ranges_iterator<Traits>>& range, 
                            const Function& m_K,
                            const SourceVectorType& source_vector,
                            const SourceParticleIterator& source_particles_begin) {
        typedef Traits::position position;
        typedef Vector<double,D> double_d;
        const size_t N = std::distance(range.begin(),range.end());
        const double_d* pbegin = &get<position>(*range.begin());
        const size_t index = pbegin - &get<position>(source_particles_begin)[0];
        double sum = 0;
        for (int i = 0; i < N; ++i) {
            const Vector<double,D>& pi = pbegin[i]; 
            sum += m_K(pi-position,position,pi)*source_vector[index+i];
        }
        return sum;
    }

    template <typename Iterator, 
              typename Function, 
              typename SourceVectorType,
                    typename Traits=Iterator::traits_type,
                    typename SourceParticleIterator=Traits::raw_pointer, 
                    unsigned int D=Traits::dimension>
    double calculate_K_direct(const Vector<double,D>& position,
                            iterator_range<Iterator>& range, 
                            const Function& m_K,
                            const SourceVectorType& source_vector,
                            const SourceParticleIterator& source_particles_begin) {
        typedef Traits::position position;
        typedef Iterator::reference particle_reference;
        typedef Vector<double,D> double_d;
        
        double sum = 0;
        for (particle_reference i: range) {
            const Vector<double,D>& pi = get<position>(i); 
            const size_t index = &pi- &get<position>(source_particles_begin)[0];
            sum += m_K(pi-position,position,pi)*source_vector[index];
        }
        return sum;
    }

    template <unsigned int D>
    struct theta_condition {
        typedef Vector<double,D> double_d;
        const double_d& m_low;
        const double_d& m_high;
        const double m_r2;
        const double m_r;
        static const double m_theta = 0.5;
        constexpr double m_theta2 = m_theta*m_theta;
        theta_condition(const double_d& low, const double_d& high):
            m_low(low),m_high(high),
            m_r2(0.25*(high-low).squaredNorm()),
            m_r(std::sqrt(m_r2))
        {}

        bool theta_condition(const double_d& low, const double_d& high) {
        double d = 0.5*(high + low - m_low - m_high).norm(); 
        double other_r2 = 0.25*(high-low).squaredNorm();
        if (other_r2 < m_r2) {
            const double other_r = std::sqrt(other_r2);
            return m_r2 > m_theta2*std::pow(d-other_r,2);
        } else {
            return other_r2 < m_theta2*std::pow(d-m_r,2);
        }
    }
}
            


template <typename Expansions, typename NeighbourQuery, 
          typename SourceVectorType>
struct calculate_P2M_and_M2M {
    typedef NeighbourQuery::traits_type traits_type;
    typedef traits_type::raw_pointer SourceParticleIterator, 
    typedef Expansions::vector_type vector_type;
    typedef typename NeighbourQuery::reference reference;
    typedef typename NeighbourQuery::particle_iterator particle_iterator;
    typedef std::vector<vector_type> StorageVectorType;
    typedef detail::bbox<dimension> box_type;

    static const unsigned int dimension = traits_type::dimension;

    const NeighbourQuery &m_search;
    const SourceVectorType &m_source_vector;
    const SourceParticleIterator &m_source_particles_begin;
    StorageVectorType &m_W;
        
    calculate_S_expansion(const SourceNeighbourQuery& search, 
                          const SourceVectorType& source_vector,
                          StorageVectorType& W):
        m_search(search),
        m_source_vector(source_vector),
        m_W(W)
    {}

    vector_type& operator()(const reference bucket) {
        box_type my_box(m_search.get_bucket_bounds_low(bucket),
                            m_search.get_bucket_bounds_high(bucket));
        size_t my_index = m_search.get_bucket_index(bucket);
        vector_type& W = m_W[my_index];

        if (search.is_leaf_node(bucket) { // leaf node
            calculate_P2M(W, my_box, 
                          search.get_bucket_particles(bucket),
                          m_source_vector,m_source_particles_begin);
        } else { // assume binary tree for now
            const reference child1 = search.get_child1(bucket);
            const reference child2 = search.get_child2(bucket);
            box_type child1_box(search.get_bucket_bounds_low(child1),
                                search.get_bucket_bounds_high(child1));
            box_type child2_box(search.get_bucket_bounds_low(child2),
                                search.get_bucket_bounds_high(child2));

            vector_type& child1_W = this->operator()(child1);
            vector_type& child2_W = this->operator()(child2);
            Expansions::M2M(W,my_box,child1_box,child1_W);
            Expansions::M2M(W,my_box,child1_box,child2_W);
        }
        return W;
    }
};

template <typename Expansions, typename NeighbourQuery, 
          typename SourceVectorType>
struct calculate_M2L_and_L2L {
    typedef typename NeighbourQuery::reference reference;
    typedef typename NeighbourQuery::pointer pointer;
    typedef NeighbourQuery::traits_type traits_type;
    typedef Expansions::vector_type vector_type;
    typedef typename NeighbourQuery::reference reference;
    typedef typename NeighbourQuery::particle_iterator particle_iterator;
    typedef std::vector<vector_type> StorageVectorType;
    static const unsigned int dimension = traits_type::dimension;
    typedef detail::bbox<dimension> box_type;

    const NeighbourQuery &m_search;
    StorageVectorType &m_g;
    StorageVectorType &m_W;
    std::vector<std::vector<pointer>>& m_connectivity;

    calculate_M2L_and_L2L(const NeighbourQuery& search, 
                StorageVectorType &W,
                  StorageVectorType& g,
                  std::vector<std::vector<pointer>>& connectivity):
        m_search(search),
        m_W(W),
        m_g(g),
        m_connectivity(connectivity)
    {}

    void calculate_dive(std::vector<pointer>& connected_buckets_parent,
                            const vector_type& g_parent, 
                            const box_type& box_parent, 
                            const reference bucket) {
        box_type target_box(m_search.get_bucket_bounds_low(bucket),
                        m_search.get_bucket_bounds_high(bucket));
        size_t target_index = m_search.get_bucket_index(bucket);
        vector_type& g = m_g[index];
        std::vector<pointer>& connected_buckets = m_connectivity[target_index];
        //connected_buckets.reserve(connected_buckets_parent.size());

        // expansion from parent
        Expansions::L2L(g,box,box_parent,g_parent);

        detail::theta_condition<dimension> theta(target_box.bmin,target_box.bmax);
        // expansions from weakly connected buckets on this level
        // and store strongly connected buckets to connectivity list
        for (const pointer source_pointer: connected_buckets_parent) {
            const reference source_bucket = *source_pointer;
            box_type source_box(m_query->get_bucket_bounds_low(source_bucket),
                                m_query->get_bucket_bounds_high(source_bucket));
            size_t source_index = m_search.get_bucket_index(other_bucket);
            if (theta(source_box.bmin,source_box.bmax)) {
                connected_buckets.push_back(&source_bucket);
            } else {
                Expansions::M2L(g,target_box,source_box,m_W[source_index],m_K)
            }
        }

        if (!search.is_leaf_node(bucket) { // leaf node
            calculate_dive(g,box,m_search.get_child1(bucket));
            calculate_dive(g,box,m_search.get_child2(bucket));
        }
    }

    // only called for root nodes, uses the "calculate_L2L_dive" function
    // to recurse down the tree
    void operator()(const reference bucket) {
        // do a single-level FMM on the root nodes (no heirarchy above this)
        size_t target_index = m_search.get_bucket_index(bucket);
        box_type target_box(m_search.get_bucket_bounds_low(bucket),
                        m_search.get_bucket_bounds_high(bucket));
        detail::theta_condition<dimension> theta(target_box.bmin,target_box.bmax);
        root_iterator_range_type root_buckets = m_query->get_root_buckets();
        std::vector<pointer>& connected_buckets = m_connectivity[target_index];
        for (const reference source_bucket: root_buckets) {
            box_type source_box(m_query->get_bucket_bounds_low(source_bucket),
                                m_query->get_bucket_bounds_high(source_bucket));
            size_t source_index = m_search.get_bucket_index(other_bucket);
            if (theta(source_box.bmin,source_box.bmax)) {
                connected_buckets.push_back(&source_bucket);
            } else {
                Expansions::M2L(m_g[target_index],target_box,source_box,m_W[source_index],m_K)
            }
        }

        // now dive into the tree and do a proper FMM
        if (!search.is_leaf_node(bucket) { 
            vector_type& g = m_g[target_index];
            calculate_dive(connected_buckets,g,box,m_search.get_child1(bucket));
            calculate_dive(connected_buckets,g,box,m_search.get_child2(bucket));
        }
    }

};

   
template <typename Expansions, typename Function, typename NeighbourQuery>
class FastMultipoleMethod {
    typedef NeighbourQuery::traits_type traits_type;
    typedef NeighbourQuery::reference reference;
    typedef NeighbourQuery::pointer pointer;
    typedef Expansions::vector_type vector_type;
    typedef traits_type::template vector_type<vector_type> storage_type;
    typedef traits_type::template vector_type<traits_type::template vector_type<pointer>> connectivity_type;
    typedef iterator_range<NeighbourQuery::root_iterator> root_iterator_range;
    typedef iterator_range<NeighbourQuery::neighbouring_iterator> neighbouring_iterator_range;
    typedef m_query::particle_iterator particle_iterator;
    typedef particle_iterator::reference particle_reference;
    static const unsigned int dimension = traits_type::dimension;
    typedef detail::bbox<dimension> box_type;
    storage_type m_W;
    storage_type m_g;
    connectivity_type m_connectivity; 
    const NeighbourQuery *m_query;
    const Function m_K;

public:

    FastMultipoleMethod(const NeighbourQuery &query, const Function& K):
        m_query(&query),m_K(K)
    {}

    template <typename VectorType>
    calculate_expansions(const VectorType& source_vector) {
        root_iterator_range_type root_buckets = m_query->get_root_buckets();

        m_W.resize(m_query->number_of_buckets());
        m_g.resize(m_query->number_of_buckets());

        // upward sweep of tree
        // calculate P2M and M2M expansions for source buckets (recursive up
        // the tree so use function object)
        std::for_each(root_buckets.begin(),
                      root_buckets.end(),
                      calculate_P2M_and_M2M<Expansions,NeighbourQuery,
                                            VectorType,storage_type>(
                                                *m_query,
                                                source_vector,
                                                m_W));


        // downward sweep of tree
        // calculate L2L translations for source buckets (recursive 
        // down the tree so use function object)
        std::for_each(root_buckets.begin(),
                      root_buckets.end(),
                      calculate_M2L_and_L2L<Expansions,NeighbourQuery,
                                            VectorType,storage_type>(
                                                *m_query,
                                                m_W,
                                                m_g,
                                                m_connectivity));

    }


    // evaluate expansions for given point
    double evaluate_R_expansion(const double_d& p, const VectorType& source_vector) {
        const reference bucket = m_query->get_bucket(p);
        const size_t index = m_query->get_bucket_index(bucket); 
        box_type box(m_query->get_bucket_bounds_low(bucket),
                     m_query->get_bucket_bounds_high(bucket));

        double sum = Expansions::L2P(p,box,m_g[index].begin());
        for (pointer source_pointer: m_connectivity[index]) { 
            sum += detail::calculate_K_direct(
                    m_query->get_bucket_particles(*source_pointer)
                    ,m_K,source_vector);
        }
        return sum;
    }
};




