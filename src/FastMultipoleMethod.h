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
#include "Log.h"

namespace Aboria {

namespace detail {

    template <typename InOutIterator, unsigned int D, unsigned int N, typename PositionIterator, typename InputIterator> 
    struct BlackBoxExpansions {
        typedef detail::bbox<D> box_type;
        const size_t ncheb = std::pow(N,D); 
        typedef std::array<double,ncheb> vector_type;


        inline 
        void L2M(InOutIterator accum_begin, 
                 box_type& box, 
                 PositionIterator position_begin, 
                 PositionIterator position_end, 
                 InputIterator source_begin) {

            Chebyshev_Rn cheb_rn;
            const size_t np = position_end - position_begin;
            cheb_rn.calculate_Sn(position_begin,np,N);

            int_d start(0);
            int_d end(N-1);

            for (int i=0; i<np; ++i) {
                lattice_iterator<D> mj(m_start,m_end,m_start);
                for (int j=0; j<ncheb; ++j,++mj) {
                    accum_begin[j] += row_Rn(*mj,i)*source_begin[i];
                }
            }
        }


        inline 
        void M2M(InOutIterator accum_begin, 
                 box_type& target_box, 
                 box_type& source_box, 
                 InputIterator source_begin) {

            int_d start(0);
            int_d end(N-1);
            std::array<double_d,ncheb> source_positions;

            lattice_iterator<dimension> mi(m_start,m_end,m_start);
            for (int i=0; i<ncheb; ++i,++mi) {
                const double_d pi_unit_box = detail::chebyshev_node_nd(*mi);
                source_positions[i] = 0.5*(pi_unit_box+1)*(source_box.bmax-source_box.bmin) 
                                        + source_box.bmin;
            }
             
            cheb_rn.calculate_Sn(source_positions.begin(),ncheb,N);

            for (int i=0; i<ncheb; ++i) {
                const double_d pi = 0.5*(cheb_node_unit_box+1)*(box.bmax-box.bmin) + box.bmin;
                lattice_iterator<D> mj(m_start,m_end,m_start);
                for (int j=0; j<ncheb; ++j,++mj) {
                    accum_begin[j] += row_Rn(*mj,i)*source_begin[i];
                }
            }
        }

        template <typename VectorType, unsigned int D, unsigned int N, typename SourceVectorType>
        void calculate_L2M(vector_type& sum, detail::bbox& my_box, iterator_range<ranges_iterator<traits_type>>& range, sou) {
            const size_t N = std::distance(range.begin(),range.end());
            const double_d* pbegin = &get<position>(*range.begin());
            const size_t pindex = pbegin - &get<position>(source_particles)[0];
            return Expansions::L2M(sum, my_box,pbegin,pbegin + N,source_vector.begin() + index);  
        }

        // assume serial processing of particles, this could be more efficient for ranges iterators
        template <typename Iterator, typename = typename
            std::enable_if<!std::is_same<Iterator,ranges_iterator<traits_type>>>
            void calculate_L2M(vector_type& sum, detail::bbox& my_box, iterator_range<Iterator>& range) {

                const size_t N = std::distance(range.begin(),range.end());
                std::vector<double_d> positions(N);
                std::transform(range.begin(),range.end(),positions.begin(),
                        [](const_particle_reference p) { return get<position>(p); });
                std::vector<double> source_values(N);
                std::transform(range.begin(),range.end(),source_values.begin(),
                        [](const_particle_reference p) { 
                        const size_t index = (&get<position>(p))-(&get<position>(source_particles)[0]); 
                        return source_vector[index];
                        });

                Expansions::L2M(sum, my_box,positions.begin(), positions.end(), source_values.begin());

            }
    }
            


template <typename Expansions, typename SourceNeighbourQuery, 
          typename SourceVectorType, typename StorageVectorType>
struct calculate_S_expansion {
    SourceNeighbourQuery &m_search;
    SourceVectorType &m_source_vector;
    StorageVectorType &m_W;
    typedef typename SourceNeighbourQuery::reference reference;
    typedef typename SourceNeighbourQuery::particle_iterator particle_iterator;
    typedef Expansions::vector_type vector_type;
    typedef SourceNeighbourQuery::traits_type traits_type;
    static const unsigned int dimension = traits_type::dimension;
    
    calculate_S_expansion(SourceNeighbourQuery& search, 
                          SourceVectorType& source_vector,
                          StorageVectorType& W):
        m_search(search),
        m_source_vector(source_vector),
        m_W(W)
    {}

    vector_type& operator()(const reference bucket) {
        detail::bbox my_box(search.get_bucket_bounds_low(bucket),
                            search.get_bucket_bounds_high(bucket));
        size_t my_index = search.get_index(bucket);
        vector_type& W = m_W[my_index];

        if (search.is_leaf_node(bucket) { // leaf node
            calculate_L2M(W, my_box, search.get_bucket_particles(bucket));
        } else { // assume binary tree for now
            const reference child1 = search.get_child1(bucket);
            const reference child2 = search.get_child2(bucket);
            detail::bbox child1_box(search.get_bucket_bounds_low(child1),
                                search.get_bucket_bounds_high(child1));
            detail::bbox child2_box(search.get_bucket_bounds_low(child2),
                                search.get_bucket_bounds_high(child2));

            vector_type& child1_W = this->operator()(child1);
            vector_type& child2_W = this->operator()(child2);
            Expansions::M2M(W.begin(),my_box,child1_box,child1_W.begin());
            Expansions::M2M(W.begin(),my_box,child1_box,child2_W.begin());
        }
        return W;
    }
};

template <typename Expansions, typename StorageVectorType
          typename SourceNeighbourQuery, typename TargetNeighbourQuery>
struct translate_S_R_accumulator {
    TargetNeighbourQuery &m_target;
    SourceNeighbourQuery &m_source;
    StorageVectorType &m_W;
    const static unsigned int dimension = TargetNeighbourQuery::dimension;
    detail::bbox<Dtarget> m_target_bbox;
    typedef typename SourceNeighbourQuery::bucket_reference reference;
    typedef Expansions::vector_type vector_type;
    vector_type &m_sum;

    evaluate_B_l(SourceNeighbourQuery& source,
                 TargetNeighbourQuery& target,
                 detail::bbox<Dtarget>& target_bbox,
                 StorageVectorType& W,
                 vector_type& sum)
        :q(q),target(target),source(source),m_W(W),m_target_bbox(target_bbox),m_sum {}

    vector_type operator()(const reference source_bucket) {
        detail::bbox<dimension> source_bbox(source.get_bucket_bounds_low(source_bucket),
                                          source.get_bucket_bounds_high(source_bucket));
        size_t source_index = m_source.get_index(source_bucket); 
        Expansions::M2L(m_sum.begin(),m_target_bbox,source_bbox,m_W[source_index].begin());
    }
};
 

template <typename Expansions, typename StorageVectorType 
          typename SourceNeighbourQuery, typename TargetNeighbourQuery>
struct translate_S_R {
    TargetNeighbourQuery &m_target;
    SourceNeighbourQuery &m_source;
    StorageVectorType& m_W;
    StorageVectorType& m_g;
    const static unsigned int dimension = SourceNeighbourQuery::dimension;
    typedef typename TargetNeighbourQuery::bucket_reference reference;
    typedef Expansions::vector_type vector_type;
    typedef iterator_range<SourceNeighbourQuery::well_separated_iterator> well_separated_iterator_range_type;

    evaluate_B_l(SourceNeighbourQuery &source,
                 TargetNeighbourQuery &target,
                 StorageVectorType& W,
                 StorageVectorType& g)
        :m_target(target),m_source(source),m_W(W),m_g(g) {}

    void operator()(const reference target_bucket) {
        detail::bbox<dimension> target_bbox(m_target.get_bucket_bounds_low(target_bucket),
                                          m_target.get_bucket_bounds_high(target_bucket));
        size_t target_index = m_target.get_index(target_bucket); 
        well_separated_iterator_range_type source_buckets 
                    = source.get_well_separated_buckets(target_bbox);
        std::for_each(source_buckets.begin(),
                      source_buckets.end(),
                      translate_S_R_accumulator(m_source,m_target,target_bbox,m_W,m_g[my_index]));
        
    }
};
   
template <typename Traits, typename Expansions,
         typename SourceNeighbourQuery, typename TargetNeighbourQuery>
class FastMultipoleMethod {
    typedef SourceNeighbourQuery::traits_type traits_type;
    typedef ExpansionsSource::vector_type vector_type;
    typedef traits_type::template vector_type<vector_type> storage_type;
    typedef iterator_range<SourceNeighbourQuery::root_iterator> source_root_iterator_range;
    storage_type m_W;
    storage_type m_g;
    SourceNeighbourQuery *m_source;

public:

    SingleLevelFMM(SourceParticlesQuery &source):
        m_source(&source)
    {}

    template <typename SourceVectorType>
    calculate_expansions(SourceVectorType& source_vector) {
        source_root_iterator_range source_root_buckets = m_source->get_root_buckets();
        source_all_iterator_range source_all_buckets = m_source.get_all_buckets();

        m_W.resize(m_source->number_of_buckets());
        m_g.resize(m_source->number_of_buckets());

        // calculate L2M and M2M expansions for source buckets
        std::for_each(source_root_buckets.begin(),
                      source_root_buckets.end(),
                      calculate_S_expansion<ExpansionsSource,SourceNeighbourQuery,
                                            SourceVectorType,source_storage_type>(
                                                *m_source,
                                                source_vector,
                                                m_W));


        // calculate M2L translations for source buckets
        std::for_each(source_all_buckets.begin(),
                      source_all_buckets.end(),
                      translate_S_R<Expansions,StorageVectorType,
                                    SourceNeighbourQuery,TargetNeighbourQuery>(
                              m_source
                              m_W,m_g));

        // calculate L2L translations for source buckets
        std::for_each(source_root_buckets.begin(),
                      source_root_buckets.end(),
                      calculate_R_expansion<ExpansionsSource,SourceNeighbourQuery,
                                            SourceVectorType,source_storage_type>(
                                                *m_source,
                                                source_vector,
                                                m_W));

    }


    // evaluate expansions for given point
    Expansion::Scalar evaluate_R_expansion(const double_d& p) {
        target_neighbouring_iterator_range neighbouring_buckets = m_source->get_neighbouring_buckets(p);
        

    }
};




