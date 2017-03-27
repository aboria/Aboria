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

            lattice_iterator<D> mi(m_start,m_end,m_start);
            for (int i=0; i<ncheb; ++i,++mi) {
                for (int j=0; j<np; ++j) {
                    accum_begin[i] += row_Rn(*mi,j)*source_begin[j];
                }
            }
        }


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
             
            cheb_rn.calculate_Sn_with_bbox(source_positions.begin(),target_box,ncheb,N);

            lattice_iterator<D> mi(m_start,m_end,m_start);
            for (int i=0; i<ncheb; ++i,++mi) {
                for (int j=0; j<ncheb; ++j) {
                    accum_begin[i] += cheb_rn(*mi,j)*source_begin[j];
                }
            }
        }


        template <typename Function>
        void M2L(InOutIterator accum_begin, 
                 box_type& target_box, 
                 box_type& source_box, 
                 InputIterator source_begin,
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
                    accum_begin[i] += K(pj-pi,pi,pj)*source_begin[j];
                }
            }


        }

     

    };


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
            


template <typename Expansions, typename NeighbourQuery, 
          typename SourceVectorType, typename StorageVectorType>
struct calculate_L2M_and_M2M {
    NeighbourQuery &m_search;
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
        detail::bbox my_box(m_search.get_bucket_bounds_low(bucket),
                            m_search.get_bucket_bounds_high(bucket));
        size_t my_index = m_search.get_index(bucket);
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

template <typename Expansions, typename NeighbourQuery, 
          typename SourceVectorType, typename StorageVectorType>
struct calculate_L2L {
    const NeighbourQuery &m_search;
    const SourceVectorType &m_source_vector;
    StorageVectorType &m_g;
    typedef typename NeighbourQuery::reference reference;
    typedef Expansions::vector_type vector_type;
    typedef SourceNeighbourQuery::traits_type traits_type;
    static const unsigned int dimension = traits_type::dimension;
    
    calculate_S_expansion(NeighbourQuery& search, 
                          SourceVectorType& source_vector,
                          StorageVectorType& g):
        m_search(search),
        m_source_vector(source_vector),
        m_g(g)
    {}

    void calculate_L2L_dive(const vector_type& g_parent, detail::bbox<dimension>& bbox_parent, const reference bucket) {
        detail::bbox box(m_search.get_bucket_bounds_low(bucket),
                        m_search.get_bucket_bounds_high(bucket));
        size_t index = m_search.get_index(bucket);
        vector_type& g = m_g[index];
        Expansions::L2L(g.begin(),box,bbox_parent,g_parent.begin());
        if (search.is_leaf_node(bucket) { // leaf node
            calculate_L2L_dive(g,box,m_search.get_child1(bucket));
            calculate_L2L_dive(g,box,m_search.get_child2(bucket));
        }
    }

    vector_type& operator()(const reference bucket) {
        if (search.is_leaf_node(bucket) { // leaf node
            detail::bbox box(m_search.get_bucket_bounds_low(bucket),
                        m_search.get_bucket_bounds_high(bucket));
            size_t my_index = m_search.get_index(bucket);
            const vector_type& g = m_g[my_index];
            calculate_L2L_dive(g,box,m_search.get_child1(bucket));
            calculate_L2L_dive(g,box,m_search.get_child2(bucket));
        }
    }

};

   
template <typename Expansions, typename NeighbourQuery>
class FastMultipoleMethod {
    typedef NeighbourQuery::traits_type traits_type;
    typedef Expansions::vector_type vector_type;
    typedef traits_type::template vector_type<vector_type> storage_type;
    typedef iterator_range<NeighbourQuery::root_iterator> root_iterator_range_type;
    typedef iterator_range<NeighbourQuery::all_iterator> all_iterator_range_type;
    storage_type m_W;
    storage_type m_g;
    NeighbourQuery *m_source;

public:

    FastMultipoleMethod(NeighbourQuery &source):
        m_source(&source)
    {}

    template <typename Function, typename VectorType>
    calculate_expansions(const Function& K, const VectorType& source_vector) {
        root_iterator_range_type root_buckets = m_source->get_root_buckets();
        all_iterator_range_type all_buckets = m_source->get_all_buckets();

        m_W.resize(m_source->number_of_buckets());
        m_g.resize(m_source->number_of_buckets());

        // calculate L2M and M2M expansions for source buckets (recursive up
        // the tree so use function object)
        std::for_each(root_buckets.begin(),
                      root_buckets.end(),
                      calculate_P2M_and_M2M<Expansions,NeighbourQuery,
                                            VectorType,storage_type>(
                                                *m_source,
                                                source_vector,
                                                m_W));


        // calculate M2L translations for source buckets (this should be the most work)
        std::for_each(all_buckets.begin(),all_buckets.end(),
                      [&](const reference target_bucket) {
            detail::bbox<dimension> target_bbox(m_source->get_bucket_bounds_low(target_bucket),
                                          m_target->get_bucket_bounds_high(target_bucket));
            size_t target_index = m_source->get_index(target_bucket); 
            well_separated_iterator_range_type source_buckets 
                    = m_source->get_well_separated_buckets(target_bucket);
            std::for_each(source_buckets.begin(), source_buckets.end(), 
                    [&](const reference source_bucket) {
                detail::bbox<dimension> source_bbox(m_source->get_bucket_bounds_low(source_bucket),
                                                  m_source->get_bucket_bounds_high(source_bucket));
                size_t source_index = m_source->get_index(source_bucket); 
                Expansions::M2L(m_sum.begin(),m_target_bbox,source_bbox,m_W[source_index].begin(),K)
            });
        });


        // calculate L2L translations for source buckets (recursive 
        // down the tree so use function object)
        std::for_each(root_buckets.begin(),
                      root_buckets.end(),
                      calculate_L2L<Expansions,NeighbourQuery,
                                            VectorType,storage_type>(
                                                *m_source,
                                                source_vector,
                                                m_W));

    }


    // evaluate expansions for given point
    Expansion::Scalar evaluate_R_expansion(const double_d& p) {
        const reference bucket = m_source->get_bucket(p);
        const size_t index = m_source->get_index(bucket); 
        detail::bbox<dimension> box(m_source->get_bucket_bounds_low(bucket),
                                    m_source->get_bucket_bounds_high(bucket));
        target_neighbouring_iterator_range neighbouring_buckets = m_source->get_neighbouring_buckets(bucket);

        Expansions::L2P(p,m_target_bbox,source_bbox,m_W[source_index].begin(),K)
        

    }
};




