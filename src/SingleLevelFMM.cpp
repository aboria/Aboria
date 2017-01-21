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

//
// Acknowledgement: This source was modified from the Thrust example bucket_sort2d.cu
//


#ifndef SINGLELEVELFMM_H_
#define SINGLELEVELFMM_H_

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

struct mq_expansions {
    double R(x,y);
    double S(x,y);
    double S_R(x,y)
        etc.
}

template <typename Expansions, unsigned int R, 
          typename SourceNeighbourQuery>
struct calculate_S_expansion {
    SourceNeighbourQuery &search;
    typedef typename SourceNeighbourQuery::bucket_reference reference;
    
    calculate_S_expansion(SourceNeighbourQuery search):
        search(search) {}

    Vector<double,R> operator()(const reference bucket) {
        Vector<double,R> sum = 0;
        for (auto i&: search.get_bucket_particles(bucket)) {
            sum += Expansions::b<R>(r,xi,centre);
        }
        return sum;
    }
};

template <typename Expansions, unsigned int P1, unsigned int P2,
          typename SourceNeighbourQuery, typename TargetNeighbourQuery>
struct translate_S_R_accumulator {
    TargetNeighbourQuery &target;
    SourceNeighbourQuery &source;
    const static unsigned int Dsource = TargetNeighbourQuery::dimension;
    const static unsigned int Dtarget = SourceNeighbourQuery::dimension;
    detail::bbox<Dtarget> target_bbox;
    typedef typename SourceNeighbourQuery::bucket_reference reference;

    evaluate_B_l(SourceNeighbourQuery &source,
                 TargetNeighbourQuery &target,
                 size_t target_centre)
        :q(q),target(target),source(source),target_centre(target_centre) {}

    double operator()(const Vector& sum, const reference source_bucket) {
        detail::bbox<Dsource> source_bbox = source.get_bucket_bbox(source_bucket);
        if (well_separated(source_bbox,target_bbox)) {
            return sum + Expansions::S_R();
        } else {
            return sum;
        }
    }
};
 

template <typename Expansions, unsigned int P1, unsigned int P2,
          typename SourceNeighbourQuery, typename TargetNeighbourQuery>
struct translate_S_R {
    SourceNeighbourQuery &search;
    TargetNeighbourQuery &target;
    Traits::template vector_type<Vector<double,P1> &B;
    const static unsigned int Dtarget = SourceNeighbourQuery::dimension;
    typedef typename TargetNeighbourQuery::bucket_reference reference;

    evaluate_B_l(SourceNeighbourQuery &source,
                 TargetNeighbourQuery &target,
                 Traits::template vector_type<Vector<double,P1> &B
                 )
        :target(target),source(source),B(B) {}

    Vector<double,P2> operator()(const reference target_bucket) {
        detail::bbox<Dtarget> target_bbox = target.get_bucket_bbox(target_bucket);
        auto source_buckets = source.get_bucket_range();
        return Aboria::accumulate(source_buckets.begin(),
                           source_buckets.end(),
                           Vector<double,P2>(0),
                           translate_S_R_accumulator(source,target,target_bbox));
        
    }
};
   
template <typename Traits, typename Expansions, unsigned int P1, unsigned int P2,
         typename SourceNeighbourSearch, typename TargetNeighbourSearch>
class SingleLevelFMM {
    SourceNeighbourQuery m_source;
    TargetNeighbourQuery m_target;
    Traits::template vector_type<Vector<double,P2> A;
    Traits::template vector_type<Vector<double,P1> B;

public:

    SingleLevelFMM(SourceParticlesQuery &source, TargetParticlesQuery &target):
        m_source(source),
        m_target(target)
    {}

    // calculate S expansions for source buckets
    calculate_S_expansions() {
        auto source_buckets = m_source.get_bucket_range();
        B.resize(source_bucket_ids.size());
        return Aboria::transform(source_buckets.begin(),
                                 source_buckets.end(),
                                 B.begin(), 
                                 calculate_S_expansion<
                                    Expansions,
                                    P1,
                                    SourceNeighbourQuery>(m_source));
        
    }

    // calculate S|R translations for target buckets
    calculate_S_R_translations() {
        auto target_buckets = m_target.get_bucket_range();
        A.resize(target_bucket_ids.size());
        return Aboria::transform(target_buckets.begin(),
                                 target_buckets.end(),
                                 A.begin(), 
                                 translate_S_R<
                                    Expansions,
                                    P1,P2,
                                    SourceNeighbourQuery,
                                    TargetNeighbourQuery>(m_source,
                                                          m_target,
                                                          B));
    }

    // evaluate R expansions for given particle set
    Expansion::Scalar evaluate_R_expansion(position &p) {


    }
};




