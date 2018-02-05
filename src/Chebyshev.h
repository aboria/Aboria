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


#ifndef CHEBYSHEV_H_
#define CHEBYSHEV_H_

#include "detail/Chebyshev.h"
#include "NeighbourSearchBase.h"

#ifdef HAVE_EIGEN
#include <Eigen/Core>
#endif

namespace Aboria {

template <unsigned int D,
         typename InputIterator,
         typename OutputIterator,
         typename PositionIterator,
         typename Kernel>
void chebyshev_interpolation(
        InputIterator input_iterator_begin, 
        InputIterator input_iterator_end, 
        OutputIterator output_iterator_begin,  
        OutputIterator output_iterator_end,  
        PositionIterator source_positions_begin,
        PositionIterator target_positions_begin,
        Kernel kernel,
        unsigned int n) {
#ifdef HAVE_EIGEN
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_type;
    typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type; 
    typedef Eigen::Map<vector_type> map_type;
    typedef Vector<double,D> double_d;
    typedef Vector<int,D> int_d;

    const size_t sourceN = std::distance(input_iterator_begin,input_iterator_end);
    const size_t targetN = std::distance(output_iterator_begin,output_iterator_end);
    const size_t ncheb = std::pow(n,D);
    const int_d start = int_d(0);
    const int_d end = int_d(n);

    map_type source_values(&(*input_iterator_begin),sourceN);
    map_type target_values(&(*output_iterator_begin),targetN);

    // Calculate Sn's at source and target positions
    detail::Chebyshev_Rn<D> source_Rn,target_Rn;
    source_Rn.calculate_Sn(source_positions_begin,sourceN,n);
    target_Rn.calculate_Sn(target_positions_begin,targetN,n);
    
    // fill source_Rn matrix
    matrix_type source_Rn_matrix(ncheb,sourceN);
    for (int i=0; i<sourceN; ++i) {
        lattice_iterator<D> mi(start,end);
        for (int j=0; j<ncheb; ++j,++mi) {
            source_Rn_matrix(j,i) = source_Rn(*mi,i);
        }
    }

    //First compute the weights at the Chebyshev nodes ym by anterpolation 
    vector_type W = source_Rn_matrix*source_values;

    // fill kernel matrix
    matrix_type kernel_matrix(ncheb,ncheb);
    lattice_iterator<D> mi(start,end);
    for (int i=0; i<ncheb; ++i,++mi) {
        const double_d pi = source_Rn.get_position(*mi);
        lattice_iterator<D> mj(start,end);
        for (int j=0; j<ncheb; ++j,++mj) {
            const double_d pj = target_Rn.get_position(*mj);
            kernel_matrix(i,j) = kernel(pi-pj,pj,pi);
        }
    }

    //Next compute f ðxÞ at the Chebyshev nodes xl:
    vector_type fcheb = kernel_matrix*W;

    // fill target_Rn matrix
    matrix_type target_Rn_matrix(targetN,ncheb);
    for (int i=0; i<targetN; ++i) {
        lattice_iterator<D> mj(start,end);
        for (int j=0; j<ncheb; ++j,++mj) {
            target_Rn_matrix(i,j) = target_Rn(*mj,i);
        }
    }

    //Last compute f ðxÞ at the observation points xi by interpolation:
    target_values = target_Rn_matrix*fcheb;
#else
    ERROR("chebyshev_interpolation requires the Eigen library");
#endif
}


}

#endif
