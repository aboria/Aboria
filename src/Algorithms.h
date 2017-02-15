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


#ifndef ALGORITHMS_H_
#define ALGORITHMS_H_

namespace Aboria {

/*
void dense_linear_operator<variable_lhs,variable_rhs>(particles_i, particle_j, kernel_functor) 
void sparse_radial_linear_operator<variable_lhs,variable_rhs>(particles_i, particle_j, kernel_functor, radius) 
void linear_operator<variable_lhs,variable_rhs>(particles_i, particle_j, kernel_functor) 
void linear_operator_solve<variable_lhs,variable_rhs>(particles_i, particle_j, kernel_functor) 
void non_linear_operator<variable_lhs>(particles_i, particle_j, kernel_functor) 
void non_linear_operator_solve<variable_lhs>(particles_i, particle_j, kernel_functor, gradient_kernel_functor) 
*/

void chebyshev_interpolation(
        input_iterator_begin, 
        output_iterator_begin,  
        source_positions,
        target_positions,
        kernel,
        unsigned int n) {
    typedef Eigen::Matrix<double,Eigen::dynamic,Eigen::dynamic> matrix_type;
    typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type; 
    typedef Eigen::Map<vector_type> map_type;

    // fill source_Rn matrix
    detail::Chebyshev_Rn<D> Rn_eval;
    Rn.calculate_Sn(source_positions,n);
    const int_d start = int_d(0);
    const int_d end = int_d(n-1);
    lattice_iterator<D> m(start,end,start);
    const size_t ncheb = std::pow(n,D);
    matrix_type source_Rn(ncheb,positions.size());
    for (int i=0; i<positions.size(); ++i) {
        for (int j=0; j<ncheb; ++j,++m) {
            source_Rn(j,i) = Rn_eval(m,i);
        }
    }

    // fill kernel matrix
    matrix_type kernel_matrix(ncheb,ncheb);
    for (int i=0; i<ncheb; ++i,++mi) {
        double_d pi = ;
        for (int j=0; j<ncheb; ++j,++m) {
            // get position
            double_d pj = ;
            kernel_matrix(i,j) = kernel(pi,pj);
        }
    }

    // fill target_Rn matrix
    Rn.calculate_Sn(target_positions,n);
    const int_d start = int_d(0);
    const int_d end = int_d(n-1);
    lattice_iterator<D> m(start,end,start);
    const size_t ncheb = std::pow(n,D);
    matrix_type target_Rn(positions.size(),ncheb);
    for (int i=0; i<positions.size(); ++i) {
        for (int j=0; j<ncheb; ++j,++m) {
            target_Rn(i,j) = Rn_eval(m,i);
        }
    }

    map_type source_values(&(*input_iterator_begin),source_size);
    map_type target_values(&(*output_iterator_begin),target_size);

    //First compute the weights at the Chebyshev nodes ym by anterpolation 
    //Next compute f ðxÞ at the Chebyshev nodes xl:
    //Last compute f ðxÞ at the observation points xi by interpolation:
    target_values = target_Rn*kernel_matrix*source_Rn*source_values;
}

}

#endif
