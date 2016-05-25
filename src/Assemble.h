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


#ifndef ASSEMBLE_H_
#define ASSEMBLE_H_

#include <type_traits>
#include "Symbolic.h"


#ifdef HAVE_EIGEN
#include <Eigen/Sparse>
#include <Eigen/Dense>
#endif


namespace Aboria {

#ifdef HAVE_EIGEN

template <typename Scalar, int Rows, int Cols, typename Expr, 
         typename = typename std::enable_if<proto::matches<R, bivariate_expr> >::type>
void assemble(Eigen::Matrix<Scalar, Rows, Cols>& matrix, const Expr& expr) {
    const R &reactants = proto::value(*this).second;
    P products = proto::value(*this).first;
        
    typedef typename std::result_of<bivariate_expr(R)>::type::first_type particles_a_type;
    typedef typename std::result_of<bivariate_expr(R)>::type::second_type particles_b_type;

    const particles_a_type& a = bivariate_expr()(expr).first;
    const particles_b_type& b = bivariate_expr()(expr).second;

    const size_t na = a.size();
    const size_t nb = b.size();

    matrix.resize(a,b);

    CHECK((matrix.rows() == na) && (matrix.cols() == nb), "matrix size is not compatible with expression. expr = ("<<na<<","<<nb<<") and matrix = ("<<matrix.rows()<<","<<matrix.cols()<<").")

    for (size_t i=0; i<na; ++i) {
        for (size_t j=0; j<na; ++j) {
            matrix(i,j) = expr.eval(a[i],b[i]);
        }
    }
}

template <typename Scalar, typename Expr, typename ifExpr,
         typename = typename std::enable_if<proto::matches<R, bivariate_expr> >::type>
void assemble(Eigen::SparseMatrix<Scalar>& matrix, const Expr& expr, const ifExpr& if_expr) {
    const R &reactants = proto::value(*this).second;
    P products = proto::value(*this).first;
        
    typedef typename std::result_of<bivariate_expr(R)>::type::first_type particles_a_type;
    typedef typename std::result_of<bivariate_expr(R)>::type::second_type particles_b_type;

    const particles_a_type& a = bivariate_expr()(expr).first;
    const particles_b_type& b = bivariate_expr()(expr).second;

    const size_t na = a.size();
    const size_t nb = b.size();

    matrix.resize(nb,na);

    CHECK((matrix.cols() == na) && (matrix.rows() == nb), "matrix size is not compatible with expression. expr = ("<<na<<","<<nb<<") and matrix = ("<<matrix.cols()<<","<<matrix.rows()<<").")

    std::map<size_t,size_t> map_id_to_index;
    for (size_t j=0; j<nb; ++j) {
        map_id_to_index[get<id>(b[j])] = j;
    }

    //std::vector<size_t> cols_sizes(na);
    typedef Eigen::Triplet<Scalar> triplet_type;
    std::vector<triplet_type> tripletList;
    tripletList.reserve(na*5); 

    for (size_t i=0; i<na; ++i) {
        //cols_sizes[i] = 0;
        const particles_a_type::value_type& ai = a[i];
        for (auto bj: b.get_neighbours(get<position>(ai))) {
            if (if_expr.eval(ai,bi)) {
                const size_t j = map_id_to_index[get<id>(bj)];
                tripletList.push_back(triplet_type(i,j,expr.eval(ai,bi)));
                //++cols_sizes[i];
            }
        }
    }
    matrix.setFromTriplets(tripletList.begin(),tripletList.end());
}

#endif

}

#endif //ASSEMBLE_H_
