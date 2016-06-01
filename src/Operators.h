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


#ifndef OPERATORS_H_
#define OPERATORS_H_

#include <type_traits>
#include "Symbolic.h"


#ifdef HAVE_EIGEN
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#endif


namespace Aboria {
    /*
    template<std::size_t... I>
    auto make_unknown(const size_t i, const size_t nb, data_type& rhs, index_sequence<I...>) 
        -> decltype(std::make_tuple(std::get<I>(rhs[i+I*nb]...))) {
        return std::make_tuple(std::get<I>(rhs[i+I*nb]...));
    }
    */

#ifdef HAVE_EIGEN
    template <typename Expr, typename IfExpr> 
    class MatrixReplacement;

    template<typename Rhs, typename Expr, typename IfExpr> class 
    MatrixReplacement_ProductReturnType;

}

namespace Eigen {
    namespace internal {
        template<typename Expr, typename IfExpr>
        struct traits<Aboria::MatrixReplacement<Expr,IfExpr>> :  Eigen::internal::traits<Eigen::SparseMatrix<double> >
        {};
        template <typename Rhs, typename Expr, typename IfExpr>
        struct traits<Aboria::MatrixReplacement_ProductReturnType<Rhs,Expr,IfExpr> > {
            // The equivalent plain objet type of the product. This type is used if the product needs to be evaluated into a temporary.
            typedef Eigen::Matrix<typename Rhs::Scalar, Eigen::Dynamic, Rhs::ColsAtCompileTime> ReturnType;
        };
    }
}

namespace Aboria {
    // Inheriting EigenBase should not be needed in the future.
    template <typename Expr, typename IfExpr> 
    class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement<Expr,IfExpr>> {
        public:
        typedef typename std::result_of<detail::bivariate_expr(Expr)>::type::first_type label_a_type_ref;
        typedef typename std::result_of<detail::bivariate_expr(Expr)>::type::second_type label_b_type_ref;
        typedef typename std::remove_reference<label_a_type_ref>::type label_a_type;
        typedef typename std::remove_reference<label_b_type_ref>::type label_b_type;
        typedef typename label_a_type::particles_type particles_a_type;
        typedef typename label_b_type::particles_type particles_b_type;
        typedef typename particles_b_type::double_d double_d;
        typedef typename particles_b_type::position position;
        const static unsigned int dimension = particles_b_type::dimension;                          \


            // Expose some compile-time information to Eigen:
            typedef double Scalar;
            typedef double RealScalar;
            typedef size_t Index;
            enum {
                ColsAtCompileTime = Eigen::Dynamic,
                RowsAtCompileTime = Eigen::Dynamic,
                MaxColsAtCompileTime = Eigen::Dynamic,
                MaxRowsAtCompileTime = Eigen::Dynamic
            };

            MatrixReplacement(const Expr& expr, const IfExpr& if_expr):m_expr(expr),m_if_expr(if_expr) {};

            Index rows() const { 
                const particles_a_type& a = detail::bivariate_expr()(m_expr).first.get_particles();
                return a.size(); 
            }

            Index cols() const { 
                const particles_b_type& b = detail::bivariate_expr()(m_expr).second.get_particles();
                return b.size();
            }
            void resize(Index a_rows, Index a_cols) {
                // This method should not be needed in the future.
                assert(a_rows==0 && a_cols==0 || a_rows==rows() && a_cols==cols());
            }
            // In the future, the return type should be Eigen::Product<MatrixReplacement,Rhs>
            template<typename Rhs>
            MatrixReplacement_ProductReturnType<Rhs,Expr,IfExpr> operator*(const Eigen::MatrixBase<Rhs>& x) const {
                return MatrixReplacement_ProductReturnType<Rhs,Expr,IfExpr>(*this, x.derived());
            }

            const Expr& m_expr;
            const IfExpr& m_if_expr;

    };
    // The proxy class representing the product of a MatrixReplacement with a MatrixBase<>
    template<typename Rhs, typename Expr, typename IfExpr>
    class MatrixReplacement_ProductReturnType : public Eigen::ReturnByValue<MatrixReplacement_ProductReturnType<Rhs,Expr,IfExpr> > {
        typedef typename MatrixReplacement<Expr,IfExpr>::particles_a_type particles_a_type;
        typedef typename MatrixReplacement<Expr,IfExpr>::particles_b_type particles_b_type;
        public:
            typedef typename MatrixReplacement<Expr,IfExpr>::Index Index;

            // The ctor store references to the matrix and right-hand-side object (usually a vector).
            MatrixReplacement_ProductReturnType(const MatrixReplacement<Expr,IfExpr>& matrix, const Rhs& rhs)
                : m_matrix(matrix), m_rhs(rhs)
            {}

            Index rows() const { return m_matrix.rows(); }
            Index cols() const { return m_rhs.cols(); }
            // This function is automatically called by Eigen. It must evaluate the product of matrix * rhs into y.
            template<typename Dest>
            void evalTo(Dest& y) const {
                const particles_a_type& a = detail::bivariate_expr()(m_matrix.m_expr).first.get_particles();
                const particles_b_type& b = detail::bivariate_expr()(m_matrix.m_expr).second.get_particles();

                const Index na = m_matrix.rows();
                const Index nb = m_matrix.cols();

                ASSERT(nb == m_rhs.rows(),"In Eigen-Aboria Matrix-Vector product: The number of columns in the  Matrix replacement and the length of the rhs vector do not match");

                y.setZero(na);
                for (size_t i=0; i<na; ++i) {
                    typename particles_a_type::const_reference ai = a[i];
                    for (auto pairj: b.get_neighbours(get<typename particles_a_type::position>(ai))) {
                        const typename particles_a_type::double_d & dx = std::get<1>(pairj);
                        typename particles_b_type::const_reference bj = std::get<0>(pairj);
                        double sum = 0;
                        if (eval(m_matrix.m_if_expr,dx,ai,bj)) {
                            sum += eval(m_matrix.expr,dx,ai,bj);
                        }
                        y(i) = sum;
                    }
                }
            }

        protected:
            const MatrixReplacement<Expr,IfExpr>& m_matrix;
            typename Rhs::Nested m_rhs;
    };






    
    template <typename Expr, typename IfExpr=detail::SymbolicExpr<typename proto::terminal<bool>::type>,
    typename = typename std::enable_if<proto::matches<Expr, detail::bivariate_expr>::value >::type>
    MatrixReplacement<Expr,IfExpr> create_eigen_operator(const Expr& expr, const IfExpr& if_expr = IfExpr(proto::terminal<bool>::type::make(true))) {
        return MatrixReplacement<Expr,IfExpr>(expr,if_expr);
    }

#endif

}

#endif //OPERATORS_H_
