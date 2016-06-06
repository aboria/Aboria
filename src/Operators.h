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
    template <unsigned int NI, unsigned int NJ, typename Blocks> 
    class MatrixReplacement;

    template<typename Rhs, unsigned int NI, unsigned int NJ, typename Blocks> class 
    MatrixReplacement_ProductReturnType;

}

namespace Eigen {
    namespace internal {
        template<unsigned int NI, unsigned int NJ, typename Blocks>
        struct traits<Aboria::MatrixReplacement<NI, NJ, Blocks>> :  Eigen::internal::traits<Eigen::SparseMatrix<double> >
        {};
        template <typename Rhs, unsigned int NI, unsigned int NJ, typename Blocks>
        struct traits<Aboria::MatrixReplacement_ProductReturnType<Rhs,NI,NJ,Blocks> > {
            // The equivalent plain objet type of the product. This type is used if the product needs to be evaluated into a temporary.
            typedef Eigen::Matrix<typename Rhs::Scalar, Eigen::Dynamic, Rhs::ColsAtCompileTime> ReturnType;
        };
    }
}

namespace Aboria {

    struct one {
        struct dummy {
            size_t size() { return 1; }
        };
        const dummy& get_particles() { return d; }
        dummy d;
    };

    int sum() {
        return 0;
    }

    template<typename T1, typename... T>
    int sum(T1 s, T... ts) {
        std::cout << "s = "<<s<<std::endl;
        return s + sum(ts...);
    }

    // Inheriting EigenBase should not be needed in the future.
    template <unsigned int NI, unsigned int NJ, typename Blocks> 
    class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement<NI,NJ,Blocks>> {

        public:
        
    
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

            MatrixReplacement(const Blocks& blocks):m_blocks(blocks) {};


            template<std::size_t... I>
            Index rows_impl(index_sequence<I...>) const {
                return sum(std::get<0>(std::get<I*NJ>(m_blocks)).size()...);
            }

            template<std::size_t... J>
            Index cols_impl(index_sequence<J...>) const {
                return sum(std::get<1>(std::get<J>(m_blocks)).size()...);
            }

            Index rows() const { 
                std::cout << "rows = " << rows_impl(make_index_sequence<NI>()) << std::endl;
                return rows_impl(make_index_sequence<NI>());
            }

            Index cols() const { 
                std::cout << "cols = " << cols_impl(make_index_sequence<NJ>())<< std::endl;
                return cols_impl(make_index_sequence<NJ>());
            }

            template <int I>
            Index start_col() const { 
                return cols_impl(make_index_sequence<I>());
            }

            template <int I>
            Index size_col() const { 
                return std::get<1>(std::get<I>(m_blocks)).size();
            }

            template <int I>
            Index start_row() const { 
                return rows_impl(make_index_sequence<I>());
            }

            template <int I>
            Index size_row() const { 
                return std::get<0>(std::get<I*NJ>(m_blocks)).size();
            }

            void resize(Index a_rows, Index a_cols) {
                // This method should not be needed in the future.
                assert(a_rows==0 && a_cols==0 || a_rows==rows() && a_cols==cols());
            }
            // In the future, the return type should be Eigen::Product<MatrixReplacement,Rhs>
            template<typename Rhs>
            MatrixReplacement_ProductReturnType<Rhs,NI,NJ,Blocks> operator*(const Eigen::MatrixBase<Rhs>& x) const {
                return MatrixReplacement_ProductReturnType<Rhs,NI,NJ,Blocks>(*this, x.derived());
            }

            const Blocks m_blocks;

    };
    // The proxy class representing the product of a MatrixReplacement with a MatrixBase<>
    template<typename Rhs, unsigned int NI, unsigned int NJ, typename Blocks>
    class MatrixReplacement_ProductReturnType : public Eigen::ReturnByValue<MatrixReplacement_ProductReturnType<Rhs,NI,NJ,Blocks> > {
        public:
            typedef typename MatrixReplacement<NI,NJ,Blocks>::Index Index;

            // The ctor store references to the matrix and right-hand-side object (usually a vector).
            MatrixReplacement_ProductReturnType(const MatrixReplacement<NI,NJ,Blocks>& matrix, const Rhs& rhs)
                : m_matrix(matrix), m_rhs(rhs)
            {}

            Index rows() const { return m_matrix.rows(); }
            Index cols() const { return m_rhs.cols(); }

            template <typename Dest, typename Source, typename particles_a_type, typename particles_b_type,
                                     typename expr_type, typename if_expr_type>
            void evalTo_block(Eigen::VectorBlock<Dest> y, const Eigen::VectorBlock<Source>& rhs, const std::tuple<const particles_a_type&,const particles_b_type&,expr_type,if_expr_type>& block) const {
                typedef typename particles_b_type::double_d double_d;
                typedef typename particles_b_type::position position;
                const static unsigned int dimension = particles_b_type::dimension;

                typedef typename std::result_of<detail::bivariate_expr(expr_type)>::type::first_type label_a_type_ref;
                typedef typename std::result_of<detail::bivariate_expr(expr_type)>::type::second_type label_b_type_ref;
                typedef typename std::remove_reference<label_a_type_ref>::type label_a_type_check;
                typedef typename std::remove_reference<label_b_type_ref>::type label_b_type_check;
                typedef typename label_a_type_check::particles_type particles_a_type_check;
                typedef typename label_b_type_check::particles_type particles_b_type_check;

                static_assert(std::is_same<particles_a_type_check,particles_a_type>::value, "particles type a in expression not the same as particles type given to block");
                static_assert(std::is_same<particles_b_type_check,particles_b_type>::value, "particles type b in expression not the same as particles type given to block");

                const particles_a_type& a = std::get<0>(block);
                const particles_b_type& b = std::get<1>(block);
                const expr_type expr = std::get<2>(block);
                const if_expr_type if_expr = std::get<3>(block);

                const Index na = a.size();
                const Index nb = b.size();

                for (size_t i=0; i<na; ++i) {
                    typename particles_a_type::const_reference ai = a[i];
                    double sum = 0;
                    if (boost::is_same<if_expr_type,detail::SymbolicExpr<typename proto::terminal<bool>::type>>::value) {
                        for (size_t j=0; j<nb; ++j) {
                            typename particles_b_type::const_reference bj = b[j];
                            const double_d dx = a.correct_dx_for_periodicity(get<position>(bj)-get<position>(ai));
                            sum += eval(expr,dx,ai,bj)*rhs(j);
                        }
                    } else {
                        for (auto pairj: b.get_neighbours(get<position>(ai))) {
                            const double_d & dx = std::get<1>(pairj);
                            typename particles_b_type::const_reference bj = std::get<0>(pairj);
                            const Index j = &get<position>(bj) - get<position>(b).data();
                            if (eval(if_expr,dx,ai,bj)) {
                                sum += eval(expr,dx,ai,bj)*rhs(j);
                            }
                        }
                    }
                    y(i) += sum;
                }
            }


            template <typename Dest, typename Source, typename particles_b_type,
                                     typename expr_type, typename if_expr_type>
            void evalTo_block(Eigen::VectorBlock<Dest> y, const Eigen::VectorBlock<Source>& rhs, const std::tuple<const one&,const particles_b_type&,expr_type,if_expr_type>& block) const {

                typedef typename particles_b_type::double_d double_d;
                typedef typename particles_b_type::position position;
                const static unsigned int dimension = particles_b_type::dimension;

                typedef typename std::result_of<detail::univariate_expr(expr_type)>::type label_b_type_ref;
                typedef typename std::remove_reference<label_b_type_ref>::type label_b_type_check;
                typedef typename label_b_type_check::particles_type particles_b_type_check;

                static_assert(std::is_same<particles_b_type_check,particles_b_type>::value, "particles type b in expression not the same as particles type given to block");

                const particles_b_type& b = std::get<1>(block).get_particles();
                const expr_type expr = std::get<2>(block);
                const if_expr_type if_expr = std::get<3>(block);

                const Index nb = b.size();

                double sum = 0;
                for (size_t j=0; j<nb; ++j) {
                    typename particles_b_type::const_reference bj = b[j];
                    if (eval(if_expr,bj)) {
                        sum += eval(expr,bj)*rhs(j);
                    }
                }
                y(0) += sum;
            }

            template <typename Dest, typename Source, typename particles_a_type,
                                     typename expr_type, typename if_expr_type>
            void evalTo_block(Eigen::VectorBlock<Dest> y, const Eigen::VectorBlock<Source>& rhs, const std::tuple<const particles_a_type&,const one&,expr_type,if_expr_type>& block) const {

                typedef typename particles_a_type::double_d double_d;
                typedef typename particles_a_type::position position;
                const static unsigned int dimension = particles_a_type::dimension;

                typedef typename std::result_of<detail::univariate_expr(expr_type)>::type label_a_type_ref;
                typedef typename std::remove_reference<label_a_type_ref>::type label_a_type_check;
                typedef typename label_a_type_check::particles_type particles_a_type_check;

                static_assert(std::is_same<particles_a_type_check,particles_a_type>::value, "particles type a in expression not the same as particles type given to block");

                const particles_a_type& a = std::get<0>(block).get_particles();
                const expr_type expr = std::get<2>(block);
                const if_expr_type if_expr = std::get<3>(block);

                const Index na = a.size();

                for (size_t i=0; i<na; ++i) {
                    typename particles_a_type::const_reference ai = a[i];
                    if (eval(if_expr,ai)) {
                        y(i) += eval(expr,ai)*rhs(0);
                    }
                }
            }


            template<typename Dest>
            void evalTo_impl2(Dest& y) const {}
             
            template<typename Dest, typename I, typename J, typename T1, typename ... T>
            void evalTo_impl2(Dest& y, const std::tuple<I,J,T1>& block, const T&... other_blocks) const {
                evalTo_block(y.segment(m_matrix.template start_col<J::value>(),m_matrix.template size_col<J::value>()),
                        m_rhs.segment(m_matrix.template start_row<I::value>(),m_matrix.template size_row<I::value>()),
                        std::get<2>(block));
                evalTo_impl2(y,other_blocks...);
            }

            template<typename Dest, std::size_t... I>
            void evalTo_impl1(Dest& y, index_sequence<I...>) const {
                evalTo_impl2(y, std::make_tuple(mpl::int_<I/NJ>(), mpl::int_<I%NJ>(), std::get<I>(m_matrix.m_blocks))...);
            }

            // This function is automatically called by Eigen. It must evaluate the product of matrix * rhs into y.
            template<typename Dest>
            void evalTo(Dest& y) const {
                y.setZero();
                evalTo_impl1(y, make_index_sequence<NI*NJ>());
            }
            
        protected:
            const MatrixReplacement<NI,NJ,Blocks>& m_matrix;
            typename Rhs::Nested m_rhs;
    };






    
    template <typename A, unsigned int A_depth, 
              typename B, unsigned int B_depth, 
              typename Expr, typename IfExpr=detail::SymbolicExpr<typename proto::terminal<bool>::type>>
    MatrixReplacement<1,1,std::tuple<std::tuple<const A&,const B&,Expr,IfExpr>>> 
    create_eigen_operator(
            const Label<A_depth,A>& a, 
            const Label<B_depth,B>& b, 
            const Expr& expr, 
            const IfExpr& if_expr = IfExpr(proto::terminal<bool>::type::make(true))) 
    {
        return MatrixReplacement<1,1,std::tuple<std::tuple<const A&, const B&, Expr,IfExpr>>>(
                std::make_tuple(
                    std::tie(proto::value(a).get_particles(),proto::value(b).get_particles(),expr,if_expr)
                ));
    }

    template <unsigned int NI, unsigned NJ, typename ... T>
    MatrixReplacement<NI,NJ,std::tuple<T...>> create_block_eigen_operator(const MatrixReplacement<1,1,std::tuple<T>>&... operators) {
        return MatrixReplacement<NI,NJ,std::tuple<T...>>(std::make_tuple(std::get<0>(operators.m_blocks)...));
    }


#endif

}

#endif //OPERATORS_H_
