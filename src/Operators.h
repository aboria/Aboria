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

namespace Aboria {
    template <unsigned int NI, unsigned int NJ, typename Blocks>
    class MatrixReplacement;
}

#ifdef HAVE_EIGEN
#include "detail/Operators.h"


namespace Aboria {

// Inheriting EigenBase should not be needed in the future.
    template <unsigned int NI, unsigned int NJ, typename Blocks>
    class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement<NI,NJ,Blocks>> {

            typedef typename tuple_ns::tuple_element<0,Blocks>::type first_block_type;
            typedef typename first_block_type::position position;
            
        public:

            // Expose some compile-time information to Eigen:
            typedef typename first_expr_type::Scalar Scalar;
            typedef Scalar RealScalar;
            typedef size_t Index;
            typedef int StorageIndex;
            enum {
                ColsAtCompileTime = Eigen::Dynamic,
                RowsAtCompileTime = Eigen::Dynamic,
                MaxColsAtCompileTime = Eigen::Dynamic,
                MaxRowsAtCompileTime = Eigen::Dynamic,
                IsRowMajor = false
            };
            typedef detail::InnerIterator<MatrixReplacement> InnerIterator;

            MatrixReplacement(const Blocks& blocks):m_blocks(blocks) {};

            CUDA_HOST_DEVICE
            Index rows() const {
                //std::cout << "rows = " << rows_impl(detail::make_index_sequence<NI>()) << std::endl;
                //
#ifdef __CUDA_ARCH__
                ERROR_CUDA("MatrixReplacement class unusable from device code");
                return 0;
#else
                return rows_impl(detail::make_index_sequence<NI>());
#endif
            }

            CUDA_HOST_DEVICE
            Index cols() const {
                //std::cout << "cols = " << cols_impl(detail::make_index_sequence<NJ>())<< std::endl;
#ifdef __CUDA_ARCH__
                ERROR_CUDA("MatrixReplacement class unusable from device code");
                return 0;
#else
                return cols_impl(detail::make_index_sequence<NJ>());
#endif
            }

            CUDA_HOST_DEVICE
            Index innerSize() const {
#ifdef __CUDA_ARCH__
                ERROR_CUDA("MatrixReplacement class unusable from device code");
                return 0;
#else
                return rows();
#endif
            }

            CUDA_HOST_DEVICE
            Index outerSize() const {
#ifdef __CUDA_ARCH__
                ERROR_CUDA("MatrixReplacement class unusable from device code");
                return 0;
#else
                return cols();
#endif
            }

            void resize(Index a_rows, Index a_cols) {
                // This method should not be needed in the future.
                assert(a_rows==0 && a_cols==0 || a_rows==rows() && a_cols==cols());
            }

            Scalar coeff(const Index i, const Index j) const {
                return coeff_impl(i,j,detail::make_index_sequence<NI*NJ>());
            }

            template<typename Rhs>
            Eigen::Product<MatrixReplacement,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
                return Eigen::Product<MatrixReplacement,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
            }

            template<typename Derived>
            void assemble(Eigen::DenseBase<Derived>& matrix) const {
                const size_t na = rows();
                const size_t nb = cols();
                matrix.resize(na,nb);
                CHECK((matrix.rows() == na) && (matrix.cols() == nb), "matrix size is not compatible with expression.");
                assemble_impl(matrix,detail::make_index_sequence<NI*NJ>());
            }

            void assemble(Eigen::SparseMatrix<Scalar>& matrix) {
                const size_t na = rows();
                const size_t nb = cols();
                //matrix.resize(na,nb);
                CHECK((matrix.rows() == na) && (matrix.cols() == nb), "matrix size is not compatible with expression.");

                typedef Eigen::Triplet<Scalar> triplet_type;
                std::vector<triplet_type> tripletList;
                // TODO: can we estimate this better?
                tripletList.reserve(na*5);

                assemble_impl(tripletList,detail::make_index_sequence<NI*NJ>());

                matrix.setFromTriplets(tripletList.begin(),tripletList.end());
            }

            

        private:
            template<std::size_t... I>
            Index rows_impl(detail::index_sequence<I...>) const {
                return detail::sum(tuple_ns::get<0>(tuple_ns::get<I*NJ>(m_blocks)).size()...);
            }

            template<std::size_t... J>
            Index cols_impl(detail::index_sequence<J...>) const {
                return detail::sum(tuple_ns::get<1>(tuple_ns::get<J>(m_blocks)).size()...);
            }

            template <int I>
            Index start_col() const {
                return cols_impl(detail::make_index_sequence<I>());
            }

            template <int I>
            Index size_col() const {
                return tuple_ns::get<1>(tuple_ns::get<I>(m_blocks)).size();
            }

            template <int I>
            Index start_row() const {
                return rows_impl(detail::make_index_sequence<I>());
            }

            template <int I>
            Index size_row() const {
                return tuple_ns::get<0>(tuple_ns::get<I*NJ>(m_blocks)).size();
            }

            template <
              typename particles_a_type, typename particles_b_type,
              typename expr_type>
            Scalar coeff_impl_block(const Index i, const Index j, const tuple_ns::tuple<const particles_a_type&,const particles_b_type&,expr_type>& block) const {

                ASSERT(i>=0, "i less than zero");
                ASSERT(j>=0, "j less than zero");
                const particles_a_type& a = tuple_ns::get<0>(block);
                const particles_b_type& b = tuple_ns::get<1>(block);
                const expr_type& expr = tuple_ns::get<2>(block);
                ASSERT(i < a.size(),"i greater than a.size()");
                ASSERT(j < b.size(),"j greater than b.size()");
                typename particles_a_type::const_reference ai = a[i];
                typename particles_b_type::const_reference bj = b[j];
                typename particles_a_type::double_d dx = get<position>(bj)-get<position>(ai);
                return expr.eval(dx,ai,bj);
            }

            template<std::size_t... I>
            Scalar coeff_impl(const Index i, const Index j, detail::index_sequence<I...>) const {

                return detail::sum(
                        ((i>=start_row<I/NJ>())&&(i<start_row<I/NJ+1>())
                         &&
                         (j>=start_col<I%NJ>())&&(j<start_col<I%NJ+1>()))?
                            (coeff_impl_block(i-start_row<I/NJ>(),
                                              j-start_col<I%NJ>(),
                                              tuple_ns::get<I>(m_blocks))):
                            (0.0)...
                        );
            }

            template <typename particles_a_type, typename particles_b_type,
                      typename expr_type>
            void assemble_block_impl(const size_t startI, const size_t startJ,
                    std::vector<Eigen::Triplet<Scalar>>& triplets,
                    const tuple_ns::tuple<
                            const particles_a_type&,
                            const particles_b_type&,
                            expr_type>& block) const {
                typedef typename particles_b_type::double_d double_d;
                typedef typename particles_b_type::position position;
                const particles_a_type& a = tuple_ns::get<0>(block);
                const particles_b_type& b = tuple_ns::get<1>(block);
                const expr_type& expr = tuple_ns::get<2>(block);
                
                expr.assemble(a,b,triplets,startI,startJ);
            }

            template <typename particles_a_type, typename particles_b_type,
                      typename expr_type, typename if_expr_type, typename Derived>
            void assemble_block_impl(
                    const Eigen::MatrixBase<Derived> &matrix,
                    const tuple_ns::tuple<
                            const particles_a_type&,
                            const particles_b_type&,
                            expr_type,if_expr_type>& block) const {
                typedef typename particles_b_type::double_d double_d;
                typedef typename particles_b_type::position position;
                const particles_a_type& a = tuple_ns::get<0>(block);
                const particles_b_type& b = tuple_ns::get<1>(block);
                expr_type expr = tuple_ns::get<2>(block);
                expr.assemble(a,b,matrix);
            }

            template<std::size_t... I>
            void assemble_impl(std::vector<Eigen::Triplet<Scalar>>& triplets, detail::index_sequence<I...>) const {
                int dummy[] = { 0, (
                        assemble_block_impl(
                            start_row<I/NJ>(),start_col<I%NJ>(),
                            triplets,tuple_ns::get<I>(m_blocks)),void(),0)... };
                static_cast<void>(dummy);
            }


            template<typename Derived, std::size_t... I>
            void assemble_impl(Eigen::DenseBase<Derived>& matrix, detail::index_sequence<I...>) const {
                int dummy[] = { 0, (
                        assemble_block_impl(
                            matrix.block(start_row<I/NJ>(),start_col<I%NJ>()
                            ,start_row<I/NJ+1>()-start_row<I/NJ>()
                            ,start_col<I%NJ+1>()-start_col<I%NJ>())
                            ,tuple_ns::get<I>(m_blocks)),void(),0)... };
                static_cast<void>(dummy);
            }

            const Blocks m_blocks;

    };

/// creates a matrix-free linear operator for use with Eigen
///
template <typename A,
          typename B,
          typename Expr,
          typename BlockType=tuple_ns::tuple<const A&,const B&,const Expr&>,
          typename TupleType=tuple_ns::tuple<BlockType>>
MatrixReplacement<1,1,TupleType>
create_operator(
        const A& a,
        const B& b,
        const Expr& expr) {
    return MatrixReplacement<1,1,TupleType>(
            tuple_ns::make_tuple(
                tuple_ns::make_tuple(a,b,expr)
            ));
}

/// creates a matrix-free linear block operator for use with Eigen
///
template <unsigned int NI, unsigned int NJ, typename ... T>
MatrixReplacement<NI,NJ,tuple_ns::tuple<typename tuple_ns::tuple_element<0,T>::type...>>
create_block_operator(const MatrixReplacement<1,1,T>&... operators) {
    typedef tuple_ns::tuple<typename tuple_ns::tuple_element<0,T>::type...> tuple_type;
    return MatrixReplacement<NI,NJ,tuple_type>(tuple_type(tuple_ns::get<0>(operators.m_blocks)...));
}



}
#endif //HAVE_EIGEN

#endif //OPERATORS_H_
