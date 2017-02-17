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
#include "detail/Evaluate.h"
#include "Symbolic.h"

#ifdef HAVE_EIGEN
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#endif

namespace Aboria {
namespace detail {
#ifdef HAVE_EIGEN
    template <unsigned int NI, unsigned int NJ, typename Blocks>
    class MatrixReplacement;
}
}

namespace Eigen {
    namespace internal {
        // MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
        template<unsigned int NI, unsigned int NJ, typename Blocks>
        struct traits<Aboria::MatrixReplacement<NI, NJ, Blocks>> :  public Eigen::internal::traits<Eigen::SparseMatrix<double> > {};
    }

template<typename T,unsigned int N>
struct NumTraits<Aboria::Vector<T,N>> {
  typedef Aboria::Vector<T,N> Scalar;
  typedef Aboria::Vector<T,N> Real;
  typedef Aboria::Vector<T,N> NonInteger;
  typedef Aboria::Vector<T,N> Literal;
  typedef Aboria::Vector<T,N> Nested;
  enum {
    IsComplex = 0,
    IsInteger = 0,
    IsSigned = 1,
    RequireInitialization = 0,
    ReadCost = N,
    AddCost = N,
    MulCost = N
  };

  inline static Real epsilon() { return Real(std::numeric_limits<T>::epsilon()); }
  inline static Real dummy_precision() { return Real(std::numeric_limits<T>::epsilon()); }
  inline static Scalar highest() { return Scalar(std::numeric_limits<T>::max()); }
  inline static Scalar lowest() { return Scalar(std::numeric_limits<T>::lowest()); }
  inline static int digits10() { return N*std::numeric_limits<T>::digits10;}
};
}

namespace Aboria {
namespace detail {

    template <typename MatrixReplacement, typename Scalar=MatrixReplacement::Scalar>
    class InnerIterator {
        public:
        typedef const Scalar* pointer;
        typedef std::forward_iterator_tag iterator_category;
        typedef const Scalar& reference;
        typedef const Scalar value_type;
        typedef std::ptrdiff_t difference_type;

        InnerIterator(const MatrixReplacement& mat, const Index row):
            m_mat(mat),m_row(row),m_col(0) {};

        InnerIterator(const InnerIterator& other):
            m_mat(other.m_mat),m_row(other.m_row),m_col(other.m_col) {};

        Scalar value() const {
            return dereference();
        }

        Index index() const {
            return m_col;
        }

        operator bool() const {
            return (m_row < m_mat.rows()) && (m_col < m_mat.cols());
        }

        bool equal(InnerIterator const& other) const {
            return (m_row == other.m_row)&&(m_col == other.m_col);
        }

        Scalar dereference() const {
            return m_mat.coeff(m_row,m_col);
        }

        void increment() {
            m_col++;
            ASSERT(m_col < m_mat.cols(),"InnerIterator outside cols range");
        }

        Scalar operator *() {
            return dereference();
        }
        Scalar operator ->() {
            return dereference();
        }
        InnerIterator& operator++() {
            increment();
            return *this;
        }
        InnerIterator operator++(int) {
            InnerIterator tmp(*this);
            operator++();
            return tmp;
        }

        size_t operator-(InnerIterator start) const {
            ASSERT(m_row == start.m_row,"Difference between InnerIterators must have identical row numbers");
            return (m_col-start.m_col);
        }

        inline bool operator==(const InnerIterator& rhs) {
            return equal(rhs);
        }

        inline bool operator!=(const InnerIterator& rhs){
            return !operator==(rhs);
        }

        private:


        friend class boost::iterator_core_access;
        const Index m_row;
        Index m_col;
        const MatrixReplacement& m_mat;
    };
    


    template<typename T1=void>
    int sum() {
        return 0;
    }

    template<typename T1, typename... T>
    T1 sum(T1 s, T... ts) {
        return s + sum(ts...);
    }


    // Inheriting EigenBase should not be needed in the future.
    template <unsigned int NI, unsigned int NJ, typename Blocks>
    class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement<NI,NJ,Blocks>> {

            typedef typename tuple_ns::tuple_element<0,Blocks>::type first_block_type;
            typedef typename tuple_ns::tuple_element<2,first_block_type>::type first_expr_type;
            
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
            typedef InnerIterator<MatrixReplacement> InnerIterator;

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

    template <typename Dest, typename Source,
              typename particles_a_type, typename particles_b_type,
              typename expr_type>
    void evalTo_block(Eigen::VectorBlock<Dest> y, const Eigen::VectorBlock<Source>& rhs, const tuple_ns::tuple<const particles_a_type&,const particles_b_type&,expr_type>& block) {
        const particles_a_type& a = tuple_ns::get<0>(block);
        const particles_b_type& b = tuple_ns::get<1>(block);
        const expr_type& expr = tuple_ns::get<2>(block);

        expr.evaluate(a,b,y,rhs);
    }

    template<typename Dest, unsigned int NI, unsigned int NJ, typename Blocks, typename Rhs>
    void evalTo_unpack_blocks(Dest& y, const MatrixReplacement<NI,NJ,Blocks>& lhs, const Rhs& rhs) {}

    template<typename Dest, unsigned int NI, unsigned int NJ, typename Blocks, typename Rhs, typename I, typename J, typename T1, typename ... T>
    void evalTo_unpack_blocks(Dest& y, const MatrixReplacement<NI,NJ,Blocks>& lhs, const Rhs& rhs, const tuple_ns::tuple<I,J,T1>& block, const T&... other_blocks) {
        evalTo_block(y.segment(lhs.template start_row<I::value>(),lhs.template size_row<I::value>()),
                rhs.segment(lhs.template start_col<J::value>(),lhs.template size_col<J::value>()),
                tuple_ns::get<2>(block));
        evalTo_unpack_blocks(y,lhs,rhs,other_blocks...);
    }

    template<typename Dest, unsigned int NI, unsigned int NJ, typename Blocks, typename Rhs, std::size_t... I>
    void evalTo_impl(Dest& y, const MatrixReplacement<NI,NJ,Blocks>& lhs, const Rhs& rhs, detail::index_sequence<I...>) {
        evalTo_unpack_blocks(y,lhs,rhs,tuple_ns::make_tuple(mpl::int_<I/NJ>(), mpl::int_<I%NJ>(), tuple_ns::get<I>(lhs.m_blocks))...);
    }
}


// Implementation of MatrixReplacement * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen {
namespace internal {
template<typename Rhs, unsigned int NI, unsigned int NJ, typename Blocks>
struct generic_product_impl<Aboria::MatrixReplacement<NI,NJ,Blocks>, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
    : generic_product_impl_base<Aboria::MatrixReplacement<NI,NJ,Blocks>,Rhs,generic_product_impl<Aboria::MatrixReplacement<NI,NJ,Blocks>,Rhs> > {

    typedef typename Product<Aboria::MatrixReplacement<NI,NJ,Blocks>,Rhs>::Scalar Scalar;
    template<typename Dest>
    CUDA_HOST_DEVICE
    static void scaleAndAddTo(Dest& y, const Aboria::MatrixReplacement<NI,NJ,Blocks>& lhs, const Rhs& rhs, const Scalar& alpha) {
        // This method should implement "y += alpha * lhs * rhs" inplace,
        // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
#ifdef __CUDA_ARCH__
        ERROR_CUDA("MatrixReplacement class unusable from device code");
#else
        assert(alpha==Scalar(1) && "scaling is not implemented");
        evalTo_impl(y, lhs, rhs, Aboria::detail::make_index_sequence<NI*NJ>());
#endif
    }
};
}
}

namespace Aboria {
    /// creates a matrix-free linear operator for use with Eigen
    ///
    template <typename A,
              typename B,
              typename Expr,
              typename BlockType=tuple_ns::tuple<const A&,const B&,const Expr&>,
              typename TupleType=tuple_ns::tuple<BlockType>
    MatrixReplacement<1,1,TupleType>
    create_operator(
            const A& a,
            const B& b,
            const Expr& expr)
    {
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


#endif

}

#endif //OPERATORS_H_
