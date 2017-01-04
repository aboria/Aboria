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
#include <unsupported/Eigen/IterativeSolvers>
#endif

namespace Aboria {
#ifdef HAVE_EIGEN
    template <unsigned int NI, unsigned int NJ, typename Blocks> 
    class MatrixReplacement;
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

    struct OneDummy {
        struct vector_tag {
            size_t size() { return 1; }
        };
        vector_tag get_particles() { return d; }

        vector_tag d;
    };

    struct SizeOne {
        size_t size() const { return 1; }
    };

    struct One {
        const SizeOne& get_size_one() const { return one; }
        SizeOne one;
    };

    namespace detail {

    int sum() {
        return 0;
    }

    template<typename T1, typename... T>
    T1 sum(T1 s, T... ts) {
        return s + sum(ts...);
    }

    }

    // Inheriting EigenBase should not be needed in the future.
    template <unsigned int NI, unsigned int NJ, typename Blocks> 
    class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement<NI,NJ,Blocks>> {

            typedef typename tuple_ns::tuple_element<0,Blocks>::type first_block_type;
            typedef typename tuple_ns::tuple_element<2,first_block_type>::type first_expr_type;
        public:
    
            // Expose some compile-time information to Eigen:
            typedef typename detail::symbolic_helper<first_expr_type>::result_base_type Scalar;
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

            MatrixReplacement(const Blocks& blocks):m_blocks(blocks) {};

            template<std::size_t... I>
            Index rows_impl(detail::index_sequence<I...>) const {
                return detail::sum(tuple_ns::get<0>(tuple_ns::get<I*NJ>(m_blocks)).size()...);
            }

            template<std::size_t... J>
            Index cols_impl(detail::index_sequence<J...>) const {
                return detail::sum(tuple_ns::get<1>(tuple_ns::get<J>(m_blocks)).size()...);
            }

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

            void resize(Index a_rows, Index a_cols) {
                // This method should not be needed in the future.
                assert(a_rows==0 && a_cols==0 || a_rows==rows() && a_cols==cols());
            }

            template < 
              typename particles_a_type, typename particles_b_type,
              typename expr_type, typename if_expr_type>
            Scalar coeff_impl_block(const Index i, const Index j, const tuple_ns::tuple<const particles_a_type&,const particles_b_type&,expr_type,if_expr_type>& block) const {
                
                //TODO: Have to copy the expressions (I think), since proto returns a non-const reference
                //to the stored constants, and if I want this to be const
                expr_type expr = tuple_ns::get<2>(block);
                if_expr_type if_expr = tuple_ns::get<3>(block);

                if (is_trivially_zero(expr) || is_trivially_false(if_expr)) {
                    return 0;
                } else {
                    typedef typename particles_b_type::double_d double_d;
                    typedef typename particles_b_type::position position;

                    const particles_a_type& a = tuple_ns::get<0>(block);
                    const particles_b_type& b = tuple_ns::get<1>(block);
                    typename particles_a_type::const_reference ai = a[i];
                    typename particles_b_type::const_reference bj = b[j];
                    const double_d dx = a.correct_dx_for_periodicity(get<position>(bj)-get<position>(ai));

                    if (eval(if_expr,dx,ai,bj)) {
                        return eval(expr,dx,ai,bj);
                    } else {
                        return 0;
                    }
                }

            }



            template < 
              typename particles_b_type,
              typename expr_type, typename if_expr_type>
            Scalar coeff_impl_block(const Index i, const Index j, const tuple_ns::tuple<const SizeOne&,const particles_b_type&,expr_type,if_expr_type>& block) const {
             
                //TODO: Have to copy the expressions (I think), since proto returns a non-const reference
                //to the stored constants, and if I want this to be const
                expr_type expr = tuple_ns::get<2>(block);
                if_expr_type if_expr = tuple_ns::get<3>(block);
                
                if (is_trivially_zero(expr) || is_trivially_false(if_expr)) {
                    return 0;
                } else {
                    typedef typename particles_b_type::double_d double_d;
                    typedef typename particles_b_type::position position;

                    const particles_b_type& b = tuple_ns::get<1>(block);
                    typename particles_b_type::const_reference bj = b[j];

                    
                    if (eval(if_expr,bj)) {
                        return eval(expr,bj);
                    } else {
                        return 0;
                    }
                }

            }

            template < 
              typename particles_a_type,
              typename expr_type, typename if_expr_type>
            Scalar coeff_impl_block(const Index i, const Index j, const tuple_ns::tuple<const particles_a_type&,const SizeOne&,expr_type,if_expr_type>& block) const {
 
                //TODO: Have to copy the expressions (I think), since proto returns a non-const reference
                //to the stored constants, and if I want this to be const
                expr_type expr = tuple_ns::get<2>(block);
                if_expr_type if_expr = tuple_ns::get<3>(block);
                             

                if (is_trivially_zero(expr) || is_trivially_false(if_expr)) {
                    return 0;
                } else {
                    typedef typename particles_a_type::double_d double_d;
                    typedef typename particles_a_type::position position;

                    const particles_a_type& a = tuple_ns::get<0>(block);
                    typename particles_a_type::const_reference ai = a[i];
       
                    if (eval(if_expr,ai)) {
                        return eval(expr,ai);
                    } else {
                        return 0;
                    }
                }

            }

            template < 
              typename expr_type, typename if_expr_type>
            Scalar coeff_impl_block(const Index i, const Index j, const tuple_ns::tuple<const SizeOne&,const SizeOne&,expr_type,if_expr_type>& block) const {
                //TODO: Have to copy the expressions (I think), since proto returns a non-const reference
                //to the stored constants, and if I want this to be const
                expr_type expr = tuple_ns::get<2>(block);
                if_expr_type if_expr = tuple_ns::get<3>(block);
                 
                if (eval(if_expr)) {
                    return eval(expr);
                } else {
                    return 0;
                }

            }

            template<std::size_t... I>
            Scalar coeff_impl(const Index i, const Index j, detail::index_sequence<I...>) const {
                return detail::sum(
                        ((i>=start_row<I>())&&(i<start_row<I+1>()))?
                            (coeff_impl_block(i,j,tuple_ns::get<I*NJ>(m_blocks))):
                            (0.0)...
                        );
            }

            Scalar coeff(const Index i, const Index j) const {
                return coeff_impl(i,j,detail::make_index_sequence<NI>());
            }

            template<typename Rhs>
            Eigen::Product<MatrixReplacement,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
                return Eigen::Product<MatrixReplacement,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
            }

            template<typename Derived>
            void assemble(Eigen::DenseBase<Derived>& matrix) {
                const size_t na = rows();
                const size_t nb = cols();
                matrix.resize(na,nb);
                CHECK((matrix.rows() == na) && (matrix.cols() == nb), "matrix size is not compatible with expression.");
                for (size_t i=0; i<na; ++i) {
                    for (size_t j=0; j<nb; ++j) {
                        matrix(i,j) = coeff(i,j);
                    }
                }
            }

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


            const Blocks m_blocks;

    };

    template <typename Dest, typename Source, 
              typename particles_a_type, typename particles_b_type,
              typename expr_type, typename if_expr_type>
    void evalTo_block(Eigen::VectorBlock<Dest> y, const Eigen::VectorBlock<Source>& rhs, const tuple_ns::tuple<const particles_a_type&,const particles_b_type&,expr_type,if_expr_type>& block) {
        typedef typename particles_b_type::double_d double_d;
        typedef typename particles_b_type::position position;
        const particles_a_type& a = tuple_ns::get<0>(block);
        const particles_b_type& b = tuple_ns::get<1>(block);

        //TODO: Have to copy the expressions (I think), since proto returns a non-const reference
        //to the stored constants, and the eval product expression in eigen assumes that your matrix
        //replacement is a constant
        expr_type expr = tuple_ns::get<2>(block);
        if_expr_type if_expr = tuple_ns::get<3>(block);

        const size_t na = a.size();
        const size_t nb = b.size();


        if (is_trivially_zero(expr) || is_trivially_false(if_expr)) {
            //std::cout << "zero a x b block" <<std::endl;
            return;
        }

        if (is_trivially_true(if_expr)) {
            //std::cout << "dense a x b block" <<std::endl;
            ASSERT(!a.get_periodic().any(),"periodic does not work with dense");
            
            const size_t parallel_size = 20;
            const size_t block_size = 20;
            if (na > parallel_size) {
                #pragma omp parallel for
                for (size_t i=0; i<na; ++i) {
                    typename particles_a_type::const_reference ai = a[i];
                    double sum = 0;
                    for (size_t j=0; j<na; ++j) {
                        typename particles_b_type::const_reference bj = b[j];
                        sum += eval(expr,get<position>(bj)-get<position>(ai),ai,bj)*rhs(j);
                    }
                    y(i) += sum;
                }
            } else {
                for (size_t i=0; i<na; ++i) {
                    typename particles_a_type::const_reference ai = a[i];
                    double sum = 0;
                    for (size_t j=0; j<na; ++j) {
                        typename particles_b_type::const_reference bj = b[j];
                        //std::cout << "a = "<<get<position>(ai)<<" b = "<<get<position>(bj)<<std::endl;
                        //std::cout << "using dx = "<<dx<<" rhs(j) = "<<rhs(j)<<" eval = "<<eval(expr,dx,ai,bj)<<std::endl;
                        sum += eval(expr,get<position>(bj)-get<position>(ai),ai,bj)*rhs(j);
                    }
                    y(i) += sum;
                }
            }
        } else {
            //std::cout << "sparse a x b block" <<std::endl;
            #pragma omp parallel for
            for (size_t i=0; i<na; ++i) {
                typename particles_a_type::const_reference ai = a[i];
                double sum = 0;
                //std::cout << "evaluating fucntion for particle at "<<get<position>(ai)<<std::endl;
                for (auto pairj: b.get_neighbours(get<position>(ai))) {
                    const double_d & dx = tuple_ns::get<1>(pairj);
                    typename particles_b_type::const_reference bj = tuple_ns::get<0>(pairj);
                    //std::cout << "looking at particle with dx = "<<dx<<std::endl;
                    const size_t j = &get<position>(bj) - get<position>(b).data();
                    if (eval(if_expr,dx,ai,bj)) {
                        //std::cout <<"if expression is true. eval = "<<eval(expr,dx,ai,bj)<<std::endl;
                        sum += eval(expr,dx,ai,bj)*rhs(j);
                    }
                }
                y(i) += sum;
            }
        }
    }


    template <typename Dest, typename Source, typename particles_b_type,
              typename expr_type, typename if_expr_type>
    void evalTo_block(Eigen::VectorBlock<Dest> y, const Eigen::VectorBlock<Source>& rhs, const tuple_ns::tuple<const SizeOne&,const particles_b_type&,expr_type,if_expr_type>& block) {

        typedef typename particles_b_type::double_d double_d;
        typedef typename particles_b_type::position position;

        const particles_b_type& b = tuple_ns::get<1>(block);
        expr_type expr = tuple_ns::get<2>(block);
        if_expr_type if_expr = tuple_ns::get<3>(block);

        const size_t nb = b.size();

        if (is_trivially_zero(expr) || is_trivially_false(if_expr)) {
            //std::cout << "zero one x b block" <<std::endl;
            return;
        }

        double sum = 0;
        //std::cout << "dens one x b block" <<std::endl;
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
    void evalTo_block(Eigen::VectorBlock<Dest> y, const Eigen::VectorBlock<Source>& rhs, const tuple_ns::tuple<const particles_a_type&,const SizeOne&,expr_type,if_expr_type>& block) {

        typedef typename particles_a_type::double_d double_d;
        typedef typename particles_a_type::position position;

        const particles_a_type& a = tuple_ns::get<0>(block);
        expr_type expr = tuple_ns::get<2>(block);
        if_expr_type if_expr = tuple_ns::get<3>(block);

        const size_t na = a.size();

        if (is_trivially_zero(expr) || is_trivially_false(if_expr)) {
            //std::cout << "zero a x one block" <<std::endl;
            return;
        }

        //std::cout << "dens a x one block" <<std::endl;
        for (size_t i=0; i<na; ++i) {
            typename particles_a_type::const_reference ai = a[i];
            if (eval(if_expr,ai)) {
                y(i) += eval(expr,ai)*rhs(0);
            }
        }
    }

    template <typename Dest, typename Source,
                             typename expr_type, typename if_expr_type>
    void evalTo_block(Eigen::VectorBlock<Dest> y, const Eigen::VectorBlock<Source>& rhs, const tuple_ns::tuple<const SizeOne&,const SizeOne&,expr_type,if_expr_type>& block) {

        expr_type expr = tuple_ns::get<2>(block);
        if_expr_type if_expr = tuple_ns::get<3>(block);

        if (is_trivially_zero(expr) || is_trivially_false(if_expr)) {
            //std::cout << "zero one x one block" <<std::endl;
            return;
        }

        proto::display_expr(expr);
        //std::cout << "dens one x one block" <<std::endl;
        y(0) += eval(expr)*rhs(0);
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
    template <typename A, unsigned int A_depth, 
              typename B, unsigned int B_depth, 
              typename Expr, typename IfExpr=detail::SymbolicExpr<typename proto::terminal<bool>::type>>
    MatrixReplacement<1,1,tuple_ns::tuple<tuple_ns::tuple<const A&,const B&,
        typename detail::symbolic_helper<Expr>::deep_copy_type,
        typename detail::symbolic_helper<IfExpr>::deep_copy_type>>>
    create_eigen_operator(
            const Label<A_depth,A>& a, 
            const Label<B_depth,B>& b, 
            const Expr& expr, 
            const IfExpr& if_expr = IfExpr(proto::terminal<bool>::type::make(true))) 
    {
        typedef tuple_ns::tuple<const A&, const B&, 
            typename detail::symbolic_helper<Expr>::deep_copy_type,
            typename detail::symbolic_helper<IfExpr>::deep_copy_type> block_type;
        return MatrixReplacement<1,1,tuple_ns::tuple<block_type>>(
                tuple_ns::make_tuple(
                    block_type(proto::value(a).get_particles(),proto::value(b).get_particles(),
                        deep_copy(expr),
                        deep_copy(if_expr))
                ));
    }

    template <typename B, unsigned int B_depth, 
              typename Expr, typename IfExpr=detail::SymbolicExpr<typename proto::terminal<bool>::type>>
    MatrixReplacement<1,1,tuple_ns::tuple<tuple_ns::tuple<const SizeOne&,const B&,
        typename detail::symbolic_helper<Expr>::deep_copy_type,
        typename detail::symbolic_helper<IfExpr>::deep_copy_type>>>
    create_eigen_operator(
            const One& a, 
            const Label<B_depth,B>& b, 
            const Expr& expr, 
            const IfExpr& if_expr = IfExpr(proto::terminal<bool>::type::make(true))) 
    {
        typedef tuple_ns::tuple<const SizeOne&, const B&, 
            typename detail::symbolic_helper<Expr>::deep_copy_type,
            typename detail::symbolic_helper<IfExpr>::deep_copy_type> block_type;
        return MatrixReplacement<1,1,tuple_ns::tuple<block_type>>(
                tuple_ns::make_tuple(
                    block_type(a.get_size_one(),proto::value(b).get_particles(),
                        deep_copy(expr),
                        deep_copy(if_expr))
                ));
    }

    template <typename A, unsigned int A_depth, 
              typename Expr, typename IfExpr=detail::SymbolicExpr<typename proto::terminal<bool>::type>>
    MatrixReplacement<1,1,tuple_ns::tuple<tuple_ns::tuple<const A&,const SizeOne&,
        typename detail::symbolic_helper<Expr>::deep_copy_type,
        typename detail::symbolic_helper<IfExpr>::deep_copy_type>>>
    create_eigen_operator(
            const Label<A_depth,A>& a, 
            const One& b, 
            const Expr& expr, 
            const IfExpr& if_expr = IfExpr(proto::terminal<bool>::type::make(true))) 
    {
        typedef tuple_ns::tuple<const A&, const SizeOne&, 
            typename detail::symbolic_helper<Expr>::deep_copy_type,
            typename detail::symbolic_helper<IfExpr>::deep_copy_type> block_type;
        return MatrixReplacement<1,1,tuple_ns::tuple<block_type>>(
                tuple_ns::make_tuple(
                    block_type(proto::value(a).get_particles(),b.get_size_one(),
                        deep_copy(expr),
                        deep_copy(if_expr))
                ));
    }

    template <typename Expr, typename IfExpr=detail::SymbolicExpr<typename proto::terminal<bool>::type>>
    MatrixReplacement<1,1,tuple_ns::tuple<tuple_ns::tuple<const SizeOne&,const SizeOne&,
        typename detail::symbolic_helper<Expr>::deep_copy_type,
        typename detail::symbolic_helper<IfExpr>::deep_copy_type>>>
    create_eigen_operator(
            const One& a, 
            const One& b, 
            const Expr& expr, 
            const IfExpr& if_expr = IfExpr(proto::terminal<bool>::type::make(true))) 
    {
        typedef tuple_ns::tuple<const SizeOne&, const SizeOne&, 
            typename detail::symbolic_helper<Expr>::deep_copy_type,
            typename detail::symbolic_helper<IfExpr>::deep_copy_type> block_type;
        return MatrixReplacement<1,1,tuple_ns::tuple<block_type>>(
                tuple_ns::make_tuple(
                    block_type(a.get_size_one(),b.get_size_one(),
                        deep_copy(expr),
                        deep_copy(if_expr))
                ));
    }
    /*
    create_block_eigen_operator_impl(const Blocks&... blocks) {
        typedef 
        return MatrixReplacement<NI,NJ,tuple_ns::tuple<Blocks...>>(
        */
        
    template <unsigned int NI, unsigned int NJ, typename ... T>
    MatrixReplacement<NI,NJ,tuple_ns::tuple<typename tuple_ns::tuple_element<0,T>::type...>> 
    create_block_eigen_operator(const MatrixReplacement<1,1,T>&... operators) {
        typedef tuple_ns::tuple<typename tuple_ns::tuple_element<0,T>::type...> tuple_type;
        return MatrixReplacement<NI,NJ,tuple_type>(tuple_type(tuple_ns::get<0>(operators.m_blocks)...));
    }


#endif

}

#endif //OPERATORS_H_
