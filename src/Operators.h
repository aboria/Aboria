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
    int sum(T1 s, T... ts) {
        return s + sum(ts...);
    }

    }

    // Inheriting EigenBase should not be needed in the future.
    template <unsigned int NI, unsigned int NJ, typename Blocks> 
    class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement<NI,NJ,Blocks>> {

        public:
        
    
            // Expose some compile-time information to Eigen:
            typedef double Scalar;
            typedef double RealScalar;
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
                return detail::sum(std::get<0>(std::get<I*NJ>(m_blocks)).size()...);
            }

            template<std::size_t... J>
            Index cols_impl(detail::index_sequence<J...>) const {
                return detail::sum(std::get<1>(std::get<J>(m_blocks)).size()...);
            }

            Index rows() const { 
                //std::cout << "rows = " << rows_impl(detail::make_index_sequence<NI>()) << std::endl;
                return rows_impl(detail::make_index_sequence<NI>());
            }

            Index cols() const { 
                //std::cout << "cols = " << cols_impl(detail::make_index_sequence<NJ>())<< std::endl;
                return cols_impl(detail::make_index_sequence<NJ>());
            }

            Index innerSize() const { 
                return rows();
            }

            Index outerSize() const { 
                return cols();
            }

            template <int I>
            Index start_col() const { 
                return cols_impl(detail::make_index_sequence<I>());
            }

            template <int I>
            Index size_col() const { 
                return std::get<1>(std::get<I>(m_blocks)).size();
            }

            template <int I>
            Index start_row() const { 
                return rows_impl(detail::make_index_sequence<I>());
            }

            template <int I>
            Index size_row() const { 
                return std::get<0>(std::get<I*NJ>(m_blocks)).size();
            }

            void resize(Index a_rows, Index a_cols) {
                // This method should not be needed in the future.
                assert(a_rows==0 && a_cols==0 || a_rows==rows() && a_cols==cols());
            }

            template < 
              typename particles_a_type, typename particles_b_type,
              typename expr_type, typename if_expr_type>
            Scalar coeff_impl_block(const Index i, const Index j, const std::tuple<const particles_a_type&,const particles_b_type&,expr_type,if_expr_type>& block) const {
                
                typedef typename particles_b_type::double_d double_d;
                typedef typename particles_b_type::position position;

                const particles_a_type& a = std::get<0>(block);
                const particles_b_type& b = std::get<1>(block);
                typename particles_a_type::const_reference ai = a[i];
                typename particles_b_type::const_reference bj = b[j];
                const double_d dx = a.correct_dx_for_periodicity(get<position>(bj)-get<position>(ai));

                //TODO: Have to copy the expressions (I think), since proto returns a non-const reference
                //to the stored constants, and if I want this to be const
                expr_type expr = std::get<2>(block);
                if_expr_type if_expr = std::get<3>(block);
                 
                if (eval(if_expr,dx,ai,bj)) {
                    return eval(expr,dx,ai,bj);
                } else {
                    return 0;
                }

            }



            template < 
              typename particles_b_type,
              typename expr_type, typename if_expr_type>
            Scalar coeff_impl_block(const Index i, const Index j, const std::tuple<const SizeOne&,const particles_b_type&,expr_type,if_expr_type>& block) const {
                
                typedef typename particles_b_type::double_d double_d;
                typedef typename particles_b_type::position position;

                const particles_b_type& b = std::get<1>(block);
                typename particles_b_type::const_reference bj = b[j];

                //TODO: Have to copy the expressions (I think), since proto returns a non-const reference
                //to the stored constants, and if I want this to be const
                expr_type expr = std::get<2>(block);
                if_expr_type if_expr = std::get<3>(block);
                 
                if (eval(if_expr,bj)) {
                    return eval(expr,bj);
                } else {
                    return 0;
                }

            }

            template < 
              typename particles_a_type,
              typename expr_type, typename if_expr_type>
            Scalar coeff_impl_block(const Index i, const Index j, const std::tuple<const particles_a_type&,const SizeOne&,expr_type,if_expr_type>& block) const {
                
                typedef typename particles_a_type::double_d double_d;
                typedef typename particles_a_type::position position;

                const particles_a_type& a = std::get<0>(block);
                typename particles_a_type::const_reference ai = a[i];

                //TODO: Have to copy the expressions (I think), since proto returns a non-const reference
                //to the stored constants, and if I want this to be const
                expr_type expr = std::get<2>(block);
                if_expr_type if_expr = std::get<3>(block);
                 
                if (eval(if_expr,ai)) {
                    return eval(expr,ai);
                } else {
                    return 0;
                }

            }

            template < 
              typename expr_type, typename if_expr_type>
            Scalar coeff_impl_block(const Index i, const Index j, const std::tuple<const SizeOne&,const SizeOne&,expr_type,if_expr_type>& block) const {
                //TODO: Have to copy the expressions (I think), since proto returns a non-const reference
                //to the stored constants, and if I want this to be const
                expr_type expr = std::get<2>(block);
                if_expr_type if_expr = std::get<3>(block);
                 
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
                            (coeff_impl_block(i,j,std::get<I*NJ>(m_blocks))):
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
    void evalTo_block(Eigen::VectorBlock<Dest> y, const Eigen::VectorBlock<Source>& rhs, const std::tuple<const particles_a_type&,const particles_b_type&,expr_type,if_expr_type>& block) {
        typedef typename particles_b_type::double_d double_d;
        typedef typename particles_b_type::position position;
        const particles_a_type& a = std::get<0>(block);
        const particles_b_type& b = std::get<1>(block);

        //TODO: Have to copy the expressions (I think), since proto returns a non-const reference
        //to the stored constants, and the eval product expression in eigen assumes that your matrix
        //replacement is a constant
        expr_type expr = std::get<2>(block);
        if_expr_type if_expr = std::get<3>(block);

        const size_t na = a.size();
        const size_t nb = b.size();

        if ((detail::is_const<expr_type>::value && (std::abs(eval(expr,double_d(),typename particles_a_type::value_type(),typename particles_b_type::value_type()))<=std::numeric_limits<double>::epsilon()))
                ||
                (detail::is_const<if_expr_type>::value && (!eval(if_expr,double_d(),typename particles_a_type::value_type(),typename particles_b_type::value_type())))
           ) {
            //std::cout << "zero a x b block" <<std::endl;
            return;
        }

        if (detail::is_const<if_expr_type>::value && eval(if_expr,double_d(),typename particles_a_type::value_type(),typename particles_b_type::value_type())==true) {
            //std::cout << "dense a x b block" <<std::endl;
            for (size_t i=0; i<na; ++i) {
                typename particles_a_type::const_reference ai = a[i];
                double sum = 0;
                for (size_t j=0; j<nb; ++j) {
                    typename particles_b_type::const_reference bj = b[j];
                    const double_d dx = a.correct_dx_for_periodicity(get<position>(bj)-get<position>(ai));
                    //std::cout << "a = "<<get<position>(ai)<<" b = "<<get<position>(bj)<<std::endl;
                    //std::cout << "using dx = "<<dx<<" rhs(j) = "<<rhs(j)<<" eval = "<<eval(expr,dx,ai,bj)<<std::endl;
                    sum += eval(expr,dx,ai,bj)*rhs(j);
                }
                y(i) += sum;
            }
        } else {
            //std::cout << "sparse a x b block" <<std::endl;
            for (size_t i=0; i<na; ++i) {
                typename particles_a_type::const_reference ai = a[i];
                double sum = 0;
                //std::cout << "evaluating fucntion for particle at "<<get<position>(ai)<<std::endl;
                for (auto pairj: b.get_neighbours(get<position>(ai))) {
                    const double_d & dx = std::get<1>(pairj);
                    typename particles_b_type::const_reference bj = std::get<0>(pairj);
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
    void evalTo_block(Eigen::VectorBlock<Dest> y, const Eigen::VectorBlock<Source>& rhs, const std::tuple<const SizeOne&,const particles_b_type&,expr_type,if_expr_type>& block) {

        typedef typename particles_b_type::double_d double_d;
        typedef typename particles_b_type::position position;

        const particles_b_type& b = std::get<1>(block);
        expr_type expr = std::get<2>(block);
        if_expr_type if_expr = std::get<3>(block);

        const size_t nb = b.size();

        if ((detail::is_const<expr_type>::value && (std::abs(eval(expr,typename particles_b_type::value_type()))<=std::numeric_limits<double>::epsilon())) ||
                (detail::is_const<if_expr_type>::value && (!eval(if_expr,typename particles_b_type::value_type())))) {
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
    void evalTo_block(Eigen::VectorBlock<Dest> y, const Eigen::VectorBlock<Source>& rhs, const std::tuple<const particles_a_type&,const SizeOne&,expr_type,if_expr_type>& block) {

        typedef typename particles_a_type::double_d double_d;
        typedef typename particles_a_type::position position;

        const particles_a_type& a = std::get<0>(block);
        expr_type expr = std::get<2>(block);
        if_expr_type if_expr = std::get<3>(block);

        const size_t na = a.size();

        if ((detail::is_const<expr_type>::value && (std::abs(eval(expr,typename particles_a_type::value_type()))<=std::numeric_limits<double>::epsilon())) ||
            (detail::is_const<if_expr_type>::value && (!eval(if_expr,typename particles_a_type::value_type())))) {
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
    void evalTo_block(Eigen::VectorBlock<Dest> y, const Eigen::VectorBlock<Source>& rhs, const std::tuple<const SizeOne&,const SizeOne&,expr_type,if_expr_type>& block) {

        expr_type expr = std::get<2>(block);
        if_expr_type if_expr = std::get<3>(block);

        //TODO: should use return type instead of double
        if ((detail::is_const<expr_type>::value && (std::abs(eval(expr))<=std::numeric_limits<double>::epsilon())) ||

            (detail::is_const<if_expr_type>::value && (!eval(if_expr)))) {
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
    void evalTo_unpack_blocks(Dest& y, const MatrixReplacement<NI,NJ,Blocks>& lhs, const Rhs& rhs, const std::tuple<I,J,T1>& block, const T&... other_blocks) {
        evalTo_block(y.segment(lhs.template start_row<I::value>(),lhs.template size_row<I::value>()),
                rhs.segment(lhs.template start_col<J::value>(),lhs.template size_col<J::value>()),
                std::get<2>(block));
        evalTo_unpack_blocks(y,lhs,rhs,other_blocks...);
    }

    template<typename Dest, unsigned int NI, unsigned int NJ, typename Blocks, typename Rhs, std::size_t... I>
    void evalTo_impl(Dest& y, const MatrixReplacement<NI,NJ,Blocks>& lhs, const Rhs& rhs, detail::index_sequence<I...>) {
        evalTo_unpack_blocks(y,lhs,rhs,std::make_tuple(mpl::int_<I/NJ>(), mpl::int_<I%NJ>(), std::get<I>(lhs.m_blocks))...);
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
    static void scaleAndAddTo(Dest& y, const Aboria::MatrixReplacement<NI,NJ,Blocks>& lhs, const Rhs& rhs, const Scalar& alpha) {
        // This method should implement "y += alpha * lhs * rhs" inplace,
        // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
        assert(alpha==Scalar(1) && "scaling is not implemented");
        evalTo_impl(y, lhs, rhs, Aboria::detail::make_index_sequence<NI*NJ>());
    }
};
}
}

namespace Aboria {    
    template <typename A, unsigned int A_depth, 
              typename B, unsigned int B_depth, 
              typename Expr, typename IfExpr=detail::SymbolicExpr<typename proto::terminal<bool>::type>>
    MatrixReplacement<1,1,std::tuple<std::tuple<const A&,const B&,
        typename detail::symbolic_helper<Expr>::deep_copy_type,
        typename detail::symbolic_helper<IfExpr>::deep_copy_type>>>
    create_eigen_operator(
            const Label<A_depth,A>& a, 
            const Label<B_depth,B>& b, 
            const Expr& expr, 
            const IfExpr& if_expr = IfExpr(proto::terminal<bool>::type::make(true))) 
    {
        typedef std::tuple<const A&, const B&, 
            typename detail::symbolic_helper<Expr>::deep_copy_type,
            typename detail::symbolic_helper<IfExpr>::deep_copy_type> block_type;
        return MatrixReplacement<1,1,std::tuple<block_type>>(
                std::make_tuple(
                    block_type(proto::value(a).get_particles(),proto::value(b).get_particles(),
                        deep_copy(expr),
                        deep_copy(if_expr))
                ));
    }

    template <typename B, unsigned int B_depth, 
              typename Expr, typename IfExpr=detail::SymbolicExpr<typename proto::terminal<bool>::type>>
    MatrixReplacement<1,1,std::tuple<std::tuple<const SizeOne&,const B&,
        typename detail::symbolic_helper<Expr>::deep_copy_type,
        typename detail::symbolic_helper<IfExpr>::deep_copy_type>>>
    create_eigen_operator(
            const One& a, 
            const Label<B_depth,B>& b, 
            const Expr& expr, 
            const IfExpr& if_expr = IfExpr(proto::terminal<bool>::type::make(true))) 
    {
        typedef std::tuple<const SizeOne&, const B&, 
            typename detail::symbolic_helper<Expr>::deep_copy_type,
            typename detail::symbolic_helper<IfExpr>::deep_copy_type> block_type;
        return MatrixReplacement<1,1,std::tuple<block_type>>(
                std::make_tuple(
                    block_type(a.get_size_one(),proto::value(b).get_particles(),
                        deep_copy(expr),
                        deep_copy(if_expr))
                ));
    }

    template <typename A, unsigned int A_depth, 
              typename Expr, typename IfExpr=detail::SymbolicExpr<typename proto::terminal<bool>::type>>
    MatrixReplacement<1,1,std::tuple<std::tuple<const A&,const SizeOne&,
        typename detail::symbolic_helper<Expr>::deep_copy_type,
        typename detail::symbolic_helper<IfExpr>::deep_copy_type>>>
    create_eigen_operator(
            const Label<A_depth,A>& a, 
            const One& b, 
            const Expr& expr, 
            const IfExpr& if_expr = IfExpr(proto::terminal<bool>::type::make(true))) 
    {
        typedef std::tuple<const A&, const SizeOne&, 
            typename detail::symbolic_helper<Expr>::deep_copy_type,
            typename detail::symbolic_helper<IfExpr>::deep_copy_type> block_type;
        return MatrixReplacement<1,1,std::tuple<block_type>>(
                std::make_tuple(
                    block_type(proto::value(a).get_particles(),b.get_size_one(),
                        deep_copy(expr),
                        deep_copy(if_expr))
                ));
    }

    template <typename Expr, typename IfExpr=detail::SymbolicExpr<typename proto::terminal<bool>::type>>
    MatrixReplacement<1,1,std::tuple<std::tuple<const SizeOne&,const SizeOne&,
        typename detail::symbolic_helper<Expr>::deep_copy_type,
        typename detail::symbolic_helper<IfExpr>::deep_copy_type>>>
    create_eigen_operator(
            const One& a, 
            const One& b, 
            const Expr& expr, 
            const IfExpr& if_expr = IfExpr(proto::terminal<bool>::type::make(true))) 
    {
        typedef std::tuple<const SizeOne&, const SizeOne&, 
            typename detail::symbolic_helper<Expr>::deep_copy_type,
            typename detail::symbolic_helper<IfExpr>::deep_copy_type> block_type;
        return MatrixReplacement<1,1,std::tuple<block_type>>(
                std::make_tuple(
                    block_type(a.get_size_one(),b.get_size_one(),
                        deep_copy(expr),
                        deep_copy(if_expr))
                ));
    }

    template <unsigned int NI, unsigned NJ, typename ... T>
    MatrixReplacement<NI,NJ,std::tuple<T...>> create_block_eigen_operator(const MatrixReplacement<1,1,std::tuple<T>>&... operators) {
        return MatrixReplacement<NI,NJ,std::tuple<T...>>(std::make_tuple(std::get<0>(operators.m_blocks)...));
    }


#endif

}

#endif //OPERATORS_H_
