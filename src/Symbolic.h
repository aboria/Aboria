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


#ifndef SYMBOLIC_H_
#define SYMBOLIC_H_

#include "Vector.h"
#include "detail/Symbolic.h"

namespace Aboria {

    using detail::SymbolicDomain;

    /// evaluate a given expression that returns a constant value (scaler or vector)
    /// \params expr the expression to evaluate. Must be an expression that returns a constant, i.e. that does not depend on a particle's variables
    /// \return the result of the expression
    //TODO: move univariate and bivariate down here to
    //TODO: check its a constant expression...
    template <class Expr>
    typename proto::result_of::eval<Expr const,  detail::ConstantCtx const>::type
    eval(Expr const &expr) {
            detail::ConstantCtx const ctx;
            return proto::eval(expr,ctx);
    }

    template<typename Expr, typename Unknown=std::tuple<>>  
    typename detail::symbolic_helper<Expr>::template result<Unknown> 
    eval(Expr const &expr, 
            const typename detail::symbolic_helper<Expr>::particle_a_reference& particle_a, 
            const Unknown& unknown_a=std::tuple<>()) {
        typename detail::symbolic_helper<Expr>::template univariate_context_type<Unknown> const ctx(particle_a,unknown_a);
        return proto::eval(expr, ctx);
    }

    template<typename Expr, typename Unknown=std::tuple<>>  
    typename detail::symbolic_helper<Expr>::template result<Unknown> 
    eval(Expr const &expr, 
            const typename detail::symbolic_helper<Expr>::double_d& dx,
            const typename detail::symbolic_helper<Expr>::particle_a_reference& particle_a, 
            const typename detail::symbolic_helper<Expr>::particle_b_reference& particle_b, 
            const Unknown& unknown_a=std::tuple<>(), 
            const Unknown& unknown_b=std::tuple<>()) {
        typename detail::symbolic_helper<Expr>::template bivariate_context_type<Unknown> const ctx(dx,particle_a,particle_b,unknown_a,unknown_b);
        return proto::eval(expr, ctx);
    }

    template<typename ParticleType, typename unknown_tuple_type, typename Expr,  
        typename=typename boost::enable_if<proto::matches<Expr, detail::univariate_expr> >::type>
    typename proto::result_of::eval<Expr const, detail::ParticleCtx<ParticleType,unknown_tuple_type> const>::type
    eval_bivariate(Expr const &expr, 
            typename ParticleType::const_reference particle, 
            const unknown_tuple_type& unknown_tuple) {
        detail::ParticleCtx<ParticleType,unknown_tuple_type> const ctx(particle,unknown_tuple);
        return proto::eval(expr, ctx);
    }

    template<typename ParticleType1, typename ParticleType2, typename Expr,  
        typename=typename boost::enable_if<proto::matches<Expr, detail::bivariate_expr>>::type>
    typename proto::result_of::eval<Expr const, detail::TwoParticleCtx<ParticleType1,ParticleType2,std::tuple<>> const>::type
    eval(Expr const &expr, const typename ParticleType1::double_d& dx, 
            typename ParticleType1::const_reference particle1, 
            typename ParticleType2::const_reference particle2) {
        detail::TwoParticleCtx<ParticleType1,ParticleType2,std::tuple<>> const ctx(dx, particle1, particle2, std::tuple<>(), std::tuple<>());
        return proto::eval(expr, ctx);
    }


    template<typename ParticleType1, typename ParticleType2, typename unknown_tuple_type, typename Expr,  
        typename=typename boost::enable_if<proto::matches<Expr, detail::bivariate_expr>>::type>
    typename proto::result_of::eval<Expr const, detail::TwoParticleCtx<ParticleType1,ParticleType2,unknown_tuple_type> const>::type
    eval(Expr const &expr, const typename ParticleType1::double_d& dx, 
            typename ParticleType1::const_reference particle1, 
            typename ParticleType2::const_reference particle2, 
            const unknown_tuple_type& unknown_tuple1,
            const unknown_tuple_type& unknown_tuple2) {
        detail::TwoParticleCtx<ParticleType1,ParticleType2,unknown_tuple_type> const ctx(dx, particle1, particle2, unknown_tuple1, unknown_tuple2);
        return proto::eval(expr, ctx);
    }


    /// define a symbol from a Variable type T. This symbol can then be used in expressions
    /// \param T the Variable type of the symbol
    template<typename T>
    struct Symbol
        : detail::SymbolicExpr<typename proto::terminal<detail::symbolic<T> >::type> {

            /// type of the (terminal) expression used to express Symbol
            typedef typename proto::terminal<detail::symbolic<T> >::type expr_type;
            /// value_type of the underlying Variable type
            typedef typename T::value_type value_type;
            /// type of internal data class used to store Symbol information (e.g. buffering)
            typedef detail::symbolic<T> data_type;

            /// constructor
            explicit Symbol()
                : detail::SymbolicExpr<expr_type>( expr_type::make(data_type()) )
            {}
    };

    template<unsigned int I>
    struct Unknown 
        : detail::SymbolicExpr<typename proto::terminal<detail::unknown<mpl::int_<I> > >::type> {

            typedef typename proto::terminal<detail::unknown<mpl::int_<I>> >::type expr_type;
            typedef detail::unknown<mpl::int_<I>> data_type;

            /// constructor
            explicit Unknown()
                : detail::SymbolicExpr<expr_type>( expr_type::make(data_type()) )
            {}
    };




    /// define a label with a given depth \p I that referres to a particle 
    /// container with type \p P
    /// \param I the label depth. If an expression contains a number of 
    /// nested sums over particles, this depth specifies which sum the label 
    /// belongs to
    /// \param P the type of the particle container used to store the 
    /// particle and variable data
    template<unsigned int I, typename P>
    struct Label 
        : detail::LabelExpr<typename proto::terminal<detail::label<mpl::int_<I>,P > >::type> {

            /// type of the (terminal) expression used to express Label
            typedef typename proto::terminal<detail::label<mpl::int_<I>,P > >::type expr_type;
            /// type of internal data class used to store Label information 
            typedef detail::label<mpl::int_<I>,P > data_type;

            /// constructor
            explicit Label(P& p)
                : detail::LabelExpr<expr_type>( expr_type::make(data_type(p)))
            {}


            //BOOST_PROTO_EXTENDS_USING_ASSIGN(Label)
    };


    template<unsigned int I, typename P>
    Label<I,P> create_label(P& p) {
        return Label<I,P>(p);
    }

    /// a symbolic class used to refer to the difference between neighbouring 
    /// particles position vectors. Note that for periodic domains this might 
    /// be different than `get<position)(a)-get<position>(b)`, and in this case 
    /// always gives the *shortest* position difference
    template<typename L1, typename L2>
    struct Dx 
        : detail::SymbolicExpr<typename proto::terminal<detail::dx<typename L1::data_type,typename L2::data_type> >::type> {

            typedef typename detail::dx<typename L1::data_type,typename L2::data_type> data_type;
            typedef typename proto::terminal<data_type>::type expr_type;

            /// constructor
            explicit Dx(L1& la, L2& lb)
                : detail::SymbolicExpr<expr_type>( expr_type::make(data_type(proto::value(la),proto::value(lb))) )
            {}
    };

    
    template<typename L1, typename L2>
    Dx<L1,L2> create_dx(L1& la, L2& lb) {
        return Dx<L1,L2>(la,lb);
    }

    /// a symbolic class used to return a normally distributed random variable. This uses
    /// the random number generator of the current particle to generate the random
    /// variable
    struct Normal
        : proto::terminal<detail::normal>::type {};

    
    /// a symbolic class that, when evaluated, returns a Vect3d class. 
    template <typename T>
    struct VectorSymbolic
        : detail::SymbolicExpr<typename proto::terminal<detail::vector<T> >::type> {

            typedef typename proto::terminal<detail::vector<T> >::type expr_type;
            typedef detail::vector<T> data_type;

            explicit VectorSymbolic()
                : detail::SymbolicExpr<expr_type>( expr_type::make(data_type()) )
            {}

    };

    /// an accumulation expression, for example a sum or product over neighbouring 
    /// particles
    /// \param T a functor class that performs the accumulation, for example `std::plus<double>`  
    template <typename T>
    struct Accumulate
        : detail::SymbolicExpr<typename proto::terminal<detail::accumulate<T> >::type> {

            typedef typename proto::terminal<detail::accumulate<T> >::type expr_type;
            typedef detail::accumulate<T> data_type;

            /// empty constructor, makes an instantiation of the functor class \p T
            explicit Accumulate()
                : detail::SymbolicExpr<expr_type>( expr_type::make(data_type()) )
            {}

            /// constructor that passes a single argument to the functor class \p T
            template<typename T2>
                explicit Accumulate(const T2& arg)
                : detail::SymbolicExpr<expr_type>( expr_type::make(data_type(arg)) )
            {}

            /// a function used to set the initial value of the accumulation
            /// \param the initial value. Must be the same type used by the 
            /// accumulation functor \p T
            void set_init(const typename data_type::init_type & arg) {
                proto::value(*this).set_init(arg);
            }

    };

    /// convenient functor to get a minumum value using the Accumulate expression
    template <typename T>
    struct min {
        typedef T result_type;
        T operator()(const T arg1, const T arg2) const {
            return std::min(arg1,arg2);
        }
    };

    /// convenient functor to get a maximum value using the Accumulate expression
    template <typename T>
    struct max {
        typedef T result_type;
        T operator()(const T arg1, const T arg2) const {
            return std::max(arg1,arg2);
        }
    };

    /// a symbolic class that refers to a Geometry class. 
    template <typename T>
    struct GeometrySymbolic
        : detail::SymbolicExpr<typename proto::terminal<detail::geometry<T> >::type> {

            typedef typename proto::terminal<detail::geometry<T> >::type expr_type;
            typedef detail::geometry<T> data_type;

            explicit GeometrySymbolic()
                : detail::SymbolicExpr<expr_type>( expr_type::make(data_type()) )
            {}

    };

    /// a symbolic class that refers to a Geometries class. 
    template <typename T>
    struct GeometriesSymbolic
        : detail::SymbolicExpr<typename proto::terminal<detail::geometries<T> >::type> {

            typedef typename proto::terminal<detail::geometries<T> >::type expr_type;
            typedef detail::geometries<T> data_type;

            explicit GeometriesSymbolic()
                : detail::SymbolicExpr<expr_type>( expr_type::make(data_type()) )
            {}
    };


}



#endif /* SYMBOLIC_H_ */
