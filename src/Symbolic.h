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

/// returns a symbolic expression with the terminals held by value, rather than
/// the normal held by reference.
template <typename Expr>
typename detail::symbolic_helper<Expr>::deep_copy_type
deep_copy(Expr const &expr) {
  return proto::deep_copy(proto::as_expr<detail::SymbolicDomain>(expr));
}

/// evaluate a given expression that returns a constant value (scaler or vector)
/// \params expr the expression to evaluate. Must be an expression that returns
/// a constant, i.e. that does not depend on a particle's variables \return the
/// result of the expression
// TODO: move univariate and bivariate down here to
template <class Expr>
typename boost::enable_if<
    detail::is_const<Expr>,
    typename detail::symbolic_helper<Expr>::result_base_type const>::type
// typename detail::symbolic_helper<Expr>::result>::type
eval(Expr &expr) {
  typename detail::symbolic_helper<Expr>::const_context_type const ctx;
  return proto::eval(expr, ctx);
}

/// evaluate a given expression that depends on a single label $i$,
/// for a single input particle
/// \params expr the expression to evaluate
/// \params particle_a the particle to be substituted for the label $i$
/// \return the result of the expression after substituting in the particle
/// values
template <typename Expr>
typename boost::enable_if<detail::is_univariate<Expr>,
                          typename detail::symbolic_helper<Expr>::result>::type
eval(Expr &expr,
     const typename detail::symbolic_helper<Expr>::particle_a_reference
         &particle_a) {
  typedef
      typename detail::symbolic_helper<Expr>::univariate_context_type ctx_type;
  typedef typename detail::symbolic_helper<Expr>::label_a_type label_type;
  ctx_type const ctx(fusion::make_map<label_type>(particle_a));
  return proto::eval(expr, ctx);
}

/// evaluate a given expression that depends on a single label $i$,
/// for a single input particle
/// \params expr the expression to evaluate
/// \params particle_a the particle to be substituted for the label $i$
/// \return the result of the expression after substituting in the particle
/// values
template <typename Expr, typename ParticleReference>
typename boost::enable_if<detail::is_univariate_with_no_label<Expr>,
                          typename detail::symbolic_helper<
                              Expr>::template result<ParticleReference>>::type
eval(Expr &expr, const ParticleReference &particle_a) {
  typedef typename detail::symbolic_helper<
      Expr>::template univariate_context_type<ParticleReference>
      ctx_type;
  typedef typename detail::symbolic_helper<Expr>::label_a_type label_type;
  ctx_type const ctx(fusion::make_map<label_type>(particle_a));
  return proto::eval(expr, ctx);
}

/// evaluate a given expression that returns a constant value (scaler or vector)
/// \params expr the expression to evaluate. Must be an expression that returns
/// a constant, i.e. that does not depend on a particle's variables \params
/// particle_a dummy arguement, not used \return the result of the expression
template <typename Expr, typename AnyRef>
typename boost::enable_if<
    detail::is_const<Expr>,
    typename detail::symbolic_helper<Expr>::result_base_type const>::type
// why do I need this?, should be
// typename detail::symbolic_helper<Expr>::result>::type
eval(Expr &expr, const AnyRef &particle_a) {
  typename detail::symbolic_helper<Expr>::const_context_type const ctx;
  return proto::eval(expr, ctx);
}

/// evaluate a given expression that depends on two labels $i$ and $j$
/// \params expr the expression to evaluate
/// \params dx the shortest vector from particle_a to particle_b
/// \params particle_a the particle to be substituted for the label $i$
/// \params particle_b the particle to be substituted for the label $j$
/// \return the result of the expression after substituting in the particle
/// values
template <typename Expr>
typename boost::enable_if<detail::is_bivariate<Expr>,
                          typename detail::symbolic_helper<Expr>::result>::type
eval(Expr &expr, const typename detail::symbolic_helper<Expr>::double_d &dx,
     const typename detail::symbolic_helper<Expr>::particle_a_reference
         &particle_a,
     const typename detail::symbolic_helper<Expr>::particle_b_reference
         &particle_b) {

  typedef
      typename detail::symbolic_helper<Expr>::bivariate_context_type ctx_type;
  typedef typename detail::symbolic_helper<Expr>::label_a_type label_a_type;
  typedef typename detail::symbolic_helper<Expr>::label_b_type label_b_type;
  ctx_type const ctx(
      fusion::make_map<label_a_type, label_b_type>(particle_a, particle_b),
      fusion::make_list(boost::cref(dx)));
  return proto::eval(expr, ctx);
}

/// evaluate a given expression that depends on a single label $i$,
/// for a single input particle
/// \params expr the expression to evaluate
/// \params dx dummy arguement not used
/// \params particle_a the particle to be substituted for the label $i$
/// \params particle_b dummy particle, not used
/// \return the result of the expression after substituting in the particle
/// values
template <typename Expr, typename AnyRef>
typename boost::enable_if<
    mpl::and_<detail::is_univariate<Expr>,
              mpl::equal_to<
                  typename detail::symbolic_helper<Expr>::label_a_type::depth,
                  mpl::int_<0>>>,
    typename detail::symbolic_helper<Expr>::result>::type
eval(Expr &expr, const typename detail::symbolic_helper<Expr>::double_d &dx,
     const typename detail::symbolic_helper<Expr>::particle_a_reference
         &particle_a,
     const AnyRef &particle_b) {

  typedef
      typename detail::symbolic_helper<Expr>::univariate_context_type ctx_type;
  typedef typename detail::symbolic_helper<Expr>::label_a_type label_type;

  ctx_type const ctx(fusion::make_map<label_type>(particle_a));
  return proto::eval(expr, ctx);
}

/// evaluate a given expression that depends on a single label $i$,
/// for a single input particle
/// \params expr the expression to evaluate
/// \params dx dummy arguement not used
/// \params particle_a dummy particle, not used
/// \params particle_b the particle to be substituted for the label $i$
/// \return the result of the expression after substituting in the particle
/// values
template <typename Expr, typename AnyRef>
typename boost::enable_if<
    mpl::and_<detail::is_univariate<Expr>,
              mpl::equal_to<
                  typename detail::symbolic_helper<Expr>::label_a_type::depth,
                  mpl::int_<1>>>,
    typename detail::symbolic_helper<Expr>::result>::type
eval(Expr &expr, const typename detail::symbolic_helper<Expr>::double_d &dx,
     const AnyRef &particle_a,
     const typename detail::symbolic_helper<Expr>::particle_a_reference
         &particle_b) {

  typedef
      typename detail::symbolic_helper<Expr>::univariate_context_type ctx_type;
  typedef typename detail::symbolic_helper<Expr>::label_a_type label_type;

  ctx_type const ctx(fusion::make_map<label_type>(particle_b));
  return proto::eval(expr, ctx);
}

/// evaluate a given expression that returns a constant value (scalar or vector)
/// \params expr the expression to evaluate. Must be an expression that returns
/// a constant, i.e. that does not depend on a particle's variables \params dx
/// dummy argument, not used \params particle_a dummy argument, not used \params
/// particle_b dummy argument, not used \return the result of the expression
template <typename Expr, typename AnyDx, typename AnyRef1, typename AnyRef2>
typename boost::enable_if<
    detail::is_const<Expr>,
    typename detail::symbolic_helper<Expr>::result_base_type const>::type
// why do I need this?, should be
// typename detail::symbolic_helper<Expr>::result>::type
eval(Expr &expr, const AnyDx &dx, const AnyRef1 &particle_a,
     const AnyRef2 &particle_b) {
  typename detail::symbolic_helper<Expr>::const_context_type const ctx;
  return proto::eval(expr, ctx);
}

/// returns `true` if the expression always evaluates to less than or equal
/// to `std::numeric_limits<double>::epsilon()`
template <class Expr>
typename boost::enable_if<detail::is_const<Expr>, bool>::type
is_trivially_zero(Expr &expr) {
  // TODO: should use return type instead of double
  return std::abs(eval(expr)) <= std::numeric_limits<double>::epsilon();
}

/// returns `false` for non constant expressions (might be not zero)
template <class Expr>
typename boost::enable_if<mpl::not_<detail::is_const<Expr>>, bool>::type
is_trivially_zero(Expr &expr) {
  return false;
}

/// returns `true` if the expression always evaluates to `false`
template <class IfExpr>
typename boost::enable_if<detail::is_const<IfExpr>, bool>::type
is_trivially_false(IfExpr &expr) {
  return eval(expr) == false;
}

/// returns `false` for non constant expressions (might be not `false`)
template <class IfExpr>
typename boost::enable_if<mpl::not_<detail::is_const<IfExpr>>, bool>::type
is_trivially_false(IfExpr &expr) {
  return false;
}

/// returns `true` if the expression always evaluates to `true`
template <class IfExpr>
typename boost::enable_if<detail::is_const<IfExpr>, bool>::type
is_trivially_true(IfExpr &expr) {
  return eval(expr) == true;
}

/// returns `false` for non constant expressions (might be not `true`)
template <class IfExpr>
typename boost::enable_if<mpl::not_<detail::is_const<IfExpr>>, bool>::type
is_trivially_true(IfExpr &expr) {
  return false;
}

/// define a symbol from a Variable type T. This symbol can then be used in
/// expressions \param T the Variable type of the symbol
template <typename T>
struct Symbol : detail::SymbolicExpr<
                    typename proto::terminal<detail::symbolic<T>>::type> {

  /// type of the (terminal) expression used to express Symbol
  typedef typename proto::terminal<detail::symbolic<T>>::type expr_type;
  /// the underlying Variable type
  typedef T variable_type;
  /// value_type of the underlying Variable type
  typedef typename T::value_type value_type;
  /// type of internal data class used to store Symbol information (e.g.
  /// buffering)
  typedef detail::symbolic<T> data_type;

  /// constructor
  explicit Symbol()
      : detail::SymbolicExpr<expr_type>(expr_type::make(data_type())) {}
};

/*
template<unsigned int I>
struct Unknown
    : detail::SymbolicExpr<typename proto::terminal<detail::unknown<mpl::int_<I>
> >::type> {

        typedef typename proto::terminal<detail::unknown<mpl::int_<I>> >::type
expr_type; typedef detail::unknown<mpl::int_<I>> data_type;

        /// constructor
        explicit Unknown()
            : detail::SymbolicExpr<expr_type>( expr_type::make(data_type()) )
        {}
};
*/

/// define a label with a given depth \p I that referres to a particle
/// container with type \p P
/// \param I the label depth. If an expression contains a number of
/// nested sums over particles, this depth specifies which sum the label
/// belongs to
/// \param P the type of the particle container used to store the
/// particle and variable data
template <unsigned int I, typename P>
struct Label
    : detail::LabelExpr<
          typename proto::terminal<detail::label<mpl::int_<I>, P>>::type> {

  /// type of the (terminal) expression used to express Label
  typedef
      typename proto::terminal<detail::label<mpl::int_<I>, P>>::type expr_type;
  /// type of internal data class used to store Label information
  typedef detail::label<mpl::int_<I>, P> data_type;

  /// constructor
  explicit Label(P &p)
      : detail::LabelExpr<expr_type>(expr_type::make(data_type(p))) {}

  template <typename Variable> void resize_buffer(const size_t n) {
    get<Variable>(proto::value(*this).get_buffers()).resize(n);
  }
  // BOOST_PROTO_EXTENDS_USING_ASSIGN(Label)
};

template <unsigned int I, typename P> Label<I, P> create_label(P &p) {
  return Label<I, P>(p);
}

/// a symbolic class used to refer to the difference between neighbouring
/// particles position vectors. Note that for periodic domains this might
/// be different than `get<position)(a)-get<position>(b)`, and in this case
/// always gives the *shortest* position difference
template <typename L1, typename L2>
struct Dx
    : detail::SymbolicExpr<typename proto::terminal<
          detail::dx<typename L1::data_type, typename L2::data_type>>::type> {

  typedef typename detail::dx<typename L1::data_type, typename L2::data_type>
      data_type;
  typedef typename proto::terminal<data_type>::type expr_type;

  /// constructor
  explicit Dx(L1 &la, L2 &lb)
      : detail::SymbolicExpr<expr_type>(
            expr_type::make(data_type(proto::value(la), proto::value(lb)))) {}
};

/// helper function to create a Dx class that refers to the shortest distance
/// between two particles given by labels la and lb
/// \see Dx
/// \see Label
template <typename L1, typename L2> Dx<L1, L2> create_dx(L1 &la, L2 &lb) {
  return Dx<L1, L2>(la, lb);
}

/// a symbolic class used to return a normally distributed random variable. This
/// uses the random number generator of the current particle to generate the
/// random variable
struct Normal : proto::terminal<detail::normal>::type {};

/// a symbolic class used to return a uniformly distributed random variable.
/// This uses the random number generator of the current particle to generate
/// the random variable
struct Uniform : proto::terminal<detail::uniform>::type {};

/// a symbolic class that, when evaluated, returns a Vect3d class.
template <typename T, unsigned int N>
struct VectorSymbolic
    : detail::SymbolicExpr<
          typename proto::terminal<detail::vector<T, N>>::type> {

  typedef typename proto::terminal<detail::vector<T, N>>::type expr_type;
  typedef detail::vector<T, N> data_type;

  explicit VectorSymbolic()
      : detail::SymbolicExpr<expr_type>(expr_type::make(data_type())) {}
};

/// an accumulation expression, for example a sum or product over neighbouring
/// particles
/// \param T a functor class that performs the accumulation, for example
/// `std::plus<double>`
template <typename T>
struct Accumulate : detail::SymbolicExpr<
                        typename proto::terminal<detail::accumulate<T>>::type> {

  typedef typename proto::terminal<detail::accumulate<T>>::type expr_type;
  typedef detail::accumulate<T> data_type;

  /// empty constructor, makes an instantiation of the functor class \p T
  explicit Accumulate()
      : detail::SymbolicExpr<expr_type>(expr_type::make(data_type())) {}

  /// constructor that passes a single argument to the functor class \p T
  template <typename T2>
  explicit Accumulate(const T2 &arg)
      : detail::SymbolicExpr<expr_type>(expr_type::make(data_type(arg))) {}

  /// a function used to set the initial value of the accumulation
  /// \param the initial value. Must be the same type used by the
  /// accumulation functor \p T
  void set_init(const typename data_type::init_type &arg) {
    proto::value(*this).set_init(arg);
  }
};

/// an accumulation expression, for example a sum or product, over neighbouring
/// particles within a given radius
/// \param T a functor class that performs the accumulation, for example
/// `std::plus<double>`
template <typename T, int LNormNumber = 2,
          typename Terminal = typename detail::accumulate_within_distance<
              T, mpl::int_<LNormNumber>>>
struct AccumulateWithinDistance
    : detail::SymbolicExpr<typename proto::terminal<Terminal>::type> {

  typedef Terminal data_type;
  typedef typename proto::terminal<data_type>::type expr_type;

  /// empty constructor, makes an instantiation of the functor class \p T
  explicit AccumulateWithinDistance(const double max_distance = 1)
      : detail::SymbolicExpr<expr_type>(
            expr_type::make(data_type(max_distance))) {}

  /// constructor that passes a single argument to the functor class \p T
  template <typename T2>
  explicit AccumulateWithinDistance(const double max_distance, const T2 &arg)
      : detail::SymbolicExpr<expr_type>(
            expr_type::make(data_type(max_distance, arg))) {}

  /// a function used to set the initial value of the accumulation
  /// \param the initial value. Must be the same type used by the
  /// accumulation functor \p T
  void set_init(const typename data_type::init_type &arg) {
    proto::value(*this).set_init(arg);
  }

  void set_max_distance(const double max_distance) {
    proto::value(*this).set_max_distance(max_distance);
  }
};

/// convenient functor to get a minumum value using the Accumulate expression
/// \code
///     Accumulate<max<double>> max;
/// \endcode
template <typename T> struct min {
  typedef T result_type;
  T operator()(const T arg1, const T arg2) const {
    // std::cout << "finding min of "<<arg1<<" and "<<arg2<<std::endl;
    return std::min(arg1, arg2);
  }
};

/// convenient functor to get a maximum value using the Accumulate expression
/// \code
///     Accumulate<min<double>> min;
/// \endcode
template <typename T> struct max {
  typedef T result_type;
  T operator()(const T arg1, const T arg2) const {
    return std::max(arg1, arg2);
  }
};

/*
/// a symbolic class that refers to a Geometry class.
template <typename T>
struct GeometrySymbolic
    : detail::SymbolicExpr<typename proto::terminal<detail::geometry<T> >::type>
{

        typedef typename proto::terminal<detail::geometry<T> >::type expr_type;
        typedef detail::geometry<T> data_type;

        explicit GeometrySymbolic()
            : detail::SymbolicExpr<expr_type>( expr_type::make(data_type()) )
        {}

};

/// a symbolic class that refers to a Geometries class.
template <typename T>
struct GeometriesSymbolic
    : detail::SymbolicExpr<typename proto::terminal<detail::geometries<T>
>::type> {

        typedef typename proto::terminal<detail::geometries<T> >::type
expr_type; typedef detail::geometries<T> data_type;

        explicit GeometriesSymbolic()
            : detail::SymbolicExpr<expr_type>( expr_type::make(data_type()) )
        {}
};
*/

} // namespace Aboria

#include "detail/Contexts.h"

#endif /* SYMBOLIC_H_ */
