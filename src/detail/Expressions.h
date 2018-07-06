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

#ifndef EXPRESSIONS_DETAIL_H_
#define EXPRESSIONS_DETAIL_H_

#include "detail/Symbolic.h"

namespace Aboria {
namespace detail {

///////////////////
/// Expressions ///
///////////////////

template <typename Expr>
struct SymbolicExpr : proto::extends<Expr, SymbolicExpr<Expr>, SymbolicDomain> {

  explicit SymbolicExpr(Expr const &expr)
      : proto::extends<Expr, SymbolicExpr<Expr>, SymbolicDomain>(expr) {}
};

template <typename Expr, typename Enable = void> struct symbolic_helper {};

template <typename Expr>
struct symbolic_helper<
    Expr, typename boost::enable_if<is_const<typename proto::result_of::as_expr<
              Expr, detail::SymbolicDomain>::type>>::type> {

  typedef typename proto::result_of::as_expr<Expr, detail::SymbolicDomain>::type
      expr_type;

  typedef typename std::remove_const<
      typename proto::result_of::deep_copy<expr_type>::type>::type
      deep_copy_type;

  typedef EvalCtx<> const_context_type;
  typedef
      typename proto::result_of::eval<expr_type, const_context_type const>::type
          result;

  typedef typename std::remove_cv<
      typename std::remove_reference<result>::type>::type result_base_type;
};

template <typename Expr>
struct symbolic_helper<
    Expr,
    typename boost::enable_if<is_univariate<typename proto::result_of::as_expr<
        Expr, detail::SymbolicDomain>::type>>::type> {

  typedef typename proto::result_of::as_expr<Expr, detail::SymbolicDomain>::type
      expr_type;

  typedef typename std::remove_const<
      typename proto::result_of::deep_copy<expr_type>::type>::type
      deep_copy_type;

  typedef typename fusion::result_of::at_c<
      typename std::result_of<get_labels(expr_type, fusion::nil)>::type,
      0>::type label_a_type_ref;
  typedef typename std::remove_const<
      typename std::remove_reference<label_a_type_ref>::type>::type
      label_a_type;
  typedef typename label_a_type::particles_type particles_a_type;
  typedef typename particles_a_type::double_d double_d;
  typedef typename particles_a_type::const_reference particle_a_reference;
  typedef EvalCtx<fusion::map<fusion::pair<label_a_type, particle_a_reference>>>
      univariate_context_type;
  typedef typename proto::result_of::eval<
      expr_type, univariate_context_type const>::type result;

  typedef typename std::remove_cv<
      typename std::remove_reference<result>::type>::type result_base_type;
};

template <typename Expr>
struct symbolic_helper<
    Expr, typename boost::enable_if<
              is_univariate_with_no_label<typename proto::result_of::as_expr<
                  Expr, detail::SymbolicDomain>::type>>::type> {

  typedef typename proto::result_of::as_expr<Expr, detail::SymbolicDomain>::type
      expr_type;

  typedef typename std::remove_const<
      typename proto::result_of::deep_copy<expr_type>::type>::type
      deep_copy_type;

  struct dummy_type {};

  typedef dummy_type label_a_type;

  template <typename ParticleReference>
  using univariate_context_type =
      EvalCtx<fusion::map<fusion::pair<label_a_type, ParticleReference>>>;

  template <typename ParticleReference>
  using result = typename proto::result_of::eval<
      expr_type, univariate_context_type<ParticleReference> const>::type;

  template <typename ParticleReference>
  using result_base_type = typename std::remove_cv<
      typename std::remove_reference<result<ParticleReference>>::type>::type;
};

template <typename Expr>
struct symbolic_helper<
    Expr,
    typename boost::enable_if<is_bivariate<typename proto::result_of::as_expr<
        Expr, detail::SymbolicDomain>::type>>::type> {
  typedef typename proto::result_of::as_expr<Expr, detail::SymbolicDomain>::type
      expr_type;

  typedef typename std::remove_const<
      typename proto::result_of::deep_copy<expr_type>::type>::type
      deep_copy_type;

  typedef typename fusion::result_of::at_c<
      typename std::result_of<get_labels(expr_type, fusion::nil)>::type,
      0>::type label_first_type_ref;
  typedef typename fusion::result_of::at_c<
      typename std::result_of<get_labels(expr_type, fusion::nil)>::type,
      1>::type label_second_type_ref;
  typedef typename std::remove_const<
      typename std::remove_reference<label_first_type_ref>::type>::type
      label_first_type;
  typedef typename std::remove_const<
      typename std::remove_reference<label_second_type_ref>::type>::type
      label_second_type;

  static_assert(mpl::not_equal_to<typename label_first_type::depth,
                                  typename label_second_type::depth>::value,
                "label a depth equal to label b");

  typedef
      typename mpl::if_<mpl::less<typename label_first_type::depth,
                                  typename label_second_type::depth>,
                        label_first_type, label_second_type>::type label_a_type;
  typedef
      typename mpl::if_<mpl::greater<typename label_first_type::depth,
                                     typename label_second_type::depth>,
                        label_first_type, label_second_type>::type label_b_type;

  typedef typename label_a_type::particles_type particles_a_type;
  typedef typename label_b_type::particles_type particles_b_type;
  typedef typename particles_a_type::const_reference particle_a_reference;
  typedef typename particles_b_type::const_reference particle_b_reference;
  typedef typename particles_a_type::double_d double_d;
  typedef EvalCtx<fusion::map<fusion::pair<label_a_type, particle_a_reference>,
                              fusion::pair<label_b_type, particle_b_reference>>,
                  fusion::list<const double_d &>>
      bivariate_context_type;
  typedef typename proto::result_of::eval<
      expr_type, bivariate_context_type const>::type result;

  typedef typename std::remove_cv<
      typename std::remove_reference<result>::type>::type result_base_type;
};

/*
template<typename Expr>
struct SymbolicExpr<Expr,typename boost::enable_if<proto::matches<Expr,
detail::univariate_expr> >::type> : proto::extends<Expr, SymbolicExpr<Expr>,
SymbolicDomain> {

    typedef typename std::result_of<detail::univariate_expr(Expr)>::type
label_a_type_ref; typedef typename std::remove_reference<label_a_type_ref>::type
label_a_type; typedef typename label_a_type::particles_type particles_a_type;

    explicit SymbolicExpr(Expr const &expr)
        : proto::extends<Expr, SymbolicExpr<Expr>, SymbolicDomain>(expr)
    {}

    template <typename unknown_tuple_type=std::tuple<>>
    typename proto::result_of::eval<Expr const,
ParticleCtx<particles_a_type,unknown_tuple_type> const>::type eval( typename
particles_a_type::const_reference particle, unknown_tuple_type unknown_tuple =
std::tuple<>()) const { ParticleCtx<particles_a_type,unknown_tuple_type> const
ctx(particle,unknown_tuple); return proto::eval(*this, ctx);
    }


};

template<typename Expr>
struct SymbolicExpr<Expr,typename boost::enable_if<proto::matches<Expr,
detail::bivariate_expr> >::type> : proto::extends<Expr, SymbolicExpr<Expr>,
SymbolicDomain> {

    typedef typename std::result_of<detail::bivariate_expr(Expr)>::type::first
label_a_type_ref; typedef typename
std::result_of<detail::bivariate_expr(Expr)>::type::second label_b_type_ref;
    typedef typename std::remove_reference<label_a_type_ref>::type label_a_type;
    typedef typename std::remove_reference<label_b_type_ref>::type label_b_type;
    typedef typename label_a_type::particles_type particles_a_type;
    typedef typename label_b_type::particles_type particles_b_type;
    typedef typename particles_a_type::double_d double_d;

    explicit SymbolicExpr(Expr const &expr)
        : proto::extends<Expr, SymbolicExpr<Expr>, SymbolicDomain>(expr)
    {}

    template <typename unknown_tuple_type=std::tuple<>>
    typename proto::result_of::eval<Expr const,
TwoParticleCtx<particles_a_type,particles_b_type,unknown_tuple_type>
const>::type eval( const double_d& dx, typename
particles_a_type::const_reference particle1, typename
particles_b_type::const_reference particle2,unknown_tuple_type unknown_tuple1 =
std::tuple<>(), unknown_tuple_type unknown_tuple2 = std::tuple<>()) const {
        TwoParticleCtx<particles_a_type,particles_b_type, unknown_tuple_type>
const ctx(dx,particle1,particle2,unknown_tuple1,unknown_tuple2); return
proto::eval(*this, ctx);
    }

};
*/

template <typename Expr>
struct GeometryExpr : proto::extends<Expr, GeometryExpr<Expr>, GeometryDomain> {
  explicit GeometryExpr(Expr const &expr)
      : proto::extends<Expr, GeometryExpr<Expr>, GeometryDomain>(expr) {}
};

template <typename Expr>
struct LabelExpr : proto::extends<Expr, LabelExpr<Expr>, LabelDomain> {
  explicit LabelExpr(Expr const &expr)
      : proto::extends<Expr, LabelExpr<Expr>, LabelDomain>(expr) {}

  BOOST_PROTO_EXTENDS_USING_ASSIGN(LabelExpr)
};

template <unsigned int I, typename P>
struct Label
    : LabelExpr<typename proto::terminal<label<mpl::int_<I>, P>>::type> {

  typedef typename proto::terminal<label<mpl::int_<I>, P>>::type expr_type;
  typedef label<mpl::int_<I>, P> data_type;

  explicit Label(P &p) : LabelExpr<expr_type>(expr_type::make(data_type(p))) {}

  // BOOST_PROTO_EXTENDS_USING_ASSIGN(Label)
};

} // namespace detail
} // namespace Aboria
#endif
