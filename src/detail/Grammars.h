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

#ifndef GRAMMARS_DETAIL_H_
#define GRAMMARS_DETAIL_H_

#include "detail/Symbolic.h"

namespace Aboria {
namespace detail {

////////////////
/// Grammars ///
////////////////

// A grammar which matches all the assignment operators,
// so we can easily disable them.
struct AssignOps : proto::switch_<struct AssignOpsCases> {};

// Here are the cases used by the switch_ above.
struct AssignOpsCases {
  template <typename Tag, int D = 0> struct case_ : proto::not_<_> {};

  template <int D> struct case_<proto::tag::plus_assign, D> : _ {};
  template <int D> struct case_<proto::tag::minus_assign, D> : _ {};
  template <int D> struct case_<proto::tag::multiplies_assign, D> : _ {};
  template <int D> struct case_<proto::tag::divides_assign, D> : _ {};
  template <int D> struct case_<proto::tag::modulus_assign, D> : _ {};
  template <int D> struct case_<proto::tag::shift_left_assign, D> : _ {};
  template <int D> struct case_<proto::tag::shift_right_assign, D> : _ {};
  template <int D> struct case_<proto::tag::bitwise_and_assign, D> : _ {};
  template <int D> struct case_<proto::tag::bitwise_or_assign, D> : _ {};
  template <int D> struct case_<proto::tag::bitwise_xor_assign, D> : _ {};
};

struct SymbolicGrammar
    : proto::or_<
          proto::terminal<_>,
          proto::and_<proto::nary_expr<_, proto::vararg<SymbolicGrammar>>,
                      proto::not_<AssignOps>>>

{};

struct get_particles : proto::callable {
  template <typename Sig> struct result;

  template <typename This, typename LabelType>
  struct result<This(const LabelType &)>
      : boost::add_reference<typename LabelType::particles_type> {};

  /*
  template<typename This, typename LabelType>
      struct result<This(LabelType)>
      : boost::add_reference<typename LabelType::particles_type>
      {};
      */

  template <typename LabelType>
  typename LabelType::particles_type &operator()(const LabelType &label) {
    return label.get_particles();
  }
};

struct LabelGrammar
    : proto::when<proto::terminal<label<_, _>>, get_particles(proto::_value)> {
};

struct SubscriptGrammar
    : proto::when<proto::terminal<label<_, _>>, proto::_value> {};

struct GeometryGrammar
    : proto::function<proto::terminal<geometry<_>>, SymbolicGrammar,
                      SymbolicGrammar, SymbolicGrammar> {};

struct GeometriesGrammar : proto::function<proto::terminal<geometries<_>>,
                                           LabelGrammar, SymbolicGrammar> {};

struct AccumulateGrammar : proto::function<proto::terminal<accumulate<_>>,
                                           LabelGrammar, SymbolicGrammar> {};

struct AccumulateWithinDistanceGrammar
    : proto::function<proto::terminal<accumulate_within_distance<_, _>>,
                      LabelGrammar, SymbolicGrammar> {};

/*
struct add_dummy_if_empty: proto::callable {
    template<typename Sig>
    struct result;

    template<typename This>
    struct result<This(const fusion::nil&)>:
        fusion::result_of::as_list<dummy_type>
    {};


    template<typename This, typename T, typename = typename
        std::enable_if<!std::is_same<T,fusion::nil>::value>::type>
    struct result<This(T)> {
        typedef T type;
    };

    template<typename This>
    typename result<add_dummy_if_empty(const fusion::nil& state)>::type
    operator()(const fusion::nil& state) {
        return fusion::make_list(dummy());
        //return typename result<add_dummy_if_empty(const T&)>::type();
    }

    template<typename This, typename T, typename = typename
        std::enable_if<!std::is_same<T,fusion::nil>::value>::type>
    typename result<add_dummy_if_empty(const T&)>::type
    operator()(const T& state) {
        return state;
        return
    }
};
*/

struct remove_label : proto::callable {
  template <typename Sig> struct result;

  template <typename This, typename T1, typename T2>
  struct result<This(T1, T2)>
      : fusion::result_of::as_list<typename fusion::result_of::remove_if<
            typename std::remove_reference<T2>::type,
            boost::is_same<
                // typename std::remove_reference<typename
                // std::remove_const<T1>::type>::type &,
                T1, mpl::_>>::type> {};

  /*
  template<typename This, typename T1, typename T2>
  struct result<This(T1, const T2&)>:
      fusion::result_of::as_list<
          typename
  fusion::result_of::remove_if<T2,boost::is_same<T1,mpl::_>>::type
      >
  {};
  */

  template <typename T1, typename T2>
  typename result<remove_label(const T1 &, const T2 &)>::type
  operator()(const T1 &label, const T2 &state) {
    return fusion::as_list(
        fusion::remove_if<boost::is_same<const T1 &, mpl::_>>(state));
  }
};

struct push_back_if_new : proto::callable {
  template <typename Sig> struct result;

  /*
  template<typename This, typename T1, typename T2>
  struct result<This(const T1&, const T2&)> {
      typedef fusion::cons<const T1&,
              typename boost::result_of<remove_label(const T1&,const
  T2&)>::type> type;
  };
  */

  template <typename This, typename T1, typename T2>
  struct result<This(T1 &, T2)> {
    typedef fusion::cons<const T1 &, typename boost::result_of<remove_label(
                                         const T1 &, const T2 &)>::type>
        type;
  };

  /*
  template<typename This, typename T1, typename T2>
  struct result<This(T1&, const T2&)> {
      typedef fusion::cons<const T1&,
              typename boost::result_of<remove_label(const T1&,const
  T2&)>::type> type;
  };
  */

  /*
  template<typename This, typename T1, typename T2>
  struct result<This(T1&, const T2&)>:
      fusion::result_of::as_list<
          typename fusion::result_of::push_back<
              typename boost::result_of<remove_label(T1&,const T2&)>::type,
  T1>::type
      >
  {};
  */

  template <typename T1, typename T2>
  typename result<push_back_if_new(const T1 &, const T2 &)>::type
  operator()(const T1 &label, const T2 &state) {
    return typename result<push_back_if_new(const T1 &, const T2 &)>::type(
        label, remove_label()(label, state));
  }
};

struct get_a_from_dx : proto::callable {
  template <typename Sig> struct result;

  template <typename This, typename DxType> struct result<This(DxType &)> {
    typedef const typename DxType::label_a_type &type;
  };

  template <typename DxType>
  const typename DxType::label_a_type &operator()(const DxType &dx) {
    return dx.get_label_a();
  }
};

struct get_b_from_dx : proto::callable {
  template <typename Sig> struct result;

  template <typename This, typename DxType> struct result<This(DxType &)> {
    typedef const typename DxType::label_b_type &type;
  };

  template <typename DxType>
  const typename DxType::label_b_type &operator()(const DxType &dx) {
    return dx.get_label_b();
  }
};

struct get_labels
    : proto::or_<

          // Put all label terminals at the head of the
          // list that we're building in the "state" parameter
          proto::when<proto::terminal<label<_, _>>,
                      push_back_if_new(proto::_value, proto::_state)>,
          proto::when<proto::terminal<dx<_, _>>,
                      push_back_if_new(
                          get_a_from_dx(proto::_value),
                          push_back_if_new(get_b_from_dx(proto::_value),
                                           proto::_state))>
          // don't add other terminals
          ,
          proto::when<proto::terminal<_>, proto::_state>,
          proto::when<proto::function<proto::terminal<accumulate<_>>, _, _>,
                      remove_label(proto::_value(proto::_child1),
                                   get_labels(proto::_child2, proto::_state))>,
          proto::when<
              proto::function<proto::terminal<accumulate_within_distance<_, _>>,
                              _, _>,
              remove_label(proto::_value(proto::_child1),
                           get_labels(proto::_child2, proto::_state))>,
          proto::otherwise<proto::fold<_, proto::_state, get_labels>>> {};

} // namespace detail
struct norm_fun;
struct inf_norm_fun;
struct dot_fun;
namespace detail {

struct norm_dx
    : proto::or_<proto::function<proto::terminal<Aboria::norm_fun>,
                                 proto::terminal<dx<_, _>>>,
                 proto::function<proto::terminal<Aboria::inf_norm_fun>,
                                 proto::terminal<dx<_, _>>>,
                 proto::function<proto::terminal<Aboria::dot_fun>,
                                 proto::terminal<dx<_, _>>,
                                 proto::terminal<dx<_, _>>>> {};

struct accumulate_within_distance_expr
    : proto::or_<
          proto::when<proto::terminal<_>, mpl::bool_<false>()>,
          proto::when<
              proto::function<proto::terminal<accumulate_within_distance<_, _>>,
                              _, _>,
              mpl::bool_<true>()>,
          proto::when<proto::function<proto::terminal<accumulate<_>>, _, _>,
                      mpl::bool_<false>()>,
          proto::when<proto::nary_expr<_, proto::vararg<_>>,
                      proto::fold<_, mpl::bool_<false>(),
                                  mpl::or_<accumulate_within_distance_expr,
                                           proto::_state>()>>

          > {};

namespace result_of {

template <typename Expr>
struct accumulate_within_distance_expr
    : boost::result_of<Aboria::detail::accumulate_within_distance_expr(Expr)> {
};

} // namespace result_of

struct range_if_expr
    : proto::or_<proto::when<proto::less<norm_dx, SymbolicGrammar>,
                             SymbolicGrammar(proto::_right)>,
                 proto::when<proto::less_equal<norm_dx, SymbolicGrammar>,
                             SymbolicGrammar(proto::_right)>,
                 proto::when<proto::greater<SymbolicGrammar, norm_dx>,
                             SymbolicGrammar(proto::_left)>,
                 proto::when<proto::greater_equal<SymbolicGrammar, norm_dx>,
                             SymbolicGrammar(proto::_left)>,
                 proto::when<proto::logical_and<range_if_expr, SymbolicGrammar>,
                             range_if_expr(proto::_left)>,
                 proto::when<proto::logical_and<SymbolicGrammar, range_if_expr>,
                             range_if_expr(proto::_right)>
                 /*
                     , proto::when<
                         proto::nary_expr<_, proto::vararg<range_if_expr> >
                         ,range_if_expr(proto::_right)
                         >
                         */
                 > {};

namespace result_of {
template <typename Expr> struct get_labels {
  typedef typename boost::result_of<Aboria::detail::get_labels(
      Expr, fusion::nil)>::type get_labels_result;
  typedef typename std::remove_reference<get_labels_result>::type type;
};

} // namespace result_of

template <typename Expr>
struct is_const
    : mpl::and_<
          mpl::equal_to<typename fusion::result_of::size<
                            typename result_of::get_labels<Expr>::type>::type,
                        mpl::int_<0>>,
          mpl::not_<typename result_of::accumulate_within_distance_expr<
              Expr>::type>> {};

template <typename Expr>
struct is_univariate
    : mpl::equal_to<typename fusion::result_of::size<
                        typename result_of::get_labels<Expr>::type>::type,
                    mpl::int_<1>> {};

template <typename Expr>
struct is_univariate_with_no_label
    : mpl::and_<
          mpl::equal_to<typename fusion::result_of::size<
                            typename result_of::get_labels<Expr>::type>::type,
                        mpl::int_<0>>,
          typename result_of::accumulate_within_distance_expr<Expr>::type> {};

template <typename Expr>
struct is_bivariate
    : mpl::equal_to<typename fusion::result_of::size<
                        typename result_of::get_labels<Expr>::type>::type,
                    mpl::int_<2>> {};

template <typename Expr, unsigned int I> struct get_label_c {
  typedef typename result_of::get_labels<Expr>::type labels_type;
  typedef typename fusion::result_of::value_at_c<get_labels, I>::type type;
};

} // namespace detail
} // namespace Aboria
#endif
