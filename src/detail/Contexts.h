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

#ifndef CONTEXTS_DETAIL_H_
#define CONTEXTS_DETAIL_H_

#include "detail/Symbolic.h"

namespace Aboria {
namespace detail {

////////////////
/// Contexts ///
////////////////

// Here is an evaluation context that indexes into a lazy vector
// expression, and combines the result.
template <typename labels_type, typename dx_type> struct EvalCtx {
  typedef typename fusion::result_of::size<labels_type>::type size_type;
  typedef typename fusion::result_of::size<dx_type>::type dx_size_type;
  static constexpr int dx_size = size_type::value * (size_type::value - 1) / 2;

  // BOOST_MPL_ASSERT_MSG(dx_size_type::value==dx_size,DX_SIZE_NOT_CONSISTENT_WITH_LABELS_SIZE,(dx_size,dx_size_type));
  static_assert(dx_size_type::value == dx_size,
                "dx size not consitent with labels_size");

  EvalCtx(labels_type labels = fusion::nil(), dx_type dx = fusion::nil())
      : m_labels(labels), m_dx(dx) {}

  template <typename Expr
            // defaulted template parameters, so we can
            // specialize on the expressions that need
            // special handling.
            ,
            typename Tag = typename proto::tag_of<Expr>::type,
            typename Enable = void>
  struct eval : proto::default_eval<Expr, EvalCtx const> {
    // terminals with variables, labels, uniform or normal should not be
    // evaluated...
    static_assert(
        mpl::not_<mpl::or_<
            proto::matches<Expr, proto::terminal<uniform>>,
            proto::matches<Expr, proto::terminal<normal>>,
            proto::matches<Expr, proto::terminal<label<_, _>>>,
            proto::matches<Expr, proto::terminal<symbolic<_>>>>>::value,
        "\nError: labels, symbols, uniform or normals can not be evaluated on "
        "their own.\n\tAlways subscript a symbol with a label, e.g. p[a] or "
        "symbol[label]");
  };

  // Handle normal and uniform subscripts here...
  template <typename Expr>
  struct eval<
      Expr, proto::tag::subscript,
      typename boost::enable_if<mpl::and_<
          proto::matches<typename proto::result_of::child_c<Expr, 1>::type,
                         proto::terminal<label<_, _>>>,
          mpl::or_<
              proto::matches<typename proto::result_of::child_c<Expr, 0>::type,
                             proto::terminal<normal>>,
              proto::matches<typename proto::result_of::child_c<Expr, 0>::type,
                             proto::terminal<uniform>>>,
          mpl::greater<size_type, mpl::int_<0>>>>::type> {

    typedef typename proto::result_of::child_c<Expr, 1>::type child1_type;
    typedef typename proto::result_of::value<child1_type>::type label_type;

    static_assert(fusion::result_of::has_key<labels_type, label_type>::value,
                  "label not in evaluation context");

    typedef double result_type;

    result_type operator()(Expr &expr, EvalCtx const &ctx) const {
      // Normal and uniform terminal types have a operator() that takes a
      // generator. Pass the random generator for the labeled particle to this
      // operator()
      return proto::value(proto::child_c<0>(expr))(
          // need to const_cast this cause everything is
          // normally held as a const &. Could cause problems???
          const_cast<generator_type &>(
              get<generator>(fusion::at_key<label_type>(ctx.m_labels))));
    }
  };

  // Handle other subscripts here...
  template <typename Expr>
  struct eval<
      Expr, proto::tag::subscript,
      typename boost::enable_if<mpl::and_<
          proto::matches<typename proto::result_of::child_c<Expr, 1>::type,
                         proto::terminal<label<_, _>>>,
          mpl::not_<mpl::or_<
              proto::matches<typename proto::result_of::child_c<Expr, 0>::type,
                             proto::terminal<normal>>,
              proto::matches<typename proto::result_of::child_c<Expr, 0>::type,
                             proto::terminal<uniform>>>>,
          mpl::greater<size_type, mpl::int_<0>>>>::type> {

    typedef typename proto::result_of::child_c<Expr, 1>::type child1_type;
    typedef typename proto::result_of::value<child1_type>::type label_type;

    /*
    typedef typename label_type::particles_type particles_type;
    typedef typename particles_type::const_reference particle_ref;
    typedef typename fusion::pair<label_type,particle_ref> search_type;

    */
    static_assert(fusion::result_of::has_key<labels_type, label_type>::value,
                  "label not in evaluation context");

    typedef typename proto::result_of::child_c<Expr, 0>::type child0_type;
    typedef typename proto::result_of::value<child0_type>::type symbolic_type;
    typedef typename symbolic_type::variable_type variable_type;
    typedef const typename variable_type::value_type &result_type;

    result_type operator()(Expr &expr, EvalCtx const &ctx) const {
      return get<variable_type>(fusion::at_key<label_type>(ctx.m_labels));
    }
  };

  // Handle dx terminals here...
  template <typename Expr>
  struct eval<Expr, proto::tag::terminal,
              typename boost::enable_if<
                  mpl::and_<proto::matches<Expr, proto::terminal<dx<_, _>>>,
                            mpl::equal<size_type, mpl::int_<2>>>>::type> {
    typedef typename fusion::result_of::front<const dx_type>::type result_type;
    typedef typename proto::result_of::value<Expr>::type expr_dx;
    typedef typename expr_dx::label_a_type expr_label_a_type;
    typedef typename expr_dx::label_b_type expr_label_b_type;

    /*
    BOOST_MPL_ASSERT_MSG((fusion::result_of::has_key<labels_type,expr_label_a_type>::value),ASDFASDFASDF,(expr_dx,labels_type,expr_label_a_type));
    BOOST_MPL_ASSERT_MSG((fusion::result_of::has_key<labels_type,expr_label_b_type>::value),ASDFASDFASDF,(expr_label_b_type));
    */

    static_assert(
        fusion::result_of::has_key<labels_type, expr_label_a_type>::value,
        "dx label a not in evaluation context");
    static_assert(
        fusion::result_of::has_key<labels_type, expr_label_b_type>::value,
        "dx label b not in evaluation context");

    result_type operator()(Expr &expr, EvalCtx const &ctx) const {
      return fusion::front(ctx.m_dx);
    }
  };

  template <typename result_type, typename label_type, typename expr_type,
            typename accumulate_type>
  static result_type dense_sum_impl(
      const label_type &label, expr_type &expr, accumulate_type &accum,
      const EvalCtx &ctx,
      mpl::int_<0>) { // note: using tag dispatching here cause I couldn't
                      // figure out how to do this via enable_if....

    result_type sum = accum.init;

    auto particles = label.get_particles();

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < particles.size(); ++i) {
      const auto &p = particles[i];
      auto new_labels = fusion::make_map<label_type>(p);
      EvalCtx<decltype(new_labels), decltype(ctx.m_dx)> const new_ctx(
          new_labels, ctx.m_dx);
      sum = accum.functor(sum, proto::eval(expr, new_ctx));
    }
    return sum;
  }

  template <typename result_type, typename label_b_type, typename expr_type,
            typename accumulate_type>
  static result_type dense_sum_impl(const label_b_type &label, expr_type &expr,
                                    accumulate_type &accum, const EvalCtx &ctx,
                                    mpl::int_<1>) {

    typedef typename label_b_type::particles_type particles_b_type;
    typedef typename particles_b_type::position position;
    typedef typename position::value_type double_d;
    typedef typename std::remove_reference<typename fusion::result_of::at_c<
        labels_type, 0>::type>::type::first_type label_a_type;
    typedef typename std::remove_reference<typename fusion::result_of::at_c<
        labels_type, 0>::type>::type::second_type const_a_reference;
    typedef typename particles_b_type::const_reference const_b_reference;
    typedef typename fusion::map<fusion::pair<label_a_type, const_a_reference>,
                                 fusion::pair<label_b_type, const_b_reference>>
        map_type;
    typedef fusion::list<const double_d &> list_type;

    const particles_b_type &particlesb = label.get_particles();
    ASSERT(!particlesb.get_periodic().any(),
           "periodic does not work with dense");
    const_a_reference ai = fusion::front(ctx.m_labels).second;
    const size_t nb = particlesb.size();

    result_type sum = accum.init;
    if (is_trivially_zero(expr)) {
      return sum;
    } else {
      for (size_t i = 0; i < nb; ++i) {
        const_b_reference bi = particlesb[i];

        EvalCtx<map_type, list_type> const new_ctx(
            fusion::make_map<label_a_type, label_b_type>(ai, bi),
            fusion::make_list(get<position>(bi) - get<position>(ai)));

        sum = accum.functor(sum, proto::eval(expr, new_ctx));
      }
    }
    return sum;
  }

  template <typename result_type, typename label_b_type, typename expr_type,
            typename accumulate_type, typename dummy = size_type>
  static result_type sparse_sum_impl(const label_b_type &label, expr_type &expr,
                                     accumulate_type &accum, const EvalCtx &ctx,
                                     mpl::int_<1>) {

    typedef typename label_b_type::particles_type particles_b_type;
    const particles_b_type &particlesb = label.get_particles();

    typedef typename std::remove_reference<typename fusion::result_of::at_c<
        labels_type, 0>::type>::type::first_type label_a_type;

    typedef typename std::remove_reference<typename fusion::result_of::at_c<
        labels_type, 0>::type>::type::second_type const_a_reference;

    typedef typename particles_b_type::position position;
    typedef typename position::value_type double_d;
    typedef typename particles_b_type::const_reference const_b_reference;

    const_a_reference ai = fusion::front(ctx.m_labels).second;

    typedef typename fusion::map<fusion::pair<label_a_type, const_a_reference>,
                                 fusion::pair<label_b_type, const_b_reference>>
        map_type;

    typedef fusion::list<const double_d &> list_type;
    const int LNormNumber = accumulate_type::norm_number_type::value;

    result_type sum = accum.init;
    // TODO: get query range and put it in box search
    for (auto b = distance_search<LNormNumber>(
             particlesb.get_query(), get<position>(ai), accum.max_distance);
         b != false; ++b) {
      EvalCtx<map_type, list_type> const new_ctx(
          // fusion::make_map<label_a_type, label_b_type>(ai, *b),
          map_type(ai, *b),
          list_type(b.dx())); // fusion::make_list(b.dx()));

      sum = accum.functor(sum, proto::eval(expr, new_ctx));
    }
    return sum;
  }

  template <typename Expr>
  struct eval<Expr, proto::tag::function,
              typename boost::enable_if<
                  proto::matches<Expr, AccumulateGrammar>>::type> {

    typedef typename proto::result_of::child_c<Expr, 0>::type child0_type;
    typedef typename proto::result_of::value<child0_type>::type
        functor_terminal_type;
    typedef typename functor_terminal_type::functor_type functor_type;
    typedef typename functor_type::result_type result_type;

    result_type operator()(Expr &expr, EvalCtx const &ctx) const {
      return dense_sum_impl<result_type>(
          proto::value(proto::child_c<1>(expr)), proto::child_c<2>(expr),
          proto::value(proto::child_c<0>(expr)), ctx, size_type());
    }
  };

  template <typename Expr>
  struct eval<Expr, proto::tag::function,
              typename boost::enable_if<mpl::and_<
                  proto::matches<Expr, AccumulateWithinDistanceGrammar>,
                  mpl::equal<size_type, mpl::int_<1>>>>::type> {

    typedef typename proto::result_of::child_c<Expr, 0>::type child0_type;
    typedef typename proto::result_of::value<child0_type>::type
        functor_terminal_type;
    typedef typename functor_terminal_type::functor_type functor_type;
    typedef typename functor_type::result_type result_type;

    result_type operator()(Expr &expr, EvalCtx const &ctx) const {
      return sparse_sum_impl<result_type>(
          proto::value(proto::child_c<1>(expr)), proto::child_c<2>(expr),
          proto::value(proto::child_c<0>(expr)), ctx, size_type());
    }
  };

  labels_type m_labels;
  dx_type m_dx;
};

} // namespace detail
} // namespace Aboria
#endif
