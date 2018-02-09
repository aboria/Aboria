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

#ifndef EVALUATE_DETAIL_H_
#define EVALUATE_DETAIL_H_

#include "Symbolic.h"

namespace Aboria {
namespace detail {

//
// Helper classes
//

struct return_second {
  template <typename T1, typename T2>
  const T2 &operator()(const T1 &ignored, const T2 &to_return) {
    return to_return;
  }
};

template <typename LabelType>
struct is_not_my_label
    : proto::and_<
          proto::terminal<label<_, _>>,
          proto::if_<mpl::not_<boost::is_same<proto::_value, LabelType>>()>> {};

template <typename VariableType>
struct is_my_symbol : proto::terminal<symbolic<VariableType>> {};

template <typename VariableType, typename LabelType>
struct is_not_aliased
    : proto::or_<
          proto::and_<proto::terminal<proto::_>,
                      proto::not_<proto::and_<
                          proto::if_<boost::is_same<
                              VariableType,
                              typename LabelType::particles_type::position>()>,
                          proto::terminal<dx<_, _>>>>,
                      proto::not_<proto::and_<
                          proto::if_<boost::is_same<
                              VariableType,
                              typename LabelType::particles_type::position>()>,
                          proto::terminal<accumulate_within_distance<_, _>>>>>,
          proto::and_<
              proto::nary_expr<proto::_, proto::vararg<is_not_aliased<
                                             VariableType, LabelType>>>,
              proto::not_<proto::subscript<is_my_symbol<VariableType>,
                                           is_not_my_label<LabelType>>>>> {};

// expose alias checking for testing in metafunctions.h
template <typename SymbolType, typename LabelType, typename ExprRHS>
typename boost::enable_if<
    mpl::not_<proto::matches<typename proto::result_of::as_expr<ExprRHS>::type,
                             is_not_aliased<typename SymbolType::variable_type,
                                            typename LabelType::data_type>>>,
    mpl::true_>::type
alias_check(SymbolType const &, LabelType const &, ExprRHS const &) {
  return mpl::true_();
}

// expose alias checking for testing in metafunctions.h
template <typename SymbolType, typename LabelType, typename ExprRHS>
typename boost::enable_if<
    proto::matches<typename proto::result_of::as_expr<ExprRHS>::type,
                   is_not_aliased<typename SymbolType::variable_type,
                                  typename LabelType::data_type>>,
    mpl::false_>::type
alias_check(SymbolType const &, LabelType const &, ExprRHS const &) {
  return mpl::false_();
}

template <typename LabelType, typename ExprRHS>
typename boost::enable_if<detail::is_univariate<ExprRHS>, void>::type
check_valid_assign_expr(const LabelType &label, ExprRHS const &expr) {
  typedef typename LabelType::particles_type particles_type;
  typedef typename detail::result_of::get_labels<ExprRHS>::type rhs_labels_type;
  typedef typename fusion::result_of::at_c<rhs_labels_type, 0>::type
      rhs_label_type_ref;
  typedef typename std::remove_cv<
      typename std::remove_reference<rhs_label_type_ref>::type>::type
      rhs_label_type;

  const rhs_label_type &rhs_label =
      fusion::at_c<0>(detail::get_labels()(expr, fusion::nil()));
  static_assert(std::is_same<rhs_label_type, LabelType>::value,
                "Labels on LHS and RHS of assign expression do not match");
  particles_type &particles = label.get_particles();
  const auto &particles_in_expr = rhs_label.get_particles();
  CHECK(&particles_in_expr == &particles,
        "labels on LHS and RHS of assign expression do not refer to the same "
        "particles container");
}

template <typename ParticlesType, typename ExprRHS>
typename boost::enable_if<
    mpl::or_<detail::is_const<ExprRHS>,
             detail::is_univariate_with_no_label<ExprRHS>>,
    void>::type
check_valid_assign_expr(const ParticlesType &particles, ExprRHS const &expr) {}

template <typename ParticlesType, typename ExprRHS>
typename boost::enable_if<detail::is_bivariate<ExprRHS>, void>::type
check_valid_assign_expr(const ParticlesType &particles, ExprRHS const &expr) {
  static_assert(!detail::is_bivariate<ExprRHS>::value,
                "asignment expression must be constant or univariate");
}

} // namespace detail
} // namespace Aboria
#endif
