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

#ifndef SYMBOLIC_ASSIGNMENT_DETAIL_H_
#define SYMBOLIC_ASSIGNMENT_DETAIL_H_

#include "detail/Evaluate.h"
#include "detail/Symbolic.h"

namespace Aboria {
namespace detail {

#define SUBSCRIPT_TYPE                                                         \
  boost::proto::exprns_::basic_expr<                                           \
      boost::proto::tagns_::tag::subscript,                                    \
      boost::proto::argsns_::list2<                                            \
          Aboria::detail::SymbolicExpr<boost::proto::exprns_::expr<            \
              boost::proto::tagns_::tag::terminal,                             \
              boost::proto::argsns_::term<                                     \
                  Aboria::detail::symbolic<VariableType>>,                     \
              0l>> &,                                                          \
          Aboria::Label<0u, ParticlesType> &>,                                 \
      2l>

// TODO: this seems a messy way to define a symbol subscripted by a label. might
// be better to put a subscript operator in the Symbol class?
// for: cleaner
// against: messes up the grammer (already existing subscript expression).
// Could have Symbol *not* be an expression, but then can't assume labels
// within expressions anymore....
template <typename VariableType, typename ParticlesType>
struct SymbolicExpr<SUBSCRIPT_TYPE>
    : proto::extends<SUBSCRIPT_TYPE, SymbolicExpr<SUBSCRIPT_TYPE>,
                     SymbolicDomain> {
  typedef SUBSCRIPT_TYPE Expr;

  typedef typename ParticlesType::position position;

  typedef typename std::remove_const<
      typename std::remove_reference<typename proto::result_of::value<
          typename proto::result_of::child_c<Expr, 1>::type>::type>::type>::type
      label_type;
  typedef typename proto::result_of::value<
      typename proto::result_of::child_c<Expr, 0>::type>::type symbol_type;

#undef SUBSCRIPT_TYPE

  explicit SymbolicExpr(Expr const &expr)
      : proto::extends<Expr, SymbolicExpr<Expr>, SymbolicDomain>(expr),
        msymbol(proto::value(proto::child_c<0>(expr))),
        mlabel(proto::value(proto::child_c<1>(expr))) {}

#define DEFINE_THE_OP(functor, the_op)                                         \
  template <typename ExprRHS>                                                  \
  const SymbolicExpr &operator the_op(ExprRHS const &expr) const {             \
    BOOST_MPL_ASSERT_NOT((boost::is_same<VariableType, id>));                  \
    evaluate_nonlinear<VariableType, functor>(                                 \
        proto::as_expr<SymbolicDomain>(expr), mlabel);                         \
    return *this;                                                              \
  }

  DEFINE_THE_OP(std::plus<void>, +=)
  DEFINE_THE_OP(std::minus<void>, -=)
  DEFINE_THE_OP(std::divides<void>, /=)
  DEFINE_THE_OP(std::multiplies<void>, *=)
  DEFINE_THE_OP(detail::return_second, =)

private:
  symbol_type &msymbol;
  label_type &mlabel;
};

} // namespace detail
} // namespace Aboria
#endif
