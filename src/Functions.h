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

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include "Symbolic.h"
#include <math.h>

namespace mpl = boost::mpl;
namespace proto = boost::proto;

namespace Aboria {

/// a macro to generate a symbolic function taking a single argument
/// using a given functor
/// \param function_name the name of the symbolic function to generate
/// \param function_to_wrap a functor that defines a typedef `result_type`
/// with the returned type, and a operator() that takes a single argument
/// and returns `result_type`
/// \param domain the domain of the symbolic function (normally SymbolicDomain)
#define ABORIA_UNARY_FUNCTION(function_name, function_to_wrap, domain)         \
  template <typename Expr>                                                     \
  typename proto::result_of::make_expr<proto::tag::function, domain,           \
                                       function_to_wrap,                       \
                                       Expr const &>::type const               \
  function_name(Expr const &arg) {                                             \
    return proto::make_expr<proto::tag::function, domain>(function_to_wrap(),  \
                                                          boost::ref(arg));    \
  }

/// a macro to generate a symbolic function taking two arguments
/// using a given functor
/// \param function_name the name of the symbolic function to generate
/// \param function_to_wrap a functor that defines a typedef `result_type`
/// with the returned type, and a operator() that takes two arguments
/// and returns `result_type`
/// \param domain the domain of the symbolic function (normally SymbolicDomain)
#define ABORIA_BINARY_FUNCTION(function_name, function_to_wrap, domain)        \
  template <typename Expr1, typename Expr2>                                    \
  typename proto::result_of::make_expr<proto::tag::function, domain,           \
                                       function_to_wrap, Expr1 const &,        \
                                       Expr2 const &>::type const              \
  function_name(Expr1 const &arg1, Expr2 const &arg2) {                        \
    return proto::make_expr<proto::tag::function, domain>(                     \
        function_to_wrap(), boost::ref(arg1), boost::ref(arg2));               \
  }

/// a macro to generate a symbolic function taking three arguments
/// using a given functor
/// \param function_name the name of the symbolic function to generate
/// \param function_to_wrap a functor that defines a typedef `result_type`
/// with the returned type, and a operator() that takes three arguments
/// and returns `result_type`
/// \param domain the domain of the symbolic function (normally SymbolicDomain)
#define ABORIA_TERNARY_FUNCTION(function_name, function_to_wrap, domain)       \
  template <typename Expr1, typename Expr2, typename Expr3>                    \
  typename proto::result_of::make_expr<                                        \
      proto::tag::function, domain, function_to_wrap, Expr1 const &,           \
      Expr2 const &, Expr3 const &>::type const                                \
  function_name(Expr1 const &arg1, Expr2 const &arg2, Expr3 const &arg3) {     \
    return proto::make_expr<proto::tag::function, domain>(                     \
        function_to_wrap(), boost::ref(arg1), boost::ref(arg2),                \
        boost::ref(arg3));                                                     \
  }

#define ABORIA_TAGGED_FUNCTION(name, tag)                                      \
  template <typename LABEL, typename CONDITIONAL, typename ARG, typename INIT> \
  typename proto::result_of::make_expr<tag, SymbolicDomain, LABEL const &,     \
                                       CONDITIONAL const &, ARG const &,       \
                                       INIT const &>::type const               \
  name(LABEL const &label, CONDITIONAL const &conditional, ARG const &arg,     \
       INIT const &init) {                                                     \
    return proto::make_expr<tag, SymbolicDomain>(                              \
        boost::ref(label), boost::ref(conditional), boost::ref(arg),           \
        boost::ref(init));                                                     \
  }

struct norm_fun {
  typedef double result_type;

  template <typename T, unsigned int N>
  result_type operator()(const Vector<T, N> &vector) const {
    return vector.norm();
  }
};

/// a symbolic norm function for Vect3d
///
ABORIA_UNARY_FUNCTION(norm, norm_fun, SymbolicDomain);

struct inf_norm_fun {
  typedef double result_type;

  template <typename T, unsigned int N>
  result_type operator()(const Vector<T, N> &vector) const {
    return vector.inf_norm();
  }
};

/// a symbolic norm function for Vect3d
///
ABORIA_UNARY_FUNCTION(inf_norm, inf_norm_fun, SymbolicDomain);

struct dot_fun {
  typedef double result_type;

  template <typename T, unsigned int N>
  result_type operator()(const Vector<T, N> &vector1,
                         const Vector<T, N> &vector2) const {
    return vector1.dot(vector2);
  }
};

/// a symbolic dot-product function for Vect3d
///
ABORIA_BINARY_FUNCTION(dot, dot_fun, SymbolicDomain);

struct exp_fun {
  typedef double result_type;

  result_type operator()(const double &arg) const { return std::exp(arg); }
};

/// a symbolic exponential function for scalars
///
ABORIA_UNARY_FUNCTION(exp, exp_fun, SymbolicDomain);

struct sqrt_fun {
  typedef double result_type;

  result_type operator()(const double &arg) const { return std::sqrt(arg); }
};

/// a symbolic square root function for scalars
///
ABORIA_UNARY_FUNCTION(sqrt, sqrt_fun, SymbolicDomain);

struct sign_fun {
  typedef double result_type;

  result_type operator()(const double &arg) const {
    return (0.0 < arg) - (arg < 0.0);
  }
};

/// a symbolic sign function for scalars
///
ABORIA_UNARY_FUNCTION(sign, sign_fun, SymbolicDomain);

struct erf_fun {
  typedef double result_type;

  result_type operator()(const double &arg) const { return std::erf(arg); }
};

/// a symbolic error function for scalars
///
ABORIA_UNARY_FUNCTION(erf, erf_fun, SymbolicDomain);

struct erfc_fun {
  typedef double result_type;

  result_type operator()(const double &arg) const { return std::erfc(arg); }
};

/// a symbolic complimentary error function for scalars
///
ABORIA_UNARY_FUNCTION(erfc, erfc_fun, SymbolicDomain);

struct log_fun {
  typedef double result_type;

  result_type operator()(const double &arg) const { return std::log(arg); }
};

/// a symbolic log function for scalars
///
ABORIA_UNARY_FUNCTION(log, log_fun, SymbolicDomain);

struct abs_fun {
  typedef double result_type;

  result_type operator()(const double &arg) const { return std::abs(arg); }
};

/// a symbolic absolute value function for scalars
///
ABORIA_UNARY_FUNCTION(abs, abs_fun, SymbolicDomain);

struct pow_fun {
  typedef double result_type;

  double operator()(const double d, const double exp) const {
    return std::pow(d, exp);
  }

  double operator()(const double d, const int exp) const {
    return std::pow(d, exp);
  }
};

/// a symbolic power function for that takes two arguements arg1
/// and arg2 and returns arg1 to the power of arg2
ABORIA_BINARY_FUNCTION(pow, pow_fun, SymbolicDomain);

template <typename T> struct reflect_fun {
  typedef Vector<double, 3> Vect3d;
  typedef Vect3d result_type;

  result_type operator()(const Vect3d &from, const Vect3d &to,
                         const T geometry) const {
    Vect3d result = to;
    reflect_once(from, result, geometry);
    return result;
  }
};

ABORIA_TERNARY_FUNCTION(reflect_, reflect_fun<Expr3>, SymbolicDomain);

} // namespace Aboria
#endif /* FUNCTIONS_H_ */
