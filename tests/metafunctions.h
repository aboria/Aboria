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

#ifndef METAFUNCTIONSTEST_H_
#define METAFUNCTIONSTEST_H_

#include <cxxtest/TestSuite.h>

#include "Aboria.h"
#include <boost/mpl/assert.hpp>

using namespace Aboria;

class MetafunctionsTest : public CxxTest::TestSuite {
public:
  void test_aliasing(void) {
    ABORIA_VARIABLE(scalar, double, "scalar")
    typedef Particles<std::tuple<scalar>> ParticlesType;
    ParticlesType particles;
    Symbol<scalar> s;
    Symbol<position_d<3>> r;
    Label<0, ParticlesType> a(particles);
    Label<1, ParticlesType> b(particles);
    auto dx = create_dx(a, b);
    Accumulate<std::plus<double>> sum;
    AccumulateWithinDistance<std::plus<vdouble3>> sparse_sum;
    Accumulate<std::plus<vdouble3>> sumv;
    VectorSymbolic<double, 3> vector;

    static_assert(decltype(detail::alias_check(s, a, sum(b, s[b])))::value,
                  "should be aliased");

    static_assert(decltype(detail::alias_check(s, a, sum(b, 1 + s[b])))::value,
                  "should be aliased");

    static_assert(decltype(detail::alias_check(
                      r, a, sparse_sum(b, vector(s[b], 0, 0))))::value,
                  "should be aliased");

    static_assert(decltype(detail::alias_check(r, a, sumv(b, r[b])))::value,
                  "should be aliased");

    static_assert(decltype(detail::alias_check(
                      r, a, sumv(b, vdouble3(1, 2, 3) + dx)))::value,
                  "should be aliased");

    static_assert(!decltype(detail::alias_check(s, a, sum(b, s[a])))::value,
                  "should not be aliased");

    static_assert(!decltype(detail::alias_check(s, a, sum(b, norm(dx))))::value,
                  "should not be aliased");

    static_assert(!decltype(detail::alias_check(s, a, s[a]))::value,
                  "should not be aliased");

    static_assert(!decltype(detail::alias_check(s, a, s[a] + r[a][1]))::value,
                  "should not be aliased");
  }

  void test_expr_matching(void) {
    ABORIA_VARIABLE(scalar, double, "scalar")
    typedef Particles<std::tuple<scalar>> ParticlesType;
    ParticlesType particles;
    Symbol<scalar> s;
    Symbol<id> id_;
    Label<0, ParticlesType> a(particles);
    Label<1, ParticlesType> b(particles);
    auto dx = create_dx(a, b);

    static_assert(
        proto::matches<decltype(norm(dx) < 3), detail::range_if_expr>::value,
        "norm(dx)<3 does not match range_if_expr");

    static_assert(proto::matches<decltype(inf_norm(dx) < 3),
                                 detail::range_if_expr>::value,
                  "inf_norm(dx)<3 does not match range_if_expr");

    static_assert(
        proto::matches<decltype(dot(dx, dx) < 3), detail::range_if_expr>::value,
        "dot(dx,dx)<3 does not match range_if_expr");

    static_assert(
        proto::matches<decltype(id_[a] != id_[b] && norm(dx) < 3),
                       detail::range_if_expr>::value,
        "id_[a] != id_b[] && norm(dx)<3 does not match range_if_expr");

    static_assert(
        proto::matches<decltype(norm(dx) < 3 && id_[a] != id_[b]),
                       detail::range_if_expr>::value,
        "norm(dx)<3 && id_[a] != id_b[] does not match range_if_expr");

    static_assert(
        proto::matches<decltype(norm(dx) <= 3), detail::range_if_expr>::value,
        "norm(dx)<=3 does not match range_if_expr");

    static_assert(proto::matches<decltype(norm(dx) <= s[a]),
                                 detail::range_if_expr>::value,
                  "norm(dx)<=s[a] does not match range_if_expr");

    static_assert(proto::matches<decltype(dot(dx, dx) <= s[a]),
                                 detail::range_if_expr>::value,
                  "dot(dx,dx)<=s[a] does not match range_if_expr");

    static_assert(
        proto::matches<decltype(s[a] > norm(dx)), detail::range_if_expr>::value,
        "s[a] > norm(dx) does not match range_if_expr");

    static_assert(proto::matches<decltype(s[a] > dot(dx, dx)),
                                 detail::range_if_expr>::value,
                  "s[a] > dot(dx,dx) does not match range_if_expr");

    static_assert(
        proto::matches<decltype(10 > norm(dx)), detail::range_if_expr>::value,
        "10 > norm(dx) does not match range_if_expr");

    static_assert(
        proto::matches<decltype(10 >= norm(dx)), detail::range_if_expr>::value,
        "10 >= norm(dx) does not match range_if_expr");

    static_assert(
        !proto::matches<decltype(norm(dx) > 3), detail::range_if_expr>::value,
        "norm(dx)>3 matchs range_if_expr");

    static_assert(
        !proto::matches<decltype(dx > 3), detail::range_if_expr>::value,
        "dx>3 matchs range_if_expr");

    static_assert(
        !proto::matches<decltype(norm(dx) >= 3), detail::range_if_expr>::value,
        "norm(dx)>=3 matchs range_if_expr");

    static_assert(
        !proto::matches<decltype(s[a] == 3), detail::range_if_expr>::value,
        "s[a]==3 matchs range_if_expr");

    static_assert(!proto::matches<decltype(proto::lit(true)),
                                  detail::range_if_expr>::value,
                  "lit(true) matchs range_if_expr");

    static_assert(
        !proto::matches<decltype(proto::lit(1)), detail::range_if_expr>::value,
        "lit(1) matchs range_if_expr");
  }

  void test_get_labels(void) {
    ABORIA_VARIABLE(scalar, double, "scalar")
    typedef Particles<std::tuple<scalar>> ParticlesType;
    ParticlesType particles;
    Symbol<scalar> s;
    Label<0, ParticlesType> a(particles);
    Label<1, ParticlesType> b(particles);
    auto dx = create_dx(a, b);
    Accumulate<std::plus<double>> sum_dense;
    AccumulateWithinDistance<std::plus<double>> sum_sparse(2);

    // test constant expressions behave as expected
    auto nil_vect = detail::get_labels()(proto::lit(1), fusion::nil_());
#if not defined(__CUDACC__)
    BOOST_MPL_ASSERT_MSG(
        (std::is_same<decltype(nil_vect), fusion::nil_>::value),
        RESULT_OF_GET_LABELS_ON_CONSTANT_EXPRESSION_IS_NOT_A_NIL_LIST,
        (decltype(nil_vect)));
#endif
    static_assert(detail::is_const<proto::literal<int>>::value,
                  "result of is_const on constant expression is not true");
    static_assert(!detail::is_univariate<proto::literal<int>>::value,
                  "result of is_univariate on constant expression is true");
    static_assert(!detail::is_bivariate<proto::literal<int>>::value,
                  "result of is_univariate on constant expression is true");

    // test univariate expressions behave as expected
    auto vect_w_a = detail::get_labels()(s[a], fusion::nil_());

    typedef const detail::label<boost::mpl::int_<0>, ParticlesType>
        &expected_type_a;
    typedef const detail::label<boost::mpl::int_<1>, ParticlesType>
        &expected_type_b;
#if not defined(__CUDACC__)
    BOOST_MPL_ASSERT_MSG(
        (std::is_same<decltype(vect_w_a),
                      fusion::cons<expected_type_a>>::value),
        RESULT_OF_GET_LABELS_ON_UNIVARIATE_EXPRESSION_S_A_IS_NOT_CORRECT,
        (decltype(vect_w_a)));
#endif
    static_assert(
        std::is_same<std::result_of<detail::get_labels(decltype(s[a]),
                                                       fusion::nil_)>::type,
                     fusion::cons<expected_type_a>>::value,
        "result of get_labels on univariate expression s[a] is not correct");
#if not defined(__CUDACC__)
    BOOST_MPL_ASSERT_MSG(
        (std::is_same<std::result_of<detail::get_labels(decltype(s[a] + s[a]),
                                                        fusion::nil_)>::type,
                      fusion::cons<expected_type_a>>::value),
        RESULT_OF_GET_LABELS_ON_UNIVARIATE_EXPRESSION_S_A_PLUS_S_A_IS_NOT_CORRECT,
        (std::result_of<detail::get_labels(decltype(s[a] + s[a]),
                                           fusion::nil_)>::type));
#endif
    static_assert(detail::is_univariate<decltype(s[a])>::value,
                  "result of is_univariate on expression s[a] is not true");
    static_assert(
        detail::is_univariate<decltype(s[a] + s[a])>::value,
        "result of is_univariate on expression s[a]+s[a] is not true");
    static_assert(!detail::is_const<decltype(s[a] + s[a])>::value,
                  "result of is_const on expression s[a]+s[a] is true");
    static_assert(!detail::is_bivariate<decltype(s[a] + s[a])>::value,
                  "result of is_bivariate on expression s[a]+s[a] is true");

    // test bivariate expressions behave as expected
    auto vect_w_a_and_b = detail::get_labels()(s[a] + s[b], fusion::nil_());
#if not defined(__CUDACC__)
    BOOST_MPL_ASSERT_MSG(
        (std::is_same<decltype(vect_w_a_and_b),
                      fusion::cons<expected_type_b,
                                   fusion::cons<expected_type_a>>>::value),
        RESULT_OF_GET_LABELS_ON_BIVARIATE_EXPRESSION_S_A_PLUS_S_B_IS_NOT_CORRECT,
        (decltype(vect_w_a_and_b)));
#endif
    static_assert(detail::is_bivariate<decltype(s[a] + s[b])>::value,
                  "result of is_bivariate on expression s[a]+s[b] is not true");
    auto vect_w_a_and_b2 =
        detail::get_labels()(s[a] + s[b] + s[a] * s[b], fusion::nil_());
    static_assert(
        std::is_same<decltype(vect_w_a_and_b2),
                     fusion::cons<expected_type_b,
                                  fusion::cons<expected_type_a>>>::value,
        "result of get_labels on bivariate expression s[a]+s[b]+s[a]*s[b] is "
        "not correct");
    static_assert(
        detail::is_bivariate<decltype(s[a] + s[b] + s[a] * s[b])>::value,
        "result of is_bivariate on expression s[a]+s[b]+s[a]*s[b] is not true");

    // check if_else
    static_assert(
        detail::is_bivariate<decltype(if_else(s[a] == 1, 2, s[b]))>::value,
        "result of is_bivariate on expression if_else(s[a]==1,2,s[b]) is not "
        "true");

    // check dx
    static_assert(detail::is_bivariate<decltype(dx)>::value,
                  "result of is_bivariate on expression dx is not true");

    // check sums
    auto vect_w_a_and_b3 =
        detail::get_labels()(sum_sparse(b, s[a] + s[a]), fusion::nil_());
#if not defined(__CUDACC__)
    BOOST_MPL_ASSERT_MSG(
        (std::is_same<decltype(vect_w_a_and_b3),
                      fusion::cons<expected_type_a>>::value),
        RESULT_OF_GET_LABELS_ON_UNIVARIATE_EXPRESSION_SUM_IS_NOT_CORRECT,
        (decltype(vect_w_a_and_b3)));
#endif

    // check sums
    auto vect_w_a_and_b4 =
        detail::get_labels()(sum_dense(a, sum_sparse(b, 1)), fusion::nil_());
#if not defined(__CUDACC__)
    BOOST_MPL_ASSERT_MSG(
        (std::is_same<decltype(vect_w_a_and_b4), fusion::nil_>::value),
        RESULT_OF_GET_LABELS_ON_CONST_EXPRESSION_SUM_IS_NOT_CORRECT,
        (decltype(vect_w_a_and_b4)));
#endif

    static_assert(
        detail::is_univariate<decltype(sum_sparse(b, s[a] + s[a]))>::value,
        "result of is_univariate on expression sum_sparse(b,s[a]+s[a]) is not "
        "true");

    static_assert(
        detail::is_univariate<decltype(sum_sparse(b, norm(dx)))>::value,
        "result of is_univariate on expression sum_sparse(b,norm(dx)) is not "
        "true");

    static_assert(detail::is_const<decltype(
                      sum_dense(b, if_else(s[b] == 1, s[b], 0)))>::value,
                  "result of is_const on expression "
                  "sum_dense(b,if_else(s[b]==1,s[b],0)) is not true");

    static_assert(
        detail::is_const<decltype(sum_dense(a, sum_sparse(b, 1)))>::value,
        "result of is_const on expression "
        "sum_dense(a,sum_sparse(b,1)) is not true");
  }

  void test_thrust(void) {
#if defined(__aboria_use_thrust_algorithms__) || defined(__CUDACC__)
    ABORIA_VARIABLE(scalar, double, "scalar")
    typedef Particles<std::tuple<scalar>> ParticlesType;
    ParticlesType particles;
    auto b = particles.begin();

    static_assert(!(thrust::detail::is_trivial_iterator<decltype(b)>::value),
                  "aboria zip iterator is a trivial iterator!");
    static_assert(!(thrust::detail::is_pointer<decltype(b)>::value),
                  "aboria zip iterator is a pointer!");
    static_assert(!(thrust::detail::is_thrust_pointer<decltype(b)>::value),
                  "aboria zip iterator is a thrust pointer!");

    thrust::device_vector<float> A0;
    thrust::device_vector<float> A1;
    auto thrust_zip =
        thrust::make_zip_iterator(thrust::make_tuple(A0.begin(), A1.begin()));
    static_assert(
        !(thrust::detail::is_trivial_iterator<decltype(thrust_zip)>::value),
        "thrust zip iterator is a trivial iterator!");
#endif
  }
};

#endif /* SYMBOLICTEST_H_ */
