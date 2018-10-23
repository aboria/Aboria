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

#ifndef UTILS_H_
#define UTILS_H_

#include <cxxtest/TestSuite.h>

#include "Aboria.h"

#ifdef HAVE_EIGEN
#include "RedSVD/RedSVD.h"
#include "detail/LowRank.h"
#endif

using namespace Aboria;

class UtilsTest : public CxxTest::TestSuite {
public:
  void test_bucket_indicies(void) {
    typedef Vector<unsigned int, 3> vect;
    detail::bucket_index<3> bi(vect(4, 7, 2));
    unsigned int index = bi.collapse_index_vector(vect(1, 2, 1));
    // index = 1 + 2*2 + 1*2*7 = 19
    TS_ASSERT_EQUALS(index, 19);
    vect vindex = bi.reassemble_index_vector(index);
    TS_ASSERT_EQUALS(vindex[0], 1);
    TS_ASSERT_EQUALS(vindex[1], 2);
    TS_ASSERT_EQUALS(vindex[2], 1);

    typedef Vector<unsigned int, 4> vect2;
    detail::bucket_index<4> bi2(vect2(4, 7, 1, 6));
    index = bi2.collapse_index_vector(vect2(1, 2, 0, 4));
    // index = 4 + 0*6 + 2*6*1 + 1*6*1*7 = 58
    TS_ASSERT_EQUALS(index, 58);
    vect2 vindex2 = bi2.reassemble_index_vector(index);
    TS_ASSERT_EQUALS(vindex2[0], 1);
    TS_ASSERT_EQUALS(vindex2[1], 2);
    TS_ASSERT_EQUALS(vindex2[2], 0);
    TS_ASSERT_EQUALS(vindex2[3], 4);
  }

  void test_point_to_bucket_indicies(void) {
    const unsigned int D = 3;
    vdouble3 min(0, 0, 0);
    vdouble3 max(1, 1, 1);
    bbox<D> bounds(min, max);
    vint3 size(5, 5, 5);
    vdouble3 side_length = (max - min) / size;
    detail::point_to_bucket_index<3> ptob(size, side_length, bounds);

    vint3 point_bucket_index_true(2, 2, 2);
    vint3 point_bucket_index =
        ptob.find_bucket_index_vector(vdouble3(0.5, 0.5, 0.5));
    for (int i = 0; i < 3; ++i) {
      TS_ASSERT_EQUALS(point_bucket_index_true[i], point_bucket_index[i]);
    }

    int index_true = 2 * 5 * 5 + 2 * 5 + 2;
    int index = ptob(vdouble3(0.5, 0.5, 0.5));
    TS_ASSERT_EQUALS(index_true, index);

    /*
    int index2_true = 2;
    int index2 = ptob.get_min_index_by_quadrant(0.5 - 1e-5, 0, true);
    TS_ASSERT_EQUALS(index2_true, index2);

    int index3_true = 3;
    int index3 = ptob.get_min_index_by_quadrant(0.5 + 1e-5, 0, true);
    TS_ASSERT_EQUALS(index3_true, index3);

    int index4_true = 1;
    int index4 = ptob.get_min_index_by_quadrant(0.5 - 1e-5, 0, false);
    TS_ASSERT_EQUALS(index4_true, index4);

    int index5_true = 2;
    int index5 = ptob.get_min_index_by_quadrant(0.5 + 1e-5, 0, false);
    TS_ASSERT_EQUALS(index5_true, index5);
    */
  }

  void test_circle_cube_intersection(void) {
    vdouble2 centre(0, 0);
    double radius = 1;
    bbox<2> cube(vdouble2(-0.1, -0.1), vdouble2(0.1, 0.1));
    TS_ASSERT_EQUALS(circle_intersect_cube(centre, radius, cube), true);

    cube.bmin = vdouble2(-0.005, 0.11);
    cube.bmax = vdouble2(0.005, 0.12);
    TS_ASSERT_EQUALS(circle_intersect_cube(centre, radius, cube), true);

    cube.bmin = vdouble2(-1.1, -5);
    cube.bmax = vdouble2(-0.9, 5);
    TS_ASSERT_EQUALS(circle_intersect_cube(centre, radius, cube), true);

    cube.bmin = vdouble2(-1.1, -5);
    cube.bmax = vdouble2(-1.01, 5);
    TS_ASSERT_EQUALS(circle_intersect_cube(centre, radius, cube), false);

    cube.bmin = vdouble2(-1.1, -1.1);
    cube.bmax = vdouble2(1.1, 1.1);
    TS_ASSERT_EQUALS(circle_intersect_cube(centre, radius, cube), true);

    cube.bmin = vdouble2(0.99, 0.99);
    cube.bmax = vdouble2(0.999, 0.999);
    TS_ASSERT_EQUALS(circle_intersect_cube(centre, radius, cube), false);
  }

  void test_linear_transform(void) {
    auto skew = [](const vdouble2 &v) {
      return vdouble2(v[0] + 0.3 * v[1], v[1]);
    };
    auto t = create_linear_transform<2>(skew);
    TS_ASSERT_EQUALS(t.get_eigen_vertices()[0], t.get_eigen_vertices()[1]);

    bbox<2> unit(vdouble2::Constant(-1), vdouble2::Constant(1));
    auto unitt = t(unit);
    TS_ASSERT_EQUALS(unitt[0], 2.6);
    TS_ASSERT_EQUALS(unitt[1], 2);

    bbox<2> unit2(vdouble2::Constant(0), vdouble2::Constant(1));

    auto unitt3 = t(unit2);
    TS_ASSERT_EQUALS(unitt3[0], 1.3);
    TS_ASSERT_EQUALS(unitt3[1], 1);

    auto skew2 = [](const vdouble2 &v) {
      return vdouble2(v[0] - 0.3 * v[1], v[1]);
    };
    auto t2 = create_linear_transform<2>(skew2);
    TS_ASSERT_EQUALS(!(t2.get_eigen_vertices()[0]), t2.get_eigen_vertices()[1]);

    auto unitt2 = t2(unit);
    TS_ASSERT_EQUALS(unitt2[0], 2.6);
    TS_ASSERT_EQUALS(unitt2[1], 2);

    auto t3 = create_scale_transform(vdouble2(1, 2));
    auto v_transformed = t3(vdouble2(2, 2));
    TS_ASSERT_EQUALS(v_transformed[0], 2.0);
    TS_ASSERT_EQUALS(v_transformed[1], 4.0);

    auto unit_transformed = t3(unit);
    TS_ASSERT_EQUALS(unit_transformed[0], 2.0);
    TS_ASSERT_EQUALS(unit_transformed[1], 4.0);
  }

  void test_low_rank(void) {
#ifdef HAVE_EIGEN
    const unsigned int D = 2;
    const size_t N = 100;
    // randomly generate a bunch of positions over a range
    const double pos_min = 0;
    const double pos_max = 1;
    std::uniform_real_distribution<double> U1(pos_min, pos_max);
    std::uniform_real_distribution<double> U2(pos_min - 2, pos_max - 2);
    generator_type generator(time(NULL));
    auto gen1 = std::bind(U1, generator);
    auto gen2 = std::bind(U2, generator);

    typedef Particles<std::tuple<>, D> ParticlesType;
    typedef ParticlesType::const_reference const_reference;
    typedef typename ParticlesType::position position;
    ParticlesType particles1(N);
    ParticlesType particles2(N);
    const double error_aim = 0.01;

    for (size_t i = 0; i < N; i++) {
      for (size_t d = 0; d < D; ++d) {
        get<position>(particles1)[i][d] = gen1();
        get<position>(particles2)[i][d] = gen2();
      }
    }

    const double c = 0.01;
    auto kernel_fun = [&c](const_reference pa, const_reference pb) {
      return std::sqrt((get<position>(pb) - get<position>(pa)).squaredNorm() +
                       c);
    };

    KernelDense<ParticlesType, ParticlesType, decltype(kernel_fun)> kernel(
        particles1, particles2, kernel_fun);

    // fill a matrix with the result
    Eigen::Matrix<double, N, N> fixed_mat;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dyn_mat(N, N);

    kernel.assemble(fixed_mat);
    kernel.assemble(dyn_mat);

    double rms_error_scale = 0;
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < N; j++) {
        double Ztilde = kernel.coeff(i, j);
        rms_error_scale += std::pow(Ztilde, 2);
      }
    }

    for (int k = 1; k < 10; ++k) {
      Eigen::Matrix<double, N, Eigen::Dynamic> U(N, k);
      Eigen::Matrix<double, Eigen::Dynamic, N> V(k, N);

      Eigen::Matrix<double, N, N> fixed_mat_copy = fixed_mat;
      size_t est_k = detail::adaptive_cross_approximation_full(
          fixed_mat_copy, k, error_aim, U, V);

      // check accuracy
      double rms_error_full_fixed = 0;
      for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
          double Ztilde = 0;
          for (size_t kk = 0; kk < est_k; kk++) {
            Ztilde += U(i, kk) * V(kk, j);
          }
          rms_error_full_fixed += std::pow(Ztilde - fixed_mat(i, j), 2);
        }
      }

      std::cout << "fixed-ACA (full): k = " << k << " est_k = " << est_k
                << " rms error = "
                << std::sqrt(rms_error_full_fixed / rms_error_scale)
                << std::endl;

      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dyn_mat_copy =
          dyn_mat;
      est_k = detail::adaptive_cross_approximation_full(dyn_mat_copy, k,
                                                        error_aim, U, V);

      // check accuracy
      double rms_error_full_dyn = 0;
      for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
          double Ztilde = 0;
          for (size_t kk = 0; kk < est_k; kk++) {
            Ztilde += U(i, kk) * V(kk, j);
          }
          rms_error_full_dyn += std::pow(Ztilde - dyn_mat(i, j), 2);
        }
      }

      std::cout << "dyn-ACA (full): k = " << k << " est_k = " << est_k
                << " rms error = "
                << std::sqrt(rms_error_full_dyn / rms_error_scale) << std::endl;

      est_k = detail::adaptive_cross_approximation_partial(fixed_mat, k,
                                                           error_aim, U, V);

      // check accuracy
      double rms_error_partial_fixed = 0;
      for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
          double Ztilde = 0;
          for (size_t kk = 0; kk < est_k; kk++) {
            Ztilde += U(i, kk) * V(kk, j);
          }
          rms_error_partial_fixed += std::pow(Ztilde - fixed_mat(i, j), 2);
        }
      }

      std::cout << "fixed-ACA (partial): k = " << k << " est_k = " << est_k
                << " rms error = "
                << std::sqrt(rms_error_partial_fixed / rms_error_scale)
                << std::endl;

      est_k = detail::adaptive_cross_approximation_partial(dyn_mat, k,
                                                           error_aim, U, V);

      // check accuracy
      double rms_error_partial_dyn = 0;
      for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
          double Ztilde = 0;
          for (size_t kk = 0; kk < est_k; kk++) {
            Ztilde += U(i, kk) * V(kk, j);
          }
          rms_error_partial_dyn += std::pow(Ztilde - dyn_mat(i, j), 2);
        }
      }

      std::cout << "dyn-ACA (partial): k = " << k << " est_k = " << est_k
                << " rms error = "
                << std::sqrt(rms_error_partial_dyn / rms_error_scale)
                << std::endl;

      est_k = detail::adaptive_cross_approximation_partial(kernel, k, error_aim,
                                                           U, V);

      // check accuracy
      double rms_error_kernel = 0;
      for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
          double Ztilde = 0;
          for (size_t kk = 0; kk < est_k; kk++) {
            Ztilde += U(i, kk) * V(kk, j);
          }
          rms_error_kernel += std::pow(Ztilde - kernel.coeff(i, j), 2);
        }
      }

      std::cout << "kernel-ACA (partial): k = " << k << " est_k = " << est_k
                << " rms error = "
                << std::sqrt(rms_error_kernel / rms_error_scale) << std::endl;

      typedef RedSVD::RedSVD<decltype(fixed_mat)> fixed_svd_type;
      fixed_svd_type svd(fixed_mat, k);
      fixed_svd_type::ScalarVector diagonals = svd.singularValues();
      decltype(fixed_mat) fixed_mat_approx =
          svd.matrixU() *
          Eigen::DiagonalWrapper<fixed_svd_type::ScalarVector>(diagonals) *
          svd.matrixV().transpose();

      double rms_error_svd_fixed = 0;
      for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
          rms_error_svd_fixed +=
              std::pow(fixed_mat_approx(i, j) - fixed_mat(i, j), 2);
        }
      }

      std::cout << "fixed-RSVD: k = " << k << " rms error = "
                << std::sqrt(rms_error_svd_fixed / rms_error_scale)
                << std::endl;

      typedef RedSVD::RedSVD<decltype(dyn_mat)> dyn_svd_type;
      dyn_svd_type svd2(dyn_mat, k);
      dyn_svd_type::ScalarVector diagonals2 = svd2.singularValues();
      decltype(dyn_mat) dyn_mat_approx =
          svd2.matrixU() *
          Eigen::DiagonalWrapper<dyn_svd_type::ScalarVector>(diagonals) *
          svd2.matrixV().transpose();

      double rms_error_svd_dyn = 0;
      for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
          rms_error_svd_dyn +=
              std::pow(dyn_mat_approx(i, j) - dyn_mat(i, j), 2);
        }
      }

      std::cout << "dyn-RSVD: k = " << k << " rms error = "
                << std::sqrt(rms_error_svd_dyn / rms_error_scale) << std::endl;
      std::cout << "--------------------" << std::endl;
      if (k == 9) {
        TS_ASSERT_LESS_THAN(std::sqrt(rms_error_full_fixed / rms_error_scale),
                            error_aim);
        TS_ASSERT_LESS_THAN(std::sqrt(rms_error_full_dyn / rms_error_scale),
                            error_aim);
        TS_ASSERT_LESS_THAN(
            std::sqrt(rms_error_partial_fixed / rms_error_scale), error_aim);
        TS_ASSERT_LESS_THAN(std::sqrt(rms_error_partial_dyn / rms_error_scale),
                            error_aim);
        TS_ASSERT_LESS_THAN(std::sqrt(rms_error_kernel / rms_error_scale),
                            error_aim);
        TS_ASSERT_LESS_THAN(std::sqrt(rms_error_svd_fixed / rms_error_scale),
                            0.001);
        TS_ASSERT_LESS_THAN(std::sqrt(rms_error_svd_dyn / rms_error_scale),
                            0.001);
      }
    }

#endif
  }
};

#endif /* CONSTRUCTORS_H_ */
