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

#ifndef BEM_TEST_H_
#define BEM_TEST_H_

#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <cxxtest/TestSuite.h>

//[bem
/*`

Here we do BEM

*/
#include "Aboria.h"
#include <random>

using namespace Aboria;

//<-
class BEMTest : public CxxTest::TestSuite {
public:
  void test_BEM(void) {
#ifdef HAVE_EIGEN &&HAVE_H2LIB
    std::cout << "---------------------\n"
              << "Running BEM test....\n"
              << "---------------------" << std::endl;
    //->
    //=int main() {
    auto funct = [](const double x, const double y) {
      return std::exp(-9 * std::pow(x - 0.5, 2) - 9 * std::pow(y - 0.25, 2));
      // return x;
    };

    typedef Eigen::Matrix<double, 2, 2> mat2x2;

    vdouble2 min(0);
    vdouble2 max(1);
    const double alpha = (max - min).prod() / (4.0 * pi);
    const double inv4alpha = 1.0 / (4.0 * alpha);
    const double inv4sqrtalpha = 1.0 / (4.0 * std::sqrt(alpha));
    const double sqrtpialpha = std::sqrt(pi * alpha);

    auto integrate_wave_space = [](const vdouble2 &K, const vdouble2 &dx_a,
                                   const vdouble2 &dx_b) {
      mat2x2 result;
      const double k2 = K.squaredNorm();
      const vdouble2 xd = dx_b - dx_a;
      const double sin_plus_sin = std::sin(dx_b[0] * K[0] + dx_b[1] * K[1]) -
                                  std::sin(dx_a[0] * K[0] + dx_a[1] * K[1]);
      const double xd_dot_K = xd[0] * K[0] + xd[1] * K[1];
      const double C = 4.0 * pi * h * (1 + alpha * k2) * std::exp(-alpha * k2) *
                       sin_plus_sin / (tau * std::pow(k2, 2) * xd_dot_K);
      result(0, 0) = C * (-k2 + K[0] * K[0]);
      result(0, 1) = C * (K[0] * K[1]);
      result(1, 0) = result(0, 1);
      result(1, 1) = C * (-k2 + K[1] * K[1]);
      return result;
    };

    auto integrate_real_space0 =
        [](const double_d &dx_a, const double_d &dx_b) {
          mat2x2 result;
          vdouble2 d = dx_b - dx_a;
          const double h = d.norm();
          d /= h;
          const double E1 =
              boost::math::expint(1, std::pow(h, 2) * 0.25 * inv4alpha);
          const double erfh = boost::math::erf(h * inv4sqrtalpha);
          const double C1 = 0.5 * h * E1;
          const double C2 = 2 * sqrtpialpha * erfh;
          result(0, 0) = C1 + C2 * d[0] * d[0];
          result(0, 1) = C2 * d[0] * d[1];
          result(1, 0) = result(1, 0);
          result(1, 1) = C1 + C2 * d[1] * d[1];
          return result;
        }

    // x = x - Llamba
    auto integrate_real_space = [](const double_d &dx_a, const double_d &dx_b) {
      const vdouble2 d = dx_b - dx_a;

      auto f =
          [&](const double t) {
            mat2x2 result;
            const vdouble2 x = d * t + dx_a;
            const double r2 = x.squaredNorm();
            const double E1 = boost::math::expint(1, r2 * inv4alpha);
            const double exp = std::exp(-r2 * inv4alpha);
            result(0, 0) = 0.5 * E1 + (x[0] * x[0] / r2 - 1) * exp;
            result(0, 1) = (x[0] * x[1] / r2) * exp;
            result(1, 0) = result(0, 1);
            result(0, 0) = 0.5 * E1 + (x[1] * x[1] / r2 - 1) * exp;
            return result;
          }

          return d.norm() *
          boost::math::quadrature::gauss<mat2x2, 7>::integrate(f, 0, 1);
    };

    ABORIA_VARIABLE(traction0, double, "traction0")
    ABORIA_VARIABLE(traction1, double, "traction1")
    ABORIA_VARIABLE(velocity0, double, "velocity0")
    ABORIA_VARIABLE(velocity1, double, "velocity1")

    typedef Particles<std::tuple<traction0, traction1, ve>, 2, std::vector,
                      SearchMethod>
        ParticlesType;
    typedef position_d<2> position;
    typedef typename ParticlesType::const_reference const_particle_reference;
    typedef typename position::value_type const &const_position_reference;
    typedef Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> map_type;
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_type;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_type;

    const int Ntest = 100;
    const int Nboundary = 100;

    ParticlesType boundary(Nboundary);
    ParticlesType test(Ntest * Ntest);

    const double dx = (max[0] - min[0]) / Ntest;
    const vdouble2 centre = 0.5 * (max - min);
    const double r = 0.2 * (max[0] - min[0]);
    for (int i = 0; i < Ntest; ++i) {
      for (int j = 0; j < Ntest; ++j) {
        get<position>(test)[i * Ntest + j] = vdouble2(i * dx, j * dx);
      }
    }
    for (int i = 0; i < Nboundary; ++i) {
      const double theta = i * 2 * pi / Nboundary;
      get<position>(test)[i * Ntest + j] =
          centre + r * vdouble2(std::cos(theta), std::sin(theta));
    }

    matrix_type A(2 * Nboundary, 2 * Nboundary);
    for (int i = 0; i < Nboundary; ++i) {
      for (int j = 0; j < Nboundary; ++j) {

        // knots.init_neighbour_search(min,max,periodic);

        augment.push_back(p);

        auto kernel = [](const_position_reference dx,
                         const_particle_reference a,
                         const_particle_reference b) {
          return std::sqrt(dx.squaredNorm() + get<constant2>(b));
        };

        auto one = [](const_position_reference dx, const_particle_reference a,
                      const_particle_reference b) { return 1.0; };

        auto G = create_dense_operator(knots, knots, kernel);
        auto P = create_dense_operator(knots, augment, one);
        auto Pt = create_dense_operator(augment, knots, one);
        auto Zero = create_zero_operator(augment, augment);

        auto W = create_block_operator<2, 2>(G, P, Pt, Zero);

        auto G_test = create_dense_operator(test, knots, kernel);
        auto W_test = create_block_operator<2, 2>(G_test, P, Pt, Zero);

        vector_type phi(N + 1), gamma(N + 1);
        for (size_t i = 0; i < knots.size(); ++i) {
          const double x = get<position>(knots[i])[0];
          const double y = get<position>(knots[i])[1];
          phi[i] = funct(x, y);
        }
        phi[knots.size()] = 0;

        matrix_type W_matrix(N + 1, N + 1);
        W.assemble(W_matrix);

        gamma = W_matrix.ldlt().solve(phi);

        Eigen::GMRES<matrix_type> gmres;
        gmres.setMaxIterations(max_iter);
        gmres.set_restart(restart);
        gmres.compute(W_matrix);
        gamma = gmres.solve(phi);
        std::cout << "GMRES:       #iterations: " << gmres.iterations()
                  << ", estimated error: " << gmres.error() << std::endl;

        phi = W * gamma;
        double rms_error = 0;
        double scale = 0;
        for (size_t i = 0; i < knots.size(); ++i) {
          const double x = get<position>(knots[i])[0];
          const double y = get<position>(knots[i])[1];
          const double truth = funct(x, y);
          const double eval_value = phi[i];
          rms_error += std::pow(eval_value - truth, 2);
          scale += std::pow(truth, 2);
          // TS_ASSERT_DELTA(eval_value,truth,2e-3);
        }

        std::cout << "rms_error for global support, at centers  = "
                  << std::sqrt(rms_error / scale) << std::endl;
        TS_ASSERT_LESS_THAN(std::sqrt(rms_error / scale), 1e-6);

        phi = W_test * gamma;
        rms_error = 0;
        scale = 0;
        for (size_t i = 0; i < test.size(); ++i) {
          const double x = get<position>(test[i])[0];
          const double y = get<position>(test[i])[1];
          const double truth = funct(x, y);
          const double eval_value = phi[i];
          rms_error += std::pow(eval_value - truth, 2);
          scale += std::pow(truth, 2);
          // TS_ASSERT_DELTA(eval_value,truth,2e-3);
        }
        std::cout << "rms_error for global support, away from centers  = "
                  << std::sqrt(rms_error / scale) << std::endl;
        TS_ASSERT_LESS_THAN(std::sqrt(rms_error / scale), 1e-5);

//=}
//]
#endif // HAVE_EIGEN
      }

      template <template <typename> class SearchMethod>
      void helper_compact(void) {
#ifdef HAVE_EIGEN
        std::cout << "---------------------\n"
                  << "Running compact test....\n"
                  << "---------------------" << std::endl;

        //[rbf_interpolation_compact
        //=#include "Aboria.h"
        //=#include <random>
        //=int main() {
        auto funct = [](const double x, const double y) {
          return std::exp(-9 * std::pow(x - 0.5, 2) -
                          9 * std::pow(y - 0.25, 2));
          // return x;
        };

        ABORIA_VARIABLE(alpha, double, "alpha value")
        ABORIA_VARIABLE(interpolated, double, "interpolated value")
        ABORIA_VARIABLE(constant2, double, "c2 value")

        typedef Particles<std::tuple<alpha, constant2, interpolated>, 2,
                          std::vector, SearchMethod>
            ParticlesType;
        typedef position_d<2> position;
        typedef
            typename ParticlesType::const_reference const_particle_reference;
        typedef typename position::value_type const &const_position_reference;
        typedef Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> map_type;
        typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_type;
        typedef Eigen::SparseMatrix<double> matrix_type;
        ParticlesType knots, augment;
        ParticlesType test;

        const double hfac = 4.0;
        vdouble2 min(0);
        vdouble2 max(1);
        vdouble2 periodic(false);

        const int N = 1000;

        const int max_iter = 100;
        const int restart = 101;
        const double delta = std::pow(double(N), -1.0 / 2.0);
        const double h = hfac * delta;
        std::cout << "using h = " << h << std::endl;
        const double RASM_size = 2 * h;
        const int RASM_n = N * std::pow(RASM_size, 2) / (max - min).prod();
        const double RASM_buffer = 0.9 * RASM_size;
        std::cout << "RASM_size = " << RASM_size << " RASM_n = " << RASM_n
                  << " RASM_buffer = " << RASM_buffer << std::endl;

        typename ParticlesType::value_type p;

        std::default_random_engine generator(123);
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        for (int i = 0; i < N; ++i) {
          get<position>(p) =
              vdouble2(distribution(generator), distribution(generator));
          knots.push_back(p);

          get<position>(p) =
              vdouble2(distribution(generator), distribution(generator));
          test.push_back(p);
        }
        augment.push_back(p);

        knots.init_neighbour_search(min, max, periodic);

        auto kernel = [h](const_position_reference dx,
                          const_particle_reference a,
                          const_particle_reference b) {
          return std::pow(2.0 - dx.norm() / h, 4) * (1.0 + 2.0 * dx.norm() / h);
        };

        auto one = [](const_position_reference dx, const_particle_reference a,
                      const_particle_reference b) { return 1.0; };

        auto G = create_sparse_operator(knots, knots, 2 * h, kernel);
        auto P = create_dense_operator(knots, augment, one);
        auto Pt = create_dense_operator(augment, knots, one);
        auto Zero = create_zero_operator(augment, augment);

        auto W = create_block_operator<2, 2>(G, P, Pt, Zero);

        auto G_test = create_sparse_operator(test, knots, 2 * h, kernel);
        auto W_test = create_block_operator<2, 2>(G_test, P, Pt, Zero);

        vector_type phi(N + 1), gamma(N + 1);
        for (size_t i = 0; i < knots.size(); ++i) {
          const double x = get<position>(knots[i])[0];
          const double y = get<position>(knots[i])[1];
          phi[i] = funct(x, y);
        }
        phi[knots.size()] = 0;

        matrix_type W_matrix(N + 1, N + 1);
        W.assemble(W_matrix);

        Eigen::ConjugateGradient<matrix_type, Eigen::Lower | Eigen::Upper,
                                 Eigen::DiagonalPreconditioner<double>>
            cg_test;
        cg_test.setMaxIterations(max_iter);
        cg_test.compute(W_matrix);
        gamma = cg_test.solve(phi);
        std::cout << "CG:          #iterations: " << cg_test.iterations()
                  << ", estimated error: " << cg_test.error() << std::endl;

        Eigen::ConjugateGradient<matrix_type, Eigen::Lower | Eigen::Upper,
                                 RASMPreconditioner<Eigen::HouseholderQR>>
            cg;
        cg.setMaxIterations(max_iter);
        cg.preconditioner().set_buffer_size(RASM_buffer);
        cg.preconditioner().set_number_of_particles_per_domain(RASM_n);
        cg.preconditioner().analyzePattern(W);
        cg.compute(W_matrix);
        gamma = cg.solve(phi);
        std::cout << "CG-RASM:     #iterations: " << cg.iterations()
                  << ", estimated error: " << cg.error() << std::endl;

        Eigen::MINRES<matrix_type, Eigen::Lower | Eigen::Upper,
                      RASMPreconditioner<Eigen::HouseholderQR>>
            minres;
        minres.setMaxIterations(max_iter);
        minres.preconditioner().set_buffer_size(RASM_buffer);
        minres.preconditioner().set_number_of_particles_per_domain(RASM_n);
        minres.preconditioner().analyzePattern(W);
        minres.compute(W_matrix);
        gamma = minres.solve(phi);
        std::cout << "MINRES-RASM: #iterations: " << minres.iterations()
                  << ", estimated error: " << minres.error() << std::endl;

        Eigen::GMRES<matrix_type, RASMPreconditioner<Eigen::HouseholderQR>>
            gmres;
        gmres.setMaxIterations(max_iter);
        gmres.preconditioner().set_buffer_size(RASM_buffer);
        gmres.preconditioner().set_number_of_particles_per_domain(RASM_n);
        gmres.preconditioner().analyzePattern(W);
        gmres.set_restart(restart);
        gmres.compute(W_matrix);
        gamma = gmres.solve(phi);
        std::cout << "GMRES-RASM:  #iterations: " << gmres.iterations()
                  << ", estimated error: " << gmres.error() << std::endl;

        Eigen::DGMRES<matrix_type, RASMPreconditioner<Eigen::HouseholderQR>>
            dgmres;
        dgmres.setMaxIterations(max_iter);
        dgmres.preconditioner().set_buffer_size(RASM_buffer);
        dgmres.preconditioner().set_number_of_particles_per_domain(RASM_n);
        dgmres.preconditioner().analyzePattern(W);
        dgmres.set_restart(restart);
        dgmres.compute(W_matrix);
        gamma = dgmres.solve(phi);
        std::cout << "DGMRES-RASM:  #iterations: " << dgmres.iterations()
                  << ", estimated error: " << dgmres.error() << std::endl;

        Eigen::BiCGSTAB<matrix_type, RASMPreconditioner<Eigen::HouseholderQR>>
            bicg;
        bicg.setMaxIterations(max_iter);
        bicg.preconditioner().set_buffer_size(RASM_buffer);
        bicg.preconditioner().set_number_of_particles_per_domain(RASM_n);
        bicg.preconditioner().analyzePattern(W);
        bicg.compute(W_matrix);
        gamma = bicg.solve(phi);
        std::cout << "BiCGSTAB-RASM:#iterations: " << bicg.iterations()
                  << ", estimated error: " << bicg.error() << std::endl;

        double rms_error = 0;
        double scale = 0;
        phi = W * gamma;
        for (size_t i = 0; i < knots.size(); ++i) {
          const double x = get<position>(knots[i])[0];
          const double y = get<position>(knots[i])[1];
          const double truth = funct(x, y);
          const double eval_value = phi[i];
          rms_error += std::pow(eval_value - truth, 2);
          scale += std::pow(truth, 2);
          // TS_ASSERT_DELTA(eval_value,truth,2e-3);
        }
        std::cout << "rms_error for compact support, at centers  = "
                  << std::sqrt(rms_error / scale) << std::endl;
        TS_ASSERT_LESS_THAN(std::sqrt(rms_error / scale), 1e-4);

        rms_error = 0;
        scale = 0;
        phi = W_test * gamma;
        for (size_t i = 0; i < test.size(); ++i) {
          const double x = get<position>(test[i])[0];
          const double y = get<position>(test[i])[1];
          const double truth = funct(x, y);
          const double eval_value = phi[i];
          rms_error += std::pow(eval_value - truth, 2);
          scale += std::pow(truth, 2);
          // TS_ASSERT_DELTA(eval_value,truth,2e-3);
        }

        std::cout << "rms_error for compact support, away from centers  = "
                  << std::sqrt(rms_error / scale) << std::endl;
        TS_ASSERT_LESS_THAN(std::sqrt(rms_error / scale), 1e-2);

//=}
//]
#endif // HAVE_EIGEN
      }

      template <template <typename> class SearchMethod> void helper_h2(void) {
#ifdef HAVE_EIGEN
        std::cout << "---------------------\n"
                  << "Running h2 test....\n"
                  << "---------------------" << std::endl;
        //[rbf_interpolation_h2
        //=#include "Aboria.h"
        //=#include <random>
        //=int main() {
        auto funct = [](const double x, const double y) {
          return std::exp(-9 * std::pow(x - 0.5, 2) -
                          9 * std::pow(y - 0.25, 2));
          // return x;
        };

        ABORIA_VARIABLE(alpha, double, "alpha value")
        ABORIA_VARIABLE(interpolated, double, "interpolated value")

        typedef Particles<std::tuple<alpha, interpolated>, 2, std::vector,
                          SearchMethod>
            ParticlesType;
        typedef position_d<2> position;
        typedef
            typename ParticlesType::const_reference const_particle_reference;
        typedef typename position::value_type const &const_position_reference;
        typedef Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> map_type;
        typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_type;
        typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
            matrix_type;
        ParticlesType knots;
        ParticlesType augment;
        ParticlesType test;

        const double c = 3.0;
        const double c2 = std::pow(c, 2);
        vdouble2 min(0);
        vdouble2 max(1);
        vdouble2 periodic(false);

        const int N = 10000;

        const double RASM_size = 0.3 / c;
        const int RASM_n = N * std::pow(RASM_size, 2) / (max - min).prod();
        const double RASM_buffer = 0.9 * RASM_size;
        // std::cout << "RASM_size = "<<RASM_size<<" RASM_n = "<<RASM_n<<"
        // RASM_buffer = "<<RASM_buffer<<std::endl;

        const int nx = 3;
        const int max_iter = 100;
        const int restart = 101;
        const double delta = 1.0 / nx;
        typename ParticlesType::value_type p;

        std::default_random_engine generator;
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        for (int i = 0; i < N; ++i) {
          get<position>(p) =
              vdouble2(distribution(generator), distribution(generator));
          knots.push_back(p);

          get<position>(p) =
              vdouble2(distribution(generator), distribution(generator));
          test.push_back(p);
        }

        knots.init_neighbour_search(min, max, periodic);

        augment.push_back(p);

        auto kernel = [&](const vdouble2 &dx, const vdouble2 &a,
                          const vdouble2 &b) {
          return std::exp(-dx.squaredNorm() * c2);
        };

        auto one = [](const_position_reference dx, const_particle_reference a,
                      const_particle_reference b) { return 1.0; };

        auto G = create_h2_operator<4>(knots, knots, kernel);
        auto P = create_dense_operator(knots, augment, one);
        auto Pt = create_dense_operator(augment, knots, one);
        auto Zero = create_zero_operator(augment, augment);

        auto W = create_block_operator<2, 2>(G, P, Pt, Zero);

        auto Gtest = create_h2_operator(G.get_first_kernel(), test);
        auto Ptest = create_dense_operator(test, augment, one);
        auto Wtest = create_block_operator<2, 2>(Gtest, Ptest, Pt, Zero);

        vector_type phi(N), gamma(N);
        for (size_t i = 0; i < knots.size(); ++i) {
          const double x = get<position>(knots[i])[0];
          const double y = get<position>(knots[i])[1];
          phi[i] = funct(x, y);
        }
        // phi[knots.size()] = 0;

        // Eigen::BiCGSTAB<decltype(W),
        // RASMPreconditioner<Eigen::HouseholderQR>> bicg;
        // Eigen::BiCGSTAB<decltype(W), ExtMatrixPreconditioner<2>> bicg;
        // bicg.preconditioner().set_buffer_size(RASM_buffer);
        // bicg.preconditioner().set_number_of_particles_per_domain(RASM_n);
        // bicg.setMaxIterations(max_iter);
        // bicg.compute(W);
        // gamma = bicg.solve(phi);
        // std::cout << "BiCGSTAB-RASM:#iterations: " << bicg.iterations() << ",
        // estimated error: " << bicg.error() << std::endl;

        // Eigen::DGMRES<decltype(W),  RASMPreconditioner<Eigen::HouseholderQR>>
        // dgmres;
        Eigen::DGMRES<decltype(G), ExtMatrixPreconditioner<2>> dgmres;
        // Eigen::DGMRES<decltype(W)> dgmres;
        // dgmres.preconditioner().set_buffer_size(RASM_buffer);
        // dgmres.preconditioner().set_number_of_particles_per_domain(RASM_n);
        // dgmres.preconditioner().analyzePattern(W);
        dgmres.setMaxIterations(max_iter);
        dgmres.set_restart(restart);
        dgmres.compute(G);
        gamma = dgmres.solve(phi);
        std::cout << "DGMRES-EXT:  #iterations: " << dgmres.iterations()
                  << ", estimated error: " << dgmres.error()
                  << " true error = " << (G * gamma - phi).norm() << std::endl;

        phi = G * gamma;
        double rms_error = 0;
        double scale = 0;
        for (size_t i = 0; i < knots.size(); ++i) {
          const double x = get<position>(knots[i])[0];
          const double y = get<position>(knots[i])[1];
          const double truth = funct(x, y);
          const double eval_value = phi[i];
          rms_error += std::pow(eval_value - truth, 2);
          scale += std::pow(truth, 2);
          // TS_ASSERT_DELTA(eval_value,truth,2e-3);
        }

        std::cout << "rms_error for global support, at centers  = "
                  << std::sqrt(rms_error / scale) << std::endl;
        TS_ASSERT_LESS_THAN(std::sqrt(rms_error / scale), 1e-4);

        rms_error = 0;
        scale = 0;
        vector_type phi_test = Gtest * gamma;
        for (size_t i = 0; i < test.size(); ++i) {
          const vdouble2 p = get<position>(test)[i];
          const double eval_value = phi_test[i];
          const double truth = funct(p[0], p[1]);
          rms_error += std::pow(eval_value - truth, 2);
          scale += std::pow(truth, 2);
          // TS_ASSERT_DELTA(eval_value,truth,2e-3);
        }
        std::cout << "rms_error for global support, away from centers  = "
                  << std::sqrt(rms_error / scale) << std::endl;
        TS_ASSERT_LESS_THAN(std::sqrt(rms_error / scale), 1e-4);

//=}
//]
#endif // HAVE_EIGEN
      }

      void test_CellListOrdered() {
        std::cout << "-------------------------------------------\n"
                  << "Running tests on CellListOrdered....\n"
                  << "------------------------------------------" << std::endl;
        helper_global<CellListOrdered>();
        helper_compact<CellListOrdered>();
      }

      void test_CellList() {
        std::cout << "-------------------------------------------\n"
                  << "Running tests on CellList....\n"
                  << "------------------------------------------" << std::endl;
        helper_global<CellList>();
        helper_compact<CellList>();
      }

      void test_kdtree() {
        std::cout << "-------------------------------------------\n"
                  << "Running tests on kdtree....\n"
                  << "------------------------------------------" << std::endl;
        helper_compact<Kdtree>();
        helper_h2<Kdtree>();
      }

      void test_HyperOctree() {
        std::cout << "-------------------------------------------\n"
                  << "Running tests on HyperOctree....\n"
                  << "------------------------------------------" << std::endl;
        helper_h2<HyperOctree>();
        helper_compact<Kdtree>();
      }
    };

#endif /* RBF_INTERPOLATION_TEST_H_ */
    //->
