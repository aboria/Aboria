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

#ifndef CHEBYSHEV_TEST_H_
#define CHEBYSHEV_TEST_H_

#include <cxxtest/TestSuite.h>

#include <chrono>
#include <random>
#include <time.h>
typedef std::chrono::system_clock Clock;
#include "Chebyshev.h"
#include "FastMultipoleMethod.h"
#include "Kernels.h"
#include "Level1.h"

using namespace Aboria;

class ChebyshevTest : public CxxTest::TestSuite {
  ABORIA_VARIABLE(source, double, "source");
  ABORIA_VARIABLE(target_cheb, double, "target chebyshev");
  ABORIA_VARIABLE(target_manual, double, "target manual");
  ABORIA_VARIABLE(target_fmm, double, "target fmm");
  ABORIA_VARIABLE(target_h2, double, "target h2");

public:
  template <unsigned int N, typename ParticlesType, typename KernelFunction>
  void helper_fast_methods_calculate(ParticlesType &particles,
                                     const KernelFunction &kernel,
                                     const double scale) {

#ifdef HAVE_EIGEN
    const unsigned int dimension = ParticlesType::dimension;

    // perform the operation using chebyshev interpolation operator
    auto t0 = Clock::now();
    auto C = create_chebyshev_operator(particles, particles, N, kernel);
    auto t1 = Clock::now();
    std::chrono::duration<double> time_cheb_setup = t1 - t0;
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_type;
    typedef Eigen::Map<vector_type> map_type;
    map_type source_vect(get<source>(particles).data(), particles.size());
    map_type target_vect(get<target_cheb>(particles).data(), particles.size());

    t0 = Clock::now();
    target_vect = C * source_vect;
    t1 = Clock::now();
    std::chrono::duration<double> time_cheb_eval = t1 - t0;

    const double L2_cheb = std::inner_product(
        std::begin(get<target_cheb>(particles)),
        std::end(get<target_cheb>(particles)),
        std::begin(get<target_manual>(particles)), 0.0,
        [](const double t1, const double t2) { return t1 + t2; },
        [](const double t1, const double t2) {
          // std::cout << "t1 = " << t1 << "t2 = " << t2 << std::endl;

          return (t1 - t2) * (t1 - t2);
        });

    std::cout << "dimension = " << dimension << ". N = " << N
              << ". L2_cheb error = " << L2_cheb
              << ". L2_cheb relative error is " << std::sqrt(L2_cheb / scale)
              << ". time_cheb_setup  = " << time_cheb_setup.count()
              << ". time_cheb_eval = " << time_cheb_eval.count() << std::endl;

    // TODO: is there a better test than this, maybe shouldn't randomly do it?
    if (dimension == 3 && N >= 9)
      TS_ASSERT_LESS_THAN(std::sqrt(L2_cheb / scale), 0.01);
#endif
  }

  template <unsigned int D, template <typename, typename> class StorageVector,
            template <typename> class SearchMethod>
  void helper_fast_methods(size_t N) {
    typedef Vector<double, D> double_d;
    typedef Vector<bool, D> bool_d;
    // randomly generate a bunch of positions over a range
    const double pos_min = 0.0;
    const double pos_max = 1.0;
    std::uniform_real_distribution<double> U(pos_min, pos_max);
    generator_type generator(time(NULL));
    auto gen = std::bind(U, generator);

    typedef Particles<
        std::tuple<source, target_cheb, target_manual, target_fmm, target_h2>,
        D, StorageVector, SearchMethod>
        ParticlesType;
    typedef typename ParticlesType::position position;
    ParticlesType particles(N);

    for (size_t i = 0; i < N; i++) {
      for (size_t d = 0; d < D; ++d) {
        get<position>(particles)[i][d] = gen();
      }
      get<source>(particles)[i] = gen();
    }
    particles.init_neighbour_search(double_d::Constant(pos_min),
                                    double_d::Constant(pos_max),
                                    bool_d::Constant(false));

    // generate a source vector using a smooth cosine
    auto source_fn = [&](const double_d &p) {
      // return (p-double_d(0)).norm();
      double ret = 1.0;
      const double scale = 2.0 * detail::PI / (pos_max - pos_min);
      for (size_t i = 0; i < D; ++i) {
        ret *= cos((p[i] - pos_min) * scale);
      }
      return ret / N;
    };
    std::transform(std::begin(get<position>(particles)),
                   std::end(get<position>(particles)),
                   std::begin(get<source>(particles)), source_fn);

    const double c = 0.01;
    auto kernel = [&c](const double_d &pa, const double_d &pb) {
      return std::sqrt((pb - pa).squaredNorm() + c);
    };

    // perform the operation manually
    std::fill(std::begin(get<target_manual>(particles)),
              std::end(get<target_manual>(particles)), 0.0);

    auto t0 = Clock::now();
    for (size_t i = 0; i < N; i++) {
      const double_d pi = get<position>(particles)[i];
      for (size_t j = 0; j < N; j++) {
        const double_d pj = get<position>(particles)[j];
        get<target_manual>(particles)[i] +=
            kernel(pi, pj) * get<source>(particles)[j];
      }
    }
    auto t1 = Clock::now();
    std::chrono::duration<double> time_manual = t1 - t0;

    const double scale = std::accumulate(
        std::begin(get<target_manual>(particles)),
        std::end(get<target_manual>(particles)), 0.0,
        [](const double t1, const double t2) { return t1 + t2 * t2; });

    std::cout << "MANUAL TIMING: dimension = " << D
              << ". number of particles = " << N
              << ". time = " << time_manual.count() << " scale = " << scale
              << std::endl;

    helper_fast_methods_calculate<3>(particles, kernel, scale);
    helper_fast_methods_calculate<4>(particles, kernel, scale);
    helper_fast_methods_calculate<5>(particles, kernel, scale);
    helper_fast_methods_calculate<6>(particles, kernel, scale);
    helper_fast_methods_calculate<7>(particles, kernel, scale);
    helper_fast_methods_calculate<8>(particles, kernel, scale);
    helper_fast_methods_calculate<9>(particles, kernel, scale);
    helper_fast_methods_calculate<10>(particles, kernel, scale);
  }

  template <unsigned int D> void helper_Rn_calculation(void) {
    const double tol = 1e-10;
    // randomly generate a bunch of positions over a range
    std::uniform_real_distribution<double> U(-10, 100);
    generator_type generator(time(NULL));
    typedef Vector<double, D> double_d;
    typedef Vector<int, D> int_d;
    const size_t N = 50;
    std::vector<double_d> positions(N);
    for (size_t i = 0; i < N; i++) {
      for (size_t d = 0; d < D; ++d) {
        positions[i][d] = U(generator);
      }
    }
    detail::Chebyshev_Rn<D> Rn;
    for (int n = 1; n < 10; ++n) {
      Rn.calculate_Sn(std::begin(positions), N, n);
      const int_d start = int_d::Constant(0);
      const int_d end = int_d::Constant(n);
      auto range = iterator_range<lattice_iterator<D>>(
          lattice_iterator<D>(start, end), lattice_iterator<D>());
      const double_d scale =
          double_d::Constant(1.0) / (Rn.box.bmax - Rn.box.bmin);
      for (size_t i = 0; i < positions.size(); ++i) {
        const double_d &x =
            (2 * positions[i] - Rn.box.bmin - Rn.box.bmax) * scale;
        for (const int_d &m : range) {
          TS_ASSERT_DELTA(Rn(m, i), detail::chebyshev_Rn_slow(x, m, n), tol);
        }
      }
    }
    const double_d scale =
        double_d::Constant(1.0) / (Rn.box.bmax - Rn.box.bmin);
    for (size_t i = 0; i < positions.size(); ++i) {
      const unsigned int n = 4;
      const double_d &x =
          (2 * positions[i] - Rn.box.bmin - Rn.box.bmax) * scale;
      detail::ChebyshevRnSingle<D, n> cheb_rn(positions[i], Rn.box);
      const int_d start = int_d::Constant(0);
      const int_d end = int_d::Constant(n);
      auto range = iterator_range<lattice_iterator<D>>(
          lattice_iterator<D>(start, end), lattice_iterator<D>());

      for (const int_d &m : range) {
        TS_ASSERT_DELTA(cheb_rn(m), detail::chebyshev_Rn_slow(x, m, n), tol);
      }
    }
  }

  void test_chebyshev_polynomial_calculation(void) {
    const double tol = 1e-10;
    // evaluate polynomial of order k at i-th root
    // of polynomial of order n
    // should by cos(k*(2*i-1)/(2*n)*PI
    std::cout << "testing polynomial calculation..." << std::endl;
    const int n = 4;
    for (int i = 0; i < n; ++i) {
      const double x = cos((2.0 * i + 1.0) / (2.0 * n) * detail::PI);
      for (int k = 0; k < n; ++k) {
        TS_ASSERT_DELTA(detail::chebyshev_polynomial(x, k),
                        cos(k * (2.0 * i + 1.0) / (2.0 * n) * detail::PI), tol);
      }
    }
  }

  void test_fast_methods_CellList(void) {
    const size_t N = 1000;
    std::cout << "BUCKET_SEARCH_SERIAL: testing 1D..." << std::endl;
    helper_fast_methods<1, std::vector, CellList>(N);
    std::cout << "BUCKET_SEARCH_SERIAL: testing 2D..." << std::endl;
    helper_fast_methods<2, std::vector, CellList>(N);
    std::cout << "BUCKET_SEARCH_SERIAL: testing 3D..." << std::endl;
    helper_fast_methods<3, std::vector, CellList>(N);
  }

  void test_fast_methods_CellListOrdered(void) {
    const size_t N = 1000;
    std::cout << "BUCKET_SEARCH_PARALLEL: testing 1D..." << std::endl;
    helper_fast_methods<1, std::vector, CellListOrdered>(N);
    std::cout << "BUCKET_SEARCH_PARALLEL: testing 2D..." << std::endl;
    helper_fast_methods<2, std::vector, CellListOrdered>(N);
    std::cout << "BUCKET_SEARCH_PARALLEL: testing 3D..." << std::endl;
    helper_fast_methods<3, std::vector, CellListOrdered>(N);
  }

  void test_fast_methods_kd_tree(void) {
    const size_t N = 1000;
    std::cout << "KD_TREE: testing 1D..." << std::endl;
    helper_fast_methods<1, std::vector, Kdtree>(N);
    std::cout << "KD_TREE: testing 2D..." << std::endl;
    helper_fast_methods<2, std::vector, Kdtree>(N);
    std::cout << "KD_TREE: testing 3D..." << std::endl;
    helper_fast_methods<3, std::vector, Kdtree>(N);
  }

  void test_fast_methods_HyperOctree(void) {
    const size_t N = 1000;
    std::cout << "OCTTREE: testing 1D..." << std::endl;
    helper_fast_methods<1, std::vector, HyperOctree>(N);
    std::cout << "OCTTREE: testing 2D..." << std::endl;
    helper_fast_methods<2, std::vector, HyperOctree>(N);
    std::cout << "OCTTREE: testing 3D..." << std::endl;
    helper_fast_methods<3, std::vector, HyperOctree>(N);
  }

  void test_Rn_calculation(void) {
    std::cout << "testing 1D..." << std::endl;
    helper_Rn_calculation<1>();
    std::cout << "testing 2D..." << std::endl;
    helper_Rn_calculation<2>();
    std::cout << "testing 3D..." << std::endl;
    helper_Rn_calculation<3>();
    std::cout << "testing 4D..." << std::endl;
    helper_Rn_calculation<4>();
  }
};

#endif
