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

#ifndef FMM_TEST_H_
#define FMM_TEST_H_

#include <cxxtest/TestSuite.h>

#include <chrono>
#include <random>
#include <time.h>
typedef std::chrono::system_clock Clock;
#include "Level1.h"
#ifdef HAVE_EIGEN
#include "Kernels.h"
#endif
#include "Chebyshev.h"
#include "FastMultipoleMethod.h"
#ifdef HAVE_GPERFTOOLS
#include <gperftools/profiler.h>
#endif

using namespace Aboria;

class FMMTest : public CxxTest::TestSuite {
  ABORIA_VARIABLE(source, double, "source");
  ABORIA_VARIABLE(target_manual, double, "target manual");
  ABORIA_VARIABLE(target_fmm, double, "target fmm");
  ABORIA_VARIABLE(vsource, vdouble2, "vector-valued source");
  ABORIA_VARIABLE(vtarget_manual, vdouble2, "vector-valued target manual");
  ABORIA_VARIABLE(vtarget_fmm, vdouble2, "vector-valued target fmm");

public:
  template <unsigned int N, typename Source, typename TargetManual,
            typename TargetFMM, typename ParticlesType, typename KernelFunction,
            typename P2PKernelFunction>
  void helper_fast_methods_calculate(ParticlesType &particles,
                                     const KernelFunction &kernel,
                                     const P2PKernelFunction &p2pkernel,
                                     const double scale) {
    const unsigned int dimension = ParticlesType::dimension;
    typedef typename TargetFMM::value_type value_type;
    typedef detail::VectorTraits<value_type> scalar_traits;

    auto t0 = Clock::now();
    auto fmm =
        make_fmm(particles, particles,
                 make_black_box_expansion<dimension, N>(kernel), p2pkernel);
    auto t1 = Clock::now();
    std::chrono::duration<double> time_fmm_setup = t1 - t0;
    std::fill(std::begin(get<TargetFMM>(particles)),
              std::end(get<TargetFMM>(particles)), scalar_traits::Zero());
    t0 = Clock::now();
    fmm.matrix_vector_multiply(get<TargetFMM>(particles),
                               get<Source>(particles));
    t1 = Clock::now();
    std::chrono::duration<double> time_fmm_eval = t1 - t0;

    double L2_fmm = std::inner_product(
        std::begin(get<TargetFMM>(particles)),
        std::end(get<TargetFMM>(particles)),
        std::begin(get<TargetManual>(particles)), 0.0,
        [](const double t1, const double t2) { return t1 + t2; },
        [](const value_type &t1, const value_type &t2) {
          return scalar_traits::squaredNorm(t1 - t2);
        });

    std::cout << "for fmm:" << std::endl;
    std::cout << "dimension = " << dimension << ". N = " << N
              << ". L2_fmm error = " << L2_fmm << ". L2_fmm relative error is "
              << std::sqrt(L2_fmm / scale)
              << ". time_fmm_setup = " << time_fmm_setup.count()
              << ". time_fmm_eval = " << time_fmm_eval.count() << std::endl;

    if (N == 3) {
      TS_ASSERT_LESS_THAN(L2_fmm / scale, 1e-2);
    }

#ifdef HAVE_EIGEN
    typedef typename ParticlesType::reference reference;
    for (reference p : particles) {
      get<TargetFMM>(p) = scalar_traits::Zero();
    }
    t0 = Clock::now();
    auto fmm_eigen =
        create_fmm_operator<N>(particles, particles, kernel, p2pkernel);
    t1 = Clock::now();
    time_fmm_setup = t1 - t0;

    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_type;
    const int length = scalar_traits::length;
    vector_type target_eigen(length * particles.size());
    vector_type source_eigen(length * particles.size());
    for (size_t i = 0; i < particles.size(); ++i) {
      for (int j = 0; j < length; ++j) {
        source_eigen[i * length + j] =
            scalar_traits::Index(get<Source>(particles)[i], j);
      }
    }
    t0 = Clock::now();
    target_eigen = fmm_eigen * source_eigen;
    t1 = Clock::now();
    time_fmm_eval = t1 - t0;
    for (size_t i = 0; i < particles.size(); ++i) {
      for (int j = 0; j < length; ++j) {
        scalar_traits::Index(get<TargetFMM>(particles)[i], j) =
            target_eigen[i * length + j];
      }
    }

    L2_fmm = std::inner_product(
        std::begin(get<TargetFMM>(particles)),
        std::end(get<TargetFMM>(particles)),
        std::begin(get<TargetManual>(particles)), 0.0,
        [](const double t1, const double t2) { return t1 + t2; },
        [](const value_type &t1, const value_type &t2) {
          return scalar_traits::squaredNorm(t1 - t2);
        });

    std::cout << "for fmm operator:" << std::endl;
    std::cout << "dimension = " << dimension << ". N = " << N
              << ". L2_fmm error = " << L2_fmm << ". L2_fmm relative error is "
              << std::sqrt(L2_fmm / scale)
              << ". time_fmm_setup = " << time_fmm_setup.count()
              << ". time_fmm_eval = " << time_fmm_eval.count() << std::endl;

    if (N == 3) {
      TS_ASSERT_LESS_THAN(L2_fmm / scale, 1e-2);
    }
#endif
  }

  template <unsigned int D, template <typename, typename> class StorageVector,
            template <typename> class SearchMethod>
  void helper_fast_methods(size_t N) {
    typedef Vector<double, D> double_d;
    typedef Vector<int, D> int_d;
    typedef Vector<bool, D> bool_d;
    // randomly generate a bunch of positions over a range
    const double pos_min = 0;
    const double pos_max = 1;
    std::uniform_real_distribution<double> U(pos_min, pos_max);
    generator_type generator(time(NULL));
    auto gen = std::bind(U, generator);
    typedef Vector<double, D> double_d;
    typedef Vector<int, D> int_d;
    const int num_particles_per_bucket = 10;

    typedef Particles<std::tuple<source, target_manual, target_fmm>, D,
                      StorageVector, SearchMethod>
        ParticlesType;
    typedef typename ParticlesType::position position;
    typedef typename ParticlesType::const_reference const_reference;
    ParticlesType particles(N);

    for (size_t i = 0; i < N; i++) {
      for (size_t d = 0; d < D; ++d) {
        get<position>(particles)[i][d] = gen();
        get<source>(particles)[i] = gen();
      }
      get<target_fmm>(particles)[i] = 0.0;
    }
    particles.init_neighbour_search(
        int_d::Constant(pos_min), int_d::Constant(pos_max),
        bool_d::Constant(false), num_particles_per_bucket);

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
    auto p2pkernel = [&c](const_reference pa, const_reference pb) {
      return std::sqrt((get<position>(pb) - get<position>(pa)).squaredNorm() +
                       c);
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

    helper_fast_methods_calculate<1, source, target_manual, target_fmm>(
        particles, kernel, p2pkernel, scale);
    helper_fast_methods_calculate<2, source, target_manual, target_fmm>(
        particles, kernel, p2pkernel, scale);
    helper_fast_methods_calculate<3, source, target_manual, target_fmm>(
        particles, kernel, p2pkernel, scale);

#ifdef HAVE_EIGEN

    if (D == 2) {
      typedef Particles<std::tuple<vsource, vtarget_manual, vtarget_fmm>, D,
                        StorageVector, SearchMethod>
          vParticlesType;
      vParticlesType vparticles(N);

      for (size_t i = 0; i < N; i++) {
        for (size_t d = 0; d < D; ++d) {
          get<position>(vparticles)[i][d] = gen();
        }
      }

      vparticles.init_neighbour_search(
          int_d::Constant(pos_min), int_d::Constant(pos_max),
          bool_d::Constant(false), num_particles_per_bucket);

      // generate a source vector using a smooth cosine
      auto vsource_fn = [&](const double_d &p) {
        // return (p-double_d(0)).norm();
        vdouble2 ret = vdouble2::Constant(1.0);
        const double scale = 2.0 * detail::PI / (pos_max - pos_min);
        for (size_t i = 0; i < D; ++i) {
          ret[0] *= cos((p[i] - pos_min) * scale);
          ret[1] *= sin((p[i] - pos_min) * scale);
        }
        return ret / N;
      };
      std::transform(std::begin(get<position>(vparticles)),
                     std::end(get<position>(vparticles)),
                     std::begin(get<vsource>(vparticles)), vsource_fn);

      const double c = 0.1;
      auto vkernel = [&c](const double_d &pa, const double_d &pb) {
        const double_d x = pb - pa;
        const double r2 = x.squaredNorm();
        const double exp = std::exp(-r2 / std::pow(c, 2));
        Eigen::Matrix<double, 2, 2> ret;
        if (r2 == 0) {
          ret(0, 0) = (-0.5) * exp;
          ret(0, 1) = 0.5 * exp;
          ret(1, 0) = ret(0, 1);
          ret(1, 1) = (-0.5) * exp;
        } else {
          ret(0, 0) = (x[0] * x[0] / r2 - 1) * exp;
          ret(0, 1) = (x[0] * x[1] / r2) * exp;
          ret(1, 0) = ret(0, 1);
          ret(1, 1) = (x[1] * x[1] / r2 - 1) * exp;
        }
        return ret;
      };
      auto vp2pkernel = [&](auto pa, auto pb) {
        return vkernel(get<position>(pa), get<position>(pb));
      };

      // perform the operation manually
      std::fill(std::begin(get<vtarget_manual>(vparticles)),
                std::end(get<vtarget_manual>(vparticles)), vdouble2::Zero());

      t0 = Clock::now();
      for (size_t i = 0; i < N; i++) {
        const double_d pi = get<position>(vparticles)[i];
        for (size_t j = 0; j < N; j++) {
          const double_d pj = get<position>(vparticles)[j];
          get<vtarget_manual>(vparticles)[i] +=
              vkernel(pi, pj) * get<vsource>(vparticles)[j];
        }
      }
      t1 = Clock::now();
      time_manual = t1 - t0;

      const double vscale = std::accumulate(
          std::begin(get<vtarget_manual>(vparticles)),
          std::end(get<vtarget_manual>(vparticles)), 0.0,
          [](const double t1, const vdouble2 t2) { return t1 + t2.dot(t2); });

      std::cout << "VECTOR-VALUED - MANUAL TIMING: dimension = " << D
                << ". number of particles = " << N
                << ". time = " << time_manual.count() << " scale = " << vscale
                << std::endl;

      helper_fast_methods_calculate<1, vsource, vtarget_manual, vtarget_fmm>(
          vparticles, vkernel, vp2pkernel, vscale);
      helper_fast_methods_calculate<2, vsource, vtarget_manual, vtarget_fmm>(
          vparticles, vkernel, vp2pkernel, vscale);
      helper_fast_methods_calculate<3, vsource, vtarget_manual, vtarget_fmm>(
          vparticles, vkernel, vp2pkernel, vscale);
    }
#endif
  }

  template <typename Expansions>
  void helper_fmm_operators(Expansions &expansions) {
    const unsigned int D = Expansions::dimension;
    typedef Vector<double, D> double_d;
    typedef typename Expansions::l_expansion_type l_expansion_type;
    typedef typename Expansions::m_expansion_type m_expansion_type;

    // unit box
    bbox<D> parent(double_d::Constant(0.0), double_d::Constant(1.0));
    bbox<D> leaf1(double_d::Constant(0.0), double_d::Constant(1.0));
    leaf1.bmax[0] = 0.5;
    bbox<D> leaf2(double_d::Constant(0.0), double_d::Constant(1.0));
    leaf2.bmin[0] = 0.5;
    std::cout << "parent = " << parent << " leaf1 = " << leaf1
              << " leaf2 = " << leaf2 << std::endl;

    // create n particles, 2 leaf boxes, 1 parent box
    std::uniform_real_distribution<double> U(0, 1);
    generator_type generator(time(NULL));
    const size_t n = 10;
    double_d particles_in_leaf1[n];
    double_d particles_in_leaf2[n];
    double source_leaf1[n];
    double field_just_self_leaf1[n];
    double field_all_leaf1[n];
    double source_leaf2[n];
    double field_just_self_leaf2[n];
    double field_all_leaf2[n];

    auto f = [](const double_d &p) { return p[0]; };

    for (size_t i = 0; i < n; ++i) {
      particles_in_leaf1[i][0] = 0.5 * U(generator);
      particles_in_leaf2[i][0] = 0.5 * U(generator) + 0.5;
      for (size_t j = 1; j < D; ++j) {
        particles_in_leaf1[i][j] = U(generator);
        particles_in_leaf2[i][j] = U(generator);
      }
      source_leaf1[i] = f(particles_in_leaf1[i]);
      source_leaf2[i] = f(particles_in_leaf2[i]);
    }

    for (size_t i = 0; i < n; ++i) {
      field_just_self_leaf1[i] = 0;
      field_just_self_leaf2[i] = 0;
      for (size_t j = 0; j < n; ++j) {
        field_just_self_leaf1[i] +=
            source_leaf1[j] *
            expansions.m_K(particles_in_leaf1[i], particles_in_leaf1[j]);
        field_just_self_leaf2[i] +=
            source_leaf2[j] *
            expansions.m_K(particles_in_leaf2[i], particles_in_leaf2[j]);
      }
      field_all_leaf1[i] = field_just_self_leaf1[i];
      field_all_leaf2[i] = field_just_self_leaf2[i];
      for (size_t j = 0; j < n; ++j) {
        field_all_leaf1[i] +=
            source_leaf2[j] *
            expansions.m_K(particles_in_leaf1[i], particles_in_leaf2[j]);
        field_all_leaf2[i] +=
            source_leaf1[j] *
            expansions.m_K(particles_in_leaf2[i], particles_in_leaf1[j]);
      }
    }

    // check P2M, and L2P
    m_expansion_type expansionM_leaf1 = {0};

    for (size_t i = 0; i < n; ++i) {
      expansions.P2M(expansionM_leaf1, leaf1, particles_in_leaf1[i],
                     source_leaf1[i]);
    }

    m_expansion_type expansionL_leaf1 = {0};
    expansions.M2L(expansionL_leaf1, leaf1, leaf1, expansionM_leaf1);

    double L2 = 0;
    double scale = 0;
    for (size_t i = 0; i < n; ++i) {
      const double check =
          expansions.L2P(particles_in_leaf1[i], leaf1, expansionL_leaf1);
      L2 += std::pow(check - field_just_self_leaf1[i], 2);
      scale += std::pow(field_just_self_leaf1[i], 2);
      TS_ASSERT_LESS_THAN(std::abs(check - field_just_self_leaf1[i]), 2e-4);
    }

    TS_ASSERT_LESS_THAN(std::sqrt(L2 / scale), 1e-4);

    m_expansion_type expansionM_leaf2 = {};
    for (size_t i = 0; i < n; ++i) {
      expansions.P2M(expansionM_leaf2, leaf2, particles_in_leaf2[i],
                     source_leaf2[i]);
    }

    l_expansion_type expansionL_leaf2 = {};
    expansions.M2L(expansionL_leaf2, leaf2, leaf2, expansionM_leaf2);

    L2 = 0;
    for (size_t i = 0; i < n; ++i) {
      const double check =
          expansions.L2P(particles_in_leaf2[i], leaf2, expansionL_leaf2);
      L2 += std::pow(check - field_just_self_leaf2[i], 2);
      scale += std::pow(field_just_self_leaf2[i], 2);
    }
    TS_ASSERT_LESS_THAN(std::sqrt(L2 / scale), 1e-4);

    // check M2M and L2L
    m_expansion_type expansionM_parent = {};
    expansions.M2M(expansionM_parent, parent, leaf1, expansionM_leaf1);
    expansions.M2M(expansionM_parent, parent, leaf2, expansionM_leaf2);
    l_expansion_type expansionL_parent = {};
    expansions.M2L(expansionL_parent, parent, parent, expansionM_parent);

    l_expansion_type reexpansionL_leaf1 = {};
    expansions.L2L(reexpansionL_leaf1, leaf1, parent, expansionL_parent);

    L2 = 0;
    scale = 0;
    for (size_t i = 0; i < n; ++i) {
      const double check =
          expansions.L2P(particles_in_leaf1[i], leaf1, reexpansionL_leaf1);
      L2 += std::pow(check - field_all_leaf1[i], 2);
      scale += std::pow(field_all_leaf1[i], 2);
    }
    TS_ASSERT_LESS_THAN(std::sqrt(L2 / scale), 1e-4);
  }

  void test_fmm_operators() {
    const unsigned int D = 2;
    typedef Vector<double, D> double_d;
    auto kernel = [](const double_d &pa, const double_d &pb) {
      return std::sqrt((pb - pa).squaredNorm() + 0.1);
    };
    detail::BlackBoxExpansions<D, 10, decltype(kernel), 1, 1> expansions(
        kernel);
    helper_fmm_operators(expansions);
  }

  void test_fast_methods_CellList(void) {
    const size_t N = 1000;
#ifdef HAVE_GPERFTOOLS
    ProfilerStart("fmm_CellList");
#endif
    std::cout << "BUCKET_SEARCH_SERIAL: testing 1D..." << std::endl;
    helper_fast_methods<1, std::vector, CellList>(N);
    std::cout << "BUCKET_SEARCH_SERIAL: testing 2D..." << std::endl;
    helper_fast_methods<2, std::vector, CellList>(N);
    std::cout << "BUCKET_SEARCH_SERIAL: testing 3D..." << std::endl;
    helper_fast_methods<3, std::vector, CellList>(N);
#ifdef HAVE_GPERFTOOLS
    ProfilerStop();
#endif
  }

  void test_fast_methods_CellListOrdered(void) {
    const size_t N = 1000;
#ifdef HAVE_GPERFTOOLS
    ProfilerStart("fmm_CellListOrdered");
#endif
    std::cout << "BUCKET_SEARCH_PARALLEL: testing 1D..." << std::endl;
    helper_fast_methods<1, std::vector, CellListOrdered>(N);
    std::cout << "BUCKET_SEARCH_PARALLEL: testing 2D..." << std::endl;
    helper_fast_methods<2, std::vector, CellListOrdered>(N);
    std::cout << "BUCKET_SEARCH_PARALLEL: testing 3D..." << std::endl;
    helper_fast_methods<3, std::vector, CellListOrdered>(N);
#ifdef HAVE_GPERFTOOLS
    ProfilerStop();
#endif
  }

  void test_fast_methods_kd_tree(void) {
    const size_t N = 1000;
#ifdef HAVE_GPERFTOOLS
    ProfilerStart("fmm_kd_tree");
#endif
    std::cout << "KD_TREE: testing 1D..." << std::endl;
    helper_fast_methods<1, std::vector, Kdtree>(N);
    std::cout << "KD_TREE: testing 2D..." << std::endl;
    helper_fast_methods<2, std::vector, Kdtree>(N);
    std::cout << "KD_TREE: testing 3D..." << std::endl;
    helper_fast_methods<3, std::vector, Kdtree>(N);
#ifdef HAVE_GPERFTOOLS
    ProfilerStop();
#endif
  }

  void test_fast_methods_kd_tree_nanoflann(void) {
    const size_t N = 1000;
#ifdef HAVE_GPERFTOOLS
    ProfilerStart("fmm_kd_tree");
#endif
    std::cout << "KD_TREE: testing 1D..." << std::endl;
    helper_fast_methods<1, std::vector, KdtreeNanoflann>(N);
    std::cout << "KD_TREE: testing 2D..." << std::endl;
    helper_fast_methods<2, std::vector, KdtreeNanoflann>(N);
    std::cout << "KD_TREE: testing 3D..." << std::endl;
    helper_fast_methods<3, std::vector, KdtreeNanoflann>(N);
#ifdef HAVE_GPERFTOOLS
    ProfilerStop();
#endif
  }

  void test_fast_methods_HyperOctree(void) {
    const size_t N = 1000;
#ifdef HAVE_GPERFTOOLS
    ProfilerStart("fmm_oct_tree");
#endif
    std::cout << "OCTTREE: testing 1D..." << std::endl;
    helper_fast_methods<1, std::vector, HyperOctree>(N);
    std::cout << "OCTTREE: testing 2D..." << std::endl;
    helper_fast_methods<2, std::vector, HyperOctree>(N);
    std::cout << "OCTTREE: testing 3D..." << std::endl;
    helper_fast_methods<3, std::vector, HyperOctree>(N);
#ifdef HAVE_GPERFTOOLS
    ProfilerStop();
#endif
  }
};

#endif
