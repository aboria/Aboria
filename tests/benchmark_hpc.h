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

#ifndef BENCHMARK_HPC_H_
#define BENCHMARK_HPC_H_

#include "Aboria.h"
#include <chrono>
#include <cxxtest/TestSuite.h>
typedef std::chrono::system_clock Clock;
#include <fstream> // std::ofstream
#include <thread>
#ifdef HAVE_GPERFTOOLS
#include <gperftools/profiler.h>
#endif

using namespace Aboria;

class BenchmarkHPC : public CxxTest::TestSuite {
public:
  template <template <typename, typename> class VectorT,
            template <typename> class SearchMethod, unsigned int D,
            typename VectorT2>
  std::tuple<double, double> md_step(const VectorT2 &position_vect,
                                     const size_t repeats) {
    typedef Vector<double, D> double_d;
    typedef Vector<bool, D> bool_d;
    const size_t n = position_vect.size();
    std::cout << "md_step: D = " << D << " N = " << N << " n = " << n
              << " repeats = " << repeats << std::endl;
    typedef Particles<std::tuple<>, D, VectorT, SearchMethod> particles_t;
    typedef typename particles_t::position position;
    particles_t particles(n);

    double_d min = double_d::Constant(-0.1);
    double_d max = double_d::Constant(1.1);
    const double diameter = 0.1;
    bool_d periodic = bool_d::Constant(false);

    auto t0 = Clock::now();
    detail::copy(position_vect.begin(), position_vect.end(),
                 get<position>(particles).begin());

    particles.init_neighbour_search(min, max, periodic);
    auto query = particles.get_query();
    auto t1 = Clock::now();
    std::chrono::duration<double> time_setup = t1 - t0;

    const double dt = 1.0;
    const double mass = 1.0;
    const double k = 1.0;
    auto kernel = [=] CUDA_HOST_DEVICE(typename particles_t::raw_reference i) {
      double_d force_sum = double_d::Constant(0);
      for (auto j = euclidean_search(query, get<position>(i), diameter);
           j != false; ++j) {
        const double r = j.dx().norm();
        if (r != 0) {
          force_sum -= k * (diameter / r - 1) * j.dx();
        }
      }
      get<position>(i) += dt * dt * force_sum / mass;
    };

    detail::for_each(particles.begin(), particles.end(), kernel);

    t0 = Clock::now();
#ifdef HAVE_GPERFTOOLS
    typedef typename particles_type::traits_type traits_type;
    if (std::is_same<SearchMethod<traits_type>,
                     Kdtree<traits_type>>::value) {
      ProfilerStart(
          ("md_step_h2_kdtree" + std::to_string(N) + "_" + std::to_string(D))
              .c_str());
    } else if (std::is_same<SearchMethod<traits_type>,
                            HyperOctree<traits_type>>::value) {

      ProfilerStart(
          ("md_step_octtree" + std::to_string(N) + "_" + std::to_string(D))
              .c_str());
    } else if (std::is_same<SearchMethod<traits_type>,
                            CellList<traits_type>>::value) {

      ProfilerStart(
          ("md_step_cell_list" + std::to_string(N) + "_" + std::to_string(D))
              .c_str());
    } else {
      ProfilerStart(("md_step_cell_list_ordered" + std::to_string(N) + "_" +
                     std::to_string(D))
                        .c_str());
    }

#endif
    for (size_t ii = 0; ii < repeats; ++ii) {
      detail::for_each(particles.begin(), particles.end(), kernel);
    }
#ifdef HAVE_GPERFTOOLS
    ProfilerStop();
#endif
    t1 = Clock::now();
    std::chrono::duration<double> time_step = t1 - t0;
    std::cout << "time setup = " << time_setup.count() << std::endl;
    std::cout << "time step = " << time_step.count() / repeats << std::endl;

    return std::make_tuple(time_setup.count(), time_step.count() / repeats);
  }

  template <unsigned int D> void helper_md_step() {
    std::ofstream file_setup, file_step;
    // const size_t base_repeatsN2 = 1e7;
    const size_t base_repeatsN = 2e5;
    const size_t ncpus = 2;
    file_setup.open("benchmark_md_setup" + std::to_string(D) + ".csv");
    file_step.open("benchmark_md_step" + std::to_string(D) + ".csv");
    file_setup << "#" << std::setw(25) << "N";
    file_step << "#" << std::setw(25) << "N";

    for (size_t i = 1; i <= ncpus; ++i) {
      file_setup << std::setw(25)
                 << "OPENMP_cell_list_ordered" + std::to_string(i);
      file_step << std::setw(25)
                << "OPENMP_cell_list_ordered" + std::to_string(i);

      file_step << std::setw(25) << "OPENMP_octtree" + std::to_string(i);
      file_setup << std::setw(25) << "OPENMP_octtree" + std::to_string(i);

#if not defined(__CUDACC__)
      file_setup << std::setw(25) << "OPENMP_cell_list" + std::to_string(i);
      file_step << std::setw(25) << "OPENMP_cell_list" + std::to_string(i);
#endif
    }

    file_setup << std::endl;
    file_step << std::endl;

    for (double i = 1000; i < 1000000; i *= 1.1) {
      //{
      //    double i = 10000;
      const size_t N = i;
      const int seed = 0;
      // randomly generate a bunch of positions over a range
      const double pos_min = 0;
      const double pos_max = 1;
      typedef Vector<double, D> double_d;

      std::vector<double_d> positions(N);

      detail::counting_iterator<int> start(0);
      detail::transform(
          start, start + N, positions.begin(), [=](const int idx) {
            generator_type generator(seed);
            detail::uniform_real_distribution<double> U(pos_min, pos_max);
            generator.discard(idx * D);
            double_d p;
            for (size_t d = 0; d < D; ++d) {
              p[d] = U(generator);
            }
            return p;
          });

      // const size_t repeatsN2 = base_repeatsN2 / std::pow(N, 2) + 1;
      const size_t repeatsN = base_repeatsN / N + 1;
      file_setup << std::setw(15) << N;
      file_step << std::setw(15) << N;

      for (size_t i = 2; i <= ncpus; ++i) {
#ifdef HAVE_OPENMP
        omp_set_num_threads(i);
#endif
        double time_setup;
        double time_step;

        std::tie(time_setup, time_step) =
            md_step<thrust::device_vector, CellListOrdered, D>(positions,
                                                               repeatsN);
        file_setup << std::setw(25) << time_setup;
        file_step << std::setw(25) << time_step;

        std::tie(time_setup, time_step) =
            md_step<thrust::device_vector, HyperOctree, D>(positions, repeatsN);
        file_setup << std::setw(25) << time_setup;
        file_step << std::setw(25) << time_step;

#if not defined(__CUDACC__)
        std::tie(time_setup, time_step) =
            md_step<thrust::device_vector, CellList, D>(positions, repeatsN);
        file_setup << std::setw(25) << time_setup;
        file_step << std::setw(25) << time_step;
#endif
      }

      file_setup << std::endl;
      file_step << std::endl;
    }
    file_setup.close();
    file_step.close();
  }

  void test_multiquadric() {
    // helper_multiquadric<1>();
    helper_md_step<2>();
    // helper_multiquadric<3>();
    // helper_multiquadric<4>();
  }
};

#endif