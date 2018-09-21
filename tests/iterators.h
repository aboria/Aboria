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

#ifndef ITERATORS_TESTS_H_
#define ITERATORS_TESTS_H_

#include <cxxtest/TestSuite.h>

#include "Aboria.h"

using namespace Aboria;

class IteratorsTest : public CxxTest::TestSuite {
public:
  typedef std::mt19937 generator_type;
  generator_type generator;

  void test_lattice_within_distance(void) {
    using particles_t = Particles<std::tuple<>, 1, std::vector, CellList>;
    using position = typename particles_t::position;
    particles_t particles(100);
    std::uniform_real_distribution<double> uni(0, 1);
    for (auto &p : get<position>(particles)) {
      p[0] = uni(generator);
    }
    particles.init_neighbour_search(
        vdouble1::Constant(0), vdouble1::Constant(1), vbool1::Constant(false));

    //   0    1    2    3    4    5    6     7    8   9
    // |----|----|----|----|----|----|----|----|----|----|
    // 0   0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9  1.0

    int count = 0;
    for (auto p = particles.get_query().get_buckets_near_point<2>(
             vdouble1::Constant(0.5), 0.05);
         p != false; ++p) {
      ++count;
      std::cout << *p;
    }
    std::cout << std::endl;
    TS_ASSERT_EQUALS(count, 2);

    count = 0;
    for (auto p = particles.get_query().get_buckets_near_point<2>(
             vdouble1::Constant(0.5), 0.15);
         p != false; ++p) {
      ++count;
      std::cout << *p;
    }
    std::cout << std::endl;
    TS_ASSERT_EQUALS(count, 4);

    count = 0;
    for (auto p = particles.get_query().get_buckets_near_point<2>(
             vdouble1::Constant(0.1), 0.15);
         p != false; ++p) {
      ++count;
      std::cout << *p;
    }
    std::cout << std::endl;
    TS_ASSERT_EQUALS(count, 3);

    count = 0;
    for (auto p = particles.get_query().get_buckets_near_point<2>(
             vdouble1::Constant(-0.1), 0.05);
         p != false; ++p) {
      ++count;
      std::cout << *p;
    }
    std::cout << std::endl;
    TS_ASSERT_EQUALS(count, 0);

    count = 0;
    for (auto p = particles.get_query().get_buckets_near_point<2>(
             vdouble1::Constant(-0.1), 0.15);
         p != false; ++p) {
      ++count;
      std::cout << *p;
    }
    std::cout << std::endl;
    TS_ASSERT_EQUALS(count, 1);
  }

  void test_lattice_within_distance_2d(void) {
    using particles_t = Particles<std::tuple<>, 2, std::vector, CellList>;
    using position = typename particles_t::position;
    particles_t particles(1000);
    std::uniform_real_distribution<double> uni(0, 1);
    for (auto &p : get<position>(particles)) {
      p = vdouble2(uni(generator), uni(generator));
    }
    particles.init_neighbour_search(vdouble2(0, 0), vdouble2(1, 1),
                                    vbool2(false, false));

    //       0    1    2    3    4    5    6     7    8   9
    // 0.0 |----|----|----|----|----|----|----|----|----|----|
    // 0.1 |----|----|----|----|----|----|----|----|----|----|
    // 0.2 |----|----|----|----|----|----|----|----|----|----|
    // 0.3 |----|----|----|----|----|----|----|----|----|----|
    // 0.4 |----|----|----|----|----|----|----|----|----|----|
    // 0.5 |----|----|----|----|----|----|----|----|----|----|
    // 0.6 |----|----|----|----|----|----|----|----|----|----|
    // 0.7 |----|----|----|----|----|----|----|----|----|----|
    // 0.8 |----|----|----|----|----|----|----|----|----|----|
    // 0.9 |----|----|----|----|----|----|----|----|----|----|
    //     0   0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9  1.0

    int count = 0;
    for (auto p = particles.get_query().get_buckets_near_point<2>(
             vdouble2(0.5, 0.5), 0.05);
         p != false; ++p) {
      ++count;
      std::cout << *p;
    }
    std::cout << std::endl;
    TS_ASSERT_EQUALS(count, 4);

    count = 0;
    for (auto p = particles.get_query().get_buckets_near_point<2>(
             vdouble2(0.5, 0.5), 0.1001);
         p != false; ++p) {
      ++count;
      std::cout << *p;
    }
    std::cout << std::endl;
    TS_ASSERT_EQUALS(count, 12);

    count = 0;
    for (auto p = particles.get_query().get_buckets_near_point<2>(
             vdouble2(0.55, 0.55), 0.1001);
         p != false; ++p) {
      ++count;
      std::cout << *p;
    }
    std::cout << std::endl;
    TS_ASSERT_EQUALS(count, 9);

    count = 0;
    for (auto p = particles.get_query().get_buckets_near_point<2>(
             vdouble2(0.55, 0.55), 0.049999);
         p != false; ++p) {
      ++count;
      std::cout << *p;
    }
    std::cout << std::endl;
    TS_ASSERT_EQUALS(count, 1);

    count = 0;
    for (auto p = particles.get_query().get_buckets_near_point<2>(
             vdouble2(-0.001, 0.001), 2.0);
         p != false; ++p) {
      ++count;
      std::cout << *p;
    }
    std::cout << std::endl;
    TS_ASSERT_EQUALS(count, 100);

    count = 0;
    for (auto p = particles.get_query().get_buckets_near_point<2>(
             vdouble2(-0.001, 0.001), 0.01);
         p != false; ++p) {
      ++count;
      std::cout << *p;
    }
    std::cout << std::endl;
    TS_ASSERT_EQUALS(count, 1);

    count = 0;
    for (auto p = particles.get_query().get_buckets_near_point<2>(
             vdouble2(1.001, 1.001), 0.01);
         p != false; ++p) {
      ++count;
      std::cout << *p;
    }
    std::cout << std::endl;
    TS_ASSERT_EQUALS(count, 1);
  }
};

#endif
