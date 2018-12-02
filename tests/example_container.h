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

#ifndef EXAMPLE_CONTAINER_H_
#define EXAMPLE_CONTAINER_H_

#include <cxxtest/TestSuite.h>

#include <random>

#include "Aboria.h"
using namespace Aboria;

class ExampleContainerTest : public CxxTest::TestSuite {
public:
  typedef std::mt19937 generator_type;
  generator_type generator;

  template <template <typename> class SearchMethod> void helper_example(void) {
    //[example_neighbour_search
    /*`
    This example creates $N$ randomly distributed particles in a periodic
    domain. For each particle, the number of neighbouring particles within a
    certain radius is counted.

    */
    //<-
    // clang-format off
    //->
//=#include <random>
//=
//=#include "Aboria.h"
//=using namespace Aboria;
//=
//=int main() {
    //<-
    // clang-format on
    //->
    // Create a 2d particle container with one additional variable which
    // will hold the count for each particle
    ABORIA_VARIABLE(count, int, "number_of_neighbours")
    //<-
    typedef Particles<std::tuple<count>, 2, std::vector, SearchMethod>
        container_type;
    //->
    //=typedef Particles<std::tuple<count>,2> container_type;

    typedef typename container_type::position position;

    // set radius and number of particles
    const double radius = 0.5;
    const int N = 100;

    // create the $N$ particles
    container_type particles(N);

    // create N particles uniformly distributed
    std::uniform_real_distribution<double> uni(0, 1);
    for (int i = 0; i < N; ++i) {
      get<position>(particles)[i] = vdouble2(uni(generator), uni(generator));
    }

    // initiate neighbour search on a periodic 2d domain of side length L
    // set average number of particles per cell to 1
    particles.init_neighbour_search(vdouble2(0, 0), vdouble2(1, 1),
                                    vbool2(true, true));

    // count neighbours around each particle (note: includes self!) and store in
    // `count` variable
    for (int i = 0; i < N; ++i) {
      get<count>(particles)[i] = 0;
      for (auto j = euclidean_search(particles.get_query(),
                                     get<position>(particles)[i], radius);
           j != false; ++j) {
        get<count>(particles)[i]++;
      }
    }
  }
  //]

  void test_CellListOrdered() { helper_example<CellListOrdered>(); }

  void test_CellList() { helper_example<CellList>(); }

  void test_Kdtree() { helper_example<Kdtree>(); }

  void test_KdtreeNanoflann() { helper_example<KdtreeNanoflann>(); }

  void test_HyperOctree() { helper_example<HyperOctree>(); }
};

#endif /* MD_TEST_H_ */
