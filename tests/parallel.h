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

#ifndef PARALLEL_H_
#define PARALLEL_H_

#include <chrono>
#include <cxxtest/TestSuite.h>
typedef std::chrono::system_clock Clock;

#include "Aboria.h"

using namespace Aboria;

class ParallelTest : public CxxTest::TestSuite {
public:
  ABORIA_VARIABLE(thrust_neighbour_count, int, "thrust_neighbour_count")

  void test_documentation(void) {
    //[parallel
    /*`
    [section Parallelism in Aboria]

    Aboria can use OpenMP or CUDA to utilise multiple cores or Nvidia GPUs that
    you might have. In general, Aboria uses the parallel algorithms and vectors
    provided with [@http://thrust.github.io Thrust] to do this. From there, you
    can either use Thrust's OpenMP or CUDA backends to provide the type of
    parallelism you wish. However, there are a few parts of Aboria that are
    OpenMP only (notably the entirety of the [link aboria.symbolic_expressions
    symbolic] and [link aboria.evaluating_and_solving_kernel_op kernel] APIs).

    [section OpenMP]

    You don't have to do much to start using OpenMP for the high level symbolic
    or kernel interfaces, all that is required is that you install OpenMP and
    Thrust and add the `HAVE_THRUST` compiler definition (see [link
    aboria.installation_and_getting_started]). For lower-level programming you
    will need to add a few pragmas to your code, a few examples of which are
    discussed below.

    First, let's create a set of `N` particles as usual

    */
    const size_t N = 100;
    ABORIA_VARIABLE(neighbour_count, int, "neighbour_count");
    typedef Particles<std::tuple<neighbour_count>, 2> particle_t;
    typedef particle_t::position position;
    particle_t particles(N);

    /*`
    Now we will loop through the particles and set their initial positions
    randomly. In order to use an OpenMP parallel for loop, we will stick to a
    simple index-based loop for loop to iterate through the particles, like so
    */

#pragma omp parallel for
    for (size_t i = 0; i < particles.size(); ++i) {
      std::uniform_real_distribution<double> uniform(0, 1);
      auto gen = get<generator>(particles)[i];
      get<position>(particles)[i] = vdouble2(uniform(gen), uniform(gen));
    }

    /*`
    Now we can initialise the neighbourhood search data structure. Note that all
    creation and updates to the spatial data structures are run in parallel.

    [note currently the only data structure that is created or updated in serial
    is [classref Aboria::KdtreeNanoflann]. All the rest are done in parallel
    using either OpenMP or CUDA]
    */

    particles.init_neighbour_search(
        vdouble2::Constant(0), vdouble2::Constant(1), vbool2::Constant(false));

    /*`

    We will use Aboria's range search to look for neighbouring pairs within a
    cutoff, and once again use OpenMPs parallel loop. All queries to the spatial
    data structures are thread-safe and can be used in parallel.
    */
    const double radius = 0.1;

#pragma omp parallel for
    for (size_t i = 0; i < particles.size(); ++i) {
      for (auto j = euclidean_search(particles.get_query(),
                                     get<position>(particles)[i], radius);
           j != false; ++j) {
        ++get<neighbour_count>(particles)[i];
      }
    }

    /*`
    In general, that is 90% of what you need to know, just add a couple of
    OpenMP pragmas to your loops and you are ready to go!

    [endsect]

    */
//<-
#ifdef HAVE_THRUST
    //->
    /*`

    [section CUDA]

    [caution CUDA support in Aboria is experimental, and is not tested
    regularly. We welcome feedback by any CUDA users if anything doesn't
    work for you]

    Writing CUDA compatible code is slightly more involved. Aboria uses the
    Thrust library for CUDA parallism, and follows similar patterns (i.e.
    STL-like).

    Most importantly, we need to make sure that all the particle data is
    contained in vectors that are stored on the GPU. To do this we use a
    `thrust::device_vector` as the base storage vector for our particles
    class

    [note we want to use the type `thrust_neighbour_count` within a device
    function, so we need to define this type outside any host functions
    (including main).]

    */

    //=ABORIA_VARIABLE(thrust_neighbour_count, int, "thrust_neighbour_count");
    //=
    //=int main() {

    typedef Particles<std::tuple<thrust_neighbour_count>, 2,
                      thrust::device_vector, CellListOrdered>
        thrust_particle_t;
    thrust_particle_t thrust_particles(N);

    /*`

    Since all our data is on the device, we cannot use raw for loops to access
    this data without copying it back to the host, an expensive operation.
    Instead, Thrust provides a wide variety of parallel algorithms to manipulate
    the data. Aboria's [classref Aboria::zip_iterator] is compatible with the
    Thrust framework, so can be used in a similar fashion to Thrust's own
    `zip_iterator` (except, unlike Thrust's `zip_iterator`, we can take
    advantage of Aboria's tagged `reference` and `value_types`).

    We can use Thrust's `tabulate` algorithm to loop through the particles and
    set their initial positions randomly.

    */

    thrust::tabulate(get<position>(thrust_particles).begin(),
                     get<position>(thrust_particles).end(),
                     [] __device__(const int i) {
                       thrust::default_random_engine gen;
                       thrust::uniform_real_distribution<float> uni(0, 1);
                       gen.discard(i);
                       return vdouble2(uni(gen), uni(gen));
                     });

    /*`
    Now we can initialise the neighbourhood search data structure. Note that we
    are using [classref Aboria::CellListOrdered] data structure, which is
    similar to [classref Aboria::CellList] but instead relies on reordering the
    particles to arrange them into cells, which is more amenable to
    parallelisation using a GPU.
    */

    thrust_particles.init_neighbour_search(
        vdouble2::Constant(0), vdouble2::Constant(1), vbool2::Constant(false));

    /*`

    We can use any of Aboria's range searches within a Thrust algorithm. Below
    we will implement a range search around each particle, counting all
    neighbours within range. Note that we need to copy all of the variables from
    the outer scope to the lambda function, since the lambda will run on the
    device, and won't be able to access any host memory.

    [note The [classref Aboria::NeighbourQueryBase] class for each spatial data
    structure is designed to be copyable to the GPU, but the [classref
    Aboria::Particles] class is not, so while the `query` variable is copyable
    to the device, the `thrust_particles` variable is not.]

    [note The type of variable `i` in the lambda will be deduced as [classref
    Aboria::Particles::raw_reference]. This is different to [classref
    Aboria::Particles::reference] when using `thrust::device_vector`, but acts
    in a similar fashion]
    */

    thrust::for_each(
        thrust_particles.begin(), thrust_particles.end(),
        [radius, query = thrust_particles.get_query()] __device__(auto i) {
          get<thrust_neighbour_count>(i) = 0;
          for (auto j = euclidean_search(query, get<position>(i), radius);
               j != false; ++j) {
            ++get<thrust_neighbour_count>(i);
          }
        });

    /*`

    While we have exclusively used `thrust::for_each` above, the iterators that
    Aboria provides for the [classref Aboria::Particles] container should work
    with all of Thrust's algorithms. For example, you might wish to restructure
    the previous code as a transform:
     */

    thrust::transform(
        thrust_particles.begin(), thrust_particles.end(),
        get<thrust_neighbour_count>(thrust_particles).begin(),
        [radius, query = thrust_particles.get_query()] __device__(auto i) {
          int sum = 0;
          for (auto j = euclidean_search(query, get<position>(i), radius);
               j != false; ++j) {
            ++sum;
          }
          return sum;
        });
//=}

/*`
[endsect]
*/
//<-
#endif // HAVE_THRUST
       //->
       /*`
       [endsect]
       */
    //]
  }
};

#endif /* PARALLEL_H_ */
