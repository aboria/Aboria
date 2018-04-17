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
    symbolic] and [link aboria evaluating_and_solving_kernel_op kernel] APIs ,
    the [link
    aboria.evaluating_and_solving_kernel_op.creating_fast_multipole_method_o
    fast multipole method] and [link
    aboria.evaluating_and_solving_kernel_op.creating_hierarchical_matrix_ope H2
    matrices]).

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
    typedef Particles<std::tuple<neighbour_count>> particle_t;
    typedef particle_t::position position;
    particle_type particles(N);

    /*`
    Now we will loop through the particles and set their initial positions
    randomly. We will use an index-based loop for this purpose so use can use an
    OpenMP parallel loop, like so
    */

#pragma omp parallel for
    for (int i = 0; i < particles.size(); ++i) {
      std::uniform_real_distribution<double> uniform(0, 1);
      auto gen = get<generator>(particles)[i];
      get<position>(particles)[i] = vdouble2(uniform(gen), uniform(gen));
    }

    /*`
    Now we can initialise the neighbourhood search data structure. Note that all
    creation and updates to the spatial data structures are run in parallel.
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
    for (int i = 0; i < particles.size(); ++i) {
      for (auto i = euclidean_search(particles.get_query(),
                                     get<position>(particles)[i], radius);
           i != false; ++i) {
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

    Firstly, we need to make sure that all the particle data is contained in
    vectors that are stored on the GPU. To do this we use a
    `thrust::device_vector` as the base storage vector for our particles
    class

    */

    typedef Particles<std::tuple<neighbour_count>, 2, thrust::device_vector,
                      CellListOrdered>
        thrust_particle_t;
    thrust_particle_t thrust_particles(N);

    /*`
    Since all our data is on the device, we cannot use raw for loops to access
    this data without copying it back to the host, an expensive operation.
    Instead, Thrust provides a wide variety of parallel algorithms to manipulate
    the data. Aboria's zip_iterators are compatible with the Thrust framework,
    so can be used in a similar fashion to Thrust's own zip_iterators (except,
    unlike Thrust's `zip_iterator`, we can take advantage of Aboria's tagged
    `reference` and `value_types`).

    We can use Thrust's `for_each` algorithm to loop through the particles and
    set their initial positions randomly.
    */

    thrust::for_each(thrust_particles.begin(), thrust_particles.end(),
                     [] __device__(auto i) {
                       /*
                        * set a random position, and initialise velocity
                        */
                       auto gen = get<generator>(i);
                       thrust::uniform_real_distribution<float> uni(0, 1);
                       get<thrust_position>(i) = vdouble2(uni(gen), uni(gen));
                     });

    /*`
    Now we can initialise the neighbourhood search data structure. Note that we
    are using Aboria's `CellListOrdered` data structure, which is similar to
    `CellList` but instead relies on reordering the particles to arrange them
    into cells, which is more amenable to parallelisation using a GPU.
    */

    thrust_particles.init_neighbour_search(
        vdouble2::Constant(0), vdouble2::Constant(1), vbool2::Constant(false));

    /*`

    We can use any of Aboria's range searches within a Thrust algorithm. Once
    again we will implement a range search around each particle, counting all
    neighbours within range. Note that we need to copy all of the variables from
    the outer scope to the lambda function, since the lambda will run on the
    device, and won't be able to access any host memory.

    [note The NeighbourQueryBase class for each spatial data structure is
    designed to be copyable to the GPU, but the Particles class is not, so while
    the `query` variable is copyable to the device, the `thrust_particles`
    variable is not.]

    [note The type of variable `i` in the lambda will be deduced as
    Particles::raw_reference. This is different to Particles::reference when
    using `thrust::device_vector`, but acts in a similar fashion]
    */

    thrust::for_each(
        thrust_particles.begin(), thrust_particles.end(),
        [radius, query = thrust_particles.get_query()] __device__(auto i) {
          get<neighbour_count>(i) = 0;
          for (auto j = euclidean_search(query, get<position>(i), radius);
               j != false; ++j) {
            ++get<neighbour_count>(i);
          }
        });

    /*`

    While we have exclusively used `thrust::for_each` above, the iterators that
    Aboria provides for the Particles container should work with all of Thrust's
    algorithms. For example, you might wish to restructure the previous code as
    a transform:
     */

    thrust::transform(
        thrust_particles.begin(), thrust_particles.end(),
        get<neighbour_count>(thrust_particles).begin(),
        [radius, query = thrust_particles.get_query()] __device__(auto i) {
          int sum =
              0 for (auto j = euclidean_search(query, get<position>(i), radius);
                     j != false; ++j) {
            ++sum;
          }
          return sum;
        });

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

  template <unsigned int D, typename Reference> struct set_random {
    double a, b;
    int N;

    set_random(double a, double b, int N) : a(a), b(b), N(N) {}

    CUDA_HOST_DEVICE
    void operator()(Reference arg) {
      Vector<double, D> &p = get<position_d<D>>(arg);
      generator_type &gen = get<generator>(arg);

#if defined(__CUDACC__)
      thrust::uniform_real_distribution<float> dist(a, b);
      thrust::uniform_int_distribution<int> dist_int(0, N - 1);
#else
      std::uniform_real_distribution<float> dist(a, b);
      std::uniform_int_distribution<int> dist_int(0, N - 1);
#endif
      for (size_t d = 0; d < D; ++d) {
        get<id_to_find>(arg) = dist_int(gen);
        p[d] = dist(gen);
      }
    }
  };

  template <typename ParticlesType> struct brute_force_check {
    typedef typename ParticlesType::raw_reference reference;
    typedef typename ParticlesType::raw_pointer pointer;
    typedef typename ParticlesType::raw_const_reference const_reference;
    typedef typename ParticlesType::double_d double_d;
    typedef typename ParticlesType::int_d int_d;
    typedef typename ParticlesType::query_type query_type;
    typedef typename ParticlesType::position position;
    static const unsigned int D = ParticlesType::dimension;

    query_type query;

    brute_force_check(ParticlesType &particles)
        : query(particles.get_query()) {}

    ABORIA_HOST_DEVICE_IGNORE_WARN
    void operator()(reference i) {
      pointer begin = query.get_particles_begin();
      get<index_brute>(i) = -1;
      for (size_t ii = 0; ii < query.number_of_particles(); ++ii) {
        const_reference j = *(begin + ii);
        if (get<id>(j) == get<id_to_find>(i)) {
          get<index_brute>(i) = ii;
        }
      }
    }
  };
  template <typename ParticlesType> struct aboria_check {
    typedef typename ParticlesType::raw_reference reference;
    typedef typename ParticlesType::raw_pointer pointer;
    typedef typename ParticlesType::raw_const_reference const_reference;
    typedef typename ParticlesType::double_d double_d;
    typedef typename ParticlesType::query_type query_type;
    typedef typename ParticlesType::position position;
    static const unsigned int D = ParticlesType::dimension;

    query_type query;

    aboria_check(ParticlesType &particles) : query(particles.get_query()) {}

    ABORIA_HOST_DEVICE_IGNORE_WARN
    void operator()(reference i) {
      auto result = query.find(get<id_to_find>(i));
      if (result == query.get_particles_begin() + query.number_of_particles()) {
        get<index_aboria>(i) = -1;
      } else {
        get<index_aboria>(i) = result - query.get_particles_begin();
      }
    }
  };

  template <unsigned int D, template <typename, typename> class VectorType,
            template <typename> class SearchMethod>
  void helper_d_random(const int N, const bool push_back_construction,
                       const bool set_domain) {
    typedef Particles<std::tuple<id_to_find, index_brute, index_aboria>, D,
                      VectorType, SearchMethod>
        particles_type;
    typedef position_d<D> position;
    typedef Vector<double, D> double_d;
    typedef Vector<bool, D> bool_d;
    double_d min = double_d::Constant(-1);
    double_d max = double_d::Constant(1);
    bool is_periodic = false;
    bool_d periodic = bool_d::Constant(is_periodic);
    particles_type particles;

    std::cout << "random test (D=" << D << " periodic= " << is_periodic
              << "  N=" << N
              << " push_back_construction = " << push_back_construction
              << ", set domain = " << set_domain << "):" << std::endl;

    unsigned seed1 =
        std::chrono::system_clock::now().time_since_epoch().count();
    seed1 = 4238735308;
    std::cout << "seed is " << seed1 << std::endl;
    particles.set_seed(seed1);
    generator_type gen(seed1);

    std::uniform_real_distribution<float> uniform(-1.0, 1.0);
    std::uniform_int_distribution<int> uniform_int(0, N - 1);

    if (push_back_construction) {
      particles.init_id_search();
      if (set_domain) {
        particles.init_neighbour_search(min, max, periodic);
      }
      typename particles_type::value_type p;
      for (int i = 0; i < N; ++i) {
        get<id_to_find>(p) = uniform_int(gen);
        for (size_t d = 0; d < D; ++d) {
          get<position>(p)[d] = uniform(gen);
        }
        particles.push_back(p);
      }
    } else {
      particles.resize(N);
      detail::for_each(
          std::begin(particles), std::end(particles),
          set_random<D, typename particles_type::raw_reference>(-1.0, 1.0, N));

      if (set_domain) {
        particles.init_neighbour_search(min, max, periodic);
      }
      particles.init_id_search();
    }

    /*
    // push 3 more particles
    for (int i=0; i<3; ++i) {
        for (size_t d = 0; d < D; ++d) {
            get<position>(p)[d] = uniform(gen);
        }
        particles.push_back(p);
    }
    */

    // delete random particle
    std::uniform_int_distribution<int> uniform_int_N(0, N - 1);
    particles.erase(particles.begin() + uniform_int_N(gen));

    // delete last particle
    particles.erase(particles.begin() + particles.size() - 1);

    // delete random particle + last and call update on this range
    std::uniform_int_distribution<int> uniform_int_N_minus_2(0, N - 3);
    const size_t random_index = uniform_int_N_minus_2(gen);
    get<alive>(particles)[random_index] = false;
    *get<alive>(particles.end() - 1) = false;
    if (particles.is_ordered()) {
      particles.update_positions(particles.begin(), particles.end());
    } else {
      particles.update_positions(particles.begin() + random_index,
                                 particles.end());
    }

    // brute force search
    auto t0 = Clock::now();
    Aboria::detail::for_each(particles.begin(), particles.end(),
                             brute_force_check<particles_type>(particles));
    auto t1 = Clock::now();
    std::chrono::duration<double> dt_brute = t1 - t0;

    // Aboria search
    t0 = Clock::now();
    Aboria::detail::for_each(particles.begin(), particles.end(),
                             aboria_check<particles_type>(particles));
    t1 = Clock::now();
    std::chrono::duration<double> dt_aboria = t1 - t0;
    for (size_t i = 0; i < particles.size(); ++i) {
      if (int(get<index_brute>(particles)[i]) !=
          int(get<index_aboria>(particles)[i])) {
        std::cout << "error in finding id pid = "
                  << static_cast<const size_t &>(get<id_to_find>(particles)[i])
                  << std::endl;
      }
      TS_ASSERT_EQUALS(int(get<index_brute>(particles)[i]),
                       int(get<index_aboria>(particles)[i]));
    }

    std::cout << "\ttiming result: Aboria = " << dt_aboria.count()
              << " versus brute force = " << dt_brute.count() << std::endl;
  }

  template <template <typename, typename> class VectorType,
            template <typename> class SearchMethod>
  void helper_d_test_list_random() {
    helper_d_random<2, VectorType, SearchMethod>(14, false, false);
    helper_d_random<2, VectorType, SearchMethod>(14, true, false);
    helper_d_random<2, VectorType, SearchMethod>(14, false, true);
    helper_d_random<2, VectorType, SearchMethod>(14, true, true);

    helper_d_random<2, VectorType, SearchMethod>(1000, false, false);
    helper_d_random<2, VectorType, SearchMethod>(1000, true, false);
    helper_d_random<2, VectorType, SearchMethod>(1000, false, true);
    helper_d_random<2, VectorType, SearchMethod>(1000, true, true);
  }

  void test_std_vector_CellList(void) {
    helper_d_test_list_random<std::vector, CellList>();
  }

  void test_std_vector_CellListOrdered(void) {
    helper_d_test_list_random<std::vector, CellListOrdered>();
  }

  void test_std_vector_Kdtree(void) {
#if not defined(__CUDACC__)
    helper_d_test_list_random<std::vector, Kdtree>();
#endif
  }

  void test_std_vector_HyperOctree(void) {
    helper_d_test_list_random<std::vector, HyperOctree>();
  }

  void test_thrust_vector_CellList(void) {
#if defined(HAVE_THRUST)
    // helper_d_test_list_random<thrust::device_vector,CellList>();
#endif
  }

  void test_thrust_vector_CellListOrdered(void) {
#if defined(HAVE_THRUST)
    helper_d_test_list_random<thrust::device_vector, CellListOrdered>();
#endif
  }

  void test_thrust_vector_HyperOctree(void) {
#if defined(HAVE_THRUST)
    helper_d_test_list_random<thrust::device_vector, HyperOctree>();
#endif
  }
};

#endif /* ID_SEARCH_H_ */
