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

#ifndef ID_SEARCH_H_
#define ID_SEARCH_H_

#include <chrono>
#include <cxxtest/TestSuite.h>
typedef std::chrono::system_clock Clock;

#include "Aboria.h"

using namespace Aboria;

class IDSearchTest : public CxxTest::TestSuite {
public:
  ABORIA_VARIABLE(id_to_find, size_t, "id_to_find")
  ABORIA_VARIABLE(index_brute, int, "found index")
  ABORIA_VARIABLE(index_aboria, int, "found index")

  void test_documentation(void) {
#if not defined(__CUDACC__)
    //[id_search
    /*`
    [section ID Searching]

    As well as neighbourhood searching, Aboria has functionality to search by
    each particle's unique id. First, let's create a set of `N` particles and
    randomly rearrange their position in the vector

    */
    const size_t N = 100;
    typedef Particles<> particle_type;
    particle_type particles(N);
    std::default_random_engine g;
    std::shuffle(particles.begin(), particles.end(), g);

    /*`

    This will create a set of particles that each have a unique id between `0`
    and `N-1` (but their positions in the vector `particles` is randomised).
    Now, we are going to turn on the id search capability for `particles`
    */

    particles.init_id_search();

    /*`

    Then we will try and find the particle with id equal to 2.

    */

    auto id_2 = particles.get_query().find(2);
    //<-
    (void)id_2;
    //->
    assert(*get<id>(id_2) == 2);

    /*`
    Note that each `find` function (e.g. [memberref
    Aboria::CellListQuery::find]) returns an iterator to the
    particle set. If we try and search for an id which doesn't exist, then this
    iterator will point to the end of the particle vector
    */

    auto id_2N = particles.get_query().find(2 * N);
    //<-
    (void)id_2N;
    //->
    assert(id_2N == iterator_to_raw_pointer(particles.end()));

/*`
Finally, a note on performance: The id search is done by internally creating
vectors of id and indicies ordered by id. Keeping these vectors ordered at each
call to [memberref Aboria::Particles::update_positions] takes O(Nlog(N)). The
call to `find` performs a binary search on the ordered vectors, with takes
O(log(N)) time.

[endsect]

*/
//]
#endif
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
