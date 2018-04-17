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

#ifndef TERMINAL_DETAIL_H_
#define TERMINAL_DETAIL_H_

#include "Vector.h"

namespace Aboria {
namespace detail {

////////////////////////
/// Terminal Classes ///
////////////////////////

template <typename I, typename P> struct label {
  typedef P particles_type;
  typedef I depth;

  label(P &p)
      : m_p(p), m_buffers(new typename P::data_type()),
        m_min(new typename P::value_type()),
        m_max(new typename P::value_type()) {}

  P &get_particles() const { return m_p; }
  typename P::data_type &get_buffers() const { return *m_buffers; }
  typename P::value_type &get_min() const { return *m_min; }
  typename P::value_type &get_max() const { return *m_max; }

  P &m_p;
  std::shared_ptr<typename P::data_type> m_buffers;
  std::shared_ptr<typename P::value_type> m_min;
  std::shared_ptr<typename P::value_type> m_max;
};

template <typename T> struct symbolic {
  typedef T variable_type;
  typedef typename T::value_type value_type;
};

/*
   template<typename I>
   struct unknown {
   typedef typename I::value value;
   };
   */

template <typename T> struct accumulate {
  typedef T functor_type;
  typedef typename T::result_type init_type;
  accumulate() : init(VectorTraits<init_type>::Zero()){};
  accumulate(const T &functor)
      : functor(functor), init(VectorTraits<init_type>::Zero()){};
  void set_init(const init_type &arg) { init = arg; }
  T functor;
  init_type init;
};

template <typename T, typename LNormNumber> struct accumulate_within_distance {
  typedef T functor_type;
  typedef LNormNumber norm_number_type;
  typedef typename T::result_type init_type;
  accumulate_within_distance(const double max_distance, const T &functor = T())
      : functor(functor), max_distance(max_distance),
        init(VectorTraits<init_type>::Zero()){};
  void set_init(const init_type &arg) { init = arg; }
  void set_max_distance(const double arg) { max_distance = arg; }
  T functor;
  double max_distance;
  init_type init;
};

template <typename T, unsigned int N> struct vector {
  typedef Vector<T, N> result_type;

  template <typename... Types>
  result_type operator()(const Types... args) const {
    return result_type(args...);
  }
};

template <typename L1, typename L2> struct dx {
  typedef L1 label_a_type;
  typedef L2 label_b_type;
  label_a_type la;
  label_b_type lb;
  dx(label_a_type &la, label_b_type &lb) : la(la), lb(lb){};
  const label_a_type &get_label_a() const { return la; }
  const label_b_type &get_label_b() const { return lb; }
};

struct normal {
  // typedef std::mt19937 generator_type;
  normal(){};
  normal(uint32_t seed) : generator(seed){};
  double operator()() {
    std::normal_distribution<double> normal_distribution;
    return normal_distribution(generator);
  }
  /*
     double operator()(const size_t& id) const {
     std::normal_distribution<double> normal_distribution;
     generator_type gen(id);
     return normal_distribution(gen);
     }
     */
  double operator()(generator_type &gen) const {
    std::normal_distribution<double> normal_distribution;
    return normal_distribution(gen);
  }
  generator_type generator;
};

struct uniform {
  // typedef std::mt19937 generator_type;
  uniform(){};
  uniform(uint32_t seed) : generator(seed){};
  double operator()() {

    std::uniform_real_distribution<double> uniform_distribution;
    return uniform_distribution(generator);
  }
  /*
     double operator()(const size_t& id) const {
     std::uniform_real_distribution<double> uniform_distribution;
     generator_type gen(id);
     return uniform_distribution(gen);
     }
     */
  double operator()(generator_type &gen) const {
    std::uniform_real_distribution<double> normal_distribution;
    return normal_distribution(gen);
  }
  generator_type generator;
};

template <typename T> struct geometry {
  typedef T result_type;

  template <typename... A> result_type operator()(const A &... args) const {
    return T(args...);
  }
};

template <typename T> struct geometries {
  typedef T result_type;

  // result_type operator()(const Vect3d& centre, const double radius, const
  // bool in) const
  template <typename... A> result_type operator()(const A &... args) const {
    return T(args...);
  }
};

} // namespace detail
} // namespace Aboria
#endif
