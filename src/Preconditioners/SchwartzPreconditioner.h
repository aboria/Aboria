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

#ifndef SCHWARTZ_PRECONDITIONER_H_
#define SCHWARTZ_PRECONDITIONER_H_

#include "Preconditioners/detail/Preconditioners.h"
#include "hilbert/hilbert.h"
#include <atomic>

#ifdef HAVE_EIGEN
namespace Aboria {

template <typename Operator, typename Query> struct schwartz_decomposition {
  static const unsigned int D = Query::dimension;
  typedef Vector<double, D> double_d;
  typedef position_d<D> position;
  typedef double Scalar;
  typedef size_t Index;

  typedef detail::storage_vector_type<size_t> storage_vector_type;
  typedef detail::storage_vector_type<double> storage_vector_double_type;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> matrix_type;

  const typename Query::child_iterator *m_nodes;
  storage_vector_type *m_domain_indicies;
  storage_vector_type *m_domain_buffer;
  matrix_type *m_domain_matrix;
  size_t m_max_buffer_n;
  size_t m_max_bucket_n;
  const Operator &m_operator;
  Query m_query;
  int m_level;
  bool is_root;

  CUDA_HOST_DEVICE
  void operator()(const int i) const {
    auto &ci = m_nodes[i];
    const auto &bounds = m_query.get_bounds(ci);
    LOG(3, "process_node with bounds "
               << bounds << " is leaf node = " << m_query.is_leaf_node(*ci)
               << " is root = " << is_root);

    double_d middle;
    size_t n = 0;
    {
      for (auto particle = m_query.get_bucket_particles(*ci); particle != false;
           ++particle, ++n) {
        const double_d &p = get<position>(*particle);
        if (n == 0) {
          middle = p;
        } else {
          middle += p;
        }
      }
      middle /= static_cast<double>(n);
    }
    // middle = 0.5 * (bounds.bmax + bounds.bmin);
    const double_d side = 0.5 * (bounds.bmax - bounds.bmin);
    const double search_radius_mult = 0.2;
    const double side_norm = side.norm();
    double search_radius = (1.0 - search_radius_mult) * side_norm;

    // skip over empty buckets on non root levels
    if (!is_root && m_query.get_bucket_particles(*ci) == false) {
      return;
    }

    // find out how many particles and how many buffer particles there are
    size_t n_buffer = 0;
    size_t n_indicies = 0;

    if (is_root) {
      n_buffer = 0;
      n_indicies = m_query.number_of_particles();
    } else {
      /*
      for (auto bucket =
               m_query.template get_buckets_near_point<-1>(middle, 2 * side);
           bucket != false; ++bucket) {
        if (bucket.get_child_iterator() == ci) {
          n_indicies += m_query.get_bucket_particles(*bucket).distance_to_end();
        } else {
          n_buffer += m_query.get_bucket_particles(*bucket).distance_to_end();
        }
      }
      */

      do {
        search_radius += search_radius_mult * side_norm;
        n_indicies = 0;
        n_buffer = 0;
        for (auto particle = euclidean_search(m_query, middle, search_radius);
             particle != false; ++particle) {
          const bool is_ci =
              (particle.get_bucket_iterator().get_child_iterator() == ci);
          if (is_ci) {
            n_indicies++;
          } else {
            n_buffer++;
          }
        }
      } while (n_buffer < m_max_buffer_n);
      ASSERT_CUDA(n_indicies == n);
    }

    // allocate for buffer + indicies
    storage_vector_type buffer_tmp;
    storage_vector_type buffer_distances;
    storage_vector_type indicies_tmp;
    // storage_vector_double_type buffer_r_tmp;
    buffer_tmp.resize(n_buffer);
    buffer_distances.resize(n_buffer);
    indicies_tmp.resize(n_indicies);

    // add particles in bucket to indicies
    // add particles in neighbouring buckets to buffer
    size_t i_indicies = 0;
    size_t i_buffer = 0;
    if (is_root) {
      for (size_t i = 0; i < indicies_tmp.size(); ++i) {
        indicies_tmp[i] = i;
      }
      i_indicies = n_indicies;
    } else {
      /*
      for (auto bucket =
               m_query.template get_buckets_near_point<-1>(middle, 2 * side);
           bucket != false; ++bucket) {
        const bool is_ci = (bucket.get_child_iterator() == ci);
        for (auto particle = m_query.get_bucket_particles(*bucket);
             particle != false; ++particle) {
          const double_d &p = get<position>(*particle);
          const size_t index =
              &p - get<position>(m_query.get_particles_begin());
          if (is_ci) {
            indicies_tmp[i_indicies++] = index;
          } else {
            // buffer_r_tmp[i_buffer] = (p - middle).squaredNorm();
            buffer_tmp[i_buffer++] = index;
          }
        }
      }
      */
      for (auto particle = euclidean_search(m_query, middle, search_radius);
           particle != false; ++particle) {
        const double_d &p = get<position>(*particle);
        const bool is_ci =
            (particle.get_bucket_iterator().get_child_iterator() == ci);
        const size_t index = &p - get<position>(m_query.get_particles_begin());
        if (is_ci) {
          indicies_tmp[i_indicies++] = index;
        } else {
          buffer_distances[i_buffer] = (p - middle).squaredNorm();
          // buffer_distances[i_buffer] = Aboria::abs(p - middle).sum();
          buffer_tmp[i_buffer] = index;
          i_buffer++;
        }
      }
    }

    LOG(3, "\tfound  " << indicies_tmp.size() << " indices in bucket and "
                       << buffer_tmp.size() << " in buffer");

    ASSERT_CUDA(i_indicies == n_indicies);
    ASSERT_CUDA(i_buffer == n_buffer);

    if (buffer_tmp.size() > m_max_buffer_n ||
        indicies_tmp.size() > m_max_bucket_n) {
      // random shuffle
#if defined(__CUDACC__)
      thrust::default_random_engine gen;
#else
      generator_type gen;
#endif
      // advance forward so no random streams intersect
      gen.discard((indicies_tmp.size() + buffer_tmp.size()) * i);

      auto shuffle = [&](storage_vector_type &v) {
        for (int i = v.size() - 1; i > 0; --i) {
#if defined(__CUDACC__)
          thrust::uniform_int_distribution<int> uniform(0, i);
          size_t tmp = v[i];
          const auto random_index = uniform(gen);
          v[i] = v[random_index];
          v[random_index] = tmp;
#else
          std::uniform_int_distribution<int> uniform(0, i);
          std::swap(v[i], v[uniform(gen)]);
#endif
        }
      };

      if (buffer_tmp.size() > m_max_buffer_n) {
        // shuffle(buffer_tmp);
        detail::sort_by_key(buffer_distances.begin(), buffer_distances.end(),
                            buffer_tmp.begin());
      }
      if (indicies_tmp.size() > m_max_bucket_n) {
        shuffle(indicies_tmp);
      }
    }
    // copy random chosen buffer_tmp to buffer along with indicies
    // copy random chosen indicies_tmp to indicies
    auto &buffer = m_domain_buffer[i];
    auto &indicies = m_domain_indicies[i];
    const size_t buffer_size = std::min(buffer_tmp.size(), m_max_buffer_n);
    const size_t indicies_size = std::min(
        indicies_tmp.size(), m_max_bucket_n + m_max_buffer_n - buffer_size);
    buffer.resize(indicies_size + buffer_size);
    /*
    std::cout << "indicies size = " << indicies_size
              << " buffer size = " << buffer_size
              << " requested indicies size is " << m_max_bucket_n
              << " requested buffer size is " << m_max_buffer_n
              << " number of particles = " << m_query.number_of_particles()
              << std::endl;
              */

    // do buffer and indicies
    indicies.resize(indicies_size);
    for (size_t i = 0; i < indicies_size; ++i) {
      buffer[i] = indicies_tmp[i];
      indicies[i] = indicies_tmp[i];
    }

    // finish off with just buffer
    for (size_t i = 0; i < buffer_size; ++i) {
      buffer[indicies_size + i] = buffer_tmp[i];
    }

    LOG(3, "\tafter filtering, found  "
               << indicies.size() << " indices in bucket and " << buffer.size()
               << " in buffer");

    // ASSERT(buffer.size() > 0, "no particles in buffer");
    // ASSERT_CUDA(indicies.size() > 0);

    // fill domain matrix
    const int size = buffer.size();
    auto &domain_matrix = m_domain_matrix[i];
    domain_matrix.resize(size, size);
    for (size_t i = 0; i < buffer.size(); ++i) {
      for (size_t j = 0; j < buffer.size(); ++j) {
        domain_matrix(i, j) = m_operator.coeff(buffer[i], buffer[j]);
      }
    }
  }
};

template <typename T, typename Query> struct schwartz_matrices_from_kernel {
  static const unsigned int D = Query::dimension;
  typedef Vector<double, D> double_d;
  typedef position_d<D> position;
  typedef double Scalar;
  typedef size_t Index;

  typedef detail::storage_vector_type<size_t> storage_vector_type;
  typedef detail::storage_vector_type<double> storage_vector_double_type;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> matrix_type;

  const typename Query::child_iterator *m_nodes;
  storage_vector_type *m_domain_buffer;
  matrix_type *m_domain_matrix;
  const T m_function;
  Query m_query;

  CUDA_HOST_DEVICE
  void operator()(const int i) const {

    auto &buffer = m_domain_buffer[i];

    // fill domain matrix
    const int size = buffer.size();
    auto &domain_matrix = m_domain_matrix[i];
    domain_matrix.resize(size, size);
    for (size_t i = 0; i < buffer.size(); ++i) {
      auto a = m_query.get_particles_begin()[buffer[i]];
      for (size_t j = 0; j < buffer.size(); ++j) {
        auto b = m_query.get_particles_begin()[buffer[j]];
        domain_matrix(i, j) = m_function(a, b);
      }
    }
  }
};

template <typename T, typename Query>
struct schwartz_matrices_from_interpolation {
  static const unsigned int D = Query::dimension;
  typedef Vector<double, D> double_d;
  typedef position_d<D> position;
  typedef double Scalar;
  typedef size_t Index;

  typedef detail::storage_vector_type<size_t> storage_vector_type;
  typedef detail::storage_vector_type<double> storage_vector_double_type;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> matrix_type;

  const typename Query::child_iterator *m_nodes;
  storage_vector_type *m_domain_buffer;
  matrix_type *m_domain_matrix;
  matrix_type &m_finer_operator;
  Query m_query;

  CUDA_HOST_DEVICE
  void operator()(const int i) const {

    auto &buffer = m_domain_buffer[i];

    // fill domain matrix
    const int size = buffer.size();
    auto &domain_matrix = m_domain_matrix[i];
    domain_matrix.resize(size, size);
    for (size_t i = 0; i < buffer.size(); ++i) {
      auto a = m_query.get_particles_begin()[buffer[i]];
      for (size_t j = 0; j < buffer.size(); ++j) {
        auto b = m_query.get_particles_begin()[buffer[j]];
        domain_matrix(i, j) = m_function(a, b);
      }
    }
  }
};

template <typename Solver> struct factorize_matrix {
  typedef double Scalar;

  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> matrix_type;

  Solver operator()(matrix_type &domain_matrix) {
    LOG(3, "factorize matrix with size (" << domain_matrix.rows() << ','
                                          << domain_matrix.cols() << ')');

    Solver solver;
    if (domain_matrix.rows() == 0 || domain_matrix.cols() == 0) {
      return solver;
    }

    solver.compute(domain_matrix);

    Eigen::VectorXd b = Eigen::VectorXd::Random(domain_matrix.rows());
    Eigen::VectorXd x = solver.solve(b);
    double relative_error = (domain_matrix * x - b).norm() / b.norm();
    if (relative_error > 1e-3 || std::isnan(relative_error)) {
      std::cout << "ERROR in domain solve: relative error = " << relative_error
                << std::endl;
    }
    return solver;
  }
};

template <typename T> struct factorize_matrix<Eigen::ColPivHouseholderQR<T>> {
  typedef double Scalar;
  typedef Eigen::ColPivHouseholderQR<T> Solver;

  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> matrix_type;

  Solver operator()(matrix_type &domain_matrix) {
    LOG(3, "factorize matrix with size (" << domain_matrix.rows() << ','
                                          << domain_matrix.cols() << ')');

    Solver solver;
    if (domain_matrix.rows() == 0 || domain_matrix.cols() == 0) {
      return solver;
    }

    solver.compute(domain_matrix);

    Eigen::VectorXd b = Eigen::VectorXd::Random(domain_matrix.rows());
    Eigen::VectorXd x = solver.solve(b);
    double relative_error = (domain_matrix * x - b).norm() / b.norm();
    if (relative_error > 1e-3 || std::isnan(relative_error)) {
      std::cout << "ERROR in domain solve: relative error = " << relative_error
                << std::endl;
      std::cout << "rank of (" << domain_matrix.rows() << ","
                << domain_matrix.cols() << ") matrix is " << solver.rank()
                << std::endl;
    }
    return solver;
  }
};

template <typename Rhs, typename Dest, typename Solver> struct solve_domain {
  typedef double Scalar;
  typedef detail::storage_vector_type<size_t> storage_vector_type;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> matrix_type;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> vector_type;
  typedef Solver solver_type;

  const storage_vector_type *m_domain_indicies;
  const storage_vector_type *m_domain_buffer;
  const solver_type *domain_factorized_matrix;
  Dest &x;
  const Rhs &b;
  const Eigen::VectorXd &m_count;
  mutable std::atomic_flag lock;

  solve_domain(const storage_vector_type *m_domain_indicies,
               const storage_vector_type *m_domain_buffer,
               const solver_type *domain_factorized_matrix, Dest &x,
               const Rhs &b, const Eigen::VectorXd &m_count)
      : m_domain_indicies(m_domain_indicies), m_domain_buffer(m_domain_buffer),
        domain_factorized_matrix(domain_factorized_matrix), x(x), b(b),
        m_count(m_count), lock(ATOMIC_FLAG_INIT) {}

  solve_domain(const solve_domain &other)
      : m_domain_indicies(other.m_domain_indicies),
        m_domain_buffer(other.m_domain_buffer),
        domain_factorized_matrix(other.domain_factorized_matrix), x(other.x),
        b(other.b), m_count(other.m_count), lock(ATOMIC_FLAG_INIT) {}

  void operator()(const int i) const {
    const storage_vector_type &buffer = m_domain_buffer[i];
    const storage_vector_type &indicies = m_domain_indicies[i];
    if (buffer.size() == 0)
      return;

    const size_t nb = buffer.size();

    vector_type domain_x;
    vector_type domain_b;
    domain_x.resize(nb);
    domain_b.resize(nb);

    // copy b values from big vector
    for (size_t j = 0; j < buffer.size(); ++j) {
      domain_b[j] = m_count[buffer[j]] * b[buffer[j]];
      // domain_b[j] = b[buffer[j]];
    }

    // solve domain
    domain_x = domain_factorized_matrix[i].solve(domain_b);

    while (lock.test_and_set(std::memory_order_acquire)) // acquire lock
      ;                                                  // spin
    // copy accumulate x values to big vector
    for (size_t j = 0; j < indicies.size(); ++j) {
      x[indicies[j]] += m_count[buffer[j]] * domain_x[j];
      // x[indicies[j]] += domain_x[j];
    }

    // non-restricted
    for (size_t j = indicies.size(); j < buffer.size(); ++j) {
      x[buffer[j]] += m_count[buffer[j]] * domain_x[j];
      // x[buffer[j]] += domain_x[j];
    }
    lock.clear(std::memory_order_release); // release lock
  }
};

template <typename Operator, typename Solver> class SchwartzDecomposition {
  typedef double Scalar;
  typedef size_t Index;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> matrix_type;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> vector_type;
  typedef Solver solver_type;

  typedef detail::storage_vector_type<size_t> storage_vector_type;
  typedef detail::storage_vector_type<double> storage_vector_double_type;

  typedef std::vector<std::vector<storage_vector_type>> connectivity_t;
  typedef std::vector<std::vector<solver_type>> solvers_t;
  typedef std::vector<std::vector<matrix_type>> matrices_t;

  typedef
      typename std::remove_cv<typename std::remove_reference<Operator>::type>::
          type::FirstBlock kernel_t;
  static_assert(std::is_same<typename kernel_t::row_elements_type,
                             typename kernel_t::col_elements_type>::value,
                "Schwartz preconditioner restricted to identical row and col "
                "particle sets");
  static_assert(kernel_t::BlockRows == 1 && kernel_t::BlockCols == 1,
                "Schwartz preconditioner not currently implemented for "
                "vector-valued kernels");
  typedef typename kernel_t::row_elements_type particles_t;
  typedef typename particles_t::query_type query_t;
  typedef MatrixReplacement<1, 1, std::tuple<kernel_t>> operator_t;
  typedef std::vector<operator_t> level_operators_t;
  typedef std::vector<matrix_type> interpolation_operators_t;
  typedef std::vector<size_t> level_connectivity_t;
  typedef std::vector<particles_t> level_particles_t;

private:
  size_t m_max_buffer_n;
  size_t m_mult_buffer;

  connectivity_t m_indicies;
  connectivity_t m_buffer;
  solvers_t m_factorized_matrix;
  matrices_t m_matrix;
  level_operators_t m_A;
  const operator_t *m_fineA;
  level_particles_t m_particles;
  interpolation_operators_t m_interpolation_matricies;
  std::default_random_engine m_generator;
  std::vector<vector_type> m_count;
  double m_smoother_weighting;

public:
  SchwartzDecomposition()
      : m_max_buffer_n(1), m_mult_buffer(1), m_fineA(nullptr),
        m_smoother_weighting(1.0) {}

  void set_max_buffer_n() {
    CHECK(m_fineA, "need to call set_operator() first");
    const size_t n = m_fineA->rows();
    const size_t nBucket =
        m_fineA->get_first_kernel().get_row_elements().get_max_bucket_size();
    m_max_buffer_n =
        m_mult_buffer * std::pow(n, 1.0 / 3.0) * std::pow(nBucket, 2.0 / 3.0);
    CHECK(m_max_buffer_n <= n,
          "coarse size is larger than number of particles");
  }
  void set_mult_buffer(size_t arg) {
    m_mult_buffer = arg;
    if (m_fineA) {
      set_max_buffer_n();
    }
  }
  void set_smoother_weighting(double arg) { m_smoother_weighting = arg; }

  void set_operator(const Operator &mat) {
    m_fineA = &mat;
    set_max_buffer_n();
  }
  const operator_t &get_operator() { return *m_fineA; }
  const operator_t &get_operator_at_level(const int i) {
    ASSERT(i > 0, "need i > 0");
    return m_A[i - 1];
  }

  size_t get_n_levels() const { return m_particles.size() + 1; }

  size_t get_size_of_level(int i) const {
    return i == 0 ? m_fineA->get_first_kernel().get_row_elements().size()
                  : m_particles[i - 1].size();
  }
  const level_particles_t &get_particles() const { return m_particles; }
  const connectivity_t &get_indicies() const { return m_indicies; }
  const connectivity_t &get_buffer() const { return m_buffer; }

  void construct_levels(int n_levels) {
    m_particles.clear();
    m_particles.reserve(n_levels - 1);

    // need to emplace back these....
    m_A.clear();
    m_interpolation_matricies.clear();

    m_indicies.resize(n_levels);
    m_buffer.resize(n_levels);
    m_matrix.resize(n_levels);

    construct_decomposition_level(0, n_levels);
    for (int i = 0; i < n_levels - 1; ++i) {

      construct_particle_level(i, n_levels);
      construct_operator_level(i);
      construct_decomposition_level(i + 1, n_levels);
    }
  }

  void construct_particle_level(int i, int n_levels) {
    LOG(2, "SchwartzDecomposition: creating particles on level " << i + 1);
    const auto &a = m_fineA->get_first_kernel().get_row_elements();
    const int n_particles =
        i == n_levels - 2 ? m_max_buffer_n : a.size() / std::pow(2, i + 1);
    const int n_finer = i == 0 ? a.size() : m_particles[i - 1].size();

    // create new particle set for this level

    m_particles.emplace_back(n_particles);
    auto &particles = m_particles.back();
    typedef typename particles_t::position position;

    if (0) {
      std::vector<bitmask_t> hindex(n_finer);
      std::vector<Vector<bitmask_t, query_t::dimension>> hcoord(n_finer);
      const int nBits =
          std::numeric_limits<bitmask_t>::digits / query_t::dimension - 1;
      std::cout << "using nBits = " << nBits << std::endl;
      bitmask_t max = static_cast<bitmask_t>(1) << nBits;
      typedef typename particles_t::double_d double_d;
      const auto &prev_level_position =
          (i == 0 ? get<position>(a) : get<position>(m_particles[i - 1]));
      detail::transform(prev_level_position.begin(), prev_level_position.end(),
                        hcoord.begin(), [&](const double_d &p) {
                          Vector<bitmask_t, query_t::dimension> ret;
                          for (size_t i = 0; i < query_t::dimension; ++i) {
                            ret[i] = max * p[i];
                          }
                          return ret;
                        });
      detail::transform(hcoord.begin(), hcoord.end(), hindex.begin(),
                        [&](const auto &c) {
                          bitmask_t tmp[query_t::dimension];
                          for (size_t i = 0; i < query_t::dimension; ++i) {
                            tmp[i] = c[i];
                          }
                          return hilbert_c2i(query_t::dimension, nBits, tmp);
                        });
      level_connectivity_t chosen_indices_tmp(n_finer);
      detail::tabulate(chosen_indices_tmp.begin(), chosen_indices_tmp.end(),
                       [](const int i) { return i; });
      detail::sort_by_key(hindex.begin(), hindex.end(),
                          chosen_indices_tmp.begin());
      level_connectivity_t chosen_indices(n_particles);
      const double ratio = static_cast<double>(n_finer) / n_particles;
      detail::tabulate(chosen_indices.begin(), chosen_indices.end(),
                       [&](const int i) {
                         return chosen_indices_tmp[static_cast<int>(ratio * i)];
                       });
    }

    // choose every second particle
    level_connectivity_t chosen_indices(n_particles);
    const double ratio = static_cast<double>(n_finer) / n_particles;
    detail::tabulate(chosen_indices.begin(), chosen_indices.end(),
                     [&](const int i) { return ratio * i; });

    if (0) {
      // count domains per particle for previous level
      level_connectivity_t count(n_finer);
      count_domains_per_particle(count, i);
      // sort by the count and get rid of larger counts
      level_connectivity_t chosen_indices(n_finer);
      detail::tabulate(chosen_indices.begin(), chosen_indices.end(),
                       [](const int i) { return i; });
      detail::sort_by_key(count.begin(), count.end(), chosen_indices.begin());
      chosen_indices.resize(n_particles);
    }

    if (0) {
      level_connectivity_t chosen_indices(n_finer);
      detail::tabulate(chosen_indices.begin(), chosen_indices.end(),
                       [](const int i) { return i; });
      while (chosen_indices.size() > n_particles) {
        double min_dist = std::numeric_limits<double>::max();
        int min_index = -1;
        for (int i = 0; i < chosen_indices.size(); ++i) {
          for (int j = 0; j < chosen_indices.size(); ++j) {
            const double dist =
                (get<position>(particles)[i] - get<position>(particles)[j])
                    .squaredNorm();
            if (dist < min_dist) {
              min_index = i;
              min_dist = dist;
            }
          }
        }
        chosen_indices.erase(chosen_indices.begin() + min_index);
      }
    }

    if (0) {
      // choose a random set of particles
      level_connectivity_t chosen_indices(n_finer);
      detail::tabulate(chosen_indices.begin(), chosen_indices.end(),
                       [](const int i) { return i; });
      detail::random_unique(chosen_indices.begin(), chosen_indices.end(),
                            n_particles, m_generator);
      chosen_indices.resize(n_particles);
    }

    // TODO: need to copy chosen_indicies to gpu if particles are there...

    detail::tabulate(
        particles.begin(), particles.end(),
        [chosen_indices =
             Aboria::iterator_to_raw_pointer(chosen_indices.begin()),
         prev_level = iterator_to_raw_pointer(
             (i == 0 ? a.cbegin() : m_particles[i - 1].cbegin()))](
            const int index) { return prev_level[chosen_indices[index]]; });

    // set particle ids to point to finer level index so that we can
    // implement the restriction operator later on...
    detail::tabulate(get<id>(particles).begin(), get<id>(particles).end(),
                     [chosen_indices = Aboria::iterator_to_raw_pointer(
                          chosen_indices.begin())](const int index) {
                       return chosen_indices[index];
                     });

    const double ratio2 = static_cast<double>(a.size()) / particles.size();
    const int max_bucket_size = std::sqrt(ratio2) * a.get_max_bucket_size();
    // const int max_bucket_size = a.get_max_bucket_size();

    LOG(2, "\t: creating new particle set with "
               << particles.size() << " particles and with max bucket size of "
               << max_bucket_size);
    particles.init_neighbour_search(a.get_min(), a.get_max(), a.get_periodic(),
                                    max_bucket_size);

    double min_dist = std::numeric_limits<double>::max();
    for (size_t i = 0; i < particles.size(); ++i) {
      for (size_t j = 0; j < particles.size(); ++j) {
        if (i != j) {
          const double dist =
              (get<position>(particles)[i] - get<position>(particles)[j])
                  .squaredNorm();
          if (dist < min_dist)
            min_dist = dist;
        }
      }
    }
    LOG(2, "\t: minimum distance is " << std::sqrt(min_dist));

    // TODO: dodgy!!!
    // m_A.back().get_first_kernel().get_matrix() *= volume_scale_factor;

    // create interpolation matrix to interpolate back to finest level
    /*
    std::vector<Eigen::Triplet<double>> triplets;
    for (int i = 0; i < a.size(); ++i) {
      for (auto j = euclidean_search(a.get_query(), radius); j != false;
           ++j) {
        triplets.push_back(Eigen::Triplet<double>(i,j_index,kernel);
      }
    }
    SpMat A(m, m);
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();
    */
  }

  struct wendland_interpolate {
    typedef typename particles_t::position position;

    interpolation_operators_t::reference m_interpolation_matricies;
    const query_t &m_query_source;
    const query_t &m_query_dest;
    const double m_interpolation_length_scale;

    void operator()(const int i) {
      const auto &pi = m_query_dest.get_particles_begin()[i];

      double sum_weight = 0;
      const double c = 1.0 / m_interpolation_length_scale;
      for (auto j = euclidean_search(m_query_source, get<position>(pi),
                                     2 * m_interpolation_length_scale);
           j != false; ++j) {
        const int j_index = &(get<position>(*j)) -
                            get<position>(m_query_source.get_particles_begin());
        const double r = j.dx().norm();
        const double weight = std::pow(2.0 - r * c, 4) * (1.0 + 2.0 * r * c);
        m_interpolation_matricies(i, j_index) = weight;
        sum_weight += weight;
      }
      const double inv_sum_weight = sum_weight == 0.0 ? 0.0 : 1.0 / sum_weight;
      for (auto j = euclidean_search(m_query_source, get<position>(pi),
                                     2 * m_interpolation_length_scale);
           j != false; ++j) {
        const int j_index = &(get<position>(*j)) -
                            get<position>(m_query_source.get_particles_begin());
        m_interpolation_matricies(i, j_index) *= inv_sum_weight;
      }
    }
  };

  double get_interpolation_length_scale(const int i) {
    ASSERT(i >= 0, "only valid for i >= 0");
    ASSERT(i < get_n_levels() - 1, "only valid for i < n_levels-1");
    const double volume =
        (m_particles[i].get_max() - m_particles[i].get_min()).prod();
    const double bucket_volume =
        m_particles[i].get_max_bucket_size() * volume / m_particles[i].size();
    return 1.0 * std::pow(bucket_volume, 1.0 / m_particles[i].dimension);
  }

  void construct_operator_level(int i) {
    LOG(2, "SchwartzDecomposition: creating operator on level " << i + 1);
    m_A.emplace_back(
        kernel_t(m_particles[i], m_particles[i], m_fineA->get_first_kernel()));
  }

  void construct_decomposition_level(int i, int n_levels) {
    LOG(2, "SchwartzDecomposition: creating decomposition level " << i);

    const auto &a = m_fineA->get_first_kernel().get_row_elements();

    // get list of leaf cells
    const query_t &query =
        i == 0 ? a.get_query() : m_particles[i - 1].get_query();
    auto df_search = query.breadth_first();
    for (; df_search != false; ++df_search) {
    }

    typename query_t::breadth_first_iterator::value_type root_node_vector = {
        query.get_root()};
    auto &nodes = i == n_levels - 1 ? root_node_vector : *df_search;

    LOG(2, "\t: processing level with " << nodes.size() << " nodes");

    auto &indicies = m_indicies[i];
    auto &buffer = m_buffer[i];
    auto &matrix = m_matrix[i];

    indicies.resize(nodes.size());
    buffer.resize(nodes.size());
    matrix.resize(nodes.size());

    // const int num_threads = omp_get_max_threads();
    const int max_bucket_size = i == 0
                                    ? a.get_max_bucket_size()
                                    : m_particles[i - 1].get_max_bucket_size();
    const size_t size = i == 0 ? a.size() : m_particles[i - 1].size();
    const size_t mult_buffer = m_mult_buffer;
    const size_t max_buffer = (nodes.size() == 1)
                                  ? m_max_buffer_n
                                  : max_bucket_size * (mult_buffer - 1);

    const size_t max_indicies = max_bucket_size;

    CHECK(max_buffer <= size,
          "number of buffer particles larger than size of level")
    CHECK(max_indicies <= size,
          "number of bucket particles larger than size of level")

    using traits_t = typename query_t::traits_type;
    if (traits_t::data_on_GPU) {
      // need to copy data from/to gpu
      typename traits_t::template vector<storage_vector_type> tmp_indicies(
          nodes.size());
      typename traits_t::template vector<storage_vector_type> tmp_buffer(
          nodes.size());
      typename traits_t::template vector<matrix_type> tmp_matrix(nodes.size());

      auto count = traits_t::make_counting_iterator(0);

      if (i == 0) {
        detail::for_each(
            count, count + nodes.size(),
            schwartz_decomposition<operator_t, query_t>{
                Aboria::iterator_to_raw_pointer(nodes.begin()),
                Aboria::iterator_to_raw_pointer(tmp_indicies.begin()),
                Aboria::iterator_to_raw_pointer(tmp_buffer.begin()),
                Aboria::iterator_to_raw_pointer(tmp_matrix.begin()), max_buffer,
                max_indicies, *m_fineA, query, i, nodes.size() == 1});
      } else {
        detail::for_each(
            count, count + nodes.size(),
            schwartz_decomposition<operator_t, query_t>{
                Aboria::iterator_to_raw_pointer(nodes.begin()),
                Aboria::iterator_to_raw_pointer(tmp_indicies.begin()),
                Aboria::iterator_to_raw_pointer(tmp_buffer.begin()),
                Aboria::iterator_to_raw_pointer(tmp_matrix.begin()), max_buffer,
                max_indicies, m_A[i - 1], query, i, nodes.size() == 1});
      }

      detail::copy(tmp_indicies.begin(), tmp_indicies.end(), indicies.begin());
      detail::copy(tmp_buffer.begin(), tmp_buffer.end(), buffer.begin());
      detail::copy(tmp_matrix.begin(), tmp_matrix.end(), matrix.begin());
    } else {
      // no copy required
      auto count = traits_t::make_counting_iterator(0);
      if (i == 0) {
        detail::for_each(count, count + nodes.size(),
                         schwartz_decomposition<operator_t, query_t>{
                             Aboria::iterator_to_raw_pointer(nodes.begin()),
                             Aboria::iterator_to_raw_pointer(indicies.begin()),
                             Aboria::iterator_to_raw_pointer(buffer.begin()),
                             Aboria::iterator_to_raw_pointer(matrix.begin()),
                             max_buffer, max_indicies, *m_fineA, query, i,
                             nodes.size() == 1});
      } else {
        detail::for_each(count, count + nodes.size(),
                         schwartz_decomposition<operator_t, query_t>{
                             Aboria::iterator_to_raw_pointer(nodes.begin()),
                             Aboria::iterator_to_raw_pointer(indicies.begin()),
                             Aboria::iterator_to_raw_pointer(buffer.begin()),
                             Aboria::iterator_to_raw_pointer(matrix.begin()),
                             max_buffer, max_indicies, m_A[i - 1], query, i,
                             nodes.size() == 1});
      }
    }
  }

  void construct_factorized_levels(int n_levels) {
    m_factorized_matrix.resize(n_levels);

    for (int i = 0; i < n_levels; ++i) {
      LOG(2, "SchwartzDecomposition: factorising domain matrix level " << i);
      auto &matrix = m_matrix[i];
      auto &factorized_matrix = m_factorized_matrix[i];
      factorized_matrix.resize(matrix.size());

      detail::transform(matrix.begin(), matrix.end(), factorized_matrix.begin(),
                        factorize_matrix<solver_type>());
    }
  }

  void count_domains_per_particle(vector_type &out, int level) {
    const int n = level == 0
                      ? m_fineA->get_first_kernel().get_row_elements().size()
                      : m_particles[level - 1].size();
    out = vector_type::Zero(n);
    for (const auto &buffer : m_buffer[level]) {
      detail::for_each(buffer.begin(), buffer.end(),
                       [&](const size_t &index) { out[index] += 1; });
    }
  }

  void count_domains_per_particle(std::vector<size_t> &out, int level) {
    detail::fill(out.begin(), out.end(), 0);
    for (const auto &buffer : m_buffer[level]) {
      detail::for_each(buffer.begin(), buffer.end(),
                       [&](const size_t &index) { out[index] += 1; });
    }
  }

  void construct_domain_particle_count(int n_levels) {
    m_count.resize(n_levels);
    for (int i = 0; i < n_levels; ++i) {
      count_domains_per_particle(m_count[i], i);
      const int max_count = m_count[i].maxCoeff();
      LOG(2, "SchwartzDecomposition: level "
                 << i << " domain count max = " << max_count)
      LOG(2, "SchwartzDecomposition: level "
                 << i << " domain count min= " << m_count[i].minCoeff())
      LOG(2, "SchwartzDecomposition: level "
                 << i << " domain count averate= " << m_count[i].mean())
      // m_count[i] = m_count[i].cwiseInverse();
      // m_count[i].setConstant(1.0 / std::sqrt(m_count[i].mean()));
      // m_count[i].setConstant(1.0);
      if (i != n_levels - 1) {
        // m_count[i] = m_count[i].cwiseInverse().cwiseSqrt();
        m_count[i].setConstant(
            std::sqrt(1.0 / (2.0 * n_levels * m_count[i].maxCoeff())));
        // m_count[i].setConstant(
        //    std::sqrt(1.0 / (n_levels * m_count[i].maxCoeff())));

        // m_count[i] *= std::sqrt(m_smoother_weighting);
        // m_count[i] *= std::sqrt(0.1);
        // m_count[i].setConstant(1.0);
      } else {
        m_count[i].setConstant(1.0);
      }
    }
  }

  void hierarchical_schwartz_decomposition(int max_levels, bool interpolate) {
    const kernel_t &kernel = m_fineA->get_first_kernel();
    const particles_t &a = kernel.get_row_elements();

    CHECK(&a == &(kernel.get_col_elements()),
          "Schwartz preconditioner restricted to identical row and col "
          "particle "
          "sets");

    const int potential_n_levels =
        std::ceil((std::log(a.size()) - std::log(m_max_buffer_n)) /
                  std::log(2)) +
        1;
    const int n_levels = std::min(max_levels, potential_n_levels);
    LOG(2, "SchwartzDecomposition: creating " << n_levels << " levels");

    construct_levels(n_levels);
    construct_domain_particle_count(n_levels);
    construct_factorized_levels(n_levels);

    for (size_t i = 0; i < m_indicies.size(); ++i) {
      auto min_indicies = detail::transform_reduce(
          m_indicies[i].begin(), m_indicies[i].end(), detail::get_size(),
          std::numeric_limits<size_t>::max(), detail::min());

      auto max_indicies = detail::transform_reduce(
          m_indicies[i].begin(), m_indicies[i].end(), detail::get_size(),
          std::numeric_limits<size_t>::min(), detail::max());

      auto min_buffer = detail::transform_reduce(
          m_buffer[i].begin(), m_buffer[i].end(), detail::get_size(),
          std::numeric_limits<size_t>::max(), detail::min());

      auto max_buffer = detail::transform_reduce(
          m_buffer[i].begin(), m_buffer[i].end(), detail::get_size(),
          std::numeric_limits<size_t>::min(), detail::max());

      auto count =
          detail::transform_reduce(m_indicies[i].begin(), m_indicies[i].end(),
                                   detail::get_size(), 0, detail::plus());

      const int rows = i == 0 ? m_fineA->rows() : m_A[i - 1].rows();
      const int cols = i == 0 ? m_fineA->cols() : m_A[i - 1].cols();
      auto &particles = i == 0 ? kernel.get_row_elements() : m_particles[i - 1];

      LOG(2, "SchwartzPreconditioner: finished factorizing, on level "
                 << i << " found " << m_indicies[i].size() << " domains, with "
                 << min_indicies << "--" << max_indicies << " particles ("
                 << count << " total), and " << min_buffer << "--" << max_buffer
                 << " buffer particles. particles size is " << particles.size()
                 << ". kernel matrix is size (" << rows << ',' << cols << ")");
      if (interpolate && i > 0) {
        LOG(2, "\tInterpolation matrix size: ("
                   << m_interpolation_matricies[i - 1].rows() << "x"
                   << m_interpolation_matricies[i - 1].cols() << ")");
      }
    }
  }

  template <typename Source, typename Dest>
  void solve_level(const Source &source, Dest &dest, const int i) const {
    ASSERT(source.size() == dest.size(), "source and dest not equal size");
    const int top_level = get_n_levels() - 1;
    ASSERT(static_cast<size_t>(source.size()) == get_size_of_level(i),
           "source and dest not the same size as particles");

    if (i != top_level) {
      auto count = boost::make_counting_iterator(0);
      detail::for_each(
          count, count + m_indicies[i].size(),
          solve_domain<Source, Dest, solver_type>(
              Aboria::iterator_to_raw_pointer(m_indicies[i].begin()),
              Aboria::iterator_to_raw_pointer(m_buffer[i].begin()),
              Aboria::iterator_to_raw_pointer(m_factorized_matrix[i].begin()),
              dest, source, m_count[i]));
    } else {
      solve_domain<Source, Dest, solver_type>(
          Aboria::iterator_to_raw_pointer(m_indicies[i].begin()),
          Aboria::iterator_to_raw_pointer(m_buffer[i].begin()),
          Aboria::iterator_to_raw_pointer(m_factorized_matrix[i].begin()), dest,
          source, m_count[i])(0);
    }
  }

  template <typename Source, typename Dest>
  void restrict_multi_level_up(const Source &source, Dest &dest,
                               const int n_levels) const {
    ASSERT(n_levels > 0, "n_levels greater than 1");
    dest.resize(m_particles[n_levels - 1].size());
    ASSERT(source.size() ==
               m_fineA->get_first_kernel().get_row_elements().size(),
           "source should be same size as finest particles layer");
    // restrict b to m_r[m_i]
    detail::tabulate(dest.data(), dest.data() + dest.size(),
                     [&source, n_levels, this](const int index) {
                       int curr_index = index;
                       for (int i = n_levels - 1; i >= 0; --i) {
                         curr_index = get<id>(m_particles[i])[curr_index];
                       }
                       return source[curr_index];
                     });
  }

  template <typename Source, typename Dest>
  void restrict_level_up(const Source &source, Dest &dest, const int i) const {
    dest.resize(m_particles[i].size());
    ASSERT(i < static_cast<int>(get_n_levels()), "already at highest level");
    ASSERT(static_cast<size_t>(source.size()) ==
               (i == 0 ? m_fineA->get_first_kernel().get_row_elements().size()
                       : m_particles[i - 1].size()),
           "source size is incorrect");

    // restrict b to m_r[m_i]
    const bool use_interpolation_mat = m_interpolation_matricies.size() != 0;
    if (use_interpolation_mat) {
      dest = m_interpolation_matricies[i].transpose() * source;
    } else {
      detail::tabulate(dest.data(), dest.data() + dest.size(),
                       [&source, id = get<id>(m_particles[i].begin())](
                           const int index) { return source[id[index]]; });
    }
  }

  template <typename Source, typename Dest>
  void interpolate_multi_level_down(const Source &source, Dest &dest,
                                    const int n_levels) const {

    ASSERT(n_levels > 0, "n_levels greater than 1");
    ASSERT(dest.size() == m_fineA->get_first_kernel().get_row_elements().size(),
           "dest should be same size as finest particles layer");
    ASSERT(source.size() == m_particles[n_levels - 1].size(),
           "source should be same size as num particles on this level");

    const bool use_interpolation_mat = m_interpolation_matricies.size() != 0;
    if (use_interpolation_mat) {
      ASSERT(m_interpolation_matricies.size() > n_levels - 1,
             "no interpolation matrix!");
      dest += m_interpolation_matricies[n_levels - 1] * source;
    } else {
      // accumulate m_u[m_i] to x
      auto count = boost::make_counting_iterator(0);
      detail::for_each(count, count + source.size(),
                       [&dest, &source, this, n_levels](const int index) {
                         int curr_index = index;
                         for (int i = n_levels - 1; i >= 0; --i) {
                           curr_index = get<id>(m_particles[i])[curr_index];
                         }
                         dest[curr_index] += source[index];
                       });
    }
  }

  template <typename Source, typename Dest>
  void interpolate_level_down(const Source &source, Dest &dest,
                              const int i) const {

    ASSERT(i > 0, "already at lowest level");
    ASSERT(static_cast<size_t>(source.size()) == m_particles[i - 1].size(),
           "source should be same size as num particles on this level");
    ASSERT(static_cast<size_t>(dest.size()) == get_size_of_level(i - 1),
           "dest should be same size as finest level particles");

    const bool use_interpolation_mat = m_interpolation_matricies.size() != 0;
    if (use_interpolation_mat) {
      ASSERT(m_interpolation_matricies.size() > static_cast<size_t>(i - 1),
             "no interpolation matrix!");
      dest += m_interpolation_matricies[i - 1] * source;
    } else {
      // accumulate m_u[m_i] to x
      auto count = boost::make_counting_iterator(0);
      detail::for_each(
          count, count + source.size(),
          [&dest, &source, id = get<id>(m_particles[i - 1].begin())](
              const int index) { dest[id[index]] += source[index]; });
    }
  }

  template <typename X_Vector, typename B_Vector>
  void gauss_seidel(const B_Vector &b, X_Vector &x_k, const int level,
                    const bool lower) const {
    const int n = b.size();
    const auto &A =
        level == 0 ? m_fineA->get_first_kernel().get_matrix() : m_A[level - 1];
    ASSERT(x_k.size() == n, "x_k has inconsistent size");
    ASSERT(A.rows() == n, "operator has inconsistent rows");
    ASSERT(A.cols() == n, "operator has inconsistent rows");

    if (lower) {
      x_k = b - A.template triangularView<Eigen::StrictlyUpper>() * x_k;
      A.template triangularView<Eigen::Lower>().solveInPlace(x_k);
    } else {
      x_k = b - A.template triangularView<Eigen::StrictlyLower>() * x_k;
      A.template triangularView<Eigen::Upper>().solveInPlace(x_k);
    }
  }

  void v_cycle(std::vector<vector_type> &m_u, std::vector<vector_type> &m_r,
               int level) const {
    LOG(3, "\trunning v_cycle on level " << level);
    const size_t i = level;
    if (level == static_cast<int>(get_n_levels()) - 1) {
      LOG(3, "\t\tdirect solve on level " << level
                                          << ". r_norm = " << m_r[i].norm()
                                          << " u_norm = " << m_u[i].norm());
      solve_level(m_r[i], m_u[i], i);
    } else {
      // pre-smoothing (u <- u + A-1 r)
      vector_type residual;
      if (i == 0) {
        residual = m_r[i] - (*m_fineA) * m_u[i];
      } else {
        residual = m_r[i] - m_A[i - 1] * m_u[i];
      }
      // for (int j = 0; j < std::pow(2, level); ++j) {
      for (int j = 0; j < 1; ++j) {
        LOG(3, "\t\tpre-smoothing on level "
                   << level << " residual_norm = " << residual.norm()
                   << " u_norm = " << m_u[i].norm());
        solve_level(residual, m_u[i], i);
        /*
        std::cout << "level " << i << " residual norm = " << residual.norm()
                  << std::endl;
        std::cout << "level " << i << " u norm = " << m_u[i].norm()
                  << std::endl;
                  */
        if (i == 0) {
          residual = m_r[i] - (*m_fineA) * m_u[i];
        } else {
          residual = m_r[i] - m_A[i - 1] * m_u[i];
        }
      }
      // restrict residual to coarser level
      // r <- R^j-1 (r - A u^j)
      // auto &r_i = m_r[i];
      restrict_level_up(residual, m_r[i + 1], i);
      LOG(3, "\t\trestricted to level "
                 << level + 1 << ". residual_norm = " << residual.norm()
                 << " m_r[level+1]_norm = " << m_r[i + 1].norm());

      // v_cycle up the hierarchy
      m_u[i + 1].setZero(m_particles[i].size());

      v_cycle(m_u, m_r, i + 1);

      // interpolate result to finer level
      // u(1) <-- u(1) + Rt*u(0)
      interpolate_level_down(m_u[i + 1], m_u[i], i + 1);
      // std::cout << "level " << i << " u norm = " << m_u[i].norm() <<
      // std::endl;

      // update residual and perform post-smoothing
      // r <- r - Au
      // u <- u + A-1 r
      // for (int j = 0; j < std::pow(2, level); ++j) {
      for (int j = 0; j < 1; ++j) {

        LOG(3, "\t\tpost-smoothing on level " << level);

        if (i == 0) {
          residual = m_r[i] - (*m_fineA) * m_u[i];
        } else {
          residual = m_r[i] - m_A[i - 1] * m_u[i];
        }

        solve_level(residual, m_u[i], i);

        /*
        std::cout << "post-smoothing level " << i
                  << " residual norm = " << residual.norm() << std::endl;
        std::cout << "level " << i << " u norm = " << m_u[i].norm()
                  << std::endl;
                  */
      }
    }
  }
};

} // namespace Aboria

#endif // HAVE_EIGEN
#endif
