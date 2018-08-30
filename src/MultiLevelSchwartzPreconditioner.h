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

#ifndef MULTILEVELSCHWARTZPRECONDITIONER_H_
#define MULTILEVELSCHWARTZPRECONDITIONER_H_

#ifdef HAVE_EIGEN

namespace Aboria {

template <typename Operator, typename Solver>
class MultiLevelSchwartzPreconditioner {
  typedef double Scalar;
  typedef size_t Index;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> matrix_type;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> vector_type;
  typedef Solver solver_type;

  typedef detail::storage_vector_type storage_vector_type;

  typedef std::vector<std::vector<storage_vector_type>> connectivity_t;
  typedef std::vector<std::vector<solver_type>> solvers_t;
  typedef std::vector<std::vector<matrix_type>> matrices_t;

  typedef typename Operator::FirstBlock kernel_t;
  static_assert(
      std::is_same<typename kernel_t::row_elements_type,
                   typename kernel_t::col_elements_type>::value,
      "Multi-Level Schwartz preconditioner restricted to identical row and col "
      "particle sets");
  static_assert(
      kernel_t::BlockRows == 1 && kernel_t::BlockCols == 1,
      "Multi-Level Schwartz preconditioner not currently implemented for "
      "vector-valued kernels");
  typedef kernel_t::row_elements_type particles_t;
  typedef MatrixReplacement<1, 1, std::tuple<kernel_t>> operator_t;
  typedef typedef std::vector<operator_t> level_operators_t;
  typedef std::vector<storage_vector_type> level_connectivity_t;
  typedef std::vector<particles_t> level_particles_t;

protected:
  bool m_isInitialized;

private:
  size_t m_max_buffer_n;
  int m_multiplicitive;

  connectivity_t m_indicies;
  connectivity_t m_buffer;
  solvers_t m_factorized_matrix;
  matrices_t m_matrix;
  level_operators_t m_A;
  const operator_t *m_fineA;
  level_connectivity_t m_restriction;
  level_particles_t m_particles;

  mutable std::vector<vector_type> m_r;
  mutable std::vector<vector_type> m_u;

  Index m_rows;
  Index m_cols;

public:
  typedef typename vector_type::StorageIndex StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic
  };

  MultiLevelSchwartzPreconditioner()
      : m_isInitialized(false), m_max_buffer_n(300), m_multiplicitive(0) {}

  template <typename MatType>
  explicit MultiLevelSchwartzPreconditioner(const MatType &mat) {
    compute(mat);
  }

  Index rows() const { return m_rows; }
  Index cols() const { return m_cols; }

  void set_max_buffer_n(size_t arg) { m_max_buffer_n = arg; }
  void set_multiplicative(int arg) { m_multiplicitive = arg; }

  template <typename MatType>
  MultiLevelSchwartzPreconditioner &analyzePattern(const MatType &mat) {
    LOG(2, "MultiLevelSchwartzPreconditioner: analyzePattern: do nothing");
    return *this;
  }

  template <typename T, typename Query> struct process_node {
    static const unsigned int D = Query::dimension;
    typedef Vector<double, D> double_d;
    typedef position_d<D> position;
    typename Query::child_iterator *m_nodes;
    storage_vector_type *m_domain_indicies;
    storage_vector_type *m_domain_buffer;
    matrix_type *m_domain_matrix;
    size_t m_max_buffer_n;
    size_t m_max_bucket_n;
    const T m_function;
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

      const double_d middle = 0.5 * (bounds.bmax + bounds.bmin);
      const double_d side = 0.9 * (bounds.bmax - bounds.bmin);

      // skip over empty buckets
      if (m_query.get_bucket_particles(*ci)) {
        return;
      }

      // find out how many particles and how many buffer particles there are
      size_t n_indicies = 0;
      size_t n_buffer = 0;

      if (is_root) {
        n_buffer = 0;
        n_indicies = m_query.number_of_particles();
      } else {
        for (auto bucket =
                 m_query.template get_buckets_near_point<-1>(middle, side);
             bucket != false; ++bucket) {
          if (bucket.get_child_iterator() == ci) {
            n_indicies +=
                m_query.get_bucket_particles(*bucket).distance_to_end();
          } else {
            n_buffer += m_query.get_bucket_particles(*bucket).distance_to_end();
          }
        }
      }

      // allocate for buffer + indicies
      storage_vector_type buffer_tmp;
      storage_vector_type indicies_tmp;
      buffer_tmp.resize(n_buffer);
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
        for (auto bucket =
                 m_query.template get_buckets_near_point<-1>(middle, side);
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
              buffer_tmp[i_buffer++] = index;
            }
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
          shuffle(buffer_tmp);
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

      LOG(3, "\tafter filtering, found  " << indicies.size()
                                          << " indices in bucket and "
                                          << buffer.size() << " in buffer");

      // ASSERT(buffer.size() > 0, "no particles in buffer");
      ASSERT_CUDA(indicies.size() > 0);

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

  struct factorize_matrix {
    solver_type operator()(matrix_type &domain_matrix) {
      LOG(3, "factorize matrix with size (" << domain_matrix.rows() << ','
                                            << domain_matrix.cols() << ')');

      solver_type solver;
      if (domain_matrix.rows() == 0 || domain_matrix.cols() == 0) {
        return solver;
      }

      solver.compute(domain_matrix);

      Eigen::VectorXd b = Eigen::VectorXd::Random(domain_matrix.rows());
      Eigen::VectorXd x = solver.solve(b);
      double relative_error = (domain_matrix * x - b).norm() / b.norm();
      if (relative_error > 1e-3 || std::isnan(relative_error)) {
        std::cout << "relative error = " << relative_error << std::endl;
      }
      return solver;
    }
  };

  void factorize_impl_block() {
    typedef typename particles_t::query_type query_type;
    typedef typename query_type::traits_type traits_type;

    const kernel_t &kernel = m_fineA->get_first_kernel();
    const particles_t &a = kernel.get_row_elements();

    CHECK(&a == &(kernel.get_col_elements()),
          "Schwartz preconditioner restricted to identical row and col "
          "particle "
          "sets");

    const int n_levels =
        (std::log(a.size()) - std::log(m_max_buffer_n)) / std::log(2) + 1;
    LOG(2,
        "MultiLevelSchwartzPreconditioner: creating " << n_levels << " levels");
    m_particles.resize(n_levels - 1);
    m_restriction.resize(n_levels - 1);

    for (int i = 0; i < n_levels - 1; ++i) {
      if (i == 1) {
        m_A.emplace_back(*m_fineA);
      } else {
        m_A.emplace_back(m_A[i - 2]);
      }
      // setup particles on this level
      const int n_particles = a.size() / std::pow(2, i);
      auto &particles = m_particles[i - 1];
      particles.resize(n_particles);
      detail::tabulate(particles.begin(), particles.end(),
                       [prev_level = iterator_to_raw_pointer(
                            (i == 1 ? a.begin() : m_particles[i - 2].begin()))](
                           const int index) { return prev_level[index * 2]; });

      const int max_bucket_size =
          i == n_levels - 2 ? m_max_buffer_n : a.get_max_bucket_size();
      particles.init_neighbour_search(a.get_min(), a.get_max(),
                                      a.get_periodic(), max_bucket_size);

      // setup A operator on this level
    }

    m_indicies.resize(n_levels);
    m_buffer.resize(n_levels);
    m_matrix.resize(n_levels);
    m_factorized_matrix.resize(n_levels);
    m_r.resize(m_levels);
    m_u.resize(m_levels);

    for (int i = 0; i < n_levels; ++i) {
      // get list of leaf cells
      const query_type &query =
          i == 0 ? a.get_query() : m_particles[i - 1].get_query();
      auto df_search = query.breadth_first();
      for (; df_search != false; ++df_search) {
      }
      auto &nodes = *df_search;
      if (i == n_levels - 1) {
        CHECK(nodes.size() == 1 && query.is_leaf_node(nodes[0]),
              "top level should have one leaf node");
      }

      auto &indicies = m_indicies[i];
      auto &buffer = m_buffer[i];
      auto &matrix = m_matrix[i];
      auto &factorized_matrix = m_factorized_matrix[i];

      indicies.resize(nodes.size());
      buffer.resize(nodes.size());
      matrix.resize(nodes.size());
      factorized_matrix.resize(nodes.size());

      if (traits_type::data_on_GPU) {
        // need to copy data from/to gpu
        typename traits_type::template vector<storage_vector_type> tmp_indicies(
            nodes.size());
        typename traits_type::template vector<storage_vector_type> tmp_buffer(
            nodes.size());
        typename traits_type::template vector<matrix_type> tmp_matrix(
            nodes.size());
        typename traits_type::template vector<matrix_type> tmp_coupling_matrix(
            coupling_matrix.size());

        auto count = traits_type::make_counting_iterator(0);
        detail::for_each(
            count, count + nodes.size(),
            process_node<typename Kernel::function_type, query_type>{
                iterator_to_raw_pointer(nodes.begin()),
                iterator_to_raw_pointer(tmp_indicies.begin()),
                iterator_to_raw_pointer(tmp_buffer.begin()),
                iterator_to_raw_pointer(tmp_matrix.begin()), m_max_buffer_n,
                a.get_max_bucket_size(), kernel.get_kernel_function(), query, i,
                nodes.size() == 1});

        detail::copy(tmp_indicies.begin(), tmp_indicies.end(),
                     indicies.begin());
        detail::copy(tmp_buffer.begin(), tmp_buffer.end(), buffer.begin());
        detail::copy(tmp_matrix.begin(), tmp_matrix.end(), matrix.begin());
      } else {
        // no copy required
        auto count = traits_type::make_counting_iterator(0);
        detail::for_each(
            count, count + nodes.size(),
            process_node<typename Kernel::function_type, query_type>{
                iterator_to_raw_pointer(nodes.begin()),
                iterator_to_raw_pointer(indicies.begin()),
                iterator_to_raw_pointer(buffer.begin()),
                iterator_to_raw_pointer(matrix.begin()), m_max_buffer_n,
                a.get_max_bucket_size(), kernel.get_kernel_function(), query, i,
                nodes.size() == 1});
      }

      // factorize domain matrices
      detail::transform(matrix.begin(), matrix.end(), factorized_matrix.begin(),
                        factorize_matrix());
    }
  }

  template <unsigned int NI, unsigned int NJ, typename Blocks>
  MultiLevelSchwartzPreconditioner &
  factorize(const MatrixReplacement<NI, NJ, Blocks> &mat) {
    LOG(2, "MultiLevelSchwartzPreconditioner: factorize");
    m_rows = mat.rows();
    m_cols = mat.cols();
    m_fineA = &mat.get_first_kernel();
    factorize_impl_block();

    for (int i = 0; i < m_indicies.size(); ++i) {
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

      auto pA = i == 0 ? m_fineA : &m_A[i - 1];
      auto &particles =
          i == 0 ? m_fineA->get_row_elements() : m_particles[i - 1];

      LOG(2, "MultiLevelSchwartzPreconditioner: finished factorizing, on level "
                 << i << " found " << m_indicies[i].size() << " domains, with "
                 << min_indicies << "--" << max_indicies << " particles ("
                 << count << " total), and " << min_buffer << "--" << max_buffer
                 << " buffer particles. particles size is " << particles.size()
                 << ". kernel matrix is size (" << pA->rows() << ','
                 << pA->cols() << ")");
    }

    m_isInitialized = true;
    return *this;
  }

  template <int _Options, typename _StorageIndex, int RefOptions,
            typename RefStrideType>
  MultiLevelSchwartzPreconditioner &
  factorize(const Eigen::Ref<
            const Eigen::SparseMatrix<Scalar, _Options, _StorageIndex>,
            RefOptions, RefStrideType> &mat) {
    CHECK(m_indicies[0].size() > 0,
          "MultiLevelSchwartzPreconditioner::factorize(): cannot factorize "
          "sparse "
          "matrix, call factorize with a Aboria MatrixReplacement class "
          "instead");
    return *this;
  }

  template <typename Derived>
  MultiLevelSchwartzPreconditioner &
  factorize(const Eigen::DenseBase<Derived> &mat) {
    CHECK(m_indicies[0].size() > 0,
          "MultiLevelSchwartzPreconditioner::analyzePattern(): cannot "
          "factorize dense "
          "matrix, call factorize with a Aboria MatrixReplacement class "
          "instead");
    return *this;
  }

  template <typename MatType>
  MultiLevelSchwartzPreconditioner &compute(const MatType &mat) {
    analyzePattern(mat);
    return factorize(mat);
  }

  template <typename Rhs, typename Dest> struct solve_domain {
    const storage_vector_type *m_domain_indicies;
    const storage_vector_type *m_domain_buffer;
    const solver_type *domain_factorized_matrix;
    const matrix_type *m_coupling_matrix;
    Dest &x;
    const Rhs &b;
    int m_level;
    bool is_root;

    void operator()(const int i) const {
      const storage_vector_type &buffer = m_domain_buffer[i];
      const storage_vector_type &indicies = m_domain_indicies[i];
      if (indicies.size() == 0)
        return;

      const size_t nb = buffer.size();

      vector_type domain_x;
      vector_type domain_b;
      domain_x.resize(nb);
      domain_b.resize(nb);

      // copy b values from big vector
      for (size_t j = 0; j < buffer.size(); ++j) {
        domain_b[j] = b[buffer[j]];
      }

      // solve domain
      domain_x = domain_factorized_matrix[i].solve(domain_b);

      // copy accumulate x values to big vector
      for (size_t j = 0; j < indicies.size(); ++j) {
        x[indicies[j]] += domain_x[j];
      }
      /*
       // non-restricted
      for (size_t j = 0; j < buffer.size(); ++j) {
        x[buffer[j]] += domain_x[sub_index++];
      }
      */
    }
  };

  /** \internal */
  template <typename Rhs, typename Dest>
  void _solve_impl(const Rhs &b, Dest &x) const {

    vector_type tmp;
    if (m_multiplicitive == 2 && m_levels > 1) {
      m_r[0] = b;
      m_u[0] = vector_type::Zero(a.size());
      // symmetric multiplicative
      for (int i = 0; i < m_indicies.size(); ++i) {
        auto count = boost::make_counting_iterator(0);
        detail::for_each(
            count, count + m_indicies[i].size(),
            solve_domain<Rhs, Dest>{
                iterator_to_raw_pointer(m_indicies[i].begin()),
                iterator_to_raw_pointer(m_buffer[i].begin()),
                iterator_to_raw_pointer(m_factorized_matrix[i].begin()), m_u[i],
                m_r[i], i, m_indicies[i].size() == 1});
        if (i < m_indicies.size() - 1) {
          if (i == 0) {
            tmp = m_r[0] - (*m_fineA) * u[0]
          } else {
            tmp = m_r[i] - m_A[i - 1] * u[i];
          }
          // restrict to coarser level
          m_r[i + 1].resize(m_particles[i].size());
          detail::tabulate(m_r[i + 1].data(),
                           m_r[i + 1].data() + m_r[i + 1].size(),
                           [&tmp, id = get<id>(m_particles[i].begin())](
                               const int index) { return tmp[2 * id[index]]; });
          m_u[i + 1] = vector_type::Zero(m_particles[i].size());
        }
      }
      for (int i = m_indicies.size() - 2; i >= 0; ++i) {
        // u(1) <-- u(1) + Rt*u(0)
        detail::for_each(count, count + m_particles[i].size(),
                         [&m_u, id = get<id>(m_particles[i].begin())](
                             const int upper_index) {
                           const int lower_index = 2 * id[upper_index];
                           m_u[i][lower_index] += m_u[i + 1][upper_index];
                         });

        m_r[i] -= m_A[i] * m_u[i];

        auto count = boost::make_counting_iterator(0);
        detail::for_each(
            count, count + m_indicies[i].size(),
            solve_domain<Rhs, Dest>{
                iterator_to_raw_pointer(m_indicies[i].begin()),
                iterator_to_raw_pointer(m_buffer[i].begin()),
                iterator_to_raw_pointer(m_factorized_matrix[i].begin()), m_u[i],
                m_r[i], i, m_indicies[i].size() == 1});
      }
      x = m_u[0];

    } else if (m_multiplicitive == 1 && m_levels > 1) {
      // non-symmetric multiplicative
      x.setZero();
      vector_type b_copy = b;

      for (int i = 0; i < std::min(m_levels, 2); ++i) {
        solve_domain<Rhs, Dest> fsolve_domain{
            iterator_to_raw_pointer(m_indicies[i].begin()),
            iterator_to_raw_pointer(m_buffer[i].begin()),
            iterator_to_raw_pointer(m_factorized_matrix[i].begin()),
            iterator_to_raw_pointer(m_coupling_matrix[i].begin()),
            x,
            b_copy,
            i,
            m_indicies[i].size() == 1,
            m_interpolate};

        auto count = boost::make_counting_iterator(0);
        detail::for_each(count, count + m_indicies[i].size(), fsolve_domain);

        if (i == 0) {
          b_copy -= m_A * x;
        }
      }
    } else {
      x.setZero();

      for (int i = 0; i < std::min(m_levels, 2); ++i) {
        solve_domain<Rhs, Dest> fsolve_domain{
            iterator_to_raw_pointer(m_indicies[i].begin()),
            iterator_to_raw_pointer(m_buffer[i].begin()),
            iterator_to_raw_pointer(m_factorized_matrix[i].begin()),
            iterator_to_raw_pointer(m_coupling_matrix[i].begin()),
            x,
            b,
            i,
            m_indicies[i].size() == 1,
            m_interpolate};

        auto count = boost::make_counting_iterator(0);
        detail::for_each(count, count + m_indicies[i].size(), fsolve_domain);
      }
    }
  }

  template <typename Rhs>
  inline const Eigen::Solve<MultiLevelSchwartzPreconditioner, Rhs>
  solve(const Eigen::MatrixBase<Rhs> &b) const {
    eigen_assert(static_cast<typename Rhs::Index>(m_rows) == b.rows() &&
                 "MultiLevelSchwartzPreconditioner::solve(): invalid number of "
                 "rows of the "
                 "right hand side matrix b");
    eigen_assert(m_isInitialized &&
                 "MultiLevelSchwartzPreconditioner is not initialized.");
    return Eigen::Solve<MultiLevelSchwartzPreconditioner, Rhs>(*this,
                                                               b.derived());
  }

  Eigen::ComputationInfo info() { return Eigen::Success; }
}; //

} // namespace Aboria

#endif // HAVE_EIGEN
#endif
