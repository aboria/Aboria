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

#ifndef REDUCED_ORDER_PRECONDITIONER_H_
#define REDUCED_ORDER_PRECONDITIONER_H_

#include "detail/Preconditioners.h"

namespace Aboria {

#ifdef HAVE_EIGEN

template <typename Solver> class NystromPreconditioner {
  typedef double Scalar;
  typedef size_t Index;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> matrix_type;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> vector_type;
  typedef Solver solver_type;
  typedef std::vector<size_t> storage_vector_type;
  typedef std::vector<storage_vector_type> connectivity_type;

protected:
  bool m_isInitialized;

private:
  size_t m_random;
  double m_lambda;

  std::vector<storage_vector_type> m_domain_indicies;
  std::vector<solver_type> m_domain_factorized_matrix;
  std::vector<matrix_type> m_domain_Kux;
  std::vector<vint2> m_domain_range;
  Index m_rows;
  Index m_cols;

public:
  typedef typename vector_type::StorageIndex StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic
  };

  NystromPreconditioner()
      : m_isInitialized(false), m_random(0), m_lambda(1e-8) {}

  template <typename MatType>
  explicit NystromPreconditioner(const MatType &mat) {
    compute(mat);
  }

  Index rows() const { return m_rows; }
  Index cols() const { return m_cols; }

  void set_number_of_random_particles(size_t n) { m_random = n; }
  void set_lambda(double val) { m_lambda = val; }

  template <typename Kernel>
  void analyze_impl_block(const Index start_row, const Kernel &kernel) {
    typedef typename Kernel::row_elements_type row_elements_type;
    typedef typename Kernel::col_elements_type col_elements_type;

    static_assert(std::is_same<row_elements_type, col_elements_type>::value,
                  "Nystrom preconditioner restricted to identical row and col "
                  "particle sets");
    const row_elements_type &a = kernel.get_row_elements();
    CHECK(&a == &(kernel.get_col_elements()),
          "Nystrom preconditioner restricted to identical row and col "
          "particle "
          "sets");

    const size_t domain_index = m_domain_indicies.size();
    m_domain_indicies.push_back(connectivity_type::value_type());
    m_domain_range.push_back(vint2(start_row, start_row + a.size()));
    m_domain_Kux.push_back(matrix_type());
    storage_vector_type &indicies = m_domain_indicies[domain_index];

    if (m_random >= a.size()) {
      // add all indicies
      indicies.resize(a.size());
      std::iota(indicies.begin(), indicies.end(), 0);
    } else {
      // add some random indicies
      std::uniform_int_distribution<int> uniform_index(0, a.size() - 1);
      std::default_random_engine generator;

      for (size_t d = 0; d < m_random; ++d) {
        bool in_indicies;
        size_t proposed_index;
        do {
          proposed_index = uniform_index(generator) + start_row;

          // check not in indicies
          in_indicies =
              indicies.end() !=
              std::find(indicies.begin(), indicies.end(), proposed_index);

        } while (in_indicies);
        indicies.push_back(proposed_index);
      }
    }
    ASSERT(indicies.size() > 0, "no particles in domain");
  }

  template <typename RowParticles, typename ColParticles>
  void
  analyze_impl_block(const Index start_row,
                     const KernelZero<RowParticles, ColParticles> &kernel) {}

  template <unsigned int NI, unsigned int NJ, typename Blocks, std::size_t... I>
  void analyze_impl(const MatrixReplacement<NI, NJ, Blocks> &mat,
                    detail::index_sequence<I...>) {
    int dummy[] = {0, (analyze_impl_block(mat.template start_row<I>(),
                                          std::get<I * NJ + I>(mat.m_blocks)),
                       0)...};
    static_cast<void>(dummy);
  }

  template <unsigned int NI, unsigned int NJ, typename Blocks>
  NystromPreconditioner &
  analyzePattern(const MatrixReplacement<NI, NJ, Blocks> &mat) {

    LOG(2, "NystromPreconditioner: analyze pattern");
    m_rows = mat.rows();
    m_cols = mat.cols();
    analyze_impl(mat, detail::make_index_sequence<NI>());

    int count = 0;
    int minsize_indicies = std::numeric_limits<int>::max();
    int maxsize_indicies = std::numeric_limits<int>::min();
    for (size_t domain_index = 0; domain_index < m_domain_indicies.size();
         ++domain_index) {
      const int size_indicies = m_domain_indicies[domain_index].size();
      count += size_indicies;
      if (size_indicies < minsize_indicies)
        minsize_indicies = size_indicies;
      if (size_indicies > maxsize_indicies)
        maxsize_indicies = size_indicies;
    }
    LOG(2, "NystromPreconditioner: finished analysis, found "
               << m_domain_indicies.size() << " domains, with "
               << minsize_indicies << "--" << maxsize_indicies << " particles ("
               << count << " total)")
    return *this;
  }

  template <int _Options, typename _StorageIndex>
  NystromPreconditioner &analyzePattern(
      const Eigen::SparseMatrix<Scalar, _Options, _StorageIndex> &mat) {
    CHECK(m_domain_indicies.size() > 0,
          "NystromPreconditioner::analyzePattern(): cannot analyze sparse "
          "matrix, "
          "call analyzePattern using a Aboria MatrixReplacement class first");
    return *this;
  }

  template <int _Options, typename _StorageIndex, int RefOptions,
            typename RefStrideType>
  NystromPreconditioner &
  analyzePattern(const Eigen::Ref<
                 const Eigen::SparseMatrix<Scalar, _Options, _StorageIndex>,
                 RefOptions, RefStrideType> &mat) {
    CHECK(m_domain_indicies.size() > 0,
          "NystromPreconditioner::analyzePattern(): cannot analyze sparse "
          "matrix, "
          "call analyzePattern using a Aboria MatrixReplacement class first");
    return *this;
  }

  template <typename Derived>
  NystromPreconditioner &analyzePattern(const Eigen::DenseBase<Derived> &mat) {
    CHECK(m_domain_indicies.size() > 0,
          "NystromPreconditioner::analyzePattern(): cannot analyze dense "
          "matrix, "
          "call analyzePattern need to pass a Aboria MatrixReplacement class "
          "first");
    return *this;
  }

  template <typename MatType>
  NystromPreconditioner &factorize(const MatType &mat) {
    LOG(2, "NystromPreconditioner: factorizing domain");
    eigen_assert(
        static_cast<typename MatType::Index>(m_rows) == mat.rows() &&
        "SchwartzPreconditioner::solve(): invalid number of rows of mat");
    eigen_assert(
        static_cast<typename MatType::Index>(m_cols) == mat.cols() &&
        "SchwartzPreconditioner::solve(): invalid number of rows of mat");

    m_domain_factorized_matrix.resize(m_domain_indicies.size());

    matrix_type Kuu;

    for (size_t domain_index = 0;
         domain_index < m_domain_factorized_matrix.size(); ++domain_index) {
      const storage_vector_type &indicies = m_domain_indicies[domain_index];
      solver_type &solver = m_domain_factorized_matrix[domain_index];
      vint2 &range = m_domain_range[domain_index];
      matrix_type &Kux = m_domain_Kux[domain_index];

      const size_t size = indicies.size();
      // std::cout << "domain "<<domain_index<<"indicies =
      // "<<indicies.size()<<" buffer =  "<<buffer.size()<<" random =
      // "<<random.size()<<std::endl;

      Kuu.resize(size, size);
      Kux.resize(size, range[1] - range[0]);

      size_t i = 0;
      for (const size_t &big_index_i : indicies) {
        size_t j = 0;
        for (const size_t &big_index_j : indicies) {
          Kuu(i, j++) = mat.coeff(big_index_i, big_index_j);
        }
        j = 0;
        for (int big_index_j = range[0]; big_index_j < range[1];
             ++big_index_j) {
          Kux(i, j++) = mat.coeff(big_index_i, big_index_j);
        }
        ++i;
      }
      Kuu += Kux * (Kux.transpose());
      solver.compute(Kuu);

      Eigen::VectorXd b = Eigen::VectorXd::Random(Kuu.rows());
      Eigen::VectorXd x = solver.solve(b);
      double relative_error = (Kuu * x - b).norm() / b.norm();
      if (relative_error > 1e-3 || std::isnan(relative_error)) {
        std::cout << "relative error = " << relative_error << std::endl;
      }
    }

    m_isInitialized = true;

    return *this;
  }

  template <typename MatType>
  NystromPreconditioner &compute(const MatType &mat) {
    analyzePattern(mat);
    return factorize(mat);
  }

  /** \internal */
  template <typename Rhs, typename Dest>
  void _solve_impl(const Rhs &b, Dest &x) const {
    x = b;
    for (size_t i = 0; i < m_domain_indicies.size(); ++i) {
      auto range = m_domain_range[i];
      auto Kux = m_domain_Kux[i];
      auto solver = m_domain_factorized_matrix[i];

      x.segment(range[0], range[1]) =
          (1.0 / m_lambda) *
          (b.segment(range[0], range[1]) -
           (Kux.transpose()) *
               solver.solve(Kux * b.segment(range[0], range[1])));
    }
  }

  template <typename Rhs>
  inline const Eigen::Solve<NystromPreconditioner, Rhs>
  solve(const Eigen::MatrixBase<Rhs> &b) const {
    eigen_assert(
        static_cast<typename Rhs::Index>(m_rows) == b.rows() &&
        "NystromPreconditioner::solve(): invalid number of rows of the "
        "right hand side matrix b");
    eigen_assert(m_isInitialized &&
                 "NystromPreconditioner is not initialized.");
    return Eigen::Solve<NystromPreconditioner, Rhs>(*this, b.derived());
  }

  Eigen::ComputationInfo info() { return Eigen::Success; }
};
#endif
} // namespace Aboria

#endif
