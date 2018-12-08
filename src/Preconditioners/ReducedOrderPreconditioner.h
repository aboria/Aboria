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

#include "H2Lib.h"

namespace Aboria {

#ifdef HAVE_H2LIB
template <typename Solver> class ReducedOrderPreconditioner {
  typedef Solver solver_type;
  typedef size_t Index;
  typedef double Scalar;
  typedef H2LibMatrix h2_matrix_type;
  Index m_rows;
  Index m_cols;
  double m_tol;
  std::vector<size_t> m_col_sizes;
  std::vector<size_t> m_row_sizes;
  std::vector<std::shared_ptr<solver_type>> m_solvers;
  std::vector<const h2_matrix_type *> m_h2mats;

public:
  typedef size_t StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic
  };

  ReducedOrderPreconditioner() : m_tol(1e-5) {}

  template <typename MatType>
  explicit ReducedOrderPreconditioner(const MatType &mat) : m_tol(1e-5) {
    compute(mat);
  }

  Index rows() const { return m_rows; }
  Index cols() const { return m_cols; }

  void set_tolerance(const double tol) { m_tol = tol; }

  template <unsigned int NI, unsigned int NJ, typename Blocks>
  ReducedOrderPreconditioner &
  analyzePattern(const MatrixReplacement<NI, NJ, Blocks> &mat) {

    LOG(2, "ExtMatrixPreconditioner: analyze pattern");
    m_rows = mat.rows();
    m_cols = mat.cols();

    return *this;
  }

  struct factorize_block {
    std::vector<size_t> &m_col_sizes;
    std::vector<size_t> &m_row_sizes;
    std::vector<std::shared_ptr<solver_type>> &m_solvers;
    std::vector<const h2_matrix_type *> &m_h2mats;
    double m_tol;
    int i;

    factorize_block(std::vector<size_t> &col_sizes,
                    std::vector<size_t> &row_sizes,
                    std::vector<std::shared_ptr<solver_type>> &solvers,
                    std::vector<const h2_matrix_type *> &h2mats, double tol)
        : m_col_sizes(col_sizes), m_row_sizes(row_sizes), m_solvers(solvers),
          m_h2mats(h2mats), m_tol(tol), i(0) {}

    template <typename Block> void operator()(const Block &block) {
      LOG(2, "ReducedOrderPreconditioner: block " << i << ": non h2 block");
      m_solvers[i] = nullptr;
      m_col_sizes[i] = block.cols();
      m_row_sizes[i] = block.rows();
      ++i;
    }

    template <typename RowParticles, typename ColParticles, typename PositionF,
              typename F>
    void operator()(
        const KernelH2<RowParticles, ColParticles, PositionF, F> &kernel) {

      m_h2mats[i] = &kernel.get_h2_matrix();
      m_col_sizes[i] = kernel.cols();
      m_row_sizes[i] = kernel.rows();

      LOG(2, "ReducedOrderPreconditioner: block "
                 << i << ": factorise h2 matrix with tolerance " << m_tol);
      m_solvers[i] = std::make_shared<solver_type>(
          m_h2mats[i]->get_ph2matrix(), m_h2mats[i]->get_pblock(), m_tol);
      std::vector<double> b(m_col_sizes[i]);
      std::vector<double> b2(m_row_sizes[i]);
      for (size_t ii = 0; ii < m_col_sizes[i]; ++ii) {
        b[ii] = 1.0;
      }
      m_h2mats[i]->matrix_vector_multiply(b2, 1, false, b);
      m_solvers[i]->solve(b2, b2);
      double sum = 0;
      double sum2 = 0;
      for (size_t ii = 0; ii < m_row_sizes[i]; ++ii) {
        sum += std::pow(b2[ii] - b[ii], 2);
        sum2 += std::pow(b[ii], 2);
      }
      LOG(2, "ReducedOrderPreconditioner: block "
                 << i << ": factorisation accuracy: " << std::sqrt(sum / sum2));

      // m_solvers[i]->setMaxIterations(m_inner_iterations);
      // LOG(2,"ExtMatrixPreconditioner: block "<<i<<": set precon");
      // m_solvers[i]->preconditioner().setDroptol(0.1);
      // m_solvers[i]->preconditioner().compute(m_str_ext_matrices[i]);
      ++i;
    }
  };

  template <unsigned int NI, unsigned int NJ, typename Blocks>
  ReducedOrderPreconditioner &
  factorize(const MatrixReplacement<NI, NJ, Blocks> &mat) {
    LOG(2, "ReducedOrderPreconditioner: factorizing domain");

    m_rows = mat.rows();
    m_cols = mat.cols();
    m_solvers.resize(NI);
    m_h2mats.resize(NI);
    m_col_sizes.resize(NI);
    m_row_sizes.resize(NI);
    detail::apply_function_to_diagonal_blocks(
        factorize_block(m_col_sizes, m_row_sizes, m_solvers, m_h2mats, m_tol),
        mat);
    m_isInitialized = true;

    return *this;
  }

  template <typename MatType>
  ReducedOrderPreconditioner &compute(const MatType &mat) {
    analyzePattern(mat);
    return factorize(mat);
  }

  template <typename Rhs, typename Dest>
  void _solve_impl(const Rhs &b, Dest &x) const {
    size_t row = 0;
    size_t col = 0;
    Eigen::Matrix<double, Eigen::Dynamic, 1> x_buffer;

    for (size_t i = 0; i < m_solvers.size(); ++i) {
      auto b_segment = b.segment(col, m_col_sizes[i]);
      auto x_segment = x.segment(row, m_row_sizes[i]);
      if (m_solvers[i] != nullptr) { // solver only exists for h2 blocks
        LOG(2, "ReducedOrderPreconditioner: block "
                   << i << " solve for " << m_row_sizes[i] << "x"
                   << m_col_sizes[i] << " matrix");
        x_buffer.resize(m_col_sizes[i]);
        m_solvers[i]->solve(b_segment, x_buffer);
        x_segment = x_buffer;
      } else {
        x_segment = b_segment;
      }
      row += m_row_sizes[i];
      col += m_col_sizes[i];
    }
    LOG(2, "ReducedOrderPreconditioner: done solve_impl");
  }

  template <typename Rhs>
  inline const Eigen::Solve<ReducedOrderPreconditioner, Rhs>
  solve(const Eigen::MatrixBase<Rhs> &b) const {
    eigen_assert(m_rows == b.rows() &&
                 "ReducedOrderPreconditioner::solve(): invalid number of rows "
                 "of the right hand side matrix b");
    eigen_assert(m_isInitialized &&
                 "ReducedOrderPreconditioneris not initialized.");
    return Eigen::Solve<ReducedOrderPreconditioner, Rhs>(*this, b.derived());
  }

  Eigen::ComputationInfo info() { return Eigen::Success; }

protected:
  bool m_isInitialized;
};
#endif // HAVE_H2LIB

} // namespace Aboria

#endif
