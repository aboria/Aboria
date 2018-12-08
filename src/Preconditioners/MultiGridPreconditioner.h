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

#ifndef MULTIGRID_PRECONDITIONER_H_
#define MULTIGRID_PRECONDITIONER_H_

#include "Preconditioners/SchwartzPreconditioner.h"

#ifdef HAVE_EIGEN
namespace Aboria {

template <typename Operator, typename Solver> class MultiGridPreconditioner {
  typedef double Scalar;
  typedef size_t Index;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> matrix_type;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> vector_type;
  typedef Solver solver_type;
  typedef SchwartzDecomposition<Operator, Solver> decomposition_t;

protected:
  bool m_isInitialized;

private:
  mutable std::vector<vector_type> m_r;
  mutable std::vector<vector_type> m_u;

  Index m_rows;
  Index m_cols;

  mutable decomposition_t m_decom;

public:
  typedef typename vector_type::StorageIndex StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic
  };

  MultiGridPreconditioner() : m_isInitialized(false), m_decom() {}

  template <typename MatType>
  explicit MultiGridPreconditioner(const MatType &mat) {
    compute(mat);
  }

  Index rows() const { return m_rows; }
  Index cols() const { return m_cols; }

  void set_max_buffer_n(size_t arg) { m_decom.set_max_buffer_n(arg); }

  void set_mult_buffer(size_t arg) { m_decom.set_mult_buffer(arg); }
  void set_smoother_weighting(double arg) {
    m_decom.set_smoother_weighting(arg);
  }

  template <typename MatType>
  MultiGridPreconditioner &analyzePattern(const MatType &mat) {
    LOG(2, "MultiGridPreconditioner: analyzePattern: do nothing");
    return *this;
  }

  template <unsigned int NI, unsigned int NJ, typename Blocks>
  MultiGridPreconditioner &
  factorize(const MatrixReplacement<NI, NJ, Blocks> &mat) {
    LOG(2, "MultiGridPreconditioner: factorize");
    m_rows = mat.rows();
    m_cols = mat.cols();

    m_decom.set_operator(mat);
    const bool set_interpolate = false;
    m_decom.hierarchical_schwartz_decomposition(std::numeric_limits<int>::max(),
                                                set_interpolate);
    // m_decom.hierarchical_schwartz_decomposition(2, set_interpolate);

    const size_t n_levels = m_decom.get_n_levels();
    m_r.resize(n_levels);
    m_u.resize(n_levels);

    // allocate r's and u's
    const auto &particles = m_decom.get_particles();
    for (size_t i = 0; i < particles.size(); ++i) {
      // r <- R^j-1 r
      m_r[i + 1].resize(particles[i].size());
      m_u[i + 1].resize(particles[i].size());
    }

    m_isInitialized = true;
    return *this;
  }

  template <int _Options, typename _StorageIndex, int RefOptions,
            typename RefStrideType>
  MultiGridPreconditioner &
  factorize(const Eigen::Ref<
            const Eigen::SparseMatrix<Scalar, _Options, _StorageIndex>,
            RefOptions, RefStrideType> &mat) {
    CHECK(m_decom.get_n_levels() > 0,
          "MultiGridPreconditioner::factorize(): cannot factorize "
          "sparse "
          "matrix, call factorize with a Aboria MatrixReplacement class "
          "instead");
    return *this;
  }

  template <typename Derived>
  MultiGridPreconditioner &factorize(const Eigen::DenseBase<Derived> &mat) {
    CHECK(m_decom.get_n_levels() > 0,
          "MultiGridPreconditioner::analyzePattern(): cannot "
          "factorize dense "
          "matrix, call factorize with a Aboria MatrixReplacement class "
          "instead");
    return *this;
  }

  template <typename MatType>
  MultiGridPreconditioner &compute(const MatType &mat) {
    analyzePattern(mat);
    return factorize(mat);
  }

  /** \internal */
  template <typename Rhs, typename Dest>
  void _solve_impl(const Rhs &b, Dest &x) const {
    LOG(3, "Solving MultiGridPreconditioner:");

    m_u[0] = vector_type::Zero(x.size());
    m_r[0] = b;

    if (0) {
      // restrict r to coarsest level
      for (size_t i = 0; i < m_decom.get_n_levels() - 1; ++i) {
        m_decom.restrict_level_up(m_r[i], m_r[i + 1], i);
      }

      // solve directly at coarsest level
      const int top_level = m_decom.get_n_levels() - 1;
      m_u[top_level].setZero();
      m_decom.v_cycle(m_u, m_r, top_level);

      // go back down doing v-cycles to get the initial condition
      // at level 0
      for (int i = top_level - 1; i >= 0; --i) {
        m_u[i].setZero();
        m_decom.interpolate_level_down(m_u[i + 1], m_u[i], i + 1);
        m_decom.v_cycle(m_u, m_r, i);
      }
    }

    m_decom.v_cycle(m_u, m_r, 0);

    // finally, retrieve result
    x = m_u[0];

    // regenerate hierarchy
    // m_decom.hierarchical_schwartz_decomposition();

    // const_cast<MultiGridPreconditioner *>(this)
    //    ->m_decom.hierarchical_schwartz_decomposition();
  }

  template <typename Rhs>
  inline const Eigen::Solve<MultiGridPreconditioner, Rhs>
  solve(const Eigen::MatrixBase<Rhs> &b) const {
    eigen_assert(static_cast<typename Rhs::Index>(m_rows) == b.rows() &&
                 "MultiGridPreconditioner::solve(): invalid number of "
                 "rows of the "
                 "right hand side matrix b");
    eigen_assert(m_isInitialized &&
                 "MultiGridPreconditioner is not initialized.");
    return Eigen::Solve<MultiGridPreconditioner, Rhs>(*this, b.derived());
  }

  Eigen::ComputationInfo info() { return Eigen::Success; }
}; // namespace Aboria

} // namespace Aboria

#endif // HAVE_EIGEN
#endif
