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

#ifndef KERNELS_H_
#define KERNELS_H_

#include <type_traits>

#include "Particles.h"
#include "detail/Chebyshev.h"

#include "FastMultipoleMethod.h"
#include "detail/Kernels.h"

#ifdef HAVE_H2LIB
#include "H2Lib.h"
#endif

#ifdef HAVE_EIGEN
#include <Eigen/Core>

namespace Aboria {

template <typename RowElements, typename ColElements, typename F>
class KernelBase {
protected:
  typedef typename detail::kernel_helper<RowElements, ColElements, F>
      kernel_helper;
  static const unsigned int dimension = RowElements::dimension;
  typedef Vector<double, dimension> double_d;
  typedef Vector<int, dimension> int_d;
  typedef position_d<dimension> position;
  typedef double_d const &const_position_reference;
  typedef typename RowElements::const_reference const_row_reference;
  typedef typename ColElements::const_reference const_col_reference;

public:
  typedef typename kernel_helper::Block Block;
  typedef typename kernel_helper::Scalar Scalar;
  static const size_t BlockRows = kernel_helper::BlockRows;
  static const size_t BlockCols = kernel_helper::BlockCols;
  typedef typename kernel_helper::BlockRHSVector BlockRHSVector;
  typedef typename kernel_helper::BlockLHSVector BlockLHSVector;

  typedef RowElements row_elements_type;
  typedef ColElements col_elements_type;
  typedef F function_type;
  typedef size_t Index;
  enum { ColsAtCompileTime = -1, RowsAtCompileTime = -1 };

  KernelBase(const RowElements &row_elements, const ColElements &col_elements,
             const F &function)
      : m_row_elements(row_elements), m_col_elements(col_elements),
        m_function(function){};

  RowElements &get_row_elements() { return m_row_elements; }

  const RowElements &get_row_elements() const { return m_row_elements; }

  ColElements &get_col_elements() { return m_col_elements; }

  const function_type &get_kernel_function() const { return m_function; }

  const ColElements &get_col_elements() const { return m_col_elements; }

  size_t rows() const { return m_row_elements.size() * BlockRows; }

  size_t cols() const { return m_col_elements.size() * BlockCols; }

  Scalar coeff(const size_t i, const size_t j) const {
    ASSERT(i < rows(), "i greater than rows()");
    ASSERT(j < cols(), "j greater than cols()");
    const int pi = std::floor(static_cast<float>(i) / BlockRows);
    const int ioffset = i - pi * BlockRows;
    const int pj = std::floor(static_cast<float>(j) / BlockCols);
    const int joffset = j - pj * BlockCols;
    const Block block =
        Block(m_function(m_row_elements[pi], m_col_elements[pj]));
    return block(ioffset, joffset);
  }

  template <typename MatrixType> void assemble(const MatrixType &matrix) const {
    MatrixType::MATRIX_ASSEMBLE_NOT_IMPLEMENTED;
  }

  template <typename Triplet>
  void assemble(std::vector<Triplet> &triplets, const size_t startI = 0,
                const size_t startJ = 0) const {
    Triplet::TRIPLET_ASSEMBLE_NOT_IMPLEMENTED;
  }

  /// Evaluates a matrix-free linear operator given by \p expr \p if_expr,
  /// and particle sets \p a and \p b on a vector rhs and
  /// accumulates the result in vector lhs
  template <typename VectorLHS, typename VectorRHS>
  void evaluate(VectorLHS &lhs, const VectorRHS &rhs) const {
    VectorLHS::EVALUATE_NOT_IMPLEMENTED;
  }

protected:
  const RowElements &m_row_elements;
  const ColElements &m_col_elements;
  const F m_function;
};

template <typename RowElements, typename ColElements, typename F>
class KernelDense : public KernelBase<RowElements, ColElements, F> {
protected:
  typedef KernelBase<RowElements, ColElements, F> base_type;
  static const unsigned int dimension = base_type::dimension;
  typedef typename base_type::double_d double_d;
  typedef typename base_type::position position;
  typedef typename base_type::int_d int_d;
  typedef typename base_type::const_row_reference const_row_reference;
  typedef typename base_type::const_col_reference const_col_reference;

public:
  typedef typename base_type::Block Block;
  typedef typename base_type::BlockLHSVector BlockLHSVector;
  typedef typename base_type::Scalar Scalar;
  static const size_t BlockRows = base_type::BlockRows;
  static const size_t BlockCols = base_type::BlockCols;

  KernelDense(const RowElements &row_elements, const ColElements &col_elements,
              const F &function)
      : base_type(row_elements, col_elements, function){};

  template <typename Derived>
  void assemble(const Eigen::DenseBase<Derived> &matrix) const {

    const RowElements &a = this->m_row_elements;
    const ColElements &b = this->m_col_elements;
    const size_t na = a.size();
    const size_t nb = b.size();

    ASSERT(matrix.rows() == static_cast<typename Derived::Index>(this->rows()),
           "matrix has incompatible row size");
    ASSERT(matrix.cols() == static_cast<typename Derived::Index>(this->cols()),
           "matrix has incompatible col size");

    for (size_t i = 0; i < na; ++i) {
      for (size_t j = 0; j < nb; ++j) {
        const_cast<Eigen::DenseBase<Derived> &>(matrix)
            .template block<BlockRows, BlockCols>(i * BlockRows,
                                                  j * BlockCols) =
            Block(this->m_function(a[i], b[j]));
      }
    }
  }

  template <typename Triplet>
  void assemble(std::vector<Triplet> &triplets, const size_t startI = 0,
                const size_t startJ = 0) const {

    const RowElements &a = this->m_row_elements;
    const ColElements &b = this->m_col_elements;

    const size_t na = a.size();
    const size_t nb = b.size();

    for (size_t i = 0; i < na; ++i) {
      const_row_reference ai = a[i];
      for (size_t j = 0; j < nb; ++j) {
        const_col_reference bj = b[j];
        const Block element = static_cast<Block>(this->m_function(ai, bj));
        for (size_t ii = 0; ii < BlockRows; ++ii) {
          for (size_t jj = 0; jj < BlockCols; ++jj) {
            triplets.push_back(Triplet(i * BlockRows + ii + startI,
                                       j * BlockCols + jj + startJ,
                                       element(ii, jj)));
          }
        }
      }
    }
  }

  /// Evaluates a matrix-free linear operator given by \p expr \p if_expr,
  /// and particle sets \p a and \p b on a vector rhs and
  /// accumulates the result in vector lhs
  template <typename DerivedLHS, typename DerivedRHS>
  void evaluate(Eigen::DenseBase<DerivedLHS> &lhs,
                const Eigen::DenseBase<DerivedRHS> &rhs) const {

    const RowElements &a = this->m_row_elements;
    const ColElements &b = this->m_col_elements;

    const size_t na = a.size();
    const size_t nb = b.size();

    CHECK(static_cast<size_t>(lhs.size()) == this->rows(),
          "lhs size is inconsistent");
    CHECK(static_cast<size_t>(rhs.size()) == this->cols(),
          "rhs size is inconsistent");

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < na; ++i) {
      const_row_reference ai = a[i];
      for (size_t j = 0; j < nb; ++j) {
        const_col_reference bj = b[j];
        lhs.template segment<BlockRows>(i * BlockRows) +=
            this->m_function(ai, bj) *
            rhs.template segment<BlockCols>(j * BlockCols);
      }
    }
  }

  /// Evaluates a matrix-free linear operator given by \p expr \p if_expr,
  /// and particle sets \p a and \p b on a vector rhs and
  /// accumulates the result in vector lhs
  template <typename LHSType, typename RHSType>
  void evaluate(std::vector<LHSType> &lhs,
                const std::vector<RHSType> &rhs) const {

    const RowElements &a = this->m_row_elements;
    const ColElements &b = this->m_col_elements;

    const size_t na = a.size();
    const size_t nb = b.size();

    CHECK(lhs.size() == na, "lhs size is inconsistent");
    CHECK(rhs.size() == nb, "rhs size is inconsistent");

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < na; ++i) {
      const_row_reference ai = a[i];
      for (size_t j = 0; j < nb; ++j) {
        const_col_reference bj = b[j];
        lhs[i] += this->m_function(ai, bj) * rhs[j];
      }
    }
  }
};

template <typename RowElements, typename ColElements, typename F>
class KernelMatrix : public KernelBase<RowElements, ColElements, F> {
protected:
  typedef KernelBase<RowElements, ColElements, F> base_type;
  static const unsigned int dimension = base_type::dimension;
  typedef typename base_type::double_d double_d;
  typedef typename base_type::int_d int_d;
  typedef typename base_type::const_row_reference const_row_reference;
  typedef typename base_type::const_col_reference const_col_reference;

  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_type;

  matrix_type m_matrix;

public:
  typedef typename base_type::Scalar Scalar;
  typedef typename base_type::Block Block;
  static const size_t BlockRows = base_type::BlockRows;
  static const size_t BlockCols = base_type::BlockCols;

  KernelMatrix(const RowElements &row_elements, const ColElements &col_elements,
               const F &function)
      : base_type(row_elements, col_elements, function) {
    assemble_matrix();
  };

  void assemble_matrix() {
    const RowElements &a = this->m_row_elements;
    const ColElements &b = this->m_col_elements;

    m_matrix.resize(this->rows(), this->cols());
    for (size_t i = 0; i < a.size(); ++i) {
      const_row_reference ai = a[i];
      for (size_t j = 0; j < b.size(); ++j) {
        const_col_reference bj = b[j];
        m_matrix.template block<BlockRows, BlockCols>(i * BlockRows,
                                                      j * BlockCols) =
            static_cast<Block>(this->m_function(ai, bj));
      }
    }
  }

  Scalar coeff(const size_t i, const size_t j) const { return m_matrix(i, j); }

  template <typename MatrixType> void assemble(const MatrixType &matrix) const {
    const_cast<MatrixType &>(matrix) = m_matrix;
  }

  /// Evaluates a matrix-free linear operator given by \p expr \p if_expr,
  /// and particle sets \p a and \p b on a vector rhs and
  /// accumulates the result in vector lhs
  template <typename LHSType, typename RHSType>
  void evaluate(std::vector<LHSType> &lhs,
                const std::vector<RHSType> &rhs) const {

    const RowElements &a = this->m_row_elements;
    const ColElements &b = this->m_col_elements;

    const size_t na = a.size();
    const size_t nb = b.size();

    ASSERT(lhs.size() == na, "lhs size is inconsistent");
    ASSERT(rhs.size() == nb, "rhs size is inconsistent");

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < na; ++i) {
      for (size_t j = 0; j < nb; ++j) {
        lhs[i] += m_matrix.template block<BlockRows, BlockCols>(i * BlockRows,
                                                                j * BlockCols) *
                  rhs[j];
      }
    }
  }

  /// Evaluates a matrix-free linear operator given by \p expr \p if_expr,
  /// and particle sets \p a and \p b on a vector rhs and
  /// accumulates the result in vector lhs
  template <typename DerivedLHS, typename DerivedRHS>
  void evaluate(Eigen::DenseBase<DerivedLHS> &lhs,
                const Eigen::DenseBase<DerivedRHS> &rhs) const {
    ASSERT(lhs.size() == this->rows(), "lhs size not consistent")
    ASSERT(rhs.size() == this->cols(), "lhs size not consistent")
    lhs += m_matrix * rhs;
  }
};

template <typename RowElements, typename ColElements, typename PositionF,
          size_t QuadratureOrder = 8,
          typename F =
              detail::position_kernel<RowElements, ColElements, PositionF>>
class KernelChebyshev : public KernelDense<RowElements, ColElements, F> {
  typedef KernelDense<RowElements, ColElements, F> base_type;
  typedef typename base_type::double_d double_d;
  typedef typename base_type::int_d int_d;
  typedef typename base_type::const_position_reference const_position_reference;
  typedef typename base_type::const_row_reference const_row_reference;
  typedef typename base_type::const_col_reference const_col_reference;

  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_type;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_type;
  typedef Eigen::Map<vector_type> map_type;

  static const unsigned int dimension = base_type::dimension;
  matrix_type m_col_Rn_matrix, m_row_Rn_matrix, m_kernel_matrix;
  unsigned int m_order;
  unsigned int m_ncheb;
  const int_d m_start;
  const int_d m_end;
  PositionF m_position_function;
  mutable vector_type m_W;
  mutable vector_type m_fcheb;

public:
  typedef typename base_type::Scalar Scalar;
  typedef typename base_type::Block Block;
  static const size_t BlockRows = base_type::BlockRows;
  static const size_t BlockCols = base_type::BlockCols;

  KernelChebyshev(const RowElements &row_elements,
                  const ColElements &col_elements, const unsigned int n,
                  const PositionF &function)
      : base_type(row_elements, col_elements, F(function)), m_order(n),
        m_ncheb(std::pow(n, dimension)), m_start(int_d::Constant(0)),
        m_end(int_d::Constant(n)), m_position_function(function) {
    set_n(n);
  };

  void set_n(const unsigned int n) {
    m_order = n;
    m_ncheb = std::pow(n, dimension);
    m_W.resize(m_ncheb * BlockCols);
    m_fcheb.resize(m_ncheb * BlockRows);

    update_row_positions();
    update_col_positions();
    update_kernel_matrix();
  }

  void update_row_positions() {
    bbox<dimension> row_box(this->m_row_elements.get_min(),
                            this->m_row_elements.get_max());
    detail::integrate_chebyshev<RowElements, BlockRows, QuadratureOrder>
        integrate(this->m_row_elements, m_order, row_box);

    // fill row_Rn matrix
    m_row_Rn_matrix.resize(this->m_row_elements.size() * BlockRows,
                           m_ncheb * BlockRows);
    integrate(m_row_Rn_matrix);
  }

  void update_kernel_matrix() {
    bbox<dimension> row_box(this->m_row_elements.get_min(),
                            this->m_row_elements.get_max());
    bbox<dimension> col_box(this->m_col_elements.get_min(),
                            this->m_col_elements.get_max());
    detail::ChebyshevRn<dimension> col_Rn(m_order, row_box);
    detail::ChebyshevRn<dimension> row_Rn(m_order, col_box);

    // fill kernel matrix
    m_kernel_matrix.resize(m_ncheb * BlockRows, m_ncheb * BlockCols);
    lattice_iterator<dimension> mi(m_start, m_end);
    for (size_t i = 0; i < m_ncheb; ++i, ++mi) {
      const double_d pi = col_Rn.get_position(*mi);
      lattice_iterator<dimension> mj(m_start, m_end);
      for (size_t j = 0; j < m_ncheb; ++j, ++mj) {
        const double_d pj = row_Rn.get_position(*mj);
        m_kernel_matrix.template block<BlockRows, BlockCols>(i * BlockRows,
                                                             j * BlockCols) =
            static_cast<Block>(m_position_function(pi, pj));
      }
    }
  }

  void update_col_positions() {
    bbox<dimension> col_box(this->m_col_elements.get_min(),
                            this->m_col_elements.get_max());
    detail::integrate_chebyshev<ColElements, BlockCols, QuadratureOrder>
        integrate(this->m_col_elements, m_order, col_box);

    // fill row_Rn matrix
    m_col_Rn_matrix.resize(this->m_col_elements.size() * BlockCols,
                           m_ncheb * BlockCols);
    integrate(m_col_Rn_matrix);
  }

  /// Evaluates a matrix-free linear operator given by \p expr \p if_expr,
  /// and particle sets \p a and \p b on a vector rhs and
  /// accumulates the result in vector lhs
  template <typename DerivedLHS, typename DerivedRHS>
  void evaluate(Eigen::DenseBase<DerivedLHS> &lhs,
                const Eigen::DenseBase<DerivedRHS> &rhs) const {

    const ColElements &b = this->m_col_elements;

    CHECK(!b.get_periodic().any(), "chebyshev operator assumes not periodic");
    ASSERT(static_cast<typename DerivedLHS::Index>(this->rows()) == lhs.rows(),
           "lhs vector has incompatible size");
    ASSERT(static_cast<typename DerivedRHS::Index>(this->cols()) == rhs.rows(),
           "rhs vector has incompatible size");

    // First compute the weights at the Chebyshev nodes ym
    // by anterpolation
    m_W = m_col_Rn_matrix.transpose() * rhs.derived();

    // Next compute f ðxÞ at the Chebyshev nodes xl:
    m_fcheb = m_kernel_matrix * m_W;

    // Last compute f ðxÞ at the observation points xi by interpolation:
    lhs = m_row_Rn_matrix * m_fcheb;
  }
};

#ifdef HAVE_H2LIB

template <typename RowElements, typename ColElements, typename PositionF,
          typename F>
class KernelH2 : public KernelDense<RowElements, ColElements, F> {
  typedef KernelDense<RowElements, ColElements, F> base_type;
  typedef typename ColElements::query_type query_type;
  static const unsigned int dimension = base_type::dimension;
  typedef typename detail::H2LibBlackBoxExpansions<
      dimension, PositionF, base_type::BlockRows, base_type::BlockCols>
      expansions_type;
  typedef H2LibMatrix h2_matrix_type;

  PositionF m_position_function;
  h2_matrix_type m_h2_matrix;

public:
  typedef typename base_type::Block Block;
  typedef typename base_type::Scalar Scalar;
  static const size_t BlockRows = base_type::BlockRows;
  static const size_t BlockCols = base_type::BlockCols;
  static_assert(detail::is_particles<RowElements>::value &&
                    detail::is_particles<ColElements>::value,
                "only implemented for particle elements");

  KernelH2(const RowElements &row_elements, const ColElements &col_elements,
           const int order, const PositionF &position_function,
           const F &function, const double eta)
      : base_type(row_elements, col_elements, function),
        m_position_function(position_function),
        m_h2_matrix(row_elements, col_elements,
                    expansions_type(order, position_function), function, eta) {}

  const h2_matrix_type &get_h2_matrix() const { return m_h2_matrix; }

  const PositionF &get_position_function() const { return m_position_function; }

  void compress(const double tol) { m_h2_matrix.compress(tol); }

  /// Evaluates a h2 matrix linear operator given by \p expr \p if_expr,
  /// and particle sets \p a and \p b on a vector rhs and
  /// accumulates the result in vector lhs
  template <typename VectorLHS, typename VectorRHS>
  void evaluate(VectorLHS &lhs, const VectorRHS &rhs) const {
    m_h2_matrix.matrix_vector_multiply(lhs, 1.0, false, rhs);
  }
};
#endif

template <typename RowElements, typename ColElements, typename PositionF,
          typename F, unsigned int N>
class KernelFMM : public KernelDense<RowElements, ColElements, F> {
  typedef KernelDense<RowElements, ColElements, F> base_type;
  typedef typename base_type::position position;
  typedef typename base_type::double_d double_d;
  typedef typename base_type::int_d int_d;
  typedef typename base_type::const_position_reference const_position_reference;
  typedef typename base_type::const_row_reference const_row_reference;
  typedef typename base_type::const_col_reference const_col_reference;
  typedef typename ColElements::query_type query_type;
  static const unsigned int dimension = base_type::dimension;
  typedef typename detail::BlackBoxExpansions<
      dimension, N, PositionF, base_type::BlockRows, base_type::BlockCols>
      expansions_type;
  typedef FastMultipoleMethod<expansions_type, F, RowElements, ColElements>
      fmm_type;

  expansions_type m_expansions;
  fmm_type m_fmm;

public:
  typedef typename base_type::Block Block;
  typedef typename base_type::Scalar Scalar;
  static const size_t BlockRows = base_type::BlockRows;
  static const size_t BlockCols = base_type::BlockCols;
  static_assert(detail::is_particles<RowElements>::value &&
                    detail::is_particles<ColElements>::value,
                "only implemented for particle elements");

  KernelFMM(const RowElements &row_elements, const ColElements &col_elements,
            const PositionF &position_function, const F &function)
      : base_type(row_elements, col_elements, function),
        m_expansions(position_function),
        m_fmm(row_elements, col_elements, m_expansions, function){};

  /// Evaluates a matrix-free linear operator given by \p expr \p if_expr,
  /// and particle sets \p a and \p b on a vector rhs and
  /// accumulates the result in vector lhs
  template <typename VectorLHS, typename VectorRHS>
  void evaluate(VectorLHS &lhs, const VectorRHS &rhs) const {
    m_fmm.matrix_vector_multiply(lhs, rhs);
  }
};

template <typename RowElements, typename ColElements, typename FRadius,
          typename FWithDx,
          typename F =
              detail::sparse_kernel<RowElements, ColElements, FRadius, FWithDx>>
class KernelSparse : public KernelBase<RowElements, ColElements, F> {
protected:
  typedef KernelBase<RowElements, ColElements, F> base_type;
  typedef typename base_type::position position;
  typedef typename base_type::double_d double_d;
  typedef typename base_type::const_position_reference const_position_reference;
  typedef typename base_type::const_row_reference const_row_reference;
  typedef typename base_type::const_col_reference const_col_reference;
  static_assert(detail::is_particles<RowElements>::value &&
                    detail::is_particles<ColElements>::value,
                "only implemented for particle elements");

  FRadius m_radius_function;
  FWithDx m_dx_function;

public:
  typedef typename base_type::Block Block;
  typedef typename base_type::Scalar Scalar;
  static const size_t BlockRows = base_type::BlockRows;
  static const size_t BlockCols = base_type::BlockCols;

  KernelSparse(const RowElements &row_elements, const ColElements &col_elements,
               const FRadius &radius_function, const FWithDx &withdx_function)
      : base_type(row_elements, col_elements,
                  F(col_elements, radius_function, withdx_function)),
        m_radius_function(radius_function), m_dx_function(withdx_function){};

  /*
   * shouldn't need this anymore....
  Block coeff(const size_t i, const size_t j) const {
      const int pi = std::floor(static_cast<float>(i)/BlockRows);
      const int ioffset = i - pi*BlockRows;
      const int pj = std::floor(static_cast<float>(j)/BlockCols);
      const int joffset = j - pj*BlockCols;
      ASSERT(pi < this->m_row_elements.size(),"pi greater than a.size()");
      ASSERT(pj < this->m_col_elements.size(),"pj greater than b.size()");
      const_row_reference ai = this->m_row_elements[pi];
      const_col_reference bj = this->m_col_elements[pj];
      const_position_reference dx = get<position>(bj)-get<position>(ai);
      if (dx.squaredNorm() < std::pow(m_radius_function(ai),2)) {
          return this->m_function(dx,ai,bj)(ioffset,joffset);
      } else {
          return 0.0;
      }
  }
  */

  template <typename MatrixType> void assemble(const MatrixType &matrix) const {

    const RowElements &a = this->m_row_elements;
    const ColElements &b = this->m_col_elements;

    const size_t na = a.size();

    const_cast<MatrixType &>(matrix).setZero();

    // sparse a x b block
    for (size_t i = 0; i < na; ++i) {
      const_row_reference ai = a[i];
      const double radius = m_radius_function(ai);
      for (auto pairj =
               euclidean_search(b.get_query(), get<position>(ai), radius);
           pairj != false; ++pairj) {
        const_col_reference bj = *pairj;
        const_position_reference dx = pairj.dx();
        const size_t j = &get<position>(bj) - get<position>(b).data();
        const_cast<MatrixType &>(matrix).template block<BlockRows, BlockCols>(
            i * BlockRows, j * BlockCols) =
            static_cast<Block>(m_dx_function(dx, ai, bj));
      }
    }
  }

  template <typename Triplet>
  void assemble(std::vector<Triplet> &triplets, const size_t startI = 0,
                const size_t startJ = 0) const {

    const RowElements &a = this->m_row_elements;
    const ColElements &b = this->m_col_elements;

    const size_t na = a.size();

    // sparse a x b block
    // std::cout << "sparse a x b block" << std::endl;
    for (size_t i = 0; i < na; ++i) {
      const_row_reference ai = a[i];
      const double radius = m_radius_function(ai);
      for (auto pairj =
               euclidean_search(b.get_query(), get<position>(ai), radius);
           pairj != false; ++pairj) {
        const_col_reference bj = *pairj;
        const_position_reference dx = pairj.dx();
        const size_t j = &get<position>(bj) - get<position>(b).data();
        const Block element = static_cast<Block>(m_dx_function(dx, ai, bj));
        for (size_t ii = 0; ii < BlockRows; ++ii) {
          for (size_t jj = 0; jj < BlockCols; ++jj) {
            triplets.push_back(Triplet(i * BlockRows + ii + startI,
                                       j * BlockCols + jj + startJ,
                                       element(ii, jj)));
          }
        }
      }
    }
  }

  /// Evaluates a matrix-free linear operator given by \p expr \p if_expr,
  /// and particle sets \p a and \p b on a vector rhs and
  /// accumulates the result in vector lhs
  ///
  template <typename LHSType, typename RHSType>
  void evaluate(std::vector<LHSType> &lhs,
                const std::vector<RHSType> &rhs) const {

    const RowElements &a = this->m_row_elements;
    const ColElements &b = this->m_col_elements;

    const size_t na = a.size();

    ASSERT(na == rhs.size(), "lhs vector has incompatible size");
    ASSERT(b.size() == lhs.size(), "rhs vector has incompatible size");

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < na; ++i) {
      const_row_reference ai = a[i];
      const double radius = m_radius_function(ai);
      for (auto pairj =
               euclidean_search(b.get_query(), get<position>(ai), radius);
           pairj != false; ++pairj) {
        const_position_reference dx = pairj.dx();
        const_col_reference bj = *pairj;
        const size_t j = &get<position>(bj) - get<position>(b).data();
        lhs[i] += m_dx_function(dx, ai, bj) * rhs[j];
      }
    }
  }

  template <typename DerivedLHS, typename DerivedRHS>
  void evaluate(Eigen::DenseBase<DerivedLHS> &lhs,
                const Eigen::DenseBase<DerivedRHS> &rhs) const {

    ASSERT(static_cast<typename DerivedLHS::Index>(this->rows()) == lhs.rows(),
           "lhs vector has incompatible size");
    ASSERT(static_cast<typename DerivedRHS::Index>(this->cols()) == rhs.rows(),
           "rhs vector has incompatible size");

    const RowElements &a = this->m_row_elements;
    const ColElements &b = this->m_col_elements;

    const size_t na = a.size();

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < na; ++i) {
      const_row_reference ai = a[i];
      const double radius = m_radius_function(ai);
      for (auto pairj =
               euclidean_search(b.get_query(), get<position>(ai), radius);
           pairj != false; ++pairj) {
        const_position_reference dx = pairj.dx();
        const_col_reference bj = *pairj;
        const size_t j = &get<position>(bj) - get<position>(b).data();
        lhs.template segment<BlockRows>(i * BlockRows) +=
            m_dx_function(dx, ai, bj) *
            rhs.template segment<BlockCols>(j * BlockCols);
      }
    }
  }
};

template <typename RowElements, typename ColElements, typename F,
          typename RadiusFunction = detail::constant_radius<RowElements>>
class KernelSparseConst
    : public KernelSparse<RowElements, ColElements, RadiusFunction, F> {
  typedef KernelSparse<RowElements, ColElements, RadiusFunction, F> base_type;
  typedef typename base_type::position position;
  typedef typename base_type::double_d double_d;
  typedef typename base_type::const_position_reference const_position_reference;
  typedef typename base_type::const_row_reference const_row_reference;
  typedef typename base_type::const_col_reference const_col_reference;

public:
  typedef typename base_type::Block Block;
  typedef typename base_type::Scalar Scalar;
  static const size_t BlockRows = base_type::BlockRows;
  static const size_t BlockCols = base_type::BlockCols;

  KernelSparseConst(const RowElements &row_elements,
                    const ColElements &col_elements, const double radius,
                    const F &function)
      : base_type(row_elements, col_elements, RadiusFunction(radius),
                  function) {}
};

template <typename RowElements, typename ColElements,
          typename F = detail::zero_kernel<RowElements, ColElements>>
class KernelZero : public KernelBase<RowElements, ColElements, F> {

  typedef KernelBase<RowElements, ColElements, F> base_type;
  typedef typename base_type::position position;
  typedef typename base_type::double_d double_d;
  typedef typename base_type::const_position_reference const_position_reference;
  typedef typename base_type::const_row_reference const_row_reference;
  typedef typename base_type::const_col_reference const_col_reference;

public:
  typedef typename base_type::Block Block;
  typedef typename base_type::Scalar Scalar;

  KernelZero(const RowElements &row_elements, const ColElements &col_elements)
      : base_type(row_elements, col_elements, F()){};

  template <typename MatrixType> void assemble(const MatrixType &matrix) const {
    const_cast<MatrixType &>(matrix).setZero();
  }

  template <typename Triplet>
  void assemble(std::vector<Triplet> &triplets, const size_t startI = 0,
                const size_t startJ = 0) const {}

  /// Evaluates a matrix-free linear operator given by \p expr \p if_expr,
  /// and particle sets \p a and \p b on a vector rhs and
  /// accumulates the result in vector lhs
  template <typename VectorLHS, typename VectorRHS>
  void evaluate(VectorLHS &lhs, const VectorRHS &rhs) const {}
};

} // namespace Aboria

#endif // HAVE_EIGEN

#endif
