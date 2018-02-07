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

#ifndef OPERATORS_DETAIL_H_
#define OPERATORS_DETAIL_H_

#include "Vector.h"
#include <limits>
#include <type_traits>

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>

namespace Eigen {
namespace internal {
// MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
template <unsigned int NI, unsigned int NJ, typename Blocks>
struct traits<Aboria::MatrixReplacement<NI, NJ, Blocks>>
    : public Eigen::internal::traits<Eigen::SparseMatrix<double>> {};
} // namespace internal

template <typename T, unsigned int N> struct NumTraits<Aboria::Vector<T, N>> {
  typedef Aboria::Vector<T, N> Scalar;
  typedef Aboria::Vector<T, N> Real;
  typedef Aboria::Vector<T, N> NonInteger;
  typedef Aboria::Vector<T, N> Literal;
  typedef Aboria::Vector<T, N> Nested;
  enum {
    IsComplex = 0,
    IsInteger = 0,
    IsSigned = 1,
    RequireInitialization = 0,
    ReadCost = N,
    AddCost = N,
    MulCost = N
  };

  inline static Real epsilon() {
    return Real(std::numeric_limits<T>::epsilon());
  }
  inline static Real dummy_precision() {
    return Real(std::numeric_limits<T>::epsilon());
  }
  inline static Scalar highest() {
    return Scalar(std::numeric_limits<T>::max());
  }
  inline static Scalar lowest() {
    return Scalar(std::numeric_limits<T>::lowest());
  }
  inline static int digits10() { return N * std::numeric_limits<T>::digits10; }
};
} // namespace Eigen

namespace Aboria {
namespace detail {

template <typename MatrixReplacement,
          typename Scalar = typename MatrixReplacement::Scalar>
class InnerIterator {
public:
  typedef const Scalar *pointer;
  typedef std::forward_iterator_tag iterator_category;
  typedef const Scalar &reference;
  typedef const Scalar value_type;
  typedef std::ptrdiff_t difference_type;
  typedef size_t Index;

  InnerIterator(const MatrixReplacement &mat, const Index row)
      : m_mat(mat), m_row(row), m_col(0){};

  InnerIterator(const InnerIterator &other)
      : m_mat(other.m_mat), m_row(other.m_row), m_col(other.m_col){};

  Scalar value() const { return dereference(); }

  Index index() const { return m_col; }

  operator bool() const {
    return (m_row < m_mat.rows()) && (m_col < m_mat.cols());
  }

  bool equal(InnerIterator const &other) const {
    return (m_row == other.m_row) && (m_col == other.m_col);
  }

  Scalar dereference() const { return m_mat.coeff(m_row, m_col); }

  void increment() {
    m_col++;
    ASSERT(m_col < m_mat.cols(), "InnerIterator outside cols range");
  }

  Scalar operator*() { return dereference(); }
  Scalar operator->() { return dereference(); }
  InnerIterator &operator++() {
    increment();
    return *this;
  }
  InnerIterator operator++(int) {
    InnerIterator tmp(*this);
    operator++();
    return tmp;
  }

  size_t operator-(InnerIterator start) const {
    ASSERT(m_row == start.m_row,
           "Difference between InnerIterators must have identical row numbers");
    return (m_col - start.m_col);
  }

  inline bool operator==(const InnerIterator &rhs) { return equal(rhs); }

  inline bool operator!=(const InnerIterator &rhs) { return !operator==(rhs); }

private:
  friend class boost::iterator_core_access;
  const MatrixReplacement &m_mat;
  const Index m_row;
  Index m_col;
};

template <typename T1 = void> int sum() { return 0; }

template <typename T1, typename... T> T1 sum(T1 s, T... ts) {
  return s + sum(ts...);
}

template <typename Dest, typename Source, typename Block>
void evalTo_block(Eigen::VectorBlock<Dest> y,
                  const Eigen::VectorBlock<Source> &rhs, const Block &block) {
  block.evaluate(y, rhs);
}

template <typename Dest, unsigned int NI, unsigned int NJ, typename Blocks,
          typename Rhs>
void evalTo_unpack_blocks(Dest &y, const MatrixReplacement<NI, NJ, Blocks> &lhs,
                          const Rhs &rhs) {}

template <typename Dest, unsigned int NI, unsigned int NJ, typename Blocks,
          typename Rhs, typename I, typename J, typename T1, typename... T>
void evalTo_unpack_blocks(Dest &y, const MatrixReplacement<NI, NJ, Blocks> &lhs,
                          const Rhs &rhs, const std::tuple<I, J, T1 &> &block,
                          const T &... other_blocks) {
  evalTo_block(y.segment(lhs.template start_row<I::value>(),
                         lhs.template size_row<I::value>()),
               rhs.segment(lhs.template start_col<J::value>(),
                           lhs.template size_col<J::value>()),
               std::get<2>(block));
  evalTo_unpack_blocks(y, lhs, rhs, other_blocks...);
}

template <typename Dest, unsigned int NI, unsigned int NJ, typename Blocks,
          typename Rhs, std::size_t... I>
void evalTo_impl(Dest &y, const MatrixReplacement<NI, NJ, Blocks> &lhs,
                 const Rhs &rhs, detail::index_sequence<I...>) {
  evalTo_unpack_blocks(
      y, lhs, rhs,
      std::tuple<mpl::int_<I / NJ>, mpl::int_<I % NJ>,
                 typename std::tuple_element<I, Blocks>::type const &>(
          mpl::int_<I / NJ>(), mpl::int_<I % NJ>(),
          std::get<I>(lhs.m_blocks))...);
}
} // namespace detail
} // namespace Aboria

// Implementation of MatrixReplacement * Eigen::DenseVector though a
// specialization of internal::generic_product_impl:
namespace Eigen {
namespace internal {
template <typename Rhs, unsigned int NI, unsigned int NJ, typename Blocks>
struct generic_product_impl<Aboria::MatrixReplacement<NI, NJ, Blocks>, Rhs,
                            SparseShape, DenseShape,
                            GemvProduct> // GEMV stands for matrix-vector
    : generic_product_impl_base<
          Aboria::MatrixReplacement<NI, NJ, Blocks>, Rhs,
          generic_product_impl<Aboria::MatrixReplacement<NI, NJ, Blocks>,
                               Rhs>> {

  typedef
      typename Product<Aboria::MatrixReplacement<NI, NJ, Blocks>, Rhs>::Scalar
          Scalar;
  template <typename Dest>
  CUDA_HOST_DEVICE static void
  scaleAndAddTo(Dest &y, const Aboria::MatrixReplacement<NI, NJ, Blocks> &lhs,
                const Rhs &rhs, const Scalar &alpha) {
    // This method should implement "y += alpha * lhs * rhs" inplace,
    // however, for iterative solvers, alpha is always equal to 1, so let's not
    // bother about it.
#ifdef __CUDA_ARCH__
    ERROR_CUDA("MatrixReplacement class unusable from device code");
#else
    assert(alpha == Scalar(1) && "scaling is not implemented");
    evalTo_impl(y, lhs, rhs, Aboria::detail::make_index_sequence<NI * NJ>());
#endif
  }
};
} // namespace internal
} // namespace Eigen

#endif // OPERATORS_H_
