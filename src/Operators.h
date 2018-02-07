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

#ifndef OPERATORS_H_
#define OPERATORS_H_

#include <type_traits>

namespace Aboria {
template <unsigned int NI, unsigned int NJ, typename Blocks>
class MatrixReplacement;
}

#ifdef HAVE_EIGEN
#include "detail/Operators.h"

namespace Aboria {

/// \brief A matrix-replacement class for use with Eigen
///
/// This provides a class that acts like a sparse matrix within
/// Eigen, at least for matrix-vector multiplication and the
/// iterative solvers (it has not been tested with anything else
/// and will not work, for example, in matrix-matrix addition or
/// multiplication).
///
/// It is templated on a set of \p NI x \p NJ blocks, each containing a
/// kernel defined in Kernels.h . This kernel performs the actual
/// operator
///
///  \tparam NI the number of kernel block rows in the operator
///  \tparam NJ the number of kernel block columns in the operator
///  \tparam Blocks a 3-tuple of (`particle_type_row`,
///  `particle_type_col`,`kernel_type`),
///         where `particle_type_row` is the type of the row particle set,
///         `particle_type_col` is the type of the column particle set,
///         and `kernel_type` is the type of kernel used by the block
///
///  \see Aboria::create_dense_operator()
///       create_sparse_operator()
///       create_zero_operator()
///       create_chebyshev_operator()
template <unsigned int NI, unsigned int NJ, typename Blocks>
class MatrixReplacement
    : public Eigen::EigenBase<MatrixReplacement<NI, NJ, Blocks>> {

  typedef typename std::tuple_element<0, Blocks>::type first_block_type;

public:
  // Expose some compile-time information to Eigen:
  typedef typename first_block_type::Scalar Scalar;
  typedef Scalar RealScalar;
  typedef size_t Index;
  typedef int StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    RowsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
    MaxRowsAtCompileTime = Eigen::Dynamic,
    IsRowMajor = false
  };
  typedef detail::InnerIterator<MatrixReplacement> InnerIterator;

  MatrixReplacement(const Blocks &blocks) : m_blocks(blocks){};
  MatrixReplacement(Blocks &&blocks) : m_blocks(std::move(blocks)){};

  CUDA_HOST_DEVICE
  Index rows() const {
    // std::cout << "rows = " << rows_impl(detail::make_index_sequence<NI>()) <<
    // std::endl;
    //
#ifdef __CUDA_ARCH__
    ERROR_CUDA("MatrixReplacement class unusable from device code");
    return 0;
#else
    return rows_impl(detail::make_index_sequence<NI>());
#endif
  }

  CUDA_HOST_DEVICE
  Index cols() const {
    // std::cout << "cols = " << cols_impl(detail::make_index_sequence<NJ>())<<
    // std::endl;
#ifdef __CUDA_ARCH__
    ERROR_CUDA("MatrixReplacement class unusable from device code");
    return 0;
#else
    return cols_impl(detail::make_index_sequence<NJ>());
#endif
  }

  CUDA_HOST_DEVICE
  Index innerSize() const {
#ifdef __CUDA_ARCH__
    ERROR_CUDA("MatrixReplacement class unusable from device code");
    return 0;
#else
    return rows();
#endif
  }

  CUDA_HOST_DEVICE
  Index outerSize() const {
#ifdef __CUDA_ARCH__
    ERROR_CUDA("MatrixReplacement class unusable from device code");
    return 0;
#else
    return cols();
#endif
  }

  void resize(Index a_rows, Index a_cols) {
    // This method should not be needed in the future.
    assert(a_rows == 0 && a_cols == 0 || a_rows == rows() && a_cols == cols());
  }

  Scalar coeff(const Index i, const Index j) const {
    return coeff_impl(i, j, detail::make_index_sequence<NI * NJ>());
  }

  template <typename Rhs>
  Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct>
  operator*(const Eigen::MatrixBase<Rhs> &x) const {
    return Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct>(
        *this, x.derived());
  }

  template <unsigned int I, unsigned int J>
  const typename std::tuple_element<I * NJ + J, Blocks>::type &
  get_kernel() const {
    return std::get<I * NJ + J>(m_blocks);
  }

  const typename std::tuple_element<0, Blocks>::type &get_first_kernel() const {
    return get_kernel<0, 0>();
  }

  template <unsigned int I, unsigned int J>
  typename std::tuple_element<I * NJ + J, Blocks>::type &get_kernel() {
    return std::get<I * NJ + J>(m_blocks);
  }

  typename std::tuple_element<0, Blocks>::type &get_first_kernel() {
    return get_kernel<0, 0>();
  }

  template <typename Derived>
  void assemble(Eigen::DenseBase<Derived> &matrix) const {
    const size_t na = rows();
    const size_t nb = cols();
    // matrix.resize(na,nb);
    CHECK((static_cast<size_t>(matrix.rows()) == na) &&
              (static_cast<size_t>(matrix.cols()) == nb),
          "matrix size is not compatible with expression.");
    assemble_impl(matrix, detail::make_index_sequence<NI * NJ>());
  }

  template <int _Options, typename _StorageIndex>
  void assemble(Eigen::SparseMatrix<Scalar, _Options, _StorageIndex> &matrix) {
    const size_t na = rows();
    const size_t nb = cols();
    // matrix.resize(na,nb);
    CHECK((static_cast<size_t>(matrix.rows()) == na) &&
              (static_cast<size_t>(matrix.cols()) == nb),
          "matrix size is not compatible with expression.");

    typedef Eigen::Triplet<Scalar> triplet_type;
    std::vector<triplet_type> tripletList;
    // TODO: can we estimate this better?
    tripletList.reserve(na * 5);

    assemble_impl(tripletList, detail::make_index_sequence<NI * NJ>());

    matrix.setFromTriplets(tripletList.begin(), tripletList.end());
  }

  template <std::size_t... I>
  Index rows_impl(detail::index_sequence<I...>) const {
    return detail::sum(std::get<I * NJ>(m_blocks).rows()...);
  }

  template <std::size_t... J>
  Index cols_impl(detail::index_sequence<J...>) const {
    return detail::sum(std::get<J>(m_blocks).cols()...);
  }

  template <int I> Index start_col() const {
    return cols_impl(detail::make_index_sequence<I>());
  }

  template <int I> Index size_col() const {
    return std::get<I>(m_blocks).cols();
  }

  template <int I> Index start_row() const {
    return rows_impl(detail::make_index_sequence<I>());
  }

  template <int I> Index size_row() const {
    return std::get<I * NJ>(m_blocks).rows();
  }

  template <typename block_type>
  Scalar coeff_impl_block(const Index i, const Index j,
                          const block_type &block) const {
    return block.coeff(i, j);
  }

  template <std::size_t... I>
  Scalar coeff_impl(const Index i, const Index j,
                    detail::index_sequence<I...>) const {

    return detail::sum(
        ((i >= start_row<I / NJ>()) && (i < start_row<I / NJ + 1>()) &&
         (j >= start_col<I % NJ>()) && (j < start_col<I % NJ + 1>()))
            ? (coeff_impl_block(i - start_row<I / NJ>(),
                                j - start_col<I % NJ>(), std::get<I>(m_blocks)))
            : (0.0)...);
  }

  template <typename Block>
  void assemble_block_impl(const size_t startI, const size_t startJ,
                           std::vector<Eigen::Triplet<Scalar>> &triplets,
                           const Block &block) const {

    block.assemble(triplets, startI, startJ);
  }

  template <typename Block, typename Derived>
  void assemble_block_impl(const Eigen::MatrixBase<Derived> &matrix,
                           const Block &block) const {
    block.assemble(matrix);
  }

  template <std::size_t... I>
  void assemble_impl(std::vector<Eigen::Triplet<Scalar>> &triplets,
                     detail::index_sequence<I...>) const {
    int dummy[] = {
        0, (assemble_block_impl(start_row<I / NJ>(), start_col<I % NJ>(),
                                triplets, std::get<I>(m_blocks)),
            void(), 0)...};
    static_cast<void>(dummy);
  }

  template <typename Derived, std::size_t... I>
  void assemble_impl(Eigen::DenseBase<Derived> &matrix,
                     detail::index_sequence<I...>) const {
    int dummy[] = {
        0, (assemble_block_impl(
                matrix.block(start_row<I / NJ>(), start_col<I % NJ>(),
                             start_row<I / NJ + 1>() - start_row<I / NJ>(),
                             start_col<I % NJ + 1>() - start_col<I % NJ>()),
                std::get<I>(m_blocks)),
            void(), 0)...};
    static_cast<void>(dummy);
  }

  Blocks m_blocks;
};

/// \brief creates a dense matrix-free linear operator for use with Eigen
///
/// This function returns a MatrixReplacement object that acts like a
/// dense linear operator (i.e. matrix) in Eigen.
///
/// \param row_particles The rows of the linear operator index this
///                      first particle set
/// \param col_particles The columns of the linear operator index this
///                      first particle set
/// \param function A function object that returns the value of the operator
///                 for a given particle pair
///
///
/// \tparam RowParticles The type of the row particle set
/// \tparam ColParticles The type of the column particle set
/// \tparam F The type of the function object
template <typename RowParticles, typename ColParticles, typename F,
          typename Kernel = KernelDense<RowParticles, ColParticles, F>,
          typename Operator = MatrixReplacement<1, 1, std::tuple<Kernel>>>
Operator create_dense_operator(const RowParticles &row_particles,
                               const ColParticles &col_particles,
                               const F &function) {
  return Operator(
      std::make_tuple(Kernel(row_particles, col_particles, function)));
}

/// \brief creates a matrix linear operator for use with Eigen
///
/// This function returns a MatrixReplacement object that acts like a
/// dense linear operator (i.e. matrix) in Eigen.
///
/// \param row_particles The rows of the linear operator index this
///                      first particle set
/// \param col_particles The columns of the linear operator index this
///                      first particle set
/// \param function A function object that returns the value of the operator
///                 for a given particle pair
///
///
/// \tparam RowParticles The type of the row particle set
/// \tparam ColParticles The type of the column particle set
/// \tparam F The type of the function object
template <typename RowParticles, typename ColParticles, typename F,
          typename Kernel = KernelMatrix<RowParticles, ColParticles, F>,
          typename Operator = MatrixReplacement<1, 1, std::tuple<Kernel>>>
Operator create_matrix_operator(const RowParticles &row_particles,
                                const ColParticles &col_particles,
                                const F &function) {
  return Operator(
      std::make_tuple(Kernel(row_particles, col_particles, function)));
}

/// \brief creates a matrix-free linear operator using chebyshev interpolation
///        for use with Eigen
///
/// This function returns a MatrixReplacement object that acts like a
/// dense linear operator (i.e. matrix) in Eigen. It uses chebyshev
/// interpolation to speed up its application to a vector.
///
/// \param row_particles The rows of the linear operator index this
///                      first particle set
/// \param col_particles The columns of the linear operator index this
///                      first particle set
/// \param n The number of chebyshev nodes in each dimension to use
/// \param function A function object that returns the value of the operator
///                 for a given particle pair
///
/// \tparam RowParticles The type of the row particle set
/// \tparam ColParticles The type of the column particle set
/// \tparam F The type of the function object
template <typename RowParticles, typename ColParticles, typename F,
          typename Kernel = KernelChebyshev<RowParticles, ColParticles, F>,
          typename Operator = MatrixReplacement<1, 1, std::tuple<Kernel>>>
Operator create_chebyshev_operator(const RowParticles &row_particles,
                                   const ColParticles &col_particles,
                                   const unsigned int n, const F &function) {
  return Operator(
      std::make_tuple(Kernel(row_particles, col_particles, n, function)));
}

/// \brief creates a matrix-free linear operator using the Black-Box
///        Fast Multipole Method (FMM)
///
/// This function returns a MatrixReplacement object that acts like a
/// dense linear operator (i.e. matrix) in Eigen. It uses FMM
/// speed up its application to a vector. For repeated uses it is better
/// to use create_h2_operator
///
/// \param row_particles The rows of the linear operator index this
///                      first particle set
/// \param col_particles The columns of the linear operator index this
///                      first particle set
/// \param position_function A function object that returns the value of the
/// operator
///                 for a given position pair (used for well-separated position
///                 pairs)
/// \param function A function object that returns the value of the operator
///                 for a given particle pair (used for close particle pairs)
///
/// \tparam N The number of chebyshev nodes in each dimension to use
/// \tparam RowParticles The type of the row particle set
/// \tparam ColParticles The type of the column particle set
/// \tparam F The type of the function object
template <
    unsigned int N, typename RowParticles, typename ColParticles,
    typename PositionF, typename F,
    typename Kernel = KernelFMM<RowParticles, ColParticles, PositionF, F, N>,
    typename Operator = MatrixReplacement<1, 1, std::tuple<Kernel>>>
Operator create_fmm_operator(const RowParticles &row_particles,
                             const ColParticles &col_particles,
                             const PositionF &position_function,
                             const F &function) {
  return Operator(std::make_tuple(
      Kernel(row_particles, col_particles, position_function, function)));
}

#ifdef HAVE_H2LIB

/// \brief creates a matrix-free linear operator using the Black-Box
///        Fast Multipole Method (FMM) to internally create a H2
///        hierarchical matrix
///
/// This function returns a MatrixReplacement object that acts like a
/// dense linear operator (i.e. matrix) in Eigen. It internally creates
/// a H2 hierarchical matrix to speed up repeated applications of the
/// operator
///
/// \param row_particles The rows of the linear operator index this
///                      first particle set
/// \param col_particles The columns of the linear operator index this
///                      first particle set
/// \param order Chebyshev order of the approximation within each cluster
/// \param position_function A function object that returns the value of the
/// operator
///                 for a given position pair (used for well-separated position
///                 pairs)
/// \param function A function object that returns the value of the operator
///                 for a given particle pair (used for close particle pairs)
///
/// \tparam N The number of chebyshev nodes in each dimension to use
/// \tparam RowParticles The type of the row particle set
/// \tparam ColParticles The type of the column particle set
/// \tparam F The type of the function object
template <typename RowParticles, typename ColParticles, typename PositionF,
          typename F,
          typename Kernel = KernelH2<RowParticles, ColParticles, PositionF, F>,
          typename Operator = MatrixReplacement<1, 1, std::tuple<Kernel>>>
Operator create_h2_operator(const RowParticles &row_particles,
                            const ColParticles &col_particles, const int order,
                            const PositionF &position_function,
                            const F &function, const double eta = 1.0) {
  return Operator(std::make_tuple(Kernel(row_particles, col_particles, order,
                                         position_function, function, eta)));
}

#endif

/// \brief creates a sparse matrix-free linear operator for use with Eigen
///
/// This function returns a MatrixReplacement object that acts like a
/// sparse linear operator (i.e. matrix) in Eigen, in that only particle
/// pairs (i.e. a row/column pair) with a separation less that a given value
/// are considered to be non-zero
///
/// \param row_particles The rows of the linear operator index this
///                      first particle set
/// \param col_particles The columns of the linear operator index this
///                      first particle set
/// \param radius_function A function object that takes a const_reference to
///                 a particle and returns a double value. It is assumed that \p
///                 function returns a zero value for all particle pairs with a
///                 separation greater than this value
/// \param function A function object that returns the value of the operator
///                 for a given particle pair
///
///
/// \tparam RowParticles The type of the row particle set
/// \tparam ColParticles The type of the column particle set
/// \tparam F The type of the function object
template <
    typename RowParticles, typename ColParticles, typename FRadius, typename F,
    typename Kernel = KernelSparse<RowParticles, ColParticles, FRadius, F>,
    typename Operator = MatrixReplacement<1, 1, std::tuple<Kernel>>,
    typename =
        typename std::enable_if<!std::is_arithmetic<FRadius>::value>::type>
Operator create_sparse_operator(const RowParticles &row_particles,
                                const ColParticles &col_particles,
                                const FRadius &radius_function,
                                const F &function) {
  return Operator(std::make_tuple(
      Kernel(row_particles, col_particles, radius_function, function)));
}

/// \brief creates a sparse matrix-free linear operator for use with Eigen
///
/// This function returns a MatrixReplacement object that acts like a
/// sparse linear operator (i.e. matrix) in Eigen, in that only particle
/// pairs (i.e. a row/column pair) with a separation less that a given value
/// are considered to be non-zero
///
/// \param row_particles The rows of the linear operator index this
///                      first particle set
/// \param col_particles The columns of the linear operator index this
///                      first particle set
/// \param radius   It is assumed that \p function
///                 returns a zero value
///                 for all particle pairs with a separation greater than
///                 this value
/// \param function A function object that returns the value of the operator
///                 for a given particle pair
///
///
/// \tparam RowParticles The type of the row particle set
/// \tparam ColParticles The type of the column particle set
/// \tparam F The type of the function object
template <typename RowParticles, typename ColParticles, typename F,
          typename Kernel = KernelSparseConst<RowParticles, ColParticles, F>,
          typename Operator = MatrixReplacement<1, 1, std::tuple<Kernel>>>
Operator create_sparse_operator(const RowParticles &row_particles,
                                const ColParticles &col_particles,
                                const double radius, const F &function) {
  return Operator(
      std::make_tuple(Kernel(row_particles, col_particles, radius, function)));
}

/// \brief creates a zero matrix-free linear operator for use with Eigen
///
/// This function returns a MatrixReplacement object that acts like a
/// $n \times m$ zero matrix. That is the application of this linear
/// operator to a vector will always return a zero vector.
///
/// \param row_particles The rows of the linear operator index this
///                      first particle set
/// \param col_particles The columns of the linear operator index this
///                      first particle set
///
/// \tparam RowParticles The type of the row particle set
/// \tparam ColParticles The type of the column particle set
template <typename RowParticles, typename ColParticles,
          typename Kernel = KernelZero<RowParticles, ColParticles>,
          typename Operator = MatrixReplacement<1, 1, std::tuple<Kernel>>>
Operator create_zero_operator(const RowParticles &row_particles,
                              const ColParticles &col_particles) {
  return Operator(std::make_tuple(Kernel(row_particles, col_particles)));
}

/// creates a matrix-free linear block operator for use with Eigen
///
template <unsigned int NI, unsigned int NJ, typename... T>
MatrixReplacement<NI, NJ,
                  std::tuple<typename std::tuple_element<0, T>::type...>>
create_block_operator(const MatrixReplacement<1, 1, T> &... operators) {
  typedef std::tuple<typename std::tuple_element<0, T>::type...> tuple_type;
  return MatrixReplacement<NI, NJ, tuple_type>(
      tuple_type(std::get<0>(operators.m_blocks)...));
}

} // namespace Aboria
#endif // HAVE_EIGEN

#endif // OPERATORS_H_
