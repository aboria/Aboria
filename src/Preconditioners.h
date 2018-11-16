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

#ifndef PRECONDITIONERS_H_
#define PRECONDITIONERS_H_

#include <algorithm>
#include <chrono>
#include <fstream>
#include <unordered_map>

#ifdef HAVE_CAIRO
#include <cairo-svg.h>
#endif

#ifdef HAVE_EIGEN
#include <unsupported/Eigen/SparseExtra>

#include "Operators.h"

namespace Aboria {
typedef std::chrono::system_clock Clock;

namespace detail {
template <typename Function, typename Dest, unsigned int NI, unsigned int NJ,
          typename Blocks, typename Rhs>
void apply_function_to_diagonal_blocks(
    Function &&f, Dest &y, const MatrixReplacement<NI, NJ, Blocks> &mat,
    const Rhs &rhs, std::integral_constant<unsigned int, NI>) {}

template <typename Function, unsigned int NI, unsigned int NJ, typename Blocks>
void apply_function_to_diagonal_blocks(
    Function &&f, const MatrixReplacement<NI, NJ, Blocks> &mat,
    std::integral_constant<unsigned int, NI>) {}

template <typename Function, typename Dest, unsigned int NI, unsigned int NJ,
          typename Blocks, typename Rhs, unsigned int I>
void apply_function_to_diagonal_blocks(
    Function &&f, Dest &x, const MatrixReplacement<NI, NJ, Blocks> &mat,
    const Rhs &b, std::integral_constant<unsigned int, I>) {
  f(x.segment(mat.template start_row<I>(), mat.template size_row<I>()),
    b.segment(mat.template start_col<I>(), mat.template size_col<I>()),
    std::get<I * NJ + I>(mat.m_blocks));
  apply_function_to_diagonal_blocks(
      std::forward<Function>(f), x, mat, b,
      std::integral_constant<unsigned int, I + 1>());
}

template <typename Function, unsigned int NI, unsigned int NJ, typename Blocks,
          unsigned int I>
void apply_function_to_diagonal_blocks(
    Function &&f, const MatrixReplacement<NI, NJ, Blocks> &mat,
    std::integral_constant<unsigned int, I>) {
  f(std::get<I * NJ + I>(mat.m_blocks));
  apply_function_to_diagonal_blocks(
      std::forward<Function>(f), mat,
      std::integral_constant<unsigned int, I + 1>());
}

template <typename Function, typename Dest, unsigned int NI, unsigned int NJ,
          typename Blocks, typename Rhs>
void apply_function_to_diagonal_blocks(
    Function &&function, Dest &x, const MatrixReplacement<NI, NJ, Blocks> &mat,
    const Rhs &b) {
  apply_function_to_diagonal_blocks(std::forward<Function>(function), x, mat, b,
                                    std::integral_constant<unsigned int, 0>());
}

template <typename Function, unsigned int NI, unsigned int NJ, typename Blocks>
void apply_function_to_diagonal_blocks(
    Function &&function, const MatrixReplacement<NI, NJ, Blocks> &mat) {
  apply_function_to_diagonal_blocks(std::forward<Function>(function), mat,
                                    std::integral_constant<unsigned int, 0>());
}

template <typename T> struct storage_vector_type {
  T *m_data;
  size_t m_size;

  storage_vector_type() : m_data(nullptr), m_size(0) {}
  storage_vector_type(const storage_vector_type &other)
      : m_data(nullptr), m_size(0) {
    resize(other.m_size);
    for (size_t i = 0; i < m_size; ++i) {
      m_data[i] = other.m_data[i];
    }
  }

  storage_vector_type(storage_vector_type &&other) = default;

  ~storage_vector_type() {
    if (m_data)
      delete[] m_data;
  }

  storage_vector_type &operator=(const storage_vector_type &other) {
    resize(other.m_size);
    for (size_t i = 0; i < m_size; ++i) {
      m_data[i] = other.m_data[i];
    }
    return *this;
  }
  CUDA_HOST_DEVICE
  T *begin() const { return m_data; }
  CUDA_HOST_DEVICE
  T *end() const { return m_data + m_size; }
  CUDA_HOST_DEVICE
  size_t size() const { return m_size; }
  CUDA_HOST_DEVICE
  void resize(const size_t n) {
    // only supports reductions in size after initial resize
    ASSERT_CUDA(!m_data || n <= m_size);
    if (!m_data && n > 0) {
      m_data = new T[n];
    }
    m_size = n;
  }
  CUDA_HOST_DEVICE
  T &operator[](const int i) { return m_data[i]; }
  CUDA_HOST_DEVICE
  const T &operator[](const int i) const { return m_data[i]; }
};
} // namespace detail

#include "MultiLevelSchwartzPreconditioner.h"

template <typename Solver> class ChebyshevPreconditioner {
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
  int m_order;
  matrix_type m_col_Rn_matrix, m_row_Rn_matrix;
  solver_type m_factorized_matrix;
  Index m_rows;
  Index m_cols;
  mutable vector_type m_W;
  mutable vector_type m_fcheb;

public:
  typedef typename vector_type::StorageIndex StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic
  };

  ChebyshevPreconditioner() : m_isInitialized(false), m_order(10) {}

  template <typename MatType>
  explicit ChebyshevPreconditioner(const MatType &mat) {
    compute(mat);
  }

  Index rows() const { return m_rows; }
  Index cols() const { return m_cols; }

  void set_order(int arg) { m_order = arg; }

  template <typename Kernel>
  void analyze_impl_block(const Index start_row, const Kernel &kernel) {
    typedef typename Kernel::row_elements_type row_elements_type;
    typedef typename row_elements_type::query_type query_type;
    typedef typename query_type::traits_type traits_type;
    typedef typename traits_type::double_d double_d;
    typedef typename traits_type::int_d int_d;
    typedef typename traits_type::position position;
    const unsigned int dimension = query_type::dimension;

    const row_elements_type &rows = kernel.get_row_elements();
    const row_elements_type &cols = kernel.get_row_elements();

    auto row_bounds = rows.get_query().get_bounds();
    auto col_bounds = rows.get_query().get_bounds();

    detail::ChebyshevRn<dimension> row_Rn(m_order, row_bounds);
    detail::ChebyshevRn<dimension> col_Rn(m_order, col_bounds);

    const size_t ncheb = std::pow(m_order, dimension);
    LOG(2,
        "ChebyshevPreconditioner:analyzing_impl_block with ncheb = " << ncheb);

    const int_d start = int_d::Constant(0);
    const int_d end = int_d::Constant(m_order);

    // fill row_Rn matrix
    m_row_Rn_matrix.resize(rows.size(), ncheb);
    for (size_t i = 0; i < rows.size(); ++i) {
      row_Rn.set_position(get<position>(rows)[i]);
      lattice_iterator<dimension> mj(start, end);
      for (size_t j = 0; j < ncheb; ++j, ++mj) {
        m_row_Rn_matrix(i, j) = row_Rn(*mj);
      }
    }

    // fill col_Rn matrix
    m_col_Rn_matrix.resize(cols.size(), ncheb);
    for (size_t i = 0; i < cols.size(); ++i) {
      col_Rn.set_position(get<position>(rows)[i]);
      lattice_iterator<dimension> mj(start, end);
      for (size_t j = 0; j < ncheb; ++j, ++mj) {
        m_col_Rn_matrix(i, j) = col_Rn(*mj);
      }
    }

    // fill kernel matrix
    matrix_type kernel_matrix(ncheb, ncheb);
    lattice_iterator<dimension> mi(start, end);
    for (size_t i = 0; i < ncheb; ++i, ++mi) {
      const double_d pi = col_Rn.get_position(*mi);
      lattice_iterator<dimension> mj(start, end);
      for (size_t j = 0; j < ncheb; ++j, ++mj) {
        const double_d pj = row_Rn.get_position(*mj);
        kernel_matrix(i, j) = kernel.get_position_function()(pi, pj);
      }
    }

    m_factorized_matrix.compute(kernel_matrix);

    Eigen::VectorXd b = Eigen::VectorXd::Random(kernel_matrix.rows());
    Eigen::VectorXd x = m_factorized_matrix.solve(b);
    double relative_error = (kernel_matrix * x - b).norm() / b.norm();
    if (relative_error > 1e-3 || std::isnan(relative_error)) {
      std::cout << "relative error = " << relative_error << std::endl;
    }
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
  ChebyshevPreconditioner &
  analyzePattern(const MatrixReplacement<NI, NJ, Blocks> &mat) {

    LOG(2, "ChebyshevPreconditioner: analyze pattern");
    m_rows = mat.rows();
    m_cols = mat.cols();
    analyze_impl(mat, detail::make_index_sequence<NI>());
    return *this;
  }

  template <int _Options, typename _StorageIndex>
  ChebyshevPreconditioner &analyzePattern(
      const Eigen::SparseMatrix<Scalar, _Options, _StorageIndex> &mat) {
    CHECK(m_row_Rn_matrix.rows() > 0,
          "ChebyshevPreconditioner::analyzePattern(): cannot analyze sparse "
          "matrix, "
          "call analyzePattern using a Aboria MatrixReplacement class first");
    return *this;
  }

  template <int _Options, typename _StorageIndex, int RefOptions,
            typename RefStrideType>
  ChebyshevPreconditioner &
  analyzePattern(const Eigen::Ref<
                 const Eigen::SparseMatrix<Scalar, _Options, _StorageIndex>,
                 RefOptions, RefStrideType> &mat) {
    CHECK(m_row_Rn_matrix.rows() > 0,
          "ChebyshevPreconditioner::analyzePattern(): cannot analyze sparse "
          "matrix, "
          "call analyzePattern using a Aboria MatrixReplacement class first");
    return *this;
  }

  template <typename Derived>
  ChebyshevPreconditioner &
  analyzePattern(const Eigen::DenseBase<Derived> &mat) {
    CHECK(m_row_Rn_matrix.rows() > 0,
          "ChebyshevPreconditioner::analyzePattern(): cannot analyze dense "
          "matrix, "
          "call analyzePattern need to pass a Aboria MatrixReplacement class "
          "first");
    return *this;
  }

  template <typename MatType>
  ChebyshevPreconditioner &factorize(const MatType &mat) {
    LOG(2, "ChebyshevPreconditioner: factorizing domain");
    eigen_assert(
        static_cast<typename MatType::Index>(m_rows) == mat.rows() &&
        "ChebyshevPreconditioner::solve(): invalid number of rows of mat");
    eigen_assert(
        static_cast<typename MatType::Index>(m_cols) == mat.cols() &&
        "ChebyshevPreconditioner::solve(): invalid number of rows of mat");

    m_isInitialized = true;

    return *this;
  }

  template <typename MatType>
  ChebyshevPreconditioner &compute(const MatType &mat) {
    analyzePattern(mat);
    return factorize(mat);
  }

  /** \internal */
  template <typename Rhs, typename Dest>
  void _solve_impl(const Rhs &b, Dest &x) const {

    // First compute the weights at the Chebyshev nodes ym
    // by anterpolation
    m_W = m_col_Rn_matrix.transpose() * b;

    // Next compute f ðxÞ at the Chebyshev nodes xl:
    m_fcheb = m_factorized_matrix.solve(m_W);

    // Last compute f ðxÞ at the observation points xi by interpolation:
    x = m_row_Rn_matrix * m_fcheb;
  }

  template <typename Rhs>
  inline const Eigen::Solve<ChebyshevPreconditioner, Rhs>
  solve(const Eigen::MatrixBase<Rhs> &b) const {
    eigen_assert(
        static_cast<typename Rhs::Index>(m_rows) == b.rows() &&
        "ChebyshevPreconditioner::solve(): invalid number of rows of the "
        "right hand side matrix b");
    eigen_assert(m_isInitialized &&
                 "ChebyshevPreconditioner is not initialized.");
    return Eigen::Solve<ChebyshevPreconditioner, Rhs>(*this, b.derived());
  }

  Eigen::ComputationInfo info() { return Eigen::Success; }
}; // namespace Aboria

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

#if 0
template <unsigned int ReducedOrder, 
          typename InnerPreconditioner=Eigen::IncompleteLUT<double>,
          typename IterativeSolver=Eigen::DGMRES<Eigen::SparseMatrix<double>,
                                                InnerPreconditioner>>
class ExtMatrixPreconditioner {
    typedef double Scalar;
    typedef size_t Index;
    typedef Eigen::SparseMatrix<Scalar> sparse_matrix_type;
    typedef Eigen::Matrix<Scalar,Eigen::Dynamic,1> vector_type;
    typedef Eigen::Matrix<int,Eigen::Dynamic,1> index_vector_type;
    typedef InnerPreconditioner solver_type;
    Index m_rows;
    Index m_cols;
    size_t m_inner_iterations;
    std::vector<index_vector_type> m_col_maps;
    std::vector<size_t> m_col_sizes;
    std::vector<size_t> m_row_sizes;
    std::vector<index_vector_type> m_row_maps;
    std::vector<std::shared_ptr<solver_type>> m_solvers;
    std::vector<sparse_matrix_type> m_ext_matrices;
    std::vector<sparse_matrix_type> m_str_ext_matrices;

  public:
    typedef typename vector_type::StorageIndex StorageIndex;
    enum {
      ColsAtCompileTime = Eigen::Dynamic,
      MaxColsAtCompileTime = Eigen::Dynamic
    };

    ExtMatrixPreconditioner():m_inner_iterations(10)
    {}

    template<typename MatType>
    explicit ExtMatrixPreconditioner(const MatType& mat):m_inner_iterations(10) {
      compute(mat);
    }

    Index rows() const { return m_rows; }
    Index cols() const { return m_cols; }

    template<unsigned int NI, unsigned int NJ, typename Blocks>
    ExtMatrixPreconditioner& analyzePattern(const MatrixReplacement<NI,NJ,Blocks>& mat)
    {

        LOG(2,"ExtMatrixPreconditioner: analyze pattern");
        m_rows = mat.rows();
        m_cols = mat.cols();

        return *this;
    }
    

    struct factorize_block {
        std::vector<size_t> &m_col_sizes; 
        std::vector<size_t> &m_row_sizes; 
        std::vector<index_vector_type> &m_col_maps; 
        std::vector<index_vector_type> &m_row_maps; 
        std::vector<std::shared_ptr<solver_type>> &m_solvers; 
        std::vector<sparse_matrix_type> &m_ext_matrices; 
        std::vector<sparse_matrix_type> &m_str_ext_matrices; 
        size_t m_inner_iterations;
        int i;

        factorize_block(std::vector<size_t>& col_sizes,
                        std::vector<size_t>& row_sizes,
                        std::vector<index_vector_type>& col_maps,
                        std::vector<index_vector_type>& row_maps,
                        std::vector<std::shared_ptr<solver_type>>& solvers,
                        std::vector<sparse_matrix_type>& ext_matrices,
                        std::vector<sparse_matrix_type>& str_ext_matrices,
                        size_t inner_iterations):
            m_col_sizes(col_sizes),m_row_sizes(row_sizes),m_col_maps(col_maps),m_row_maps(row_maps),m_solvers(solvers),m_ext_matrices(ext_matrices),m_str_ext_matrices(str_ext_matrices),m_inner_iterations(inner_iterations),i(0) {}


        template <typename Block>
        void operator()(const Block& block) {
            LOG(2,"ExtMatrixPreconditioner: block "<<i<<": non h2 block");
            m_solvers[i] = nullptr;
            m_col_sizes[i] = block.cols();
            m_row_sizes[i] = block.rows();
            ++i;
        }

        template <typename RowParticles, typename ColParticles, typename PositionF>
        void operator()(const KernelH2<RowParticles,ColParticles,PositionF>& kernel) {
            static const unsigned int dimension = RowParticles::dimension;

            auto h2 = make_h2_matrix(kernel.get_row_particles(),
                                     kernel.get_col_particles(),
                       make_black_box_expansion<dimension,ReducedOrder>(
                                    kernel.get_position_function()));
            LOG(2,"ExtMatrixPreconditioner: block "<<i<<": generate extended matrix");
            m_ext_matrices[i] = h2.gen_extended_matrix();
            m_str_ext_matrices[i] = h2.gen_stripped_extended_matrix();
            /*
            Eigen::saveMarket(m_ext_matrices[i],"ext_matrix.mat");
            Eigen::saveMarket(m_str_ext_matrices[i],"str_ext_matrix.mat");
            std::ofstream myfile;
            myfile.open ("ext_matrix.csv");
            myfile << Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>(m_ext_matrices[i]);
            myfile.close();
            myfile.open ("str_ext_matrix.csv");
            myfile << Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>(m_str_ext_matrices[i]);
            myfile.close();
            */

            m_col_sizes[i] = kernel.cols();
            m_row_sizes[i] = kernel.rows();
            m_col_maps[i] = h2.gen_column_map();
            m_row_maps[i] = h2.gen_row_map();
            LOG(2,"ExtMatrixPreconditioner: block "<<i<<": set precon");
            //LOG(2,"ExtMatrixPreconditioner: block "<<i<<": create solver");
            //m_solvers[i] = std::make_shared<solver_type>(m_ext_matrices[i]);
            m_solvers[i] = std::make_shared<solver_type>();
            m_solvers[i]->setDroptol(0.1);
            m_solvers[i]->compute(m_str_ext_matrices[i]);
            if (m_solvers[i]->info() != Eigen::Success) {
                ERROR("ExtMatrixPreconditioner inner preconditioner could not factorize");
            }
            //m_solvers[i]->setMaxIterations(m_inner_iterations);
            //LOG(2,"ExtMatrixPreconditioner: block "<<i<<": set precon");
            //m_solvers[i]->preconditioner().setDroptol(0.1);
            //m_solvers[i]->preconditioner().compute(m_str_ext_matrices[i]);
            LOG(2,"ExtMatrixPreconditioner: block "<<i<<": factorization complete");
            ++i;
        }
    };


    template <unsigned int NI, unsigned int NJ, typename Blocks>
    ExtMatrixPreconditioner& factorize(const MatrixReplacement<NI,NJ,Blocks>& mat)
    {
        LOG(2,"ExtMatrixPreconditioner: factorizing domain");

        m_rows = mat.rows();
        m_cols = mat.cols();
        m_solvers.resize(NI);
        m_ext_matrices.resize(NI);
        m_str_ext_matrices.resize(NI);
        m_col_sizes.resize(NI);
        m_row_sizes.resize(NI);
        m_col_maps.resize(NI);
        m_row_maps.resize(NI);
        detail::apply_function_to_diagonal_blocks(
                factorize_block(m_col_sizes,m_row_sizes,m_col_maps,m_row_maps,
                                m_solvers,m_ext_matrices,m_str_ext_matrices,m_inner_iterations), 
                mat);
        m_isInitialized = true;

        return *this;
    }
    
    template<typename MatType>
    ExtMatrixPreconditioner& compute(const MatType& mat)
    {
        analyzePattern(mat);
        return factorize(mat);
    }

    template<typename Rhs, typename Dest>
    void _solve_impl(const Rhs& b, Dest& x) const {
        vector_type m_ext_b;
        vector_type m_ext_x;
        size_t row = 0;
        size_t col = 0;

        for (size_t i = 0; i < m_solvers.size(); ++i) {
            if (m_solvers[i] != nullptr) { // solver only exists for h2 blocks
                LOG(2,"ExtMatrixPreconditioner: block "<<i<<" solve");
                // construct b
                m_ext_b.resize(m_ext_matrices[i].rows());
                for (size_t j = 0; j < m_row_maps[i].size(); ++j) {
                    m_ext_b[m_row_maps[i][j]] = b[row + j];
                }
                for (int j = m_row_maps[i].size(); j < m_ext_matrices[i].rows(); ++j) {
                    m_ext_b[j] = 0;
                }
                m_ext_x = m_ext_b;
                //for (int j = 0; j < m_ext_matrices[i].cols(); ++j) {
                //    m_ext_x[j] = 0;
                //}
                //for (int j = 0; j < m_ext_matrices[i].rows(); ++j) {
                //    std::cout << "ext_b["<<j<<"] = "<< m_ext_b[j] << std::endl;
                //}

                //TODO: solve with guess?
                //m_ext_x = m_solvers[i]->solveWithGuess(m_ext_b,m_ext_x);
                //m_ext_x = m_solvers[i]->solve(m_ext_b);
                Eigen::Index iters = m_inner_iterations; 
                double tol_error = 1e-10;
                std::cout << "m_ext_b norm = "<<m_ext_b.norm() << std::endl;
                std::cout << "m_ext_x norm = "<<m_ext_x.norm() << std::endl;
                std::cout << "residual norm = "<<(m_ext_matrices[i]*m_ext_x-m_ext_b).norm() << std::endl;
                Eigen::internal::gmres(m_ext_matrices[i],m_ext_b,m_ext_x,*m_solvers[i],
                        iters,2*m_inner_iterations,tol_error);
                //LOG(2,"ExtMatrixPreconditioner: solve complete: #iterations: " << m_solvers[i]->iterations() << ", estimated error: " << m_solvers[i]->error() << " true error = "<<(m_ext_matrices[i]*m_ext_x-m_ext_b).norm());
                LOG(2,"ExtMatrixPreconditioner: solve complete: #iterations: " << iters << " true error = "<<(m_ext_matrices[i]*m_ext_x-m_ext_b).norm());

                // filter to x
                for (size_t j = 0; j < m_col_maps[i].size(); ++j) {
                    x[col + j] = m_ext_x[m_col_maps[i][j]];
                }

                // increment row/col by number of particles (NOT size of ext vectors)
                row += m_row_maps[i].size();
                col += m_col_maps[i].size();
            } else {
                LOG(2,"ExtMatrixPreconditioner: block "<<i<<" non h2 block");
                for (int j = 0; j < m_col_sizes[i]; ++j) {
                    x[col + j] = b[row + j];
                }
                // increment row/col by the size of the block
                row += m_row_sizes[i];
                col += m_col_sizes[i];
            }
            
        }
        LOG(2,"ExtMatrixPreconditioner: done solve_impl");
    }
    

    template<typename Rhs> 
    inline const Eigen::Solve<ExtMatrixPreconditioner, Rhs>
    solve(const Eigen::MatrixBase<Rhs>& b) const {
        eigen_assert(m_rows==b.rows()
                && "ExtMatrixPreconditioner::solve(): invalid number of rows of the right hand side matrix b");
        eigen_assert(m_isInitialized 
                && "ExtMatrixPreconditioner is not initialized.");
        return Eigen::Solve<ExtMatrixPreconditioner, Rhs>(*this, b.derived());
    }
    
    Eigen::ComputationInfo info() { return Eigen::Success; }

  protected:
    bool m_isInitialized;
};
#endif
class CardinalFunctionsPreconditioner {
  typedef double Scalar;
  typedef size_t Index;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> matrix_type;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> vector_type;
  typedef std::vector<size_t> storage_vector_type;
  typedef std::vector<storage_vector_type> connectivity_type;

protected:
  bool m_isInitialized;

private:
  size_t m_random;
  double m_sigma;
  double m_M;

  connectivity_type m_domain_buffer;
  std::vector<vector_type> m_weights;
  Index m_rows;
  Index m_cols;

public:
  typedef typename vector_type::StorageIndex StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic
  };

  CardinalFunctionsPreconditioner()
      : m_isInitialized(false), m_random(0), m_sigma(-1), m_M(1.0) {}

  template <typename MatType>
  explicit CardinalFunctionsPreconditioner(const MatType &mat) {
    compute(mat);
  }

  Index rows() const { return m_rows; }
  Index cols() const { return m_cols; }

  void set_number_of_random_particles(size_t n) { m_random = n; }
  void set_sigma(double value) { m_sigma = value; }
  void set_rejection_sampling_scale(double value) { m_M = value; }

  template <typename Kernel>
  void analyze_impl_block(const Index start_row, const Kernel &kernel) {
    typedef typename Kernel::row_elements_type row_elements_type;
    typedef typename Kernel::col_elements_type col_elements_type;
    typedef typename row_elements_type::query_type query_type;
    typedef typename query_type::traits_type traits_type;
    typedef typename query_type::child_iterator child_iterator;
    typedef typename traits_type::double_d double_d;
    typedef typename traits_type::int_d int_d;
    typedef typename traits_type::position position;

    static_assert(
        std::is_same<row_elements_type, col_elements_type>::value,
        "Cardinal Functions preconditioner restricted to identical row and col "
        "particle sets");
    const row_elements_type &a = kernel.get_row_elements();
    CHECK(
        &a == &(kernel.get_col_elements()),
        "Cardinal Functions preconditioner restricted to identical row and col "
        "particle "
        "sets");
    const query_type &query = a.get_query();
    std::default_random_engine generator;
    m_domain_buffer.resize(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
      storage_vector_type &buffer = m_domain_buffer[i];

      // add buffer particles through random sampling
      int nspecial = std::pow(3, query_type::dimension);
      // const int nspecial = 0;
      buffer.resize(m_random + nspecial);
      std::vector<child_iterator> buckets(m_random + nspecial);
      std::uniform_real_distribution<double> uniform(0, 1);
      std::normal_distribution<double> normal(0, m_sigma);

      // add special points
      // updates nspecial with actual number of special points
      lattice_iterator<query_type::dimension> special_it(
          int_d::Constant(0), int_d::Constant(3), int_d::Constant(0));
      for (nspecial = 0; special_it != false; ++special_it, ++nspecial) {
        const double_d &bmin = query.get_bounds().bmin;
        const double_d &bmax = query.get_bounds().bmax;
        const double_d pos =
            (*special_it) *
                (0.5 * (bmax - bmin) - std::numeric_limits<double>::epsilon()) +
            bmin;
        // std::cout <<"adding special point at pos = "<<pos<<std::endl;
        buckets[nspecial] = query.get_bucket(pos);
      }

      const double_d middle = get<position>(a)[i];
      if (m_sigma > 0) {
        const double scale2 = 1.0 / std::pow(m_sigma, 2);
        auto gaussianf = [&](const double_d &x) {
          return std::exp(-(x - middle).squaredNorm() * scale2);
        };
        std::generate(buckets.begin() + nspecial, buckets.end(), [&]() {
          double_d sp;
          bool accepted;
          do {
            for (size_t i = 0; i < query_type::dimension; i++) {
              sp[i] = normal(generator) + middle[i];
            }
            const bool in_domain =
                (sp < a.get_max()).all() && (sp >= a.get_min()).all();
            accepted =
                in_domain && uniform(generator) <
                                 kernel.get_position_function()(middle, sp) /
                                     (gaussianf(sp) * m_M);
          } while (!accepted);
          return query.get_bucket(sp);
        });
      } else {
        const double volume = (a.get_max() - a.get_min()).prod();
        std::generate(buckets.begin() + nspecial, buckets.end(), [&]() {
          double_d sp;
          bool accepted;
          do {
            for (size_t i = 0; i < query_type::dimension; i++) {
              sp[i] =
                  0.5 * (a.get_max()[i] - a.get_min()[i]) * uniform(generator) +
                  a.get_min()[i];
            }
            const bool in_domain =
                (sp < a.get_max()).all() && (sp >= a.get_min()).all();
            accepted =
                in_domain &&
                uniform(generator) <
                    kernel.get_position_function()(middle, sp) * volume / m_M;
          } while (!accepted);
          return query.get_bucket(sp);
        });
      }

      std::unordered_map<size_t, std::pair<child_iterator, size_t>> counts;
      for (int i = 0; i < buckets.size(); ++i) {
        auto bucket_index = query.get_bucket_index(*(buckets[i]));
        auto it = counts.find(bucket_index);
        if (it != counts.end()) {
          it->second.second++;
        } else {
          counts[bucket_index] = std::make_pair(buckets[i], 1);
        }
      }
      // for (auto i : counts) {
      // std::cout << "bucket index " << i.first << " with bounds "
      //<< query.get_bounds(i.second.first) << " has " << i.second.second
      //<< " counts" << std::endl;
      //}

      int out_index = 0;
      std::for_each(counts.begin(), counts.end(), [&](auto i) {
        auto ci = i.second.first;
        size_t count = i.second.second;
        auto pit = query.get_bucket_particles(*ci);
        auto num_particles = pit.distance_to_end();
        std::vector<int> bucket_indices(num_particles);
        std::iota(bucket_indices.begin(), bucket_indices.end(), 0);
        std::shuffle(bucket_indices.begin(), bucket_indices.end(), generator);
        const int trunc_count = std::min(count, bucket_indices.size());
        std::transform(
            bucket_indices.begin(), bucket_indices.begin() + trunc_count,
            buffer.begin() + out_index, [&](const int i) {
              return (&get<position>(*(pit + i)) - &get<position>(a)[0]) +
                     start_row;
            });
        // std::cout << "looking for " << count
        //<< " samples in buffer. Found at indicies ";
        // for (size_t i = out_index; i < out_index + trunc_count; ++i) {
        //  std::cout << buffer[i] << " ";
        //}
        // std::cout << std::endl;
        out_index += trunc_count;
      });
      buffer.resize(out_index);

      // ensure that cardinal index isn't in buffer
      std::remove_if(buffer.begin(), buffer.end(),
                     [&i](const int j) { return j == i; });

#ifdef HAVE_CAIRO
      const int image_size = 512;
      cairo_surface_t *surface = cairo_svg_surface_create(
          ("sampler" + std::to_string(i) + ".svg").c_str(), image_size,
          image_size);
      cairo_svg_surface_restrict_to_version(surface, CAIRO_SVG_VERSION_1_2);
      cairo_t *cr = cairo_create(surface);
      const double lw = 0.01;

      cairo_scale(cr, image_size, image_size);
      cairo_set_line_width(cr, lw);
      const double PI = boost::math::constants::pi<double>();
      cairo_set_source_rgba(cr, 0.5, 0, 0, 1.0);

      auto &pos = get<position>(a)[i];
      cairo_arc(cr, pos[0], pos[1], lw, 0, 2 * PI);
      cairo_fill(cr);

      cairo_set_source_rgba(cr, 0, 0, 0.5, 1.0);
      for (auto i : buffer) {
        auto &pos = get<position>(a)[i];
        cairo_arc(cr, pos[0], pos[1], lw, 0, 2 * PI);
        cairo_fill(cr);
      }
      cairo_destroy(cr);
      cairo_surface_destroy(surface);
#endif

      ASSERT(buffer.size() > 0, "no particles in buffer");
    }
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
  CardinalFunctionsPreconditioner &
  analyzePattern(const MatrixReplacement<NI, NJ, Blocks> &mat) {

    LOG(2, "CardinalFunctionsPreconditioner: analyze pattern");
    m_rows = mat.rows();
    m_cols = mat.cols();
    analyze_impl(mat, detail::make_index_sequence<NI>());

    int minsize_buffer = 1000;
    int maxsize_buffer = 0;
    for (size_t domain_index = 0; domain_index < m_domain_buffer.size();
         ++domain_index) {
      const int size_buffer = m_domain_buffer[domain_index].size();
      if (size_buffer < minsize_buffer)
        minsize_buffer = size_buffer;
      if (size_buffer > maxsize_buffer)
        maxsize_buffer = size_buffer;
    }
    LOG(2, "CardinalFunctionsPreconditioner: finished analysis, found "
               << m_domain_buffer.size() << " domains, with " << minsize_buffer
               << "--" << maxsize_buffer << " buffer particles")
    return *this;
  }

  template <int _Options, typename _StorageIndex>
  CardinalFunctionsPreconditioner &analyzePattern(
      const Eigen::SparseMatrix<Scalar, _Options, _StorageIndex> &mat) {
    CHECK(m_domain_buffer.size() > 0,
          "CardinalFunctionsPreconditioner::analyzePattern(): cannot analyze "
          "sparse "
          "matrix, "
          "call analyzePattern using a Aboria MatrixReplacement class first");
    return *this;
  }

  template <int _Options, typename _StorageIndex, int RefOptions,
            typename RefStrideType>
  CardinalFunctionsPreconditioner &
  analyzePattern(const Eigen::Ref<
                 const Eigen::SparseMatrix<Scalar, _Options, _StorageIndex>,
                 RefOptions, RefStrideType> &mat) {
    CHECK(m_domain_buffer.size() > 0,
          "CardinalFunctionsPreconditioner::analyzePattern(): cannot analyze "
          "sparse "
          "matrix, "
          "call analyzePattern using a Aboria MatrixReplacement class first");
    return *this;
  }

  template <typename Derived>
  CardinalFunctionsPreconditioner &
  analyzePattern(const Eigen::DenseBase<Derived> &mat) {
    CHECK(m_domain_buffer.size() > 0,
          "CardinalFunctionsPreconditioner::analyzePattern(): cannot analyze "
          "dense "
          "matrix, "
          "call analyzePattern need to pass a Aboria MatrixReplacement class "
          "first");
    return *this;
  }

  template <typename MatType>
  CardinalFunctionsPreconditioner &factorize(const MatType &mat) {
    LOG(2, "CardinalFunctionsPreconditioner: factorizing domain");
    eigen_assert(static_cast<typename MatType::Index>(m_rows) == mat.rows() &&
                 "CardinalFunctionsPreconditioner::solve(): invalid number of "
                 "rows of mat");
    eigen_assert(static_cast<typename MatType::Index>(m_cols) == mat.cols() &&
                 "CardinalFunctionsPreconditioner::solve(): invalid number of "
                 "rows of mat");

    matrix_type domain_matrix;
    const size_t N = mat.rows();
    m_weights.resize(m_domain_buffer.size());
    for (size_t domain_index = 0; domain_index < m_domain_buffer.size();
         ++domain_index) {
      const storage_vector_type &buffer = m_domain_buffer[domain_index];
      vector_type &weights = m_weights[domain_index];

      const size_t size = 1 + buffer.size();
      // std::cout << "domain "<<domain_index<<"indicies =
      // "<<indicies.size()<<" buffer =  "<<buffer.size()<<" random =
      // "<<random.size()<<std::endl;

      domain_matrix.resize(size, N);
      weights.resize(size);

      size_t i = 0;
      for (const size_t &big_index_i : {domain_index}) {
        for (size_t j = 0; j < N; ++j) {
          domain_matrix(i, j) = mat.coeff(big_index_i, j);
        }
        ++i;
      }
      for (const size_t &big_index_i : buffer) {
        for (size_t j = 0; j < N; ++j) {
          domain_matrix(i, j) = mat.coeff(big_index_i, j);
        }
        ++i;
      }

      vector_type b = vector_type::Zero(N);
      b[domain_index] = 1;

      weights = domain_matrix.transpose().colPivHouseholderQr().solve(b);
      // weights = domain_matrix.transpose()
      //              .bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
      //              .solve(b);
      double relative_error =
          (domain_matrix.transpose() * weights - b).norm() / b.norm();
      if (relative_error > 1e-3 || std::isnan(relative_error)) {
        std::cout << "domain index = " << domain_index
                  << ": relative error = " << relative_error << std::endl;
      }
    }

    m_isInitialized = true;

    return *this;
  }

  template <typename MatType>
  CardinalFunctionsPreconditioner &compute(const MatType &mat) {
    analyzePattern(mat);
    return factorize(mat);
  }

  /** \internal */
  template <typename Rhs, typename Dest>
  void _solve_impl(const Rhs &b, Dest &x) const {
    for (size_t i = 0; i < m_domain_buffer.size(); ++i) {
      const storage_vector_type &buffer = m_domain_buffer[i];
      const vector_type &weights = m_weights[i];

      // x = W * b
      x[i] = weights[0] * b[i];
      for (size_t j = 0; j < buffer.size(); ++j) {
        x[i] += weights[j] * b[buffer[j]];
      }
    }
  }

  template <typename Rhs>
  inline const Eigen::Solve<CardinalFunctionsPreconditioner, Rhs>
  solve(const Eigen::MatrixBase<Rhs> &b) const {
    eigen_assert(static_cast<typename Rhs::Index>(m_rows) == b.rows() &&
                 "CardinalFunctionsPreconditioner::solve(): invalid number of "
                 "rows of the "
                 "right hand side matrix b");
    eigen_assert(m_isInitialized &&
                 "CardinalFunctionsPreconditioner is not initialized.");
    return Eigen::Solve<CardinalFunctionsPreconditioner, Rhs>(*this,
                                                              b.derived());
  }

  Eigen::ComputationInfo info() { return Eigen::Success; }
}; // namespace Aboria

template <typename Solver> class SchwartzPreconditioner {
  typedef double Scalar;
  typedef size_t Index;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> matrix_type;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> vector_type;
  typedef Solver solver_type;

  typedef detail::storage_vector_type<size_t> storage_vector_type;

  typedef std::vector<storage_vector_type> connectivity_type;
  typedef std::vector<solver_type> solvers_type;
  typedef std::vector<matrix_type> matrices_type;

protected:
  bool m_isInitialized;

private:
  size_t m_max_buffer_n;
  mutable std::array<double, 2> m_timing;

  connectivity_type m_indicies[2];
  connectivity_type m_buffer[2];
  solvers_type m_factorized_matrix[2];
  matrices_type m_matrix[2];
  matrices_type m_coupling_matrix[2];
  matrix_type m_A;
  int m_multiplicitive;
  int m_levels;
  bool m_interpolate;
  bool m_use_root;

  Index m_rows;
  Index m_cols;

public:
  typedef typename vector_type::StorageIndex StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic
  };

  SchwartzPreconditioner()
      : m_isInitialized(false), m_max_buffer_n(300), m_timing{0, 0},
        m_multiplicitive(0), m_levels(1), m_interpolate(false),
        m_use_root(true) {}

  template <typename MatType>
  explicit SchwartzPreconditioner(const MatType &mat)
      : m_isInitialized(false), m_max_buffer_n(300), m_timing{0, 0},
        m_multiplicitive(0), m_levels(1), m_interpolate(false),
        m_use_root(true) {
    compute(mat);
  }

  Index rows() const { return m_rows; }
  Index cols() const { return m_cols; }

  void set_max_buffer_n(size_t arg) { m_max_buffer_n = arg; }
  void set_multiplicative(int arg) { m_multiplicitive = arg; }
  void set_interpolate(bool arg) { m_interpolate = arg; }
  void set_levels(int arg) { m_levels = arg; }
  void set_use_root(int arg) { m_use_root = arg; }
  const std::array<double, 2> &get_timing() { return m_timing; }

  template <typename MatType>
  SchwartzPreconditioner &analyzePattern(const MatType &mat) {
    LOG(2, "SchwartzPreconditioner: analyzePattern: do nothing");
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
    matrix_type *m_coupling_matrix;
    size_t m_max_buffer_n;
    size_t m_max_bucket_n;
    size_t m_start_row;
    const T m_function;
    Query m_query;
    int m_level;
    bool is_root;
    bool m_interpolate;

    CUDA_HOST_DEVICE
    void operator()(const int i) const {
      auto &ci = m_nodes[i];
      const auto &bounds = m_query.get_bounds(ci);
      LOG(3, "process_node with bounds "
                 << bounds << " is leaf node = " << m_query.is_leaf_node(*ci)
                 << " is root = " << is_root);

      const double_d middle = 0.5 * (bounds.bmax + bounds.bmin);
      const double_d side = 0.9 * (bounds.bmax - bounds.bmin);

      // skip over empty leaf buckets on leaf level, or leaf buckets on non leaf
      // level
      if ((m_level == 0 && m_query.get_bucket_particles(*ci) == false) ||
          (m_level > 0 && m_query.is_leaf_node(*ci))) {
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
          const auto &bucket_bounds =
              m_query.get_bounds(bucket.get_child_iterator());
          if (bucket_bounds <= bounds) {
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
          indicies_tmp[i] = m_start_row + i;
        }
        i_indicies = n_indicies;
      } else {
        for (auto bucket =
                 m_query.template get_buckets_near_point<-1>(middle, side);
             bucket != false; ++bucket) {
          const auto &bucket_bounds =
              m_query.get_bounds(bucket.get_child_iterator());
          const bool is_in_bounds = (bucket_bounds <= bounds);
          for (auto particle = m_query.get_bucket_particles(*bucket);
               particle != false; ++particle) {
            const double_d &p = get<position>(*particle);
            const size_t index =
                &p - get<position>(m_query.get_particles_begin());
            if (is_in_bounds) {
              indicies_tmp[i_indicies++] = m_start_row + index;
            } else {
              buffer_tmp[i_buffer++] = m_start_row + index;
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

      if (m_interpolate) {
        // just do buffer
        for (size_t i = 0; i < indicies_size; ++i) {
          buffer[i] = indicies_tmp[i];
        }
      } else {
        // do buffer and indicies
        indicies.resize(indicies_size);
        for (size_t i = 0; i < indicies_size; ++i) {
          buffer[i] = indicies_tmp[i];
          indicies[i] = indicies_tmp[i];
        }
      }

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

      if (m_level > 0 && m_interpolate) {
        LOG(3, "doing interpolate coupling");

        // fill coupling matrix
        int count = 0;
        if (is_root) {
          count = m_query.number_of_particles();
        } else {
          for (auto child = m_query.get_subtree(ci); child != false; ++child) {
            if (m_query.is_leaf_node(*child)) {
              count += m_query.get_bucket_particles(*child).distance_to_end();
            }
          }
        }

        auto &coupling_matrix = m_coupling_matrix[i];
        indicies.resize(count);
        coupling_matrix.resize(count, buffer.size());

        if (is_root) {
          for (int i = 0; i < count; ++i) {
            indicies[i] = i;
            auto a = m_query.get_particles_begin()[i];
            int mj = 0;
            for (const size_t &big_index_j : buffer) {
              auto b = m_query.get_particles_begin()[big_index_j];
              coupling_matrix(i, mj++) = m_function(a, b);
            }
          }
        } else {
          int mi = 0;
          for (auto child = m_query.get_subtree(ci); child != false; ++child) {
            if (m_query.is_leaf_node(*child)) {
              for (auto a = m_query.get_bucket_particles(*child); a != false;
                   ++a, ++mi) {
                const double_d &p = get<position>(*a);
                indicies[mi] =
                    &p - get<position>(m_query.get_particles_begin());
                int mj = 0;
                for (const size_t &big_index_j : buffer) {
                  auto b = m_query.get_particles_begin()[big_index_j];
                  coupling_matrix(mi, mj++) = m_function(*a, b);
                }
              }
            }
          }
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

  template <typename Kernel>
  void factorize_fill_indices_and_matrices(
      const Index start_row, const Kernel &kernel,
      const typename Kernel::row_elements_type &a, connectivity_type &indicies,
      connectivity_type &buffer, matrices_type &matrix,
      const typename Kernel::row_elements_type &a_fine,
      matrices_type &coupling_matrix, connectivity_type &coupling_indicies) {}

  template <typename Kernel>
  void factorize_impl_block(const Index start_row, const Kernel &kernel) {
    typedef typename Kernel::row_elements_type row_elements_type;
    typedef typename Kernel::col_elements_type col_elements_type;
    typedef typename row_elements_type::query_type query_type;
    typedef typename query_type::traits_type traits_type;

    static_assert(std::is_same<row_elements_type, col_elements_type>::value,
                  "Schwartz preconditioner restricted to identical row and col "
                  "particle sets");
    static_assert(Kernel::BlockRows == 1 && Kernel::BlockCols == 1,
                  "Schwartz preconditioner not currently implemented for "
                  "vector-valued kernels");

    const row_elements_type &a = kernel.get_row_elements();

    CHECK(&a == &(kernel.get_col_elements()),
          "Schwartz preconditioner restricted to identical row and col "
          "particle "
          "sets");

    const query_type &query = a.get_query();

    // get list of leaf cells
    auto df_search = query.breadth_first();
    typename traits_type::template vector<
        typename std::remove_cv<typename decltype(df_search)::value_type>::type>
        tree = {*df_search};
    for (; df_search != false; ++df_search) {
      tree.push_back(*df_search);
    }
    auto leaf_nodes = *df_search;
    const size_t leaf_size = leaf_nodes.size();
    decltype(leaf_nodes) coarse_nodes;
    if (m_use_root) {
      coarse_nodes.push_back(query.get_root());
    } else {
      int level = tree.size() - 1;
      while (tree[level].size() > leaf_size / 1.9 && level > 0) {
        --level;
      }
      coarse_nodes = tree[level];
    }

    decltype(leaf_nodes) both_nodes[2] = {leaf_nodes, coarse_nodes};

    for (int i = 0; i < std::min(m_levels, 2); ++i) {
      auto &indicies = m_indicies[i];
      auto &buffer = m_buffer[i];
      auto &matrix = m_matrix[i];
      auto &coupling_matrix = m_coupling_matrix[i];
      auto &factorized_matrix = m_factorized_matrix[i];
      auto &nodes = both_nodes[i];

      indicies.resize(nodes.size());
      buffer.resize(nodes.size());
      matrix.resize(nodes.size());
      factorized_matrix.resize(nodes.size());
      if (i > 0) {
        coupling_matrix.resize(nodes.size());
      }
      const int num_threads = omp_get_max_threads();
      const int mult_buffer = m_max_buffer_n *
                              std::pow(a.get_max_bucket_size(), -2.0 / 3.0) *
                              std::pow(a.size() / num_threads, -1.0 / 3.0);
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
                iterator_to_raw_pointer(tmp_matrix.begin()),
                iterator_to_raw_pointer(tmp_coupling_matrix.begin()),
                (nodes.size() == 1)
                    ? m_max_buffer_n
                    : a.get_max_bucket_size() * (mult_buffer - 1),
                a.get_max_bucket_size(), start_row,
                kernel.get_kernel_function(), query, i, nodes.size() == 1,
                m_interpolate});

        detail::copy(tmp_indicies.begin(), tmp_indicies.end(),
                     indicies.begin());
        detail::copy(tmp_buffer.begin(), tmp_buffer.end(), buffer.begin());
        detail::copy(tmp_matrix.begin(), tmp_matrix.end(), matrix.begin());
        detail::copy(tmp_coupling_matrix.begin(), tmp_coupling_matrix.end(),
                     coupling_matrix.begin());
      } else {
        // no copy required
        auto count = traits_type::make_counting_iterator(0);
        detail::for_each(
            count, count + nodes.size(),
            process_node<typename Kernel::function_type, query_type>{
                iterator_to_raw_pointer(nodes.begin()),
                iterator_to_raw_pointer(indicies.begin()),
                iterator_to_raw_pointer(buffer.begin()),
                iterator_to_raw_pointer(matrix.begin()),
                iterator_to_raw_pointer(coupling_matrix.begin()),
                (nodes.size() == 1)
                    ? m_max_buffer_n
                    : a.get_max_bucket_size() * (mult_buffer - 1),
                a.get_max_bucket_size(), start_row,
                kernel.get_kernel_function(), query, i, nodes.size() == 1,
                m_interpolate});
      }

      // factorize domain matrices
      detail::transform(matrix.begin(), matrix.end(), factorized_matrix.begin(),
                        factorize_matrix());
    }
  }

  template <typename RowParticles, typename ColParticles>
  void
  factorize_impl_block(const Index start_row,
                       const KernelZero<RowParticles, ColParticles> &kernel) {}

  template <unsigned int NI, unsigned int NJ, typename Blocks, std::size_t... I>
  void factorize_impl(const MatrixReplacement<NI, NJ, Blocks> &mat,
                      detail::index_sequence<I...>) {
    int dummy[] = {0, (factorize_impl_block(mat.template start_row<I>(),
                                            std::get<I * NJ + I>(mat.m_blocks)),
                       0)...};
    static_cast<void>(dummy);
  }

  template <unsigned int NI, unsigned int NJ, typename Blocks>
  SchwartzPreconditioner &
  factorize(const MatrixReplacement<NI, NJ, Blocks> &mat) {
    LOG(2, "SchwartzPreconditioner: factorize");
    m_rows = mat.rows();
    m_cols = mat.cols();
    factorize_impl(mat, detail::make_index_sequence<NI>());

    if (m_multiplicitive != 0) {
      m_A.resize(mat.rows(), mat.cols());
      mat.assemble(m_A);
    }

    for (int i = 0; i < std::min(m_levels, 2); ++i) {
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

      LOG(2, "SchwartzPreconditioner: finished factorizing, on level "
                 << i << " found " << m_indicies[i].size() << " domains, with "
                 << min_indicies << "--" << max_indicies << " particles ("
                 << count << " total), and " << min_buffer << "--" << max_buffer
                 << " buffer particles. kernel matrix is size (" << m_A.rows()
                 << ',' << m_A.cols() << ')');
    }

    m_isInitialized = true;
    return *this;
  }

  template <int _Options, typename _StorageIndex, int RefOptions,
            typename RefStrideType>
  SchwartzPreconditioner &
  factorize(const Eigen::Ref<
            const Eigen::SparseMatrix<Scalar, _Options, _StorageIndex>,
            RefOptions, RefStrideType> &mat) {
    CHECK(m_indicies[0].size() > 0,
          "SchwartzPreconditioner::factorize(): cannot factorize sparse "
          "matrix, call factorize with a Aboria MatrixReplacement class "
          "instead");
    return *this;
  }

  template <typename Derived>
  SchwartzPreconditioner &factorize(const Eigen::DenseBase<Derived> &mat) {
    CHECK(m_indicies[0].size() > 0,
          "SchwartzPreconditioner::analyzePattern(): cannot factorize dense "
          "matrix, call factorize with a Aboria MatrixReplacement class "
          "instead");
    return *this;
  }

  template <typename MatType>
  SchwartzPreconditioner &compute(const MatType &mat) {
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
    bool m_interpolate;

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

      if (m_level == 0 || !m_interpolate) {
        // copy accumulate x values to big vector
        for (size_t j = 0; j < indicies.size(); ++j) {
          x[indicies[j]] += domain_x[j];
        }
      } else {
        // interpolate x values to big vector
        if (is_root) {
          domain_b = domain_factorized_matrix[i].solve(domain_x);
          // Eigen::Map<vector_type> big_x(x, m_coupling_matrix[i].rows());
          x += m_coupling_matrix[i] * domain_b;
        } else {
          domain_b = domain_factorized_matrix[i].solve(domain_x);
          vector_type big_x = m_coupling_matrix[i] * domain_b;
          for (size_t j = 0; j < indicies.size(); ++j) {
            x[indicies[j]] += big_x[j];
          }
        }
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

    if (m_multiplicitive == 2 && m_levels > 1) {
      // symmetric multiplicative
      // v1 <- C r
      x.setZero();
      vector_type b_copy = b;
#if ABORIA_LOG_LEVEL > 1
      auto t0 = Clock::now();
#endif
      int i = 1;
      auto count = boost::make_counting_iterator(0);
      detail::for_each(
          count, count + m_indicies[i].size(),
          solve_domain<vector_type, Dest>{
              iterator_to_raw_pointer(m_indicies[i].begin()),
              iterator_to_raw_pointer(m_buffer[i].begin()),
              iterator_to_raw_pointer(m_factorized_matrix[i].begin()),
              iterator_to_raw_pointer(m_coupling_matrix[i].begin()), x, b_copy,
              i, m_indicies[i].size() == 1, m_interpolate});

      // v2 <- v1 + M-1 ( r - A v1)
      b_copy = b - m_A * x;
#if ABORIA_LOG_LEVEL > 1
      auto t1 = Clock::now();
      m_timing[i] +=
          std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
              .count();

      t0 = Clock::now();
#endif
      i = 0;
      detail::for_each(
          count, count + m_indicies[i].size(),
          solve_domain<vector_type, Dest>{
              iterator_to_raw_pointer(m_indicies[i].begin()),
              iterator_to_raw_pointer(m_buffer[i].begin()),
              iterator_to_raw_pointer(m_factorized_matrix[i].begin()),
              iterator_to_raw_pointer(m_coupling_matrix[i].begin()), x, b_copy,
              i, m_indicies[i].size() == 1, m_interpolate});

      // v3 <- v2 + C ( r - (r - A v1) - A v2)
      b_copy = b - b_copy - m_A * x;
#if ABORIA_LOG_LEVEL > 1
      t1 = Clock::now();
      m_timing[i] +=
          std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
              .count();

      t0 = Clock::now();
#endif
      i = 1;
      detail::for_each(
          count, count + m_indicies[i].size(),
          solve_domain<vector_type, Dest>{
              iterator_to_raw_pointer(m_indicies[i].begin()),
              iterator_to_raw_pointer(m_buffer[i].begin()),
              iterator_to_raw_pointer(m_factorized_matrix[i].begin()),
              iterator_to_raw_pointer(m_coupling_matrix[i].begin()), x, b_copy,
              i, m_indicies[i].size() == 1, m_interpolate});

#if ABORIA_LOG_LEVEL > 1
      t1 = Clock::now();
      m_timing[i] +=
          std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
              .count();

#endif

    } else if (m_multiplicitive == 1 && m_levels > 1) {
      // non-symmetric multiplicative
      x.setZero();
      vector_type b_copy = b;

      for (int i = 0; i < std::min(m_levels, 2); ++i) {
#if ABORIA_LOG_LEVEL > 1
        auto t0 = Clock::now();
#endif
        solve_domain<vector_type, Dest> fsolve_domain{
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
#if ABORIA_LOG_LEVEL > 1
        auto t1 = Clock::now();
        m_timing[i] +=
            std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
                .count();
#endif
      }
    } else {
      x.setZero();

      for (int i = 0; i < std::min(m_levels, 2); ++i) {
#if ABORIA_LOG_LEVEL > 1
        auto t0 = Clock::now();
#endif

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
#if ABORIA_LOG_LEVEL > 1
        auto t1 = Clock::now();
        m_timing[i] +=
            std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
                .count();
#endif
      }
    }
  }

  template <typename Rhs>
  inline const Eigen::Solve<SchwartzPreconditioner, Rhs>
  solve(const Eigen::MatrixBase<Rhs> &b) const {
    eigen_assert(
        static_cast<typename Rhs::Index>(m_rows) == b.rows() &&
        "SchwartzPreconditioner::solve(): invalid number of rows of the "
        "right hand side matrix b");
    eigen_assert(m_isInitialized &&
                 "SchwartzPreconditioner is not initialized.");
    return Eigen::Solve<SchwartzPreconditioner, Rhs>(*this, b.derived());
  }

  Eigen::ComputationInfo info() { return Eigen::Success; }
};

template <typename Solver> class SchwartzSamplingPreconditioner {
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
  double m_sigma;
  double m_M;

  connectivity_type m_domain_indicies;
  connectivity_type m_domain_buffer;
  std::vector<solver_type> m_domain_factorized_matrix;
  Index m_rows;
  Index m_cols;

public:
  typedef typename vector_type::StorageIndex StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic
  };

  SchwartzSamplingPreconditioner()
      : m_isInitialized(false), m_random(0), m_sigma(-1), m_M(1.0) {}

  template <typename MatType>
  explicit SchwartzSamplingPreconditioner(const MatType &mat) {
    compute(mat);
  }

  Index rows() const { return m_rows; }
  Index cols() const { return m_cols; }

  void set_number_of_random_particles(size_t n) { m_random = n; }
  void set_sigma(double value) { m_sigma = value; }
  void set_rejection_sampling_scale(double value) { m_M = value; }

  template <typename Kernel>
  void analyze_impl_block(const Index start_row, const Kernel &kernel) {
    typedef typename Kernel::row_elements_type row_elements_type;
    typedef typename Kernel::col_elements_type col_elements_type;
    typedef typename row_elements_type::query_type query_type;
    typedef typename query_type::traits_type traits_type;
    typedef typename query_type::child_iterator child_iterator;
    typedef typename traits_type::double_d double_d;
    typedef typename traits_type::position position;

    static_assert(std::is_same<row_elements_type, col_elements_type>::value,
                  "Schwartz preconditioner restricted to identical row and col "
                  "particle sets");
    const row_elements_type &a = kernel.get_row_elements();
    CHECK(&a == &(kernel.get_col_elements()),
          "Schwartz preconditioner restricted to identical row and col "
          "particle "
          "sets");
    std::default_random_engine generator;
    const query_type &query = a.get_query();
    for (auto i = query.get_subtree(); i != false; ++i) {
      if (query.is_leaf_node(*i)) {
        auto ci = i.get_child_iterator();
        auto bounds = query.get_bounds(ci);
        // std::cout << "looking at bucket index " <<
        // query.get_bucket_index(*ci)
        //<< " with bounds " << bounds << std::endl;

        // skip over empty buckets
        if (query.get_bucket_particles(*ci) == false)
          continue;

        const size_t domain_index = m_domain_indicies.size();
        m_domain_indicies.push_back(connectivity_type::value_type());
        m_domain_buffer.push_back(connectivity_type::value_type());
        storage_vector_type &buffer = m_domain_buffer[domain_index];
        storage_vector_type &indicies = m_domain_indicies[domain_index];

        // add particles in bucket to indicies
        for (auto particle = query.get_bucket_particles(*i); particle != false;
             ++particle) {
          const size_t index = &(get<position>(*particle)) -
                               get<position>(query.get_particles_begin());
          indicies.push_back(start_row + index);
        }

        // add buffer particles through random sampling
        // int nspecial = std::pow(3, query_type::dimension);
        int nspecial = 0;
        // const int nspecial = 0;
        buffer.resize(m_random + nspecial);
        std::vector<child_iterator> buckets(m_random + nspecial);
        std::uniform_real_distribution<double> uniform(0, 1);
        std::normal_distribution<double> normal(0, m_sigma);

        // add special points
        // updates nspecial with actual number of special points
        /*
        lattice_iterator<query_type::dimension> special_it(
            int_d::Constant(0), int_d::Constant(3), int_d::Constant(0));
        for (nspecial = 0; special_it != false; ++special_it, ++nspecial) {
          const double_d &bmin = query.get_bounds().bmin;
          const double_d &bmax = query.get_bounds().bmax;
          const double_d pos =
              (*special_it) * (0.5 * (bmax - bmin) -
                               std::numeric_limits<double>::epsilon()) +
              bmin;
          // std::cout <<"adding special point at pos = "<<pos<<std::endl;
          const bool not_in_bucket =
              (pos >= bounds.bmax).any() || (pos < bounds.bmin).any();
          if (not_in_bucket) {
            buckets[nspecial] = query.get_bucket(pos);
          } else {
            --nspecial;
          }
        }
        buffer.resize(m_random + nspecial);
        */

        const double_d middle = 0.5 * (bounds.bmin + bounds.bmax);
        if (m_sigma > 0) {
          const double_d w = bounds.bmax - bounds.bmin;
          if (m_sigma < 0.2 * w.maxCoeff()) {
            std::vector<child_iterator> pot_buckets;
            for (auto it =
                     query.template get_buckets_near_point<-1>(middle, 0.6 * w);
                 it != false; ++it) {
              // fill pot buckets (not self)
              const auto ci = it.get_child_iterator();
              auto pot_bounds = query.get_bounds(ci);
              const bool middle_not_in_pot =
                  (middle >= pot_bounds.bmax).any() ||
                  (middle < pot_bounds.bmin).any();
              if (middle_not_in_pot) {
                pot_buckets.push_back(ci);
              }
            }
            // uniformly sample pot buckets
            std::generate(buckets.begin() + nspecial, buckets.end(), [&]() {
              const int sampled_index =
                  std::floor(uniform(generator) * pot_buckets.size());
              return pot_buckets[sampled_index];
            });
          } else {
            const double scale2 = 1.0 / std::pow(m_sigma, 2);
            auto gaussianf = [&](const double_d &x) {
              return std::exp(-(x - middle).squaredNorm() * scale2);
            };
            std::generate(buckets.begin() + nspecial, buckets.end(), [&]() {
              double_d sp;
              bool accepted;
              do {
                for (size_t i = 0; i < query_type::dimension; i++) {
                  sp[i] = normal(generator) + middle[i];
                }
                const bool not_in_bucket =
                    (sp >= bounds.bmax).any() || (sp < bounds.bmin).any();
                const bool in_domain =
                    (sp < a.get_max()).all() && (sp >= a.get_min()).all();
                accepted = not_in_bucket && in_domain &&
                           uniform(generator) <
                               kernel.get_position_function()(middle, sp) /
                                   (gaussianf(sp) * m_M);
              } while (!accepted);
              return query.get_bucket(sp);
            });
          }
        } else {
          const double volume = (a.get_max() - a.get_min()).prod();
          std::generate(buckets.begin() + nspecial, buckets.end(), [&]() {
            double_d sp;
            bool accepted;
            do {
              for (size_t i = 0; i < query_type::dimension; i++) {
                sp[i] = 0.5 * (a.get_max()[i] - a.get_min()[i]) *
                            uniform(generator) +
                        a.get_min()[i];
              }
              const bool not_in_bucket =
                  (sp >= bounds.bmax).any() || (sp < bounds.bmin).any();
              const bool in_domain =
                  (sp < a.get_max()).all() && (sp >= a.get_min()).all();
              accepted =
                  not_in_bucket && in_domain &&
                  uniform(generator) <
                      kernel.get_position_function()(middle, sp) * volume / m_M;
            } while (!accepted);
            return query.get_bucket(sp);
          });
        }

        std::unordered_map<size_t, std::pair<child_iterator, size_t>> counts;
        for (size_t i = 0; i < buckets.size(); ++i) {
          auto bucket_index = query.get_bucket_index(*(buckets[i]));
          auto it = counts.find(bucket_index);
          if (it != counts.end()) {
            it->second.second++;
          } else {
            counts[bucket_index] = std::make_pair(buckets[i], 1);
          }
        }
        // for (auto i : counts) {
        // std::cout << "bucket index " << i.first << " with bounds "
        //<< query.get_bounds(i.second.first) << " has " << i.second.second
        //<< " counts" << std::endl;
        //}

        int out_index = 0;
        std::for_each(counts.begin(), counts.end(), [&](auto i) {
          auto ci = i.second.first;
          size_t count = i.second.second;
          auto pit = query.get_bucket_particles(*ci);
          auto num_particles = pit.distance_to_end();
          std::vector<int> bucket_indices(num_particles);
          std::iota(bucket_indices.begin(), bucket_indices.end(), 0);
          std::shuffle(bucket_indices.begin(), bucket_indices.end(), generator);
          const int trunc_count = std::min(count, bucket_indices.size());
          std::transform(
              bucket_indices.begin(), bucket_indices.begin() + trunc_count,
              buffer.begin() + out_index, [&](const int i) {
                return (&get<position>(*(pit + i)) - &get<position>(a)[0]) +
                       start_row;
              });
          // std::cout << "looking for " << count
          //<< " samples in buffer. Found at indicies ";
          // for (size_t i = out_index; i < out_index + trunc_count; ++i) {
          //  std::cout << buffer[i] << " ";
          //}
          // std::cout << std::endl;
          out_index += trunc_count;
        });
        buffer.resize(out_index);

#ifdef HAVE_CAIRO_TURN_OFF
        const int image_size = 512;
        cairo_surface_t *surface = cairo_svg_surface_create(
            ("sampler" + std::to_string(domain_index) + ".svg").c_str(),
            image_size, image_size);
        cairo_svg_surface_restrict_to_version(surface, CAIRO_SVG_VERSION_1_2);
        cairo_t *cr = cairo_create(surface);
        const double lw = 0.007;

        cairo_scale(cr, image_size, image_size);
        cairo_set_line_width(cr, lw);
        cairo_set_source_rgba(cr, 0, 0, 0, 0.5);
        cairo_move_to(cr, bounds.bmin[0], bounds.bmin[1]);
        cairo_line_to(cr, bounds.bmax[0], bounds.bmin[1]);
        cairo_line_to(cr, bounds.bmax[0], bounds.bmax[1]);
        cairo_line_to(cr, bounds.bmin[0], bounds.bmax[1]);
        cairo_close_path(cr);
        cairo_stroke(cr);

        const double PI = boost::math::constants::pi<double>();
        cairo_set_source_rgba(cr, 0.5, 0, 0, 0.5);
        for (auto i : indicies) {
          auto &pos = get<position>(a)[i];
          cairo_arc(cr, pos[0], pos[1], lw, 0, 2 * PI);
          cairo_fill(cr);
        }
        cairo_set_source_rgba(cr, 0, 0, 0.5, 0.5);
        for (auto i : buffer) {
          auto &pos = get<position>(a)[i];
          cairo_arc(cr, pos[0], pos[1], lw, 0, 2 * PI);
          cairo_fill(cr);
        }
        cairo_destroy(cr);
        cairo_surface_destroy(surface);
#endif

        // ASSERT(buffer.size() > 0, "no particles in buffer");
        ASSERT(indicies.size() > 0, "no particles in domain");
      }
    }
  } // namespace Aboria

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
  SchwartzSamplingPreconditioner &
  analyzePattern(const MatrixReplacement<NI, NJ, Blocks> &mat) {

    LOG(2, "SchwartzSamplingPreconditioner: analyze pattern");
    m_rows = mat.rows();
    m_cols = mat.cols();
    analyze_impl(mat, detail::make_index_sequence<NI>());

    int count = 0;
    int minsize_buffer = 1000;
    int maxsize_buffer = 0;
    int minsize_indicies = 1000;
    int maxsize_indicies = 0;
    for (size_t domain_index = 0; domain_index < m_domain_indicies.size();
         ++domain_index) {
      const int size_indicies = m_domain_indicies[domain_index].size();
      const int size_buffer = m_domain_buffer[domain_index].size();
      count += size_indicies;
      if (size_buffer < minsize_buffer)
        minsize_buffer = size_buffer;
      if (size_buffer > maxsize_buffer)
        maxsize_buffer = size_buffer;
      if (size_indicies < minsize_indicies)
        minsize_indicies = size_indicies;
      if (size_indicies > maxsize_indicies)
        maxsize_indicies = size_indicies;
    }
    LOG(2, "SchwartzSamplingPreconditioner: finished analysis, found "
               << m_domain_indicies.size() << " domains, with "
               << minsize_indicies << "--" << maxsize_indicies << " particles ("
               << count << " total), and " << minsize_buffer << "--"
               << maxsize_buffer << " buffer particles")
    return *this;
  }

  template <int _Options, typename _StorageIndex>
  SchwartzSamplingPreconditioner &analyzePattern(
      const Eigen::SparseMatrix<Scalar, _Options, _StorageIndex> &mat) {
    CHECK(m_domain_indicies.size() > 0,
          "SchwartzSamplingPreconditioner::analyzePattern(): cannot analyze "
          "sparse "
          "matrix, "
          "call analyzePattern using a Aboria MatrixReplacement class first");
    return *this;
  }

  template <int _Options, typename _StorageIndex, int RefOptions,
            typename RefStrideType>
  SchwartzSamplingPreconditioner &
  analyzePattern(const Eigen::Ref<
                 const Eigen::SparseMatrix<Scalar, _Options, _StorageIndex>,
                 RefOptions, RefStrideType> &mat) {
    CHECK(m_domain_indicies.size() > 0,
          "SchwartzSamplingPreconditioner::analyzePattern(): cannot analyze "
          "sparse "
          "matrix, "
          "call analyzePattern using a Aboria MatrixReplacement class first");
    return *this;
  }

  template <typename Derived>
  SchwartzSamplingPreconditioner &
  analyzePattern(const Eigen::DenseBase<Derived> &mat) {
    CHECK(m_domain_indicies.size() > 0,
          "SchwartzSamplingPreconditioner::analyzePattern(): cannot analyze "
          "dense "
          "matrix, "
          "call analyzePattern need to pass a Aboria MatrixReplacement class "
          "first");
    return *this;
  }

  template <typename MatType>
  SchwartzSamplingPreconditioner &factorize(const MatType &mat) {
    LOG(2, "SchwartzSamplingPreconditioner: factorizing domain");
    eigen_assert(static_cast<typename MatType::Index>(m_rows) == mat.rows() &&
                 "SchwartzSamplingPreconditioner::solve(): invalid number of "
                 "rows of mat");
    eigen_assert(static_cast<typename MatType::Index>(m_cols) == mat.cols() &&
                 "SchwartzSamplingPreconditioner::solve(): invalid number of "
                 "rows of mat");

    m_domain_factorized_matrix.resize(m_domain_indicies.size());

    matrix_type domain_matrix;

    for (size_t domain_index = 0;
         domain_index < m_domain_factorized_matrix.size(); ++domain_index) {
      const storage_vector_type &buffer = m_domain_buffer[domain_index];
      const storage_vector_type &indicies = m_domain_indicies[domain_index];
      solver_type &solver = m_domain_factorized_matrix[domain_index];

      const size_t size = indicies.size() + buffer.size();
      // std::cout << "domain "<<domain_index<<"indicies =
      // "<<indicies.size()<<" buffer =  "<<buffer.size()<<" random =
      // "<<random.size()<<std::endl;

      domain_matrix.resize(size, size);

      size_t i = 0;
      for (const size_t &big_index_i : indicies) {
        size_t j = 0;
        for (const size_t &big_index_j : indicies) {
          domain_matrix(i, j++) = mat.coeff(big_index_i, big_index_j);
        }
        for (const size_t &big_index_j : buffer) {
          domain_matrix(i, j++) = mat.coeff(big_index_i, big_index_j);
        }
        ++i;
      }
      for (const size_t &big_index_i : buffer) {
        size_t j = 0;
        for (const size_t &big_index_j : indicies) {
          domain_matrix(i, j++) = mat.coeff(big_index_i, big_index_j);
        }
        for (const size_t &big_index_j : buffer) {
          domain_matrix(i, j++) = mat.coeff(big_index_i, big_index_j);
        }
        ++i;
      }

      solver.compute(domain_matrix);

      Eigen::VectorXd b = Eigen::VectorXd::Random(domain_matrix.rows());
      Eigen::VectorXd x = solver.solve(b);
      double relative_error = (domain_matrix * x - b).norm() / b.norm();
      if (relative_error > 1e-3 || std::isnan(relative_error)) {
        std::cout << "domain index = " << domain_index
                  << ": relative error = " << relative_error << std::endl;
      }
    }

    m_isInitialized = true;

    return *this;
  }

  template <typename MatType>
  SchwartzSamplingPreconditioner &compute(const MatType &mat) {
    analyzePattern(mat);
    return factorize(mat);
  }

  /** \internal */
  template <typename Rhs, typename Dest>
  void _solve_impl(const Rhs &b, Dest &x) const {
    // loop over domains and invert relevent sub-matricies in
    // mat
    vector_type domain_x;
    vector_type domain_b;
    x = b;
    for (size_t i = 0; i < m_domain_indicies.size(); ++i) {
      if (m_domain_indicies.size() == 0)
        continue;

      const storage_vector_type &buffer = m_domain_buffer[i];
      const storage_vector_type &indicies = m_domain_indicies[i];

      const size_t nb = indicies.size() + buffer.size();
      domain_x.resize(nb);
      domain_b.resize(nb);

      // copy x values from big vector
      size_t sub_index = 0;
      for (size_t j = 0; j < indicies.size(); ++j) {
        domain_b[sub_index++] = b[indicies[j]];
      }
      for (size_t j = 0; j < buffer.size(); ++j) {
        domain_b[sub_index++] = b[buffer[j]];
      }

      // solve domain
      domain_x = m_domain_factorized_matrix[i].solve(domain_b);

      // copy accumulate b values to big vector
      sub_index = 0;
      for (size_t j = 0; j < indicies.size(); ++j) {
        x[indicies[j]] = domain_x[sub_index++];
      }
    }
  }

  template <typename Rhs>
  inline const Eigen::Solve<SchwartzSamplingPreconditioner, Rhs>
  solve(const Eigen::MatrixBase<Rhs> &b) const {
    eigen_assert(static_cast<typename Rhs::Index>(m_rows) == b.rows() &&
                 "SchwartzSamplingPreconditioner::solve(): invalid number of "
                 "rows of the "
                 "right hand side matrix b");
    eigen_assert(m_isInitialized &&
                 "SchwartzSamplingPreconditioner is not initialized.");
    return Eigen::Solve<SchwartzSamplingPreconditioner, Rhs>(*this,
                                                             b.derived());
  }

  Eigen::ComputationInfo info() { return Eigen::Success; }
}; // namespace Aboria

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

template <typename Solver> class NystromSwartzPreconditioner {
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
  NystromPreconditioner<Solver> m_nystrom;
  SchwartzPreconditioner<Solver> m_swartz;

  Index m_rows;
  Index m_cols;

public:
  typedef typename vector_type::StorageIndex StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic
  };

  NystromSwartzPreconditioner() : m_isInitialized(false) {}

  template <typename MatType>
  explicit NystromSwartzPreconditioner(const MatType &mat) {
    compute(mat);
  }

  Index rows() const { return m_rows; }
  Index cols() const { return m_cols; }

  NystromPreconditioner<Solver> &nystrom() { return m_nystrom; }
  SchwartzPreconditioner<Solver> &swartz() { return m_swartz; }

  template <typename MatType>
  NystromSwartzPreconditioner &analyzePattern(const MatType &mat) {
    m_rows = mat.rows();
    m_cols = mat.cols();
    m_nystrom.analyzePattern(mat);
    m_swartz.analyzePattern(mat);
    return *this;
  }

  template <typename MatType>
  NystromSwartzPreconditioner &factorize(const MatType &mat) {
    LOG(2, "NystromSwartzPreconditioner: factorizing domain");
    m_nystrom.factorize(mat);
    m_swartz.factorize(mat);

    m_isInitialized = true;

    return *this;
  }

  template <typename MatType>
  NystromSwartzPreconditioner &compute(const MatType &mat) {
    analyzePattern(mat);
    return factorize(mat);
  }

  /** \internal */
  template <typename Rhs, typename Dest>
  void _solve_impl(const Rhs &b, Dest &x) const {
    m_nystrom._solve_impl(b, x);
    vector_type x_tmp(x.size());
    m_swartz._solve_impl(b, x_tmp);
    x += x_tmp;
  }

  template <typename Rhs>
  inline const Eigen::Solve<NystromSwartzPreconditioner, Rhs>
  solve(const Eigen::MatrixBase<Rhs> &b) const {
    eigen_assert(
        static_cast<typename Rhs::Index>(m_rows) == b.rows() &&
        "NystromSwartzPreconditioner::solve(): invalid number of rows of the "
        "right hand side matrix b");
    eigen_assert(m_isInitialized &&
                 "NystromSwartzPreconditioner is not initialized.");
    return Eigen::Solve<NystromSwartzPreconditioner, Rhs>(*this, b.derived());
  }

  Eigen::ComputationInfo info() { return Eigen::Success; }
};

} // namespace Aboria

#endif // HAVE_EIGEN
#endif
