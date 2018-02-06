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

#ifndef LOWRANK_H
#define LOWRANK_H
#ifdef HAVE_EIGEN

#include <Eigen/Core>

namespace Aboria {
namespace detail {

/*
 * Implementation of partially-pivoted ACA, using reference:
 *
 * Kezhong Zhao, M. N. Vouvakis and Jin-Fa Lee, "The adaptive cross
 * approximation algorithm for accelerated method of moments computations of EMC
 * problems," in IEEE Transactions on Electromagnetic Compatibility, vol. 47,
 * no. 4, pp. 763-773, Nov. 2005. URL:
 * http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1580747&isnumber=33389
 */
template <typename Kernel, typename DerivedU, typename DerivedV>
size_t adaptive_cross_approximation_partial(const Kernel &Z, const size_t max_k,
                                            const double epsilon,
                                            Eigen::MatrixBase<DerivedU> &U,
                                            Eigen::MatrixBase<DerivedV> &V) {
  ASSERT(U.rows() == static_cast<typename DerivedU::Index>(Z.rows()),
         "number of U rows not equal to number of Z rows");
  ASSERT(V.cols() == static_cast<typename DerivedV::Index>(Z.cols()),
         "number of V cols not equal to number of Z cols");
  typedef Kernel matrix_type;
  typedef Eigen::Matrix<double, matrix_type::RowsAtCompileTime, 1>
      col_vector_type;
  typedef Eigen::Matrix<double, 1, matrix_type::ColsAtCompileTime>
      row_vector_type;
  typedef size_t index_type;
  // init
  // init I_1 = 1 and set Z = 0
  const index_type rows = Z.rows();
  const index_type cols = Z.cols();
  const double epsilon2 = std::pow(epsilon, 2);
  std::list<index_type> row_pivots(rows - 1);
  std::list<index_type> col_pivots(cols);
  typedef typename std::list<index_type>::iterator pivot_iterator;

  auto rowi = row_pivots.begin();
  auto coli = col_pivots.begin();
  for (size_t i = 0; i < cols; ++i) {
    *coli = i;
    ++coli;
  }
  for (size_t i = 1; i < rows; ++i) {
    *rowi = i;
    ++rowi;
  }

  LOG(3, "adaptive_cross_approximation (parial) (Z = ["
             << rows << "," << cols << "] max_k = " << max_k
             << " epsilon = " << epsilon << ")")

  index_type k = 0;
  index_type i = 0;
  index_type j = 0;
  double uv_norm2 = 1;
  double Z_norm2 = 0;
  while (uv_norm2 > epsilon2 * Z_norm2 && k < max_k) {
    row_vector_type Rrow(cols);
    for (size_t jj = 0; jj < cols; ++jj) {
      Rrow(jj) = Z.coeff(i, jj);
    }
    for (size_t l = 0; l < k; ++l) {
      Rrow -= U.col(l)(i) * V.row(l);
    }

    // find next col index
    double max = 0;
    pivot_iterator it_max;
    for (auto it = col_pivots.begin(); it != col_pivots.end(); ++it) {
      const double aij = std::abs(Rrow(*it));
      if (max < aij) {
        max = aij;
        j = *it;
        it_max = it;
      }
    }
    col_pivots.erase(it_max);

    LOG(4, "\tcol pivot = " << j);

    V.row(k) = Rrow / Rrow(j);

    col_vector_type Rcol(rows);
    for (size_t ii = 0; ii < rows; ++ii) {
      Rcol(ii) = Z.coeff(ii, j);
    }
    for (size_t l = 0; l < k; ++l) {
      Rcol -= V.row(l)(j) * U.col(l);
    }

    U.col(k) = Rcol;
    uv_norm2 = V.row(k).norm() * U.col(k).norm();
    double sum = 0;
    for (size_t l = 0; l < k; ++l) {
      const double innerU = U.col(l).dot(U.col(k));
      const double innerV = V.row(l).dot(V.row(k));
      sum += std::abs(innerU) * std::abs(innerV);
    }
    Z_norm2 += 2 * sum + uv_norm2;
    LOG(4, "\tuv_norm2 = " << uv_norm2 << " Z_norm2 = " << Z_norm2);

    // increment k
    k++;

    // find row index for next iteration
    if (k < max_k) {
      max = 0;
      for (auto it = row_pivots.begin(); it != row_pivots.end(); ++it) {
        const double aij = std::abs(Rcol(*it));
        if (max < aij) {
          max = aij;
          i = *it;
          it_max = it;
        }
      }
      row_pivots.erase(it_max);
      LOG(4, "\trow pivot = " << i);
    }
  }
  return k;
}

/*
 * Implementation of fully-pivoted ACA (same as fully-pivoted LU decomp),
 * using reference:
 *
 * Note: overwrites input matrix Z
 *
 * Bebendorf, Mario, and Sergej Rjasanow. "Adaptive low-rank approximation of
 * collocation matrices." Computing 70.1 (2003): 1-24.
 */
template <typename DerivedZ, typename DerivedU, typename DerivedV>
size_t adaptive_cross_approximation_full(Eigen::MatrixBase<DerivedZ> &Z,
                                         const size_t max_k,
                                         const double epsilon,
                                         Eigen::MatrixBase<DerivedU> &U,
                                         Eigen::MatrixBase<DerivedV> &V) {
  ASSERT(U.rows() == Z.rows(),
         "number of U rows not equal to number of Z rows");
  ASSERT(V.cols() == Z.cols(),
         "number of V cols not equal to number of Z cols");
  typedef size_t index_type;
  // init
  // init I_1 = 1 and set Z = 0
  const index_type rows = Z.rows();
  const index_type cols = Z.cols();
  const double epsilon2 = std::pow(epsilon, 2);

  LOG(3, "adaptive_cross_approximation (full) (Z = ["
             << rows << "," << cols << "] max_k = " << max_k
             << " epsilon = " << epsilon << ")")

  index_type k = 0;
  index_type i = 0;
  index_type j = 0;
  double B_norm2 = Z.squaredNorm();
  double Z_norm2 = epsilon2 * B_norm2 + 1.0;
  while (Z_norm2 > epsilon2 * B_norm2 && k < max_k) {
    // find max element i,j
    double max = 0;
    for (size_t ii = 0; ii < rows; ++ii) {
      for (size_t jj = 0; jj < cols; ++jj) {
        const double aij = std::abs(Z(ii, jj));
        if (max < aij) {
          max = aij;
          i = ii;
          j = jj;
        }
      }
    }

    LOG(4, "\trow pivot = " << i);
    LOG(4, "\tcol pivot = " << j);

    // store u and v vectors
    //
    V.row(k) = Z.row(i) / Z(i, j);
    U.col(k) = Z.col(j);

    // update Z
    Z -= U.col(k) * V.row(k);

    // stopping criteria
    Z_norm2 = Z.squaredNorm();

    // increment k
    ++k;
  }
  return k;
}

} // namespace detail
} // namespace Aboria

#endif
#endif /* ifndef LOWRANK_H */
