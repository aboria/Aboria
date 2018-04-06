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

#ifndef FAST_MULTIPOLE_METHOD_DETAIL_H_
#define FAST_MULTIPOLE_METHOD_DETAIL_H_

#include "../Get.h"
#include "CudaInclude.h"
#include "Log.h"
#include "NeighbourSearchBase.h"
#include "SpatialUtil.h"
#include "Traits.h"
#include "Vector.h"
#include "detail/Chebyshev.h"
#include "detail/Kernels.h"
#include <iostream>

#ifdef HAVE_H2LIB
extern "C" {
#include <amatrix.h>
#include <cluster.h>
#undef I
}
#endif

namespace Aboria {
namespace detail {

template <unsigned int D, unsigned int N> struct MultiquadricExpansions {
  typedef bbox<D> box_type;
  static const size_t ncheb = std::pow(N, D);
  typedef std::array<double, ncheb> expansion_type;
#ifdef HAVE_EIGEN
  typedef Eigen::Matrix<double, ncheb, ncheb> matrix_type;
  typedef Eigen::Matrix<double, ncheb, Eigen::Dynamic> p2m_matrix_type;
  typedef Eigen::Matrix<double, Eigen::Dynamic, ncheb> m2p_matrix_type;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> dynamic_vector_type;
#endif
  typedef Vector<double, D> double_d;
  typedef Vector<int, D> int_d;
  static const unsigned int dimension = D;
  const double m_c2;

  MultiquadricExpansions(const double c) : m_c2(c * c) {}

  static void P2M(expansion_type &accum, const box_type &box,
                  const double_d &position, const double &source) {}

  static void M2M(expansion_type &accum, const box_type &target_box,
                  const box_type &source_box, const expansion_type &source) {}

  void M2L(expansion_type &accum, const box_type &target_box,
           const box_type &source_box, const expansion_type &source) {}

  static void L2L(expansion_type &accum, const box_type &target_box,
                  const box_type &source_box, const expansion_type &source) {}

  static double L2P(const double_d &p, const box_type &box,
                    const expansion_type &source) {
    return 0.0;
  }
};

template <size_t D, size_t N, typename Function, size_t BlockRows,
          size_t BlockCols>
struct BlackBoxExpansions {
  typedef bbox<D> box_type;
  static constexpr size_t ncheb = ipow(N, D);

  static constexpr size_t block_rows = BlockRows;
  static constexpr size_t block_cols = BlockCols;

  typedef typename std::conditional<
      BlockRows == 1, double, Vector<double, BlockRows>>::type row_scalar_type;
  typedef typename std::conditional<
      BlockCols == 1, double, Vector<double, BlockCols>>::type col_scalar_type;
  typedef std::array<col_scalar_type, ncheb> m_expansion_type;
  typedef std::array<row_scalar_type, ncheb> l_expansion_type;

#ifdef HAVE_EIGEN
  typedef Eigen::Matrix<double, ncheb * BlockRows, ncheb * BlockCols>
      m2l_matrix_type;
  typedef Eigen::Matrix<double, ncheb * BlockRows, ncheb * BlockRows>
      l2l_matrix_type;
  typedef Eigen::Matrix<double, ncheb * BlockCols, ncheb * BlockCols>
      m2m_matrix_type;
  typedef Eigen::Matrix<double, ncheb * BlockCols, Eigen::Dynamic>
      p2m_matrix_type;
  typedef Eigen::Matrix<double, Eigen::Dynamic, ncheb * BlockRows>
      l2p_matrix_type;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> p2p_matrix_type;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> p_vector_type;
  typedef Eigen::Matrix<double, ncheb * BlockCols, 1> m_vector_type;
  typedef Eigen::Matrix<double, ncheb * BlockRows, 1> l_vector_type;
  typedef Eigen::Matrix<double, BlockCols, BlockCols> blockcol_type;
  typedef Eigen::Matrix<double, BlockRows, BlockRows> blockrow_type;
#endif

  typedef Vector<double, D> double_d;
  typedef Vector<int, D> int_d;
  static const unsigned int dimension = D;
  static const unsigned int number_of_nodes_in_each_direction = N;
  Function m_K;
  std::array<double_d, ncheb> m_cheb_points;

  BlackBoxExpansions(const Function &K) : m_K(K) {
    // precalculate cheb_points
    lattice_iterator<dimension> mi(int_d::Constant(0), int_d::Constant(N));
    for (size_t i = 0; i < ncheb; ++i, ++mi) {
      m_cheb_points[i] = detail::chebyshev_node_nd(*mi, N);
    }
  }

  static void P2M(m_expansion_type &accum, const box_type &box,
                  const double_d &position, const col_scalar_type &source) {

    detail::ChebyshevRnSingle<D, N> cheb_rn(position, box);
    lattice_iterator<dimension> mj(int_d::Constant(0), int_d::Constant(N));
    for (size_t j = 0; j < ncheb; ++j, ++mj) {
      // std::cout << "accumulating P2M from "<<position<<" to node "<<*mj<<"
      // with Rn = "<<cheb_rn(*mj)<<std::endl;
      accum[j] += cheb_rn(*mj) * source;
    }
  }

#ifdef HAVE_EIGEN

  /*
  template <typename Derived>
  static void P2M(m_expansion_type& accum,
           const box_type& box,
           const double_d& position,
           const Eigen::DenseBase<Derived>& source) {

      detail::ChebyshevRnSingle<D,N> cheb_rn(position,box);
      lattice_iterator<dimension> mj(int_d::Constant(0),int_d::Constant(N));
      for (int j=0; j<ncheb; ++j,++mj) {
          //std::cout << "accumulating P2M from "<<position<<" to node "<<*mj<<"
  with Rn = "<<cheb_rn(*mj)<<std::endl; accum[j] += cheb_rn(*mj)*source;
      }
  }
  */

  template <typename ParticlesType>
  static void P2M_matrix(p2m_matrix_type &matrix, const box_type &box,
                         const std::vector<size_t> &indicies,
                         const ParticlesType &particles) {
    typedef typename ParticlesType::position position;
    matrix.resize(ncheb * BlockCols, indicies.size() * BlockCols);
    for (size_t i = 0; i < indicies.size(); ++i) {
      const double_d &p = get<position>(particles)[indicies[i]];
      detail::ChebyshevRnSingle<D, N> cheb_rn(p, box);
      lattice_iterator<dimension> mj(int_d::Constant(0), int_d::Constant(N));
      for (size_t j = 0; j < ncheb; ++j, ++mj) {
        // check ij
        matrix.template block<BlockCols, BlockCols>(i * BlockCols,
                                                    j * BlockCols) =
            cheb_rn(*mj) * blockcol_type::Identity();
      }
    }
  }
#endif

  void M2M(m_expansion_type &accum, const box_type &target_box,
           const box_type &source_box, const m_expansion_type &source) const {

    for (size_t j = 0; j < ncheb; ++j) {
      const double_d &pj_unit_box = m_cheb_points[j];
      const double_d pj =
          0.5 * (pj_unit_box + 1) * (source_box.bmax - source_box.bmin) +
          source_box.bmin;
      detail::ChebyshevRnSingle<D, N> cheb_rn(pj, target_box);

      lattice_iterator<D> mi(int_d::Constant(0), int_d::Constant(N));
      for (size_t i = 0; i < ncheb; ++i, ++mi) {
        accum[i] += cheb_rn(*mi) * source[j];
      }
    }
  }

#ifdef HAVE_EIGEN
  void M2M_matrix(m2m_matrix_type &matrix, const box_type &target_box,
                  const box_type &source_box) const {
    for (size_t j = 0; j < ncheb; ++j) {
      const double_d &pj_unit_box = m_cheb_points[j];
      const double_d pj =
          0.5 * (pj_unit_box + 1) * (source_box.bmax - source_box.bmin) +
          source_box.bmin;
      detail::ChebyshevRnSingle<D, N> cheb_rn(pj, target_box);

      lattice_iterator<D> mi(int_d::Constant(0), int_d::Constant(N));
      for (size_t i = 0; i < ncheb; ++i, ++mi) {
        matrix.template block<BlockCols, BlockCols>(i * BlockCols,
                                                    j * BlockCols) =
            cheb_rn(*mi) * blockcol_type::Identity();
      }
    }
  }
#endif

  void M2L(l_expansion_type &accum, const box_type &target_box,
           const box_type &source_box, const m_expansion_type &source) const {

    for (size_t i = 0; i < ncheb; ++i) {
      const double_d &pi_unit_box = m_cheb_points[i];
      const double_d pi =
          0.5 * (pi_unit_box + 1) * (target_box.bmax - target_box.bmin) +
          target_box.bmin;

      for (size_t j = 0; j < ncheb; ++j) {
        const double_d &pj_unit_box = m_cheb_points[j];
        const double_d pj =
            0.5 * (pj_unit_box + 1) * (source_box.bmax - source_box.bmin) +
            source_box.bmin;
        accum[i] += m_K(pi, pj) * source[j];
      }
    }
  }

#ifdef HAVE_EIGEN
  void M2L_matrix(m2l_matrix_type &matrix, const box_type &target_box,
                  const box_type &source_box) const {
    for (size_t i = 0; i < ncheb; ++i) {
      const double_d &pi_unit_box = m_cheb_points[i];
      const double_d pi =
          0.5 * (pi_unit_box + 1) * (target_box.bmax - target_box.bmin) +
          target_box.bmin;
      for (size_t j = 0; j < ncheb; ++j) {
        const double_d &pj_unit_box = m_cheb_points[j];
        const double_d pj =
            0.5 * (pj_unit_box + 1) * (source_box.bmax - source_box.bmin) +
            source_box.bmin;
        matrix.template block<BlockRows, BlockCols>(
            i * BlockRows, j * BlockCols) = m_K(pi, pj);
      }
    }
  }
#endif

  void L2L(l_expansion_type &accum, const box_type &target_box,
           const box_type &source_box, const l_expansion_type &source) const {
    // M2M(accum,target_box,source_box,source);
    for (size_t i = 0; i < ncheb; ++i) {
      const double_d &pi_unit_box = m_cheb_points[i];
      const double_d pi =
          0.5 * (pi_unit_box + 1) * (target_box.bmax - target_box.bmin) +
          target_box.bmin;
      detail::ChebyshevRnSingle<D, N> cheb_rn(pi, source_box);

      lattice_iterator<D> mj(int_d::Constant(0), int_d::Constant(N));
      for (size_t j = 0; j < ncheb; ++j, ++mj) {
        accum[i] += cheb_rn(*mj) * source[j];
      }
    }
  }

#ifdef HAVE_EIGEN
  void L2L_matrix(l2l_matrix_type &matrix, const box_type &target_box,
                  const box_type &source_box) const {
    for (size_t i = 0; i < ncheb; ++i) {
      const double_d &pi_unit_box = m_cheb_points[i];
      const double_d pi =
          0.5 * (pi_unit_box + 1) * (target_box.bmax - target_box.bmin) +
          target_box.bmin;
      detail::ChebyshevRnSingle<D, N> cheb_rn(pi, source_box);

      lattice_iterator<D> mj(int_d::Constant(0), int_d::Constant(N));
      for (size_t j = 0; j < ncheb; ++j, ++mj) {
        matrix.template block<BlockRows, BlockRows>(i * BlockRows,
                                                    j * BlockRows) =
            cheb_rn(*mj) * blockrow_type::Identity();
      }
    }
  }
#endif

  static row_scalar_type L2P(const double_d &p, const box_type &box,
                             const l_expansion_type &source) {
    row_scalar_type sum = detail::VectorTraits<row_scalar_type>::Zero();
    detail::ChebyshevRnSingle<D, N> cheb_rn(p, box);
    lattice_iterator<dimension> mj(int_d::Constant(0), int_d::Constant(N));
    for (size_t j = 0; j < ncheb; ++j, ++mj) {
      sum += cheb_rn(*mj) * source[j];
    }
    return sum;
  }

#ifdef HAVE_EIGEN

  /*
  template <typename Derived>
  static void L2P(const double_d& p,
             const box_type& box,
             const l_expansion_type& source,
             const Eigen::DenseBase<Derived>& sum
             ) {
      detail::ChebyshevRnSingle<D,N> cheb_rn(p,box);
      lattice_iterator<dimension> mj(int_d::Constant(0),int_d::Constant(N));
      for (int j=0; j<ncheb; ++j,++mj) {
          const_cast<Eigen::DenseBase<Derived>&>(sum) += cheb_rn(*mj)*source[j];
      }
  }

  static row_scalar_type L2P(const double_d& p,
             const box_type& box,
             const l_vector_type& source) {
      detail::ChebyshevRnSingle<D,N> cheb_rn(p,box);
      lattice_iterator<dimension> mj(int_d::Constant(0),int_d::Constant(N));
      row_scalar_type sum = detail::VectorTraits<row_scalar_type>::Zero();
      for (int j=0; j<ncheb; ++j,++mj) {
          sum += cheb_rn(*mj)*source.segment<BlockRows>(j*BlockRows);
      }
      return sum;
  }
  */

  template <typename ParticlesType>
  static void L2P_matrix(l2p_matrix_type &matrix, const box_type &box,
                         const std::vector<size_t> &indicies,
                         const ParticlesType &particles) {
    typedef typename ParticlesType::position position;
    matrix.resize(indicies.size() * BlockRows, ncheb * BlockRows);
    for (size_t i = 0; i < indicies.size(); ++i) {
      const double_d &p = get<position>(particles)[indicies[i]];
      detail::ChebyshevRnSingle<D, N> cheb_rn(p, box);
      lattice_iterator<dimension> mj(int_d::Constant(0), int_d::Constant(N));
      for (size_t j = 0; j < ncheb; ++j, ++mj) {
        matrix.template block<BlockRows, BlockRows>(i * BlockRows,
                                                    j * BlockRows) =
            cheb_rn(*mj) * blockrow_type::Identity();
      }
    }
  }

#endif
};

#ifdef HAVE_H2LIB
template <size_t D, typename Function, size_t BlockRows, size_t BlockCols>
struct H2LibBlackBoxExpansions {
  typedef bbox<D> box_type;

  static constexpr size_t block_rows = BlockRows;
  static constexpr size_t block_cols = BlockCols;

  typedef Vector<double, BlockRows> row_scalar_type;
  typedef Vector<double, BlockCols> col_scalar_type;
#ifndef HAVE_EIGEN
  static_assert(BlockRows == BlockCols == 1,
                "Need to define HAVE_EIGEN to use matrix-valued kernel");
#endif

  typedef Vector<double, D> double_d;
  typedef Vector<int, D> int_d;
  static const unsigned int dimension = D;
  Function m_K;
  std::vector<std::vector<double_d>> m_cheb_points;
  size_t m_order;
  size_t m_beta;
  size_t m_max_tree_depth;

  H2LibBlackBoxExpansions(const size_t order, const Function &K,
                          const size_t beta = 1,
                          const size_t max_tree_depth = 30)
      : m_K(K), m_order(order), m_beta(beta),
        m_max_tree_depth(beta == 0 ? 1 : max_tree_depth) {

    // precalculate cheb_points
    m_cheb_points.resize(m_max_tree_depth);
    size_t curr_order = m_order;
    for (auto &cheb_points : m_cheb_points) {
      cheb_points.resize(std::pow(curr_order, D));
      lattice_iterator<dimension> mi(int_d::Constant(0),
                                     int_d::Constant(curr_order));
      for (size_t i = 0; i < cheb_points.size(); ++i, ++mi) {
        cheb_points[i] = detail::chebyshev_node_nd(*mi, curr_order);
      }
      curr_order += m_beta;
    }
  }

  const std::vector<double_d> &get_cheb_points(const size_t order) const {
    const int level = (m_beta == 0) ? 0 : (order - m_order) / m_beta;
    ASSERT_CUDA(level >= 0);
    ASSERT_CUDA(level < static_cast<int>(m_cheb_points.size()));
    ASSERT_CUDA(m_cheb_points[level].size() == std::pow(order, D));
    return m_cheb_points[level];
  }

  template <typename ParticlesType>
  void P2M_trans_amatrix(pamatrix matrix, const pccluster t,
                         const uint *indicies, const uint indicies_size,
                         const ParticlesType &particles) const {
    typedef typename ParticlesType::position position;
    box_type box;
    for (size_t i = 0; i < D; ++i) {
      box.bmin[i] = t->bmin[i];
      box.bmax[i] = t->bmax[i];
    }
    const size_t ncheb = std::pow(m_order, D);
    ASSERT_CUDA(matrix->rows == indicies_size);
    ASSERT_CUDA(matrix->cols == ncheb * BlockCols);
    // resize_amatrix(matrix,ncheb,indicies_size);
    clear_amatrix(matrix);
    detail::ChebyshevRn<D> cheb_rn(m_order, box);
    for (size_t i = 0; i < indicies_size; i += BlockCols) {
      const double_d &p = get<position>(particles)[indicies[i] / BlockCols];
      cheb_rn.set_position(p);
      lattice_iterator<dimension> mj(int_d::Constant(0),
                                     int_d::Constant(m_order));
      for (size_t j = 0; j < ncheb; ++j, ++mj) {
        const double tmp = cheb_rn(*mj);
        for (size_t ii = 0; ii < BlockCols; ++ii) {
          setentry_amatrix(matrix, i + ii, j * BlockCols + ii, tmp);
        }
      }
    }
    /*
    std::cout << "P2M_trans = " << std::endl;
    for (size_t i = 0; i < indicies_size; ++i) {
        std::cout << "| ";
        for (size_t j = 0; j < m_ncheb * BlockCols; ++j) {
            std::cout << getentry_amatrix(matrix, i, j) << " ";
        }
        std::cout << "|" << std::endl;
    }
    */
  }

  void M2M_amatrix(pamatrix matrix, pccluster target_t, pccluster source_t,
                   const size_t target_order, const size_t source_order) const {
    // resize_amatrix(matrix,m_ncheb,m_ncheb);
    box_type target_box, source_box;
    for (size_t i = 0; i < D; ++i) {
      target_box.bmin[i] = target_t->bmin[i];
      target_box.bmax[i] = target_t->bmax[i];
      source_box.bmin[i] = source_t->bmin[i];
      source_box.bmax[i] = source_t->bmax[i];
    }
    detail::ChebyshevRn<D> cheb_rn(target_order, target_box);
    const auto &cheb_points = get_cheb_points(source_order);
    const size_t target_ncheb = std::pow(target_order, D);
    const size_t source_ncheb = std::pow(source_order, D);
    ASSERT_CUDA(matrix->rows == target_ncheb * BlockCols);
    ASSERT_CUDA(matrix->cols == source_ncheb * BlockCols);
    clear_amatrix(matrix);
    for (size_t j = 0; j < source_ncheb; ++j) {
      const double_d &pj_unit_box = cheb_points[j];
      const double_d pj =
          0.5 * (pj_unit_box + 1) * (source_box.bmax - source_box.bmin) +
          source_box.bmin;
      cheb_rn.set_position(pj);

      lattice_iterator<D> mi(int_d::Constant(0), int_d::Constant(target_order));
      for (size_t i = 0; i < target_ncheb; ++i, ++mi) {
        const double tmp = cheb_rn(*mi);
        for (size_t ii = 0; ii < BlockCols; ++ii) {
          setentry_amatrix(matrix, i * BlockCols + ii, j * BlockCols + ii, tmp);
        }
      }
    }
  }

  void M2M_trans_amatrix(pamatrix matrix, pccluster target_t,
                         pccluster source_t, const size_t target_order,
                         const size_t source_order) const {

    // resize_amatrix(matrix,m_ncheb,m_ncheb);
    box_type target_box, source_box;
    for (size_t i = 0; i < D; ++i) {
      target_box.bmin[i] = target_t->bmin[i];
      target_box.bmax[i] = target_t->bmax[i];
      source_box.bmin[i] = source_t->bmin[i];
      source_box.bmax[i] = source_t->bmax[i];
    }
    detail::ChebyshevRn<D> cheb_rn(target_order, target_box);

    const auto &cheb_points = get_cheb_points(source_order);
    const size_t target_ncheb = std::pow(target_order, D);
    const size_t source_ncheb = std::pow(source_order, D);
    ASSERT_CUDA(matrix->cols == target_ncheb * BlockCols);
    ASSERT_CUDA(matrix->rows == source_ncheb * BlockCols);
    clear_amatrix(matrix);
    for (size_t j = 0; j < source_ncheb; ++j) {
      const double_d &pj_unit_box = cheb_points[j];
      const double_d pj =
          0.5 * (pj_unit_box + 1) * (source_box.bmax - source_box.bmin) +
          source_box.bmin;
      cheb_rn.set_position(pj);

      lattice_iterator<D> mi(int_d::Constant(0), int_d::Constant(target_order));
      for (size_t i = 0; i < target_ncheb; ++i, ++mi) {
        const double tmp = cheb_rn(*mi);
        for (size_t ii = 0; ii < BlockCols; ++ii) {
          setentry_amatrix(matrix, j * BlockCols + ii, i * BlockCols + ii, tmp);
        }
      }
    }
    /*
    std::cout << "M2M_trans = " << std::endl;
    for (size_t i = 0; i < m_ncheb * BlockCols; ++i) {
        std::cout << "| ";
        for (size_t j = 0; j < m_ncheb * BlockCols; ++j) {
            std::cout << getentry_amatrix(matrix, i, j) << " ";
        }
        std::cout << "|" << std::endl;
    }
    */
  }

  void M2L_amatrix(pamatrix matrix, pccluster target_t, pccluster source_t,
                   const size_t target_order, const size_t source_order) const {
    // don't resize, already done in new_uniform
    // resize_amatrix(matrix,m_ncheb,m_ncheb);
    box_type target_box, source_box;
    for (size_t i = 0; i < D; ++i) {
      target_box.bmin[i] = target_t->bmin[i];
      target_box.bmax[i] = target_t->bmax[i];
      source_box.bmin[i] = source_t->bmin[i];
      source_box.bmax[i] = source_t->bmax[i];
    }
    const auto &source_cheb_points = get_cheb_points(source_order);
    const auto &target_cheb_points = get_cheb_points(target_order);
    const size_t target_ncheb = std::pow(target_order, D);
    const size_t source_ncheb = std::pow(source_order, D);
    ASSERT_CUDA(matrix->rows == target_ncheb * BlockRows);
    ASSERT_CUDA(matrix->cols == source_ncheb * BlockCols);
    clear_amatrix(matrix);
    for (size_t i = 0; i < target_ncheb; ++i) {
      const double_d &pi_unit_box = target_cheb_points[i];
      const double_d pi =
          0.5 * (pi_unit_box + 1) * (target_box.bmax - target_box.bmin) +
          target_box.bmin;
      for (size_t j = 0; j < source_ncheb; ++j) {
        const double_d &pj_unit_box = source_cheb_points[j];
        const double_d pj =
            0.5 * (pj_unit_box + 1) * (source_box.bmax - source_box.bmin) +
            source_box.bmin;
        // std::cout << "pi = "<<pi<<" pj = "<<pj<< "m_k =
        // "<<m_K(pi,pj)<<std::endl;
#ifdef HAVE_EIGEN
        const Eigen::Matrix<double, BlockRows, BlockCols> tmp(m_K(pi, pj));
        // std::cout << "M2L ij = "<<tmp<<" versus "<<m_K(pi,pj) << std::endl;
        for (size_t ii = 0; ii < BlockRows; ++ii) {
          for (size_t jj = 0; jj < BlockCols; ++jj) {
            setentry_amatrix(matrix, i * BlockRows + ii, j * BlockCols + jj,
                             tmp(ii, jj));
          }
        }
#else
        setentry_amatrix(matrix, i, j, m_K(pi, pj));
#endif
      }
    }
    /*
    std::cout << "M2L = " << std::endl;
    for (int i=0; i<m_ncheb*BlockRows; ++i) {
        std::cout << "| ";
        for (int j=0; j<m_ncheb*BlockCols; ++j) {
            std::cout << getentry_amatrix(matrix,i,j) << " ";
        }
        std::cout << "|" <<std::endl;
    }
    */
  }

  void L2L_amatrix(pamatrix matrix, pccluster target_t, pccluster source_t,
                   const size_t target_order, const size_t source_order) const {
    // resize_amatrix(matrix,m_ncheb,m_ncheb);
    box_type target_box, source_box;
    for (size_t i = 0; i < D; ++i) {
      target_box.bmin[i] = target_t->bmin[i];
      target_box.bmax[i] = target_t->bmax[i];
      source_box.bmin[i] = source_t->bmin[i];
      source_box.bmax[i] = source_t->bmax[i];
    }
    detail::ChebyshevRn<D> cheb_rn(source_order, source_box);
    const auto &cheb_points = get_cheb_points(target_order);
    const size_t target_ncheb = std::pow(target_order, D);
    const size_t source_ncheb = std::pow(source_order, D);
    ASSERT_CUDA(matrix->rows == target_ncheb * BlockRows);
    ASSERT_CUDA(matrix->cols == source_ncheb * BlockRows);
    clear_amatrix(matrix);

    // std:cout << "L2L: D ="<<D<<" order = "<<m_order<<" target_box =
    // "<<target_box<<" source_box = "<<source_box<<std::endl;
    for (size_t i = 0; i < target_ncheb; ++i) {
      const double_d &pi_unit_box = cheb_points[i];
      const double_d pi =
          0.5 * (pi_unit_box + 1) * (target_box.bmax - target_box.bmin) +
          target_box.bmin;
      cheb_rn.set_position(pi);
      lattice_iterator<D> mj(int_d::Constant(0), int_d::Constant(source_order));
      for (size_t j = 0; j < source_ncheb; ++j, ++mj) {
        const double tmp = cheb_rn(*mj);
        // std::cout << "mj = "<<*mj<<" tmp  = "<<tmp<<std::endl;
        for (size_t ii = 0; ii < BlockRows; ++ii) {
          setentry_amatrix(matrix, i * BlockRows + ii, j * BlockRows + ii, tmp);
        }
      }
    }
    /*
    std::cout << "L2L = " << std::endl;
    for (int i=0; i<m_ncheb*BlockRows; ++i) {
        std::cout << "| ";
        for (int j=0; j<m_ncheb*BlockRows; ++j) {
            std::cout << getentry_amatrix(matrix,i,j) << " ";
        }
        std::cout << "|" <<std::endl;
    }
    */
  }

  template <typename ParticlesType>
  void L2P_amatrix(pamatrix matrix, const pccluster t, const uint *indicies,
                   const uint indicies_size,
                   const ParticlesType &particles) const {
    typedef typename ParticlesType::position position;
    // resize_amatrix(matrix,indicies_size,m_ncheb);
    box_type box;
    for (size_t i = 0; i < D; ++i) {
      box.bmin[i] = t->bmin[i];
      box.bmax[i] = t->bmax[i];
    }
    const size_t ncheb = std::pow(m_order, D);
    ASSERT_CUDA(matrix->rows == indicies_size);
    ASSERT_CUDA(matrix->cols == ncheb * BlockRows);
    clear_amatrix(matrix);
    // std:cout << "L2P: D ="<<D<<" order = "<<m_order<<" box =
    // "<<box<<std::endl;
    detail::ChebyshevRn<D> cheb_rn(m_order, box);
    for (size_t i = 0; i < indicies_size; i += BlockRows) {
      const double_d &p = get<position>(particles)[indicies[i] / BlockRows];
      // std::cout << "p = "<<p<<std::endl;
      cheb_rn.set_position(p);
      lattice_iterator<dimension> mj(int_d::Constant(0),
                                     int_d::Constant(m_order));
      for (size_t j = 0; j < ncheb; ++j, ++mj) {
        const double tmp = cheb_rn(*mj);
        for (size_t ii = 0; ii < BlockRows; ++ii) {
          setentry_amatrix(matrix, i + ii, j * BlockRows + ii, tmp);
        }
      }
    }
    /*
    std::cout << "L2P = " << std::endl;
    for (int i=0; i<indicies_size; ++i) {
        std::cout << "| ";
        for (int j=0; j<m_ncheb*BlockRows; ++j) {
            std::cout << getentry_amatrix(matrix,i,j) << " ";
        }
        std::cout << "|" <<std::endl;
    }
    */
  }

  template <typename ParticlesType>
  void P2M_amatrix(pamatrix matrix, const pccluster t, const uint *indicies,
                   const uint indicies_size,
                   const ParticlesType &particles) const {
    typedef typename ParticlesType::position position;
    box_type box;
    for (size_t i = 0; i < D; ++i) {
      box.bmin[i] = t->bmin[i];
      box.bmax[i] = t->bmax[i];
    }
    const size_t ncheb = std::pow(m_order, D);
    ASSERT_CUDA(matrix->rows == ncheb * BlockCols);
    ASSERT_CUDA(matrix->cols == indicies_size);
    clear_amatrix(matrix);
    detail::ChebyshevRn<D> cheb_rn(m_order, box);
    for (size_t i = 0; i < indicies_size; i += BlockCols) {
      const double_d &p = get<position>(particles)[indicies[i] / BlockCols];
      cheb_rn.set_position(p);
      lattice_iterator<dimension> mj(int_d::Constant(0),
                                     int_d::Constant(m_order));
      for (size_t j = 0; j < ncheb; ++j, ++mj) {
        const double tmp = cheb_rn(*mj);
        for (size_t ii = 0; ii < BlockCols; ++ii) {
          setentry_amatrix(matrix, j * BlockCols + ii, i + ii, tmp);
        }
      }
    }
  }
};
#endif

template <typename Expansions, typename Traits, typename T,
          typename VectorType = typename Expansions::expansion_type,
          typename SourceParticleIterator = typename Traits::raw_pointer,
          unsigned int D = Traits::dimension>
void calculate_P2M(VectorType &sum, const bbox<D> &box,
                   const ranges_iterator<Traits> &range,
                   const std::vector<T> &source_vector,
                   const SourceParticleIterator &source_particles_begin,
                   const Expansions &expansions) {
  typedef typename Traits::position position;
  const size_t N = range.distance_to_end();
  const Vector<double, D> *pbegin = &get<position>(*range);
  const size_t index = pbegin - &get<position>(source_particles_begin)[0];
  for (size_t i = 0; i < N; ++i) {
    const Vector<double, D> &pi = pbegin[i];
    expansions.P2M(sum, box, pi, source_vector[index + i]);
  }
}

// assume serial processing of particles, this could be more efficient for
// ranges iterators
template <typename Expansions, typename Iterator, typename T,
          typename VectorType = typename Expansions::expansion_type,
          typename Traits = typename Iterator::traits_type,
          typename SourceParticleIterator = typename Traits::raw_pointer,
          unsigned int D = Traits::dimension,
          typename = typename std::enable_if<
              !std::is_same<Iterator, ranges_iterator<Traits>>::value>>
void calculate_P2M(VectorType &sum, const bbox<D> &box, const Iterator &range,
                   const std::vector<T> &source_vector,
                   const SourceParticleIterator &source_particles_begin,
                   const Expansions &expansions) {

  typedef typename Traits::position position;
  for (auto i = range; i != false; ++i) {
    const Vector<double, D> &pi = get<position>(*i);
    const size_t index = &pi - &get<position>(source_particles_begin)[0];
    expansions.P2M(sum, box, pi, source_vector[index]);
  }
}

#ifdef HAVE_EIGEN

// TODO: move this somewhere sensible
template <size_t Size> struct ConvertToDoubleOrVector {
  typedef
      typename std::conditional<Size == 1, double, Vector<double, Size>>::type
          return_type;

  template <typename Derived>
  static return_type convert(const Eigen::DenseBase<Derived> &vector) {
    return vector;
  }
};

template <> struct ConvertToDoubleOrVector<1> {
  typedef double return_type;

  template <typename Derived>
  static return_type convert(const Eigen::DenseBase<Derived> &vector) {
    return vector[0];
  }
};

template <typename Expansions, typename Traits, typename Derived,
          typename VectorType = typename Expansions::expansion_type,
          typename SourceParticleIterator = typename Traits::raw_pointer,
          unsigned int D = Traits::dimension>
void calculate_P2M(VectorType &sum, const bbox<D> &box,
                   const ranges_iterator<Traits> &range,
                   const Eigen::DenseBase<Derived> &source_vector,
                   const SourceParticleIterator &source_particles_begin,
                   const Expansions &expansions) {
  typedef typename Traits::position position;
  const size_t N = range.distance_to_end();
  const Vector<double, D> *pbegin = &get<position>(*range);
  const size_t index = pbegin - &get<position>(source_particles_begin)[0];
  const size_t block_size = Expansions::block_cols;
  for (size_t i = 0; i < N; ++i) {
    const Vector<double, D> &pi = pbegin[i];
    expansions.P2M(sum, box, pi,
                   ConvertToDoubleOrVector<block_size>::convert(
                       source_vector.template segment<block_size>((index + i) *
                                                                  block_size)));
  }
}

// assume serial processing of particles, this could be more efficient for
// ranges iterators
template <typename Expansions, typename Iterator, typename Derived,
          typename VectorType = typename Expansions::expansion_type,
          typename Traits = typename Iterator::traits_type,
          typename SourceParticleIterator = typename Traits::raw_pointer,
          unsigned int D = Traits::dimension,
          typename = typename std::enable_if<
              !std::is_same<Iterator, ranges_iterator<Traits>>::value>>
void calculate_P2M(VectorType &sum, const bbox<D> &box, const Iterator &range,
                   const Eigen::DenseBase<Derived> &source_vector,
                   const SourceParticleIterator &source_particles_begin,
                   const Expansions &expansions) {

  typedef typename Traits::position position;
  constexpr size_t block_size = Expansions::block_cols;
  for (auto i = range; i != false; ++i) {
    const Vector<double, D> &pi = get<position>(*i);
    const size_t index = &pi - &get<position>(source_particles_begin)[0];
    expansions.P2M(
        sum, box, pi,
        ConvertToDoubleOrVector<block_size>::convert(
            source_vector.template segment<block_size>(index * block_size)));
  }
}
#endif

template <typename Expansions, typename Traits, typename T,
          typename VectorType = typename Expansions::expansion_type,
          typename ParticleIterator = typename Traits::raw_pointer,
          unsigned int D = Traits::dimension>
void calculate_L2P(std::vector<T> &target_vector, const VectorType &source,
                   const bbox<D> &box, const ranges_iterator<Traits> &range,
                   const ParticleIterator &target_particles_begin,
                   const Expansions &expansions) {
  typedef typename Traits::position position;
  LOG(3, "calculate_L2P (range): box = " << box);
  const size_t N = range.distance_to_end();
  const Vector<double, D> *pbegin_range = &get<position>(*range);
  const Vector<double, D> *pbegin = &get<position>(target_particles_begin)[0];
  const size_t index = pbegin_range - pbegin;
  for (size_t i = index; i < index + N; ++i) {
    const Vector<double, D> &pi = pbegin[i];
    target_vector[i] += expansions.L2P(pi, box, source);
  }
}

// assume serial processing of particles, this could be more efficient for
// ranges iterators
template <typename Expansions, typename Iterator, typename T,
          typename VectorType = typename Expansions::expansion_type,
          typename Traits = typename Iterator::traits_type,
          typename ParticleIterator = typename Traits::raw_pointer,
          unsigned int D = Traits::dimension,
          typename = typename std::enable_if<
              !std::is_same<Iterator, ranges_iterator<Traits>>::value>>
void calculate_L2P(std::vector<T> &target_vector, const VectorType &source,
                   const bbox<D> &box, const Iterator &range,
                   const ParticleIterator &target_particles_begin,
                   const Expansions &expansions) {

  LOG(3, "calculate_L2P: box = " << box);
  typedef typename Traits::position position;
  for (auto i = range; i != false; ++i) {
    const Vector<double, D> &pi = get<position>(*i);
    const size_t index = &pi - &get<position>(target_particles_begin)[0];
    target_vector[index] += expansions.L2P(pi, box, source);
  }
}

#ifdef HAVE_EIGEN
template <typename Expansions, typename Traits, typename Derived,
          typename VectorType = typename Expansions::expansion_type,
          typename ParticleIterator = typename Traits::raw_pointer,
          unsigned int D = Traits::dimension>
void calculate_L2P(const Eigen::DenseBase<Derived> &target_vector,
                   const VectorType &source, const bbox<D> &box,
                   const ranges_iterator<Traits> &range,
                   const ParticleIterator &target_particles_begin,
                   const Expansions &expansions) {
  typedef typename Traits::position position;
  LOG(3, "calculate_L2P (range): box = " << box);
  const size_t N = range.distance_to_end();
  const Vector<double, D> *pbegin_range = &get<position>(*range);
  const Vector<double, D> *pbegin = &get<position>(target_particles_begin)[0];
  const size_t index = pbegin_range - pbegin;
  const size_t block_size = Expansions::block_rows;
  for (size_t i = index; i < index + N; ++i) {
    const Vector<double, D> &pi = pbegin[i];
    const auto l2p = expansions.L2P(pi, box, source);
    for (size_t ii = 0; ii < block_size; ++ii) {
      const_cast<Eigen::DenseBase<Derived> &>(
          target_vector)[i * block_size + ii] +=
          detail::VectorTraits<decltype(l2p)>::Index(l2p, ii);
    }
  }
}

// assume serial processing of particles, this could be more efficient for
// ranges iterators
template <typename Expansions, typename Iterator, typename Derived,
          typename VectorType = typename Expansions::expansion_type,
          typename Traits = typename Iterator::traits_type,
          typename ParticleIterator = typename Traits::raw_pointer,
          unsigned int D = Traits::dimension,
          typename = typename std::enable_if<
              !std::is_same<Iterator, ranges_iterator<Traits>>::value>>
void calculate_L2P(const Eigen::DenseBase<Derived> &target_vector,
                   const VectorType &source, const bbox<D> &box,
                   const Iterator &range,
                   const ParticleIterator &target_particles_begin,
                   const Expansions &expansions) {

  LOG(3, "calculate_L2P: box = " << box);
  typedef typename Traits::position position;
  const size_t block_size = Expansions::block_rows;
  for (auto i = range; i != false; ++i) {
    const Vector<double, D> &pi = get<position>(*i);
    const size_t index = &pi - &get<position>(target_particles_begin)[0];
    const auto &l2p = expansions.L2P(pi, box, source);
    for (size_t ii = 0; ii < block_size; ++ii) {
      const_cast<Eigen::DenseBase<Derived> &>(
          target_vector)[index * block_size + ii] +=
          detail::VectorTraits<decltype(l2p)>::Index(l2p, ii);
    }
  }
}
#endif

template <typename Kernel, typename Traits, typename TargetType,
          typename SourceType,
          typename ParticleIterator = typename Traits::raw_pointer,
          unsigned int D = Traits::dimension>
void calculate_P2P(std::vector<TargetType> &target_vector,
                   const std::vector<SourceType> &source_vector,
                   const ranges_iterator<Traits> &target_range,
                   const ranges_iterator<Traits> &source_range,
                   const ParticleIterator &target_particles_begin,
                   const ParticleIterator &source_particles_begin,
                   const Kernel &kernel) {
  typedef typename Traits::position position;

  const size_t n_target = target_range.distance_to_end();
  const size_t n_source = source_range.distance_to_end();

  const Vector<double, D> *pbegin_target_range = &get<position>(*target_range);
  const Vector<double, D> *pbegin_target =
      &get<position>(target_particles_begin)[0];
  const size_t index_target = pbegin_target_range - pbegin_target;

  const Vector<double, D> *pbegin_source_range = &get<position>(*source_range);
  const Vector<double, D> *pbegin_source =
      &get<position>(source_particles_begin)[0];

  const size_t index_source = pbegin_source_range - pbegin_source;

  auto pi = target_range;
  for (size_t i = index_target; i < index_target + n_target; ++i, ++pi) {
    auto pj = source_range;
    for (size_t j = index_source; j < index_source + n_source; ++j, ++pj) {
      target_vector[i] += kernel(*pi, *pj) * source_vector[j];
    }
  }
}

// assume serial processing of particles, this could be more efficient for
// ranges iterators
template <typename Kernel, typename TargetIterator, typename SourceIterator,
          typename TargetType, typename SourceType,
          typename Traits = typename TargetIterator::traits_type,
          typename ParticleIterator = typename Traits::raw_pointer,
          unsigned int D = Traits::dimension,
          typename = typename std::enable_if<
              !(std::is_same<TargetIterator, ranges_iterator<Traits>>::value &&
                std::is_same<SourceIterator, ranges_iterator<Traits>>::value)>>
void calculate_P2P(std::vector<TargetType> &target_vector,
                   const std::vector<SourceType> &source_vector,
                   const TargetIterator &target_range,
                   const SourceIterator &source_range,
                   const ParticleIterator &target_particles_begin,
                   const ParticleIterator &source_particles_begin,
                   const Kernel &kernel) {

  typedef typename Traits::position position;
  for (auto i = target_range; i != false; ++i) {
    const size_t target_index =
        &get<position>(*i) - &get<position>(target_particles_begin)[0];
    for (auto j = source_range; j != false; ++j) {
      const size_t source_index =
          &get<position>(*j) - &get<position>(source_particles_begin)[0];
      LOG(4, "calculate_P2P: i = " << target_index << " j = " << source_index);
      target_vector[target_index] +=
          kernel(*i, *j) * source_vector[source_index];
    }
  }
}

#ifdef HAVE_EIGEN
template <typename Kernel, typename Traits, typename TargetDerived,
          typename SourceDerived,
          typename ParticleIterator = typename Traits::raw_pointer,
          unsigned int D = Traits::dimension>
void calculate_P2P(const Eigen::DenseBase<TargetDerived> &target_vector,
                   const Eigen::DenseBase<SourceDerived> &source_vector,
                   const ranges_iterator<Traits> &target_range,
                   const ranges_iterator<Traits> &source_range,
                   const ParticleIterator &target_particles_begin,
                   const ParticleIterator &source_particles_begin,
                   const Kernel &kernel) {
  typedef typename Traits::position position;
  typedef typename Traits::raw_const_reference const_row_reference;
  typedef typename Traits::raw_const_reference const_col_reference;

  const size_t n_target = target_range.distance_to_end();
  const size_t n_source = source_range.distance_to_end();

  const Vector<double, D> *pbegin_target_range = &get<position>(*target_range);
  const Vector<double, D> *pbegin_target =
      &get<position>(target_particles_begin)[0];
  const size_t index_target = pbegin_target_range - pbegin_target;

  const Vector<double, D> *pbegin_source_range = &get<position>(*source_range);
  const Vector<double, D> *pbegin_source =
      &get<position>(source_particles_begin)[0];

  const size_t index_source = pbegin_source_range - pbegin_source;

  typedef detail::kernel_helper_ref<const_row_reference, const_col_reference,
                                    Kernel>
      helper;

  auto pi = target_range;
  for (size_t i = index_target; i < index_target + n_target; ++i, ++pi) {
    auto pj = source_range;
    for (size_t j = index_source; j < index_source + n_source; ++j, ++pj) {
      const_cast<Eigen::DenseBase<TargetDerived> &>(target_vector)
          .template segment<helper::block_rows>(i * helper::block_rows) +=
          kernel(*pi, *pj) * source_vector.template segment<helper::block_cols>(
                                 j * helper::block_cols);
    }
  }
}

// assume serial processing of particles, this could be more efficient for
// ranges iterators
template <typename Kernel, typename TargetIterator, typename SourceIterator,
          typename TargetDerived, typename SourceDerived,
          typename Traits = typename TargetIterator::traits_type,
          typename ParticleIterator = typename Traits::raw_pointer,
          unsigned int D = Traits::dimension,
          typename = typename std::enable_if<
              !(std::is_same<TargetIterator, ranges_iterator<Traits>>::value &&
                std::is_same<SourceIterator, ranges_iterator<Traits>>::value)>>
void calculate_P2P(const Eigen::DenseBase<TargetDerived> &target_vector,
                   const Eigen::DenseBase<SourceDerived> &source_vector,
                   const TargetIterator &target_range,
                   const SourceIterator &source_range,
                   const ParticleIterator &target_particles_begin,
                   const ParticleIterator &source_particles_begin,
                   const Kernel &kernel) {

  typedef typename Traits::position position;
  typedef typename Traits::raw_const_reference const_row_reference;
  typedef typename Traits::raw_const_reference const_col_reference;

  // TODO: lots of common code here to consolodate
  typedef typename std::result_of<Kernel(
      const_row_reference, const_col_reference)>::type FunctionReturn;
  typedef typename std::conditional<std::is_arithmetic<FunctionReturn>::value,
                                    Eigen::Matrix<FunctionReturn, 1, 1>,
                                    FunctionReturn>::type Block;

  const int block_rows = Block::RowsAtCompileTime;
  const int block_cols = Block::ColsAtCompileTime;

  static_assert(block_rows > 0,
                "kernel function must return fixed size matrix");
  static_assert(block_cols > 0,
                "kernel function must return fixed size matrix");

  for (auto i = target_range; i != false; ++i) {
    const size_t target_index =
        &get<position>(*i) - &get<position>(target_particles_begin)[0];
    for (auto j = source_range; j != false; ++j) {
      const size_t source_index =
          &get<position>(*j) - &get<position>(source_particles_begin)[0];
      LOG(4, "calculate_P2P: i = " << target_index << " j = " << source_index);
      const_cast<Eigen::DenseBase<TargetDerived> &>(target_vector)
          .template segment<block_rows>(target_index * block_rows) +=
          kernel(*i, *j) *
          source_vector.template segment<block_cols>(source_index * block_cols);
    }
  }
}
#endif

/*
template <typename Traits,
          typename Kernel,
          typename SourceVectorType,
          typename SourceParticleIterator=typename Traits::raw_pointer,
          unsigned int D=Traits::dimension>
double calculate_P2P_single(const Vector<double,D>& p,
                        const iterator_range<ranges_iterator<Traits>>& range,
                        const Kernel& kernel,
                        const SourceVectorType& source_vector,
                        const SourceParticleIterator& source_particles_begin) {
    typedef typename Traits::position position;
    typedef Vector<double,D> double_d;
    const size_t N = std::distance(range.begin(),range.end());
    const double_d* pbegin = &get<position>(*range.begin());
    const size_t index = pbegin - &get<position>(source_particles_begin)[0];
    double sum = 0;
    for (size_t i = 0; i < N; ++i) {
        const Vector<double,D>& pi = pbegin[i];
        sum += kernel(p,pi)*source_vector[index+i];
    }
    return sum;
}

template <typename Iterator,
          typename Kernel,
          typename SourceVectorType,
                typename Traits=typename Iterator::traits_type,
                typename SourceParticleIterator=typename Traits::raw_pointer,
                unsigned int D=Traits::dimension,
                typename = typename
    std::enable_if<!std::is_same<Iterator,ranges_iterator<Traits>>::value>>
double calculate_P2P_position(const Vector<double,D>& p,
                        const iterator_range<Iterator>& range,
                        const Kernel& kernel,
                        const SourceVectorType& source_vector,
                        const SourceParticleIterator& source_particles_begin) {
    typedef typename Traits::position position;
    typedef typename Iterator::reference reference;

    double sum = 0;
    for (reference i: range) {
        const Vector<double,D>& pi = get<position>(i);
        const size_t index = &pi- &get<position>(source_particles_begin)[0];
        sum += kernel(p,pi)*source_vector[index];
    }
    return sum;
}
*/

#ifdef HAVE_EIGEN
template <typename RowParticlesType, typename ColParticlesType, typename Kernel>
void P2P_matrix(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &matrix,
                const std::vector<size_t> &row_indicies,
                const std::vector<size_t> &col_indicies,
                const RowParticlesType &row_particles,
                const ColParticlesType &col_particles, const Kernel &kernel) {

  typedef detail::kernel_helper_ref<typename RowParticlesType::const_reference,
                                    typename ColParticlesType::const_reference,
                                    Kernel>
      helper;

  matrix.resize(row_indicies.size() * helper::block_rows,
                col_indicies.size() * helper::block_cols);
  for (size_t i = 0; i < row_indicies.size(); ++i) {
    for (size_t j = 0; j < col_indicies.size(); ++j) {
      matrix.block<helper::block_rows, helper::block_cols>(
          i * helper::block_rows, j * helper::block_cols) =
          kernel(row_particles[i], col_particles[j]);
    }
  }
}
#endif

#ifdef HAVE_H2LIB
template <typename RowParticlesType, typename ColParticlesType, typename Kernel>
void P2P_amatrix(pamatrix matrix, const uint *row_indicies,
                 const uint row_indicies_size, const uint *col_indicies,
                 const uint col_indicies_size,
                 const RowParticlesType &row_particles,
                 const ColParticlesType &col_particles, const Kernel &kernel) {

  typedef detail::kernel_helper_ref<typename RowParticlesType::const_reference,
                                    typename ColParticlesType::const_reference,
                                    Kernel>
      helper;

  ASSERT_CUDA(matrix->rows == row_indicies_size);
  ASSERT_CUDA(matrix->cols == col_indicies_size);

  // resize_amatrix(matrix,row_indicies_size,col_indicies_size);
  for (size_t i = 0; i < row_indicies_size; i += helper::block_rows) {
    const auto &pi = row_particles[row_indicies[i] / helper::block_rows];
    for (size_t j = 0; j < col_indicies_size; j += helper::block_cols) {
      const auto &pj = col_particles[col_indicies[j] / helper::block_cols];
      const typename helper::Block tmp(kernel(pi, pj));
      for (size_t ii = 0; ii < helper::block_rows; ++ii) {
        for (size_t jj = 0; jj < helper::block_cols; ++jj) {
          setentry_amatrix(matrix, i + ii, j + jj, tmp(ii, jj));
        }
      }
    }
  }
  /*
  std::cout << "P2P = " << std::endl;
  for (int i = 0; i < row_indicies_size; ++i) {
    std::cout << "| ";
    for (int j = 0; j < col_indicies_size; ++j) {
      std::cout << getentry_amatrix(matrix, i, j) << " ";
    }
    std::cout << "|" << std::endl;
  }
  */
}
#endif

/*
template <unsigned int D>
struct theta_condition {
    typedef Vector<double,D> double_d;
    const double_d& m_low;
    const double_d& m_high;
    const double m_r2;
    const double m_r;
    static constexpr double m_theta = 0.5;
    static constexpr double m_theta2 = m_theta*m_theta;
    theta_condition(const double_d& low, const double_d& high):
        m_low(low),m_high(high),
        m_r2(0.25*(high-low).squaredNorm()),
        m_r(std::sqrt(m_r2))
    {}

    bool check(const double_d& low, const double_d& high) const {
        double d = 0.5*(high + low - m_low - m_high).norm();
        double other_r2 = 0.25*(high-low).squaredNorm();
        if (other_r2 < m_r2) {
            const double other_r = std::sqrt(other_r2);
            return m_r2 > m_theta2*std::pow(d-other_r,2);
        } else {
            return other_r2 > m_theta2*std::pow(d-m_r,2);
        }
    }
};
*/

template <unsigned int D> struct theta_condition {
  typedef Vector<double, D> double_d;
  const double_d &m_low;
  const double_d &m_high;
  double m_max_diam;
  static constexpr double m_theta = 0.5;
  static constexpr double m_theta2 = m_theta * m_theta;
  theta_condition(const double_d &low, const double_d &high)
      : m_low(low), m_high(high),
        m_max_diam(std::numeric_limits<double>::min()) {
    for (size_t i = 0; i < D; ++i) {
      if (high[i] - low[i] > m_max_diam)
        m_max_diam = high[i] - low[i];
    }
  }

  bool check(const double_d &low, const double_d &high) const {
    double_d dist = double_d::Zero();
    double max_diam = m_max_diam;
    for (size_t i = 0; i < D; ++i) {
      if (m_high[i] < low[i]) {
        dist[i] = low[i] - m_high[i];
      } else if (high[i] < m_low[i]) {
        dist[i] = m_low[i] - high[i];
      }
      if (high[i] - low[i] > max_diam)
        max_diam = high[i] - low[i];
    }
    return std::pow(max_diam, 2) > 1.0 * dist.squaredNorm();
  }
};

} // namespace detail
} // namespace Aboria

#endif
