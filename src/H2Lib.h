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

#ifndef H2_LIB_H_
#define H2_LIB_H_

#include "detail/FastMultipoleMethod.h"

#ifdef HAVE_H2LIB

#include "detail/H2Lib.h"

namespace Aboria {

class H2Lib_LR_Decomposition {
  std::unique_ptr<h2matrix, decltype(&del_h2matrix)> A;
  std::unique_ptr<h2matrix, decltype(&del_h2matrix)> L;
  std::unique_ptr<h2matrix, decltype(&del_h2matrix)> R;
  ptruncmode tm;

public:
  H2Lib_LR_Decomposition(const ph2matrix h2mat, const pblock broot,
                         const double tol)
      : A(nullptr, del_h2matrix), L(nullptr, del_h2matrix),
        R(nullptr, del_h2matrix), tm(new_releucl_truncmode()) {

    pclusterbasis rbcopy = clone_clusterbasis(h2mat->rb);
    pclusterbasis cbcopy = clone_clusterbasis(h2mat->cb);
    A = std::unique_ptr<h2matrix, decltype(&del_h2matrix)>(
        clone_h2matrix(h2mat, rbcopy, cbcopy), del_h2matrix);

    pclusteroperator Arwf = prepare_row_clusteroperator(A->rb, A->cb, tm);
    pclusteroperator Acwf = prepare_col_clusteroperator(A->rb, A->cb, tm);

    pclusterbasis Lrcb = build_from_cluster_clusterbasis(A->rb->t);
    pclusterbasis Lccb = build_from_cluster_clusterbasis(A->cb->t);
    pclusterbasis Rrcb = build_from_cluster_clusterbasis(A->rb->t);
    pclusterbasis Rccb = build_from_cluster_clusterbasis(A->cb->t);

    L = std::unique_ptr<h2matrix, decltype(&del_h2matrix)>(
        build_from_block_lower_h2matrix(broot, Lrcb, Lccb), del_h2matrix);
    R = std::unique_ptr<h2matrix, decltype(&del_h2matrix)>(
        build_from_block_upper_h2matrix(broot, Rrcb, Rccb), del_h2matrix);

    pclusteroperator Lrwf = prepare_row_clusteroperator(L->rb, L->cb, tm);
    pclusteroperator Lcwf = prepare_col_clusteroperator(L->rb, L->cb, tm);
    pclusteroperator Rrwf = prepare_row_clusteroperator(R->rb, R->cb, tm);
    pclusteroperator Rcwf = prepare_col_clusteroperator(R->rb, R->cb, tm);

    lrdecomp_h2matrix(A.get(), Arwf, Acwf, L.get(), Lrwf, Lcwf, R.get(), Rrwf,
                      Rcwf, tm, tol);
  }

  // TODO: match eigen's interface for solver
  /*
  template <typename DerivedRHS>
  void solve(const Eigen::DenseBase<DerivedRHS> &source,
                   Eigen::Matrix<double,Eigen::Dynamic,1> &dest) {
      ASSERT(source.cols() == 1 || source.rows() == 1,"solve must take
  vectors");

      dest = source;
      pavector x = new_pointer_avector(
                         dest.data(),
                         dest.size());
      lrsolve_h2matrix_avector(L,R,x);
  }
  */

  template <typename T1, typename T2>
  void solve(const std::vector<T1> &source, std::vector<T2> &dest) {
    static_assert(std::is_standard_layout<T2>::value,
                  "T2 must be standard layout");
    typedef typename detail::VectorTraits<T2> traitsT2;

    detail::copy(std::begin(source), std::end(source), std::begin(dest));

    pavector dest_avector =
        new_pointer_avector(reinterpret_cast<double *>(dest.data()),
                            dest.size() * traitsT2::length);

    lrsolve_h2matrix_avector(L.get(), R.get(), dest_avector);
  }
};

class H2LibCholeskyDecomposition {
  std::unique_ptr<h2matrix, decltype(&del_h2matrix)> A;
  std::unique_ptr<h2matrix, decltype(&del_h2matrix)> L;
  ptruncmode tm;

public:
  H2LibCholeskyDecomposition(const ph2matrix h2mat, const pblock broot,
                             const double tol)
      : A(nullptr, del_h2matrix), L(nullptr, del_h2matrix),
        tm(new_releucl_truncmode()) {

    pclusterbasis rbcopy = clone_clusterbasis(h2mat->rb);
    pclusterbasis cbcopy = clone_clusterbasis(h2mat->cb);
    A = std::unique_ptr<h2matrix, decltype(&del_h2matrix)>(
        clone_h2matrix(h2mat, rbcopy, cbcopy), del_h2matrix);

    pclusterbasis Lrcb = build_from_cluster_clusterbasis(h2mat->rb->t);
    pclusterbasis Lccb = build_from_cluster_clusterbasis(h2mat->rb->t);
    ASSERT_CUDA(h2mat->rb->t == h2mat->cb->t);

    L = std::unique_ptr<h2matrix, decltype(&del_h2matrix)>(
        build_from_block_lower_h2matrix(broot, Lrcb, Lccb), del_h2matrix);

    pclusteroperator Lrwf = prepare_row_clusteroperator(L->rb, L->cb, tm);
    pclusteroperator Lcwf = prepare_col_clusteroperator(L->rb, L->cb, tm);

    pclusteroperator Arwf = nullptr;
    pclusteroperator Acwf = nullptr;
    init_cholesky_h2matrix(A.get(), &Arwf, &Acwf, tm);

    choldecomp_h2matrix(A.get(), Arwf, Acwf, L.get(), Lrwf, Lcwf, tm, tol);

    /*
    pclusteroperator rwfh2 = prepare_row_clusteroperator(A->rb, A->cb, tm);
    pclusteroperator cwfh2 = prepare_col_clusteroperator(A->rb, A->cb, tm);
    const double error = norm2_h2matrix(A);
    addmul_h2matrix(-1.0, L, true, L, A, rwfh2, cwfh2, tm, tol);
    std::cout << "error in factorisation is
    "<<norm2_h2matrix(A)/error<<std::endl;
    */
  }

  // TODO: match eigen's interface for solver
  template <typename DerivedRHS>
  void solve(const Eigen::DenseBase<DerivedRHS> &source,
                   Eigen::Matrix<double,Eigen::Dynamic,1> &dest) {
      ASSERT(source.cols() == 1 || source.rows() == 1,"solve must take vectors"); 
      dest = source; 
      pavector dest_avector = new_pointer_avector(reinterpret_cast<double *>(dest.data()), dest.size());

      cholsolve_h2matrix_avector(L.get(),dest_avector);
  }

  template <typename T1, typename T2>
  void solve(const std::vector<T1> &source, std::vector<T2> &dest) {
    static_assert(std::is_standard_layout<T2>::value,
                  "T2 must be standard layout");
    typedef typename detail::VectorTraits<T2> traitsT2;

    detail::copy(std::begin(source), std::end(source), std::begin(dest));

    pavector dest_avector =
        new_pointer_avector(reinterpret_cast<double *>(dest.data()),
                            dest.size() * traitsT2::length);

    cholsolve_h2matrix_avector(L.get(), dest_avector);
  }

  ///
  /// @brief returns the determinant of the H2 matrix that was used to make this
  /// decomposition
  ///
  /// @return double
  /// @return NaN if negative eigenvalues found (i.e. H2 matrix not positive
  /// definite)
  double determinant() const { return std::exp(log_determinant()); }

  ///
  /// @brief returns the log determinant of the H2 matrix that was used to make
  /// this decomposition.
  ///
  /// @return double
  /// @return NaN if negative eigenvalues found (i.e. H2 matrix not positive
  /// definite)
  /// @return -HUGE_VAL if zero eigenvalues found
  double log_determinant() const { return 2 * log_determinant_sum(L.get()); }

private:
  double log_determinant_sum(ph2matrix h2) const {
    ASSERT(h2->rsons == h2->csons, "L matrix not square!");
    // sum up in log space
    double ld = 0;
    // inaccessible leaf
    if (h2->f) {
      ASSERT(h2->f->rows == h2->f->cols, "diagonal f matrix not square");
      for (size_t i = 0; i < h2->f->rows; ++i) {
        ld += std::log(getentry_amatrix(h2->f, i, i));
      }
    } else {
      ASSERT(!h2->u, "found accessible leaf on diagonal!");
      for (size_t i = 0; i < h2->rsons; ++i) {
        ld += log_determinant_sum(h2->son[i + i * h2->rsons]);
      }
    }
    return ld;
  }
};

class HLib_LR_Decomposition {
  phmatrix A;
  ptruncmode tm;

public:
  HLib_LR_Decomposition(const phmatrix hmat, const double tol)
      : tm(new_releucl_truncmode()) {
    A = clone_hmatrix(hmat);
    lrdecomp_hmatrix(A, tm, tol);
  }

  ~HLib_LR_Decomposition() { del_hmatrix(A); }

  // TODO: match eigen's interface for solver
  /*
  template <typename DerivedRHS>
  void solve(const Eigen::DenseBase<DerivedRHS> &source,
                   Eigen::Matrix<double,Eigen::Dynamic,1> &dest) {
      ASSERT(source.cols() == 1 || source.rows() == 1,"solve must take
  vectors"); dest = source; pavector x = new_pointer_avector( dest.data(),
                         dest.size());
      lrsolve_hmatrix_avector(false,A,x);
  }
  */

  template <typename T1, typename T2>
  void solve(const std::vector<T1> &source, std::vector<T2> &dest) {
    static_assert(std::is_standard_layout<T2>::value,
                  "T2 must be standard layout");
    typedef typename detail::VectorTraits<T2> traitsT2;

    detail::copy(std::begin(source), std::end(source), std::begin(dest));

    pavector dest_avector =
        new_pointer_avector(reinterpret_cast<double *>(dest.data()),
                            dest.size() * traitsT2::length);

    lrsolve_hmatrix_avector(false, A, dest_avector);
  }
};

class HLibCholeskyDecomposition {
  phmatrix A;
  ptruncmode tm;

public:
  HLibCholeskyDecomposition(const phmatrix hmat, const double tol)
      : tm(new_abseucl_truncmode()) {

    A = clone_hmatrix(hmat);
    choldecomp_hmatrix(A, tm, tol);
  }

  ~HLibCholeskyDecomposition() { del_hmatrix(A); }

  // TODO: match eigen's interface for solver
  template <typename DerivedRHS>
  void solve(const Eigen::DenseBase<DerivedRHS> &source,
             Eigen::Matrix<double, Eigen::Dynamic, 1> &dest) {
    ASSERT(source.cols() == 1 || source.rows() == 1, "solve must take vectors");

    dest = source;
    pavector x = new_pointer_avector(dest.data(), dest.size());
    cholsolve_hmatrix_avector(A, x);
  }

  template <typename T>
  void solve(const std::vector<T> &source, std::vector<double> &dest) {
    detail::copy(std::begin(source), std::end(source), std::begin(dest));
    pavector x = new_pointer_avector(dest.data(), dest.size());
    cholsolve_hmatrix_avector(A, x);
  }
};

class H2LibMatrix {

  std::unique_ptr<h2matrix, decltype(&del_h2matrix)> m_h2;
  std::unique_ptr<block, decltype(&del_block)> m_block;
  std::vector<uint> m_row_idx;
  std::vector<uint> m_col_idx;


public:
  template <typename RowParticles, typename ColParticles, typename Expansions,
            typename Kernel>
  H2LibMatrix(const RowParticles &row_particles,
              const ColParticles &col_particles, const Expansions &expansions,
              const Kernel &kernel, const double eta)
      : m_h2(nullptr, del_h2matrix), m_block(nullptr, del_block) {

    // generate h2 matrix
    LOG(2, "H2LibMatrix: creating h2 matrix using "
               << row_particles.size() << " row particles and "
               << col_particles.size()
               << " column particles and eta = " << eta);

    const bool row_equals_col = static_cast<const void *>(&row_particles) ==
                                static_cast<const void *>(&col_particles);

    //
    // Create row clusters and clusterbasis
    //
    m_row_idx.resize(row_particles.size() * expansions.block_rows);
    pcluster row_t =
        set_root_idx(m_row_idx.data(), row_particles, expansions.block_rows);
    pclusterbasis row_cb = build_from_cluster_clusterbasis(row_t);
    auto data_row = std::make_tuple(&expansions, &row_particles);
    iterate_parallel_clusterbasis(
        row_cb, 0, max_pardepth, NULL,
        detail::assemble_h2matrix_row_clusterbasis<RowParticles, Expansions>,
        &data_row);

    //
    // Create col clusters and clusterbasis
    //
    pcluster col_t;
    if (row_equals_col && expansions.block_rows == expansions.block_cols) {
      LOG(2, "H2LibMatrix: row clusters the same as column clusters");
      col_t = row_t;
    } else {
      m_col_idx.resize(col_particles.size() * expansions.block_cols);
      col_t =
          set_root_idx(m_col_idx.data(), col_particles, expansions.block_cols);
    }
    pclusterbasis col_cb = build_from_cluster_clusterbasis(col_t);
    auto data_col = std::make_tuple(&expansions, &col_particles);
    iterate_parallel_clusterbasis(
        col_cb, 0, max_pardepth, NULL,
        detail::assemble_h2matrix_col_clusterbasis<ColParticles, Expansions>,
        &data_col);

    //
    // create h2 block
    //
    double eta_copy = eta;
    m_block = std::unique_ptr<block, decltype(&del_block)>(
        build_strict_block(row_t, col_t, &eta_copy, detail::admissible_max_cluster),
        del_block);

    /*
    int argc = 1;
    char name[10] = "test";
    char *argv[1];
    argv[0] = name;
    glutInit(&argc, argv);
    view_block(m_block.get());
    */

    m_h2 = std::unique_ptr<h2matrix, decltype(&del_h2matrix)>(
        build_from_block_h2matrix(m_block.get(), row_cb, col_cb), del_h2matrix);
    ph2matrix *enum_h2mat = enumerate_h2matrix(m_block.get(), m_h2.get());
    auto data_h2 = std::make_tuple(&expansions, &kernel, &row_particles,
                                   &col_particles, enum_h2mat);
    iterate_byrow_block(
        m_block.get(), 0, 0, 0, max_pardepth, NULL,
        detail::assemble_block_h2matrix<RowParticles, ColParticles, Expansions, Kernel>,
        &data_h2);
    freemem(enum_h2mat);
  }

  void compress(const double tol) {
    // ptruncmode tm = new_releucl_truncmode();
    ptruncmode tm = new_blockreleucl_truncmode();
    recompress_inplace_h2matrix(m_h2.get(), tm, tol);
  }

private:
  template <typename Particles>
  pcluster set_root_idx(uint *idx, const Particles &particles,
                        const size_t block_size) {
    uint *old_idx = idx;
    const uint dim = Particles::dimension;
    size_t sons = 0;
    for (auto ci_child = particles.get_query().get_children();
         ci_child != false; ++ci_child, ++sons) {
    }
    pcluster t = new_cluster(0, old_idx, sons, dim);
    size_t i = 0;
    for (auto ci_child = particles.get_query().get_children();
         ci_child != false; ++ci_child, ++i) {
      t->son[i] = set_idx(ci_child, idx, particles, block_size);
      idx += t->son[i]->size;
    }
    t->size = idx - old_idx;
    // bounding box
    for (size_t i = 0; i < dim; ++i) {
      t->bmin[i] = particles.get_min()[i];
      t->bmax[i] = particles.get_max()[i];
    }
    update_cluster(t);
    return t;
  }

  template <typename ChildIterator, typename Particles>
  pcluster set_idx(const ChildIterator ci, uint *idx,
                   const Particles &particles, const size_t block_size) {
    uint *old_idx = idx;
    const uint dim = Particles::dimension;
    pcluster t;
    if (particles.get_query().is_leaf_node(*ci)) { // leaf node
      for (auto pit = particles.get_query().get_bucket_particles(*ci);
           pit != false; ++pit) {
        const size_t pi = &(get<typename Particles::position>(*pit)) -
                          &get<typename Particles::position>(particles)[0];
        for (size_t i = 0; i < block_size; ++i, ++idx) {
          *idx = pi * block_size + i;
        }
      }
      // if (idx-old_idx == 0) std::cout << "HAVE EMPTY LEAF" <<std::endl;
      t = new_cluster(idx - old_idx, old_idx, 0, dim);
    } else {
      size_t sons = 0;
      for (auto ci_child = particles.get_query().get_children(ci);
           ci_child != false; ++ci_child, ++sons) {
      }
      t = new_cluster(0, old_idx, sons, dim);
      size_t i = 0;
      for (auto ci_child = particles.get_query().get_children(ci);
           ci_child != false; ++ci_child, ++i) {
        pcluster child_t = set_idx(ci_child, idx, particles, block_size);
        if (child_t->size > 0) {
          t->son[i] = child_t;
          idx += child_t->size;
        } else {
          del_cluster(child_t);
          --i;
          --(t->sons);
        }
      }
      t->size = idx - old_idx;
    }
    // bounding box
    const auto &bbox = particles.get_query().get_bounds(ci);
    for (size_t i = 0; i < dim; ++i) {
      t->bmin[i] = bbox.bmin[i];
      t->bmax[i] = bbox.bmax[i];
    }
    update_cluster(t);
    return t;
  }

public:
  ph2matrix get_ph2matrix() const { return m_h2.get(); }

  pblock get_pblock() const { return m_block.get(); }

  H2Lib_LR_Decomposition lr(const double tol) const {
    return H2Lib_LR_Decomposition(get_ph2matrix(), m_block.get(), tol);
  }

  H2LibCholeskyDecomposition chol(const double tol) const {
    return H2LibCholeskyDecomposition(get_ph2matrix(), m_block.get(), tol);
  }

  // target_vector += alpha*A*source_vector or alpha*A'*source_vector
  template <typename T1, typename T2>
  void matrix_vector_multiply(std::vector<T1> &target_vector,
                              const double alpha, const bool h2trans,
                              const std::vector<T2> &source_vector) const {
    static_assert(std::is_standard_layout<T1>::value,
                  "T1 must be standard layout");
    static_assert(std::is_standard_layout<T2>::value,
                  "T2 must be standard layout");
    typedef typename detail::VectorTraits<T1> traitsT1;
    typedef typename detail::VectorTraits<T2> traitsT2;

    pavector source_avector = new_pointer_avector(
        const_cast<double *>(
            reinterpret_cast<const double *>(source_vector.data())),
        source_vector.size() * traitsT1::length);
    // std::cout << "source: "<<source_vector[0]<<" versus
    // "<<getentry_avector(source_avector,0)<<std::endl;
    pavector target_avector =
        new_pointer_avector(reinterpret_cast<double *>(target_vector.data()),
                            target_vector.size() * traitsT2::length);
    /*
    for(int i=0; i<5; ++i) {
        std::cout << "target: "<<target_vector[i]<<" versus
    "<<getentry_avector(target_avector,i*2)<<std::endl;
        //std::cout << "target: "<<target_vector[i][1]<<" versus
    "<<getentry_avector(target_avector,i*2+1)<<std::endl;
    }
    */
    mvm_h2matrix_avector(alpha, h2trans, m_h2.get(), source_avector,
                         target_avector);
    /*
    for(int i=0; i<5; ++i) {
        std::cout << "after target: "<<target_vector[i]<<" versus
    "<<getentry_avector(target_avector,i*2)<<std::endl;
        //std::cout << "after target: "<<target_vector[i][1]<<" versus
    "<<getentry_avector(target_avector,i*2+1)<<std::endl;
    }
    */
  }

  template <typename T1, typename T2>
  void matrix_vector_multiply(Eigen::DenseBase<T1> &target_vector,
                              const double alpha, const bool h2trans,
                              const Eigen::DenseBase<T2> &source_vector) const {

    pavector source_avector = new_pointer_avector(
        const_cast<double *>(source_vector.derived().data()),
        source_vector.size());

    pavector target_avector = new_pointer_avector(
        // const_cast<double*>(
        target_vector.derived().data(), target_vector.size());

    mvm_h2matrix_avector(alpha, h2trans, m_h2.get(), source_avector,
                         target_avector);
  }
};



class HLibMatrix {

  std::unique_ptr<hmatrix, decltype(&del_hmatrix)> m_h;
  std::unique_ptr<block, decltype(&del_block)> m_block;
  std::vector<uint> m_row_idx;
  std::vector<uint> m_col_idx;

  template <typename RowParticles, typename ColParticles, typename Expansions,
            typename Kernel>
  static void assemble_block_hmatrix(pcblock b, uint bname, uint rname,
                                     uint cname, uint pardepth, void *data) {
    auto &data_cast =
        *static_cast<std::tuple<Expansions *, Kernel *, RowParticles *,
                                ColParticles *, phmatrix *> *>(data);

    const Expansions &expansions = *std::get<0>(data_cast);
    const Kernel &kernel = *std::get<1>(data_cast);
    const RowParticles &row_particles = *std::get<2>(data_cast);
    const ColParticles &col_particles = *std::get<3>(data_cast);
    phmatrix *enum_h = std::get<4>(data_cast);
    phmatrix h = enum_h[bname];

    if (h->r) {
      expansions.L2P_amatrix(&h->r->A, h->rc, h->rc->idx, h->rc->size,
                             row_particles);
      expansions.L2P_amatrix(&h->r->B, h->cc, h->cc->idx, h->cc->size,
                             col_particles);
    } else if (h->f) {
      detail::P2P_amatrix(h->f, h->rc->idx, h->rc->size, h->cc->idx,
                          h->cc->size, row_particles, col_particles, kernel);
    }
  }

  static void truncate_block_hmatrix(pcblock b, uint bname, uint rname,
                                     uint cname, uint pardepth, void *data) {
    auto &data_cast = *static_cast<std::tuple<double *, phmatrix *> *>(data);

    const double &tol = *std::get<0>(data_cast);
    phmatrix *enum_h = std::get<1>(data_cast);
    phmatrix h = enum_h[bname];

    if (h->r) {
      trunc_rkmatrix(NULL, tol, h->r);
    }
  }

public:
  template <typename RowParticles, typename ColParticles, typename Expansions,
            typename Kernel>
  HLibMatrix(const RowParticles &row_particles,
             const ColParticles &col_particles, const Expansions &expansions,
             const Kernel &kernel)
      : m_h(nullptr, del_hmatrix), m_block(nullptr, del_block) {

    // generate h2 matrix
    LOG(2, "H2LibMatrix: creating h matrix using "
               << row_particles.size() << " row particles and "
               << col_particles.size() << " column particles");

    const bool row_equals_col = static_cast<const void *>(&row_particles) ==
                                static_cast<const void *>(&col_particles);

    //
    // Create row clusters and clusterbasis
    //
    m_row_idx.resize(row_particles.size() * expansions.block_rows);
    pcluster row_t =
        set_root_idx(m_row_idx.data(), row_particles, expansions.block_rows);

    //
    // Create col clusters and clusterbasis
    //
    pcluster col_t;
    if (row_equals_col && expansions.block_rows == expansions.block_cols) {
      col_t = row_t;
    } else {
      m_col_idx.resize(col_particles.size() * expansions.block_cols);
      col_t =
          set_root_idx(m_col_idx.data(), col_particles, expansions.block_cols);
    }

    //
    // create h block
    //
    //  eta = 1.0;

    double eta = 1.0;
    m_block = std::unique_ptr<block, decltype(&del_block)>(
        build_nonstrict_block(row_t, col_t, &eta, admissible_max_cluster),
        del_block);

    static_assert(Expansions::block_rows == Expansions::block_cols,
                  "H Matrix only implemented for square kernel matrices");
    m_h = std::unique_ptr<hmatrix, decltype(&del_hmatrix)>(
        build_from_block_hmatrix(m_block.get(),
                                 expansions.m_ncheb * expansions.block_rows),
        del_hmatrix);
    phmatrix *enum_hmat = enumerate_hmatrix(m_block.get(), m_h.get());
    auto data_h = std::make_tuple(&expansions, &kernel, &row_particles,
                                  &col_particles, enum_hmat);
    iterate_byrow_block(
        m_block.get(), 0, 0, 0, max_pardepth, NULL,
        assemble_block_hmatrix<RowParticles, ColParticles, Expansions, Kernel>,
        &data_h);
    freemem(enum_hmat);
  }

private:
  template <typename Particles>
  pcluster set_root_idx(uint *idx, const Particles &particles,
                        const size_t block_size) {
    uint *old_idx = idx;
    const uint dim = Particles::dimension;
    size_t sons = 0;
    for (auto ci_child = particles.get_query().get_children();
         ci_child != false; ++ci_child, ++sons) {
    }
    pcluster t = new_cluster(0, old_idx, sons, dim);
    size_t i = 0;
    for (auto ci_child = particles.get_query().get_children();
         ci_child != false; ++ci_child, ++i) {
      t->son[i] = set_idx(ci_child, idx, particles, block_size);
      idx += t->son[i]->size;
    }
    t->size = idx - old_idx;
    // bounding box
    for (size_t i = 0; i < dim; ++i) {
      t->bmin[i] = particles.get_min()[i];
      t->bmax[i] = particles.get_max()[i];
    }
    update_cluster(t);
    return t;
  }

  template <typename ChildIterator, typename Particles>
  pcluster set_idx(const ChildIterator ci, uint *idx,
                   const Particles &particles, const size_t block_size) {
    uint *old_idx = idx;
    const uint dim = Particles::dimension;
    pcluster t;
    if (particles.get_query().is_leaf_node(*ci)) { // leaf node
      for (auto pit = particles.get_query().get_bucket_particles(*ci);
           pit != false; ++pit) {
        const size_t pi = &(get<typename Particles::position>(*pit)) -
                          &get<typename Particles::position>(particles)[0];
        for (size_t i = 0; i < block_size; ++i, ++idx) {
          *idx = pi * block_size + i;
        }
      }
      // if (idx-old_idx == 0) std::cout << "HAVE EMPTY LEAF" <<std::endl;
      t = new_cluster(idx - old_idx, old_idx, 0, dim);
    } else {
      size_t sons = 0;
      for (auto ci_child = particles.get_query().get_children(ci);
           ci_child != false; ++ci_child, ++sons) {
      }
      t = new_cluster(0, old_idx, sons, dim);
      size_t i = 0;
      for (auto ci_child = particles.get_query().get_children(ci);
           ci_child != false; ++ci_child, ++i) {
        pcluster child_t = set_idx(ci_child, idx, particles, block_size);
        if (child_t->size > 0) {
          t->son[i] = child_t;
          idx += child_t->size;
        } else {
          del_cluster(child_t);
          --i;
          --(t->sons);
        }
      }
      t->size = idx - old_idx;
    }
    // bounding box
    const auto &bbox = particles.get_query().get_bounds(ci);
    for (size_t i = 0; i < dim; ++i) {
      t->bmin[i] = bbox.bmin[i];
      t->bmax[i] = bbox.bmax[i];
    }
    update_cluster(t);
    return t;
  }

public:
  void compress(const double tol) {
    // ptruncmode tm = new_releucl_truncmode();

    phmatrix *enum_hmat = enumerate_hmatrix(m_block.get(), m_h.get());
    auto data_h = std::make_tuple(&tol, enum_hmat);
    iterate_byrow_block(m_block.get(), 0, 0, 0, max_pardepth, NULL,
                        truncate_block_hmatrix, &data_h);
    freemem(enum_hmat);

    coarsen_hmatrix(m_h.get(), NULL, tol, true);
  }

  phmatrix get_phmatrix() const { return m_h.get(); }

  HLib_LR_Decomposition lr(const double tol) const {
    return HLib_LR_Decomposition(get_phmatrix(), tol);
  }

  HLibCholeskyDecomposition chol(const double tol) const {
    return HLibCholeskyDecomposition(get_phmatrix(), tol);
  }

  // target_vector += alpha*A*source_vector or alpha*A'*source_vector
  template <typename T1, typename T2>
  void matrix_vector_multiply(std::vector<T1> &target_vector,
                              const double alpha, const bool h2trans,
                              const std::vector<T2> &source_vector) const {
    static_assert(std::is_standard_layout<T1>::value,
                  "T1 must be standard layout");
    static_assert(std::is_standard_layout<T2>::value,
                  "T2 must be standard layout");
    typedef typename detail::VectorTraits<T1> traitsT1;
    typedef typename detail::VectorTraits<T2> traitsT2;

    pavector source_avector = new_pointer_avector(
        const_cast<double *>(
            reinterpret_cast<const double *>(source_vector.data())),
        source_vector.size() * traitsT1::length);
    pavector target_avector =
        new_pointer_avector(reinterpret_cast<double *>(target_vector.data()),
                            target_vector.size() * traitsT2::length);
    mvm_hmatrix_avector(alpha, h2trans, m_h.get(), source_avector,
                        target_avector);
  }
};

template <unsigned int D, typename Function,
          typename KernelHelper = detail::position_kernel_helper<D, Function>,
          typename Block = typename KernelHelper::Block>
detail::H2LibBlackBoxExpansions<D, Function, Block::RowsAtCompileTime,
                                Block::ColsAtCompileTime>
make_h2lib_black_box_expansion(size_t order, const Function &function) {
  return detail::H2LibBlackBoxExpansions<D, Function, Block::RowsAtCompileTime,
                                         Block::ColsAtCompileTime>(order,
                                                                   function);
}

template <typename Expansions, typename Kernel, typename RowParticlesType,
          typename ColParticlesType>
HLibMatrix make_h2lib_h_matrix(const RowParticlesType &row_particles,
                               const ColParticlesType &col_particles,
                               const Expansions &expansions,
                               const Kernel &kernel) {
  return HLibMatrix(row_particles, col_particles, expansions, kernel);
}

template <typename Expansions, typename Kernel, typename RowParticlesType,
          typename ColParticlesType>
H2LibMatrix make_h2lib_matrix(const RowParticlesType &row_particles,
                              const ColParticlesType &col_particles,
                              const Expansions &expansions,
                              const Kernel &kernel, const double eta) {
  return H2LibMatrix(row_particles, col_particles, expansions, kernel, eta);
}

} // namespace Aboria

#endif // HAVE_H2LIB

#endif
