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

extern "C" {
#include <h2matrix.h>
#include <h2update.h>
#include <h2arith.h>
#include <harith.h>
#include <hcoarsen.h>
#undef I
}

namespace Aboria {

class H2Lib_LR_Decomposition {
    ph2matrix A;
    ph2matrix L;
    ph2matrix R;
    ptruncmode tm;
    double tol;

public:
    H2Lib_LR_Decomposition(const ph2matrix h2mat, const pblock broot, const double tol):
        tm(new_releucl_truncmode()),tol(tol) {

        pclusterbasis rbcopy = clone_clusterbasis(h2mat->rb);
        pclusterbasis cbcopy = clone_clusterbasis(h2mat->cb);
        A = clone_h2matrix(h2mat,rbcopy,cbcopy);

        pclusteroperator Arwf = prepare_row_clusteroperator(A->rb, A->cb, tm);
        pclusteroperator Acwf = prepare_col_clusteroperator(A->rb, A->cb, tm);

        pclusterbasis Lrcb = build_from_cluster_clusterbasis(A->rb->t);
        pclusterbasis Lccb = build_from_cluster_clusterbasis(A->cb->t);
        pclusterbasis Rrcb = build_from_cluster_clusterbasis(A->rb->t);
        pclusterbasis Rccb = build_from_cluster_clusterbasis(A->cb->t);

        L = build_from_block_lower_h2matrix(broot,Lrcb,Lccb);
        R = build_from_block_upper_h2matrix(broot,Rrcb,Rccb);

        pclusteroperator Lrwf = prepare_row_clusteroperator(L->rb, L->cb, tm);
        pclusteroperator Lcwf = prepare_col_clusteroperator(L->rb, L->cb, tm);
        pclusteroperator Rrwf = prepare_row_clusteroperator(R->rb, R->cb, tm);
        pclusteroperator Rcwf = prepare_col_clusteroperator(R->rb, R->cb, tm);

        lrdecomp_h2matrix(A, Arwf, Acwf, L, Lrwf, Lcwf, R, Rrwf, Rcwf, tm, tol);
    }

    ~H2Lib_LR_Decomposition() {
        del_h2matrix(L);
        del_h2matrix(R);
        del_h2matrix(A);
    }

    template <typename VectorType>
    void solve(VectorType& vector) {
        pavector x = new_pointer_avector(
                           vector.data(),
                           vector.size());
        lrsolve_h2matrix_avector(L,R,x);
    }
};

class H2LibCholeskyDecomposition {
    ph2matrix A;
    ph2matrix L;
    ptruncmode tm;
    double tol;

public:
    H2LibCholeskyDecomposition(const ph2matrix h2mat, const pblock broot, const double tol):
        tm(new_abseucl_truncmode()),tol(tol) {

        pclusterbasis rbcopy = clone_clusterbasis(h2mat->rb);
        pclusterbasis cbcopy = clone_clusterbasis(h2mat->cb);
        A = clone_h2matrix(h2mat,rbcopy,cbcopy);

        pclusterbasis Lrcb = build_from_cluster_clusterbasis(A->rb->t);
        pclusterbasis Lccb = build_from_cluster_clusterbasis(A->cb->t);
        ASSERT_CUDA(A->rb->t == A->cb->t);

        L = build_from_block_lower_h2matrix(broot,Lrcb,Lccb);

        pclusteroperator Lrwf = prepare_row_clusteroperator(L->rb, L->cb, tm);
        pclusteroperator Lcwf = prepare_col_clusteroperator(L->rb, L->cb, tm);

        pclusteroperator Arwf = nullptr;
        pclusteroperator Acwf = nullptr;
        init_cholesky_h2matrix(A, &Arwf, &Acwf, tm);

        choldecomp_h2matrix(A, Arwf, Acwf, L, Lrwf, Lcwf, tm, tol);

        /*
        pclusteroperator rwfh2 = prepare_row_clusteroperator(h2mat->rb, h2mat->cb, tm);
        pclusteroperator cwfh2 = prepare_col_clusteroperator(h2mat->rb, h2mat->cb, tm);
        addmul_h2matrix(-1.0, L, true, L, h2mat, rwfh2, cwfh2, tm, tol);
        const double error = norm2_h2matrix(h2mat);
        std::cout << "error in factorisation is "<<norm2_h2matrix(h2mat)/error<<std::endl;
        */
    }

    ~H2LibCholeskyDecomposition() {
        del_h2matrix(L);
        del_h2matrix(A);
    }

    template <typename VectorType>
    void solve(VectorType& vector) {
        pavector x = new_pointer_avector(
                           vector.data(),
                           vector.size());
        cholsolve_h2matrix_avector(L,x);
    }
};

class HLib_LR_Decomposition {
    phmatrix A;
    ptruncmode tm;
    double tol;

public:
    HLib_LR_Decomposition(const phmatrix hmat, const double tol):
        tm(new_releucl_truncmode()),tol(tol) {
        A = clone_hmatrix(hmat);
        lrdecomp_hmatrix(A, tm, tol);
    }

    ~HLib_LR_Decomposition() {
        del_hmatrix(A);
    }

    template <typename VectorType>
    void solve(VectorType& vector) {
        pavector x = new_pointer_avector(
                           vector.data(),
                           vector.size());
        lrsolve_hmatrix_avector(false,A,x);
    }
};

class HLibCholeskyDecomposition {
    phmatrix A;
    ptruncmode tm;
    double tol;

public:
    HLibCholeskyDecomposition(const phmatrix hmat, const double tol):
        tm(new_abseucl_truncmode()),tol(tol) {

        A = clone_hmatrix(hmat);
        choldecomp_hmatrix(A, tm, tol);
    }

    ~HLibCholeskyDecomposition() {
        del_hmatrix(A);
    }

    template <typename VectorType>
    void solve(VectorType& vector) {
        pavector x = new_pointer_avector(
                           vector.data(),
                           vector.size());
        cholsolve_hmatrix_avector(A,x);
    }
};





class H2LibMatrix {

    std::unique_ptr<h2matrix,decltype(&del_h2matrix)> m_h2;
    std::unique_ptr<block,decltype(&del_block)> m_block;
    std::vector<uint> m_row_idx;
    std::vector<uint> m_col_idx;

template <typename Particles,typename Expansions>
static void
assemble_h2matrix_row_clusterbasis(pcclusterbasis rbc, uint rname, void *data)
{
    auto& data_cast = *static_cast<
        std::tuple<Expansions*,Particles*>*>(data);

    const Expansions& expansions = *std::get<0>(data_cast);
    const Particles& particles = *std::get<1>(data_cast);
    pclusterbasis rb = (pclusterbasis) rbc;

    const uint k = expansions.m_ncheb;
    if (rb->sons > 0) {
        resize_clusterbasis(rb,k);
        for (int i = 0; i < rb->sons; ++i) {
            expansions.L2L_amatrix(&rb->son[i]->E,rb->son[i]->t,rb->t);
        }
    } else {
        resize_amatrix(&rb->V,rb->t->size,k);
        expansions.L2P_amatrix(&rb->V,rb->t,rb->t->idx,rb->t->size,particles);
        rb->k = k;
        update_clusterbasis(rb);
    }
}


template <typename RowParticles,typename ColParticles,typename Expansions,typename Kernel>
static void
assemble_block_h2matrix(pcblock b, uint bname,
				    uint rname, uint cname, uint pardepth,
				    void *data)
{
    auto& data_cast = *static_cast<
        std::tuple<Expansions*,Kernel*,RowParticles*,ColParticles*,ph2matrix*>*>(data);

    const Expansions& expansions = *std::get<0>(data_cast);
    const Kernel& kernel = *std::get<1>(data_cast);
    const RowParticles& row_particles = *std::get<2>(data_cast);
    const ColParticles& col_particles = *std::get<3>(data_cast);
    ph2matrix* enum_h2 = std::get<3>(data_cast);
    ph2matrix h2 = enum_h2[bname];

    if (h2->u) {
        const uint kr = h2->u->rb->k;
        const uint kc = h2->u->cb->k;
        resize_amatrix(&h2->u->S,kr,kc);
        expansions.M2L_amatrix(&h2->u->S,h2->u->rb->t,h2->u->cb->t);

    } else if (h2->f) {
        detail::P2P_amatrix(h2->f,
                h2->rb->t->idx,h2->rb->t->size,
                h2->cb->t->idx,h2->cb->t->size,
                row_particles,col_particles,kernel);
    }

  }

public:

    template <typename RowParticles, typename ColParticles, typename Expansions, typename Kernel>
    H2LibMatrix(const RowParticles &row_particles, 
                const ColParticles &col_particles, 
                const Expansions& expansions,const Kernel& kernel, const double eta = 1.0):
            m_h2(nullptr,del_h2matrix),
            m_block(nullptr,del_block) {

        //generate h2 matrix 
        LOG(2,"H2LibMatrix: creating h2 matrix using "<<row_particles.size()<<" row particles and "<<col_particles.size()<<" column particles");

        const bool row_equals_col = static_cast<const void*>(&row_particles) 
                                        == static_cast<const void*>(&col_particles);

        //
        // Create row clusters and clusterbasis
        //
        m_row_idx.resize(row_particles.size());
        pcluster row_t = set_root_idx(m_row_idx.data(),row_particles);
        pclusterbasis row_cb = build_from_cluster_clusterbasis(row_t);
        auto data_row = std::make_tuple(&expansions,&row_particles);
        iterate_parallel_clusterbasis(row_cb, 0, max_pardepth, NULL,
                assemble_h2matrix_row_clusterbasis<RowParticles,Expansions>,
				&data_row);


        //
        // Create col clusters and clusterbasis
        //
        pcluster col_t;
        if (row_equals_col) {
            col_t = row_t;
        } else {
            m_col_idx.resize(col_particles.size());
            col_t = set_root_idx(m_col_idx.data(),col_particles);
        }
        pclusterbasis col_cb = build_from_cluster_clusterbasis(col_t);
        auto data_col = std::make_tuple(&expansions,&col_particles);
        iterate_parallel_clusterbasis(col_cb, 0, max_pardepth, NULL,
                assemble_h2matrix_row_clusterbasis<ColParticles,Expansions>,
				&data_col);

        // 
        // create h2 block
        //
        //  eta = 1.0;

        double eta_copy = eta;
        m_block = std::unique_ptr<block,decltype(&del_block)>(
                build_strict_block(row_t,col_t,&eta_copy,admissible_max_cluster),
                del_block);

        m_h2 = std::unique_ptr<h2matrix,decltype(&del_h2matrix)>(
                    build_from_block_h2matrix(m_block.get(),row_cb,col_cb),
                    del_h2matrix);
        ph2matrix* enum_h2mat = enumerate_h2matrix(m_h2.get());
        auto data_h2 = std::make_tuple(&expansions,&kernel,&row_particles,
                                    &col_particles,enum_h2mat);
        iterate_byrow_block(m_block.get(), 0, 0, 0, max_pardepth, NULL,
		      assemble_block_h2matrix<RowParticles,ColParticles,Expansions,Kernel>, 
              &data_h2);
        freemem(enum_h2mat);
    }

    void compress(const double tol) {
        //ptruncmode tm = new_releucl_truncmode();
        ptruncmode tm = new_blockreleucl_truncmode();
        recompress_inplace_h2matrix(m_h2.get(),tm,tol);
    }


private:

    template <typename Particles>
    pcluster set_root_idx(uint* idx, const Particles& particles) {
        uint *old_idx = idx;
        const uint dim = Particles::dimension;
        size_t sons = 0;
        for (auto ci_child = particles.get_query().get_children(); 
                ci_child != false; ++ci_child,++sons) {
        }
        pcluster t = new_cluster(0,old_idx,sons,dim);
        size_t i = 0;
        for (auto ci_child = particles.get_query().get_children(); 
                ci_child != false; ++ci_child,++i) {
            t->son[i] = set_idx(ci_child,idx,particles);
            idx += t->son[i]->size;
        }
        t->size = idx-old_idx;
        // bounding box
        for (int i = 0; i < dim; ++i) {
            t->bmin[i] = particles.get_min()[i];
            t->bmax[i] = particles.get_max()[i];
        }
        update_cluster(t);
        return t;
    }

    template <typename ChildIterator, typename Particles>
    pcluster set_idx(const ChildIterator ci, uint* idx,const Particles& particles) {
        uint *old_idx = idx;
        const uint dim = Particles::dimension;
        pcluster t;
        if (particles.get_query().is_leaf_node(*ci)) { // leaf node
            const auto& prange = particles.get_query().get_bucket_particles(*ci);
            auto pit = prange.begin();
            for (auto pit = prange.begin(); pit != prange.end(); ++pit,++idx) {
                const size_t pi = &(get<typename Particles::position>(*pit))
                    - &get<typename Particles::position>(particles)[0];
                *idx = pi;
                
            }
            if (idx-old_idx == 0) std::cout << "HAVE EMPTY LEAF" <<std::endl;
            t = new_cluster(idx-old_idx,old_idx,0,dim);
        } else {
            size_t sons = 0;
            for (auto ci_child = particles.get_query().get_children(ci); 
                      ci_child != false; ++ci_child,++sons) {
            }
            t = new_cluster(0,old_idx,sons,dim);
            size_t i = 0;
            for (auto ci_child = particles.get_query().get_children(ci); 
                      ci_child != false; ++ci_child,++i) {
                pcluster child_t = set_idx(ci_child,idx,particles);
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
        const auto& bbox = particles.get_query().get_bounds(ci);
        for (int i = 0; i < dim; ++i) {
            t->bmin[i] = bbox.bmin[i];
            t->bmax[i] = bbox.bmax[i];
        }
        update_cluster(t);
        return t;
    }


public:

    const ph2matrix get_ph2matrix() const {
        return m_h2.get();
    }

    H2Lib_LR_Decomposition lr(const double tol) const {
        return H2Lib_LR_Decomposition(get_ph2matrix(),m_block.get(),tol);
    }

    H2LibCholeskyDecomposition chol(const double tol) const {
        return H2LibCholeskyDecomposition(get_ph2matrix(),m_block.get(),tol);
    }

    // target_vector += alpha*A*source_vector or alpha*A'*source_vector
    template <typename VectorTypeTarget, typename VectorTypeSource>
    void matrix_vector_multiply(VectorTypeTarget& target_vector, 
                                const double alpha, const bool h2trans,
                                const VectorTypeSource& source_vector) const {
        pavector source_avector = new_pointer_avector(
                                    const_cast<double*>(source_vector.data()),
                                    source_vector.size());
        pavector target_avector = new_pointer_avector(
                                    target_vector.data(),
                                    target_vector.size());
        mvm_h2matrix_avector(alpha,h2trans,m_h2.get(),source_avector,target_avector);
    }
};

class HLibMatrix {

    std::unique_ptr<hmatrix,decltype(&del_hmatrix)> m_h;
    std::unique_ptr<block,decltype(&del_block)> m_block;
    std::vector<uint> m_row_idx;
    std::vector<uint> m_col_idx;

template <typename RowParticles,typename ColParticles,typename Expansions,typename Kernel>
static void
assemble_block_hmatrix(pcblock b, uint bname,
				    uint rname, uint cname, uint pardepth,
				    void *data)
{
    auto& data_cast = *static_cast<
        std::tuple<Expansions*,Kernel*,RowParticles*,ColParticles*,phmatrix*>*>(data);

    const Expansions& expansions = *std::get<0>(data_cast);
    const Kernel& kernel = *std::get<1>(data_cast);
    const RowParticles& row_particles = *std::get<2>(data_cast);
    const ColParticles& col_particles = *std::get<3>(data_cast);
    phmatrix* enum_h = std::get<3>(data_cast);
    phmatrix h = enum_h[bname];

    if (h->r) {
        expansions.L2P_amatrix(&h->r->A,
                h->rc,h->rc->idx,h->rc->size,
                row_particles);
        expansions.M2P_amatrix(&h->r->B,
                h->cc,h->cc->idx,h->cc->size,
                col_particles);
    } else if (h->f) {
        detail::P2P_amatrix(h->f,
                h->rc->idx,h->rc->size,
                h->cc->idx,h->cc->size,
                row_particles,col_particles,kernel);
    }

}

static void
truncate_block_hmatrix(pcblock b, uint bname,
				    uint rname, uint cname, uint pardepth,
				    void *data)
{
    auto& data_cast = *static_cast<
        std::tuple<double*,phmatrix*>*>(data);

    const double& tol= *std::get<0>(data_cast);
    phmatrix* enum_h = std::get<1>(data_cast);
    phmatrix h = enum_h[bname];

    if (h->r) {
        trunc_rkmatrix(NULL, tol, h->r);
    }

  }

public:

    template <typename RowParticles, typename ColParticles, typename Expansions, typename Kernel>
    HLibMatrix(const RowParticles &row_particles, 
                const ColParticles &col_particles, 
                const Expansions& expansions,
                const Kernel& kernel):
            m_h(nullptr,del_hmatrix),
            m_block(nullptr,del_block) {

        //generate h2 matrix 
        LOG(2,"H2LibMatrix: creating h matrix using "<<row_particles.size()<<" row particles and "<<col_particles.size()<<" column particles");

        const bool row_equals_col = static_cast<const void*>(&row_particles) 
                                        == static_cast<const void*>(&col_particles);

        //
        // Create row clusters and clusterbasis
        //
        m_row_idx.resize(row_particles.size());
        pcluster row_t = set_root_idx(m_row_idx.data(),row_particles);
        
        //
        // Create col clusters and clusterbasis
        //
        pcluster col_t;
        if (row_equals_col) {
            col_t = row_t;
        } else {
            m_col_idx.resize(col_particles.size());
            col_t = set_root_idx(m_col_idx.data(),col_particles);
        }

        // 
        // create h block
        //
        //  eta = 1.0;

        double eta = 1.0;
        m_block = std::unique_ptr<block,decltype(&del_block)>(
                build_nonstrict_block(row_t,col_t,&eta,admissible_max_cluster),
                del_block);

        m_h = std::unique_ptr<hmatrix,decltype(&del_hmatrix)>(
                    build_from_block_hmatrix(m_block.get(),expansions.m_ncheb),
                    del_hmatrix);
        phmatrix* enum_hmat = enumerate_hmatrix(m_block.get(),m_h.get());
        auto data_h = std::make_tuple(&expansions,&kernel,&row_particles,
                                    &col_particles,enum_hmat);
        iterate_byrow_block(m_block.get(), 0, 0, 0, max_pardepth, NULL,
		      assemble_block_hmatrix<RowParticles,ColParticles,Expansions,Kernel>, 
              &data_h);
        freemem(enum_hmat);
    }


private:

    template <typename Particles>
    pcluster set_root_idx(uint* idx, const Particles& particles) {
        uint *old_idx = idx;
        const uint dim = Particles::dimension;
        size_t sons = 0;
        for (auto ci_child = particles.get_query().get_children(); 
                ci_child != false; ++ci_child,++sons) {
        }
        pcluster t = new_cluster(0,old_idx,sons,dim);
        size_t i = 0;
        for (auto ci_child = particles.get_query().get_children(); 
                ci_child != false; ++ci_child,++i) {
            t->son[i] = set_idx(ci_child,idx,particles);
            idx += t->son[i]->size;
        }
        t->size = idx-old_idx;
        // bounding box
        for (int i = 0; i < dim; ++i) {
            t->bmin[i] = particles.get_min()[i];
            t->bmax[i] = particles.get_max()[i];
        }
        update_cluster(t);
        return t;
    }

    template <typename ChildIterator, typename Particles>
    pcluster set_idx(const ChildIterator ci, uint* idx,const Particles& particles) {
        uint *old_idx = idx;
        const uint dim = Particles::dimension;
        pcluster t;
        if (particles.get_query().is_leaf_node(*ci)) { // leaf node
            const auto& prange = particles.get_query().get_bucket_particles(*ci);
            auto pit = prange.begin();
            for (auto pit = prange.begin(); pit != prange.end(); ++pit,++idx) {
                const size_t pi = &(get<typename Particles::position>(*pit))
                    - &get<typename Particles::position>(particles)[0];
                *idx = pi;
                
            }
            if (idx-old_idx == 0) std::cout << "HAVE EMPTY LEAF" <<std::endl;
            t = new_cluster(idx-old_idx,old_idx,0,dim);
        } else {
            size_t sons = 0;
            for (auto ci_child = particles.get_query().get_children(ci); 
                      ci_child != false; ++ci_child,++sons) {
            }
            t = new_cluster(0,old_idx,sons,dim);
            size_t i = 0;
            for (auto ci_child = particles.get_query().get_children(ci); 
                      ci_child != false; ++ci_child,++i) {
                pcluster child_t = set_idx(ci_child,idx,particles);
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
        const auto& bbox = particles.get_query().get_bounds(ci);
        for (int i = 0; i < dim; ++i) {
            t->bmin[i] = bbox.bmin[i];
            t->bmax[i] = bbox.bmax[i];
        }
        update_cluster(t);
        return t;
    }
    
public:

    void compress(const double tol) {
        //ptruncmode tm = new_releucl_truncmode();

        phmatrix* enum_hmat = enumerate_hmatrix(m_block.get(),m_h.get());
        auto data_h = std::make_tuple(&tol,enum_hmat);
        iterate_byrow_block(m_block.get(), 0, 0, 0, max_pardepth, NULL,
		      truncate_block_hmatrix, &data_h);
        freemem(enum_hmat);

        coarsen_hmatrix(m_h.get(),NULL,tol,true);
    }

    const phmatrix get_phmatrix() const {
        return m_h.get();
    }

    HLib_LR_Decomposition lr(const double tol) const {
        return HLib_LR_Decomposition(get_phmatrix(),tol);
    }

    HLibCholeskyDecomposition chol(const double tol) const {
        return HLibCholeskyDecomposition(get_phmatrix(),tol);
    }

    // target_vector += alpha*A*source_vector or alpha*A'*source_vector
    template <typename VectorTypeTarget, typename VectorTypeSource>
    void matrix_vector_multiply(VectorTypeTarget& target_vector, 
                                const double alpha, const bool h2trans,
                                const VectorTypeSource& source_vector) const {
        pavector source_avector = new_pointer_avector(
                                    const_cast<double*>(source_vector.data()),
                                    source_vector.size());
        pavector target_avector = new_pointer_avector(
                                    target_vector.data(),
                                    target_vector.size());
        mvm_hmatrix_avector(alpha,h2trans,m_h.get(),source_avector,target_avector);
    }
};

template <unsigned int D, typename Function> 
detail::H2LibBlackBoxExpansions<D,Function> make_h2lib_black_box_expansion(size_t order, const Function& function) {
    return detail::H2LibBlackBoxExpansions<D,Function>(order,function);
}

template <typename Expansions, typename Kernel, typename RowParticlesType, typename ColParticlesType>
HLibMatrix make_h2lib_h_matrix(const RowParticlesType& row_particles, const ColParticlesType& col_particles, const Expansions& expansions, const Kernel& kernel) {
    return HLibMatrix(row_particles,col_particles,expansions,kernel);
}
    

template <typename Expansions, typename Kernel, typename RowParticlesType, typename ColParticlesType>
H2LibMatrix make_h2lib_matrix(const RowParticlesType& row_particles, const ColParticlesType& col_particles, const Expansions& expansions, const Kernel& kernel) {
    return H2LibMatrix(row_particles,col_particles,expansions,kernel);
}


}

#endif

