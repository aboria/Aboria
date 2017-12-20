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
    H2Lib_LR_Decomposition(const ph2matrix h2mat, const pblock broot):
        tm(new_releucl_truncmode()),tol(1e-10) {

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
    H2LibCholeskyDecomposition(const ph2matrix h2mat, const pblock broot):
        tm(new_abseucl_truncmode()),tol(1e-8) {

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


template <typename RowParticles,typename ColParticles,typename Expansions>
static void
assemble_block_h2matrix(pcblock b, uint bname,
				    uint rname, uint cname, uint pardepth,
				    void *data)
{
    auto& data_cast = *static_cast<
        std::tuple<Expansions*,RowParticles*,ColParticles*,ph2matrix*>*>(data);

    const Expansions& expansions = *std::get<0>(data_cast);
    const RowParticles& row_particles = *std::get<1>(data_cast);
    const ColParticles& col_particles = *std::get<2>(data_cast);
    ph2matrix* enum_h2 = std::get<3>(data_cast);
    ph2matrix h2 = enum_h2[bname];

    if (h2->u) {
        const uint kr = h2->u->rb->k;
        const uint kc = h2->u->cb->k;
        resize_amatrix(&h2->u->S,kr,kc);
        expansions.M2L_amatrix(&h2->u->S,h2->u->rb->t,h2->u->cb->t);
    } else if (h2->f) {
        expansions.P2P_amatrix(h2->f,
                h2->rb->t->idx,h2->rb->t->size,
                h2->cb->t->idx,h2->cb->t->size,
                row_particles,col_particles);
    }

  }

public:

    template <typename RowParticles, typename ColParticles, typename Expansions>
    H2LibMatrix(const RowParticles &row_particles, 
                const ColParticles &col_particles, 
                const Expansions& expansions):
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
        m_block = std::unique_ptr<block,decltype(&del_block)>(
                new_block(row_t,col_t,false,row_t->sons,col_t->sons),
                del_block);
        size_t j = 0;
        for (auto cj = col_particles.get_query().get_children(); 
                cj != false; ++cj,++j) {
            size_t i = 0;
            for (auto ci = row_particles.get_query().get_children(); 
                    ci != false; ++ci,++i) {
                pblock child_block = generate_blocks(ci,cj,
                              row_t->son[i],
                              col_t->son[j],
                              row_particles,col_particles);
                m_block->son[i + j*row_t->sons] = child_block;
            }
        }
        update_block(m_block.get());

        m_h2 = std::unique_ptr<h2matrix,decltype(&del_h2matrix)>(
                    build_from_block_h2matrix(m_block.get(),row_cb,col_cb),
                    del_h2matrix);
        ph2matrix* enum_h2mat = enumerate_h2matrix(m_h2.get());
        auto data_h2 = std::make_tuple(&expansions,&row_particles,
                                    &col_particles,enum_h2mat);
        iterate_byrow_block(m_block.get(), 0, 0, 0, max_pardepth, NULL,
		      assemble_block_h2matrix<RowParticles,ColParticles,Expansions>, 
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
                t->son[i] = child_t;
                idx += child_t->size;
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


    template <typename RowChildIterator, typename ColChildIterator,
             typename RowParticles, typename ColParticles>
    pblock generate_blocks(
            const RowChildIterator& ci, 
            const ColChildIterator& cj,
            pcluster ci_t,
            pcluster cj_t,
            const RowParticles &row_particles,
            const ColParticles &col_particles) {

        pblock block;

        const unsigned int dimension = RowParticles::dimension;
        const auto& ci_box = row_particles.get_query().get_bounds(ci);
        const auto& cj_box = col_particles.get_query().get_bounds(cj);

        LOG(3,"generate_blocks with buckets "<<ci_box<<" and "<<cj_box);
        
        detail::theta_condition<dimension> theta(ci_box.bmin,ci_box.bmax);

        const bool ci_leaf = row_particles.get_query().is_leaf_node(*ci);
        const bool cj_leaf = col_particles.get_query().is_leaf_node(*cj);

        ASSERT_CUDA(ci_leaf ? ci_t->sons == 0 : ci_t->sons > 0);
        ASSERT_CUDA(cj_leaf ? cj_t->sons == 0 : cj_t->sons > 0);

        // each (i,j) is non-admissible leaf, admissible leaf, or otherwise
        if (!theta.check(cj_box.bmin,cj_box.bmax)) {
            block = new_block(ci_t,cj_t,true,0,0);
        } else if (ci_leaf && cj_leaf) {
            block = new_block(ci_t,cj_t,false,0,0);
        } else if (!ci_leaf && !cj_leaf) {
            block = new_block(ci_t,cj_t,false,ci_t->sons,cj_t->sons);
            size_t j = 0;
            for (auto cj_child = col_particles.get_query().get_children(cj); 
                    cj_child != false; ++cj_child,++j) {
                size_t i = 0;
                for (auto ci_child = row_particles.get_query().get_children(ci); 
                        ci_child != false; ++ci_child,++i) {
                    pblock child_block = generate_blocks(
                            ci_child,cj_child,
                            ci_t->son[i],cj_t->son[j],
                            row_particles,col_particles);
                    block->son[i + j*ci_t->sons] = child_block;
                }
            }
        } else if (!ci_leaf) {
            block = new_block(ci_t,cj_t,false,ci_t->sons,1);
            size_t i = 0;
            for (auto ci_child = row_particles.get_query().get_children(ci); 
                    ci_child != false; ++ci_child,++i) {
                pblock child_block = generate_blocks(
                        ci_child,cj,
                        ci_t->son[i],cj_t,
                        row_particles,col_particles);
                block->son[i] = child_block;
            }
        } else {
            block = new_block(ci_t,cj_t,false,1,cj_t->sons);
            size_t j = 0;
            for (auto cj_child = col_particles.get_query().get_children(cj); 
                    cj_child != false; ++cj_child,++j) {
                pblock child_block = generate_blocks(
                        ci,cj_child,
                        ci_t,cj_t->son[j],
                        row_particles,col_particles);
                block->son[j] = child_block;
            }
        }
        update_block(block);
        return block;
    }


public:

    const ph2matrix get_ph2matrix() const {
        return m_h2.get();
    }

    H2Lib_LR_Decomposition lr() const {
        return H2Lib_LR_Decomposition(get_ph2matrix(),m_block.get());
    }

    H2LibCholeskyDecomposition chol() const {
        return H2LibCholeskyDecomposition(get_ph2matrix(),m_block.get());
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

    

template <typename Expansions, typename RowParticlesType, typename ColParticlesType>
H2LibMatrix make_h2lib_matrix(const RowParticlesType& row_particles, const ColParticlesType& col_particles, const Expansions& expansions) {
    return H2LibMatrix(row_particles,col_particles,expansions);
}


}

#endif

