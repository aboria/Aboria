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
    H2Lib_LR_Decomposition(const ph2matrix h2mat):
        tm(new_abseucl_truncmode()),tol(1e-8) {

        pclusterbasis rbcopy = clone_clusterbasis(h2mat->rb);
        pclusterbasis cbcopy = clone_clusterbasis(h2mat->cb);
        A = clone_h2matrix(h2mat,rbcopy,cbcopy);

        pclusteroperator Arwf = prepare_row_clusteroperator(A->rb, A->cb, tm);
        pclusteroperator Acwf = prepare_col_clusteroperator(A->rb, A->cb, tm);

        pblock broot = build_from_h2matrix_block(A);

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
    H2LibCholeskyDecomposition(const ph2matrix h2mat):
        tm(new_abseucl_truncmode()),tol(1e-8) {

        A = clone_h2matrix(h2mat,h2mat->rb,h2mat->cb);

        pblock broot = build_from_h2matrix_block(A);

        tolower_h2matrix(A);

        pclusteroperator Arwf = prepare_row_clusteroperator(A->rb, A->cb, tm);
        pclusteroperator Acwf = prepare_col_clusteroperator(A->rb, A->cb, tm);


        pclusterbasis Lrcb = build_from_cluster_clusterbasis(A->rb->t);
        pclusterbasis Lccb = build_from_cluster_clusterbasis(A->cb->t);

        L = build_from_block_lower_h2matrix(broot,Lrcb,Lccb);

        pclusteroperator Lrwf = prepare_row_clusteroperator(Lrcb, Lccb, tm);
        pclusteroperator Lcwf = prepare_col_clusteroperator(Lrcb, Lccb, tm);

        choldecomp_h2matrix(A, Arwf, Acwf, L, Lrwf, Lcwf, tm, tol);
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

    ph2matrix m_root_h2mat;

public:

    template <typename RowParticles, typename ColParticles, typename Expansions>
    H2LibMatrix(const RowParticles &row_particles, 
                const ColParticles &col_particles, 
                const Expansions& expansions) {

        //generate h2 matrix 
        LOG(2,"H2LibMatrix: creating h2 matrix using "<<row_particles.size()<<" row particles and "<<col_particles.size()<<" column particles");


        pclusterbasis row_clusterbasis = create_root_clusterbasis(row_particles,expansions);

        const bool row_equals_col = static_cast<const void*>(&row_particles) 
                                        == static_cast<const void*>(&col_particles);
        pclusterbasis col_clusterbasis;
        if (row_equals_col) {
            //col_clusterbasis = build_from_cluster_clusterbasis(row_clusterbasis->t);
            col_clusterbasis = clone_clusterbasis(row_clusterbasis);
        } else {
            col_clusterbasis = create_root_clusterbasis(col_particles,expansions);
        }


        m_root_h2mat = new_h2matrix(row_clusterbasis,col_clusterbasis);
        m_root_h2mat->rsons = row_clusterbasis->sons;
        if (m_root_h2mat->rsons == 0) m_root_h2mat->rsons = 1;
        m_root_h2mat->csons = col_clusterbasis->sons;
        if (m_root_h2mat->csons == 0) m_root_h2mat->csons = 1;
        m_root_h2mat->son = new ph2matrix[m_root_h2mat->rsons*m_root_h2mat->csons];
        for (int i = 0; i < m_root_h2mat->rsons*m_root_h2mat->csons; ++i) {
            m_root_h2mat->son[i] = nullptr;
        }

        //loop through h2 matrix and fill that out
        size_t j = 0;
        for (auto cj = col_particles.get_query().get_children(); 
                cj != false; ++cj,++j) {
            size_t i = 0;
            for (auto ci = row_particles.get_query().get_children(); 
                    ci != false; ++ci,++i) {
                ph2matrix child_h2 = generate_matrices(ci,cj,
                              row_clusterbasis->son[i],
                              col_clusterbasis->son[j],
                              row_particles,col_particles,
                              expansions);
                ref_h2matrix(m_root_h2mat->son + i+j*m_root_h2mat->rsons,child_h2);
            }
        }
        update_h2matrix(m_root_h2mat);
    }


    ~H2LibMatrix() {
        /*
        pclusterbasis row_clusterbasis = m_root_h2mat->rb;
        pclusterbasis col_clusterbasis = m_root_h2mat->cb;
        pcluster row_cluster = row_clusterbasis->t;
        pcluster col_cluster = col_clusterbasis->t;
        */
        del_h2matrix(m_root_h2mat);
        /*
        del_clusterbasis(row_clusterbasis);
        del_cluster(cluster);
        */
    }

private:

    template <typename Particles, typename Expansions>
    pclusterbasis create_root_clusterbasis(const Particles& particles,
                               const Expansions& expansions) {

        const size_t dim = Particles::dimension;
        pcluster t;
        std::vector<std::tuple<pclusterbasis,pcluster>> cb_children;
        uint *idx;
        size_t rank_sum = 0;
        size_t max_rank = std::numeric_limits<size_t>::min();

        // create cluster 
        size_t nsons = 0;
        
        // create cluster
        for (auto ci_child = particles.get_query().get_children(); 
                ci_child != false; ++ci_child) {
            ++nsons;
        }

        // TODO: do I need to set idx for non-leafs?
        t = new_cluster(0,nullptr,nsons,dim);
        auto bbox = detail::bbox<dim>(particles.get_min(),particles.get_max());
        // bounding box
        for (int i = 0; i < dim; ++i) {
            t->bmin[i] = particles.get_min()[i];
            t->bmax[i] = particles.get_max()[i];
        }
        cb_children.resize(nsons);

        // create child clusters, link them together and calculate descendents
        size_t i = 0;
        t->desc = 1;
        size_t np = 0;
        for (auto ci_child = particles.get_query().get_children(); 
                ci_child != false; ++ci_child,++i) {
            cb_children[i] = finalise_clusterbasis(ci_child,bbox,particles,expansions);
            rank_sum += std::get<0>(cb_children[i])->ktree;
            if (std::get<0>(cb_children[i])->kbranch > max_rank) {
                max_rank = std::get<0>(cb_children[i])->kbranch;
            }
            np += std::get<1>(cb_children[i])->size;
            t->son[i] = std::get<1>(cb_children[i]);
            t->desc += std::get<1>(cb_children[i])->desc;
        }
        t->size = np;

        // create clusterbasis
        pclusterbasis cb;

        cb = new_clusterbasis(t);

        for (int i = 0; i < nsons; ++i) {
            ref_clusterbasis(cb->son + i,std::get<0>(cb_children[i]));
        }

        cb->k = Expansions::ncheb;
        cb->ktree = rank_sum + cb->k;
        //cb->kbranch = max_rank < cb->k ? cb->k : max_rank;
        cb->kbranch = max_rank + cb->k;

        return cb;
    }

    template <typename ChildIterator, typename Particles, typename Expansions>
    std::tuple<pclusterbasis,pcluster>
    finalise_clusterbasis(const ChildIterator& ci, 
                            const detail::bbox<Particles::dimension>& parent_bbox,
                               const Particles& particles,
                               const Expansions& expansions) {

        const size_t index = particles.get_query().get_bucket_index(*ci); 
        const size_t dim = Particles::dimension;
        pcluster t;
        std::vector<std::tuple<pclusterbasis,pcluster>> cb_children;
        uint *idx;
        size_t rank_sum = 0;
        size_t max_rank = std::numeric_limits<size_t>::min();
        const auto& bbox = particles.get_query().get_bounds(ci);

        // create cluster 
        size_t nsons = 0;
        if (particles.get_query().is_leaf_node(*ci)) { // leaf node
            //get target_index
            const auto& prange = particles.get_query().get_bucket_particles(*ci);
            const size_t np = std::distance(prange.begin(),prange.end());
            idx = new uint[np];
            auto pit = prange.begin();
            for (int i = 0; i < np; ++i,++pit) {
                const size_t pi = &(get<typename Particles::position>(*pit))
                    - &get<typename Particles::position>(particles)[0];
                idx[i] = pi;
            }

            t = new_cluster(np,idx,0,dim);
            t->desc = 1;

        } else {

            // create cluster
            for (auto ci_child = particles.get_query().get_children(ci); 
                    ci_child != false; ++ci_child) {
                ++nsons;
            }

            // TODO: do I need to set idx for non-leafs?
            t = new_cluster(0,nullptr,nsons,dim);
            cb_children.resize(nsons);

            // create child clusters, link them together and calculate descendents
            size_t i = 0;
            t->desc = 1;
            size_t np = 0;
            for (auto ci_child = particles.get_query().get_children(ci); 
                    ci_child != false; ++ci_child,++i) {
                cb_children[i] = finalise_clusterbasis(ci_child,bbox,particles,expansions);
                rank_sum += std::get<0>(cb_children[i])->ktree;
                if (std::get<0>(cb_children[i])->kbranch > max_rank) {
                    max_rank = std::get<0>(cb_children[i])->kbranch;
                }
                np += std::get<1>(cb_children[i])->size;
                t->son[i] = std::get<1>(cb_children[i]);
                t->desc += std::get<1>(cb_children[i])->desc;
            }
            t->size = np;
        }

        // bounding box
        for (int i = 0; i < dim; ++i) {
            t->bmin[i] = bbox.bmin[i];
            t->bmax[i] = bbox.bmax[i];
        }

        // create clusterbasis
        pclusterbasis cb;

        if (particles.get_query().is_leaf_node(*ci)) { // leaf node
            cb = new_clusterbasis(t);
            cb->k = Expansions::ncheb;
            cb->ktree = cb->k + t->size;
            cb->kbranch = cb->ktree;
            // create L2P amatrix
            resize_amatrix(&cb->V,t->size,cb->k);
            expansions.L2P_amatrix(&cb->V,bbox,idx,t->size,particles);
        } else {
            cb = new_clusterbasis(t);
            cb->k = Expansions::ncheb;
            cb->ktree = rank_sum + cb->k;
            cb->kbranch = max_rank + cb->k;
            /*
            if (max_rank < cb->k) {
                max_rank = cb->k;
                cb->kbranch = cb->k;
            } else {
                cb->kbranch = max_rank;
            }
            */
        }

        for (int i = 0; i < nsons; ++i) {
            ref_clusterbasis(cb->son + i,std::get<0>(cb_children[i]));
        }

        resize_amatrix(&cb->E,cb->k,cb->k);
        expansions.L2L_amatrix(&cb->E,bbox,parent_bbox);

        return std::make_tuple(cb,t);
    }

    
    template <typename RowChildIterator, typename ColChildIterator,
             typename RowParticles, typename ColParticles, typename Expansions>
    ph2matrix generate_matrices(
            const RowChildIterator& ci, 
            const ColChildIterator& cj,
            pclusterbasis ci_cb,
            pclusterbasis cj_cb,
            const RowParticles &row_particles,
            const ColParticles &col_particles,
            const Expansions& expansions
            ) {

        const unsigned int dimension = RowParticles::dimension;
        const auto& ci_box = row_particles.get_query().get_bounds(ci);
        const auto& cj_box = col_particles.get_query().get_bounds(cj);

        size_t ci_index = row_particles.get_query().get_bucket_index(*ci);
        size_t cj_index = col_particles.get_query().get_bucket_index(*cj);

        LOG(3,"generate_matrices with buckets "<<ci_box<<" and "<<cj_box);

        ph2matrix h2mat = new_h2matrix(ci_cb,cj_cb);
        
        detail::theta_condition<dimension> theta(ci_box.bmin,ci_box.bmax);

        const bool ci_leaf = row_particles.get_query().is_leaf_node(*ci);
        const bool cj_leaf = row_particles.get_query().is_leaf_node(*cj);
        
        // each (i,j) is non-admissible leaf, admissible leaf, or otherwise
        if (!theta.check(cj_box.bmin,cj_box.bmax)) {
            // admissible leaf - M2L
            // create uniform 
            h2mat->u = new_uniform(ci_cb,cj_cb);
            expansions.M2L_amatrix(
                    &h2mat->u->S,ci_box,cj_box);
            h2mat->rsons = 0;
            h2mat->csons = 0;
            h2mat->son = NULL;
            h2mat->f = NULL;
        } else {
            // TODO: get_children should return empty range if leaf
            // non-leaf
            h2mat->rsons = ci_cb->sons;
            h2mat->csons = cj_cb->sons;
            if (h2mat->rsons == 0 && h2mat->csons == 0) {
                h2mat->rsons = 0;
                h2mat->csons = 0;
                h2mat->son = nullptr;
            } else {
                if (h2mat->rsons == 0) h2mat->rsons = 1;
                if (h2mat->csons == 0) h2mat->csons = 1;
                h2mat->son = new ph2matrix[h2mat->rsons*h2mat->csons];
                for (int i = 0; i < h2mat->rsons*h2mat->csons; ++i) {
                    h2mat->son[i] = nullptr;
                }
            }

            if (ci_leaf && cj_leaf) {
                //non-admissible leaf - P2P
                h2mat->f = new_amatrix(ci_cb->t->size,cj_cb->t->size);
                expansions.P2P_amatrix(
                        h2mat->f,
                        ci_cb->t->idx,ci_cb->t->size,
                        cj_cb->t->idx,cj_cb->t->size,
                        row_particles,col_particles);
            } else if (!ci_leaf && !cj_leaf) {
                size_t j = 0;
                for (auto cj_child = col_particles.get_query().get_children(cj); 
                        cj_child != false; ++cj_child,++j) {
                    size_t i = 0;
                    for (auto ci_child = row_particles.get_query().get_children(ci); 
                            ci_child != false; ++ci_child,++i) {
                        ph2matrix child_h2 = generate_matrices(ci_child,cj_child,
                                          ci_cb->son[i],cj_cb->son[j],
                                          row_particles,col_particles,
                                          expansions);
                        ref_h2matrix(h2mat->son + i+j*h2mat->rsons,child_h2);
                    }
                }
            } else if (!ci_leaf) {
                size_t i = 0;
                for (auto ci_child = row_particles.get_query().get_children(ci); 
                        ci_child != false; ++ci_child,++i) {
                    ph2matrix child_h2 = generate_matrices(ci_child,cj,
                                      ci_cb->son[i],cj_cb,
                                      row_particles,col_particles,
                                      expansions);
                    ref_h2matrix(h2mat->son + i,child_h2);
                }
            } else {
                size_t j = 0;
                for (auto cj_child = col_particles.get_query().get_children(cj); 
                        cj_child != false; ++cj_child,++j) {
                    ph2matrix child_h2 = generate_matrices(ci,cj_child,
                                      ci_cb,cj_cb->son[j],
                                      row_particles,col_particles,
                                      expansions);
                    ref_h2matrix(h2mat->son + j,child_h2);
                }
            }
            // link up rows and cols for uniform blocks
            /*
            std::vector<puniform> prev_row_u(h2mat->rsons,nullptr);
            for (int j = 0; j < h2mat->csons; ++j) {
                puniform prev_col_u = nullptr;
                for (int i = 0; i < h2mat->rsons; ++i) {
                    puniform child_u = h2mat->son[i + j*h2mat->rsons]->u;
                    if (child_u != nullptr) {
                        // update prev entries
                        child_u->rprev = prev_row_u[i];
                        child_u->cprev = prev_col_u;
                        // update next entries of prev u's
                        if (child_u->rprev != nullptr) {
                            child_u->rprev->rnext = child_u;
                        }
                        if (child_u->cprev != nullptr) {
                            child_u->cprev->cnext = child_u;
                        }
                    }
                }
            }
            */
        }
        update_h2matrix(h2mat);
        return h2mat;
    }

public:

    const ph2matrix get_ph2matrix() const {
        return m_root_h2mat;
    }

    H2Lib_LR_Decomposition lr() const {
        return H2Lib_LR_Decomposition(m_root_h2mat);
    }

    H2LibCholeskyDecomposition chol() const {
        return H2LibCholeskyDecomposition(m_root_h2mat);
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
        mvm_h2matrix_avector(alpha,h2trans,m_root_h2mat,source_avector,target_avector);
    }
};

    

template <typename Expansions, typename RowParticlesType, typename ColParticlesType>
H2LibMatrix make_h2lib_matrix(const RowParticlesType& row_particles, const ColParticlesType& col_particles, const Expansions& expansions) {
    return H2LibMatrix(row_particles,col_particles,expansions);
}


}

#endif

