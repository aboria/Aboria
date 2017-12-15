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
#undef I
}

namespace Aboria {

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
        pclusterbasis col_clusterbasis = create_root_clusterbasis(col_particles,expansions);

        m_root_h2mat = new_h2matrix(row_clusterbasis,col_clusterbasis);
        m_root_h2mat->rsons = row_clusterbasis->sons;
        if (m_root_h2mat->rsons == 0) m_root_h2mat->rsons = 1;
        m_root_h2mat->csons = col_clusterbasis->sons;
        if (m_root_h2mat->csons == 0) m_root_h2mat->csons = 1;
        m_root_h2mat->son = new ph2matrix[m_root_h2mat->rsons*m_root_h2mat->csons];

        //loop through h2 matrix and fill that out
        size_t j = 0;
        for (auto cj = col_particles.get_query().get_children(); 
                cj != false; ++cj,++j) {
            size_t i = 0;
            for (auto ci = row_particles.get_query().get_children(); 
                    ci != false; ++ci,++i) {
                ph2matrix child_h2 = generate_matrices(ci,cj,
                              row_clusterbasis,row_clusterbasis->son[i],
                              col_clusterbasis,col_clusterbasis->son[j],
                              row_particles,col_particles,
                              expansions);
                ref_h2matrix(m_root_h2mat->son + i+j*m_root_h2mat->rsons,child_h2);
            }
        }
        update_h2matrix(m_root_h2mat);
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
        cb_children.resize(nsons);

        // create child clusters, link them together and calculate descendents
        size_t i = 0;
        t->desc = 1;
        for (auto ci_child = particles.get_query().get_children(); 
                ci_child != false; ++ci_child,++i) {
            cb_children[i] = finalise_clusterbasis(ci_child,particles,expansions);
            rank_sum += std::get<0>(cb_children[i])->ktree;
            if (std::get<0>(cb_children[i])->kbranch > max_rank) {
                max_rank = std::get<0>(cb_children[i])->kbranch;
            }
            t->son[i] = std::get<1>(cb_children[i]);
            t->desc += std::get<1>(cb_children[i])->desc;
        }

        // bounding box
        for (int i = 0; i < dim; ++i) {
            t->bmin[i] = particles.get_min()[i];
            t->bmax[i] = particles.get_max()[i];
        }

        // create clusterbasis
        pclusterbasis cb;

        cb = new_clusterbasis(t);

        for (int i = 0; i < nsons; ++i) {
            ref_clusterbasis(cb->son + i,std::get<0>(cb_children[i]));
        }

        cb->k = Expansions::ncheb;
        cb->ktree = rank_sum + cb->k;
        cb->kbranch = max_rank + cb->k;

        return cb;
    }

    template <typename ChildIterator, typename Particles, typename Expansions>
    std::tuple<pclusterbasis,pcluster>
    finalise_clusterbasis(const ChildIterator& ci, 
                               const Particles& particles,
                               const Expansions& expansions) {

        const size_t index = particles.get_query()->get_bucket_index(*ci); 
        const size_t dim = Particles::dimension;
        pcluster t;
        std::vector<pclusterbasis> cb_children;
        uint *idx;
        size_t rank_sum = 0;
        size_t max_rank = std::numeric_limits<size_t>::min();

        // create cluster 
        size_t nsons = 0;
        if (particles.get_query()->is_leaf_node(*ci)) { // leaf node
            //get target_index
            const auto& prange = particles.get_query()->get_bucket_particles(*ci);
            const size_t np = std::distance(prange.begin(),prange.end());
            idx = new uint[np];
            auto pit = prange.begin();
            for (int i = 0; i < np; ++i,++pit) {
                const size_t pi = &(*get<Particles::position>(pit))
                    - &get<Particles::position>(particles)[0];
                idx[i] = pi;
            }

            t = new_cluster(np,idx,0,dim);

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
            for (auto ci_child = particles.get_query().get_children(ci); 
                    ci_child != false; ++ci_child,++i) {
                cb_children[i] = finalise_clusterbasis(ci_child,particles,expansions);
                rank_sum += cb_children[i]->ktree;
                if (cb_children[i]->kbranch > max_rank) {
                    max_rank = cb_children[i]->kbranch;
                }
                t->son[i] = cb_children[i]->t;
                t->desc += cb_children[i]->t->desc;
            }
        }

        // bounding box
        const auto& bbox = particles.get_query()->get_bounds(ci);
        for (int i = 0; i < dim; ++i) {
            t->bmin[i] = bbox.bmin[i];
            t->bmax[i] = bbox.bmax[i];
        }

        // create clusterbasis
        pclusterbasis cb;

        if (particles.get_query()->is_leaf_node(*ci)) { // leaf node
            cb = new_leaf_clusterbasis(t);
            cb->k = Expansions::ncheb;
            // create P2M amatrix
            // TODO: P2M or L2P here?
            cb->V = new_amatrix(cb->k,t->size);
            expansions.P2M_amatrix(cb->V,bbox,idx,particles);
        } else {
            cb = new_clusterbasis(t);
        }

        for (int i = 0; i < nsons; ++i) {
            ref_clusterbasis(cb->son + i,cb_children[i]);
        }

        cb->ktree = rank_sum + cb->k;
        cb->kbranch = max_rank + cb->k;

        cb->E = new_amatrix(cb->k,cb->k);
        expansions.L2L_amatrix(cb->E,bbox,ci.bounds);

        return std::make_tuple(cb,t);
    }

    
    template <typename RowChildIterator, typename ColChildIterator,
             typename RowParticles, typename ColParticles, typename Expansions>
    ph2matrix generate_matrices(
            const RowChildIterator& ci, 
            const ColChildIterator& cj,
            pclusterbasis parent_ci_cb,
            pclusterbasis ci_cb,
            pclusterbasis parent_cj_cb,
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
        if (ci_leaf && cj_leaf) {
            //non-admissible leaf - P2P
            h2mat->f = new_amatrix(ci_cb->t->size,cj_cb->t->size);
            expansions.P2P_amatrix(
                    h2mat->f,ci_cb->t->idx,ci_cb->t->size,cj_cb->t->idx,cj_cb->t->size,
                    row_particles,col_particles);

        } else if (!theta.check(cj_box.bmin,cj_box.bmax)) {
            // admissible leaf - M2L
            // create uniform 
            h2mat->u = new_uniform(ci_cb,cj_cb);
            expansions.M2L_amatrix(
                    h2mat->u,ci_box,cj_box);
        } else {
            // TODO: get_children should return empty range if leaf
            // non-leaf
            h2mat->rsons = ci_cb->sons;
            if (h2mat->rsons == 0) h2mat->rsons = 1;
            h2mat->csons = cj_cb->sons;
            if (h2mat->csons == 0) h2mat->csons = 1;
            h2mat->son = new ph2matrix[h2mat->rsons*h2mat->csons];
            
            if (!ci_leaf && !cj_leaf) {
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



    // target_vector += A*source_vector
    template <typename VectorTypeTarget, typename VectorTypeSource>
    void matrix_vector_multiply(VectorTypeTarget& target_vector, 
                          const VectorTypeSource& source_vector) const {
        avector source_avector;
        source_avector.v = source_vector.data();
        source_avector.dim = source_vector.size();
        avector target_avector;
        target_vector.v = target_vector.data();
        target_vector.dim = target_vector.size();
        mvm_h2matrix_avector(1,false,m_root_h2mat,source_vector,target_vector);
    }
};


template <typename Expansions, typename RowParticlesType, typename ColParticlesType>
H2LibMatrix make_h2lib_matrix(const RowParticlesType& row_particles, const ColParticlesType& col_particles, const Expansions& expansions) {
    return H2LibMatrix(row_particles,col_particles,expansions);
}


}

#endif

