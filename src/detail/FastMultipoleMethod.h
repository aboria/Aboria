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

#include "detail/SpatialUtil.h"
#include "detail/Kernels.h"
#include "detail/Chebyshev.h"
#include "NeighbourSearchBase.h"
#include "Traits.h"
#include "CudaInclude.h"
#include "Vector.h"
#include "Get.h"
#include "Log.h"
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

    template <unsigned int D, unsigned int N> 
    struct MultiquadricExpansions {
        typedef detail::bbox<D> box_type;
        static const size_t ncheb = std::pow(N,D); 
        typedef std::array<double,ncheb> expansion_type;
#ifdef HAVE_EIGEN
        typedef Eigen::Matrix<double,ncheb,ncheb> matrix_type;
        typedef Eigen::Matrix<double,ncheb,Eigen::Dynamic> p2m_matrix_type;
        typedef Eigen::Matrix<double,Eigen::Dynamic,ncheb> m2p_matrix_type;
        typedef Eigen::Matrix<double,Eigen::Dynamic,1> dynamic_vector_type;
#endif
        typedef Vector<double,D> double_d;
        typedef Vector<int,D> int_d;
        static const unsigned int dimension = D;
        const double m_c2;

        MultiquadricExpansions(const double c):m_c2(c*c) 
        {}

        static void P2M(expansion_type& accum, 
                 const box_type& box, 
                 const double_d& position,
                 const double& source ) {

        }

        static void M2M(expansion_type& accum, 
                 const box_type& target_box, 
                 const box_type& source_box, 
                 const expansion_type& source) {

        }

        void M2L(expansion_type& accum, 
                 const box_type& target_box, 
                 const box_type& source_box, 
                 const expansion_type& source) {

          
        }

        static void L2L(expansion_type& accum, 
                 const box_type& target_box, 
                 const box_type& source_box, 
                 const expansion_type& source) {
         
        }

        static double L2P(const double_d& p,
                   const box_type& box, 
                   const expansion_type& source) {
            return 0.0;
        }

    };


    template <size_t D, size_t N, typename Function, size_t BlockRows, size_t BlockCols> 
    struct BlackBoxExpansions {
        typedef detail::bbox<D> box_type;
        static constexpr size_t ncheb = ipow(N,D); 

        static const size_t block_rows = BlockRows;
        static const size_t block_cols = BlockCols;
        
        typedef typename std::conditional<BlockRows==1,
                double,
                Vector<double,BlockRows>>::type row_scalar_type;
        typedef typename std::conditional<BlockCols==1,
                double,
                Vector<double,BlockCols>>::type col_scalar_type;
        typedef std::array<col_scalar_type,ncheb> m_expansion_type;
        typedef std::array<row_scalar_type,ncheb> l_expansion_type;

#ifdef HAVE_EIGEN
        typedef Eigen::Matrix<double,ncheb*BlockRows,ncheb*BlockCols> m2l_matrix_type;
        typedef Eigen::Matrix<double,ncheb*BlockRows,ncheb*BlockRows> l2l_matrix_type;
        typedef Eigen::Matrix<double,ncheb*BlockCols,ncheb*BlockCols> m2m_matrix_type;
        typedef Eigen::Matrix<double,ncheb*BlockCols,Eigen::Dynamic> p2m_matrix_type;
        typedef Eigen::Matrix<double,Eigen::Dynamic,ncheb*BlockRows> l2p_matrix_type;
        typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> p2p_matrix_type;
        typedef Eigen::Matrix<double,Eigen::Dynamic,1> p_vector_type;
        typedef Eigen::Matrix<double,ncheb*BlockCols,1> m_vector_type;
        typedef Eigen::Matrix<double,ncheb*BlockRows,1> l_vector_type;
        typedef Eigen::Matrix<double,BlockCols,BlockCols> blockcol_type;
        typedef Eigen::Matrix<double,BlockRows,BlockRows> blockrow_type;
#endif

        typedef Vector<double,D> double_d;
        typedef Vector<int,D> int_d;
        static const unsigned int dimension = D;
        static const unsigned int number_of_nodes_in_each_direction = N;
        Function m_K;
        std::array<double_d,ncheb> m_cheb_points; 

        BlackBoxExpansions(const Function &K):m_K(K) 
        {
            //precalculate cheb_points
            lattice_iterator<dimension> mi(int_d::Constant(0),int_d::Constant(N));
            for (int i=0; i<ncheb; ++i,++mi) {
                m_cheb_points[i] = detail::chebyshev_node_nd(*mi,N);
            }
        }


        static void P2M(m_expansion_type& accum,
                 const box_type& box,
                 const double_d& position,
                 const col_scalar_type& source) {

            detail::ChebyshevRnSingle<D,N> cheb_rn(position,box);
            lattice_iterator<dimension> mj(int_d::Constant(0),int_d::Constant(N));
            for (int j=0; j<ncheb; ++j,++mj) {
                //std::cout << "accumulating P2M from "<<position<<" to node "<<*mj<<" with Rn = "<<cheb_rn(*mj)<<std::endl;
                accum[j] += cheb_rn(*mj)*source;
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
                //std::cout << "accumulating P2M from "<<position<<" to node "<<*mj<<" with Rn = "<<cheb_rn(*mj)<<std::endl;
                accum[j] += cheb_rn(*mj)*source;
            }
        }
        */


        template <typename ParticlesType>
        static void P2M_matrix(p2m_matrix_type& matrix, 
                    const box_type& box,
                    const std::vector<size_t>& indicies,
                    const ParticlesType& particles) {
            typedef typename ParticlesType::position position;
            matrix.resize(ncheb*BlockCols,indicies.size()*BlockCols);
            for (int i = 0; i < indicies.size(); ++i) {
                const double_d& p = get<position>(particles)[indicies[i]];
                detail::ChebyshevRnSingle<D,N> cheb_rn(p,box);
                lattice_iterator<dimension> mj(int_d::Constant(0),int_d::Constant(N));
                for (int j=0; j<ncheb; ++j,++mj) {
                    // check ij
                    matrix.block<BlockCols,BlockCols>(i*BlockCols,j*BlockCols) = 
                                                        cheb_rn(*mj)*blockcol_type::Identity();
                }
                
            }
        }
#endif

        void M2M(m_expansion_type& accum, 
                 const box_type& target_box, 
                 const box_type& source_box, 
                 const m_expansion_type& source) const {

            for (int j=0; j<ncheb; ++j) {
                const double_d& pj_unit_box = m_cheb_points[j];
                const double_d pj = 0.5*(pj_unit_box+1)*(source_box.bmax-source_box.bmin) 
                    + source_box.bmin;
                detail::ChebyshevRnSingle<D,N> cheb_rn(pj,target_box);

                lattice_iterator<D> mi(int_d::Constant(0),int_d::Constant(N));
                for (int i=0; i<ncheb; ++i,++mi) {
                    accum[i] += cheb_rn(*mi)*source[j];
                }
            }
        }

#ifdef HAVE_EIGEN
        void M2M_matrix(m2m_matrix_type& matrix, 
                 const box_type& target_box, 
                 const box_type& source_box) const {
            for (int j=0; j<ncheb; ++j) {
                const double_d& pj_unit_box = m_cheb_points[j];
                const double_d pj = 0.5*(pj_unit_box+1)*(source_box.bmax-source_box.bmin) 
                    + source_box.bmin;
                detail::ChebyshevRnSingle<D,N> cheb_rn(pj,target_box);

                lattice_iterator<D> mi(int_d::Constant(0),int_d::Constant(N));
                for (int i=0; i<ncheb; ++i,++mi) {
                    matrix.block<BlockCols,BlockCols>(i*BlockCols,j*BlockCols) = 
                                                        cheb_rn(*mi)*blockcol_type::Identity();
                }
            }
        }
#endif

        void M2L(l_expansion_type& accum, 
                 const box_type& target_box, 
                 const box_type& source_box, 
                 const m_expansion_type& source) const {

            for (int i=0; i<ncheb; ++i) {
                const double_d& pi_unit_box = m_cheb_points[i];
                const double_d pi = 0.5*(pi_unit_box+1)*(target_box.bmax-target_box.bmin) 
                                                                    + target_box.bmin;

                for (int j=0; j<ncheb; ++j) {
                    const double_d& pj_unit_box = m_cheb_points[j];
                    const double_d pj = 0.5*(pj_unit_box+1)*(source_box.bmax-source_box.bmin) 
                                                                    + source_box.bmin;
                    accum[i] += m_K(pi,pj)*source[j];
                }
            }
        }

#ifdef HAVE_EIGEN
        void M2L_matrix(m2l_matrix_type& matrix, 
                 const box_type& target_box, 
                 const box_type& source_box) const {
            for (int i=0; i<ncheb; ++i) {
                const double_d& pi_unit_box = m_cheb_points[i];
                const double_d pi = 0.5*(pi_unit_box+1)*(target_box.bmax-target_box.bmin) 
                                                                    + target_box.bmin;
                for (int j=0; j<ncheb; ++j) {
                    const double_d& pj_unit_box = m_cheb_points[j];
                    const double_d pj = 0.5*(pj_unit_box+1)*(source_box.bmax-source_box.bmin) 
                                                                    + source_box.bmin;
                    matrix.block<BlockRows,BlockCols>(i*BlockRows,j*BlockCols) = m_K(pi,pj);
                }
            }


        }
#endif



        void L2L(l_expansion_type& accum, 
                 const box_type& target_box, 
                 const box_type& source_box, 
                 const l_expansion_type& source) const {
            //M2M(accum,target_box,source_box,source);
            for (int i=0; i<ncheb; ++i) {
                const double_d& pi_unit_box = m_cheb_points[i];
                const double_d pi = 0.5*(pi_unit_box+1)*(target_box.bmax-target_box.bmin) 
                    + target_box.bmin;
                detail::ChebyshevRnSingle<D,N> cheb_rn(pi,source_box);

                lattice_iterator<D> mj(int_d::Constant(0),int_d::Constant(N));
                for (int j=0; j<ncheb; ++j,++mj) {
                    accum[i] += cheb_rn(*mj)*source[j];
                }
            }
        }


#ifdef HAVE_EIGEN
        void L2L_matrix(l2l_matrix_type& matrix, 
                 const box_type& target_box, 
                 const box_type& source_box) const {
            for (int i=0; i<ncheb; ++i) {
                const double_d& pi_unit_box = m_cheb_points[i];
                const double_d pi = 0.5*(pi_unit_box+1)*(target_box.bmax-target_box.bmin) 
                    + target_box.bmin;
                detail::ChebyshevRnSingle<D,N> cheb_rn(pi,source_box);

                lattice_iterator<D> mj(int_d::Constant(0),int_d::Constant(N));
                for (int j=0; j<ncheb; ++j,++mj) {
                    matrix.block<BlockRows,BlockRows>(i*BlockRows,j*BlockRows) = 
                                                    cheb_rn(*mj)*blockrow_type::Identity();
                }
            }
        }
#endif


        static row_scalar_type L2P(const double_d& p,
                   const box_type& box, 
                   const l_expansion_type& source) {
            row_scalar_type sum = detail::VectorTraits<row_scalar_type>::Zero();
            detail::ChebyshevRnSingle<D,N> cheb_rn(p,box);
            lattice_iterator<dimension> mj(int_d::Constant(0),int_d::Constant(N));
            for (int j=0; j<ncheb; ++j,++mj) {
                sum += cheb_rn(*mj)*source[j];
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
        static void L2P_matrix(l2p_matrix_type& matrix, 
                    const box_type& box,
                    const std::vector<size_t>& indicies,
                    const ParticlesType& particles) {
            typedef typename ParticlesType::position position;
            matrix.resize(indicies.size()*BlockRows,ncheb*BlockRows);
            for (int i = 0; i < indicies.size(); ++i) {
                const double_d& p = get<position>(particles)[indicies[i]];
                detail::ChebyshevRnSingle<D,N> cheb_rn(p,box);
                lattice_iterator<dimension> mj(int_d::Constant(0),int_d::Constant(N));
                for (int j=0; j<ncheb; ++j,++mj) {
                    matrix.block<BlockRows,BlockRows>(i*BlockRows,j*BlockRows) = 
                                                cheb_rn(*mj)*blockrow_type::Identity();
                }
                
            }
        }

        
#endif


    };


#ifdef HAVE_H2LIB
    template <size_t D, typename Function, size_t BlockRows, size_t BlockCols> 
    struct H2LibBlackBoxExpansions {
        typedef detail::bbox<D> box_type;

        static const size_t block_rows = BlockRows;
        static const size_t block_cols = BlockCols;

        typedef Vector<double,BlockRows> row_scalar_type;
        typedef Vector<double,BlockCols> col_scalar_type;
#ifndef HAVE_EIGEN
        static_assert(BlockRows == BlockCols == 1,"Need to define HAVE_EIGEN to use matrix-valued kernel");
#endif

        typedef Vector<double,D> double_d;
        typedef Vector<int,D> int_d;
        static const unsigned int dimension = D;
        Function m_K;
        std::vector<double_d> m_cheb_points; 
        size_t m_order;
        size_t m_ncheb;

        H2LibBlackBoxExpansions(const size_t order, const Function &K):
            m_K(K),m_order(order),m_ncheb(std::pow(m_order,D))
        {
            m_cheb_points.resize(m_ncheb);
            //precalculate cheb_points
            lattice_iterator<dimension> mi(int_d::Constant(0),int_d::Constant(m_order));
            for (int i=0; i<m_cheb_points.size(); ++i,++mi) {
                m_cheb_points[i] = detail::chebyshev_node_nd(*mi,m_order);
            }
        }

        template <typename ParticlesType>
        void P2M_trans_amatrix(pamatrix matrix, 
                    const pccluster t,
                    const uint* indicies,
                    const uint indicies_size,
                    const ParticlesType& particles) const {
            typedef typename ParticlesType::position position;
            box_type box;
            for (int i = 0; i < D; ++i) {
                box.bmin[i] = t->bmin[i];
                box.bmax[i] = t->bmax[i];
            }
            ASSERT_CUDA(matrix->rows == m_ncheb);
            ASSERT_CUDA(matrix->cols == indicies_size);
            //resize_amatrix(matrix,m_ncheb,indicies_size);
            clear_amatrix(matrix);
            detail::ChebyshevRn<D> cheb_rn(m_order,box);
            for (int i = 0; i < indicies_size; ++i) {
                const double_d& p = get<position>(particles)[indicies[i]];
                cheb_rn.set_position(p);
                lattice_iterator<dimension> mj(int_d::Constant(0),int_d::Constant(m_order));
                for (int j=0; j<m_ncheb; ++j,++mj) {
                    const double tmp = cheb_rn(*mj);
                    for (int ii = 0; ii < BlockCols; ++ii) {
                        setentry_amatrix(matrix,j*BlockCols+ii,i*BlockCols+ii,tmp);
                    }
                }
            }
        }

        void M2M_amatrix(pamatrix matrix, 
                 pccluster target_t, 
                 pccluster source_t) const {
            //resize_amatrix(matrix,m_ncheb,m_ncheb);
            ASSERT_CUDA(matrix->rows == m_ncheb);
            ASSERT_CUDA(matrix->cols == m_ncheb);
            clear_amatrix(matrix);
            box_type target_box,source_box;
            for (int i = 0; i < D; ++i) {
                target_box.bmin[i] = target_t->bmin[i];
                target_box.bmax[i] = target_t->bmax[i];
                source_box.bmin[i] = source_t->bmin[i];
                source_box.bmax[i] = source_t->bmax[i];
            }
            detail::ChebyshevRn<D> cheb_rn(m_order,target_box);
            for (int j=0; j<m_ncheb; ++j) {
                const double_d& pj_unit_box = m_cheb_points[j];
                const double_d pj = 0.5*(pj_unit_box+1)*(source_box.bmax-source_box.bmin) 
                    + source_box.bmin;
                cheb_rn.set_position(pj);

                lattice_iterator<D> mi(int_d::Constant(0),int_d::Constant(m_order));
                for (int i=0; i<m_ncheb; ++i,++mi) {
                    const double tmp = cheb_rn(*mi);
                    for (int ii = 0; ii < BlockCols; ++ii) {
                        setentry_amatrix(matrix,i*BlockCols+ii,j*BlockCols+ii,tmp);
                    }
                }
            }
        }

        void M2L_amatrix(pamatrix matrix, 
                 pccluster target_t, 
                 pccluster source_t) const {
            // don't resize, already done in new_uniform
            //resize_amatrix(matrix,m_ncheb,m_ncheb);
            ASSERT_CUDA(matrix->rows == m_ncheb);
            ASSERT_CUDA(matrix->cols == m_ncheb);
            box_type target_box,source_box;
            for (int i = 0; i < D; ++i) {
                target_box.bmin[i] = target_t->bmin[i];
                target_box.bmax[i] = target_t->bmax[i];
                source_box.bmin[i] = source_t->bmin[i];
                source_box.bmax[i] = source_t->bmax[i];
            }
            for (int i=0; i<m_ncheb; ++i) {
                const double_d& pi_unit_box = m_cheb_points[i];
                const double_d pi = 0.5*(pi_unit_box+1)*(target_box.bmax-target_box.bmin) 
                                                                    + target_box.bmin;
                for (int j=0; j<m_ncheb; ++j) {
                    const double_d& pj_unit_box = m_cheb_points[j];
                    const double_d pj = 0.5*(pj_unit_box+1)*(source_box.bmax-source_box.bmin) 
                                                                    + source_box.bmin;
#ifdef HAVE_EIGEN
                    const Eigen::Matrix<double,BlockRows,BlockCols> tmp(m_K(pi,pj));
                    for (int ii = 0; ii < BlockRows; ++ii) {
                        for (int jj = 0; jj < BlockCols; ++jj) {
                            setentry_amatrix(matrix,i*BlockRows+ii,j*BlockCols+jj,tmp(ii,jj));
                        }
                    }
#else
                    setentry_amatrix(matrix,i,j,m_K(pi,pj));
#endif
                }
            }
        }

        void L2L_amatrix(pamatrix matrix, 
                 pccluster target_t, 
                 pccluster source_t) const {
            //resize_amatrix(matrix,m_ncheb,m_ncheb);
            ASSERT_CUDA(matrix->rows == m_ncheb);
            ASSERT_CUDA(matrix->cols == m_ncheb);
            clear_amatrix(matrix);
            box_type target_box,source_box;
            for (int i = 0; i < D; ++i) {
                target_box.bmin[i] = target_t->bmin[i];
                target_box.bmax[i] = target_t->bmax[i];
                source_box.bmin[i] = source_t->bmin[i];
                source_box.bmax[i] = source_t->bmax[i];
            }
            detail::ChebyshevRn<D> cheb_rn(m_order,source_box);
            for (int i=0; i<m_ncheb; ++i) {
                const double_d& pi_unit_box = m_cheb_points[i];
                const double_d pi = 0.5*(pi_unit_box+1)*(target_box.bmax-target_box.bmin) 
                    + target_box.bmin;
                cheb_rn.set_position(pi);
                lattice_iterator<D> mj(int_d::Constant(0),int_d::Constant(m_order));
                for (int j=0; j<m_ncheb; ++j,++mj) {
                    const double tmp = cheb_rn(*mj);
                    for (int ii = 0; ii < BlockRows; ++ii) {
                        setentry_amatrix(matrix,i*BlockRows+ii,j*BlockRows+ii,tmp);
                    }
                }
            }
        }

        template <typename ParticlesType>
        void L2P_amatrix(pamatrix matrix, 
                    const pccluster t,
                    const uint* indicies,
                    const uint indicies_size,
                    const ParticlesType& particles) const {
            typedef typename ParticlesType::position position;
            //resize_amatrix(matrix,indicies_size,m_ncheb);
            ASSERT_CUDA(matrix->rows == indicies_size);
            ASSERT_CUDA(matrix->cols == m_ncheb);
            clear_amatrix(matrix);
            box_type box;
            for (int i = 0; i < D; ++i) {
                box.bmin[i] = t->bmin[i];
                box.bmax[i] = t->bmax[i];
            }
            detail::ChebyshevRn<D> cheb_rn(m_order,box);
            for (int i = 0; i < indicies_size; ++i) {
                const double_d& p = get<position>(particles)[indicies[i]];
                cheb_rn.set_position(p);
                lattice_iterator<dimension> mj(int_d::Constant(0),int_d::Constant(m_order));
                for (int j=0; j<m_ncheb; ++j,++mj) {
                    const double tmp = cheb_rn(*mj);
                    for (int ii = 0; ii < BlockRows; ++ii) {
                        setentry_amatrix(matrix,i*BlockRows+ii,j*BlockRows+ii,tmp);
                    }
                }
                
            }
        }

        template <typename ParticlesType>
        void P2M_amatrix(pamatrix matrix, 
                    const pccluster t,
                    const uint* indicies,
                    const uint indicies_size,
                    const ParticlesType& particles) const {
            typedef typename ParticlesType::position position;
            ASSERT_CUDA(matrix->rows == m_ncheb);
            ASSERT_CUDA(matrix->cols == indicies_size);
            clear_amatrix(matrix);
            box_type box;
            for (int i = 0; i < D; ++i) {
                box.bmin[i] = t->bmin[i];
                box.bmax[i] = t->bmax[i];
            }
            detail::ChebyshevRn<D> cheb_rn(m_order,box);
            for (int i = 0; i < indicies_size; ++i) {
                const double_d& p = get<position>(particles)[indicies[i]];
                cheb_rn.set_position(p);
                lattice_iterator<dimension> mj(int_d::Constant(0),int_d::Constant(m_order));
                for (int j=0; j<m_ncheb; ++j,++mj) {
                    const double tmp = cheb_rn(*mj);
                    for (int ii = 0; ii < BlockCols; ++ii) {
                        setentry_amatrix(matrix,j*BlockCols+ii,i*BlockCols+ii,tmp);
                    }
                }
                
            }
        }

    };
#endif

    template <typename Expansions,
              typename Traits, 
              typename T,
                    typename VectorType=typename Expansions::expansion_type,
                    typename SourceParticleIterator=typename Traits::raw_pointer, 
                    unsigned int D=Traits::dimension>
    void calculate_P2M(VectorType& sum, 
                        const detail::bbox<D>& box, 
                        const iterator_range<ranges_iterator<Traits>>& range, 
                        const std::vector<T>& source_vector,
                        const SourceParticleIterator& source_particles_begin,
                        const Expansions& expansions) {
        typedef typename Traits::position position;
        const size_t N = std::distance(range.begin(),range.end());
        const Vector<double,D>* pbegin = &get<position>(*range.begin());
        const size_t index = pbegin - &get<position>(source_particles_begin)[0];
        for (int i = 0; i < N; ++i) {
            const Vector<double,D>& pi = pbegin[i]; 
            expansions.P2M(sum,box,pi,source_vector[index+i]);  
        }
    }

    // assume serial processing of particles, this could be more efficient for ranges iterators
    template <typename Expansions,
                typename Iterator, 
                typename T,
                 typename VectorType=typename Expansions::expansion_type, 
                 typename Traits=typename Iterator::traits_type,
                 typename SourceParticleIterator=typename Traits::raw_pointer, 
                 unsigned int D=Traits::dimension,
                 typename = typename
        std::enable_if<!std::is_same<Iterator,ranges_iterator<Traits>>::value>>
    void calculate_P2M(VectorType& sum, 
                        const detail::bbox<D>& box, 
                        const iterator_range<Iterator>& range, 
                        const std::vector<T>& source_vector,
                        const SourceParticleIterator& source_particles_begin,
                        const Expansions &expansions) {

        typedef typename Traits::position position;
        typedef typename Iterator::reference reference;
        for (reference i: range) {
            const Vector<double,D>& pi = get<position>(i); 
            const size_t index = &pi- &get<position>(source_particles_begin)[0];
            expansions.P2M(sum,box,pi,source_vector[index]);  
        }

    }

#ifdef HAVE_EIGEN

    //TODO: move this somewhere sensible
    template <size_t Size>
    struct ConvertToDoubleOrVector {
        typedef typename std::conditional<Size==1,
                         double,
                         Vector<double,Size>>::type return_type;

        template <typename Derived>
        static return_type convert(const Eigen::DenseBase<Derived>& vector) {
            return vector;
        }
    };

    template <>
    struct ConvertToDoubleOrVector<1> {
        typedef double return_type;

        template <typename Derived>
        static return_type convert(const Eigen::DenseBase<Derived>& vector) {
            return vector[0];
        }
    };

    template <typename Expansions,
              typename Traits, 
              typename Derived,
                    typename VectorType=typename Expansions::expansion_type,
                    typename SourceParticleIterator=typename Traits::raw_pointer, 
                    unsigned int D=Traits::dimension>
    void calculate_P2M(VectorType& sum, 
                        const detail::bbox<D>& box, 
                        const iterator_range<ranges_iterator<Traits>>& range, 
                        const Eigen::DenseBase<Derived>& source_vector,
                        const SourceParticleIterator& source_particles_begin,
                        const Expansions& expansions) {
        typedef typename Traits::position position;
        const size_t N = std::distance(range.begin(),range.end());
        const Vector<double,D>* pbegin = &get<position>(*range.begin());
        const size_t index = pbegin - &get<position>(source_particles_begin)[0];
        const size_t block_size = expansions.block_cols;
        for (int i = 0; i < N; ++i) {
            const Vector<double,D>& pi = pbegin[i]; 
            expansions.P2M(sum,box,pi,
                    ConvertToDoubleOrVector<block_size>::convert(
                        source_vector.template segment<block_size>((index+i)*block_size)));  
        }
    }

    

    // assume serial processing of particles, this could be more efficient for ranges iterators
    template <typename Expansions,
                typename Iterator, 
                typename Derived,
                 typename VectorType=typename Expansions::expansion_type, 
                 typename Traits=typename Iterator::traits_type,
                 typename SourceParticleIterator=typename Traits::raw_pointer, 
                 unsigned int D=Traits::dimension,
                 typename = typename
        std::enable_if<!std::is_same<Iterator,ranges_iterator<Traits>>::value>>
    void calculate_P2M(VectorType& sum, 
                        const detail::bbox<D>& box, 
                        const iterator_range<Iterator>& range, 
                        const Eigen::DenseBase<Derived>& source_vector,
                        const SourceParticleIterator& source_particles_begin,
                        const Expansions &expansions) {

        typedef typename Traits::position position;
        typedef typename Iterator::reference reference;
        const size_t block_size = expansions.block_cols;
        for (reference i: range) {
            const Vector<double,D>& pi = get<position>(i); 
            const size_t index = &pi- &get<position>(source_particles_begin)[0];
            expansions.P2M(sum,box,pi,
                    ConvertToDoubleOrVector<block_size>::convert(
                    source_vector.template segment<block_size>(index*block_size)));  
        }

    }
#endif

    template <typename Expansions,
              typename Traits, 
              typename T, 
                    typename VectorType=typename Expansions::expansion_type,
                    typename ParticleIterator=typename Traits::raw_pointer, 
                    unsigned int D=Traits::dimension>
    void calculate_L2P(
                        std::vector<T>& target_vector,
                        const VectorType& source, 
                        const detail::bbox<D>& box, 
                        const iterator_range<ranges_iterator<Traits>>& range, 
                        const ParticleIterator& target_particles_begin,
                        const Expansions& expansions) {
        typedef typename Traits::position position;
        LOG(3,"calculate_L2P (range): box = "<<box);
        const size_t N = std::distance(range.begin(),range.end());
        const Vector<double,D>* pbegin_range = &get<position>(*range.begin());
        const Vector<double,D>* pbegin = &get<position>(target_particles_begin)[0];
        const size_t index = pbegin_range - pbegin;
        for (int i = index; i < index+N; ++i) {
            const Vector<double,D>& pi = pbegin[i]; 
            target_vector[i] += expansions.L2P(pi,box,source);  
        }
    }

    // assume serial processing of particles, this could be more efficient for ranges iterators
    template <typename Expansions,
                typename Iterator, 
              typename T, 
                 typename VectorType=typename Expansions::expansion_type, 
                 typename Traits=typename Iterator::traits_type,
                 typename ParticleIterator=typename Traits::raw_pointer, 
                 unsigned int D=Traits::dimension,
                 typename = typename
        std::enable_if<!std::is_same<Iterator,ranges_iterator<Traits>>::value>>
    void calculate_L2P(
                        std::vector<T>& target_vector,
                        const VectorType& source, 
                        const detail::bbox<D>& box, 
                        const iterator_range<Iterator>& range, 
                        const ParticleIterator& target_particles_begin,
                        const Expansions &expansions) {

        LOG(3,"calculate_L2P: box = "<<box);
        typedef typename Traits::position position;
        typedef typename Iterator::reference reference;
        for (reference i: range) {
            const Vector<double,D>& pi = get<position>(i); 
            const size_t index = &pi- &get<position>(target_particles_begin)[0];
            target_vector[index] += expansions.L2P(pi,box,source);  
        }

    }

#ifdef HAVE_EIGEN
    template <typename Expansions,
              typename Traits, 
                typename Derived, 
                    typename VectorType=typename Expansions::expansion_type,
                    typename ParticleIterator=typename Traits::raw_pointer, 
                    unsigned int D=Traits::dimension>
    void calculate_L2P(
                        const Eigen::DenseBase<Derived>& target_vector,
                        const VectorType& source, 
                        const detail::bbox<D>& box, 
                        const iterator_range<ranges_iterator<Traits>>& range, 
                        const ParticleIterator& target_particles_begin,
                        const Expansions& expansions) {
        typedef typename Traits::position position;
        LOG(3,"calculate_L2P (range): box = "<<box);
        const size_t N = std::distance(range.begin(),range.end());
        const Vector<double,D>* pbegin_range = &get<position>(*range.begin());
        const Vector<double,D>* pbegin = &get<position>(target_particles_begin)[0];
        const size_t index = pbegin_range - pbegin;
        const size_t block_size = expansions.block_rows;
        for (int i = index; i < index+N; ++i) {
            const Vector<double,D>& pi = pbegin[i]; 
            const auto l2p = expansions.L2P(pi,box,source);
            for (int ii = 0; ii < block_size; ++ii) {
                const_cast<Eigen::DenseBase<Derived>&>(target_vector)[i*block_size+ii] += 
                    detail::VectorTraits<decltype(l2p)>::Index(l2p,ii);
            }
        }
    }

    // assume serial processing of particles, this could be more efficient for ranges iterators
    template <typename Expansions,
                typename Iterator, 
                typename Derived, 
                 typename VectorType=typename Expansions::expansion_type, 
                 typename Traits=typename Iterator::traits_type,
                 typename ParticleIterator=typename Traits::raw_pointer, 
                 unsigned int D=Traits::dimension,
                 typename = typename
        std::enable_if<!std::is_same<Iterator,ranges_iterator<Traits>>::value>>
    void calculate_L2P(
                        const Eigen::DenseBase<Derived>& target_vector,
                        const VectorType& source, 
                        const detail::bbox<D>& box, 
                        const iterator_range<Iterator>& range, 
                        const ParticleIterator& target_particles_begin,
                        const Expansions &expansions) {

        LOG(3,"calculate_L2P: box = "<<box);
        typedef typename Traits::position position;
        typedef typename Iterator::reference reference;
        const size_t block_size = expansions.block_rows;
        for (reference i: range) {
            const Vector<double,D>& pi = get<position>(i); 
            const size_t index = &pi- &get<position>(target_particles_begin)[0];
            const auto& l2p = expansions.L2P(pi,box,source);
            for (int ii = 0; ii < block_size; ++ii) {
                const_cast<Eigen::DenseBase<Derived>&>(target_vector)[index*block_size+ii] += 
                    detail::VectorTraits<decltype(l2p)>::Index(l2p,ii);
            }
        }
    }
#endif

    template <typename Kernel,
              typename Traits, 
              typename TargetType, 
              typename SourceType, 
               typename ParticleIterator=typename Traits::raw_pointer, 
               unsigned int D=Traits::dimension>
    void calculate_P2P(
                        std::vector<TargetType>& target_vector,
                        const std::vector<SourceType>& source_vector,
                        const iterator_range<ranges_iterator<Traits>>& target_range, 
                        const iterator_range<ranges_iterator<Traits>>& source_range, 
                        const ParticleIterator& target_particles_begin,
                        const ParticleIterator& source_particles_begin,
                        const Kernel& kernel) {
        typedef typename Traits::position position;

        const size_t n_target = std::distance(target_range.begin(),target_range.end());
        const size_t n_source = std::distance(source_range.begin(),source_range.end());

        const Vector<double,D>* pbegin_target_range = &get<position>(*target_range.begin());
        const Vector<double,D>* pbegin_target = &get<position>(target_particles_begin)[0];
        const size_t index_target = pbegin_target_range - pbegin_target;

        const Vector<double,D>* pbegin_source_range = &get<position>(*source_range.begin());
        const Vector<double,D>* pbegin_source = &get<position>(source_particles_begin)[0];

        const size_t index_source = pbegin_source_range - pbegin_source;

        auto pi = target_range.begin();
        for (int i = index_target; i < index_target+n_target; ++i,++pi) {
            auto pj = source_range.begin();
            for (int j = index_source; j < index_source+n_source; ++j,++pj) {
                target_vector[i] += kernel(*pi,*pj)*source_vector[j];
            }
        }
    }

    // assume serial processing of particles, this could be more efficient for ranges iterators
    template <typename Kernel,
                typename TargetIterator, 
                typename SourceIterator, 
              typename TargetType, 
              typename SourceType, 
                 typename Traits=typename TargetIterator::traits_type,
                 typename ParticleIterator=typename Traits::raw_pointer, 
                 unsigned int D=Traits::dimension,
                 typename = typename
        std::enable_if<!(std::is_same<TargetIterator,ranges_iterator<Traits>>::value
                            && std::is_same<SourceIterator,ranges_iterator<Traits>>::value)>>
    void calculate_P2P(
                        std::vector<TargetType>& target_vector,
                        const std::vector<SourceType>& source_vector,
                        const iterator_range<TargetIterator>& target_range, 
                        const iterator_range<SourceIterator>& source_range, 
                        const ParticleIterator& target_particles_begin,
                        const ParticleIterator& source_particles_begin,
                        const Kernel &kernel) {

        typedef typename Traits::position position;
        for (auto& i: target_range) {
            const size_t target_index = &get<position>(i)
                                      - &get<position>(target_particles_begin)[0];
            for (auto& j: source_range) {
                const size_t source_index = &get<position>(j) 
                                          - &get<position>(source_particles_begin)[0];
                LOG(4,"calculate_P2P: i = "<<target_index<<" j = "<<source_index);
                target_vector[target_index] += kernel(i,j)*source_vector[source_index];
            }
        }

    }

#ifdef HAVE_EIGEN
    template <typename Kernel,
              typename Traits, 
              typename TargetDerived, 
              typename SourceDerived, 
               typename ParticleIterator=typename Traits::raw_pointer, 
               unsigned int D=Traits::dimension>
    void calculate_P2P(
                        const Eigen::DenseBase<TargetDerived>& target_vector,
                        const Eigen::DenseBase<SourceDerived>& source_vector,
                        const iterator_range<ranges_iterator<Traits>>& target_range, 
                        const iterator_range<ranges_iterator<Traits>>& source_range, 
                        const ParticleIterator& target_particles_begin,
                        const ParticleIterator& source_particles_begin,
                        const Kernel& kernel) {
        typedef typename Traits::position position;
        typedef typename Traits::raw_const_reference const_row_reference;
        typedef typename Traits::raw_const_reference const_col_reference;

        const size_t n_target = std::distance(target_range.begin(),target_range.end());
        const size_t n_source = std::distance(source_range.begin(),source_range.end());

        const Vector<double,D>* pbegin_target_range = &get<position>(*target_range.begin());
        const Vector<double,D>* pbegin_target = &get<position>(target_particles_begin)[0];
        const size_t index_target = pbegin_target_range - pbegin_target;

        const Vector<double,D>* pbegin_source_range = &get<position>(*source_range.begin());
        const Vector<double,D>* pbegin_source = &get<position>(source_particles_begin)[0];

        const size_t index_source = pbegin_source_range - pbegin_source;

        typedef detail::kernel_helper_ref<const_row_reference,
                                           const_col_reference,Kernel> helper;

        auto pi = target_range.begin();
        for (int i = index_target; i < index_target+n_target; ++i,++pi) {
            auto pj = source_range.begin();
            for (int j = index_source; j < index_source+n_source; ++j,++pj) {
                const_cast<Eigen::DenseBase<TargetDerived>&>(target_vector)
                    .template segment<helper::block_rows>(i*helper::block_rows) += 
                    kernel(*pi,*pj)*source_vector.template segment<helper::block_cols>(j*helper::block_cols);
            }
        }
    }

    // assume serial processing of particles, this could be more efficient for ranges iterators
    template <typename Kernel,
                typename TargetIterator, 
                typename SourceIterator, 
              typename TargetDerived, 
              typename SourceDerived, 
                 typename Traits=typename TargetIterator::traits_type,
                 typename ParticleIterator=typename Traits::raw_pointer, 
                 unsigned int D=Traits::dimension,
                 typename = typename
        std::enable_if<!(std::is_same<TargetIterator,ranges_iterator<Traits>>::value
                            && std::is_same<SourceIterator,ranges_iterator<Traits>>::value)>>
    void calculate_P2P(
                        const Eigen::DenseBase<TargetDerived>& target_vector,
                        const Eigen::DenseBase<SourceDerived>& source_vector,
                        const iterator_range<TargetIterator>& target_range, 
                        const iterator_range<SourceIterator>& source_range, 
                        const ParticleIterator& target_particles_begin,
                        const ParticleIterator& source_particles_begin,
                        const Kernel &kernel) {

        typedef typename Traits::position position;
        typedef typename Traits::raw_const_reference const_row_reference;
        typedef typename Traits::raw_const_reference const_col_reference;

        // TODO: lots of common code here to consolodate
        typedef typename std::result_of<Kernel(const_row_reference, 
                                               const_col_reference)>::type FunctionReturn;
        typedef typename std::conditional<
                            std::is_arithmetic<FunctionReturn>::value,
                            Eigen::Matrix<FunctionReturn,1,1>,
                            FunctionReturn>::type Block;

        const int block_rows = Block::RowsAtCompileTime; 
        const int block_cols = Block::ColsAtCompileTime; 

        static_assert(block_rows > 0,"kernel function must return fixed size matrix");
        static_assert(block_cols > 0,"kernel function must return fixed size matrix");

        for (auto& i: target_range) {
            const size_t target_index = &get<position>(i)
                                      - &get<position>(target_particles_begin)[0];
            for (auto& j: source_range) {
                const size_t source_index = &get<position>(j) 
                                          - &get<position>(source_particles_begin)[0];
                LOG(4,"calculate_P2P: i = "<<target_index<<" j = "<<source_index);
                const_cast<Eigen::DenseBase<TargetDerived>&>(target_vector)
                    .template segment<block_rows>(target_index*block_rows) += 
                    kernel(i,j)*source_vector.template segment<block_cols>(source_index*block_cols);
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
        for (int i = 0; i < N; ++i) {
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
    void P2P_matrix(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& matrix, 
                const std::vector<size_t>& row_indicies,
                const std::vector<size_t>& col_indicies,
                const RowParticlesType& row_particles,
                const ColParticlesType& col_particles,
                const Kernel& kernel) {
       
        typedef detail::kernel_helper_ref<typename RowParticlesType::const_reference,
                                          typename ColParticlesType::const_reference,
                                          Kernel> helper;

        typedef typename ColParticlesType::position position;
        matrix.resize(row_indicies.size(),col_indicies.size());
        for (int i = 0; i < row_indicies.size(); ++i) {
            for (int j = 0; j < col_indicies.size(); ++j) {
                matrix.block<helper::block_rows,helper::block_cols>(
                        i*helper::block_rows,j*helper::block_cols) 
                                = kernel(row_particles[i],col_particles[j]);
            }
        }
    }
#endif

#ifdef HAVE_H2LIB
    template <typename RowParticlesType, typename ColParticlesType, typename Kernel>
    void P2P_amatrix(pamatrix matrix, 
                const uint* row_indicies,
                const uint row_indicies_size,
                const uint* col_indicies,
                const uint col_indicies_size,
                const RowParticlesType& row_particles,
                const ColParticlesType& col_particles,
                const Kernel& kernel) {
        typedef typename ColParticlesType::position position;
        ASSERT_CUDA(matrix->rows == row_indicies_size);
        ASSERT_CUDA(matrix->cols == col_indicies_size);

        typedef detail::kernel_helper_ref<typename RowParticlesType::const_reference,
                                          typename ColParticlesType::const_reference,
                                          Kernel> helper;

        //resize_amatrix(matrix,row_indicies_size,col_indicies_size);
        for (int i = 0; i < row_indicies_size; ++i) {
            const auto& pi = row_particles[row_indicies[i]];
            for (int j = 0; j < col_indicies_size; ++j) {
                const auto& pj = col_particles[col_indicies[j]];
                const typename helper::Block tmp(kernel(pi,pj));
                for (int ii = 0; ii < helper::block_rows; ++ii) {
                    for (int jj = 0; jj < helper::block_cols; ++jj) {
                        setentry_amatrix(matrix,i,j,tmp(ii,jj));
                    }
                    
                }
            }
            
        }
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

    template <unsigned int D>
    struct theta_condition {
        typedef Vector<double,D> double_d;
        const double_d& m_low;
        const double_d& m_high;
        double m_max_diam; 
        static constexpr double m_theta = 0.5;
        static constexpr double m_theta2 = m_theta*m_theta;
        theta_condition(const double_d& low, const double_d& high):
            m_low(low),m_high(high),m_max_diam(std::numeric_limits<double>::min())
        {
            for (int i = 0; i < D; ++i) {
                if (high[i]-low[i] > m_max_diam) m_max_diam = high[i]-low[i];
            }
        }

        bool check(const double_d& low, const double_d& high) const {
            double_d dist = double_d::Zero();
            double max_diam = m_max_diam;
            for (int i = 0; i < D; ++i) {
                if (m_high[i] < low[i]) {
                    dist[i] = low[i]-m_high[i];
                } else if (high[i] < m_low[i]) {
                    dist[i] = m_low[i] - high[i];
                }
                if (high[i]-low[i] > max_diam) max_diam = high[i]-low[i];
            }
            return std::pow(max_diam,2) > 1.0*dist.squaredNorm();
        }
    };

}
}

#endif
