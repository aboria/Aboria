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


#ifndef KERNELS_H_
#define KERNELS_H_

#include <type_traits>

#include "detail/Chebyshev.h"
#include <Eigen/Core>

#include "FastMultipoleMethod.h"
#include "detail/Kernels.h"

#ifdef HAVE_H2LIB
#include "H2Lib.h"
#endif



namespace Aboria {


    template<typename RowElements, typename ColElements, typename F>
    class KernelBase {
    protected:
        typedef typename detail::kernel_helper<RowElements,ColElements,F> kernel_helper;
        static const unsigned int dimension = kernel_helper::dimension;
        typedef Vector<double,dimension> double_d;
        typedef Vector<int,dimension> int_d;
        typedef typename kernel_helper::const_position_reference const_position_reference;
        typedef typename kernel_helper::const_row_reference const_row_reference;
        typedef typename kernel_helper::const_col_reference const_col_reference;

    public:

        typedef typename kernel_helper::Block Block;
        typedef typename Block::Scalar Scalar;
        static const size_t BlockRows = Block::RowsAtCompileTime;
        static const size_t BlockCols = Block::ColsAtCompileTime;
        typedef typename Eigen::Matrix<Scalar,BlockCols,1> BlockRHSVector;
        typedef typename Eigen::Matrix<Scalar,BlockRows,1> BlockLHSVector;

        typedef RowElements row_particles_type;
        typedef ColElements col_particles_type;
        typedef F function_type;
        typedef size_t Index;
        enum {
            ColsAtCompileTime = -1,
            RowsAtCompileTime = -1
        };

        KernelBase(const RowElements& row_particles,
                   const ColElements& col_particles,
                   const F& function): 
            m_function(function),
            m_row_particles(row_particles), 
            m_col_particles(col_particles)
        {};

        RowElements& get_row_particles() {
            return m_row_particles;
        }

        const RowElements& get_row_particles() const {
            return m_row_particles;
        }

        ColElements& get_col_particles() {
            return m_col_particles;
        }

        const function_type& get_kernel_function() const {
            return m_function;
        }

        const ColElements& get_col_particles() const {
            return m_col_particles;
        }

        size_t rows() const {
            return m_row_particles.size()*BlockRows;
        }

        size_t cols() const {
            return m_col_particles.size()*BlockCols;
        }


        Scalar coeff(const size_t i, const size_t j) const {
            ASSERT(i < rows(),"i greater than rows()");
            ASSERT(j < cols(),"j greater than cols()");
            const int pi = std::floor(static_cast<float>(i)/BlockRows);
            const int ioffset = i - pi*BlockRows;
            const int pj = std::floor(static_cast<float>(j)/BlockCols);
            const int joffset = j - pj*BlockCols;
            const_row_reference ai = m_row_particles[pi];
            const_col_reference bj = m_col_particles[pj];
            return m_function(ai,bj);
        }

        template<typename MatrixType>
        void assemble(const MatrixType &matrix) const {
                    MatrixType::MATRIX_ASSEMBLE_NOT_IMPLEMENTED; 
        }

        template<typename Triplet>
        void assemble(std::vector<Triplet>& triplets,
                      const size_t startI=0, const size_t startJ=0
                      ) const {
            Triplet::TRIPLET_ASSEMBLE_NOT_IMPLEMENTED; 
        }

        /// Evaluates a matrix-free linear operator given by \p expr \p if_expr,
        /// and particle sets \p a and \p b on a vector rhs and
        /// accumulates the result in vector lhs
        template<typename VectorLHS,typename VectorRHS>
        void evaluate(VectorLHS &lhs, const VectorRHS &rhs) const {
            VectorLHS::EVALUATE_NOT_IMPLEMENTED;
       }

    protected:
        const RowElements &m_row_particles;
        const ColElements &m_col_particles;
        const F& m_function;
    };

    template<typename RowElements, typename ColElements, typename F>
    class KernelDense: public KernelBase<RowElements,ColElements,F> {
    protected:
        typedef KernelBase<RowElements,ColElements,F> base_type;
        static const unsigned int dimension = base_type::dimension;
        typedef typename base_type::double_d double_d;
        typedef typename base_type::int_d int_d;
        typedef typename base_type::const_row_reference const_row_reference;
        typedef typename base_type::const_col_reference const_col_reference;
    public:
        typedef typename base_type::Block Block;
        typedef typename base_type::BlockLHSVector BlockLHSVector;
        typedef typename base_type::Scalar Scalar;
        static const size_t BlockRows = base_type::BlockRows;
        static const size_t BlockCols = base_type::BlockCols;

        KernelDense(const RowElements& row_particles,
                    const ColElements& col_particles,
                    const F& function): base_type(row_particles,
                                                  col_particles,
                                                  function) 
        {};

        template<typename MatrixType>
        void assemble(const MatrixType &matrix) const {

            const RowElements& a = this->m_row_particles;
            const ColElements& b = this->m_col_particles;
            const size_t na = a.size();
            const size_t nb = b.size();

            ASSERT(matrix.rows() == this->rows(),"matrix has incompatible row size");
            ASSERT(matrix.cols() == this->cols(),"matrix has incompatible col size");

            const bool is_periodic = !a.get_periodic().any();

            for (size_t i=0; i<na; ++i) {
                const_row_reference ai = a[i];
                for (size_t j=0; j<nb; ++j) {
                    const_col_reference  bj = b[j];
                    const_cast< MatrixType& >(matrix).template block<BlockRows,BlockCols>(i*BlockRows,j*BlockCols) = Block(this->m_function(ai,bj));
                }
            }
        }

        template<typename Triplet>
        void assemble(std::vector<Triplet>& triplets,
                      const size_t startI=0, const size_t startJ=0
                      ) const {

            const RowElements& a = this->m_row_particles;
            const ColElements& b = this->m_col_particles;

            const size_t na = a.size();
            const size_t nb = b.size();

            const bool is_periodic = !a.get_periodic().any();

            for (size_t i=0; i<na; ++i) {
                const_row_reference ai = a[i];
                for (size_t j=0; j<nb; ++j) {
                    const_col_reference bj = b[j];
                    const Block element = this->m_function(ai,bj);
                    for (int ii = 0; ii < BlockRows; ++ii) {
                        for (int jj = 0; jj < BlockCols; ++jj) {
                            triplets.push_back(Triplet(i*BlockRows+ii,
                                                     j*BlockCols+jj,
                                                     element(ii,jj)));
                        }
                    }
                }
            }
        }

        /// Evaluates a matrix-free linear operator given by \p expr \p if_expr,
        /// and particle sets \p a and \p b on a vector rhs and
        /// accumulates the result in vector lhs
        template <typename DerivedLHS, typename DerivedRHS>
        void evaluate(Eigen::DenseBase<DerivedLHS> &lhs, 
                const Eigen::DenseBase<DerivedRHS> &rhs) const {
            typedef Eigen::Matrix<Scalar,BlockRows,1> row_vector;

            const RowElements& a = this->m_row_particles;
            const ColElements& b = this->m_col_particles;

            const size_t na = a.size();
            const size_t nb = b.size();

            CHECK(lhs.size() == this->rows(),"lhs size is inconsistent");
            CHECK(rhs.size() == this->cols(),"rhs size is inconsistent");

            const bool is_periodic = !a.get_periodic().any();

            #pragma omp parallel for
            for (size_t i=0; i<na; ++i) {
                const_row_reference ai = a[i];
                for (size_t j=0; j<nb; ++j) {
                    const_col_reference bj = b[j];
                    lhs.template segment<BlockRows>(i*BlockRows) += 
                                this->m_function(ai,bj)
                                    *rhs.template segment<BlockCols>(j*BlockCols);
                }
            }
       }

        /// Evaluates a matrix-free linear operator given by \p expr \p if_expr,
        /// and particle sets \p a and \p b on a vector rhs and
        /// accumulates the result in vector lhs
        template <typename LHSType, typename RHSType>
        void evaluate(std::vector<LHSType> &lhs, 
                const std::vector<RHSType> &rhs) const {

            const RowElements& a = this->m_row_particles;
            const ColElements& b = this->m_col_particles;

            const size_t na = a.size();
            const size_t nb = b.size();

            CHECK(lhs.size() == na,"lhs size is inconsistent");
            CHECK(rhs.size() == nb,"rhs size is inconsistent");

            const bool is_periodic = !a.get_periodic().any();

            #pragma omp parallel for
            for (size_t i=0; i<na; ++i) {
                const_row_reference ai = a[i];
                for (size_t j=0; j<nb; ++j) {
                    const_col_reference bj = b[j];
                    lhs[i] += this->m_function(ai,bj)*rhs[j];
                }
            }
       }

    };

    template<typename RowElements, typename ColElements, typename F>
    class KernelMatrix: public KernelBase<RowElements,ColElements,F> {
    protected:
        typedef KernelBase<RowElements,ColElements,F> base_type;
        static const unsigned int dimension = base_type::dimension;
        typedef typename base_type::double_d double_d;
        typedef typename base_type::int_d int_d;
        typedef typename base_type::const_row_reference const_row_reference;
        typedef typename base_type::const_col_reference const_col_reference;

        typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_type;

        matrix_type m_matrix;
    public:
        typedef typename base_type::Scalar Scalar;
        typedef typename base_type::Block Block;
        static const size_t BlockRows = base_type::BlockRows;
        static const size_t BlockCols = base_type::BlockCols;

        KernelMatrix(const RowElements& row_particles,
                    const ColElements& col_particles,
                    const F& function): base_type(row_particles,
                                                  col_particles,
                                                  function) 
        {
            assemble_matrix(); 
        };

        void assemble_matrix() {
            const RowElements& a = this->m_row_particles;
            const ColElements& b = this->m_col_particles;

            const bool is_periodic = !a.get_periodic().any();

            m_matrix.resize(this->rows(),this->cols());
            for (size_t i=0; i<a.size(); ++i) {
                const_row_reference ai = a[i];
                for (size_t j=0; j<b.size(); ++j) {
                    const_col_reference bj = b[j];
                    m_matrix.template block<BlockRows,BlockCols>(
                                           i*BlockRows,j*BlockCols) = 
                                                Block(this->m_function(ai,bj));
                }
            }
        }

        Scalar coeff(const size_t i, const size_t j) const {
            return m_matrix(i,j);
        }

        template<typename MatrixType>
        void assemble(const MatrixType &matrix) const {
            const_cast< MatrixType& >(matrix) = m_matrix;
        }

        /// Evaluates a matrix-free linear operator given by \p expr \p if_expr,
        /// and particle sets \p a and \p b on a vector rhs and
        /// accumulates the result in vector lhs
        template <typename LHSType, typename RHSType>
        void evaluate(std::vector<LHSType> &lhs, 
                const std::vector<RHSType> &rhs) const {


            const RowElements& a = this->m_row_particles;
            const ColElements& b = this->m_col_particles;

            const size_t na = a.size();
            const size_t nb = b.size();

            ASSERT(lhs.size() == na,"lhs size is inconsistent");
            ASSERT(rhs.size() == nb,"rhs size is inconsistent");

            const bool is_periodic = !a.get_periodic().any();

            #pragma omp parallel for
            for (size_t i=0; i<na; ++i) {
                for (size_t j=0; j<nb; ++j) {
                    lhs[i] += m_matrix.template block<BlockRows,BlockCols>(
                                            i*BlockRows,j*BlockCols)
                                         *rhs[j];
                }
            }
        }

        /// Evaluates a matrix-free linear operator given by \p expr \p if_expr,
        /// and particle sets \p a and \p b on a vector rhs and
        /// accumulates the result in vector lhs
        template <typename DerivedLHS, typename DerivedRHS>
        void evaluate(Eigen::DenseBase<DerivedLHS> &lhs, 
                const Eigen::DenseBase<DerivedRHS> &rhs) const {
            ASSERT(lhs.size() == this->rows(),"lhs size not consistent")
            ASSERT(rhs.size() == this->cols(),"lhs size not consistent")
            lhs += m_matrix*rhs;
        }
 
    };


    template<typename RowElements, typename ColElements, typename PositionF,
        size_t QuadratureOrder = 8,
        typename ElementF=detail::position_kernel<RowElements,ColElements,PositionF>,
        typename F=detail::integrate_over_element<RowElements,ColElements,ElementF,QuadratureOrder>>
    class KernelChebyshev: public KernelDense<RowElements,ColElements,F> {
        typedef KernelDense<RowElements,ColElements,F> base_type;
        typedef typename base_type::double_d double_d;
        typedef typename base_type::int_d int_d;
        typedef typename base_type::const_position_reference const_position_reference;
        typedef typename base_type::const_row_reference const_row_reference;
        typedef typename base_type::const_col_reference const_col_reference;

        typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_type;
        typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type; 
        typedef Eigen::Map<vector_type> map_type;

        static const unsigned int dimension = base_type::dimension;
        detail::Chebyshev_Rn<dimension> col_Rn,row_Rn;
        matrix_type m_col_Rn_matrix,m_row_Rn_matrix,m_kernel_matrix;
        unsigned int m_order;
        unsigned int m_ncheb;
        const int_d m_start;
        const int_d m_end;
        mutable vector_type m_W;
        mutable vector_type m_fcheb;
        PositionF m_position_function;

    public:
        typedef typename base_type::Scalar Scalar;
        typedef typename base_type::Block Block;
        static const size_t BlockRows = base_type::BlockRows;
        static const size_t BlockCols = base_type::BlockCols;

        KernelChebyshev(const RowElements& row_particles,
                        const ColElements& col_particles,
                        const unsigned int n,
                        const PositionF& function): m_order(n),
                                            m_ncheb(std::pow(n,dimension)),
                                            m_start(0),
                                            m_end(n),
                                            m_position_function(function),
                                            base_type(row_particles,
                                                  col_particles,
                                                  F(function)) {
            set_n(n);
        };

        void set_n(const unsigned int n) { 
            m_order = n; 
            m_ncheb = std::pow(n,dimension); 
            m_W.resize(m_ncheb*BlockCols);
            m_fcheb.resize(m_ncheb*BlockRows);

            update_row_positions();
            update_col_positions();
            update_kernel_matrix();
        }

        void update_row_positions() {
            detail::integrate_chebyshev<RowElements,BlockRows,QuadratureOrder> integrate(
                                                this->row_elements,
                                                m_order,
                                                detail::bbox<dimension>(
                                                    this->row_elements.get_min(),
                                                    this->row_elements.get_max()));
                                                       
            // fill row_Rn matrix
            m_row_Rn_matrix.resize(this->row_elements.size()*BlockRows,m_ncheb*BlockRows);
            integrate(m_row_Rn_matrix);
        }

        void update_kernel_matrix() {
            // fill kernel matrix
            m_kernel_matrix.resize(m_ncheb*BlockRows,m_ncheb*BlockCols);
            lattice_iterator<dimension> mi(m_start,m_end);
            for (int i=0; i<m_ncheb; ++i,++mi) {
                const double_d pi = col_Rn.get_position(*mi);
                lattice_iterator<dimension> mj(m_start,m_end);
                for (int j=0; j<m_ncheb; ++j,++mj) {
                    const double_d pj = row_Rn.get_position(*mj);
                    m_kernel_matrix.template block<BlockRows,BlockCols>(
                                        i*BlockRows,j*BlockCols) =
                                            Block(m_position_function(pi-pj));
                }
            }
        }

        void update_col_positions() {
            detail::integrate_chebyshev<ColElements,BlockCols,QuadratureOrder> integrate(
                                            this->col_elements,
                                            m_order,
                                            detail::bbox<dimension>(
                                                    this->col_elements.get_min(),
                                                    this->col_elements.get_max()));
             
            // fill row_Rn matrix
            m_col_Rn_matrix.resize(m_ncheb*BlockCols,this->col_elements.size()*BlockCols);
            integrate(m_col_Rn_matrix.transpose());
            
        }

        /// Evaluates a matrix-free linear operator given by \p expr \p if_expr,
        /// and particle sets \p a and \p b on a vector rhs and
        /// accumulates the result in vector lhs
        template <typename DerivedLHS, typename DerivedRHS>
        void evaluate(Eigen::DenseBase<DerivedLHS> &lhs, 
                const Eigen::DenseBase<DerivedRHS> &rhs) const {

            const RowElements& a = this->m_row_particles;
            const ColElements& b = this->m_col_particles;

            CHECK(!b.get_periodic().any(),"chebyshev operator assumes not periodic");
            ASSERT(this->rows() == lhs.rows(),"lhs vector has incompatible size");
            ASSERT(this->cols() == rhs.rows(),"rhs vector has incompatible size");

            //First compute the weights at the Chebyshev nodes ym 
            //by anterpolation 
            m_W = m_col_Rn_matrix*rhs;

            //Next compute f ðxÞ at the Chebyshev nodes xl:
            m_fcheb = m_kernel_matrix*m_W;

            //Last compute f ðxÞ at the observation points xi by interpolation:
            lhs = m_row_Rn_matrix*m_fcheb;
        }
    };

#ifdef HAVE_H2LIB

    template<typename RowElements, typename ColElements, typename PositionF,
        typename F=detail::position_lambda<RowElements,ColElements,PositionF>>
    class KernelH2: public KernelDense<RowElements,ColElements,F> {
        typedef KernelDense<RowElements,ColElements,F> base_type;
        typedef typename ColElements::query_type query_type;
        static const unsigned int dimension = base_type::dimension;
        typedef typename detail::H2LibBlackBoxExpansions<dimension,PositionF> expansions_type;
        typedef H2LibMatrix h2_matrix_type;

        h2_matrix_type m_h2_matrix;
        PositionF m_position_function;

    public:
        typedef typename base_type::Block Block;
        typedef typename base_type::Scalar Scalar;
        static const size_t BlockRows = base_type::BlockRows;
        static const size_t BlockCols = base_type::BlockCols;
        static_assert(BlockRows==1, "only implemented for scalar kernels");
        static_assert(BlockCols==1, "only implemented for scalar kernels");
        static_assert(is_particles<RowElements>::value && 
                      is_particles<ColElements>::value,
                      "only implemented for particle elements");

        typedef PositionF position_function_type;

        KernelH2(const RowElements& row_particles,
                        const ColElements& col_particles,
                        const PositionF& function,
                        const int order, 
                        const double eta = 1.0):
                          m_h2_matrix(row_particles,col_particles,
                                  expansions_type(
                                      order,
                                      function),eta),
                          m_position_function(function),
                          base_type(row_particles,
                                  col_particles,
                                  F(function)) {
                          }

        const h2_matrix_type& get_h2_matrix() const {
            return m_h2_matrix;
        }

        const PositionF& get_position_function() const {
            return m_position_function;
        }


        /// Evaluates a h2 matrix linear operator given by \p expr \p if_expr,
        /// and particle sets \p a and \p b on a vector rhs and
        /// accumulates the result in vector lhs
        template<typename VectorLHS,typename VectorRHS>
        void evaluate(VectorLHS &lhs, const VectorRHS &rhs) const {
            m_h2_matrix.matrix_vector_multiply(lhs,1.0,false,rhs);
        }
    };
#endif


    template<typename RowElements, typename ColElements, typename PositionF,
        unsigned int N,
        typename F=detail::position_kernel<RowElements,ColElements,PositionF>>
    class KernelFMM: public KernelDense<RowElements,ColElements,F> {
        typedef KernelDense<RowElements,ColElements,F> base_type;
        typedef typename base_type::position position;
        typedef typename base_type::double_d double_d;
        typedef typename base_type::int_d int_d;
        typedef typename base_type::const_position_reference const_position_reference;
        typedef typename base_type::const_row_reference const_row_reference;
        typedef typename base_type::const_col_reference const_col_reference;
        typedef typename ColElements::query_type query_type;
        static const unsigned int dimension = base_type::dimension;
        typedef typename detail::BlackBoxExpansions<dimension,N,PositionF> expansions_type;
        typedef FastMultipoleMethod<expansions_type,ColElements> fmm_type;


        expansions_type m_expansions;
        fmm_type m_fmm;

    public:
        typedef typename base_type::Block Block;
        typedef typename base_type::Scalar Scalar;
        static const size_t BlockRows = base_type::BlockRows;
        static const size_t BlockCols = base_type::BlockCols;
        static_assert(BlockRows==1, "only implemented for scalar kernels");
        static_assert(BlockCols==1, "only implemented for scalar kernels");
        static_assert(detail::is_particles<RowElements>::value && 
                      detail::is_particles<ColElements>::value,
                      "only implemented for particle elements");


        KernelFMM(const RowElements& row_particles,
                        const ColElements& col_particles,
                        const PositionF& function): 
                                            m_expansions(function),
                                            m_fmm(col_particles,
                                                  m_expansions),
                                            base_type(row_particles,
                                                  col_particles,
                                                  F(function)) {
        };

        /// Evaluates a matrix-free linear operator given by \p expr \p if_expr,
        /// and particle sets \p a and \p b on a vector rhs and
        /// accumulates the result in vector lhs
        template<typename VectorLHS,typename VectorRHS>
        void evaluate(VectorLHS &lhs, const VectorRHS &rhs) const {
            m_fmm.matrix_vector_multiply(this->m_row_particles,lhs,rhs);
        }
    };

    template<typename RowElements, typename ColElements, typename FRadius, typename F>
    class KernelSparse: public KernelBase<RowElements,ColElements,F> {
    protected:
        typedef KernelBase<RowElements,ColElements,F> base_type;
        typedef typename base_type::position position;
        typedef typename base_type::double_d double_d;
        typedef typename base_type::const_position_reference const_position_reference;
        typedef typename base_type::const_row_reference const_row_reference;
        typedef typename base_type::const_col_reference const_col_reference;
        static_assert(detail::is_particles<RowElements>::value && 
                      detail::is_particles<ColElements>::value,
                      "only implemented for particle elements");
    public:
        typedef typename base_type::Block Block;
        typedef typename base_type::Scalar Scalar;
        static const size_t BlockRows = base_type::BlockRows;
        static const size_t BlockCols = base_type::BlockCols;

        KernelSparse(const RowElements& row_particles,
                    const ColElements& col_particles,
                    const FRadius& radius_function,
                    const F& function): m_radius_function(radius_function),
                                        base_type(row_particles,
                                                  col_particles,
                                                  function) 
        {};

        Block coeff(const size_t i, const size_t j) const {
            const int pi = std::floor(static_cast<float>(i)/BlockRows);
            const int ioffset = i - pi*BlockRows;
            const int pj = std::floor(static_cast<float>(j)/BlockCols);
            const int joffset = j - pj*BlockCols;
            ASSERT(pi < this->m_row_particles.size(),"pi greater than a.size()");
            ASSERT(pj < this->m_col_particles.size(),"pj greater than b.size()");
            const_row_reference ai = this->m_row_particles[pi];
            const_col_reference bj = this->m_col_particles[pj];
            const_position_reference dx = get<position>(bj)-get<position>(ai);
            if (dx.squaredNorm() < std::pow(m_radius_function(ai),2)) {
                return this->m_function(dx,ai,bj)(ioffset,joffset);
            } else {
                return 0.0;
            }
        }

        template<typename MatrixType>
        void assemble(const MatrixType &matrix) const {

            const RowElements& a = this->m_row_particles;
            const ColElements& b = this->m_col_particles;

            const size_t na = a.size();
            const size_t nb = b.size();

            const_cast< MatrixType& >(matrix).setZero();

            //sparse a x b block
            for (size_t i=0; i<na; ++i) {
                const_row_reference ai = a[i];
                const double radius = m_radius_function(ai);
                for (auto pairj: euclidean_search(b.get_query(),get<position>(ai),radius)) {
                    const_col_reference bj = detail::get_impl<0>(pairj);
                    const_position_reference dx = detail::get_impl<1>(pairj);
                    const size_t j = &get<position>(bj) - get<position>(b).data();
                    const_cast< MatrixType& >(matrix).block<BlockRows,BlockCols>(
                            i*BlockRows,j*BlockCols) = 
                                                Block(this->m_function(dx,ai,bj));

                }
            }
        }

        template<typename Triplet>
        void assemble(std::vector<Triplet>& triplets,
                      const size_t startI=0, const size_t startJ=0
                      ) const {

            const RowElements& a = this->m_row_particles;
            const ColElements& b = this->m_col_particles;

            const size_t na = a.size();
            const size_t nb = b.size();

            //sparse a x b block
            //std::cout << "sparse a x b block" << std::endl;
            for (size_t i=0; i<na; ++i) {
                const_row_reference ai = a[i];
                const double radius = m_radius_function(ai);
                for (auto pairj: euclidean_search(b.get_query(),get<position>(ai),radius)) {
                    const_col_reference bj = detail::get_impl<0>(pairj);
                    const_position_reference dx = detail::get_impl<1>(pairj);
                    const size_t j = &get<position>(bj) - get<position>(b).data();
                    const Block element = this->m_function(dx,ai,bj);
                    for (int ii = 0; ii < BlockRows; ++ii) {
                        for (int jj = 0; jj < BlockCols; ++jj) {
                            triplets.push_back(Triplet(i*BlockRows+ii+startI,
                                                       j*BlockCols+jj+startJ,
                                                       element(ii,jj)));
                        }
                    }
                }
            }
        }

        /// Evaluates a matrix-free linear operator given by \p expr \p if_expr,
        /// and particle sets \p a and \p b on a vector rhs and
        /// accumulates the result in vector lhs
        ///
        template <typename LHSType, typename RHSType>
        void evaluate(std::vector<LHSType> &lhs, 
                const std::vector<RHSType> &rhs) const {

            const RowElements& a = this->m_row_particles;
            const ColElements& b = this->m_col_particles;

            const size_t na = a.size();
            const size_t nb = b.size();

            ASSERT(na == a.size(),"lhs vector has incompatible size");
            ASSERT(nb == b.size(),"rhs vector has incompatible size");

            #pragma omp parallel for
            for (size_t i=0; i<na; ++i) {
                const_row_reference ai = a[i];
                const double radius = m_radius_function(ai);
                for (auto pairj: euclidean_search(b.get_query(),get<position>(ai),radius)) {
                    const_position_reference dx = detail::get_impl<1>(pairj);
                    const_col_reference bj = detail::get_impl<0>(pairj);
                    const size_t j = &get<position>(bj) - get<position>(b).data();
                    lhs[i] += this->m_function(dx,ai,bj)*rhs[j];
                }
            }
       }

        template <typename DerivedLHS, typename DerivedRHS>
        void evaluate(Eigen::DenseBase<DerivedLHS> &lhs, 
                const Eigen::DenseBase<DerivedRHS> &rhs) const {
            typedef Eigen::Matrix<Scalar,BlockRows,1> row_vector;

            ASSERT(this->rows() == lhs.rows(),"lhs vector has incompatible size");
            ASSERT(this->cols() == rhs.rows(),"rhs vector has incompatible size");

            const RowElements& a = this->m_row_particles;
            const ColElements& b = this->m_col_particles;

            const size_t na = a.size();
            const size_t nb = b.size();

            #pragma omp parallel for
            for (size_t i=0; i<na; ++i) {
                const_row_reference ai = a[i];
                const double radius = m_radius_function(ai);
                for (auto pairj: euclidean_search(b.get_query(),get<position>(ai),radius)) {
                    const_position_reference dx = detail::get_impl<1>(pairj);
                    const_col_reference bj = detail::get_impl<0>(pairj);
                    const size_t j = &get<position>(bj) - get<position>(b).data();
                    lhs.template segment<BlockRows>(i*BlockRows) += 
                        this->m_function(dx,ai,bj)
                            *rhs.template segment<BlockCols>(j*BlockCols);
                }
            }
       }

    private:
        FRadius m_radius_function;
    };

    namespace detail {
        template <typename RowElements,
                 typename const_row_reference=
                     typename RowElements::const_reference>
        struct constant_radius {
            const double m_radius;
            constant_radius(const double radius):m_radius(radius)
            {}
            double operator()(const_row_reference a) const {
                return m_radius; 
            }
        };
    }

    template<typename RowElements, typename ColElements, typename F,
        typename RadiusFunction=detail::constant_radius<RowElements>>
    class KernelSparseConst: public KernelSparse<RowElements,ColElements,RadiusFunction,F> {
        typedef KernelSparse<RowElements,ColElements,RadiusFunction,F> base_type;
        typedef typename base_type::position position;
        typedef typename base_type::double_d double_d;
        typedef typename base_type::const_position_reference const_position_reference;
        typedef typename base_type::const_row_reference const_row_reference;
        typedef typename base_type::const_col_reference const_col_reference;
    public:
        typedef typename base_type::Block Block;
        typedef typename base_type::Scalar Scalar;
        static const size_t BlockRows = base_type::BlockRows;
        static const size_t BlockCols = base_type::BlockCols;

        KernelSparseConst(const RowElements& row_particles,
                    const ColElements& col_particles,
                    const double radius,
                    const F& function): base_type(row_particles,
                                                  col_particles,
                                                  RadiusFunction(radius),
                                                  function) 
        {}
    };

    namespace detail {
        template <typename RowElements, typename ColElements,
                 typename double_d=
                     typename RowElements::double_d,
                 typename const_row_reference=
                     typename RowElements::const_reference,
                 typename const_col_reference=
                     typename ColElements::const_reference>
        struct zero_lambda {
            double operator()(const double_d& dx,
                            const_row_reference a,
                            const_col_reference b) const {
                                return 0.0;
            }
        };
    }

    template<typename RowElements, typename ColElements,
        typename F=detail::zero_lambda<RowElements,ColElements>>
    class KernelZero: public KernelBase<RowElements,ColElements,F> {

        typedef KernelBase<RowElements,ColElements,F> base_type;
        typedef typename base_type::position position;
        typedef typename base_type::double_d double_d;
        typedef typename base_type::const_position_reference const_position_reference;
        typedef typename base_type::const_row_reference const_row_reference;
        typedef typename base_type::const_col_reference const_col_reference;

    public:
        typedef typename base_type::Block Block;
        typedef typename base_type::Scalar Scalar;

        KernelZero(const RowElements& row_particles,
                   const ColElements& col_particles): 
            base_type(row_particles,
                      col_particles,
                      F()) 
        {};

        template<typename MatrixType>
        void assemble(const MatrixType &matrix) const {
            const_cast< MatrixType& >(matrix).setZero();
        }
        

        template<typename Triplet>
        void assemble(std::vector<Triplet>& triplets,
                      const size_t startI=0, const size_t startJ=0
                      ) const {
        }

        /// Evaluates a matrix-free linear operator given by \p expr \p if_expr,
        /// and particle sets \p a and \p b on a vector rhs and
        /// accumulates the result in vector lhs
        template<typename VectorLHS,typename VectorRHS>
        void evaluate(VectorLHS &lhs, const VectorRHS &rhs) const {
        }
    };




}

#endif
