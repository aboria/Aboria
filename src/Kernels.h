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

namespace Aboria {

    template<typename RowParticles, typename ColParticles, typename F>
    class KernelBase {
    protected:
        typedef typename RowParticles::position position;
        typedef typename RowParticles::double_d double_d;
        typedef typename position::value_type const & const_position_reference;
        typedef typename RowParticles::const_reference const_row_reference;
        typedef typename ColParticles::const_reference const_col_reference;
    public:
        typedef typename std::result_of<F(const_position_reference, 
                                  const_row_reference, 
                                  const_col_reference)>::type Scalar;

        KernelBase(const RowParticles& row_particles,
                   const ColParticles& col_particles,
                   const F& function): 
            m_function(function),
            m_row_particles(row_particles), 
            m_col_particles(col_particles)
        {};


        size_t size_row() const {
            return m_row_particles.size();
        }

        size_t size_col() const {
            return m_row_particles.size();
        }

        Scalar eval(const_position_reference dx, 
                    const_row_reference a, 
                    const_col_reference b) const {
            return m_function(dx,a,b);
        }

        Scalar coeff(const size_t i, const size_t j) const {
            ASSERT(i>=0, "i less than zero");
            ASSERT(j>=0, "j less than zero");
            ASSERT(i < m_row_particles.size(),"i greater than a.size()");
            ASSERT(j < m_col_particles.size(),"j greater than b.size()");
            const_row_reference ai = m_row_particles[i];
            const_col_reference bj = m_col_particles[i];
            const_position_reference dx = get<position>(bj)-get<position>(ai);
            return eval(dx,ai,bj);
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
        const RowParticles &m_row_particles;
        const ColParticles &m_col_particles;
        F m_function;
    };

    template<typename RowParticles, typename ColParticles, typename F>
    class KernelDense: public KernelBase<RowParticles,ColParticles,F> {
        typedef KernelBase<RowParticles,ColParticles,F> base_type;
        typedef typename base_type::position position;
        typedef typename base_type::double_d double_d;
        typedef typename base_type::const_position_reference const_position_reference;
        typedef typename base_type::const_row_reference const_row_reference;
        typedef typename base_type::const_col_reference const_col_reference;
    public:
        typename base_type::Scalar Scalar;

        KernelDense(const RowParticles& row_particles,
                    const ColParticles& col_particles,
                    const F& function): base_type(row_particles,
                                                  col_particles,
                                                  function) 
        {};

        template<typename MatrixType>
        void assemble(const MatrixType &matrix) const {

            const RowParticles& a = this->m_row_particles;
            const ColParticles& b = this->m_col_particles;
            const size_t na = a.size();
            const size_t nb = b.size();

            ASSERT(!a.get_periodic().any(),"periodic does not work with dense");

            for (size_t i=0; i<na; ++i) {
                const_row_reference ai = a[i];
                for (size_t j=0; j<nb; ++j) {
                    const_col_reference  bj = b[j];
                    const_position_reference dx = 
                        get<position>(bj)-get<position>(ai);
                    const_cast< MatrixType& >(matrix)(i,j) = this->eval(dx,ai,bj);
                }
            }
        }

        template<typename Triplet>
        void assemble(std::vector<Triplet>& triplets,
                      const size_t startI=0, const size_t startJ=0
                      ) const {

            const RowParticles& a = this->m_row_particles;
            const ColParticles& b = this->m_col_particles;

            const size_t na = a.size();
            const size_t nb = b.size();

            ASSERT(!a.get_periodic().any(),"periodic does not work with dense");

            for (size_t i=0; i<na; ++i) {
                const_row_reference ai = a[i];
                for (size_t j=0; j<nb; ++j) {
                    const_col_reference bj = b[j];
                    const_position_reference dx = 
                        get<position>(bj)-get<position>(ai);
                    triplets.push_back(Triplet(i+startI,j+startJ,
                                this->eval(dx,ai,bj)));
                }
            }
        }

        /// Evaluates a matrix-free linear operator given by \p expr \p if_expr,
        /// and particle sets \p a and \p b on a vector rhs and
        /// accumulates the result in vector lhs
        template<typename VectorLHS,typename VectorRHS>
        void evaluate(VectorLHS &lhs, const VectorRHS &rhs) const {

            const RowParticles& a = this->m_row_particles;
            const ColParticles& b = this->m_col_particles;

            const size_t na = a.size();
            const size_t nb = b.size();

            ASSERT(!a.get_periodic().any(),"periodic does not work with dense");

            const size_t parallel_size = 20;
            if (na > parallel_size) {
                #pragma omp parallel for
                for (size_t i=0; i<na; ++i) {
                    const_row_reference ai = a[i];
                    double sum = 0;
                    for (size_t j=0; j<nb; ++j) {
                        const_col_reference bj = b[j];
                        const_position_reference dx = 
                            get<position>(bj)-get<position>(ai);
                        sum += this->eval(dx,ai,bj)*rhs(j);
                    }
                    lhs[i] += sum;
                }
            } else {
                for (size_t i=0; i<na; ++i) {
                    const_row_reference ai = a[i];
                    double sum = 0;
                    for (size_t j=0; j<nb; ++j) {
                        const_col_reference bj = b[j];
                        const_position_reference dx = 
                            get<position>(bj)-get<position>(ai);
                        sum += this->eval(dx,ai,bj)*rhs[j];
                    }
                    lhs[i] += sum;
                }
            }
       }
    };

    namespace detail {
        template <typename RowParticles, typename ColParticles,
                 typename double_d=
                     typename RowParticles::double_d,
                 typename const_row_reference=
                     typename RowParticles::const_reference,
                 typename const_col_reference=
                     typename ColParticles::const_reference>
        struct zero_lambda {
            double operator()(const double_d& dx,
                            const_row_reference a,
                            const_col_reference b) const {
                                return 0.0;
            }
        };
    }

    template<typename RowParticles, typename ColParticles,
        typename F=detail::zero_lambda<RowParticles,ColParticles>>
    class KernelZero: public KernelBase<RowParticles,ColParticles,F> {

        typedef KernelBase<RowParticles,ColParticles,F> base_type;
        typedef typename base_type::position position;
        typedef typename base_type::double_d double_d;
        typedef typename base_type::const_position_reference const_position_reference;
        typedef typename base_type::const_row_reference const_row_reference;
        typedef typename base_type::const_col_reference const_col_reference;

    public:
        typename base_type::Scalar Scalar;

        KernelZero(const RowParticles& row_particles,
                    const ColParticles& col_particles): 
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
