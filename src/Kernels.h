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
        typedef typename RowParticles::position position;
        typedef typename RowParticles::double_d double_d;
        typedef const typename position::value_type& const_position_reference;
        typedef typename RowParticles::const_reference const_row_reference;
        typedef typename ColParticles::const_reference const_col_reference;
    public:
        typedef typename std::result_of<F(const_position_reference, 
                                  const_row_reference, 
                                  const_col_reference)>::type Scalar;

        KernelBase(const F& function): m_function(function) {};

        Scalar eval(const_position_reference dx, 
                    const_row_reference a, 
                    const_col_reference b) {
            return m_function(dx,a,b);
        }

        template<typename MatrixType>
        void assemble(
                const RowParticles &a, 
                const ColParticles &b, 
                const MatrixType &matrix) {
                    MatrixType::MATRIX_ASSEMBLE_NOT_IMPLEMENTED; 
        }

        template<typename Triplet>
        void assemble(const RowParticles &a,
                      const ColParticles &b,
                      std::vector<Triplet>& triplets,
                      const size_t startI=0, const size_t startJ=0
                      ) {
            Triplet::TRIPLET_ASSEMBLE_NOT_IMPLEMENTED; 
        }

        /// Evaluates a matrix-free linear operator given by \p expr \p if_expr,
        /// and particle sets \p a and \p b on a vector rhs and
        /// accumulates the result in vector lhs
        template<typename VectorLHS,typename VectorRHS>
        void evaluate(
                const RowParticles &a, 
                const ColParticles &b, 
                VectorLHS &lhs, const VectorRHS &rhs) {

            VectorLHS::EVALUATE_NOT_IMPLEMENTED;
       }

    private:
        F m_function;
    };

    template<typename RowParticles, typename ColParticles, typename F>
    class KernelDense: public KernelBase<RowParticles,ColParticles,F> {
        typedef KernelBase<RowParticles,ColParticles,F> base_type;
        typedef typename base_type::position position;
        typedef typename base_type::double_d double_d;
        typedef typename base_type::const_position_reference const_position_reference;
        typedef typename base_type::const_reference const_row_reference;
        typedef typename base_type::const_reference const_col_reference;
    public:
        typename base_type::Scalar Scalar;

        KernelDense(const F& function): base_type(function) {};

        template<typename MatrixType>
        void assemble(
                const RowParticles &a, 
                const ColParticles &b, 
                const MatrixType &matrix) {

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
        void assemble(const RowParticles &a,
                      const ColParticles &b,
                      std::vector<Triplet>& triplets,
                      const size_t startI=0, const size_t startJ=0
                      ) {

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
        void evaluate(
                const RowParticles &a, 
                const ColParticles &b, 
                VectorLHS &lhs, const VectorRHS &rhs) {

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



}

#endif
