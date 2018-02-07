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

#ifndef EVALUATE_H_
#define EVALUATE_H_

#include "Symbolic.h"
#include "detail/Evaluate.h"

namespace Aboria {

/// Evaluates a non-linear operator \p expr over a set of particles
/// given by label \p label and stores the result, using the functor
/// \p Functor, in variable with type \p VariableType
template <typename VariableType, typename Functor, typename ExprRHS,
          typename LabelType>
void evaluate_nonlinear(ExprRHS const &expr, LabelType &label) {
  typedef typename VariableType::value_type value_type;
  typedef typename LabelType::particles_type particles_type;
  typedef typename particles_type::position position;

  typedef
      typename proto::matches<ExprRHS,
                              detail::is_not_aliased<VariableType, LabelType>>
          not_aliased;
  typedef typename LabelType::particles_type particles_type;

  particles_type &particles = label.get_particles();

  // check expr is a univariate expression and that it refers to the same
  // particles container
  check_valid_assign_expr(label, expr);

  // if aliased then need to copy to a tempory buffer first
  std::vector<value_type> &buffer =
      (not_aliased::value) ? get<VariableType>(particles)
                           : get<VariableType>(label.get_buffers());
  buffer.resize(particles.size());

  // evaluate expression for all particles and store in buffer
  const size_t n = particles.size();
  Functor functor;
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
  for (size_t i = 0; i < n; i++) {
    buffer[i] =
        functor(get<VariableType>(particles)[i], eval(expr, particles[i]));
  }

  // if aliased then copy back from the buffer
  if (not_aliased::value == false) {
    const size_t n = particles.size();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < n; i++) {
      get<VariableType>(particles[i]) = buffer[i];
    }
  }

  if (boost::is_same<VariableType, position>::value) {
    particles.update_positions();
  }

  if (boost::is_same<VariableType, alive>::value) {
    particles.update_positions();
  }
}

/*
/// Evaluates a matrix-free linear operator given by \p expr \p if_expr,
/// and particle sets \p a and \p b on a vector rhs and
/// accumulates the result in vector lhs
template<typename Expr,
         typename IfExpr,
         typename ParticlesTypeA,
         typename ParticlesTypeB,
         typename VectorLHS,
         typename VectorRHS
         >
void evaluate_linear(Expr &expr,
              IfExpr &if_expr,
        const ParticlesTypeA &a,
        const ParticlesTypeB &b,
        VectorLHS &lhs, const VectorRHS &rhs) {

    typedef typename ParticlesTypeA::double_d double_d;
    typedef typename ParticlesTypeA::position position;

    const size_t na = a.size();
    const size_t nb = b.size();

    if (is_trivially_zero(expr) || is_trivially_false(if_expr)) {
        //std::cout << "zero a x b block" <<std::endl;
        return;
    }

    if (is_trivially_true(if_expr)) {
        //std::cout << "dense "<<na<<" x "<<nb<<" block" <<std::endl;
        ASSERT(!a.get_periodic().any(),"periodic does not work with dense");

        const size_t parallel_size = 20;
        const size_t block_size = 20;
        if (na > parallel_size) {
            #pragma omp parallel for
            for (size_t i=0; i<na; ++i) {
                typename ParticlesTypeA::const_reference ai = a[i];
                double sum = 0;
                for (size_t j=0; j<nb; ++j) {
                    typename ParticlesTypeB::const_reference bj = b[j];
                    sum +=
eval(expr,get<position>(bj)-get<position>(ai),ai,bj)*rhs(j);
                }
                lhs[i] += sum;
            }
        } else {
            for (size_t i=0; i<na; ++i) {
                typename ParticlesTypeA::const_reference ai = a[i];
                double sum = 0;
                for (size_t j=0; j<nb; ++j) {
                    typename ParticlesTypeB::const_reference bj = b[j];
                    //std::cout << "a = "<<get<position>(ai)<<" b =
"<<get<position>(bj)<<std::endl;
                    //std::cout << "using dx =
"<<get<position>(bj)-get<position>(ai)<<" rhs(j) = "<<rhs(j)<<" eval =
"<<eval(expr,get<position>(bj)-get<position>(ai),ai,bj)<<std::endl; sum +=
eval(expr,get<position>(bj)-get<position>(ai),ai,bj)*rhs[j];
                }
                lhs[i] += sum;
            }
        }
    } else {
        //std::cout << "sparse a x b block" <<std::endl;
        #pragma omp parallel for
        for (size_t i=0; i<na; ++i) {
            typename ParticlesTypeA::const_reference ai = a[i];
            double sum = 0;
            //std::cout << "evaluating fucntion for particle at
"<<get<position>(ai)<<std::endl; for (auto pairj:
box_search(b.get_query(),get<position>(ai))) { const double_d & dx =
tuple_ns::get<1>(pairj); typename ParticlesTypeB::const_reference bj =
tuple_ns::get<0>(pairj);
                //std::cout << "looking at particle with dx = "<<dx<<std::endl;
                const size_t j = &get<position>(bj) - get<position>(b).data();
                if (eval(if_expr,dx,ai,bj)) {
                    //std::cout <<"if expression is true. eval =
"<<eval(expr,dx,ai,bj)<<std::endl; sum += eval(expr,dx,ai,bj)*rhs(j);
                }
            }
            lhs[i] += sum;
        }
    }
}


template<typename Expr,
         typename IfExpr,
         typename ParticlesTypeA,
         typename ParticlesTypeB,
         typename Triplet
         >
void assemble(Expr &expr,
        IfExpr &if_expr,
        const ParticlesTypeA &a,
        const ParticlesTypeB &b,
        std::vector<Triplet>& triplets,
        const size_t startI=0, const size_t startJ=0) {

    typedef typename ParticlesTypeB::double_d double_d;
    typedef typename ParticlesTypeB::position position;

    const size_t na = a.size();
    const size_t nb = b.size();

    if (is_trivially_zero(expr) || is_trivially_false(if_expr)) {
        //zero a x b block
        return;
    }

    if (is_trivially_true(if_expr)) {
        //dense a x b block
        //std::cout << "dense a x b block" << std::endl;
        ASSERT(!a.get_periodic().any(),"periodic does not work with dense");

        for (size_t i=0; i<na; ++i) {
            typename ParticlesTypeA::const_reference ai = a[i];
            for (size_t j=0; j<nb; ++j) {
                typename ParticlesTypeB::const_reference bj = b[j];
                triplets.push_back(Triplet(i+startI,j+startJ,eval(expr,get<position>(bj)-get<position>(ai),ai,bj)));
            }
        }
    } else {
        //sparse a x b block
        //std::cout << "sparse a x b block" << std::endl;
        for (size_t i=0; i<na; ++i) {
            typename ParticlesTypeA::const_reference ai = a[i];
            for (auto pairj: box_search(b.get_query(),get<position>(ai))) {
                const double_d & dx = tuple_ns::get<1>(pairj);
                typename ParticlesTypeB::const_reference bj =
tuple_ns::get<0>(pairj); const size_t j = &get<position>(bj) -
get<position>(b).data(); if (eval(if_expr,dx,ai,bj)) {
                    triplets.push_back(Triplet(i+startI,j+startJ,eval(expr,dx,ai,bj)));
                }
            }
        }
    }
}

template<typename Expr,
         typename IfExpr,
         typename ParticlesTypeA,
         typename ParticlesTypeB,
         typename MatrixType
         >
void assemble(Expr &expr,
        IfExpr &if_expr,
        const ParticlesTypeA &a,
        const ParticlesTypeB &b,
        const MatrixType &matrix) {

    typedef typename ParticlesTypeB::double_d double_d;
    typedef typename ParticlesTypeB::position position;

    const size_t na = a.size();
    const size_t nb = b.size();

    if (is_trivially_zero(expr) || is_trivially_false(if_expr)) {
        //zero a x b block
        // hack from
https://eigen.tuxfamily.org/dox/TopicFunctionTakingEigenTypes.html const_cast<
MatrixType& >(matrix).setZero(); return;
    }

    if (is_trivially_true(if_expr)) {
        //dense a x b block
        ASSERT(!a.get_periodic().any(),"periodic does not work with dense");

        for (size_t i=0; i<na; ++i) {
            typename ParticlesTypeA::const_reference ai = a[i];
            for (size_t j=0; j<nb; ++j) {
                typename ParticlesTypeB::const_reference bj = b[j];
                const double_d dx = get<position>(bj)-get<position>(ai);
                const_cast< MatrixType& >(matrix)(i,j) = eval(expr,dx,ai,bj);
            }
        }
    } else {
        //sparse a x b block
        for (size_t i=0; i<na; ++i) {
            typename ParticlesTypeA::const_reference ai = a[i];
            for (auto pairj: box_search(b.get_query(),get<position>(ai))) {
                const double_d & dx = tuple_ns::get<1>(pairj);
                typename ParticlesTypeB::const_reference bj =
tuple_ns::get<0>(pairj); const size_t j = &get<position>(bj) -
get<position>(b).data(); if (eval(if_expr,dx,ai,bj)) { const_cast< MatrixType&
>(matrix)(i,j) = eval(expr,dx,ai,bj); } else { const_cast< MatrixType&
>(matrix)(i,j) = 0;
                }
            }
        }
    }
}
*/

} // namespace Aboria
#endif
