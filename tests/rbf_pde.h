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


#ifndef RBF_PDE_TEST_H_
#define RBF_PDE_TEST_H_


#include <cxxtest/TestSuite.h>

//[rbf_pde
#include "Aboria.h"

using namespace Aboria;

//<-
class RbfPdeTest : public CxxTest::TestSuite {
public:

    template<template <typename> class SearchMethod>
    void helper_Eigen(void) {
#ifdef HAVE_EIGEN
//->
//=int main() {
        auto funct = [](const double x, const double y) { 
            return std::cos(4*x+4*y);
        };

        auto laplace_funct = [](const double x, const double y) { 
            return -32*std::cos(4*x+4*y);
        };

        ABORIA_VARIABLE(boundary,uint8_t,"is boundary knot")
        ABORIA_VARIABLE(interpolated,double,"interpolated value")
        ABORIA_VARIABLE(constant2,double,"c2 value")
        ABORIA_VARIABLE(alpha,double,"alpha value")

//<-
    	typedef Particles<std::tuple<alpha,boundary,constant2,interpolated>,2,std::vector,SearchMethod> ParticlesType;
//->
//=     typedef Particles<std::tuple<alpha,boundary,constant2,interpolated>,2> ParticlesType;

        typedef position_d<2> position;
        typedef typename ParticlesType::const_reference const_particle_reference;
        typedef Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>> map_type; 
        typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type; 
       	ParticlesType knots,augment;

       	const double c = 0.5;
        const int max_iter = 100;
        const int restart = 100;
        
        const int nx = 7;
        constexpr int N = (nx+1)*(nx+1);
        const double delta = 1.0/nx;
        typename ParticlesType::value_type p;
        for (int i=0; i<=nx; ++i) {
            for (int j=0; j<=nx; ++j) {
                get<position>(p) = vdouble2(i*delta,j*delta);
                if ((i==0)||(i==nx)||(j==0)||(j==nx)) {
                    get<boundary>(p) = true;
                } else {
                    get<boundary>(p) = false;
                }
                get<constant2>(p) = std::pow(c,2);
                knots.push_back(p);
            }
        }
        augment.push_back(p);

        auto kernel = [](
                         const_particle_reference a,
                         const_particle_reference b) {
             const vdouble2 dx = get<position>(b) - get<position>(a);
             return std::exp(-dx.squaredNorm()/get<constant2>(b));
                        };

        auto laplace_kernel = [](
                         const_particle_reference a,
                         const_particle_reference b) {
             const vdouble2 dx = get<position>(b) - get<position>(a);
             return 4.0*(dx.squaredNorm()-get<constant2>(b))*
                        std::exp(-dx.squaredNorm()/get<constant2>(b))/
                            (get<constant2>(b)*get<constant2>(b));
                        };

        auto G = create_dense_operator(knots,knots,
                [kernel,laplace_kernel](
                         const_particle_reference a,
                         const_particle_reference b) {
                    if (get<boundary>(a)) {
                        return kernel(a,b);
                    } else {
                        return laplace_kernel(a,b);
                    }
                    });

        auto P = create_dense_operator(knots,augment,
                [](
                         const_particle_reference a,
                         const_particle_reference b) {
                    if (get<boundary>(a)) {
                        return 1.0;
                    } else {
                        return 0.0;
                    }
                    });

        auto Pt = create_dense_operator(augment,knots,
                [](
                         const_particle_reference a,
                         const_particle_reference b) {
                    if (get<boundary>(b)) {
                        return 1.0;
                    } else {
                        return 0.0;
                    }
                    });

        auto Zero = create_zero_operator(augment,augment);

        auto W = create_block_operator<2,2>(G, P,
                                            Pt,Zero);


        vector_type phi(N+1), gamma(N+1);
        for (int i=0; i<N; ++i) {
            const double x = get<position>(knots[i])[0];
            const double y = get<position>(knots[i])[1];
            if (get<boundary>(knots[i])) {
                phi[i] = funct(x,y);
            } else {
                phi[i] = laplace_funct(x,y);
            }
        }
        phi[N] = 0;

        std::cout << std::endl;

        Eigen::GMRES<decltype(W), Eigen::DiagonalPreconditioner<double>> gmres;
        gmres.set_restart(restart);
        gmres.setMaxIterations(max_iter);
        gmres.compute(W);
        gamma = gmres.solve(phi);
        std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;

        phi = W*gamma;
        for (int i=0; i<N; ++i) {
            const double x = get<position>(knots[i])[0];
            const double y = get<position>(knots[i])[1];
            if (get<boundary>(knots[i])) {
                TS_ASSERT_DELTA(phi[i],funct(x,y),2e-3); 
            } else {
                TS_ASSERT_DELTA(phi[i],laplace_funct(x,y),2e-3); 
            }
        }
        TS_ASSERT_DELTA(phi[N],0,2e-3); 


        map_type alpha_wrap(get<alpha>(knots).data(),N);
        map_type interp_wrap(get<interpolated>(knots).data(),N);
        alpha_wrap = gamma.segment<N>(0);

        vector_type beta = vector_type::Constant(N,gamma[N]);

        auto K = create_dense_operator(knots,knots,kernel);
        interp_wrap = K*alpha_wrap + beta;

        double rms_error = 0;
        double scale = 0;
        for (int i=0; i<N; ++i) {
            const double x = get<position>(knots[i])[0];
            const double y = get<position>(knots[i])[1];
            const double truth = funct(x,y);
            const double eval_value = get<interpolated>(knots[i]);
            rms_error += std::pow(eval_value-truth,2);
            scale += std::pow(truth,2);
            TS_ASSERT_DELTA(eval_value,truth,1e-2); 
        }
        std::cout << "rms_error for global support, at centers  = "<<std::sqrt(rms_error/scale)<<std::endl;
        TS_ASSERT_LESS_THAN(std::sqrt(rms_error/scale),1e-3);

//=}
//]
#endif // HAVE_EIGEN
    }


    void test_CellListOrdered() {
        helper_Eigen<CellListOrdered>();
    }

    void test_CellList() {
        helper_Eigen<CellList>();
    }


};

#endif /* RBF_PDE_TEST_H_ */
