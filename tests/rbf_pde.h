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

#define LOG_LEVEL 1
#include "Aboria.h"

using namespace Aboria;


class OperatorsTest : public CxxTest::TestSuite {
public:

    void test_Eigen(void) {
#ifdef HAVE_EIGEN
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

    	typedef Particles<std::tuple<alpha,boundary,constant2,interpolated>,2> ParticlesType;
        typedef position_d<2> position;
       	ParticlesType knots;

       	const double c = 0.5;
        const int max_iter = 100;
        const int restart = 100;
        double2 periodic(false);
        
        const int nx = 7;
        constexpr int N = (nx+1)*(nx+1);
        const double delta = 1.0/nx;
        ParticlesType::value_type p;
        for (int i=0; i<=nx; ++i) {
            for (int j=0; j<=nx; ++j) {
                get<position>(p) = double2(i*delta,j*delta);
                if ((i==0)||(i==nx)||(j==0)||(j==nx)) {
                    get<boundary>(p) = true;
                } else {
                    get<boundary>(p) = false;
                }
                get<constant2>(p) = std::pow(c,2);
                knots.push_back(p);
            }
        }

        Symbol<boundary> is_b;
        Symbol<position> r;
        Symbol<constant2> c2;
        Symbol<alpha> al;
        Symbol<interpolated> interp;
        Label<0,ParticlesType> a(knots);
        Label<1,ParticlesType> b(knots);
        One one;
        auto dx = create_dx(a,b);
        Accumulate<std::plus<double> > sum;

        auto kernel = deep_copy(
                exp(-pow(norm(dx),2)/c2[b])
                //sqrt(pow(norm(dx),2) + c2[b])
                );

        auto laplace_kernel = deep_copy(
                //(2*c2[b] + pow(norm(dx),2)) / pow(pow(norm(dx),2) + c2[b],1.5)
                4*(pow(norm(dx),2) - c2[b]) * exp(-pow(norm(dx),2)/c2[b])/pow(c2[a],2)
                );

        auto G = create_eigen_operator(a,b, 
                    if_else(is_b[a],
                        kernel,
                        laplace_kernel
                    )
                );
        auto P = create_eigen_operator(a,one,
                    if_else(is_b[a],
                        1.0,
                        0.0
                    )
                );
        auto Pt = create_eigen_operator(one,b,
                    if_else(is_b[b],
                        1.0,
                        0.0
                    )
                );
        auto Zero = create_eigen_operator(one,one, 0.);

        auto W = create_block_eigen_operator<2,2>(G, P,
                                                  Pt,Zero);


        Eigen::VectorXd phi(knots.size()+1), gamma;
        for (int i=0; i<knots.size(); ++i) {
            const double x = get<position>(knots[i])[0];
            const double y = get<position>(knots[i])[1];
            if (get<boundary>(knots[i])) {
                phi[i] = funct(x,y);
            } else {
                phi[i] = laplace_funct(x,y);
            }
        }
        phi[knots.size()] = 0;

        std::cout << std::endl;
       

        Eigen::GMRES<decltype(W), Eigen::DiagonalPreconditioner<double>> gmres;
        gmres.set_restart(restart);
        gmres.setMaxIterations(max_iter);
        gmres.compute(W);
        gamma = gmres.solve(phi);
        std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;


        /*
        Eigen::DGMRES<decltype(W), Eigen::DiagonalPreconditioner<double>> dgmres;
        dgmres.setMaxIterations(max_iter);
        dgmres.compute(W);
        gamma = dgmres.solve(phi);
        std::cout << "DGMRES:   #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;
        */

        phi = W*gamma;
        for (int i=0; i<knots.size(); ++i) {
            const double x = get<position>(knots[i])[0];
            const double y = get<position>(knots[i])[1];
            if (get<boundary>(knots[i])) {
                TS_ASSERT_DELTA(phi[i],funct(x,y),2e-3); 
            } else {
                TS_ASSERT_DELTA(phi[i],laplace_funct(x,y),2e-3); 
            }
        }
        TS_ASSERT_DELTA(phi[knots.size()],0,2e-3); 


        // This could be more intuitive....
        Eigen::Map<Eigen::Matrix<double,N,1>> alpha_wrap(get<alpha>(knots).data());
        alpha_wrap = gamma.segment<N>(0);

        const double beta = gamma(knots.size());

        interp[a] = sum(b,true,al[b]*kernel) + beta;

        for (int i=0; i<knots.size(); ++i) {
            const double x = get<position>(knots[i])[0];
            const double y = get<position>(knots[i])[1];
            const double truth = funct(x,y);
            const double eval_value = get<interpolated>(knots[i]);
            //TODO: bad point error, can we improve?
            TS_ASSERT_DELTA(eval_value,truth,1e-2); 
        }
#endif // HAVE_EIGEN
    }


};

#endif /* RBF_PDE_TEST_H_ */
