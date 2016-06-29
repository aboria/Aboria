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
        ABORIA_VARIABLE(boundary,uint8_t,"is boundary knot")
        ABORIA_VARIABLE(truth,double,"truth value")
        ABORIA_VARIABLE(forcing,double,"forcing value")
        ABORIA_VARIABLE(estimated,double,"estimated value")

    	typedef Particles<std::tuple<boundary,truth,forcing,estimated>,2> ParticlesType;
        typedef position_d<2> position;
       	ParticlesType knots;

       	const double c2 = std::pow(0.11,2);
        double2 min(0);
        double2 max(1);
        double2 periodic(false);
        
        const int nx = 10;
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
                knots.push_back(p);
            }
        }

        Symbol<boundary> is_b;
        Symbol<position> r;
        Label<0,ParticlesType> a(knots);
        Label<1,ParticlesType> b(knots);
        One one;
        auto dx = create_dx(a,b);
        Accumulate<std::plus<double> > sum;

        auto H = create_eigen_operator(a,b, 
                    sqrt(pow(norm(dx),2) + c2)
                );
        auto G = create_eigen_operator(a,b, 
                    if_else(is_b[a],
                        sqrt(pow(norm(dx),2) + c2),
                        (2*c2 + pow(norm(dx),2)) / pow(pow(norm(dx),2) + c2,1.5)
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
            const double F = std::cos(4*x+4*y);
            if (get<boundary>(knots[i])) {
                phi[i] = F;
            } else {
                phi[i] = -32*F;
            }
        }
        phi[knots.size()] = 0;

        Eigen::ConjugateGradient<decltype(W), 
            Eigen::Lower|Eigen::Upper,  Eigen::DiagonalPreconditioner<double>> cg;
        cg.compute(W);
        gamma = cg.solve(phi);
        std::cout << std::endl << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;

        const double beta = gamma(knots.size());
        for (int i=0; i<knots.size(); ++i) {
            const double2 p = get<position>(knots[i]);
            const double truth = std::cos(4*p[0]+4*p[1]);
            const double gamma_i = gamma(i);
            const double eval_value = eval(sum(a,true,gamma_i*sqrt(pow(norm(r[a]-p),2) + c2))) + beta;
            TS_ASSERT_DELTA(eval_value,truth,2e-3); 
        }

        
        
#endif // HAVE_EIGEN
    }


};

#endif /* RBF_PDE_TEST_H_ */
