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


#ifndef OPERATORSTEST_H_
#define OPERATORSTEST_H_


#include <cxxtest/TestSuite.h>

#define LOG_LEVEL 1
#include "Aboria.h"

using namespace Aboria;


class OperatorsTest : public CxxTest::TestSuite {
public:

    void test_Eigen(void) {
#ifdef HAVE_EIGEN
        ABORIA_VARIABLE(scalar1,double,"scalar1")
        ABORIA_VARIABLE(scalar2,double,"scalar2")

    	typedef Particles<std::tuple<scalar1,scalar2>> ParticlesType;
        typedef position_d<3> position;
       	ParticlesType particles;

       	double diameter = 0.1;
        double3 min(-1);
        double3 max(1);
        double3 periodic(true);
        
        double s_init1 = 1.0;
        double s_init2 = 2.0;
        ParticlesType::value_type p;
        get<position>(p) = double3(0,0,0);
        get<scalar1>(p) = s_init1;
        get<scalar2>(p) = s_init2;
       	particles.push_back(p);
        get<position>(p) = double3(diameter*0.9,0,0);
       	particles.push_back(p);
        get<position>(p) = double3(diameter*1.8,0,0);
       	particles.push_back(p);

        const size_t n = 3;

        particles.init_neighbour_search(min,max,diameter,periodic);

        Symbol<scalar1> s1;
        Symbol<scalar2> s2;
        Label<0,ParticlesType> a(particles);
        Label<1,ParticlesType> b(particles);
        auto dx = create_dx(a,b);

        auto A = create_eigen_operator(s1[a] + s2[b]);
        Eigen::VectorXd v(3);
        v << 1, 2, 3;
        Eigen::VectorXd ans(3);
        ans = A*v;
        for (int i=0; i<n; i++) {
            TS_ASSERT_EQUALS(ans[i],(s_init1+s_init2)*v.sum()); 
        }
        v << 0, 2, 1;
        ans = A*v;
        for (int i=0; i<n; i++) {
            TS_ASSERT_EQUALS(ans[i],(s_init1+s_init2)*v.sum()); 
        }


        auto C = create_eigen_operator(s1[a] + s2[b], norm(dx) < diameter);
        v << 1, 2, 3;
        ans = C*v;
        
        for (int i=0; i<n; i++) {
            double sum = 0;
            for (int j=0; j<n; j++) {
                if ((get<id>(particles[i]) == 0) && (get<id>(particles[j]) == 2)) {
                    sum += 0;
                } else if ((get<id>(particles[i]) == 2) && (get<id>(particles[j]) == 0)) {
                    sum += 0;
                } else {
                    sum += (s_init1+s_init2)*v[j];
                }
            }
            TS_ASSERT_EQUALS(ans[i],sum); 
        }


#endif // HAVE_EIGEN
    }




};

#endif /* OPERATORSTEST_H_ */
