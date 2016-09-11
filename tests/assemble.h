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


#ifndef ASSEMBLETEST_H_
#define ASSEMBLETEST_H_


#include <cxxtest/TestSuite.h>

#include "Aboria.h"

using namespace Aboria;


class AssembleTest : public CxxTest::TestSuite {
public:

    void test_Eigen(void) {
#ifdef HAVE_EIGEN
        ABORIA_VARIABLE(scalar,double,"scalar")

    	typedef Particles<std::tuple<scalar>> ParticlesType;
        typedef position_d<3> position;
       	ParticlesType particles;

       	double diameter = 0.1;
        double3 min(-1);
        double3 max(1);
        double3 periodic(true);
        
        double s_init = 1.0;
        ParticlesType::value_type p;
        get<position>(p) = double3(0,0,0);
        get<scalar>(p) = s_init;
       	particles.push_back(p);
        get<position>(p) = double3(diameter*0.9,0,0);
       	particles.push_back(p);
        get<position>(p) = double3(diameter*1.8,0,0);
       	particles.push_back(p);

        const size_t n = 3;

        particles.init_neighbour_search(min,max,diameter,periodic);

        Symbol<scalar> s;
        auto a = create_label<0>(particles);
        auto b = create_label<1>(particles);
        auto dx = create_dx(a,b);

        Eigen::Map<Eigen::Matrix<double,n,1> > s_vect(get<scalar>(particles).data());

        Eigen::Matrix<double,n,n> A;
        assemble(A, s[a] + s[b]);
        s_vect = A*s_vect;
        for (int i=0; i<n; i++) {
            TS_ASSERT_EQUALS(s_vect[i],(2*s_init)*particles.size()); 
            TS_ASSERT_EQUALS(get<scalar>(particles[i]),(2*s_init)*particles.size()); 
            for (int j; j<n; j++) {
                TS_ASSERT_EQUALS(A(i,j),(2*s_init)); 
            }
        }

        s_vect.setOnes();
        for (int i; i<n; i++) {
            TS_ASSERT_EQUALS(get<scalar>(particles[i]),1.0); 
        }

        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> B(n,n);
        assemble(B, s[a] + s[b]);
        s_vect = B*s_vect;
        for (int i=0; i<n; i++) {
            TS_ASSERT_EQUALS(get<scalar>(particles[i]),(2*s_init)*particles.size()); 
            for (int j; j<n; j++) {
                TS_ASSERT_EQUALS(B(i,j),(2*s_init)); 
            }
        }

        Eigen::SparseMatrix<double> C(n,n); 
        assemble(C, s[a] + s[b], norm(dx) < diameter);
        s_vect = C*s_vect;
        const double exp = (2*s_init)*particles.size();

        for (int i=0; i<n; i++) {
            if ((get<id>(particles[i]) == 0) || (get<id>(particles[i]) == 2)) {
                TS_ASSERT_EQUALS(get<scalar>(particles[i]),2*2*exp*exp); 
            }
            if (get<id>(particles[i]) == 1) {
                TS_ASSERT_EQUALS(get<scalar>(particles[i]),3*2*exp*exp); 
            }
            for (int j=0; j<n; j++) {
                if ((get<id>(particles[i]) == 1) || (get<id>(particles[j]) == 1)) {
                    TS_ASSERT_EQUALS(C.coeff(i,j),2*exp); 
                }
                else if (get<id>(particles[i]) == 0) {
                    if ((get<id>(particles[j]) == 0) || (get<id>(particles[j]) == 1)) {
                        TS_ASSERT_EQUALS(C.coeff(i,j),2*exp); 
                    } else {
                        TS_ASSERT_EQUALS(C.coeff(i,j),0); 
                    }
                }
                else if (get<id>(particles[i]) == 2) {
                    if ((get<id>(particles[j]) == 2) || (get<id>(particles[j]) == 1)) {
                        TS_ASSERT_EQUALS(C.coeff(i,j),2*exp); 
                    } else {
                        TS_ASSERT_EQUALS(C.coeff(i,j),0); 
                    }
                }
            }
        }

        s_vect.setOnes();
        for (int i=0; i<n; i++) {
            TS_ASSERT_EQUALS(get<scalar>(particles[i]),1.0); 
        }

        Eigen::SparseMatrix<double> D(n,n); 
        assemble(D, s[a] + s[b], norm(dx) < diameter);

        s_vect = D*s_vect;

        for (int i=0; i<n; i++) {
            if ((get<id>(particles[i]) == 0) || (get<id>(particles[i]) == 2)) {
                TS_ASSERT_EQUALS(get<scalar>(particles[i]),2*2); 
            }
            if (get<id>(particles[i]) == 1) {
                TS_ASSERT_EQUALS(get<scalar>(particles[i]),3*2); 
            }
        }
#endif // HAVE_EIGEN
    }




};

#endif /* OPERATORSTEST_H_ */
