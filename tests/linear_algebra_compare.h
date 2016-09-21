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


#ifndef LINEAR_ALGEBRA_COMPARE_TEST_H_
#define LINEAR_ALGEBRA_COMPARE_TEST_H_


#include <cxxtest/TestSuite.h>

#include "Aboria.h"

using namespace Aboria;


class LinearAlgebraCompareTest : public CxxTest::TestSuite {
public:

    void finite_difference_eigen(const size_t N, const size_t timesteps) {
#ifdef HAVE_EIGEN
        const size_t N3 = N*N*N;
        Eigen::SparseMatrix<double> A(N3,N3);
        typedef Eigen::Triplet<double> triplet_type;
        std::vector<triplet_type> tripletList;
        tripletList.reserve(6*N3); 
        Eigen::VectorXd s(N3);
        const double h = 1.0/N; 
        const double invh2 = 1.0/(h*h);
        const double dt = 0.1;

        printf("eigen: start setup.\n");
        for (size_t i=0; i<N; ++i) {
            for (size_t j=0; j<N; ++j) {
                for (size_t k=0; k<N; ++k) {
                    const size_t index = i*N*N+j*N+k;
                    assert(index < N3);
                    tripletList.push_back(triplet_type(index,index,-6*invh2));
                    if (index>=1) tripletList.push_back(triplet_type(index,index-1,invh2));
                    if (index+1<N3) tripletList.push_back(triplet_type(index,index+1,invh2));
                    if (index>=N) tripletList.push_back(triplet_type(index,index-N,invh2));
                    if (index+N<N3) tripletList.push_back(triplet_type(index,index+N,invh2));
                    if (index>=N*N) tripletList.push_back(triplet_type(index,index-N*N,invh2));
                    if (index+N*N<N3) tripletList.push_back(triplet_type(index,index+N*N,invh2));
                    s(index) = std::exp((double3(i*h,j*h,k*h)-double3(0.5,0.5,0.5)).squaredNorm());
                }
            }
        }
        A.setFromTriplets(tripletList.begin(),tripletList.end());

        printf("eigen: end setup.\n");


        for (size_t i=0; i<timesteps; ++i) {
            std::cout << "." <<std::flush;
            s += dt*A*s;
        }

        printf("eigen: end run.\n");
#endif
    }

    void finite_difference_aboria(const size_t N, const size_t timesteps) {
        ABORIA_VARIABLE(scalar,double,"scalar")

    	typedef Particles<std::tuple<scalar>,3> nodes_type;
        typedef position_d<3> position;
       	nodes_type nodes;

        const double h = 1.0/N; 
        double3 min(-h/2);
        double3 max(1+h/2);
        double3 periodic(false);
        const double htol = 1.01*h;
        const double invh2 = 1.0/(h*h);
        const double dt = 0.1;
        
        printf("aboria: start setup.\n");
        nodes_type::value_type p;
        for (size_t i=0; i<N; ++i) {
            for (size_t j=0; j<N; ++j) {
                for (size_t k=0; k<N; ++k) {
                    get<position>(p) = double3(i*h,j*h,k*h);
                    get<scalar>(p) = std::exp((get<position>(p)-double3(0.5,0.5,0.5)).squaredNorm());
       	            nodes.push_back(p);
                }
            }
        }

        nodes.init_neighbour_search(min,max,h,periodic);

        printf("aboria: end setup.\n");
        Symbol<scalar> s;
        Symbol<id> id_;
        Label<0,nodes_type> a(nodes);
        Label<1,nodes_type> b(nodes);
        auto dx = create_dx(a,b);
        Accumulate<std::plus<double> > sum;

        for (size_t i=0; i<timesteps; ++i) {
            std::cout << "." <<std::flush;
            //s[a] += dt*invh2*sum(b, norm(dx)<htol, if_else(id_[a]==id_[b],-6,1)*s[b]);
            s[a] += dt*invh2*sum(b, norm(dx)<htol, 6*s[b]);
        }
        printf("aboria: end run.\n");

        /*
        auto A = create_eigen_operator(a,b, 
                if_else(id_[a]==id_[b],6,-1), norm(dx)<htol );

        Eigen::VectorXd u(N);
        v << 1, 2, 3;
        Eigen::VectorXd ans(3);
        ans = A*v;

        */


    }

    void test_eigen() {
        finite_difference_eigen(20,100); 
    }

    void test_aboria() {
        finite_difference_aboria(20,100); 
    }


};

#endif /* OPERATORSTEST_H_ */
