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


#ifndef UTILS_H_
#define UTILS_H_

#include <cxxtest/TestSuite.h>

#include "Aboria.h"
#include "detail/LowRank.h"

using namespace Aboria;


class UtilsTest : public CxxTest::TestSuite {
public:
    void test_bucket_indicies(void) {
        typedef Vector<unsigned int,3> vect;
        detail::bucket_index<3> bi(vect(4,7,2));
        unsigned int index = bi.collapse_index_vector(vect(1,2,1));
        // index = 1 + 2*2 + 1*2*7 = 19 
    	TS_ASSERT_EQUALS(index,19);
        vect vindex = bi.reassemble_index_vector(index);
    	TS_ASSERT_EQUALS(vindex[0],1);
    	TS_ASSERT_EQUALS(vindex[1],2);
    	TS_ASSERT_EQUALS(vindex[2],1);

        typedef Vector<unsigned int,4> vect2;
        detail::bucket_index<4> bi2(vect2(4,7,1,6));
        index = bi2.collapse_index_vector(vect2(1,2,0,4));
        // index = 4 + 0*6 + 2*6*1 + 1*6*1*7 = 58
    	TS_ASSERT_EQUALS(index,58);
        vect2 vindex2 = bi2.reassemble_index_vector(index);
    	TS_ASSERT_EQUALS(vindex2[0],1);
    	TS_ASSERT_EQUALS(vindex2[1],2);
    	TS_ASSERT_EQUALS(vindex2[2],0);
    	TS_ASSERT_EQUALS(vindex2[3],4);

    }

    void test_low_rank(void) {
#ifdef HAVE_EIGEN
        const unsigned int D = 2;
        const size_t N = 100;
        typedef Vector<double,D> double_d;
        typedef Vector<int,D> int_d;
        typedef Vector<bool,D> bool_d;
        // randomly generate a bunch of positions over a range 
        const double pos_min = 0;
        const double pos_max = 1;
        std::uniform_real_distribution<double> U(pos_min,pos_max);
        generator_type generator(time(NULL));
        auto gen = std::bind(U, generator);
        typedef Vector<double,D> double_d;
        typedef Vector<int,D> int_d;

        typedef Particles<std::tuple<>,D> ParticlesType;
        typedef ParticlesType::const_reference const_reference;
        typedef typename ParticlesType::position position;
        ParticlesType particles(N);

        for (int i=0; i<N; i++) {
            for (int d=0; d<D; ++d) {
                get<position>(particles)[i][d] = gen();
            }
        }

        const double c = 0.01;
        auto kernel_fun = [&c](const double_d &dx, const_reference pa, const_reference  pb) {
            return std::sqrt(dx.squaredNorm() + c); 
        };

        KernelDense<ParticlesType,ParticlesType,decltype(kernel_fun)> kernel(particles,particles,kernel_fun);

        // fill a matrix with the result
        Eigen::Matrix<double,N,N> fixed_mat;
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> dyn_mat(N,N);

        kernel.assemble(fixed_mat);
        kernel.assemble(dyn_mat);

        double rms_error_scale = 0;
        for (int i=0; i<N; i++) {
            for (int j=0; j<N; j++) {
                double Ztilde = kernel.coeff(i,j);
                rms_error_scale += std::pow(Ztilde,2);
            }
        }

        for (int k = 1; k < N; ++k) {
            Eigen::Matrix<double,N,Eigen::Dynamic> U(N,k);
            Eigen::Matrix<double,Eigen::Dynamic,N> V(k,N);

            detail::adaptive_cross_approximation(fixed_mat,k,0.01,U,V);

            // check accuracy
            double rms_error_fixed = 0;
            for (int i=0; i<N; i++) {
                for (int j=0; j<N; j++) {
                    double Ztilde = 0;
                    for (int kk=0; kk<k; kk++) {
                        Ztilde += U(i,kk)*V(kk,j);
                    }
                    rms_error_fixed += std::pow(Ztilde-fixed_mat(i,j),2);
                }
            }

            std::cout << "fixed-ACA: k = "<<k<<" rms error = "<<std::sqrt(rms_error_fixed/rms_error_scale)<<std::endl;

            detail::adaptive_cross_approximation(dyn_mat,k,0.01,U,V);
            
            // check accuracy
            double rms_error_dyn = 0;
            for (int i=0; i<N; i++) {
                for (int j=0; j<N; j++) {
                    double Ztilde = 0;
                    for (int kk=0; kk<k; kk++) {
                        Ztilde += U(i,kk)*V(kk,j);
                    }
                    rms_error_dyn += std::pow(Ztilde-dyn_mat(i,j),2);
                }
            }

            std::cout << "dyn-ACA: k = "<<k<<" rms error = "<<std::sqrt(rms_error_dyn/rms_error_scale)<<std::endl;
            
            //detail::adaptive_cross_approximation(kernel,k,0.01,U,V);

            // check accuracy
            double rms_error_kernel = 0;
            for (int i=0; i<N; i++) {
                for (int j=0; j<N; j++) {
                    double Ztilde = 0;
                    for (int kk=0; kk<k; kk++) {
                        Ztilde += U(i,kk)*V(kk,j);
                    }
                    rms_error_kernel += std::pow(Ztilde-kernel.coeff(i,j),2);
                }
            }

            std::cout << "kernel-ACA: k = "<<k<<" rms error = "<<std::sqrt(rms_error_kernel/rms_error_scale)<<std::endl;
        }




#endif
    }




};


#endif /* CONSTRUCTORS_H_ */
