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


#ifndef SPEED_TEST_H_
#define SPEED_TEST_H_


#include <cxxtest/TestSuite.h>
#ifdef HAVE_EIGEN
#include <Eigen/Core>
#endif
#include "Aboria.h"
#include <chrono>
typedef std::chrono::system_clock Clock;
#include <fstream>      // std::ofstream
#ifdef HAVE_GPERFTOOLS
#include <gperftools/profiler.h>
#endif

using namespace Aboria;


class SpeedTest : public CxxTest::TestSuite {
public:
    // BLAS Level 1
    // Dense Vector Addition (c = a + b)
    double vector_addition_aboria_level1(const size_t N) {
        std::cout << "vector_addition_aboria_level1: N = "<<N<<std::endl;
        ABORIA_VARIABLE(a,double,"a")
        ABORIA_VARIABLE(b,double,"b")
        ABORIA_VARIABLE(c,double,"c")
    	typedef Particles<std::tuple<a,b,c>,3> nodes_type;
       	nodes_type nodes;
        nodes_type::value_type node;
        for (int i=0; i<N; i++) {
            get<a>(node) = i;
            get<b>(node) = i*2;
            nodes.push_back(node);
        }
        auto t0 = Clock::now();
        for (int i=0; i<N; i++) {
            get<c>(nodes)[i] = get<a>(nodes)[i] + get<b>(nodes)[i];
            //get<c>(nodes[i]) = get<a>(nodes[i]) + get<b>(nodes[i]);
        }
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        return dt.count();
    }
    
    double vector_addition_aboria_level2(const size_t N) {
        std::cout << "vector_addition_aboria_level2: N = "<<N<<std::endl;
        ABORIA_VARIABLE(a,double,"a")
        ABORIA_VARIABLE(b,double,"b")
        ABORIA_VARIABLE(c,double,"c")
    	typedef Particles<std::tuple<a,b,c>,3> nodes_type;
       	nodes_type nodes;
        nodes_type::value_type node;
        for (int i=0; i<N; i++) {
            get<a>(node) = i;
            get<b>(node) = i*2;
            nodes.push_back(node);
        }
        Symbol<a> A;
        Symbol<b> B;
        Symbol<c> C;
        Label<0,nodes_type> i(nodes);
        auto t0 = Clock::now();
        C[i] = A[i] + B[i];
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        return dt.count();
    }

    double vector_addition_eigen(const size_t N) {
        std::cout << "vector_addition_eigen: N = "<<N<<std::endl;
#ifdef HAVE_EIGEN
        typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type; 
        vector_type A(N),B(N),C(N);
        for (int i=0; i<N; i++) {
            A[i] = i;
            B[i] = 2*i;
        }
        auto t0 = Clock::now();
        C = A + B;
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        return dt.count();
#endif
    }

    // Daxpy (b += a*0.001) 
    double daxpy_aboria_level1(const size_t N) {
        std::cout << "daxpy_aboria_level1: N = "<<N<<std::endl;
        ABORIA_VARIABLE(a,double,"a")
        ABORIA_VARIABLE(b,double,"b")
    	typedef Particles<std::tuple<a,b>,3> nodes_type;
       	nodes_type nodes;
        nodes_type::value_type node;
        for (int i=0; i<N; i++) {
            get<a>(node) = i;
            get<b>(node) = i*2;
            nodes.push_back(node);
        }
        auto t0 = Clock::now();
        for (int i=0; i<N; i++) {
            get<b>(nodes)[i] += get<a>(nodes)[i]*0.001;
        }
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        return dt.count();
    }
    
    double daxpy_aboria_level2(const size_t N) {
        std::cout << "daxpy_aboria_level2: N = "<<N<<std::endl;
        ABORIA_VARIABLE(a,double,"a")
        ABORIA_VARIABLE(b,double,"b")
    	typedef Particles<std::tuple<a,b>,3> nodes_type;
       	nodes_type nodes;
        nodes_type::value_type node;
        for (int i=0; i<N; i++) {
            get<a>(node) = i;
            get<b>(node) = i*2;
            nodes.push_back(node);
        }
        Symbol<a> A;
        Symbol<b> B;
        Label<0,nodes_type> i(nodes);
        auto t0 = Clock::now();
        B[i] += A[i]*0.001;
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        return dt.count();
    }

    double daxpy_eigen(const size_t N) {
        std::cout << "daxpy_eigen: N = "<<N<<std::endl;
#ifdef HAVE_EIGEN
        typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type; 
        vector_type A(N),B(N);
        for (int i=0; i<N; i++) {
            A[i] = i;
            B[i] = 2*i;
        }
        auto t0 = Clock::now();
        B += A*0.001;
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        return dt.count();
#endif
    }


    double finite_difference_eigen(const size_t N) {
        std::cout << "finite_difference_eigen: N = "<<N<<std::endl;
#ifdef HAVE_EIGEN
        const size_t N3 = N*N*N;
        Eigen::SparseMatrix<double> A(N3,N3);
        typedef Eigen::Triplet<double> triplet_type;
        std::vector<triplet_type> tripletList;
        tripletList.reserve(6*N3); 
        Eigen::VectorXd s(N3);
        const double h = 1.0/N; 
        const double invh2 = 1.0/(h*h);
        const double delta_t = 0.1;

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

        auto t0 = Clock::now();
        s += A*s;
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        return dt.count();
    
#endif
    }

    double finite_difference_aboria(const size_t N) {
        std::cout << "finite_difference_aboria: N = "<<N<<std::endl;
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
        const double delta_t = 0.1;
        
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

        Symbol<scalar> s;
        Symbol<id> id_;
        Label<0,nodes_type> a(nodes);
        Label<1,nodes_type> b(nodes);
        auto dx = create_dx(a,b);
        Accumulate<std::plus<double> > sum;
        s.resize_buffer(nodes);

        auto t0 = Clock::now();
        s[a] += invh2*sum(b, norm(dx)<htol, if_else(id_[a]==id_[b],-6,1)*s[b]);
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        return dt.count();
    }

    double finite_difference_aboria_eigen(const size_t N) {
#ifdef HAVE_EIGEN
        std::cout << "finite_difference_aboria_eigen: N = "<<N<<std::endl;
        ABORIA_VARIABLE(scalar,double,"scalar")

    	typedef Particles<std::tuple<scalar>,3> nodes_type;
        typedef position_d<3> position;
       	nodes_type nodes;

        const double h = 1.0/N; 
        double3 min(-h/2);
        double3 max(1+h/2);
        double3 periodic(false);
        Eigen::VectorXd s(N*N*N);
        const double htol = 1.01*h;
        const double invh2 = 1.0/(h*h);
        const double delta_t = 0.1;
        
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

        Symbol<id> id_;
        Label<0,nodes_type> a(nodes);
        Label<1,nodes_type> b(nodes);
        auto dx = create_dx(a,b);
        Accumulate<std::plus<double> > sum;

        auto A = create_eigen_operator(a,b, 
                if_else(id_[a]==id_[b],6,-1), norm(dx)<htol );

        auto t0 = Clock::now();
        s += A*s;
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        return dt.count();
#endif
    }

    double multiquadric_aboria(const size_t N) {
        std::cout << "multiquadric_aboria: N = "<<N<<std::endl;
        ABORIA_VARIABLE(scalar,double,"scalar")
        ABORIA_VARIABLE(kernel_constant,double,"kernel constant")

    	typedef Particles<std::tuple<scalar,kernel_constant>,3> nodes_type;
        typedef position_d<3> position;
       	nodes_type nodes;

        const double h = 1.0/N; 
        double3 min(-h/2);
        double3 max(1+h/2);
        double3 periodic(false);
        const double htol = 1.01*h;
        const double invh2 = 1.0/(h*h);
        const double delta_t = 0.1;
        
        nodes_type::value_type p;
        for (size_t i=0; i<N; ++i) {
            for (size_t j=0; j<N; ++j) {
                for (size_t k=0; k<N; ++k) {
                    get<position>(p) = double3(i*h,j*h,k*h);
                    get<scalar>(p) = 1.5;
                    get<kernel_constant>(p) = 0.1;
       	            nodes.push_back(p);
                }
            }
        }

        nodes.init_neighbour_search(min,max,h,periodic);

        Symbol<scalar> s;
        Symbol<kernel_constant> c;
        Label<0,nodes_type> a(nodes);
        Label<1,nodes_type> b(nodes);
        auto dx = create_dx(a,b);
        Accumulate<std::plus<double> > sum;
        s.resize_buffer(nodes);
        auto t0 = Clock::now();
#ifdef HAVE_GPERFTOOLS
        ProfilerStart("multiquadric_aboria");
#endif
        s[a] += sum(b,true,sqrt(dot(dx,dx)+pow(c[b],2))*s[b]);
#ifdef HAVE_GPERFTOOLS
        ProfilerStop();
#endif
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        return dt.count();
    }

    double multiquadric_aboria_eigen(const size_t N) {
        std::cout << "multiquadric_aboria_eigen: N = "<<N<<std::endl;
        ABORIA_VARIABLE(scalar,double,"scalar")
        ABORIA_VARIABLE(kernel_constant,double,"kernel constant")

    	typedef Particles<std::tuple<scalar,kernel_constant>,3> nodes_type;
        typedef position_d<3> position;
       	nodes_type nodes;

        const double h = 1.0/N; 
        double3 min(-h/2);
        double3 max(1+h/2);
        double3 periodic(false);
        const double htol = 1.01*h;
        const double invh2 = 1.0/(h*h);
        const double delta_t = 0.1;
        Eigen::VectorXd s(N*N*N);
        
        nodes_type::value_type p;
        for (size_t i=0; i<N; ++i) {
            for (size_t j=0; j<N; ++j) {
                for (size_t k=0; k<N; ++k) {
                    const size_t index = i*N*N + j*N + k;
                    get<position>(p) = double3(i*h,j*h,k*h);
                    s[index] = 1.5;
                    get<kernel_constant>(p) = 0.1;
       	            nodes.push_back(p);
                }
            }
        }

        nodes.init_neighbour_search(min,max,h,periodic);

        Symbol<kernel_constant> c;
        Label<0,nodes_type> a(nodes);
        Label<1,nodes_type> b(nodes);
        auto dx = create_dx(a,b);
        Accumulate<std::plus<double> > sum;

        auto A = create_eigen_operator(a,b, 
                sqrt(dot(dx,dx)+pow(c[b],2))
                );

        auto t0 = Clock::now();
        s += A*s;
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        return dt.count();
    }

    double multiquadric_eigen(const size_t N) {
#ifdef HAVE_EIGEN
        std::cout << "multiquadric_eigen: N = "<<N<<std::endl;

        const double h = 1.0/N; 
        double3 min(-h/2);
        double3 max(1+h/2);
        double3 periodic(false);
        const double htol = 1.01*h;
        const double invh2 = 1.0/(h*h);
        const double delta_t = 0.1;
        const size_t N3 = N*N*N;
        Eigen::VectorXd s(N3);
        Eigen::MatrixXd A(N3,N3);
        
        for (size_t i=0; i<N; ++i) {
            for (size_t j=0; j<N; ++j) {
                for (size_t k=0; k<N; ++k) {
                    const size_t index = i*N*N + j*N + k;
                    const double3 r = double3(i*h,j*h,k*h);
                    const double c = 1.5;
                    for (size_t ii=0; ii<N; ++ii) {
                        for (size_t jj=0; jj<N; ++jj) {
                            for (size_t kk=0; kk<N; ++kk) {
                                const double3 r2 = double3(ii*h,jj*h,kk*h);
                                const double3 dx = r2-r;
                                const size_t index2 = ii*N*N + jj*N + kk;
                                A(index,index2) = std::sqrt(dx.dot(dx)+std::pow(c,2));
                            }
                        }
                    }
                    s[index] = 0.1;
                }
            }
        }

        auto t0 = Clock::now();
        s += A*s;
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        return dt.count();
#endif
    }


    void test_multiquadric() {
        std::ofstream file;
        file.open("multiquadric.csv");
        file <<"#"<< std::setw(14) << "N" 
             << std::setw(15) << "aboria" 
             << std::setw(15) << "aboria_eigen" 
             << std::setw(15) << "eigen" << std::endl;
        for (int i = 2; i < 25; ++i) {
            const size_t N = i;
            file << std::setw(15) << std::pow(N,6)
                 << std::setw(15) << std::pow(N,6)/multiquadric_aboria(N)
                 << std::setw(15) << std::pow(N,6)/multiquadric_aboria_eigen(N)
                 << std::setw(15) << std::pow(N,6)/multiquadric_eigen(N)
                 << std::endl;
        }
        file.close();
    }



    void test_finite_difference() {
        std::ofstream file;
        file.open("finite_difference.csv");
        file <<"#"<< std::setw(14) << "N" 
             << std::setw(15) << "aboria" 
             << std::setw(15) << "aboria_eigen" 
             << std::setw(15) << "eigen" << std::endl;
        for (int i = 2; i < 25; ++i) {
            const size_t N = i;
            file << std::setw(15) << std::pow(N,3)
                 << std::setw(15) << std::pow(N,3)/finite_difference_aboria(N)
                 << std::setw(15) << std::pow(N,3)/finite_difference_aboria_eigen(N)
                 << std::setw(15) << std::pow(N,3)/finite_difference_eigen(N)
                 << std::endl;
        }
        file.close();
    }

    void test_vector_addition() {
        std::ofstream file;
        file.open("vector_addition.csv");
        file <<"#"<< std::setw(14) << "N" 
             << std::setw(15) << "aboria_level1" 
             << std::setw(15) << "aboria_level2" 
             << std::setw(15) << "eigen" << std::endl;
        for (int i = 10; i < 8e6; i*=1.2) {
            const size_t N = i;
            file << std::setw(15) << N
                 << std::setw(15) << N/vector_addition_aboria_level1(N)
                 << std::setw(15) << N/vector_addition_aboria_level2(N)
                 << std::setw(15) << N/vector_addition_eigen(N)
                 << std::endl;
        }
        file.close();
    }

    void test_daxpy() {
        std::ofstream file;
        file.open("daxpy.csv");
        file <<"#"<< std::setw(14) << "N" 
             << std::setw(15) << "aboria_level1" 
             << std::setw(15) << "aboria_level2" 
             << std::setw(15) << "eigen" << std::endl;
        for (int i = 10; i < 8e6; i*=1.2) {
            const size_t N = i;
            file << std::setw(15) << N
                 << std::setw(15) << N/daxpy_aboria_level1(N)
                 << std::setw(15) << N/daxpy_aboria_level2(N)
                 << std::setw(15) << N/daxpy_eigen(N)
                 << std::endl;
        }
        file.close();

    }




};

#endif /* OPERATORSTEST_H_ */
