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
#include <thread>
#ifdef HAVE_GPERFTOOLS
#include <gperftools/profiler.h>
#endif

#ifdef HAVE_GROMACS
#include <gromacs/selection/nbsearch.h>
#include <gromacs/pbcutil/pbc.h>

#endif

#include <boost/math/constants/constants.hpp>
const double PI = boost::math::constants::pi<double>();

using namespace Aboria;

class SpeedTest : public CxxTest::TestSuite {
public:
    // BLAS Level 1
    // Dense Vector Addition (c = a + b)
    
    
    
    double vector_addition_aboria_level2(const size_t N, const size_t repeats) {
        std::cout << "vector_addition_aboria_level2: N = "<<N<<std::endl;
        //[vector_addition
        /*`
        Here we aim to compute a simple vector addition operation

        $$
        a_i = b_i + c_i \text{ for } i = 0...N.
        $$

        A particle set containing the variables $a$, $b$ and $c$ can be defined in 
        __aboria__ like so

        */
        ABORIA_VARIABLE(a_var,double,"a")
        ABORIA_VARIABLE(b_var,double,"b")
        ABORIA_VARIABLE(c_var,double,"c")
    	typedef Particles<std::tuple<a_var,b_var,c_var>,3> nodes_type;
       	nodes_type nodes(N);
        //<- 
        for (size_t i=0; i<N; i++) {
            get<a_var>(nodes)[i] = i;
            get<b_var>(nodes)[i] = i*2;
        }
        //->
        /*`        
        The vector addition operation can then be calculated using the Level 3 layer 
        like so
        */
        Symbol<a_var> a;
        Symbol<b_var> b;
        Symbol<c_var> c;
        Label<0,nodes_type> i(nodes);
        a[i] = b[i] + c[i];
        //<- 
        auto t0 = Clock::now();
#ifdef HAVE_GPERFTOOLS
        ProfilerStart("vector_addition_aboria_level2");
#endif
        for (int r=0; r<repeats; ++r) {
            a[i] = b[i] + c[i];
        }
#ifdef HAVE_GPERFTOOLS
        ProfilerStop();
#endif
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        return dt.count()/repeats;
    }

    double vector_addition_aboria_level1(const size_t N, const size_t repeats) {
        std::cout << "vector_addition_aboria_level1: N = "<<N<<std::endl;
        ABORIA_VARIABLE(a_var,double,"a")
        ABORIA_VARIABLE(b_var,double,"b")
        ABORIA_VARIABLE(c_var,double,"c")
    	typedef Particles<std::tuple<a_var,b_var,c_var>,3> nodes_type;
       	nodes_type nodes(N);
        for (size_t i=0; i<N; i++) {
            get<a_var>(nodes)[i] = i;
            get<b_var>(nodes)[i] = i*2;
        }
        //->
        /*`
        We compare this with Level 1 Aboria using the `get` functions and looping
        through the container
        */
        for (size_t i=0; i<N; i++) {
            get<a_var>(nodes)[i] = get<b_var>(nodes)[i] + get<c_var>(nodes)[i];
        }
        //<-
        auto t0 = Clock::now();
        for (int r=0; r<repeats; ++r) {
            for (size_t i=0; i<N; i++) {
                get<a_var>(nodes)[i] = get<b_var>(nodes)[i] + get<c_var>(nodes)[i];
            }
        }
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        return dt.count()/repeats;
    }

    double vector_addition_stdvector(const size_t N, const size_t repeats) {
        std::cout << "vector_addition_stdvector: N = "<<N<<std::endl;
        //->
        /*`
        We also compare against a plain `std::vector` implementation like so
        */
        std::vector<double> a(N),b(N),c(N);
        //<-
        for (size_t i=0; i<N; i++) {
            a[i] = i;
            b[i] = i*2;
        }
        //->
        for (size_t i=0; i<N; i++) {
            a[i] = b[i] + c[i];
        }
        //<-
        auto t0 = Clock::now();
        for (int r=0; r<repeats; ++r) {
            for (size_t i=0; i<N; i++) {
                a[i] = b[i] + c[i];
            }
        }
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        return dt.count()/repeats;
    }

    double vector_addition_eigen(const size_t N, const size_t repeats) {
        std::cout << "vector_addition_eigen: N = "<<N<<std::endl;
#ifdef HAVE_EIGEN
        //->
        /*`
        Finally we compare against an Eigen implementation:
        */
        typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type; 
        vector_type a(N),b(N),c(N);
        //<-
        for (size_t i=0; i<N; i++) {
            a[i] = i;
            b[i] = 2*i;
        }
        //->
        a = b + c;
        /*`
        We can measure the time taken by the last line in the code segment above for 
        varying $N$, and compare the four different implementations

        The resultant benchmarks are shown in the Figure below, where it can be seen that 
        the four approaches are very similar in speed, confirming that [Aboria][] can 
        achieve zero-cost abstraction, at least in this simple case. More complicated 
        cases are explored below.

        [$images/benchmarks/vector_addition.svg]
        */
        //]
        auto t0 = Clock::now();
        for (int r=0; r<repeats; ++r) {
            a = b + c;
        }
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        return dt.count()/repeats;
#endif
    }

    // Daxpy (b += a*0.001) 
    double daxpy_aboria_level2(const size_t N, const size_t repeats) {
        std::cout << "daxpy_aboria_level2: N = "<<N<<std::endl;
        ABORIA_VARIABLE(a_var,double,"a")
        ABORIA_VARIABLE(b_var,double,"b")
    	typedef Particles<std::tuple<a_var,b_var>,3> nodes_type;
       	nodes_type nodes(N);
        for (size_t i=0; i<N; i++) {
            get<a_var>(nodes)[i] = i;
            get<b_var>(nodes)[i] = i*2;
        }
        Symbol<a_var> a;
        Symbol<b_var> b;
        Label<0,nodes_type> i(nodes);
        //[daxpy_run
        /*`
        This benchmark is for the BLAS DAXPY operation, given by

        $$
        a_i = a_i + 0.1*b_i  \text{ for } i = 0...N.
        $$

        This is implemented in __aboria__ using 
        */
        a[i] += 0.1*b[i];
        //<-
        auto t0 = Clock::now();
        for (int r=0; r<repeats; ++r) {
            a[i] += 0.1*b[i];
        }
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        return dt.count()/repeats;
    }

    double daxpy_aboria_level1(const size_t N, const size_t repeats) {
        std::cout << "daxpy_aboria_level1: N = "<<N<<std::endl;
        ABORIA_VARIABLE(a_var,double,"a")
        ABORIA_VARIABLE(b_var,double,"b")
    	typedef Particles<std::tuple<a_var,b_var>,3> nodes_type;
       	nodes_type nodes(N);
        for (size_t i=0; i<N; i++) {
            get<a_var>(nodes)[i] = i;
            get<b_var>(nodes)[i] = i*2;
        }
        //->
        /*`
        We compare against a Level 1 implementation like so 
        */
        for (size_t i=0; i<N; i++) {
            get<a_var>(nodes)[i] += 0.1*get<b_var>(nodes)[i];
        }
        //<-
        auto t0 = Clock::now();
        for (int r=0; r<repeats; ++r) {
            for (size_t i=0; i<N; i++) {
                get<a_var>(nodes)[i] += 0.1*get<b_var>(nodes)[i];
            }
        }
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        return dt.count()/repeats;
    }

    double daxpy_stdvector(const size_t N, const size_t repeats) {
        std::cout << "daxpy_aboria_level1: N = "<<N<<std::endl;
        //->
        /*`
        and a `std::vector` implementation like so
        */
        std::vector<double> a(N),b(N);
        //<-
        for (size_t i=0; i<N; i++) {
            a[i] = i;
            b[i] = i*2;
        }
        //->
        for (size_t i=0; i<N; i++) {
            a[i] += 0.1*b[i];
        }
        //<-
        auto t0 = Clock::now();
        for (int r=0; r<repeats; ++r) {
            for (size_t i=0; i<N; i++) {
                a[i] += 0.1*b[i];
            }
        }
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        return dt.count()/repeats;
    }
    

    double daxpy_eigen(const size_t N, const size_t repeats) {
        std::cout << "daxpy_eigen: N = "<<N<<std::endl;
#ifdef HAVE_EIGEN
        //->
        /*`
        and an Eigen implementation
        */
        typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type; 
        vector_type a(N),b(N);
        //<-
        for (size_t i=0; i<N; i++) {
            a[i] = i;
            b[i] = 2*i;
        }
        //->
        a += 0.1*b;
        /*`
        The comarison benchmarks for varying $N$ are shown below

        [$images/benchmarks/daxpy.svg]
        */
        //]
        auto t0 = Clock::now();
        for (int r=0; r<repeats; ++r) {
            a += 0.1*b;
        }
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        return dt.count()/repeats;
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
                    s(index) = std::exp((vdouble3(i*h,j*h,k*h)-vdouble3(0.5,0.5,0.5)).squaredNorm());
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
        vdouble3 min(-h/2);
        vdouble3 max(1+h/2);
        vdouble3 periodic(false);
        const double htol = 1.01*h;
        const double invh2 = 1.0/(h*h);
        const double delta_t = 0.1;
        
        nodes_type::value_type p;
        for (size_t i=0; i<N; ++i) {
            for (size_t j=0; j<N; ++j) {
                for (size_t k=0; k<N; ++k) {
                    get<position>(p) = vdouble3(i*h,j*h,k*h);
                    get<scalar>(p) = std::exp((get<position>(p)-vdouble3(0.5,0.5,0.5)).squaredNorm());
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
        a.template resize_buffer<scalar>(nodes.size());

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
        vdouble3 min(-h/2);
        vdouble3 max(1+h/2);
        vdouble3 periodic(false);
        Eigen::VectorXd s(N*N*N);
        const double htol = 1.01*h;
        const double invh2 = 1.0/(h*h);
        const double delta_t = 0.1;
        
        nodes_type::value_type p;
        for (size_t i=0; i<N; ++i) {
            for (size_t j=0; j<N; ++j) {
                for (size_t k=0; k<N; ++k) {
                    get<position>(p) = vdouble3(i*h,j*h,k*h);
                    get<scalar>(p) = std::exp((get<position>(p)-vdouble3(0.5,0.5,0.5)).squaredNorm());
       	            nodes.push_back(p);
                }
            }
        }

        nodes.init_neighbour_search(min,max,htol,periodic);

        auto A = create_sparse_operator(nodes,nodes, 
                htol,
                [](const vdouble3 &dx,
                    nodes_type::const_reference a,
                    nodes_type::const_reference b) {
                    if (get<id>(a)==get<id>(b)) {
                        return 6.0;
                    } else {
                        return -1.0;
                    }
                }
                );

        auto t0 = Clock::now();
        s += A*s;
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        return dt.count();
#endif
    }

    


    template <unsigned int Dim>
    double multiquadric_aboria(const size_t N, const size_t repeats) {
        //typedef Vector<double,Dim> double_d;
        typedef double double_d;
        std::cout << "multiquadric_aboria: N = "<<N<<std::endl;
        ABORIA_VARIABLE(a_var,double_d,"a")
        ABORIA_VARIABLE(b_var,double_d,"b")
    	typedef Particles<std::tuple<a_var,b_var>,2> nodes_type;
        typedef position_d<2> position;
       	nodes_type nodes(N*N);

        const double h = 1.0/N; 
        vdouble2 min(-h/2);
        vdouble2 max(1+h/2);
        vdouble2 periodic(false);
        const double htol = 1.01*h;
        const double invh2 = 1.0/(h*h);
        const double delta_t = 0.1;
        
        for (size_t i=0; i<N; ++i) {
            for (size_t j=0; j<N; ++j) {
                const size_t index = i*N + j;
                get<position>(nodes)[index] = vdouble2(i*h,j*h);
                get<a_var>(nodes)[index] = double_d(1.5);
                get<b_var>(nodes)[index] = double_d(2.5);
            }
        }

        nodes.init_neighbour_search(min,max,h,periodic);

        Symbol<a_var> a;
        Symbol<b_var> b;
        Label<0,nodes_type> i(nodes);
        i.template resize_buffer<a_var>(nodes.size());
        //[multiquadric
        /*`
        Here we move onto a dense, $N^2$ operation, given by the non-linear operator

        $$
        a_i = a_i + \sum_j^N a_j \sqrt{\mathbf{dx}\_{ij} \cdot \mathbf{dx}\_{ij} + b_j^2}    
        \text{ for } i = 0...N.
        $$

        where $\mathbf{dx}\_{ij}$ is the shortest vector from particle $i$ to $i$. This 
        is implemented in Level 3 __aboria__ like so
        */
        Label<1,nodes_type> j(nodes);
        auto dx = create_dx(i,j);
        Accumulate<std::plus<double> > sum;
        a[i] += sum(j,true,a[j]*sqrt(dot(dx,dx)+b[j]*b[j]));
        //<-
        
        auto t0 = Clock::now();
#ifdef HAVE_GPERFTOOLS
        ProfilerStart("multiquadric_aboria");
#endif
        for (int ii=0; ii<repeats; ++ii) {
            a[i] += sum(j,true,a[j]*sqrt(dot(dx,dx)+b[j]*b[j]));
        }
#ifdef HAVE_GPERFTOOLS
        ProfilerStop();
#endif
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        std::cout << "time = "<<dt.count()/repeats<<std::endl;
        return dt.count()/repeats;
    }

    template <unsigned int Dim>
    double multiquadric_vector(const size_t inN, const size_t repeats) {
        const size_t N = inN*inN;
        //typedef Vector<double,Dim> double_d;
        std::cout << "multiquadric_vector: N = "<<N<<std::endl;
        //->
        /*`
        This is compared against a `std::vector` implementation. Note that
        this operator involves aliasing, in that the update variable $a$ 
        appears within the sum, so we need to accumulate the update to a 
        temporary buffer before we assign to $a_i$.

        The implementation is shown below (note the openMP parallel loops
        are turned off for the plot below)
        */

        std::vector<double> x(N), y(N), b(N), a(N), a_buffer(N);

        //<-
        const double h = 1.0/N; 
        vdouble2 min(-h/2);
        vdouble2 max(1+h/2);
        vdouble2 periodic(false);
        const double htol = 1.01*h;
        const double invh2 = 1.0/(h*h);
        const double delta_t = 0.1;
        
        for (size_t i=0; i<inN; ++i) {
            for (size_t j=0; j<inN; ++j) {
                const size_t index = i*inN + j;
                x[index] = i*h;
                y[index] = j*h;
                b[index] = double(0.1);
                a[index] = double(1.5);
            }
        }

        //->
        #pragma omp parallel for
        for (size_t i = 0; i < N; ++i) {
            a_buffer[i] = a[i];
            for (size_t j = 0; j < N; ++j) {
                const double dx_x = x[j]-x[i];
                const double dx_y = y[j]-x[i];
                a_buffer[i] += a[j]*std::sqrt(dx_x*dx_x+dx_y*dx_y+b[j]*b[j]);
            }
        }
        #pragma omp parallel for
        for (size_t i = 0; i < N; ++i) {
            a[i] = a_buffer[i];
        }

        /*`
        The benchmarks are shown below. 

        [$images/benchmarks/multiquadric.svg]
         */
        //]
        auto t0 = Clock::now();
        for (int r=0; r<repeats; ++r) {
            #pragma omp parallel for
            for (size_t i = 0; i < N; ++i) {
                a_buffer[i] = a[i];
                for (size_t j = 0; j < N; ++j) {
                    const double dx_x = x[j]-x[i];
                    const double dx_y = y[j]-x[i];
                    a_buffer[i] += a[j]*std::sqrt(dx_x*dx_x+dx_y*dx_y+b[j]*b[j]);
                }
            }
            #pragma omp parallel for
            for (size_t i = 0; i < N; ++i) {
                a[i] = a_buffer[i];
            }
        }
        auto t1 = Clock::now();

        std::chrono::duration<double> dt = t1 - t0;
        std::cout << "time = "<<dt.count()/repeats<<std::endl;
        return dt.count()/repeats;
    }


    template <unsigned int Dim>
    double multiquadric_aboria_eigen(const size_t N, const size_t repeats) {
#ifdef HAVE_EIGEN
        //typedef Vector<double,Dim> double_d;
        typedef double double_d;
        std::cout << "multiquadric_aboria_eigen: N = "<<N<<std::endl;
        ABORIA_VARIABLE(scalar,double_d,"scalar")
        ABORIA_VARIABLE(kernel_constant,double_d,"kernel constant")

    	typedef Particles<std::tuple<scalar,kernel_constant>,2> nodes_type;
        typedef position_d<2> position;
       	nodes_type nodes(N*N);

        const double h = 1.0/N; 
        vdouble2 min(-h/2);
        vdouble2 max(1+h/2);
        vdouble2 periodic(false);
        const double htol = 1.01*h;
        const double invh2 = 1.0/(h*h);
        const double delta_t = 0.1;
        Eigen::Matrix<double_d,Eigen::Dynamic,1> s(N*N);
        
        for (size_t i=0; i<N; ++i) {
            for (size_t j=0; j<N; ++j) {
                const size_t index = i*N + j;
                get<position>(nodes)[index] = vdouble2(i*h,j*h);
                s[index] = double_d(1.5);
                get<kernel_constant>(nodes)[index] = double_d(0.1);
            }
        }

        nodes.init_neighbour_search(min,max,h,periodic);

        auto A = create_dense_operator(nodes,nodes, 
                [](const vdouble2 &dx,
                    typename nodes_type::const_reference a,
                    typename nodes_type::const_reference b) {
                    return std::sqrt(dx.squaredNorm()+
                                get<kernel_constant>(b)*get<kernel_constant>(b));
                }
                );

        s += A*s;
        auto t0 = Clock::now();
#ifdef HAVE_GPERFTOOLS
#ifdef HAVE_OPENMP
        if (omp_get_max_threads() > 1) {
            ProfilerStart("multiquadric_aboria_eigen_parallel");
        } else {
            ProfilerStart("multiquadric_aboria_eigen_serial");
        }
#else
        ProfilerStart("multiquadric_aboria_eigen_serial");
#endif
#endif
        for (int i=0; i<repeats; ++i) {
            s += A*s;
        }
#ifdef HAVE_GPERFTOOLS
        ProfilerStop();
#endif
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        std::cout << "time = "<<dt.count()/repeats<<std::endl;
        return dt.count()/repeats;
#else
        return 0.0;
#endif
    }

    template <unsigned int Dim>
    double multiquadric_eigen(const size_t N, const size_t repeats) {
        //typedef Vector<double,Dim> double_d;
        typedef double double_d;
#ifdef HAVE_EIGEN
        std::cout << "multiquadric_eigen: N = "<<N<<" with "<<Eigen::nbThreads()<<" threads"<<std::endl;

        const double h = 1.0/N; 
        vdouble2 min(-h/2);
        vdouble2 max(1+h/2);
        vdouble2 periodic(false);
        const double htol = 1.01*h;
        const double invh2 = 1.0/(h*h);
        const double delta_t = 0.1;
        const size_t N2 = N*N;
        Eigen::Matrix<double_d,Eigen::Dynamic,1> s(N2);
        Eigen::Matrix<double_d,Eigen::Dynamic,1> c(N2);
        Eigen::Matrix<double_d,Eigen::Dynamic,Eigen::Dynamic> A(N2,N2);

        for (size_t i=0; i<N; ++i) {
            for (size_t j=0; j<N; ++j) {
                const size_t index = i*N + j;
                s[index] = double_d(0.1);
                c[index] = double_d(1.0);
            }
        }
         
        for (size_t i=0; i<N; ++i) {
            for (size_t j=0; j<N; ++j) {
                const size_t index = i*N + j;
                const vdouble2 r = vdouble2(i*h,j*h);
                for (size_t ii=0; ii<N; ++ii) {
                    for (size_t jj=0; jj<N; ++jj) {
                        const vdouble2 r2 = vdouble2(ii*h,jj*h);
                        const vdouble2 dx = r2-r;
                        const size_t index2 = ii*N + jj;
                        A(index,index2) = std::sqrt(dx.dot(dx)+c[index2]*c[index2]);
                    }
                }
            }
        }

        auto t0 = Clock::now();
        for (int i=0; i<repeats; ++i) {
            s += A*s;
        }
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        std::cout << "time = "<<dt.count()/repeats<<std::endl;
        return dt.count()/repeats;
#else
        return 0.0;
#endif
    }

    template <template <typename> class SearchMethod>
    double linear_spring_aboria(const size_t N, const double radius, const size_t repeats) {
        std::cout << "linear_spring_aboria: N = "<<N<<std::endl;

        const double r = radius;

        ABORIA_VARIABLE(a_var,vdouble3,"a")
    	typedef Particles<std::tuple<a_var>,3,std::vector,SearchMethod> nodes_type;
        typedef position_d<3> position;
       	nodes_type nodes(N*N*N);

        const double h = 1.0/N; 
        vdouble3 min(-h/2);
        vdouble3 max(1+h/2);
        vdouble3 periodic(false);
        const double htol = 1.01*h;
        const double invh2 = 1.0/(h*h);
        const double delta_t = 0.1;
        
        for (size_t i=0; i<N; ++i) {
            for (size_t j=0; j<N; ++j) {
                for (size_t k=0; k<N; ++k) {
                    const size_t index = i*N*N + j*N + k;
                    get<position>(nodes)[index] = vdouble3(i*h,j*h,k*h);
                }
            }
        }

        nodes.init_neighbour_search(min,max,r,periodic);

        Symbol<position> p;
        Symbol<a_var> a;
        Label<0,nodes_type> i(nodes);
        Label<1,nodes_type> j(nodes);
        auto dx = create_dx(i,j);
        Accumulate<std::plus<vdouble3> > sum;
        Accumulate<std::plus<int> > sumi;

        //[linear_spring
        /*`
        Finally we implement a non-linear operator involving a neighbour search, common 
        in particle-based methods. This is given by

        $$
        a_i = \sum_j^N \begin{cases}
                    \frac{r-|\mathbf{dx}\_{ij}|}{|\mathbf{dx}\_{ij}|}\mathbf{dx}\_{ij} , & 
                    \text{for } 
                      |\mathbf{dx}\_{ij}|<r \\\\
                    0 & \text{otherwise},
                    \end{cases}   \text{ for } i = 0...N.
        $$

        where $r$ is a given constant.

        */
        a[i] = sum(j,norm(dx)<r,(r-norm(dx))/norm(dx)*dx);
        /*`
        The benchmarks are shown below. 

        [$images/benchmarks/linear_spring.svg]
         */
        //]
        
        const int count = eval(sumi(i,true,sumi(j,norm(dx)<r,1)));
        std::cout << "found "<<count<<" pairs"<<std::endl;
        auto t0 = Clock::now();
#ifdef HAVE_GPERFTOOLS
        ProfilerStart("linear_spring_aboria");
#endif
        for (int ii=0; ii<repeats; ++ii) {
            a[i] = sum(j,norm(dx)<r,(r-norm(dx))/norm(dx)*dx);
        }
#ifdef HAVE_GPERFTOOLS
        ProfilerStop();
#endif
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        std::cout << "time = "<<dt.count()/repeats<<std::endl;
        return dt.count()/repeats;
    }

    double linear_spring_gromacs(const size_t N, const double radius, const size_t repeats) {
#ifdef HAVE_GROMACS
        std::cout << "linear_spring_gromacs: N = "<<N<<std::endl;

        const double r = radius;

    	typedef Particles<std::tuple<>,2> nodes_type;
        typedef position_d<2> position;
       	nodes_type nodes(N*N*N);

        const double h = 1.0/N; 
        vdouble2 min(-h/2);
        vdouble2 max(1+h/2);
        vdouble2 periodic(false);
        const double htol = 1.01*h;
        const double invh2 = 1.0/(h*h);
        const double delta_t = 0.1;

        std::vector<gmx::RVec> positions(N*N*N);
        std::vector<gmx::RVec> results(N*N*N);
        
        for (size_t i=0; i<N; ++i) {
            for (size_t j=0; j<N; ++j) {
                for (size_t k=0; k<N; ++k) {
                    const size_t index = i*N*N + j*N + k;
                    positions[index][0] = i*h;
                    positions[index][1] = j*h;
                    positions[index][2] = k*h;
                }
            }
        }

        t_pbc pbc;
        matrix box;
        // not sure on format of box, but since using epbcNONE I
        // think I can ignore it
        set_pbc(&pbc, epbcNONE, box);
        gmx::AnalysisNeighborhoodPositions analysis_positions(positions);
        gmx::AnalysisNeighborhood neighborhood;
        neighborhood.setCutoff(r);
        //neighborhood.setMode(gmx::eSearchMode_Grid)
        gmx::AnalysisNeighborhoodSearch nbsearch = neighborhood.initSearch(&pbc,analysis_positions);

        for (int a=0; a<N*N*N; ++a) {
            results[a][0] = 0;
            results[a][1] = 0;
            results[a][2] = 0;
        }
        gmx::AnalysisNeighborhoodPairSearch pairSearch = nbsearch.startPairSearch(analysis_positions);
        gmx::AnalysisNeighborhoodPair pair;
        int count = 0;
        while (pairSearch.findNextPair(&pair)) {
            ++count;
            const double d2 = pair.distance2();
            ASSERT(d2<=r*r,"bad search? ");
            const int a = pair.refIndex();
            const int b = pair.testIndex();
            const rvec& dx = pair.dx();
            const double norm_dx = std::sqrt(d2);
            const double scale = (r-norm_dx)/norm_dx;
            results[a][0] += scale*dx[0]; 
            results[a][1] += scale*dx[1]; 
            results[a][2] += scale*dx[2]; 
        }
        
        std::cout << "found "<<count<<" pairs"<<std::endl;

        auto t0 = Clock::now();
        for (int i=0; i<repeats; ++i) {
            for (int a=0; a<N*N*N; ++a) {
                results[a][0] = 0;
                results[a][1] = 0;
                results[a][2] = 0;
            }
            gmx::AnalysisNeighborhoodPairSearch pairSearch = nbsearch.startPairSearch(analysis_positions);
            gmx::AnalysisNeighborhoodPair pair;
            while (pairSearch.findNextPair(&pair)) {
                const double d2 = pair.distance2();
                ASSERT(d2<=r*r,"bad search?");
                const int a = pair.refIndex();
                const int b = pair.testIndex();
                const rvec& dx = pair.dx();
                const double norm_dx = std::sqrt(d2);
                const double scale = (r-norm_dx)/norm_dx;
                results[a][0] += scale*dx[0]; 
                results[a][1] += scale*dx[1]; 
                results[a][2] += scale*dx[2]; 
            }
        }
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        std::cout << "time = "<<dt.count()/repeats<<std::endl;
        return dt.count()/repeats;
#else
        return 0;
#endif
    }

    void test_linear_spring() {
#ifdef HAVE_OPENMP
            omp_set_num_threads(1);
#endif
        std::ofstream file;
        const size_t base_repeats = 5e6;
        for (double radius_div_h = 1.1; radius_div_h < 5; radius_div_h += 1) {
            char buffer[100];
            sprintf(buffer,"linear_spring%4.4f.csv",radius_div_h);
            file.open(buffer);
            file <<"#"<< std::setw(14) << "N" 
                << std::setw(15) << "aboria_serial" 
                << std::setw(15) << "aboria_parallel" 
                << std::setw(15) << "gromacs" << std::endl;


            for (double i = 2; i < 30; i *= 1.05) {
                const size_t N = i;
                const size_t matrix_size = std::pow(N,3)*(4.0/3.0)*PI*std::pow(radius_div_h,3);
                const size_t repeats = base_repeats/matrix_size + 1;
                const double h = 1.0/N; 
                const double radius = radius_div_h*h;
                file << std::setw(15) << std::pow(N,3);
                file << std::setw(15) << std::pow(N,6)/linear_spring_aboria<CellList>(N,radius,repeats);
                file << std::setw(15) << std::pow(N,6)/linear_spring_aboria<CellListOrdered>(N,radius,repeats);
                file << std::setw(15) << std::pow(N,6)/linear_spring_gromacs(N,radius,repeats);
                file << std::endl;
            }
            file.close();
        }
    }

    void test_multiquadric() {
        std::ofstream file;
        const size_t base_repeats = 5e6;
        file.open("multiquadric.csv");
        file <<"#"<< std::setw(14) << "N" 
             << std::setw(15) << "aboria" 
             << std::setw(15) << "aboria_eigen" 
             << std::setw(15) << "eigen" << std::endl;
        for (double i = 2; i < 100; i *= 1.05) {
            const size_t N = i;
            const size_t repeats = base_repeats/std::pow(N,4) + 1;
            file << std::setw(15) << std::pow(N,2);
#ifdef HAVE_OPENMP
            omp_set_num_threads(1);
#ifdef HAVE_EIGEN
            Eigen::setNbThreads(1);
#endif
#endif
            file << std::setw(15) << std::pow(N,4)/multiquadric_aboria<1>(N,repeats);
            file << std::setw(15) << std::pow(N,4)/multiquadric_aboria_eigen<1>(N,repeats);
            file << std::setw(15) << std::pow(N,4)/multiquadric_eigen<1>(N,repeats);
            file << std::setw(15) << std::pow(N,4)/multiquadric_vector<1>(N,repeats);
#ifdef HAVE_OPENMP
            omp_set_num_threads(4);
#ifdef HAVE_EIGEN
            Eigen::setNbThreads(4);
#endif
#endif
            file << std::setw(15) << std::pow(N,4)/multiquadric_aboria<1>(N,repeats);
            file << std::setw(15) << std::pow(N,4)/multiquadric_aboria_eigen<1>(N,repeats);
            file << std::setw(15) << std::pow(N,4)/multiquadric_eigen<1>(N,repeats);
            file << std::setw(15) << std::pow(N,4)/multiquadric_vector<1>(N,repeats);
            file << std::endl;
        }
        file.close();
    }

    void test_multiquadric_scaling() {
#ifdef HAVE_OPENMP
        std::ofstream file;
        const size_t repeats = 10;
        const size_t N = 100;
        file.open("multiquadric_scaling.csv");
        
        file <<"#"<< std::setw(14) << "N" 
             << std::setw(15) << "aboria" 
             << std::setw(15) << "aboria_eigen" 
             << std::setw(15) << "eigen" << std::endl;
        int num_cores  = std::thread::hardware_concurrency();
        for (int i = 1; i <= num_cores; ++i) {
            omp_set_num_threads(i);
#ifdef HAVE_EIGEN
            Eigen::setNbThreads(i);
#endif
            file << std::setw(15) << i;
            file << std::setw(15) << std::pow(N,4)/multiquadric_aboria<1>(N,repeats);
            file << std::setw(15) << std::pow(N,4)/multiquadric_aboria_eigen<1>(N,repeats);
            file << std::setw(15) << std::pow(N,4)/multiquadric_eigen<1>(N,repeats);
            file << std::setw(15) << std::pow(N,4)/multiquadric_vector<1>(N,repeats);
            file << std::endl;
        }
        file.close();
#endif //HAVE_OPENMP
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
        const size_t base_repeats = 1e7;
        file.open("vector_addition.csv");
        file <<"#"<< std::setw(14) << "N" 
             << std::setw(15) << "aboria_level1" 
             << std::setw(15) << "aboria_level2" 
             << std::setw(15) << "eigen"
             << std::setw(15) << "stdvector" << std::endl;

#ifdef HAVE_OPENMP
            omp_set_num_threads(1);
#endif
        for (int i = 10; i < 8e6; i*=1.2) {
            const size_t N = i;
            const size_t repeats = base_repeats/N + 1;
            file << std::setw(15) << N
                 << std::setw(15) << N/vector_addition_aboria_level1(N,repeats)
                 << std::setw(15) << N/vector_addition_aboria_level2(N,repeats)
                 << std::setw(15) << N/vector_addition_eigen(N,repeats)
                 << std::setw(15) << N/vector_addition_stdvector(N,repeats)
                 << std::endl;
        }
        file.close();
    }

    void test_daxpy() {
        std::ofstream file;
        file.open("daxpy.csv");
        const size_t base_repeats = 1e7;
        file <<"#"<< std::setw(14) << "N" 
             << std::setw(15) << "aboria_level1" 
             << std::setw(15) << "aboria_level2" 
             << std::setw(15) << "eigen"
             << std::setw(15) << "stdvector" << std::endl;
#ifdef HAVE_OPENMP
            omp_set_num_threads(1);
#endif
        for (int i = 10; i < 8e6; i*=1.2) {
            const size_t N = i;
            const size_t repeats = base_repeats/N + 1;
            file << std::setw(15) << N
                 << std::setw(15) << N/daxpy_aboria_level1(N,repeats)
                 << std::setw(15) << N/daxpy_aboria_level2(N,repeats)
                 << std::setw(15) << N/daxpy_eigen(N,repeats)
                 << std::setw(15) << N/daxpy_stdvector(N,repeats)
                 << std::endl;
        }
        file.close();

    }





};

#endif /* OPERATORSTEST_H_ */
