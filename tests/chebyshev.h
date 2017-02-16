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


#ifndef CHEBYSHEV_TEST_H_
#define CHEBYSHEV_TEST_H_

#include <cxxtest/TestSuite.h>

#include <random>
#include <time.h>
#include "Level1.h"
#include "Chebyshev.h"

using namespace Aboria;

class ChebyshevTest : public CxxTest::TestSuite {
public:
    template <unsigned int D>
    void helper_Rn_calculation(void) {
        const double tol = 1e-10;
        // randomly generate a bunch of positions over a range 
        std::uniform_real_distribution<double> U(-10,100);
        generator_type generator(time(NULL));
        typedef Vector<double,D> double_d;
        typedef Vector<int,D> int_d;
        const size_t N = 50;
        std::vector<double_d> positions(N);
        for (int i=0; i<N; i++) {
            for (int d=0; d<D; ++d) {
                positions[i][d] = U(generator);
            }
        }
        detail::Chebyshev_Rn<D> Rn;
        for (int n=1; n<10; ++n) {
            Rn.calculate_Sn(std::begin(positions),N,n);
            const int_d start = int_d(0);
            const int_d end = int_d(n-1);
            auto range = iterator_range<lattice_iterator<D>>(
                lattice_iterator<D>(start,end,start)
                ,++lattice_iterator<D>(start,end,end)
                );
            const double_d scale = double_d(1.0)/(Rn.box.bmax-Rn.box.bmin);
            for (int i=0; i<positions.size(); ++i) {
                const double_d &x =  (2*positions[i]-Rn.box.bmin-Rn.box.bmax)*scale;
                for (const int_d& m: range) {
                    TS_ASSERT_DELTA(Rn(m,i),detail::chebyshev_Rn_slow(x,m,n),tol);
                }
            }
        }
    }

    template <unsigned int D>
    void helper_chebyshev_interpolation(void) {
        const double tol = 1e-10;
        // randomly generate a bunch of positions over a range 
        std::uniform_real_distribution<double> U(0,1);
        generator_type generator(time(NULL));
        auto gen = std::bind(U, generator);
        typedef Vector<double,D> double_d;
        typedef Vector<int,D> int_d;
        const size_t N = 50;
        std::vector<double_d> positions(N);
        for (int i=0; i<N; i++) {
            for (int d=0; d<D; ++d) {
                positions[i][d] = gen();
            }
        }

        // randomly generate a source vector
        std::vector<double> source(N);
        std::generate(source.begin(), source.end(), gen);

        const double c = 1.0;
        auto kernel = [&c](const double_d &i, const double_d &j) {
            return std::sqrt((i-j).norm() + c); 
        };

        // perform the operation manually
        std::vector<double> target_manual(N,0.0);
        for (int i=0; i<N; i++) {
            const double_d pi = positions[i];
            for (int j=0; j<N; j++) {
                const double_d pj = positions[j];
                target_manual[i] += kernel(pi,pj)*source[j];
            }
        }

        std::vector<double> target(N);
        for (size_t n = 10; n < 20; ++n) {
            // perform the operation using chebyshev interpolation 
            chebyshev_interpolation<D>(
                                std::begin(source),std::end(source),
                                std::begin(target),std::end(target),
                                std::begin(positions),std::begin(positions),
                                kernel,n);
            const double L2 = std::sqrt(std::inner_product(
                        std::begin(target), std::end(target),
                        std::begin(target_manual),
                        0.0,
                        [](const double t1, const double t2) { return t1 + t2; },
                        [](const double t1, const double t2) { return (t1-t2)*(t1-t2); }
                       ));
            std::cout << "dimension = "<<D<<". n = "<<n<<". L2 error is "<<L2<<std::endl;
        }
    }
 

        


    void test_chebyshev_polynomial_calculation(void) {
        const double tol = 1e-10;
        // evaluate polynomial of order k at i-th root
        // of polynomial of order n
        // should by cos(k*(2*i-1)/(2*n)*PI
        std::cout << "testing polynomial calculation..." << std::endl;
        const int n = 4;
        for (int i=0; i<n; ++i) {
            const double x = cos((2.0*i+1.0)/(2.0*n)*detail::PI);
            for (int k=0; k<n; ++k) {
                    TS_ASSERT_DELTA(detail::chebyshev_polynomial(x,k),
                                    cos(k*(2.0*i+1.0)/(2.0*n)*detail::PI),tol);
            }
        }

    }


    void test_chebyshev_interpolation(void) {
        /*
        std::cout << "testing 1D..." << std::endl;
        helper_chebyshev_interpolation<1>();
        */
        std::cout << "testing 2D..." << std::endl;
        helper_chebyshev_interpolation<2>();
        /*
        std::cout << "testing 3D..." << std::endl;
        helper_chebyshev_interpolation<3>();
        std::cout << "testing 4D..." << std::endl;
        helper_chebyshev_interpolation<4>();
        */
    }


    void test_Rn_calculation(void) {
        std::cout << "testing 1D..." << std::endl;
        helper_Rn_calculation<1>();
        std::cout << "testing 2D..." << std::endl;
        helper_Rn_calculation<2>();
        std::cout << "testing 3D..." << std::endl;
        helper_Rn_calculation<3>();
        std::cout << "testing 4D..." << std::endl;
        helper_Rn_calculation<4>();
    }
};


#endif
