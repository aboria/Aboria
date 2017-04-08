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


#ifndef RBF_STOKES_TEST_H_
#define RBF_STOKES_TEST_H_


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
        auto u_sol= [](const double2& x) { 
            return double2(20*x[0]*std::pow(x[1],3),5*std::pow(x[0],4)-5*std::pow(x[1],4));
        };

        auto p_sol= [](const double2& x) { 
            60*std::pow(x[0],2)*x[1] - 20*std::pow(x[1],3) + c;
        }

        

        auto laplace_funct = [](const double x, const double y) { 
            return -32*std::cos(4*x+4*y);
        };

        ABORIA_VARIABLE(velocity_u,double,"velocity u")
        ABORIA_VARIABLE(velocity_v,double,"velocity v")
        ABORIA_VARIABLE(pressure,double,"pressure p")
        ABORIA_VARIABLE(alpha_1,double,"alpha value")
        ABORIA_VARIABLE(alpha_2,double,"alpha value")
        ABORIA_VARIABLE(alpha_3,double,"alpha value")
        ABORIA_VARIABLE(boundary,uint8_t,"is boundary knot")

    	typedef Particles<std::tuple<alpha_1,alpha_2,alpha_3,velocity_u,velocity_v,pressure,boundary>,2,std::vector,SearchMethod> particles_type;
        typedef particles_type::position position;
        typedef typename ParticlesType::const_reference const_reference;
        typedef Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>> map_type; 
        typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type; 
       	particles_type knots;

       	const double h = 0.5;
        const int max_iter = 100;
        const int restart = 100;
        double2 periodic(false);
        
        const int nx = 7;
        constexpr int N = (nx+1)*(nx+1);
        const double delta = 1.0/nx;
        typename ParticlesType::value_type p;
        for (int i=0; i<=nx; ++i) {
            for (int j=0; j<=nx; ++j) {
                get<position>(p) = double2(i*delta,j*delta);
                if ((i==0)||(i==nx)||(j==0)||(j==nx)) {
                    get<boundary>(p) = true;
                } else {
                    get<boundary>(p) = false;
                }
            }
            knots.push_back(p);
        }
        
        auto kernel_11 = [](const double2& dx, const_reference i, const_reference j) {
            const double r2 = dx.squaredNorm();
            if (r2==0) {
                return 5*26;
            } else {
                const double r4 = r2*r2;
                const double r = std::sqrt(r2);

                return 26*std::pow(r-1,8)*(40*r2 + r*(5+114*r2) + r4*(72-231*r))/r; 
            }
        };

        auto kernel_22 = [](const double2& dx, const_reference i, const_reference j) {
            const double r2 = dx.squaredNorm();
            if (r2==0) {
                return 5*26;
            } else {
                const double r4 = r2*r2;
                const double r = std::sqrt(r2);

                return 26*std::pow(r-1,8)*(40*r2 + r*(-5+18*r2) + r4*(984-3003*r))/r; 
            }
        };

        auto laplace_11 = [](const double2& dx, const_reference i, const_reference j) {
            const double r2 = dx.squaredNorm();
            if (r2==0) {
                return -2*6864;
            } else {
                const double r4 = r2*r2;
                const double r = std::sqrt(r2);

                return 6864*std::pow(r-1,6)*(-12*r2 + r*(-2+3*r2) + r4*(158-147*r))/r; 
            }
        };

        auto laplace_22 = [](const double2& dx, const_reference i, const_reference j) {
            const double r2 = dx.squaredNorm();
            if (r2==0) {
                return -2*6864;
            } else {
                const double r4 = r2*r2;
                const double r = std::sqrt(r2);

                return 6864*std::pow(r-1,6)*(-12*r2 + r*(-2+93*r2) + r4*(698-1617*r))/r; 
            }
        };

        auto grad_wendland_1 = [](const double2& dx, const_reference i, const_reference j) {
            const double r2 = dx.squaredNorm();
            if (r2==0) {
                return 0;
            } else {
                const double r4 = r2*r2;
                const double r = std::sqrt(r2);

                return 26*dx[0]*std::pow(r-1,9)*(45*r2 + 231*r4 + r*(5+159*r2))/r; 
            }
        };

        auto grad_wendland_2 = [](const double2& dx, const_reference i, const_reference j) {
            const double r2 = dx.squaredNorm();
            if (r2==0) {
                return 0;
            } else {
                const double r4 = r2*r2;
                const double r = std::sqrt(r2);

                return 26*dx[1]*std::pow(r-1,9)*(45*r2 + 231*r4 + r*(5+159*r2))/r; 
            }
        };

        const double search_radius = h;

        auto I11 = create_sparse_operator(knots,knots,search_radius,
                [&](const double2& dx,
                    const_reference i,
                    const_reference j) {
                    if (get<boundary>(a)) {
                        return kernel_11(dx,a,b);
                    } else {
                        return -laplace_11(dx,a,b);
                    }
                    });
        auto I12 = create_zero_operator(knots,knots);

        auto I13 = create_zero_operator(knots,knots);

        auto I21 = create_zero_operator(knots,knots);

        auto I22 = create_sparse_operator(knots,knots,search_radius,
                [&](const double2& dx,
                    const_reference i,
                    const_reference j) {
                    if (get<boundary>(a)) {
                        return kernel_22(dx,a,b);
                    } else {
                        return -laplace_22(dx,a,b);
                    }
                    });

        auto I23 = create_zero_operator(knots,knots);

        auto I31 = create_sparse_operator(knots,knots,search_radius,
                [&](const double2& dx,
                    const_reference i,
                    const_reference j) {
                    if (get<boundary>(a)) {
                        return 0;
                    } else {
                        return grad_wendland_1(dx,a,b);
                    }
                    });
        auto I32 = create_sparse_operator(knots,knots,search_radius,
                [&](const double2& dx,
                    const_reference i,
                    const_reference j) {
                    if (get<boundary>(a)) {
                        return 0;
                    } else {
                        return grad_wendland_2(dx,a,b);
                    }
                    });

        auto I33 = create_sparse_operator(knots,knots,search_radius,
                [&](const double2& dx,
                    const_reference i,
                    const_reference j) {
                    if (get<boundary>(a)) {
                        return 0;
                    } else {
                        return grad_wendland_2(dx,a,b);
                    }
                    });



create_zero_operator(knots,knots);

        auto B11 = create_sparse_operator(knots,knots,search_radius,
                [&](const double2& dx,
                    const_reference i,
                    const_reference j) {
                    if (get<boundary>(a)) {
                        return kernel_22(dx,a,b);
                    } else {
                        return laplace_22(dx,a,b);
                    }
                    });




        auto P = create_dense_operator(knots,augment,
                [](const_position_reference dx,
                         const_particle_reference a,
                         const_particle_reference b) {
                    if (get<boundary>(a)) {
                        return 1.0;
                    } else {
                        return 0.0;
                    }
                    });

        auto Pt = create_dense_operator(augment,knots,
                [](const_position_reference dx,
                         const_particle_reference a,
                         const_particle_reference b) {
                    if (get<boundary>(b)) {
                        return 1.0;
                    } else {
                        return 0.0;
                    }
                    });

        auto Zero = create_zero_operator(augment,augment);

        auto W = create_block_operator<6,2>(I11, I12,
                                            I21, I22,
                                            I31, I32,
                                            I11, I12,
                                            I21, I22,
                                            I31, I32);


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


    void test_bucket_search_parallel() {
        helper_Eigen<bucket_search_parallel>();
    }

    void test_bucket_search_serial() {
        helper_Eigen<bucket_search_serial>();
    }


};

#endif /* RBF_PDE_TEST_H_ */
