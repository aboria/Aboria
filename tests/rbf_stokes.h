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


//[rbf_stokes
#include "Aboria.h"

using namespace Aboria;

//<-
class RbfStokesTest : public CxxTest::TestSuite {
public:

    template<template <typename> class SearchMethod>
    void helper_Eigen(void) {
#ifdef HAVE_EIGEN
//->
//=int main() {
        auto u_sol= [](const vdouble2& x) { 
            return vdouble2(20*x[0]*std::pow(x[1],3),5*std::pow(x[0],4)-5*std::pow(x[1],4));
        };

        auto p_sol = [](const vdouble2& x) { 
            return 60*std::pow(x[0],2)*x[1] - 20*std::pow(x[1],3);
        };

        auto grad_p_sol = [](const vdouble2& x) { 
            return vdouble2(120*x[0]*x[1] , 60*std::pow(x[0],2) - 60*std::pow(x[1],2));
        };

        ABORIA_VARIABLE(velocity_u,double,"velocity u")
        ABORIA_VARIABLE(velocity_v,double,"velocity v")
        ABORIA_VARIABLE(pressure_x,double,"pressure gradient x")
        ABORIA_VARIABLE(pressure_y,double,"pressure gradient y")
        ABORIA_VARIABLE(boundary,uint8_t,"is boundary knot")

    	typedef Particles<std::tuple<velocity_u,velocity_v,pressure_x,pressure_y,boundary>,2,std::vector,SearchMethod> particles_type;
        typedef typename particles_type::position position;
        typedef typename particles_type::const_reference const_reference;
        typedef Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>> map_type; 
        typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type; 
        typedef Eigen::SparseMatrix<double> matrix_type; 
       	particles_type knots;

        const double mu = 1;
       	const double h = 10.0;
        const double invh = 1.0/h;
        vdouble2 periodic(false);
        vdouble2 low(-0.5);
        vdouble2 high(1.5);
        
        const int nx = 16;
        constexpr int N = (nx+1)*(nx+1);
        const double delta = 1.0/nx;
        typename particles_type::value_type p;
        for (int i=0; i<=nx; ++i) {
            for (int j=0; j<=nx; ++j) {
                get<position>(p) = vdouble2(i*delta,j*delta);
                if ((i==0)||(i==nx)||(j==0)||(j==nx)) {
                    get<boundary>(p) = true;
                } else {
                    get<boundary>(p) = false;
                }
                knots.push_back(p);
            }
        }

        knots.init_neighbour_search(low,high,periodic);
        std::cout << "total number of knots = "<<knots.size()<<std::endl;
        
        auto kernel_x = [&](const vdouble2& dx) {
            const double r = dx.norm()*invh;
            return -26.0*dx[0]*std::pow(invh,2)*std::pow(r-1.0,9)*(5.0 + 3.0*r*(15.0 + r*(53.0+77.0*r))); 
        };
        auto kernel_y = [&](const vdouble2& dx) {
            const double r = dx.norm()*invh;
            return -26.0*dx[1]*std::pow(invh,2)*std::pow(r-1.0,9)*(5.0 + 3.0*r*(15.0 + r*(53.0+77.0*r))); 
        };
        auto kernel_xx = [&](const vdouble2& dx) {
            const double r = dx.norm()*invh;
            return 26.0*std::pow(invh,2)*std::pow(r-1.0,8)*((-1.0 + r)*(5.0+3.0*r*(15.0+r*(53.0+77.0*r))) 
                                        + 132.0*(1.0+r*(8.0+21.0*r))*dx[0]*dx[0]*invh*invh); 
        };
        auto kernel_yy = [&](const vdouble2& dx) {
            const double r = dx.norm()*invh;
            return 26.0*std::pow(invh,2)*std::pow(r-1.0,8)*((-1.0 + r)*(5.0+3.0*r*(15.0+r*(53.0+77.0*r))) 
                                        + 132.0*(1.0+r*(8.0+21.0*r))*dx[1]*dx[1]*invh*invh); 
        };
        auto kernel_xy = [&](const vdouble2& dx) {
            const double r = dx.norm()*invh;
            return 3432.0*dx[0]*dx[1]*std::pow(invh,4)*std::pow(r-1.0,8)*(1.0 + r*(8.0 + 21.0*r)); 
        };
        auto laplace_xx = [&](const vdouble2& dx) {
            const double r = dx.norm()*invh;
            const double r2 = r*r;
            const double r3 = r2*r;
            const double r4 = r2*r2;
            return 6864.0*std::pow(invh,4)*std::pow(r-1.0,6)*(2.0+12.0*r-3.0*r2-158.0*r3+147.0*r4
                                    +30.0*(-3.0+r*(-18.0+49.0*r))*dx[0]*dx[0]*invh*invh);
        };
        auto laplace_yy = [&](const vdouble2& dx) {
            const double r = dx.norm()*invh;
            const double r2 = r*r;
            const double r3 = r2*r;
            const double r4 = r2*r2;
            return 6864.0*std::pow(invh,4)*std::pow(r-1.0,6)*(2.0+12.0*r-3.0*r2-158.0*r3+147.0*r4
                                    +30.0*(-3.0+r*(-18.0+49.0*r))*dx[1]*dx[1]*invh*invh);
        };
        auto laplace_xy = [&](const vdouble2& dx) {
            const double r = dx.norm()*invh;
            return 205920.0*dx[0]*dx[1]*std::pow(invh,6)*std::pow(r-1.0,6)*(-3.0+r*(-18.0+49.0*r));
        };
        auto laplace2_xx = [&](const vdouble2& dx) {
            const double r = dx.norm()*invh;
            const double r2 = r*r;
            const double r3 = r2*r;
            const double r4 = r2*r2;
            return 2471040.0*std::pow(invh,6)*std::pow(r-1.0,4)*(-1.0-4.0*r+46.0*r2-90.0*r3+49.0*r4
                                +14.0*(8.0+r*(-31.0+28.0*r))*dx[0]*dx[0]*invh*invh);
        };
        auto laplace2_yy = [&](const vdouble2& dx) {
            const double r = dx.norm()*invh;
            const double r2 = r*r;
            const double r3 = r2*r;
            const double r4 = r2*r2;
            return 2471040.0*std::pow(invh,6)*std::pow(r-1.0,4)*(-1.0-4.0*r+46.0*r2-90.0*r3+49.0*r4
                                +14.0*(8.0+r*(-31.0+28.0*r))*dx[1]*dx[1]*invh*invh);
        };
        auto laplace2_xy = [&](const vdouble2& dx) {
            const double r = dx.norm()*invh;
            return 34594560.0*dx[0]*dx[1]*std::pow(invh,8)*std::pow(r-1.0,4)*(8.0+r*(-31.0+28.0*r));
        };

        const double search_radius = h;

        auto A11 = create_sparse_operator(knots,knots,search_radius,
                [&](const vdouble2& dx,
                    const_reference i,
                    const_reference j) {
                    if (get<boundary>(i)) {
                        if (get<boundary>(j)) {
                            return -kernel_yy(dx);
                        } else {
                            return mu*laplace_yy(dx);
                        }
                    } else {
                        if (get<boundary>(j)) {
                            return mu*laplace_yy(dx);
                        } else {
                            return -mu*mu*laplace2_yy(dx) - kernel_xx(dx);
                        }
                    }
               });

        auto A12 = create_sparse_operator(knots,knots,search_radius,
                [&](const vdouble2& dx,
                    const_reference i,
                    const_reference j) {
                    if (get<boundary>(i)) {
                        if (get<boundary>(j)) {
                            return kernel_xy(dx);
                        } else {
                            return -mu*laplace_xy(dx);
                        }
                    } else {
                        if (get<boundary>(j)) {
                            return -mu*laplace_xy(dx);
                        } else {
                            return mu*mu*laplace2_xy(dx) - kernel_xy(dx);
                        }
                    }
               });

        auto A22 = create_sparse_operator(knots,knots,search_radius,
                [&](const vdouble2& dx,
                    const_reference i,
                    const_reference j) {
                    if (get<boundary>(i)) {
                        if (get<boundary>(j)) {
                            return -kernel_xx(dx);
                        } else {
                            return mu*laplace_xx(dx);
                        }
                    } else {
                        if (get<boundary>(j)) {
                            return mu*laplace_xx(dx);
                        } else {
                            return -mu*mu*laplace2_xx(dx) - kernel_yy(dx);
                        }
                    }
               });


        auto A = create_block_operator<2,2>(A11, A12,
                                            A12, A22);


        vector_type source(2*N), alpha(2*N);
        for (int i=0; i<N; ++i) {
            if (get<boundary>(knots[i])) {
                const vdouble2& velocity_solution = u_sol(get<position>(knots[i]));
                source[i] = velocity_solution[0];
                source[i+N] = velocity_solution[1];
            } else {
                source[i] = 0;
                source[i+N] = 0;
            }
        }
        std::cout << std::endl;

        // assemble matrix
        matrix_type A_matrix(2*N,2*N);
        A.assemble(A_matrix);
        std::cout << "finished assembling nonzeros = "<<A_matrix.nonZeros()<<std::endl;

        // calculate velocity weights
        Eigen::SparseLU<matrix_type> solver;
        solver.compute(A_matrix);
        if(solver.info()!=Eigen::Success) {
            std::cout << "decomp failed" <<std::endl;
            return;
        }
        alpha = solver.solve(source);
        if(solver.info()!=Eigen::Success) {
            std::cout << "solve failed" <<std::endl;
            return;
        }

        std::cout << "The relative error is:\n" 
             << (A_matrix*alpha - source).norm() / source.norm() << std::endl;

        // evaluate solution
        auto B11 = create_sparse_operator(knots,knots,search_radius,
                [&](const vdouble2& dx,
                    const_reference i,
                    const_reference j) {
                    if (get<boundary>(j)) {
                        return -kernel_yy(dx);
                    } else {
                        return mu*laplace_yy(dx);
                    }
               });

        auto B12 = create_sparse_operator(knots,knots,search_radius,
                [&](const vdouble2& dx,
                    const_reference i,
                    const_reference j) {
                    if (get<boundary>(j)) {
                        return kernel_xy(dx);
                    } else {
                        return -mu*laplace_xy(dx);
                    }
               });

        auto B22 = create_sparse_operator(knots,knots,search_radius,
                [&](const vdouble2& dx,
                    const_reference i,
                    const_reference j) {
                    if (get<boundary>(j)) {
                        return -kernel_xx(dx);
                    } else {
                        return mu*laplace_xx(dx);
                    }
               });

        auto B31 = create_sparse_operator(knots,knots,search_radius,
                [&](const vdouble2& dx,
                    const_reference i,
                    const_reference j) {
                    if (get<boundary>(j)) {
                        return 0.0;
                    } else {
                        return -kernel_xx(dx);
                    }
               });

        auto B32 = create_sparse_operator(knots,knots,search_radius,
                [&](const vdouble2& dx,
                    const_reference i,
                    const_reference j) {
                    if (get<boundary>(j)) {
                        return 0.0;
                    } else {
                        return -kernel_xy(dx);
                    }
               });

        auto B42 = create_sparse_operator(knots,knots,search_radius,
                [&](const vdouble2& dx,
                    const_reference i,
                    const_reference j) {
                    if (get<boundary>(j)) {
                        return 0.0;
                    } else {
                        return -kernel_yy(dx);
                    }
               });

        map_type u(get<velocity_u>(knots).data(),N);
        map_type v(get<velocity_v>(knots).data(),N);
        map_type prx(get<pressure_x>(knots).data(),N);
        map_type pry(get<pressure_y>(knots).data(),N);

        u = B11*alpha.head(N) + B12*alpha.tail(N);
        v = B12*alpha.head(N) + B22*alpha.tail(N);
        prx = B31*alpha.head(N) + B32*alpha.tail(N);
        pry = B32*alpha.head(N) + B42*alpha.tail(N);

        vdouble2 L2(0);
        vdouble2 scale(0);
        for (int i=0; i<N; ++i) {
            const double x = get<position>(knots[i])[0];
            const double y = get<position>(knots[i])[1];
            const vdouble2& velocity_solution = u_sol(get<position>(knots[i]));
            const vdouble2& pressure_solution = grad_p_sol(get<position>(knots[i]));
            L2[0] += (vdouble2(u[i],v[i])-velocity_solution).squaredNorm();
            L2[1] += (vdouble2(prx[i],pry[i])-pressure_solution).squaredNorm();
            scale[0] += velocity_solution.squaredNorm();
            scale[1] += pressure_solution.squaredNorm();
        }
        TS_ASSERT_DELTA(std::sqrt(L2[0]/scale[0]),0,1e-5); 
        TS_ASSERT_DELTA(std::sqrt(L2[1]/scale[1]),0,3e-2); 
        std::cout << "rms errors (u,p) = ("
                  <<std::sqrt(L2[0]/scale[0])<<","
                  <<std::sqrt(L2[1]/scale[1])<<")"<<std::endl;

#ifdef HAVE_VTK
        vtkWriteGrid("rbf_stokes",0,knots.get_grid(true));
#endif


//=}
//]
#endif // HAVE_EIGEN
    }


    void test_CellListOrdered() {
        helper_Eigen<CellListOrdered>();
    }

    void test_CellList() {
        helper_Eigen<CellList>();
    }


};

#endif /* RBF_PDE_TEST_H_ */
