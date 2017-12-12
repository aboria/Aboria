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


#ifndef RBF_INTERPOLATION_TEST_H_
#define RBF_INTERPOLATION_TEST_H_


#include <cxxtest/TestSuite.h>


//[rbf_interpolation_global
/*`

Here we interpolate a function $f(x,y)$ over the two-dimensional unit cube from $(0,0)$ to $(1,1)$, using the multiquadric Radial Basis Function (RBF). 

The function to be interpolated is

\begin{equation}
f(x,y) = \exp(-9(x-\frac{1}{2})^2 - 9(y-\frac{1}{4})^2), 
\end{equation}

and the multiquadric RBF $K(r)$ is given by

\begin{equation}
K(r) = \sqrt{r^2 + c^2}.
\end{equation}

We create a set of basis functions around a set of $N$ particles with positions $\mathbf{x}_i$. In this case the radial coordinate around each point is $r = ||\mathbf{x}_i - \mathbf{x}||$. The function $f$ is approximated using the sum of these basis functions, with the addition of a constant factor $\beta$

\begin{equation}
\overline{f}(\mathbf{x}) = \beta + \sum_i \alpha_j K(||\mathbf{x}_i-\mathbf{x}||).
\end{equation}

We exactly interpolate the function at $\mathbf{x}_i$ (i.e. set $\overline{f}(\mathbf{x}_i)=f(\mathbf{x}_i)$), leading to

\begin{equation}
f(\mathbf{x}_i) = \beta + \sum_j \alpha_j K(||\mathbf{x}_j-\mathbf{x}_i||).
\end{equation}

Note that the sum $j$ is over the same particle set.

We also need one additional constraint to find $\beta$ 

\begin{equation}
0 = \sum_j \alpha_j.
\end{equation}

We rewrite the two previous equations as a linear algebra expression

\begin{equation}
\mathbf{\Phi} = \begin{bmatrix} \mathbf{G}&\mathbf{P} \\\\ \mathbf{P}^T & \mathbf{0} \end{bmatrix} \begin{bmatrix} \mathbf{\alpha} \\\\ \mathbf{\beta} \end{bmatrix} = \mathbf{W} \mathbf{\gamma},
\end{equation}

where $\mathbf{\Phi} = [f(\mathbf{x}_1),f(\mathbf{x}_2),...,f(\mathbf{x}_N),0]$, $\mathbf{G}$ is an $N \times N$ matrix with elements $G_{ij} = K(||\mathbf{x}_j-\mathbf{x}_i||)$ and $\mathbf{P}$ is a $N \times 1$ vector of ones.

The basis function coefficients $\mathbf{\gamma}$ are found by solving this equation.
To do this, we use the iterative solvers found in the Eigen package and Aboria's ability to wrap C++ function objects as Eigen matricies.

*/
#include "Aboria.h"
#include <random>

using namespace Aboria;


//<-
class RbfInterpolationTest : public CxxTest::TestSuite {
public:

    template<template <typename> class SearchMethod>
    void helper_global(void) {
#ifdef HAVE_EIGEN
        std::cout << "---------------------\n" <<
                      "Running global test....\n" <<
                      "---------------------" << std::endl;
//->
//=int main() {
        auto funct = [](const double x, const double y) { 
            return std::exp(-9*std::pow(x-0.5,2) - 9*std::pow(y-0.25,2)); 
            //return x; 
        };

        ABORIA_VARIABLE(alpha,double,"alpha value")
        ABORIA_VARIABLE(interpolated,double,"interpolated value")
        ABORIA_VARIABLE(constant2,double,"c2 value")

//<-
    	typedef Particles<std::tuple<alpha,constant2,interpolated>,2,std::vector,SearchMethod> ParticlesType;
//->
//=     typedef Particles<std::tuple<alpha,constant2,interpolated>,2> ParticlesType;


        typedef position_d<2> position;
        typedef typename ParticlesType::const_reference const_particle_reference;
        typedef typename position::value_type const & const_position_reference;
        typedef Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>> map_type; 
        typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type; 
        typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_type; 
       	ParticlesType knots;
       	ParticlesType augment;
       	ParticlesType test;

       	const double c = 0.5;
        vdouble2 min(0);
        vdouble2 max(1);
        vdouble2 periodic(false);

        
        const int N = 1000;

        const int nx = 3;
        const int max_iter = 100;
        const int restart = 101;
        const double delta = 1.0/nx;
        typename ParticlesType::value_type p;

        std::default_random_engine generator;
        std::uniform_real_distribution<double> distribution(0.0,1.0);
        for (int i=0; i<N; ++i) {
            get<position>(p) = vdouble2(distribution(generator),
                                       distribution(generator));
            get<constant2>(p) = std::pow(c,2);  
            knots.push_back(p);

            get<position>(p) = vdouble2(distribution(generator),
                                       distribution(generator));
            get<constant2>(p) = std::pow(c,2);  
            test.push_back(p);
        }

        //knots.init_neighbour_search(min,max,periodic);

        augment.push_back(p);

        auto kernel = [](const_position_reference dx,
                         const_particle_reference a,
                         const_particle_reference b) {
                            return std::sqrt(dx.squaredNorm() + get<constant2>(b));
                        };

        auto one = [](const_position_reference dx,
                      const_particle_reference a,
                      const_particle_reference b) {
                            return 1.0;
                        };

        auto G = create_dense_operator(knots,knots,kernel);
        auto P = create_dense_operator(knots,augment,one);
        auto Pt = create_dense_operator(augment,knots,one);
        auto Zero = create_zero_operator(augment,augment);

        auto W = create_block_operator<2,2>(G, P,
                                            Pt,Zero);

        auto G_test = create_dense_operator(test,knots,kernel);
        auto W_test = create_block_operator<2,2>(G_test, P,
                                                 Pt,     Zero);

        vector_type phi(N+1), gamma(N+1);
        for (int i=0; i<knots.size(); ++i) {
            const double x = get<position>(knots[i])[0];
            const double y = get<position>(knots[i])[1];
            phi[i] = funct(x,y);
        }
        phi[knots.size()] = 0;

        matrix_type W_matrix(N+1,N+1);
        W.assemble(W_matrix);

        gamma = W_matrix.ldlt().solve(phi);

        Eigen::GMRES<matrix_type> gmres;
        gmres.setMaxIterations(max_iter);
        gmres.set_restart(restart);
        gmres.compute(W_matrix);
        gamma = gmres.solve(phi);
        std::cout << "GMRES:       #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;

        phi = W*gamma;
        double rms_error = 0;
        double scale = 0;
        for (int i=0; i<knots.size(); ++i) {
            const double x = get<position>(knots[i])[0];
            const double y = get<position>(knots[i])[1];
            const double truth = funct(x,y);
            const double eval_value = phi[i];
            rms_error += std::pow(eval_value-truth,2);
            scale += std::pow(truth,2);
            //TS_ASSERT_DELTA(eval_value,truth,2e-3); 
        }

        std::cout << "rms_error for global support, at centers  = "<<std::sqrt(rms_error/scale)<<std::endl;
        TS_ASSERT_LESS_THAN(std::sqrt(rms_error/scale),1e-6);
        
        phi = W_test*gamma;
        rms_error = 0;
        scale = 0;
        for (int i=0; i<test.size(); ++i) {
            const double x = get<position>(test[i])[0];
            const double y = get<position>(test[i])[1];
            const double truth = funct(x,y);
            const double eval_value = phi[i];
            rms_error += std::pow(eval_value-truth,2);
            scale += std::pow(truth,2);
            //TS_ASSERT_DELTA(eval_value,truth,2e-3); 
        }
        std::cout << "rms_error for global support, away from centers  = "<<std::sqrt(rms_error/scale)<<std::endl;
        TS_ASSERT_LESS_THAN(std::sqrt(rms_error/scale),1e-5);


//=}
//]
#endif // HAVE_EIGEN
    }


template<template <typename> class SearchMethod>
void helper_compact(void) {
#ifdef HAVE_EIGEN
        std::cout << "---------------------\n" <<
                     "Running compact test....\n" <<
                     "---------------------" << std::endl;

//[rbf_interpolation_compact
//=#include "Aboria.h"
//=#include <random>
//=int main() {
        auto funct = [](const double x, const double y) { 
            return std::exp(-9*std::pow(x-0.5,2) - 9*std::pow(y-0.25,2)); 
            //return x; 
        };

        ABORIA_VARIABLE(alpha,double,"alpha value")
        ABORIA_VARIABLE(interpolated,double,"interpolated value")
        ABORIA_VARIABLE(constant2,double,"c2 value")

    	typedef Particles<std::tuple<alpha,constant2,interpolated>,2,std::vector,SearchMethod> ParticlesType;
        typedef position_d<2> position;
        typedef typename ParticlesType::const_reference const_particle_reference;
        typedef typename position::value_type const & const_position_reference;
        typedef Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>> map_type; 
        typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type; 
        typedef Eigen::SparseMatrix<double> matrix_type; 
       	ParticlesType knots,augment;
       	ParticlesType test;

        const double hfac = 4.0;
        vdouble2 min(0);
        vdouble2 max(1);
        vdouble2 periodic(false);
        
        const int N = 1000;
        
        const int max_iter = 100;
        const int restart = 101;
        const double delta = std::pow(double(N),-1.0/2.0);
       	const double h = hfac*delta;
        std::cout << "using h = "<<h<<std::endl;
        const double RASM_size = 2*h;
        const int RASM_n = N*std::pow(RASM_size,2)/(max-min).prod();
        const double RASM_buffer = 0.9*RASM_size;
        std::cout << "RASM_size = "<<RASM_size<<" RASM_n = "<<RASM_n<<" RASM_buffer = "<<RASM_buffer<<std::endl;

        typename ParticlesType::value_type p;

        std::default_random_engine generator(123);
        std::uniform_real_distribution<double> distribution(0.0,1.0);
        for (int i=0; i<N; ++i) {
            get<position>(p) = vdouble2(distribution(generator),
                                       distribution(generator));
            knots.push_back(p);

            get<position>(p) = vdouble2(distribution(generator),
                                       distribution(generator));
            test.push_back(p);
        }
        augment.push_back(p);

	    knots.init_neighbour_search(min,max,periodic);

        auto kernel = [h](const_position_reference dx,
                         const_particle_reference a,
                         const_particle_reference b) {
                            return std::pow(2.0-dx.norm()/h,4)*(1.0+2.0*dx.norm()/h);
                        };

        auto one = [](const_position_reference dx,
                      const_particle_reference a,
                      const_particle_reference b) {
                            return 1.0;
                        };

        auto G = create_sparse_operator(knots,knots,2*h,kernel);
        auto P = create_dense_operator(knots,augment,one);
        auto Pt = create_dense_operator(augment,knots,one);
        auto Zero = create_zero_operator(augment,augment);

        auto W = create_block_operator<2,2>(G, P,
                                            Pt,Zero);

        auto G_test = create_sparse_operator(test,knots,2*h,kernel);
        auto W_test = create_block_operator<2,2>(G_test, P,
                                                 Pt,     Zero);


        vector_type phi(N+1), gamma(N+1);
        for (size_t i=0; i<knots.size(); ++i) {
            const double x = get<position>(knots[i])[0];
            const double y = get<position>(knots[i])[1];
            phi[i] = funct(x,y);
        }
        phi[knots.size()] = 0;

        matrix_type W_matrix(N+1,N+1);
        W.assemble(W_matrix);

        Eigen::ConjugateGradient<matrix_type, Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double>> cg_test;
        cg_test.setMaxIterations(max_iter);
        cg_test.compute(W_matrix);
        gamma = cg_test.solve(phi);
        std::cout << "CG:          #iterations: " << cg_test.iterations() << ", estimated error: " << cg_test.error() << std::endl;

        Eigen::ConjugateGradient<matrix_type, 
            Eigen::Lower|Eigen::Upper, RASMPreconditioner<Eigen::HouseholderQR>> cg;
        cg.setMaxIterations(max_iter);
        cg.preconditioner().set_buffer_size(RASM_buffer);
        cg.preconditioner().set_number_of_particles_per_domain(RASM_n);
        cg.preconditioner().analyzePattern(W);
        cg.compute(W_matrix);
        gamma = cg.solve(phi);
        std::cout << "CG-RASM:     #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;

        
        Eigen::MINRES<matrix_type, Eigen::Lower|Eigen::Upper,  RASMPreconditioner<Eigen::HouseholderQR>> minres;
        minres.setMaxIterations(max_iter);
        minres.preconditioner().set_buffer_size(RASM_buffer);
        minres.preconditioner().set_number_of_particles_per_domain(RASM_n);
        minres.preconditioner().analyzePattern(W);
        minres.compute(W_matrix);
        gamma = minres.solve(phi);
        std::cout << "MINRES-RASM: #iterations: " << minres.iterations() << ", estimated error: " << minres.error() << std::endl;

        Eigen::GMRES<matrix_type,  RASMPreconditioner<Eigen::HouseholderQR>> gmres;
        gmres.setMaxIterations(max_iter);
        gmres.preconditioner().set_buffer_size(RASM_buffer);
        gmres.preconditioner().set_number_of_particles_per_domain(RASM_n);
        gmres.preconditioner().analyzePattern(W);
        gmres.set_restart(restart);
        gmres.compute(W_matrix);
        gamma = gmres.solve(phi);
        std::cout << "GMRES-RASM:  #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;

        Eigen::DGMRES<matrix_type,  RASMPreconditioner<Eigen::HouseholderQR>> dgmres;
        dgmres.setMaxIterations(max_iter);
        dgmres.preconditioner().set_buffer_size(RASM_buffer);
        dgmres.preconditioner().set_number_of_particles_per_domain(RASM_n);
        dgmres.preconditioner().analyzePattern(W);
        dgmres.set_restart(restart);
        dgmres.compute(W_matrix);
        gamma = dgmres.solve(phi);
        std::cout << "DGMRES-RASM:  #iterations: " << dgmres.iterations() << ", estimated error: " << dgmres.error() << std::endl;

        Eigen::BiCGSTAB<matrix_type, RASMPreconditioner<Eigen::HouseholderQR>> bicg;
        bicg.setMaxIterations(max_iter);
        bicg.preconditioner().set_buffer_size(RASM_buffer);
        bicg.preconditioner().set_number_of_particles_per_domain(RASM_n);
        bicg.preconditioner().analyzePattern(W);
        bicg.compute(W_matrix);
        gamma = bicg.solve(phi);
        std::cout << "BiCGSTAB-RASM:#iterations: " << bicg.iterations() << ", estimated error: " << bicg.error() << std::endl;


        double rms_error = 0;
        double scale = 0;
        phi = W*gamma;
        for (size_t i=0; i<knots.size(); ++i) {
            const double x = get<position>(knots[i])[0];
            const double y = get<position>(knots[i])[1];
            const double truth = funct(x,y);
            const double eval_value = phi[i];
            rms_error += std::pow(eval_value-truth,2);
            scale += std::pow(truth,2);
            //TS_ASSERT_DELTA(eval_value,truth,2e-3); 
        }
        std::cout << "rms_error for compact support, at centers  = "<<std::sqrt(rms_error/scale)<<std::endl;
        TS_ASSERT_LESS_THAN(std::sqrt(rms_error/scale),1e-4);

        rms_error = 0;
        scale = 0;
        phi = W_test*gamma;
        for (size_t i=0; i<test.size(); ++i) {
            const double x = get<position>(test[i])[0];
            const double y = get<position>(test[i])[1];
            const double truth = funct(x,y);
            const double eval_value = phi[i];
            rms_error += std::pow(eval_value-truth,2);
            scale += std::pow(truth,2);
            //TS_ASSERT_DELTA(eval_value,truth,2e-3); 
        }

        std::cout << "rms_error for compact support, away from centers  = "<<std::sqrt(rms_error/scale)<<std::endl;
        TS_ASSERT_LESS_THAN(std::sqrt(rms_error/scale),1e-2);

//=}
//]
#endif // HAVE_EIGEN
    }

template<template <typename> class SearchMethod>
    void helper_h2(void) {
#ifdef HAVE_EIGEN
        std::cout << "---------------------\n"<<
                     "Running h2 test....\n" <<
                     "---------------------" << std::endl;
//[rbf_interpolation_h2
//=#include "Aboria.h"
//=#include <random>
//=int main() {
        auto funct = [](const double x, const double y) { 
            return std::exp(-9*std::pow(x-0.5,2) - 9*std::pow(y-0.25,2)); 
            //return x; 
        };

        ABORIA_VARIABLE(alpha,double,"alpha value")
        ABORIA_VARIABLE(interpolated,double,"interpolated value")

    	typedef Particles<std::tuple<alpha,interpolated>,2,std::vector,SearchMethod> ParticlesType;
        typedef position_d<2> position;
        typedef typename ParticlesType::const_reference const_particle_reference;
        typedef typename position::value_type const & const_position_reference;
        typedef Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>> map_type; 
        typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type; 
        typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_type; 
       	ParticlesType knots;
       	ParticlesType augment;
       	ParticlesType test;

       	const double c = 3.0;
        const double c2 = std::pow(c,2);
        vdouble2 min(0);
        vdouble2 max(1);
        vdouble2 periodic(false);

        const int N = 10000;

        const double RASM_size = 0.3/c;
        const int RASM_n = N*std::pow(RASM_size,2)/(max-min).prod();
        const double RASM_buffer = 0.9*RASM_size;
        //std::cout << "RASM_size = "<<RASM_size<<" RASM_n = "<<RASM_n<<" RASM_buffer = "<<RASM_buffer<<std::endl;

        const int nx = 3;
        const int max_iter = 100;
        const int restart = 101;
        const double delta = 1.0/nx;
        typename ParticlesType::value_type p;

        std::default_random_engine generator;
        std::uniform_real_distribution<double> distribution(0.0,1.0);
        for (int i=0; i<N; ++i) {
            get<position>(p) = vdouble2(distribution(generator),
                                       distribution(generator));
            knots.push_back(p);

            get<position>(p) = vdouble2(distribution(generator),
                                       distribution(generator));
            test.push_back(p);
        }

        knots.init_neighbour_search(min,max,periodic);

        augment.push_back(p);

        auto kernel = [&](const vdouble2& dx,
                         const vdouble2& a,
                         const vdouble2& b) {
                            return std::exp(-dx.squaredNorm()*c2);
                        };

        auto one = [](const_position_reference dx,
                      const_particle_reference a,
                      const_particle_reference b) {
                            return 1.0;
                        };

        auto G = create_h2_operator<4>(knots,knots,kernel);
        auto P = create_dense_operator(knots,augment,one);
        auto Pt = create_dense_operator(augment,knots,one);
        auto Zero = create_zero_operator(augment,augment);

        auto W = create_block_operator<2,2>(G, P,
                                            Pt,Zero);

        auto Gtest = create_h2_operator(G.get_first_kernel(),test);
        auto Ptest = create_dense_operator(test,augment,one);
        auto Wtest = create_block_operator<2,2>(Gtest, Ptest,
                                                Pt,Zero);

        vector_type phi(N+1), gamma(N+1);
        for (int i=0; i<knots.size(); ++i) {
            const double x = get<position>(knots[i])[0];
            const double y = get<position>(knots[i])[1];
            phi[i] = funct(x,y);
        }
        phi[knots.size()] = 0;

        //Eigen::DGMRES<decltype(W),  RASMPreconditioner<Eigen::HouseholderQR>> dgmres;
        Eigen::DGMRES<decltype(W)> dgmres;
        dgmres.setMaxIterations(max_iter);
        //dgmres.preconditioner().set_buffer_size(RASM_buffer);
        //dgmres.preconditioner().set_number_of_particles_per_domain(RASM_n);
        dgmres.set_restart(restart);
        dgmres.compute(W);
        gamma = dgmres.solve(phi);
        std::cout << "DGMRES-RASM:  #iterations: " << dgmres.iterations() << ", estimated error: " << dgmres.error() << std::endl;

        phi = W*gamma;
        double rms_error = 0;
        double scale = 0;
        for (int i=0; i<knots.size(); ++i) {
            const double x = get<position>(knots[i])[0];
            const double y = get<position>(knots[i])[1];
            const double truth = funct(x,y);
            const double eval_value = phi[i];
            rms_error += std::pow(eval_value-truth,2);
            scale += std::pow(truth,2);
            //TS_ASSERT_DELTA(eval_value,truth,2e-3); 
        }

        std::cout << "rms_error for global support, at centers  = "<<std::sqrt(rms_error/scale)<<std::endl;
        TS_ASSERT_LESS_THAN(std::sqrt(rms_error/scale),1e-4);
        
        rms_error = 0;
        scale = 0;
        vector_type phi_test = Wtest*gamma;
        for (int i=0; i<test.size(); ++i) {
            const vdouble2 p = get<position>(test)[i]; 
            const double eval_value = phi_test[i];
            const double truth = funct(p[0],p[1]);
            rms_error += std::pow(eval_value-truth,2);
            scale += std::pow(truth,2);
            //TS_ASSERT_DELTA(eval_value,truth,2e-3); 
        }
        std::cout << "rms_error for global support, away from centers  = "<<std::sqrt(rms_error/scale)<<std::endl;
        TS_ASSERT_LESS_THAN(std::sqrt(rms_error/scale),1e-4);


//=}
//]
#endif // HAVE_EIGEN
    }



    void test_bucket_search_parallel() {
        std::cout << "-------------------------------------------\n"<<
                     "Running tests on bucket_search_parallel....\n" <<
                     "------------------------------------------" << std::endl;
        helper_global<bucket_search_parallel>();
        helper_compact<bucket_search_parallel>();
    }

    void test_bucket_search_serial() {
        std::cout << "-------------------------------------------\n"<<
                     "Running tests on bucket_search_serial....\n" <<
                     "------------------------------------------" << std::endl;
        helper_global<bucket_search_serial>();
        helper_compact<bucket_search_serial>();
    }

    void test_kdtree() {
        std::cout << "-------------------------------------------\n"<<
                     "Running tests on kdtree....\n" <<
                     "------------------------------------------" << std::endl;
        helper_compact<nanoflann_adaptor>();
        helper_h2<nanoflann_adaptor>();
    }

    void test_octtree() {
        std::cout << "-------------------------------------------\n"<<
                     "Running tests on octtree....\n" <<
                     "------------------------------------------" << std::endl;
        helper_compact<nanoflann_adaptor>();
        helper_h2<octtree>();
    }

};

#endif /* RBF_INTERPOLATION_TEST_H_ */
//->
