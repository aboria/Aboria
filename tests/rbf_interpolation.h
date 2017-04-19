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
#include "Aboria.h"
#include <random>

using namespace Aboria;


//<-
class RbfInterpolationTest : public CxxTest::TestSuite {
public:

    template<template <typename> class SearchMethod>
    void helper_global(void) {
#ifdef HAVE_EIGEN
//->
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
        typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_type; 
       	ParticlesType knots;
       	ParticlesType augment;
       	ParticlesType test;

       	const double c = 0.5;
        double2 min(0);
        double2 max(1);
        double2 periodic(false);
        
        const int N = 100;
        const int nx = 3;
        const int max_iter = 100;
        const int restart = 100;
        const double delta = 1.0/nx;
        typename ParticlesType::value_type p;

        std::default_random_engine generator;
        std::uniform_real_distribution<double> distribution(0.0,1.0);
        for (int i=0; i<N; ++i) {
            get<position>(p) = double2(distribution(generator),
                                       distribution(generator));
            get<constant2>(p) = std::pow(c,2);  
            knots.push_back(p);

            get<position>(p) = double2(distribution(generator),
                                       distribution(generator));
            get<constant2>(p) = std::pow(c,2);  
            test.push_back(p);
        }

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

        vector_type phi(N+1), gamma(N+1);
        for (int i=0; i<knots.size(); ++i) {
            const double x = get<position>(knots[i])[0];
            const double y = get<position>(knots[i])[1];
            phi[i] = funct(x,y);
        }
        phi[knots.size()] = 0;

        matrix_type W_matrix(N+1,N+1);
        W.assemble(W_matrix);

        Eigen::ConjugateGradient<matrix_type, 
            Eigen::Lower|Eigen::Upper,  RASMPreconditioner> cg;
        cg.setMaxIterations(max_iter);
        cg.analyzePattern(W);
        cg.factorize(W_matrix);
        gamma = cg.solve(phi);
        std::cout << std::endl << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;

        Eigen::BiCGSTAB<matrix_type,RASMPreconditioner> bicg;
        bicg.setMaxIterations(max_iter);
        bicg.analyzePattern(W);
        bicg.factorize(W_matrix);
        gamma = bicg.solve(phi);
        std::cout << "BiCGSTAB: #iterations: " << bicg.iterations() << ", estimated error: " << bicg.error() << std::endl;

        Eigen::MINRES<matrix_type, Eigen::Lower|Eigen::Upper, RASMPreconditioner> minres;
        minres.setMaxIterations(max_iter);
        minres.analyzePattern(W);
        minres.factorize(W_matrix);
        gamma = minres.solve(phi);
        std::cout << "MINRES:   #iterations: " << minres.iterations() << ", estimated error: " << minres.error() << std::endl;

        Eigen::GMRES<matrix_type, RASMPreconditioner> gmres;
        gmres.setMaxIterations(max_iter);
        gmres.analyzePattern(W);
        gmres.set_restart(restart);
        gmres.factorize(W_matrix);
        gamma = gmres.solve(phi);
        std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;

        Eigen::DGMRES<matrix_type, RASMPreconditioner> dgmres;
        dgmres.setMaxIterations(max_iter);
        dgmres.set_restart(restart);
        dgmres.analyzePattern(W);
        dgmres.factorize(W_matrix);
        gamma = dgmres.solve(phi);
        std::cout << "DGMRES:   #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;


        phi = W*gamma;
        for (int i=0; i<knots.size(); ++i) {
            const double x = get<position>(knots[i])[0];
            const double y = get<position>(knots[i])[1];
            phi[i] = funct(x,y);
            TS_ASSERT_DELTA(phi[i],funct(x,y),2e-3); 
        }
        TS_ASSERT_DELTA(phi[knots.size()],0,2e-3); 


        // wrap Aboria vectors as Eigen vectors
        map_type alpha_wrap(get<alpha>(knots).data(),N);
        map_type interp_wrap(get<interpolated>(knots).data(),N);

        alpha_wrap = gamma.segment<N>(0);
        vector_type beta = vector_type::Constant(N,gamma[N]);
        interp_wrap = G*alpha_wrap + beta;

        double rms_error = 0;
        double scale = 0;
        for (int i=0; i<knots.size(); ++i) {
            const double x = get<position>(knots[i])[0];
            const double y = get<position>(knots[i])[1];
            const double truth = funct(x,y);
            const double eval_value = get<interpolated>(knots[i]);
            rms_error += std::pow(eval_value-truth,2);
            scale += std::pow(truth,2);
            TS_ASSERT_DELTA(eval_value,truth,2e-3); 
        }

        std::cout << "rms_error for global support, at centers  = "<<std::sqrt(rms_error/scale)<<std::endl;
        TS_ASSERT_LESS_THAN(std::sqrt(rms_error/scale),1e-8);

        map_type interp_wrap_test(get<interpolated>(test).data(),N);
        auto G_test = create_dense_operator(test,knots,kernel);
        interp_wrap_test = G_test*alpha_wrap + beta;
                
        rms_error = 0;
        scale = 0;
        for (int i=0; i<test.size(); ++i) {
            const double x = get<position>(test[i])[0];
            const double y = get<position>(test[i])[1];
            const double truth = funct(x,y);
            const double eval_value = get<interpolated>(test[i]);
            rms_error += std::pow(eval_value-truth,2);
            scale += std::pow(truth,2);
            TS_ASSERT_DELTA(eval_value,truth,2e-3); 
        }
        std::cout << "rms_error for global support, away from centers  = "<<std::sqrt(rms_error/scale)<<std::endl;
        TS_ASSERT_LESS_THAN(std::sqrt(rms_error/scale),1.5e-4);

//=}
//]
#endif // HAVE_EIGEN
    }

template<template <typename> class SearchMethod>
void helper_compact(void) {
#ifdef HAVE_EIGEN
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
       	ParticlesType knots,augment;
       	ParticlesType test;

        const double hfac = 2.0;
        double2 min(0);
        double2 max(1);
        double2 periodic(false);
        
        const int N = 1000;
        const int max_iter = 100;
        const int restart = 100;
        const double delta = std::pow(double(N),-1.0/2.0);
       	const double h = hfac*delta;
        std::cout << "using h = "<<h<<std::endl;
        typename ParticlesType::value_type p;

        std::default_random_engine generator(123);
        std::uniform_real_distribution<double> distribution(0.0,1.0);
        for (int i=0; i<N; ++i) {
            get<position>(p) = double2(distribution(generator),
                                       distribution(generator));
            knots.push_back(p);

            get<position>(p) = double2(distribution(generator),
                                       distribution(generator));
            test.push_back(p);
        }
        augment.push_back(p);

	    knots.init_neighbour_search(min,max,periodic);

        Symbol<alpha> al;
        Symbol<interpolated> interp;
        Symbol<position> r;
        Symbol<constant2> c2;
        Label<0,ParticlesType> a(knots);
        Label<1,ParticlesType> b(knots);
        Label<0,ParticlesType> i(augment);
        Label<1,ParticlesType> j(augment);
        Label<0,ParticlesType> k(test);
        auto dx = create_dx(a,b);
        auto dx2 = create_dx(k,b);
        Accumulate<std::plus<double> > sum;

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

        vector_type phi(N+1), gamma(N+1);
        for (int i=0; i<knots.size(); ++i) {
            const double x = get<position>(knots[i])[0];
            const double y = get<position>(knots[i])[1];
            phi[i] = funct(x,y);
        }
        phi[knots.size()] = 0;

        Eigen::ConjugateGradient<decltype(W), 
            Eigen::Lower|Eigen::Upper,  Eigen::DiagonalPreconditioner<double>> cg;
        cg.setMaxIterations(max_iter);
        cg.compute(W);
        gamma = cg.solve(phi);
        std::cout << std::endl << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;

        Eigen::BiCGSTAB<decltype(W), Eigen::DiagonalPreconditioner<double>> bicg;
        bicg.setMaxIterations(max_iter);
        bicg.compute(W);
        gamma = bicg.solve(phi);
        std::cout << "BiCGSTAB: #iterations: " << bicg.iterations() << ", estimated error: " << bicg.error() << std::endl;

        Eigen::MINRES<decltype(W), Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double>> minres;
        minres.setMaxIterations(max_iter);
        minres.compute(W);
        gamma = minres.solve(phi);
        std::cout << "MINRES:   #iterations: " << minres.iterations() << ", estimated error: " << minres.error() << std::endl;

        Eigen::GMRES<decltype(W), Eigen::DiagonalPreconditioner<double>> gmres;
        gmres.setMaxIterations(max_iter);
        gmres.set_restart(restart);
        gmres.compute(W);
        gamma = gmres.solve(phi);
        std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;

        Eigen::DGMRES<decltype(W), Eigen::DiagonalPreconditioner<double>> dgmres;
        dgmres.setMaxIterations(max_iter);
        dgmres.set_restart(restart);
        dgmres.compute(W);
        gamma = dgmres.solve(phi);
        std::cout << "DGMRES:   #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;

        phi = W*gamma;
        for (int i=0; i<knots.size(); ++i) {
            const double x = get<position>(knots[i])[0];
            const double y = get<position>(knots[i])[1];
            phi[i] = funct(x,y);
            TS_ASSERT_DELTA(phi[i],funct(x,y),2e-3); 
        }
        TS_ASSERT_DELTA(phi[knots.size()],0,2e-3); 


        map_type alpha_wrap(get<alpha>(knots).data(),N);
        map_type interp_wrap(get<interpolated>(knots).data(),N);

        alpha_wrap = gamma.segment<N>(0);
        vector_type beta = vector_type::Constant(N,gamma[N]);
        interp_wrap = G*alpha_wrap + beta;

        double rms_error = 0;
        double scale = 0;
        for (int i=0; i<knots.size(); ++i) {
            const double x = get<position>(knots[i])[0];
            const double y = get<position>(knots[i])[1];
            const double truth = funct(x,y);
            const double eval_value = get<interpolated>(knots[i]);
            rms_error += std::pow(eval_value-truth,2);
            scale += std::pow(truth,2);
            TS_ASSERT_DELTA(eval_value,truth,2e-3); 
        }
        std::cout << "rms_error for compact support, at centers  = "<<std::sqrt(rms_error/scale)<<std::endl;
        TS_ASSERT_LESS_THAN(std::sqrt(rms_error/scale),1e-4);

        map_type interp_wrap_test(get<interpolated>(test).data(),N);
        auto G_test = create_sparse_operator(test,knots,2*h,kernel);
        interp_wrap_test = G_test*alpha_wrap + beta;
                
        rms_error = 0;
        scale = 0;
        for (int i=0; i<test.size(); ++i) {
            const double x = get<position>(test[i])[0];
            const double y = get<position>(test[i])[1];
            const double truth = funct(x,y);
            const double eval_value = get<interpolated>(test[i]);
            rms_error += std::pow(eval_value-truth,2);
            scale += std::pow(truth,2);
            TS_ASSERT_DELTA(eval_value,truth,1e-1/truth); 
        }
        std::cout << "rms_error for compact support, away from centers  = "<<std::sqrt(rms_error/scale)<<std::endl;
        TS_ASSERT_LESS_THAN(std::sqrt(rms_error/scale),4e-2);

        
#endif // HAVE_EIGEN
    }

    void test_bucket_search_parallel() {
        helper_global<bucket_search_parallel>();
        helper_compact<bucket_search_parallel>();
    }

    void test_bucket_search_serial() {
        helper_global<bucket_search_serial>();
        helper_compact<bucket_search_serial>();
    }



};

#endif /* RBF_INTERPOLATION_TEST_H_ */
//->
