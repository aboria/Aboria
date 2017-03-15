/*
 * sph.cpp
 *
 *  Created on: 30 Jan 2014
 *      Author: mrobins
 */
#ifndef SPH_TEST_H_
#define SPH_TEST_H_

#include <cxxtest/TestSuite.h>

//[sph
#include "Aboria.h"
using namespace Aboria;
#include <boost/math/constants/constants.hpp>

const double PI = boost::math::constants::pi<double>();
const double WCON_WENDLAND = 21.0/(256.0*PI);
const unsigned int NDIM = 3;

struct F_fun {
    typedef double result_type;
    double operator()(const double r,const double h) const {
        if (r==0) return 0;
        const double q = r/h;
        if (q<=2.0) {
            return (1/std::pow(h,NDIM+2))*WCON_WENDLAND*(-4*std::pow(2-q,3)*(1+2*q) + 2*std::pow(2-q,4))/q;
        } else {
            return 0.0;
        }
    }
};

struct W_fun {
    typedef double result_type;
    double operator()(const double r,const double h) const {
        const double q = r/h;
        if (q<=2.0) {
            return (1/std::pow(h,NDIM))*WCON_WENDLAND*std::pow(2.0-q,4)*(1.0+2.0*q);
        } else {
            return 0.0;
        }
    }
};

ABORIA_BINARY_FUNCTION(F, F_fun, SymbolicDomain);
ABORIA_BINARY_FUNCTION(W, W_fun, SymbolicDomain);



//<-
class SPHTest : public CxxTest::TestSuite {
public:

    template<template <typename> class SearchMethod>
    void helper_sph(void) {
//->
//=int main() {
        ABORIA_VARIABLE(kernel_radius,double,"kernel radius");
        ABORIA_VARIABLE(velocity,double3,"velocity");
        ABORIA_VARIABLE(velocity_tmp,double3,"temp velocity");
        ABORIA_VARIABLE(varh_omega,double,"varh omega");
        ABORIA_VARIABLE(density,double,"density");
        ABORIA_VARIABLE(total_force,double3,"total force");
        ABORIA_VARIABLE(is_fixed,uint8_t,"fixed boundary");
        ABORIA_VARIABLE(pressure_div_density2,double,"pressure div density2");

        typedef Particles<std::tuple<kernel_radius,velocity,velocity_tmp,varh_omega,density,total_force,is_fixed,pressure_div_density2>,3,std::vector,SearchMethod> sph_type;
        typedef position_d<3> position;
        sph_type sph;

        Symbol<position> p;
        Symbol<id> id_;
        Symbol<velocity> v;
        Symbol<velocity_tmp> v0;
        Symbol<density> rho;
        Symbol<total_force> dvdt;
        Symbol<is_fixed> fixed;
        Symbol<varh_omega> omega;
        Symbol<kernel_radius> h;
        Symbol<pressure_div_density2> pdr2;
        Label<0,sph_type> a(sph);
        Label<1,sph_type> b(sph);
     
        const int timesteps = 500;
        const int nout = 10;
        const int timesteps_per_out = timesteps/nout;
        const double L = 31.0/1000.0;
        const int nx = 5;

        /*
         * sph parameters
         */
        const double hfac = 1.5;
        const double visc = 8.9e-07;
        const double refd = 1000.0;
        const double dens = 1000.0;
        const double gamma = 7;
        const double VMAX = 2.0*sqrt(2*9.81*L);
        const double CSFAC = 10.0;
        const double spsound = CSFAC*VMAX;
        const double prb = std::pow(refd/dens,gamma-1.0)*std::pow(spsound,2)*refd/gamma;
        const double psep = L/nx;
        double dt = std::min(0.25*hfac*psep/spsound,0.125*std::pow(hfac*psep,2)/visc);
        const double mass = dens*std::pow(psep,NDIM);

        std::cout << "h = "<<hfac*psep<<" vmax = "<<VMAX<<std::endl;

        const double time_damping = dt*500;
        double t = 0;

        const double3 low(0,0,-3.0*psep);
        const double3 high(L,L,L);
        const bool3 periodic(true,true,false);

        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < nx; j++) {
                for (int k = 0; k < nx+3; k++) {
                    typename sph_type::value_type p;
                    get<position>(p) = low + double3((i+0.5)*psep,(j+0.5)*psep,(k+0.5)*psep);
                    get<kernel_radius>(p) = hfac*psep;
                    get<velocity>(p) = double3(0,0,0);
                    get<velocity_tmp>(p) = double3(0,0,0);
                    get<varh_omega>(p) = 1.0;
                    get<density>(p) = dens;
                    get<total_force>(p) = double3(0,0,0);
                    get<is_fixed>(p) = get<position>(p)[2]<0;
                    sph.push_back(p);
                }
            }
        }

        std::cout << "starting...."<<std::endl;
        sph.init_neighbour_search(low,high,periodic);


        auto dx = create_dx(a,b);
        AccumulateWithinDistance<std::plus<double3> > sum_vect(2*hfac*psep);
        AccumulateWithinDistance<std::plus<double> > sum(2*hfac*psep);
        Accumulate<Aboria::max<double> > max;
        max.set_init(0);
        Accumulate<Aboria::min<double> > min;
        min.set_init(L);


        for (int i = 0; i < nout; ++i) {
#ifdef HAVE_VTK
            vtkWriteGrid("sph",i,sph.get_grid(true));
#endif
            for (int k = 0; k < timesteps_per_out; ++k) {
                /*
                 * 0 -> 1/2 step
                 */
                v[a] += if_else(fixed[a]==false,dt/2 * dvdt[a],0);
                if (t < time_damping) v[a] *= 0.98;
                p[a] += dt/2 * v[a];

                /*
                 * Calculate omega
                 */
                omega[a] = 1.0 - (mass/(rho[a]*NDIM))
                        *sum(b,if_else(norm(dx)<2*h[a]
                                    ,pow(norm(dx),2)*F(norm(dx),h[a])+NDIM*W(norm(dx),h[a])
                                    ,0));

                /*
                 * 1/2 -> 1 step
                 */

                /* 
                 * calculate change in density and kernel radius
                 */
                rho[a] += dt*(mass/omega[a]) 
                            * sum(b,if_else(norm(dx)<2*h[a]
                                        ,dot(v[b]-v[a],dx*F(norm(dx),h[a]))
                                        ,0));
                h[a] = pow(mass/rho[a],1.0/NDIM);
                
                /* 
                 * reset neighbour search for new kernel radius
                 */
                const double maxh = eval(max(a,h[a]));
                const double minh = eval(min(a,h[a]));
                sum_vect.set_max_distance(2*maxh);
                sum.set_max_distance(2*maxh);

                /* 
                 * advance velocity, position and calculate pdr2
                 */
                v0[a] = v[a];
                v[a] += if_else(fixed[a]==false,dt/2 * dvdt[a],0);
                pdr2[a] = prb*(pow(rho[a]/refd,gamma) - 1.0)/pow(rho[a],2);
                p[a] += dt/2 * v0[a];


                /*
                 * acceleration due to gravity
                 */
                dvdt[a] = double3(0,0,-9.81);

                /*
                 * acceleration on SPH calculation
                 */
                dvdt[a] += mass*sum_vect(b,if_else(norm(dx)<2*h[a]
                        
                        // pressure force 
                        ,((1.0/omega[a])*pdr2[a]*F(norm(dx),h[a]) + (1.0/omega[b])*pdr2[b]*F(norm(dx),h[b]))*dx 

                        // viscous force 
                        + (v[a]-v[b])*visc*(rho[a]+rho[b])/(rho[a]*rho[b])*F(norm(dx),h[a])

                        ,0)
                        );

                
                /*
                 * 1/2 -> 1 step for velocity
                 */
                v[a] = if_else(fixed[a]==false,v0[a] + dt/2 * dvdt[a],0);

                t += dt;
                    
                /*
                 * sph timestep condition
                 */
                dt = std::min(0.25*minh/spsound,0.125*std::pow(minh,2)/visc);

            }
            std::cout <<"iteration "<<i<<std::endl;
        }
    }
//]


    void test_bucket_search_parallel() {
        helper_sph<bucket_search_parallel>();
    }

    void test_bucket_search_serial() {
        helper_sph<bucket_search_serial>();
    }

    void test_nanoflann_adaptor() {
        helper_sph<nanoflann_adaptor>();
    }



};

#endif /* SPH_TEST_H_ */

