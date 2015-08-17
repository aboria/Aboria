/*
 * sph.cpp
 *
 *  Created on: 30 Jan 2014
 *      Author: mrobins
 */

#include "Aboria.h"
using namespace Aboria;

#include "Visualisation.h"


const double WCON_WENDLAND = 21.0/(256.0*PI);

struct F_fun {
    double operator()(const double r,const double h) {
        const double q = r/h;
        if (q<=2.0) {
            return (1/std::pow(h,NDIM+2))*WCON_WENDLAND*(-4*std::pow(2-q,3)*(1+2*q) + 2*std::pow(2-q,4))/q;
        } else {
            return 0.0;
        }
    }
};

struct W_fun {
    double operator()(const double r,const double h) {
        const double q = r/h;
        if (q<=2.0) {
            return (1/std::pow(h,NDIM))*WCON_WENDLAND*std::pow(2.0-q,4)*(1.0+2.0*q);
        } else {
            return 0.0;
        }
    }
};

ABORIA_BINARY_FUNCTION(F, F_fun, DataVectorDomain);
ABORIA_BINARY_FUNCTION(W, W_fun, DataVectorDomain);

int main(int argc, char **argv) {
    ABORIA_VARIABLE(kernel_radius,double,"kernel radius");
    ABORIA_VARIABLE(velocity,Vect3d,"velocity");
    ABORIA_VARIABLE(velocity_tmp,Vect3d,"temp velocity");
    ABORIA_VARIABLE(varh_omega,double,"varh omega");
    ABORIA_VARIABLE(density,double,"density");
    ABORIA_VARIABLE(total_force,Vect3d,"total force");
    ABORIA_VARIABLE(is_fixed,bool,"fixed boundary");
    ABORIA_VARIABLE(pressure_div_density2,double,"pressure div density2");

    typedef Particles<kernel_radius,velocity,velocity_tmp,varh_omega,density,total_force,is_fixed,pressure_div_density2> sph_type;
    sph_type sph;

    auto p = get_vector<position>(sph);
    auto v = get_vector<velocity>(sph);
    auto v0 = get_vector<velocity_tmp>(sph);
    auto rho = get_vector<density>(sph);
    auto dvdt = get_vector<total_force>(sph);
    auto fixed = get_vector<is_fixed>(sph);
    auto omega = get_vector<varh_omega>(sph);
    auto h = get_vector<kernel_radius>(sph);
    auto pdr2 = get_vector<pressure_div_density2>(sph);

	const int timesteps = 1000;
	const int nout = 1000;
	const int timesteps_per_out = timesteps/nout;
	const double L = 31.0/1000.0;
	const int nx = 10;

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
	const double dt = std::min(0.25*hfac*psep/spsound,0.125*std::pow(hfac*psep,2)/visc);
	const double mass = dens*std::pow(psep,NDIM);

	std::cout << "h = "<<hfac*psep<<" vmax = "<<VMAX<<std::endl;

	const double time_damping = dt*500;
	const double t = 0;

	const Vect3d low(0,0,-3.0*psep);
	const Vect3d high(L,L,L);
	const Vect3b periodic(true,true,false);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < nx; j++) {
            for (int k = -3; k < nx; k++) {
                sph_type::value_type p(low + Vect3d((i+0.5)*psep,(j+0.5)*psep,k*psep));
                set<kernel_radius>(p,hfac*psep);
                set<velocity>(p,Vect3d(0,0,0));
                set<velocity_tmp>(p,Vect3d(0,0,0));
                set<varh_omega>(p,1.0);
                set<density>(p,dens);
                set<total_force>(p,Vect3d(0,0,0));
                set<is_fixed>(p,k<0);
                sph.push_back(p);
            }
        }
    }

	/*
	 * setup output stuff
	 */
	auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	sph.copy_to_vtk_grid(grid);
	Visualisation::vtkWriteGrid("at_start_sph",0,grid);

	std::cout << "starting...."<<std::endl;
	sph.init_neighbour_search(low,high,2*hfac*psep,periodic);

    Label<0> a;
    Label<1> b;
    Dx dx;
    Accumulate<std::plus<Vect3d> > sum_vect;
    Accumulate<std::plus<double> > sum;
    Accumulate<Aboria::max<double> > max;
    max.set_init(0);
    Accumulate<Aboria::min<double> > min;
    min.set_init(L);


	for (int i = 0; i < nout; ++i) {
		for (int k = 0; k < timesteps_per_out; ++k) {
	        /*
	         * 0 -> 1/2 step
	         */
            v += dt/2 + dvdt;
            if (t < time_damping) v *= 0.98;
            p += dt/2 + v;

	        /*
	         * Calculate omega
	         */
            omega = 1.0 - (mass/rho*NDIM)*sum(b=sph,pow(norm(dx),2)*F(norm(dx),h)+NDIM*W(norm(dx),h));

	        /*
	         * 1/2 -> 1 step
	         */
            rho += dt*(mass/omega) * sum(b=sph,dot(v[a]-v[b],dx*F(norm(dx),h)));
            h = pow(mass/rho,1.0/NDIM);
            
            const double maxh = eval(max(a=sph,true,h));
            const double minh = eval(min(a=sph,true,h));

            v0 = v;
            v += if_else(fixed==false,dt/2 * dvdt,0);
            pdr2 = prb*(pow(rho/refd,gamma) - 1.0);
            p += dt/2 * v0;

	        /*
	         * acceleration on SPH calculation (pressure and viscosity force)
	         */
            dvdt += mass*sum_vect(b=sph,norm(dx)<2*h,
                    
                    ((1.0/omega[a])*pdr2[a]*F(norm(dx),h[a]) + (1.0/omega[b])*pdr2[b]*F(norm(dx),h[b]))*dx 
                    + (v[a]-v[b])*visc*(rho[a]+rho[b])/(rho[a]*rho[b])*F(norm(dx),h[a])
                    );
            
	        /*
	         * 1/2 -> 1 step for velocity
	         */
            v = if_else(fixed==false,v0 + dt/2 * dvdt,0);

            t += dt;
                
	        /*
	         * sph timestep condition
	         */
	        dt = std::min(0.25*minh/spsound,0.125*std::pow(minh,2)/visc);

		}
		std::cout <<"iteration "<<i<<std::endl;

		sph.copy_to_vtk_grid(grid);
		Visualisation::vtkWriteGrid("sph",i,grid);
	}


}
