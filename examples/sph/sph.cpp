/*
 * sphdem.cpp
 *
 *  Created on: 30 Jan 2014
 *      Author: mrobins
 */

#include "sph.h"
#include "Visualisation.h"
#include <chrono>
#include <random>


int main(int argc, char **argv) {

	auto sph = SphType::New();
	auto params = ptr<Params>(new Params());

	const int timesteps = 1000;
	const int nout = 1000;
	const int timesteps_per_out = timesteps/nout;
	const double L = 31.0/1000.0;
	const int nx = 10;

	/*
	 * sph parameters
	 */
	params->sph_hfac = 1.5;
	params->sph_visc = 8.9e-07;
	params->sph_refd = 1000.0;
	params->sph_dens = 1000.0;
	params->sph_gamma = 7;
	const double VMAX = 2.0*sqrt(2*9.81*L);
	const double CSFAC = 10.0;
	params->sph_spsound = CSFAC*VMAX;
	params->sph_prb = pow(params->sph_refd/params->sph_dens,params->sph_gamma-1.0)*pow(params->sph_spsound,2)*params->sph_refd/params->sph_gamma;
	const double psep = L/nx;
	params->sph_dt = std::min(0.25*params->sph_hfac*psep/params->sph_spsound,0.125*pow(params->sph_hfac*psep,2)/params->sph_visc);
	params->sph_mass = params->sph_dens*pow(psep,NDIM);

	std::cout << "h = "<<params->sph_hfac*psep<<" vmax = "<<VMAX<<std::endl;

	params->sph_time_damping = params->sph_dt*500;
	params->time = 0;



	auto sph_geometry = [](SphType::Value& i) {
		Vect3d acceleration(0,0,-9.8);
		return acceleration;
	};

	const Vect3d min(0,0,-3.0*psep);
	const Vect3d max(L,L,L);
	const Vect3b periodic(true,true,false);


	/*
	 * create sph and dem particles
	 */
	sph->create_particles_grid(min,max,Vect3i(nx,nx,nx+3),[psep,params](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);
		h = params->sph_hfac*psep;
		omega = 1.0;
		v  = Vect3d(0,0,0);
		v0 = Vect3d(0,0,0);
		dddt = 0;
		rho = params->sph_dens;
		f = Vect3d(0,0,0);
		if (r[2]<0) {
			fixed = true;
		} else {
			fixed = false;
		}
	});

	/*
	 * setup output stuff
	 */
	auto sph_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	sph->copy_to_vtk_grid(sph_grid);
	Visualisation::vtkWriteGrid("at_start_sph",0,sph_grid);


	std::cout << "starting...."<<std::endl;
	sph->init_neighbour_search(min,max,2*params->sph_hfac*psep,periodic);

	for (int i = 0; i < nout; ++i) {
		for (int k = 0; k < timesteps_per_out; ++k) {
			//std::this_thread::sleep_for(std::chrono::seconds(1));
			sph_timestep(sph,params,sph_geometry);
		}
		std::cout <<"iteration "<<i<<std::endl;

		sph->copy_to_vtk_grid(sph_grid);
		Visualisation::vtkWriteGrid("sph",i,sph_grid);
	}


}
