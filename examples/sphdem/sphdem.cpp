/*
 * sphdem.cpp
 * 
 * Copyright 2014 Martin Robinson
 *
 * This file is part of Aboria.
 *
 * Aboria is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Aboria is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Aboria.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 7 Feb 2014
 *      Author: robinsonm
 */

#include "sphdem.h"
#include "Visualisation.h"

#include <vtkFloatArray.h>

int calc_part_num(double diam, double L){
	int num = int(L/diam);
	return num;
}

int main(int argc, char **argv) {
	auto dem = DemType::New();
	auto sph = SphType::New();
	auto params = ptr<Params>(new Params());

	const int timesteps = 100;
	const int nout = 100;
	const int timesteps_per_out = timesteps/nout;
	const double L = 31.0/1000.0;
	const int nx = 10;


	 /* dem parameters
	 */
	params->dem_diameter = 0.0011;
	params->dem_gamma = 0.0004;
	params->dem_k = 1.0e01;
	const double dem_vol = (1.0/6.0)*PI*pow(params->dem_diameter,3);
	const double dem_dens = 1160.0;
	params->dem_mass = dem_vol*dem_dens;
	const double dem_min_reduced_mass = 0.5*params->dem_mass;
	params->dem_dt = (1.0/50.0)*PI/sqrt(params->dem_k/dem_min_reduced_mass-pow(0.5*params->dem_gamma/dem_min_reduced_mass,2));

	/*

	/*
	 * sph parameters
	 */
	params->sph_hfac = 1.3;
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


	/*
	 * define domain / geometry
	 */
	auto dem_geometry = [params](DemType::Value& i) {
		Vect3d acceleration;
		acceleration << 0,0,-9.8;
		REGISTER_DEM_PARTICLE(i);
		const double dem_diameter = params->dem_diameter;
		const double dem_k = params->dem_k;
		const double dem_gamma = params->dem_gamma;
		const double dem_mass = params->dem_mass;
		const double dem_vol = params->dem_vol;

		const double overlap = dem_diameter/2.0-r[2];
		if (overlap>0) {
			const double overlap_dot = -v[2];
			const Vect3d normal(0,0,1);
			acceleration += (dem_k*overlap + dem_gamma*overlap_dot)*normal/dem_mass;
		}
		
		return acceleration;
	};

	auto sph_geometry = [](SphType::Value& i) {
		Vect3d acceleration;
		acceleration << 0,0,-9.8;
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
		v << 0,0,0;
		v0 << 0,0,0;
		dddt = 0;
		e = 1;
		rho = params->sph_dens;
		f << 0,0,0;
		if ((r[1]<2) || (r[1]>nx-2)){
			fixed = true;
		} else {
			fixed = false;
		}
		if (r[2]<0) {
			fixed = true;
		} else {
			fixed = false;
		}
	});

	dem->create_particles(1,[L](DemType::Value& i) {
		REGISTER_DEM_PARTICLE(i);
		v << 0,0,0;
		v0 << 0,0,0;
		f << 0,0,0;
		return Vect3d(L/2,L/2,0.75*L);
	});


	/*
	 * setup output stuff
	 */
	auto sph_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	auto dem_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	sph->copy_to_vtk_grid(sph_grid);
	dem->copy_to_vtk_grid(dem_grid);
	Visualisation::vtkWriteGrid("vis/at_start_sph",0,sph_grid);
	Visualisation::vtkWriteGrid("vis/at_start_dem",0,dem_grid);


	std::cout << "starting...."<<std::endl;
	sph->init_neighbour_search(min,max,2*params->sph_hfac*psep,periodic);
	dem->init_neighbour_search(min,max,params->dem_diameter,periodic);

	for (int i = 0; i < nout; ++i) {
		for (int k = 0; k < timesteps_per_out; ++k) {
			//std::this_thread::sleep_for(std::chrono::seconds(1));
			sphdem(sph,dem,params,sph_geometry,dem_geometry);
		}
		std::cout <<"iteration "<<i<<std::endl;
		
		sph->copy_to_vtk_grid(sph_grid);
		dem->copy_to_vtk_grid(dem_grid);
		Visualisation::vtkWriteGrid("vis/sph",i,sph_grid);
		Visualisation::vtkWriteGrid("vis/dem",i,dem_grid);
	}
	
	
}
