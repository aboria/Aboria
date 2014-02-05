/*
 * sphdem.cpp
 *
 *  Created on: 30 Jan 2014
 *      Author: mrobins
 */

#include "sphdem.h"
#include "Visualisation.h"
#include <chrono>
#include <random>

int main(int argc, char **argv) {

	auto sph = SphType::New();
	auto dem = DemType::New();
	auto params = ptr<Params>(new Params());

	const double pi = 3.14;
	const int timesteps = 100000;
	const int nout = 1000;
	const int timesteps_per_out = timesteps/nout;
	const double L = 31.0/1000.0;
	const int ndem = 100;
	params->dem_diameter = 0.0011;
	params->dem_gamma = 0.0004;
	params->dem_k = 1.0e01;
	const double dem_vol = (1.0/6.0)*pi*pow(params->dem_diameter,3);
	const double dem_dens = 1160.0;
	params->dem_mass = dem_vol*dem_dens;
	const double dem_min_reduced_mass = 0.5*params->dem_mass;
	params->dem_dt = (1.0/50.0)*PI/sqrt(params->dem_k/dem_min_reduced_mass-pow(0.5*params->dem_gamma/dem_min_reduced_mass,2));


	auto geometry = [params](DemType::Value& i) {
		Vect3d acceleration;
		acceleration << 0,0,-9.8;
		const Vect3d& r = i.get_position();
		Vect3d& v = std::get<DEM_VELOCITY>(i.get_data());
		const double dem_diameter = params->dem_diameter;
		const double dem_k = params->dem_k;
		const double dem_gamma = params->dem_gamma;
		const double dem_mass = params->dem_mass;

		const double overlap = dem_diameter/2.0-r[2];
		if (overlap>0) {
			const double overlap_dot = -v[2];
			const Vect3d normal(0,0,1);
			const Vect3d gravity(0,0,-9.8);
			acceleration += (dem_k*overlap + dem_gamma*overlap_dot)*normal/dem_mass;
		}
		return acceleration;
	};

	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0,1.0);
	auto dice = std::bind ( distribution, generator );
	dem->create_particles(ndem,[ndem,L,&dice](DemType::Value& i) {
		Vect3d& v = std::get<DEM_VELOCITY>(i.get_data());
		Vect3d& f = std::get<DEM_FORCE>(i.get_data());

		v << 0,0,0;
		f << 0,0,0;
		Vect3d position(dice()*L,dice()*L,dice()*L);
		std::cout << "creating particle at "<<position<<std::endl;
		return position;
	});

	const Vect3d min(-2*params->dem_diameter,-2*params->dem_diameter,-2*params->dem_diameter);
	const Vect3d max(L+2*params->dem_diameter,L+2*params->dem_diameter,L+2*params->dem_diameter);
	const Vect3b periodic(true,true,false);
	Visualisation vis(min,max);
	auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	dem->copy_to_vtk_grid(grid);
	vis.glyph_points(grid);
	vis.start_render_loop();

	dem->init_neighbour_search(min,max,params->dem_diameter,periodic);
	for (int i = 0; i < nout; ++i) {
		for (int k = 0; k < timesteps_per_out; ++k) {
			dem_start(dem,params,geometry);
			dem_end(dem,params,geometry);
		}
		std::cout <<"iteration "<<i<<std::endl;

		vis.stop_render_loop();
		dem->copy_to_vtk_grid(grid);
		vis.restart_render_loop();
	}


}
