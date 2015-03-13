/*
 * sphdem.cpp
 *
 *  Created on: 30 Jan 2014
 *      Author: mrobins
 */

#include "dem.h"
#include "Visualisation.h"
#include <chrono>
#include <random>

int calc_part_num(double diam, double L){
	int num = int(L/diam);
	return num;
}

int main(int argc, char **argv) {

	auto dem = DemType::New();
	auto params = ptr<Params>(new Params());

	const int timesteps = 10000;
	const int ndem = 100;
	const int nout = 100;
	const int timesteps_per_out = timesteps/nout;
	const double L = 31.0/1000.0;
	const Vect3d min(-L/4,-L/4,0);
	const Vect3d max(L/4,L/4,L*2);
	const Vect3b periodic(true,true,false);


	params->dem_diameter = 0.0011;
	params->dem_gamma = 0.0004;
	params->dem_k = 1.0e01;
	const double dem_vol = (1.0/6.0)*PI*pow(params->dem_diameter,3);
	const double dem_dens = 1160.0;
	params->dem_mass = dem_vol*dem_dens;
	const double dem_min_reduced_mass = 0.5*params->dem_mass;
	params->dem_dt = (1.0/50.0)*PI/sqrt(params->dem_k/dem_min_reduced_mass-pow(0.5*params->dem_gamma/dem_min_reduced_mass,2));


	auto geometry = [params](DemType::Value& i) {
		Vect3d acceleration(0,0,-9.8);
		const Vect3d& r = i.get_position();
		Vect3d& v = i.get_data_elem<DEM_VELOCITY>();

		const double dem_diameter = params->dem_diameter;
		const double dem_k = params->dem_k;
		const double dem_gamma = params->dem_gamma;
		const double dem_mass = params->dem_mass;

		const double overlap = dem_diameter/2.0-r[2];
		if (overlap>0) {
			const double overlap_dot = -v[2];
			const Vect3d normal(0,0,1);
			acceleration += (dem_k*overlap + dem_gamma*overlap_dot)*normal/dem_mass;
		}
		return acceleration;
	};

	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(-L/4,L/4);

	auto dice = std::bind ( distribution, generator );

	dem->init_neighbour_search(min,max,params->dem_diameter,periodic);
	int i = 0;
	while (i < ndem) {
		const Vect3d position(dice(),dice(),2*dice()+L/2+params->dem_diameter/2);
		auto neigh = dem->get_neighbours(position);
		if (neigh.begin() == neigh.end()) {
			DemType::value_type particle(position);
			particle.get_data_elem<DEM_FORCE>() = Vect3d(0,0,0);
			particle.get_data_elem<DEM_VELOCITY>() = Vect3d(0,0,0);
			particle.get_data_elem<DEM_VELOCITY0>() = Vect3d(0,0,0);
			dem->push_back(particle);
			i++;
		}
	}
	

	auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	dem->copy_to_vtk_grid(grid);

	std::cout <<"starting.... "<<std::endl;

	for (int i = 0; i < nout; ++i) {
		for (int k = 0; k < timesteps_per_out; ++k) {

			dem_start(dem,params,geometry);
			dem_end(dem,params,geometry);
		}
		std::cout <<"iteration "<<i<<std::endl;

		dem->copy_to_vtk_grid(grid);
		Visualisation::vtkWriteGrid("dem",i,grid);
	}


}
