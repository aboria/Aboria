/*
 * sphdem.cpp
 *
 *  Created on: 30 Jan 2014
 *      Author: mrobins
 */

#include "sphdem.h"

int main(int argc, char **argv) {

	auto sph = SphType::New();
	auto dem = DemType::New();
	auto params = ptr<ParamTuple>(new ParamTuple());

	const double pi = 3.14;
	const int n = 1000;
	const double L = 31.0/1000.0;
	const int ndem = 2;
	const double dem_diameter = 0.0011;
	const double dem_gamma = 0.0004;
	const double dem_k = 1.0e01;
	const double dem_vol = (1.0/6.0)*pi*pow(dem_diameter,3);
	const double dem_dens = 1160.0;
	const double dem_mass = dem_vol*dem_dens;
	const double dem_min_reduced_mass = 0.5*dem_mass;
	const double dem_dt = (1.0/50.0)*PI/sqrt(dem_k/dem_min_reduced_mass-pow(0.5*dem_gamma/dem_min_reduced_mass,2));

	std::get<PARAMS_DEM_DT>(*params) = dem_dt;
	std::get<PARAMS_DEM_DIAMETER>(*params) = dem_diameter;
	std::get<PARAMS_DEM_GAMMA>(*params) = dem_gamma;
	std::get<PARAMS_DEM_K>(*params) = dem_k;
	std::get<PARAMS_DEM_MASS>(*params) = dem_mass;


	auto geometry = [params](DemType::Value& i) {
		Vect3d acceleration;
		acceleration << 0,0,-9.8;
		const Vect3d& r = i.get_position();
		Vect3d& v = std::get<DEM_VELOCITY>(i.get_data());
		const double dem_diameter = std::get<PARAMS_DEM_DIAMETER>(*params);
		const double dem_k = std::get<PARAMS_DEM_K>(*params);
		const double dem_gamma = std::get<PARAMS_DEM_GAMMA>(*params);
		const double dem_mass = std::get<PARAMS_DEM_MASS>(*params);

		const double overlap = dem_diameter/2.0-r[2];
		if (overlap>0) {
			const double overlap_dot = -v[2];
			const Vect3d normal(0,0,1);
			const Vect3d gravity(0,0,-9.8);
			acceleration += (dem_k*overlap + dem_gamma*overlap_dot)*normal/dem_mass;
		}
		return acceleration;
	};

	int c = 0;
	dem->create_particles(ndem,[ndem,L,&c](DemType::Value& i) {
		Vect3d& v = std::get<DEM_VELOCITY>(i.get_data());
		Vect3d& f = std::get<DEM_FORCE>(i.get_data());

		v << 0,0,0;
		f << 0,0,0;
		Vect3d position(c*L/ndem,0,L);
		std::cout << "creating particle at "<<position<<std::endl;
		c++;
		return position;
	});

	const Vect3d min(-2*dem_diameter,-2*dem_diameter,-2*dem_diameter);
	const Vect3d max(L+2*dem_diameter,L+2*dem_diameter,L+2*dem_diameter);
	dem->init_neighbour_search(min,max,dem_diameter);
	std::cout << "params are "<<*params<<std::endl;
	for (int i = 0; i < n; ++i) {
		std::cout <<"iteration "<<i<<std::endl;
		dem_start(dem,params,geometry);
		dem_end(dem,params,geometry);
	}


}
