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

int calc_part_num(double diam, double L){
	int num = int(L/diam);
	return num;
}

int main(int argc, char **argv) {

	auto dem = DemType::New();
	auto params = ptr<Params>(new Params());

	int timesteps;
	cin>> timesteps;
	const int nout = 1000;
	const int timesteps_per_out = timesteps/nout;
	const double L = 31.0/1000.0;
	const double Lx=L;
	const double Ly=0.0024;
	const double Lz=0.0012;
//	const int ndem = 1;
/*	const int ndemx = 7;
	const int ndemy = 7;
	const int ndemz = 7;
	const int ndem=ndemx*ndemy*ndemz;
*/



	params->dem_diameter = 0.0011;
	params->dem_gamma = 0.0004;
	params->dem_k = 1.0e01;
	const double dem_vol = (1.0/6.0)*PI*pow(params->dem_diameter,3);
	const double dem_dens = 1160.0;
	params->dem_mass = dem_vol*dem_dens;
	const double dem_min_reduced_mass = 0.5*params->dem_mass;
	params->dem_dt = (1.0/50.0)*PI/sqrt(params->dem_k/dem_min_reduced_mass-pow(0.5*params->dem_gamma/dem_min_reduced_mass,2));


	const int ndemx = calc_part_num(params->dem_diameter,Lx);
	const int ndemy = calc_part_num(params->dem_diameter,Ly);
	const int ndemz = calc_part_num(params->dem_diameter,Lz);
	const int ndem=ndemx*ndemy*ndemz;



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
/*
//	std::default_random_engine generator;
//	std::uniform_real_distribution<double> distribution(0.0,1.0);

	std::uniform_real_distribution<double> distribution(0.0, 1.0);
	std::random_device rd;
	std::default_random_engine generator( rd() );

	auto dice = std::bind ( distribution, generator );

	

	dem->create_particles(ndem,[params,ndem,L,&dice](DemType::Value& i) {
		Vect3d& v = std::get<DEM_VELOCITY>(i.get_data());
		Vect3d& f = std::get<DEM_FORCE>(i.get_data());

		v << 0,2,1;
		f << 0,0,0;
		
		Vect3d position(dice()*L-params->dem_diameter/3,dice()*L-params->dem_diameter/3,L/2);
		return position;
	});

		dem->create_particles(ndem,[params,ndem,L,&dice](DemType::Value& i) {
		Vect3d& v = std::get<DEM_VELOCITY>(dem->i.get_data());
		Vect3d& f = std::get<DEM_FORCE>(dem->i.get_data());

		v << 0,0,0;
		f << 0,0,0;
		Vect3d position(L/2-params->dem_diameter/3,L/2-params->dem_diameter/3,L/2);
		return position;
	});
*/

	const Vect3d min(0,0,0);
	const Vect3d max(L,L,L);
	const Vect3b periodic(true,true,false);

	dem->create_particles_grid(min, max, Vect3i (ndemx,ndemy,ndemz),[](DemType::Value& i) {
		Vect3d& v = std::get<DEM_VELOCITY>(i.get_data());
		Vect3d& f = std::get<DEM_FORCE>(i.get_data());

		v << 0,0,0.4;
		f << 0,0,0;
	});
	

	auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	dem->copy_to_vtk_grid(grid);
	dem->init_neighbour_search(min,max,params->dem_diameter,periodic);

	//Visualisation vis(min,max);
	//vis.glyph_points(grid);
	//vis.start_render_loop();

	for (int i = 0; i < nout; ++i) {
		for (int k = 0; k < timesteps_per_out; ++k) {
			//std::this_thread::sleep_for(std::chrono::seconds(1));
			dem_start(dem,params,geometry);
			dem_end(dem,params,geometry);
		}
		//std::cout <<"iteration "<<i<<std::endl;

		//vis.stop_render_loop();
		dem->copy_to_vtk_grid(grid);
		Visualisation::vtkWriteGrid("vis/dem",i,grid);
		//vis.restart_render_loop();
	}


}
