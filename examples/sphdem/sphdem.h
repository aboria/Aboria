/*
 * sphdem.h
 *
 *  Created on: 30 Jan 2014
 *      Author: mrobins
 */

#include "Aboria.h"
#include <tuple>
using namespace Aboria;

enum {SPH_FORCE, SPH_H, SPH_DENS, SPH_VELOCITY, SPH_VELOCITY0};
typedef std::tuple<Vect3d,double,double,Vect3d,Vect3d> SphTuple;
typedef Particles<SphTuple> SphType;

enum {DEM_FORCE, DEM_VELOCITY, DEM_VELOCITY0};
typedef std::tuple<Vect3d,Vect3d,Vect3d> DemTuple;
typedef Particles<DemTuple> DemType;

enum {PARAMS_SPH_DT,PARAMS_DEM_DT,PARAMS_DEM_MASS,PARAMS_DEM_DIAMETER,PARAMS_DEM_K,PARAMS_DEM_GAMMA};
typedef std::tuple<double,double,double,double,double,double> ParamTuple;

template<typename GeometryType>
void dem_start(ptr<DemType> dem,
		ptr<ParamTuple> params,
		GeometryType geometry) {

	const double dt = std::get<PARAMS_DEM_DT>(*params);
	const double dem_diameter = std::get<PARAMS_DEM_DIAMETER>(*params);
	const double dem_k = std::get<PARAMS_DEM_K>(*params);
	const double dem_gamma = std::get<PARAMS_DEM_GAMMA>(*params);
	const double dem_mass = std::get<PARAMS_DEM_MASS>(*params);


	dem->update_positions(dem->begin(),dem->end(),[dt](DemType::Value& i) {
		const Vect3d& r = i.get_position();
		Vect3d& f = std::get<DEM_FORCE>(i.get_data());
		Vect3d& v = std::get<DEM_VELOCITY>(i.get_data());
		Vect3d& v0 = std::get<DEM_VELOCITY0>(i.get_data());


		v0 = v + dt/2*f;
		v += dt * f;
		return r + dt * v;
	});

	std::for_each(dem->begin(),dem->end(),[&geometry,dem,dem_k,dem_gamma,dem_mass,dem_diameter](DemType::Value& i) {
		const Vect3d& r = i.get_position();
		Vect3d& f = std::get<DEM_FORCE>(i.get_data());
		Vect3d& v = std::get<DEM_VELOCITY>(i.get_data());

		f << 0,0,0;
		f = f + geometry(i);

		for (const DemType::Value& j: i.get_in_radius(dem,dem_diameter)) {
			const Vect3d& rj = j.get_position();
			const Vect3d& vj = std::get<DEM_VELOCITY>(j.get_data());
			if (i.get_id()==j.get_id()) return;

			const Vect3d dx = r-rj;
			const double r = dx.norm();
			const double overlap = dem_diameter-r;
			if (overlap>0) {
				const Vect3d normal = dx/r;
				const Vect3d dv = v-vj;
				const double overlap_dot = -dv.dot(normal);
				f += (dem_k*overlap + dem_gamma*overlap_dot)*normal/dem_mass;
			}
		}

	});

}

template<typename GeometryType>
void dem_end(ptr<DemType> dem,
		ptr<ParamTuple> params,
		GeometryType geometry) {

	const double dt = std::get<PARAMS_DEM_DT>(*params);
	const double dem_diameter = std::get<PARAMS_DEM_DIAMETER>(*params);
	const double dem_k = std::get<PARAMS_DEM_K>(*params);
	const double dem_gamma = std::get<PARAMS_DEM_GAMMA>(*params);
	const double dem_mass = std::get<PARAMS_DEM_MASS>(*params);

	std::for_each(dem->begin(),dem->end(),[dt](DemType::Value& i) {
		Vect3d& f = std::get<DEM_FORCE>(i.get_data());
		Vect3d& v = std::get<DEM_VELOCITY>(i.get_data());
		Vect3d& v0 = std::get<DEM_VELOCITY0>(i.get_data());

		v = v0 + dt/2 * f;
	});
}

