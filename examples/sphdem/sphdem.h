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

enum {DEM_FORCE, DEM_VELOCITY};
typedef std::tuple<Vect3d,Vect3d> DemTuple;
typedef Particles<DemTuple> DemType;

enum {PARAMS_SPH_DT,PARAMS_DEM_DT,PARAMS_DEM_MASS,PARAMS_DEM_DIAMETER,PARAMS_DEM_K,PARAMS_DEM_GAMMA};
typedef std::tuple<double,double,double,double,double> ParamTuple;

template<typename GeometryType>
void dem_start(ptr<DemType> dem,
		ptr<ParamTuple> params,
		ptr<GeometryType> geometry) {

	const double dt = params->get<PARAMS_DEM_DT>();
	const double dem_diameter = params->get<PARAMS_DEM_DIAMETER>();
	const double dem_k = params->get<PARAMS_DEM_K>();
	const double dem_gamma = params->get<PARAMS_DEM_GAMMA>();
	const double dem_mass = params->get<PARAMS_DEM_MASS>();


	dem.update_positions(dem.begin(),dem.end(),[](DemType::Value& i) {
		const Vect3d& r = i.get_position();
		Vect3d& f = i.get_data().get<DEM_FORCE>();
		Vect3d& v = i.get_data().get<DEM_VELOCITY>();
		Vect3d& v0 = i.get_data().get<DEM_VELOCITY0>();


		v0 = v + dt/2*f;
		v += dt * f;
		return r + dt * v;
	});

	std::for_each(dem.begin(),dem.end(),[params](DemType::Value& i) {
		Vect3d& r = i.get_position();
		Vect3d& f = i.get_data().get<DEM_FORCE>();
		Vect3d& v = i.get_data().get<DEM_VELOCITY>();

		f = 0;
		for (const DemType::Value& j: i.get_in_radius(r,dem_diameter)) {
			Vect3d& rj = j.get_position();
			Vect3d& vj = j.get_data().get<DEM_VELOCITY>();

			const Vect3d dx = r-rj;
			const double r = dx.norm();
			const double overlap = dem_diameter-r;
			if (overlap>0) {
				const Vect3d normal = dx/r;
				const Vect3d dv = v-vj;
				const double overlap_dot = -dot(dv,normal);
				f += (dem_k*overlap + dem_gamma*overlap_dot)*normal/dem_mass;
			}
		}

	});

}

template<typename GeometryType>
void dem_end(ptr<DemType> dem,
		ptr<ParamTuple> params,
		ptr<GeometryType> geometry) {

	const double dt = params->get<PARAMS_DEM_DT>();
	const double dem_diameter = params->get<PARAMS_DEM_DIAMETER>();
	const double dem_k = params->get<PARAMS_DEM_K>();
	const double dem_gamma = params->get<PARAMS_DEM_GAMMA>();
	const double dem_mass = params->get<PARAMS_DEM_MASS>();

	std::for_each(dem.begin(),dem.end(),[](DemType::Value& i) {
		Vect3d& f = i.get_data().get<DEM_FORCE>();
		Vect3d& v = i.get_data().get<DEM_VELOCITY>();
		Vect3d& v0 = i.get_data().get<DEM_VELOCITY0>();

		v = v0 + dt/2 * f;
	});
}

