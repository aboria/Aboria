/*
 * sphdem.h
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

#ifndef SPHDEM_H_
#define SPHDEM_H_

#include "Aboria.h"
#include "sph_common.h"
#include <tuple>
using namespace Aboria;


enum {DEM_FORCE, DEM_VELOCITY, DEM_VELOCITY0};
typedef std::tuple<Vect3d,Vect3d,Vect3d> DemTuple;
typedef Particles<DemTuple> DemType;

enum {SPH_FORCE, SPH_VELOCITY, SPH_VELOCITY0, SPH_DENS, SPH_POROSITY, SPH_H, SPH_DDDT,SPH_PDR2,SPH_OMEGA};
typedef std::tuple<Vect3d,Vect3d,Vect3d,double,double,double,double,double,double> SphTuple;
typedef Particles<SphTuple> SphType;

struct Params {
	double sph_dt,dem_dt,dem_mass,dem_diameter,dem_k,dem_gamma;
};

template<typename GeometryType>
void dem_start(ptr<DemType> dem,
		ptr<Params> params,
		GeometryType geometry) {

	const double dt = params->dem_dt;
	const double dem_diameter = params->dem_diameter;
	const double dem_k = params->dem_k;
	const double dem_gamma = params->dem_gamma;
	const double dem_mass = params->dem_mass;


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

		for (auto tpl: i.get_neighbours(dem)) {
			const Vect3d& dx = std::get<1>(tpl);
			const DemType::Value& j = std::get<0>(tpl);
			const Vect3d& vj = std::get<DEM_VELOCITY>(j.get_data());
			if (i.get_id()==j.get_id()) continue;

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
		ptr<Params> params,
		GeometryType geometry) {

	const double dt = params->dem_dt;
	const double dem_diameter = params->dem_diameter;
	const double dem_k = params->dem_k;
	const double dem_gamma = params->dem_gamma;
	const double dem_mass = params->dem_mass;

	std::for_each(dem->begin(),dem->end(),[dt](DemType::Value& i) {
		Vect3d& f = std::get<DEM_FORCE>(i.get_data());
		Vect3d& v = std::get<DEM_VELOCITY>(i.get_data());
		Vect3d& v0 = std::get<DEM_VELOCITY0>(i.get_data());

		v = v0 + dt/2 * f;
	});
}

template<typename GeometryType>
void integrate_dem(const double dt,ptr<DemType> dem,
		ptr<Params> params,
		GeometryType geometry) {
	const int num_it = dt/params->dem_dt;
	for (int i = 0; i < num_it; ++i) {
		dem_start(dem,params,geometry);
		dem_end(dem,params,geometry);
	}
	const double save_dem_dt = params->dem_dt;
	params->dem_dt = (dt%params->dem_dt)*params->dem_dt;
	if (params->dem_dt>0) {
		dem_start(dem,params,geometry);
		dem_end(dem,params,geometry);
	}
}

template<typename SphGeometryType,typename DemGeometryType>
void sphdem_start(ptr<SphType> sph,ptr<DemType> dem,
		ptr<Params> params,
		SphGeometryType sph_geometry,DemGeometryType dem_geometry) {

	const double dt = params->sph_dt;


	sph->update_positions(sph->begin(),sph->end(),[dt](SphType::Value& i) {
		const Vect3d& r = i.get_position();
		Vect3d& f = std::get<SPH_FORCE>(i.get_data());
		Vect3d& v = std::get<SPH_VELOCITY>(i.get_data());
		Vect3d& v0 = std::get<SPH_VELOCITY0>(i.get_data());


		v0 = v + dt/2*f;
		v += dt * f;
		return r + dt/2 * v;
	});

	std::for_each(sph->begin(),sph->end(),[dem](SphType::Value& i) {
		const Vect3d& r = i.get_position();
		double& e = std::get<SPH_POROSITY>(i.get_data());
		double& h = std::get<SPH_H>(i.get_data());
		e = 0;

		for (auto tpl: i.get_neighbours(dem)) {
			const Vect3d& dx = std::get<1>(tpl);
			const double r = dx.norm();
			e += dem_vol*W(r/h,h);
		}
	});

	std::for_each(dem->begin(),dem->end(),[sph](DemType::Value& i) {
		const Vect3d& r = i.get_position();

		for (auto tpl: i.get_neighbours(sph)) {
			const Vect3d& dx = std::get<1>(tpl);
			const double r = dx.norm();

			//interpolate needed sph values here
		}
	});


	integrate_dem(dt,dem,params,dem_geometry);


	std::for_each(sph->begin(),sph->end(),[dem](SphType::Value& i) {
		const Vect3d& r = i.get_position();
		double& e = std::get<SPH_POROSITY>(i.get_data());
		double& h = std::get<SPH_H>(i.get_data());
		e = 0;

		for (auto tpl: i.get_neighbours(dem)) {
			const Vect3d& dx = std::get<1>(tpl);
			const double r = dx.norm();
			e += dem_vol*W(r/h,h);
		}
	});

	std::for_each(sph->begin(),sph->end(),[&sph_geometry](SphType::Value& i) {
			const Vect3d& r = i.get_position();
			Vect3d& f = std::get<SPH_FORCE>(i.get_data());
			Vect3d& v = std::get<SPH_VELOCITY>(i.get_data());
			double& e = std::get<SPH_POROSITY>(i.get_data());
			double& h = std::get<SPH_H>(i.get_data());
			double& pdr2 = std::get<SPH_PDR2>(i.get_data());
			double& omega = std::get<SPH_OMEGA>(i.get_data());

			f << 0,0,0;
			f = f + sph_geometry(i);

			for (auto tpl: i.get_in_radius(sph,2*h)) {
				const Vect3d& dx = std::get<1>(tpl);
				const SphType::Value& j = std::get<0>(tpl);
				if (i.get_id()==j.get_id()) continue;
				const Vect3d& vj = std::get<SPH_VELOCITY>(j.get_data());
				double& hj = std::get<SPH_H>(j.get_data());
				double& pdr2j = std::get<SPH_PDR2>(j.get_data());
				double& omegaj = std::get<SPH_OMEGA>(j.get_data());

				const double r = dx.norm();
				const Vect3d gradW = (1.0/omega)*dx*F(r/h,h);
				const Vect3d gradWj = (1.0/omegaj)*dx*F(r/hj,hj);
				f += -sph_mass*(pdr2*gradW + pdr2j*gradWj);
			}

		});

}

template<typename SphGeometryType,typename DemGeometryType>
void sphdem_end(ptr<SphType> sph,ptr<DemType> dem,
		ptr<Params> params,
		SphGeometryType sph_geometry,DemGeometryType dem_geometry) {

	const double dt = params->sph_dt;

	std::for_each(sph->begin(),sph->end(),[dt](SphType::Value& i) {
		Vect3d& f = std::get<SPH_FORCE>(i.get_data());
		Vect3d& v = std::get<SPH_VELOCITY>(i.get_data());
		Vect3d& v0 = std::get<SPH_VELOCITY0>(i.get_data());

		v = v0 + dt/2 * f;
	});
}




#endif /* SPHDEM_H_ */
