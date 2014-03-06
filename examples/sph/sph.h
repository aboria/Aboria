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

enum {SPH_FORCE, SPH_VELOCITY, SPH_VELOCITY0, SPH_DENS, SPH_POROSITY, SPH_H, SPH_DDDT,SPH_PDR2,SPH_OMEGA};
typedef std::tuple<Vect3d,Vect3d,Vect3d,double,double,double,double,double,double> SphTuple;
typedef Particles<SphTuple> SphType;

#define REGISTER_SPH_PARTICLE(particle) \
				const Vect3d& r = particle.get_position(); \
				GET_TUPLE(Vect3d,f,SPH_FORCE,particle); \
				GET_TUPLE(Vect3d,v,SPH_VELOCITY,particle); \
				GET_TUPLE(Vect3d,v0,SPH_VELOCITY0,particle); \
				GET_TUPLE(double,rho,SPH_DENS,particle); \
				GET_TUPLE(double,e,SPH_POROSITY,particle); \
				GET_TUPLE(double,h,SPH_H,particle); \
				GET_TUPLE(double,dddt,SPH_DDDT,particle); \
				GET_TUPLE(double,pdr2,SPH_PDR2,particle); \
				GET_TUPLE(double,omega,SPH_OMEGA,particle)

#define REGISTER_NEIGHBOUR_SPH_PARTICLE(tuple) \
				const Vect3d& dx = std::get<1>(tuple); \
				const DemType::Value& j = std::get<0>(tuple) \
				const Vect3d& rj = j.get_position(); \
				GET_TUPLE(Vect3d,fj,SPH_FORCE,j); \
				GET_TUPLE(Vect3d,vj,SPH_VELOCITY,j); \
				GET_TUPLE(Vect3d,v0j,SPH_VELOCITY0,j); \
				GET_TUPLE(double,rhoj,SPH_DENS,j); \
				GET_TUPLE(double,ej,SPH_POROSITY,j); \
				GET_TUPLE(double,hj,SPH_H,j); \
				GET_TUPLE(double,dddtj,SPH_DDDT,j); \
				GET_TUPLE(double,pdr2j,SPH_PDR2,j); \
				GET_TUPLE(double,omegaj,SPH_OMEGA,j)

struct Params {
	double sph_dt,dem_dt,dem_mass,dem_diameter,dem_k,dem_gamma;
};


template<typename SphGeometryType>
void sph_start(ptr<SphType> sph,
		ptr<Params> params,
		SphGeometryType sph_geometry) {

	const double dt = params->sph_dt;


	sph->update_positions(sph->begin(),sph->end(),[dt](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);

		v0 = v + dt/2*f;
		v += dt * f;
		return r + dt/2 * v;
	});


	std::for_each(sph->begin(),sph->end(),[&sph_geometry](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);

		f << 0,0,0;
		f = f + sph_geometry(i);

		for (auto tpl: i.get_in_radius(sph,2*h)) {
			REGISTER_NEIGHBOUR_SPH_PARTICLE(tpl);

			const double r = dx.norm();
			const Vect3d gradW = (1.0/omega)*dx*F(r/h,h);
			const Vect3d gradWj = (1.0/omegaj)*dx*F(r/hj,hj);
			f += -sph_mass*(pdr2*gradW + pdr2j*gradWj);
			
		}

	});

}

template<typename SphGeometryType>
void sph_end(ptr<SphType> sph,
		ptr<Params> params,
		SphGeometryType sph_geometry) {

	const double dt = params->sph_dt;

	std::for_each(sph->begin(),sph->end(),[dt](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);

		v = v0 + dt/2 * f;
	});
}




#endif /* SPHDEM_H_ */
