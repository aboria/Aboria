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
#include <tuple>
using namespace Aboria;
#include "sph_common.h"

enum {SPH_FORCE, SPH_VELOCITY, SPH_VELOCITY0,SPH_DENS,SPH_H, SPH_DDDT,SPH_PDR2,SPH_OMEGA,SPH_FIXED};
typedef std::tuple<Vect3d,Vect3d,Vect3d,double,double,double,double,double,bool> SphTuple;
typedef Particles<SphTuple> SphType;

#define GET_TUPLE(type,name,position,particle) type& name = std::get<position>(particle.get_data())
#define REGISTER_SPH_PARTICLE(particle) \
				const Vect3d& r = particle.get_position(); \
				GET_TUPLE(bool,fixed,SPH_FIXED,particle); \
				GET_TUPLE(Vect3d,f,SPH_FORCE,particle); \
				GET_TUPLE(Vect3d,v,SPH_VELOCITY,particle); \
				GET_TUPLE(Vect3d,v0,SPH_VELOCITY0,particle); \
				GET_TUPLE(double,rho,SPH_DENS,particle); \
				GET_TUPLE(double,h,SPH_H,particle); \
				GET_TUPLE(double,dddt,SPH_DDDT,particle); \
				GET_TUPLE(double,pdr2,SPH_PDR2,particle); \
				GET_TUPLE(double,omega,SPH_OMEGA,particle)

#define REGISTER_NEIGHBOUR_SPH_PARTICLE(tuple) \
				const Vect3d& dx = std::get<1>(tuple); \
				const SphType::Value& j = std::get<0>(tuple); \
				const Vect3d& rj = j.get_position(); \
				const GET_TUPLE(bool,fixedj,SPH_FIXED,j); \
				const GET_TUPLE(Vect3d,fj,SPH_FORCE,j); \
				const GET_TUPLE(Vect3d,vj,SPH_VELOCITY,j); \
				const GET_TUPLE(Vect3d,v0j,SPH_VELOCITY0,j); \
				const GET_TUPLE(double,rhoj,SPH_DENS,j); \
				const GET_TUPLE(double,hj,SPH_H,j); \
				const GET_TUPLE(double,dddtj,SPH_DDDT,j); \
				const GET_TUPLE(double,pdr2j,SPH_PDR2,j); \
				const GET_TUPLE(double,omegaj,SPH_OMEGA,j)

struct Params {
	double sph_dt,sph_mass,sph_hfac,sph_visc,sph_refd,sph_gamma,sph_spsound,sph_prb,sph_dens,sph_time_damping;
	double time;
};


template<typename SphGeometryType>
void sph_timestep(ptr<SphType> sph,
		ptr<Params> params,
		SphGeometryType sph_geometry) {

	const double dt = params->sph_dt;
	const double sph_mass = params->sph_mass;
	const double sph_prb = params->sph_prb;
	const double sph_refd = params->sph_refd;
	const double sph_dens = params->sph_dens;
	const double sph_gamma = params->sph_gamma;
	const double sph_visc = params->sph_visc;


	/*
	 * 0 -> 1/2 step
	 */
	//std::cout << "0 -> 1/2 step"<<std::endl;

	sph->update_positions(sph->begin(),sph->end(),[params,dt,sph_mass](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);
		if (!fixed) v += dt/2 * f;
		if (params->time < params->sph_time_damping) v *= 0.98;
		return r + dt/2 * v;
	});


	/*
	 * Calculate omega
	 */
	//std::cout << "calculate omega"<<std::endl;

	std::for_each(sph->begin(),sph->end(),[sph,sph_mass](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);
		omega = -sph_mass*(0+NDIM/pow(h,NDIM)*K(0,h));
		for (auto tpl: i.get_neighbours(sph)) {
			REGISTER_NEIGHBOUR_SPH_PARTICLE(tpl);
			const double r2 = dx.squaredNorm();
			if (r2 > 4.0*h*h) continue;
			if (r2 == 0) continue;
			const double r = sqrt(r2);
			const double q = r/h;
			omega -= sph_mass*(r2*F(q,h)+NDIM/pow(h,NDIM)*K(q,h));
		}
		omega = 1.0 + omega/(rho*NDIM);
		//std::cout << "omega = "<<omega<<" rho = "<<rho<<" ndim = "<<NDIM<<std::endl;
	});

	/*
	 * Calculate change in density
	 */
	//std::cout << "calculate change in density"<<std::endl;

	std::for_each(sph->begin(),sph->end(),[sph,&sph_geometry,sph_mass](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);

		dddt = 0;
		for (auto tpl: i.get_neighbours(sph)) {
			REGISTER_NEIGHBOUR_SPH_PARTICLE(tpl);
			const double r2 = dx.squaredNorm();
			if (r2 > 4.0*h*h) continue;
			if (r2 == 0) continue;
			const double r = sqrt(r2);

			dddt += sph_mass*(v-vj).dot(dx*F(r/h,h));
		}
		dddt *= 1.0/omega;

	});


	/*
	 * 1/2 -> 1 step
	 */
	//std::cout << "1/2 -> 1 step"<<std::endl;

	std::for_each(sph->begin(),sph->end(),[dt,sph_mass](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);
		rho += dt * dddt;
		h = pow(sph_mass/rho,1.0/NDIM);
	});
	auto iterator_to_maxh =
			std::max_element(sph->begin(),sph->end(),[](SphType::Value& i, SphType::Value& j){
		GET_TUPLE(double,hi,SPH_H,i);
		GET_TUPLE(double,hj,SPH_H,j);
		return hi<hj;
	});
	const double maxh = std::get<SPH_H>(iterator_to_maxh->get_data());
	auto iterator_to_minh =
			std::min_element(sph->begin(),sph->end(),[](SphType::Value& i, SphType::Value& j){
		GET_TUPLE(double,hi,SPH_H,i);
		GET_TUPLE(double,hj,SPH_H,j);
		return hi<hj;
	});
	const double minh = std::get<SPH_H>(iterator_to_minh->get_data());

	sph->reset_neighbour_search(2.0*maxh,[dt,sph_mass,sph_prb,sph_refd,sph_gamma](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);
		v0 = v;
		if (!fixed) v += dt/2 * f;
		const double press = sph_prb*(pow(rho/sph_refd,sph_gamma) - 1.0);
		pdr2 = press/pow(rho,2);
		return r + dt/2 * v0;
	});


	/*
	 * acceleration on SPH calculation
	 */
	//std::cout << "acceleration on SPH"<<std::endl;

	std::for_each(sph->begin(),sph->end(),[sph,&sph_geometry,sph_mass,sph_visc](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);

		f = Vect3d(0,0,0);
		f += sph_geometry(i);

		for (auto tpl: i.get_neighbours(sph)) {
			REGISTER_NEIGHBOUR_SPH_PARTICLE(tpl);
			const double r2 = dx.squaredNorm();
			if (r2 > 4.0*h*h) continue;
			if (r2 == 0) continue;
			const double r = sqrt(r2);
			const double fdash = F(r/h,h);
			const double fdashj = F(r/hj,hj);

			/*
			 * pressure gradient
			 */
			f += -sph_mass*((1.0/omega)*pdr2*fdash + (1.0/omegaj)*pdr2j*fdashj)*dx;

			const Vect3d dv = v-vj;
			const double visc = sph_visc*(rho+rhoj)/(rho*rhoj);

			/*
			 * viscosity (morris)
			 */
			f += dv*sph_mass*visc*fdash;
		}

	});


	/*
	 * 1/2 -> 1 step for velocity
	 */
	//std::cout << "1/2 -> 1 step for velocity"<<std::endl;

	std::for_each(sph->begin(),sph->end(),[dt](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);

		if (!fixed) v = v0 + dt/2 * f;
	});

	/*
	 * sph timestep condition
	 */
	params->time += params->sph_dt;
	params->sph_dt = std::min(0.25*minh/params->sph_spsound,0.125*pow(minh,2)/params->sph_visc);

}

#endif /* SPHDEM_H_ */
