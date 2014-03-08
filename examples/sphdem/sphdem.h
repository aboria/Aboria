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
using namespace Aboria;

#include "sph_common.h"
#include <tuple>



enum {DEM_FORCE, DEM_VELOCITY, DEM_VELOCITY0, DEM_FORCE0,DEM_SHEPSUM};
typedef std::tuple<Vect3d,Vect3d,Vect3d,Vect3d,double> DemTuple;
typedef Particles<DemTuple> DemType;

#define GET_TUPLE(type,name,position,particle) type& name = std::get<position>(particle.get_data())
#define REGISTER_DEM_PARTICLE(particle) \
				const Vect3d& r = particle.get_position(); \
				GET_TUPLE(Vect3d,f,DEM_FORCE,particle); \
				GET_TUPLE(Vect3d,v,DEM_VELOCITY,particle); \
				GET_TUPLE(Vect3d,f0,DEM_FORCE0,particle); \
				GET_TUPLE(double,s,DEM_SHEPSUM,particle); \
				GET_TUPLE(Vect3d,v0,DEM_VELOCITY0,particle)
#define REGISTER_NEIGHBOUR_DEM_PARTICLE(tuple) \
				const Vect3d& dx = std::get<1>(tuple); \
				const DemType::Value& j = std::get<0>(tuple); \
				const Vect3d& rj = j.get_position(); \
				const GET_TUPLE(Vect3d,fj,DEM_FORCE,j); \
				const GET_TUPLE(Vect3d,vj,DEM_VELOCITY,j); \
				const GET_TUPLE(Vect3d,f0j,DEM_FORCE0,j); \
				const GET_TUPLE(double,sj,DEM_SHEPSUM,j); \
				const GET_TUPLE(Vect3d,v0j,DEM_VELOCITY0,j)


enum {SPH_FORCE, SPH_VELOCITY, SPH_VELOCITY0,SPH_FORCE0,SPH_DENS, SPH_POROSITY, SPH_H, SPH_DDDT,SPH_PDR2,SPH_OMEGA,SPH_FIXED,SPH_FORCE_EXT};
typedef std::tuple<Vect3d,Vect3d,Vect3d,Vect3d,double,double,double,double,double,double,bool,Vect3d> SphTuple;
typedef Particles<SphTuple> SphType;

#define REGISTER_SPH_PARTICLE(particle) \
				const Vect3d& r = particle.get_position(); \
				GET_TUPLE(bool,fixed,SPH_FIXED,particle); \
				GET_TUPLE(Vect3d,f,SPH_FORCE,particle); \
				GET_TUPLE(Vect3d,v,SPH_VELOCITY,particle); \
				GET_TUPLE(Vect3d,v0,SPH_VELOCITY0,particle); \
				GET_TUPLE(double,rho,SPH_DENS,particle); \
				GET_TUPLE(double,e,SPH_POROSITY,particle); \
				GET_TUPLE(double,h,SPH_H,particle); \
				GET_TUPLE(double,dddt,SPH_DDDT,particle); \
				GET_TUPLE(double,pdr2,SPH_PDR2,particle); \
				GET_TUPLE(Vect3d,f0,SPH_FORCE0,particle); \
				GET_TUPLE(Vect3d,fext,SPH_FORCE_EXT,particle); \
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
				const GET_TUPLE(double,ej,SPH_POROSITY,j); \
				const GET_TUPLE(double,hj,SPH_H,j); \
				const GET_TUPLE(double,dddtj,SPH_DDDT,j); \
				const GET_TUPLE(double,pdr2j,SPH_PDR2,j); \
				const GET_TUPLE(Vect3d,f0j,SPH_FORCE0,j); \
				const GET_TUPLE(Vect3d,fextj,SPH_FORCE_EXT,j); \
				const GET_TUPLE(double,omegaj,SPH_OMEGA,j)

struct Params {
	double sph_dt,sph_mass,sph_hfac,sph_visc,sph_refd,sph_gamma,sph_spsound,sph_prb,sph_dens;
	double dem_dt,dem_mass,dem_diameter,dem_k,dem_gamma, dem_vol, dem_time_drop;
	double time;
};


template<typename GeometryType>
void dem_timestep(ptr<DemType> dem,
		ptr<Params> params,
		GeometryType geometry) {

	const double dt = params->dem_dt;
	const double dem_diameter = params->dem_diameter;
	const double dem_k = params->dem_k;
	const double dem_gamma = params->dem_gamma;
	const double dem_mass = params->dem_mass;
	const double dem_vol = params->dem_vol;

	dem->update_positions(dem->begin(),dem->end(),[dt](DemType::Value& i) {
		REGISTER_DEM_PARTICLE(i);

		v0 = v + dt/2*(f+f0);
		v += dt * (f+f0);
		return r + dt * v0;
	});

	std::for_each(dem->begin(),dem->end(),[&geometry,dem,dem_k,dem_gamma,dem_mass,dem_diameter](DemType::Value& i) {
		REGISTER_DEM_PARTICLE(i);

		f << 0,0,0;
		f = f + geometry(i);

		for (auto tpl: i.get_neighbours(dem)) {
			REGISTER_NEIGHBOUR_DEM_PARTICLE(tpl);

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

	std::for_each(dem->begin(),dem->end(),[dt](DemType::Value& i) {
		REGISTER_DEM_PARTICLE(i);

		v = v0 + dt/2 * (f+f0);
	});

}


template<typename GeometryType>
void integrate_dem(const double dt,ptr<DemType> dem,
		ptr<Params> params,
		GeometryType geometry) {
	if (params->time < params->dem_time_drop) return;

	const int num_it = floor(dt/params->dem_dt);
	//std::cout << "going for "<<num_it<<" dem timesteps"<<std::endl;
	for (int i = 0; i < num_it; ++i) {
		dem_timestep(dem,params,geometry);
	}
	const double save_dem_dt = params->dem_dt;
	params->dem_dt = dt - num_it*params->dem_dt;
	if (params->dem_dt>0) {
		dem_timestep(dem,params,geometry);
	}
	params->dem_dt = save_dem_dt;
}

template<typename SphGeometryType,typename DemGeometryType>
void sphdem(ptr<SphType> sph,ptr<DemType> dem,
		ptr<Params> params,
		SphGeometryType sph_geometry,DemGeometryType dem_geometry) {

	const double dt = params->sph_dt;
	const double sph_mass = params->sph_mass;
	const double sph_prb = params->sph_prb;
	const double sph_refd = params->sph_refd;
	const double sph_dens = params->sph_dens;
	const double sph_gamma = params->sph_gamma;
	const double sph_visc = params->sph_visc;

	const double dem_vol = params->dem_vol;
	const double dem_diameter = params->dem_diameter;
	const double dem_mass = params->dem_mass;

	/*
	 * 0 -> 1/2 step
	 */
	//std::cout << "0 -> 1/2 step"<<std::endl;
	integrate_dem(dt/2,dem,params,dem_geometry); //update dem positions

	sph->update_positions(sph->begin(),sph->end(),[dt,sph_mass](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);
		if (!fixed) v += dt/2 * (f+f0+fext);
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
	 *  Calculate coupling force on DEM
	 */
	//std::cout << "calculate coupling force on DEM"<<std::endl;

	std::for_each(dem->begin(),dem->end(),[sph,dem_vol,sph_mass,sph_visc,dem_diameter,sph_dens,dem_mass](DemType::Value& i) {
		REGISTER_DEM_PARTICLE(i);

		f0 << 0,0,0;
		s = 0;
		Vect3d vf(0,0,0);
		double ef = 0;
		for (auto tpl: i.get_neighbours(sph)) {
			REGISTER_NEIGHBOUR_SPH_PARTICLE(tpl);
			const double r2 = dx.squaredNorm();
			if (r2 > 4.0*hj*hj) continue;
			const double r = sqrt(r2);
			const double dvWab = sph_mass*W(r/hj,hj)/rhoj;
			f0 += fj*dvWab;
			vf += vj*dvWab;
			ef += ej*dvWab;
			s += dvWab;
		}
		f0 /= s;
		vf /= s;
		ef /= s;

		if (s > 0.5) {
			const Vect3d fdrag = 3.0*PI*sph_visc*sph_dens*dem_diameter*ef*(vf-v);
			//std::cout <<"v = "<<v<<" vf = "<<vf<<" f0 = "<<f0<<" fdrag = "<<fdrag<<std::endl;

			f0 = (dem_vol*f0 + fdrag)/dem_mass;
		} else {
			f0 << 0,0,0;
		}
	});

	/*
	 * 1/2 -> 1 step
	 */
	//std::cout << "1/2 -> 1 step"<<std::endl;

	integrate_dem(dt/2,dem,params,dem_geometry); //update dem positions
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
		if (!fixed) v += dt/2 * (f+f0+fext);
		const double press = sph_prb*(pow(rho/(e*sph_refd),sph_gamma) - 1.0);
		pdr2 = press/pow(rho,2);
		return r + dt/2 * v0;
	});
	dem->reset_neighbour_search(2.0*maxh);


	/*
	 * Calculate porosity on SPH
	 */
	//std::cout << "calcaulte porosity on sph"<<std::endl;

	std::for_each(sph->begin(),sph->end(),[dem,dem_vol](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);
		e=1;
		for (auto tpl: i.get_neighbours(dem)) {
			REGISTER_NEIGHBOUR_DEM_PARTICLE(tpl);
			const double r2 = dx.squaredNorm();
			if (r2 > 4.0*h*h) continue;
			const double r = sqrt(r2);
			e -= dem_vol*W(r/h,h);
		}
	});

	/*
	 *  Calculate shep sum on DEM
	 */
	//std::cout << "calcaulte shep sum on dem"<<std::endl;

	std::for_each(dem->begin(),dem->end(),[sph,dem_vol,sph_mass](DemType::Value& i) {
		REGISTER_DEM_PARTICLE(i);
		s = 0;
		for (auto tpl: i.get_neighbours(sph)) {
			REGISTER_NEIGHBOUR_SPH_PARTICLE(tpl);
			const double r2 = dx.squaredNorm();
			if (r2 > 4.0*hj*hj) continue;
			const double r = sqrt(r2);
			const double Wab = W(r/hj,hj);
			s += sph_mass*Wab/rhoj;
		}
	});
	
	/*
	 * Calculate coupling force on SPH
	 */
	//std::cout << "calculate coupling force on sph"<<std::endl;

	std::for_each(sph->begin(),sph->end(),[dem,dem_mass](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);

		f0 << 0,0,0;
		for (auto tpl: i.get_neighbours(dem)) {
			REGISTER_NEIGHBOUR_DEM_PARTICLE(tpl);
			const double r2 = dx.squaredNorm();
			if (r2 > 4.0*h*h) continue;
			const double r = sqrt(r2);
			f0 -= dem_mass*f0j*W(r/h,h)/sj;
		}
		f0 /= rho;
	});

	dem->reset_neighbour_search(dem_diameter);

	/*
	 * acceleration on SPH calculation
	 */
	//std::cout << "acceleration on SPH"<<std::endl;

	std::for_each(sph->begin(),sph->end(),[sph,&sph_geometry,sph_mass,sph_visc](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);

		f << 0,0,0;
		fext = sph_geometry(i);
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

		if (!fixed) v = v0 + dt/2 * (f+f0+fext);
	});

	/*
	 * sph timestep condition
	 */
	params->time += params->sph_dt;
	params->sph_dt = std::min(0.25*minh/params->sph_spsound,0.125*pow(minh,2)/params->sph_visc);
}


#endif /* SPHDEM_H_ */
