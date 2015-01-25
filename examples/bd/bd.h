/*
 * ld.h
 *
 *  Created on: 25 Jan 2015
 *      Author: mrobins
 */

#ifndef BD_H_
#define BD_H_

#include "Aboria.h"
using namespace Aboria;

#include <tuple>


enum {SPECIES_VELOCITY};
typedef std::tuple<Vect3d,double> SpeciesTuple;
typedef Particles<SpeciesTuple> SpeciesType;

#define GET_TUPLE(type,name,position,particle) type& name = std::get<position>(particle.get_data())
#define REGISTER_SPECIES_PARTICLE(particle) \
				const Vect3d& r = particle.get_position(); \
				const bool alive = particle.is_alive(); \
				GET_TUPLE(Vect3d,v,SPECIES_VELOCITY,particle);
#define REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tuple) \
				const Vect3d& dx = std::get<1>(tuple); \
				const SpeciesType::Value& j = std::get<0>(tuple); \
				const Vect3d& rj = j.get_position(); \
				const bool alivej = j.is_alive(); \
				const GET_TUPLE(Vect3d,vj,SPECIES_VELOCITY,j);

struct Params {
	Params() {};
	double D,dt;
	double time;
};


void bd_timestep(ptr<SpeciesType> A, ptr<Sphere> sphere,
		ptr<Params> params) {

	const double dt = params->dt;
	const double D = params->D;

	A->update_positions(A->begin(),A->end(),[A,dt,sphere,D](SpeciesType::Value& i) {
		REGISTER_SPECIES_PARTICLE(i);

		v = sqrt(2.0*D*dt)*Vect3d(i.rand_normal(),i.rand_normal(),i.rand_normal());

		Vect3d new_position = r+v;

		if (reflect_once(r,new_position,sphere)) {
			//reflected
		}

		return new_position;
	});

}


#endif /* BD_H_ */
