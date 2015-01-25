/*
 * bd.cpp
 *
 *  Created on: 25 Jan 2015
 *      Author: mrobins
 */

#ifndef BD_CPP_
#define BD_CPP_

#include "bd.h"
#include <assert.h>     /* assert */


int main(int argc, char **argv) {

	auto A = SpeciesType::New();
	auto params = ptr<Params>(new Params());

	const double radius = 1.0;
	auto sphere = Sphere::New(Vect3d(0,0,0),radius,false);


	params->D = 1.0;
	params->dt = 0.0001;

	const double buffer = 0.1;
	const Vect3d min(-radius-buffer,-radius-buffer,-radius-buffer);
	const Vect3d max(radius+buffer,radius+buffer,radius+buffer);
	const Vect3b periodic(false,false,false);

	SpeciesType::value_type p;
	const int N = 100;
	for (int i = 0; i < N; ++i) {
		p.set_position(Vect3d(0,0,0));
		A->push_back(p);
	}

	const int timesteps = 1000;
	for (int i = 0; i < timesteps; ++i) {
		bd_timestep(A,sphere,params);
	}

	std::for_each(A->begin(),A->end(),[radius](SpeciesType::Value& i) {
		assert(i.get_position().norm() < radius);
	});


}


#endif /* BD_CPP_ */
