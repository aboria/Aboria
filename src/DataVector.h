/*
 * DataVector.h
 *
 *  Created on: 23 Feb 2015
 *      Author: robinsonm
 */

#ifndef DATAVECTOR_H_
#define DATAVECTOR_H_

#include "Particles.h"

namespace Aboria {

template<typename I, typename ParticlesType>
class DataVector {
public:
	DataVector(ParticlesType& particles):
		particles(particles)
	{};
	ParticlesType &get_particles() {
		return particles;
	}
	const ParticlesType &get_particles() const {
		return particles;
	}
	std::size_t size() const {
		return particles.size();
	}
protected:
	ParticlesType &particles;
};


}
#endif /* DATAVECTOR_H_ */
