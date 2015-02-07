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

template<typename ParticlesType>
class DataVectorBase {
public:
	DataVectorBase(ParticlesType& particles):
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

template<int I, typename ParticlesType>
class DataVector: public DataVectorBase<ParticlesType> {
public:
	typedef typename Elem<I, ParticlesType>::type value_type;

	DataVector(ParticlesType &particles):
		DataVectorBase<ParticlesType>(particles)
	{};
	value_type&
	operator []( std::size_t i ) {
		return this->particles[i].template get_elem<I>();
	}
	const value_type&
	operator []( std::size_t i ) const {
		return this->particles[i].template get_elem<I>();
	}
	void set( std::size_t i, const value_type& arg) {
		return this->particles[i].template set_elem<I>(arg);
	}
};

template<typename ParticlesType>
class DataVector<POSITION,ParticlesType>: public DataVectorBase<ParticlesType> {
public:
	typedef typename Elem<POSITION, ParticlesType>::type value_type;

	DataVector(ParticlesType &particles):
		DataVectorBase<ParticlesType>(particles)
	{};
	const value_type&
	operator []( std::size_t i ) const {
		return this->particles[i].get_position();
	}
	void set( std::size_t i, const value_type& arg) {
		return this->particles[i].set_position(arg);
	}
};

template<typename ParticlesType>
class DataVector<ID,ParticlesType>: public DataVectorBase<ParticlesType> {
public:
	typedef typename Elem<ID, ParticlesType>::type value_type;

	DataVector(ParticlesType &particles):
		DataVectorBase<ParticlesType>(particles)
	{};
	const value_type&
	operator []( std::size_t i ) const {
		return this->particles[i].get_id();
	}
};

template<typename ParticlesType>
class DataVector<ALIVE,ParticlesType>: public DataVectorBase<ParticlesType> {
public:
	typedef typename Elem<ALIVE, ParticlesType>::type value_type;

	DataVector(ParticlesType &particles):
		DataVectorBase<ParticlesType>(particles)
	{};
	const value_type&
	operator []( std::size_t i ) const {
		return this->particles[i].get_alive();
	}
	void set( std::size_t i, const value_type& arg) {
		return this->particles[i].set_alive(arg);
	}
};




}
#endif /* DATAVECTOR_H_ */
