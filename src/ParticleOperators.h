/*
 * ParticleOperators.h
 *
 *  Created on: 29 Jan 2014
 *      Author: mrobins
 */

#ifndef PARTICLEOPERATORS_H_
#define PARTICLEOPERATORS_H_

#include "Particles.h"
#include "Operator.h"
#include "Ptr.h"
#include <memory>
//#include <parallel/algorithm>

namespace Aboria {

template<typename ParticlesType, typename FunctionType>
class ForEachParticle {
public:
	ForEachParticle(ParticlesType output, FunctionType function, const std::string name):
		output(output),function(function),name(name) {}

	void print(std::ostream& out) const {
		out << "Particle for_each with function = " <<name;
	}
	void execute() {
		std::for_each(output->begin(),output->end(),function);
	}
	void reset() {}

private:
	const std::string name;
	ParticlesType output;
	FunctionType function;
};


template<typename ParticlesType, typename FunctionType>
static std::shared_ptr<Operator> for_each(ParticlesType output, FunctionType function) {
	return std::shared_ptr<Operator>(new Operator(
			ForEachParticle<ParticlesType,FunctionType>(output,function,typeid(function).name())
			));
}

template<typename ParticlesType, typename FunctionType>
static std::shared_ptr<Operator> for_each(ParticlesType output, FunctionType function, std::string name) {
	return std::shared_ptr<Operator>(new Operator(
			ForEachParticle<ParticlesType,FunctionType>(output,function,name)
	));
}

template<typename ParticlesType, typename ParticlesTypeNeighbours, typename FunctionType>
class ForEachParticleNeighbourSearch {
public:
	ForEachParticleNeighbourSearch(ParticlesType output, ParticlesTypeNeighbours neighbours, FunctionType function, const std::string name):
		output(output),neighbours(neighbours),function(function),name(name) {}

	void print(std::ostream& out) const {
		out << "Particle for_each_neighbours with function = " <<name;
	}
	void execute() {
		neighbours->refresh_neighbour_search();
		std::for_each(output->begin(),output->end(),function);
	}
	void reset() {}

private:
	const std::string name;
	ParticlesType output;
	ParticlesTypeNeighbours neighbours;
	FunctionType function;
};

template<typename ParticlesType, typename ParticlesTypeNeighbours, typename FunctionType>
static ptr<Operator> for_each_neigh(ParticlesType output, ParticlesTypeNeighbours neighbours, FunctionType function) {
	return ptr<Operator>(new Operator(
			ForEachParticleNeighbourSearch<ParticlesType,FunctionType>(output,neighbours,function,typeid(function).name())
	));
}

template<typename ParticlesType, typename ParticlesTypeNeighbours, typename FunctionType>
static ptr<Operator> for_each_neigh(ParticlesType output, ParticlesTypeNeighbours neighbours, FunctionType function, const std::string name) {
	return ptr<Operator>(new Operator(
			ForEachParticleNeighbourSearch<ParticlesType,FunctionType>(output,neighbours,function,name)
	));
}

} /* namespace Aboria */
#endif /* PARTICLEOPERATORS_H_ */
