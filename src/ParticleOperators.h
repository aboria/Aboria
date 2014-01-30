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
#include <memory>
//#include <parallel/algorithm>

namespace Aboria {

template<int N, typename FunctionType>
class ForEachParticle {
public:
	ForEachParticle(Particles<N>& output, FunctionType function, const std::string name):
		output(output),function(function),name(name) {}

	void print(std::ostream& out) const {
		out << "Particle for_each with function = " <<name;
	}
	void execute() {
		std::for_each(output.begin(),output.end(),function);
	}
	void reset() {}

private:
	const std::string name;
	Particles<N>& output;
	FunctionType function;
};

template<int N, typename FunctionType>
static std::shared_ptr<Operator> for_each(Particles<N>& output, FunctionType function) {
	return std::shared_ptr<Operator>(new Operator(
			ForEachParticle<N,FunctionType>(output,function,typeid(function).name())
			));
}

} /* namespace Aboria */
#endif /* PARTICLEOPERATORS_H_ */
