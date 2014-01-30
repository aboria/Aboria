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

namespace Aboria {

template<typename InputType, typename OutputType, typename FunctionType>
class ForEachParticle {
public:
	ForEachParticle(InputType input, OutputType output, FunctionType function):
		input(input),output(output),function(function) {}

	void print(std::ostream& out) const {
		out << "Particle for_each with function = " <<typeid(FunctionType).name();
	}
	void execute() {
		std::for_each(output.begin(),output.end(),function);
	}
	void reset() {}

private:
	InputType input;
	OutputType output;
	FunctionType function;
};

template<int N, typename FunctionType>
static std::shared_ptr<Operator> for_each(Particles<N>& output, FunctionType function) {
	return std::shared_ptr<Operator>(new Operator(ForEachParticle<int,Particles<N>,FunctionType>(0,output,function)));
}

template<typename InputType, int N, typename FunctionType>
static std::auto_ptr<Operator> for_each(InputType input, Particles<N> output, FunctionType function) {
	return std::auto_ptr<Operator>(new ForEachParticle<InputType,Particles<N>,FunctionType>(input,output,function));
}

} /* namespace Aboria */
#endif /* PARTICLEOPERATORS_H_ */
